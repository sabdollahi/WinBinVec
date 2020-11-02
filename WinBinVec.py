import torch.nn as nn
import torch.nn.functional as F
import os
import pickle
import pandas as pd
from sklearn.utils import shuffle
from sklearn.model_selection import KFold
from torch.autograd import Variable
import torch
import math
import numpy as np

class WinBinVec(nn.Module):
    def __init__(self):
        super(WinBinVec, self).__init__()
        self.fc1 = torch.nn.Linear(430, 256)
        self.fc2 = torch.nn.Linear(256, 430)
        #Kernel size and stride = Window size
        #new_dim = math.floor((dim - kernel + (2*padding))/stride) + 1
        self.conv = torch.nn.Conv1d(430, 430, kernel_size=5, stride=5, padding=0)
        #new_dim*430 = 2150
        self.fc3 = torch.nn.Linear(2150, 1024)
        self.fc4 = torch.nn.Linear(1024, 512)
        self.fc5 = torch.nn.Linear(512, 256)
        self.fc6 = torch.nn.Linear(256, 64)
        self.fc7 = torch.nn.Linear(64, 32)
        self.fc8 = torch.nn.Linear(32, 2)
        
        self.bn3 = nn.BatchNorm1d(num_features=1024)
        self.drop3 = torch.nn.Dropout(0.5)
        self.bn4 = nn.BatchNorm1d(num_features=512)
        self.drop4 = torch.nn.Dropout(0.5)
        self.bn5 = nn.BatchNorm1d(num_features=256)
        self.drop5 = torch.nn.Dropout(0.5)

    def forward(self, X):
        #Number of PPIs = 430
        #Number of Amino Acids in PPIs = 25	
        #(batch_size, 430, 25)
        #Firstly, we prefer to use an one-dimensional convolutional layer to extract useful features of each Window
        #Since Window_size=5, we set kernel (filter) = 5 and stride = 5
        X = self.conv(X)
        #For each PPI mutation binary vector, we calculate its mean
        mX = torch.mean(X, dim=2)
        mX = F.relu(self.fc1(mX))
        #It assigns a kind of weight to each PPI (Sigmoid)
        mX = F.sigmoid(self.fc2(mX))
        #We mutiply these weights to their corresponding PPI 
        mX = mX.view(X.size(0), -1, 1).expand_as(X)
        X = X * mX
        X = X.view(X.size(0),-1)
        X = self.drop3(self.bn3(self.fc3(X)))
        X = self.drop4(self.bn4(self.fc4(X)))
        X = self.drop5(self.bn5(self.fc5(X)))
        X = self.fc6(X)
        X = self.fc7(X)
        return self.fc8(X)

#Each Cancer Organ Prediction
main_folder = "MutSigPPI/"
tcga_clinical_dataframe = pickle.load(open(main_folder + "TCGA_clinical_dataframe.pickle","rb"))
classes = {"adrenal gland":0, "bladder":1, "breast":2, "GYN":3, "bile duct":4, "CRC":5, "bone marrow":6, "esophagus":7, "brain":8, "head and neck":9, "Kidney":10, "liver":11, "Lung":12, "pleura":13, "pancreatic":14, "male reproductive system":15, "other":16, "Melanoma":17, "stomach":18, "thyroid":19, "thymus":20}
which_clinicals = ['cancer_class']
tcga_clinical_dataframe = tcga_clinical_dataframe[which_clinicals]
for cancer_class in classes:
    folds_accuracy = []
    replace_statement = {}
    for cl in classes:
        if(cl != cancer_class):
            replace_statement[cl] = 0
        else:
            replace_statement[cl] = 1
    specific_cancer_patients = tcga_clinical_dataframe[tcga_clinical_dataframe["cancer_class"] == cancer_class]
    specific_cancer_patients = specific_cancer_patients.replace({'cancer_class': replace_statement})
    other_cancer_patients = tcga_clinical_dataframe[tcga_clinical_dataframe["cancer_class"] != cancer_class]
    other_cancer_patients = shuffle(other_cancer_patients).sample(n = len(specific_cancer_patients))
    other_cancer_patients = other_cancer_patients.replace({'cancer_class': replace_statement})
    K = 10 #Kfold (number of parts = K)
    kf_other = KFold(n_splits = K, shuffle = True)
    kf_specific = KFold(n_splits = K, shuffle = True)
    parts_specific = kf_specific.split(specific_cancer_patients)
    parts_other = kf_other.split(other_cancer_patients)
    indices_specific = next(parts_specific, None)
    indices_other = next(parts_other, None)

    while(indices_specific):
        lr = 0.5
        #Define the model
        model = WinBinVec()
        criterion = nn.CrossEntropyLoss()
        optimizer = torch.optim.SGD(model.parameters(), lr=lr)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, 1.0, gamma=0.95)
        batch_size = 20
        for shuffled_epoch in range(10):
            print("Cancer Class:" + cancer_class + " Shuffled Epoch=" + str(shuffled_epoch))
            training = specific_cancer_patients.iloc[indices_specific[0]]
            training_other = other_cancer_patients.iloc[indices_other[0]]
            training = shuffle(training.append(training_other))
            Y = training[['cancer_class']].values
            Y = Variable(torch.LongTensor(Y.flatten()), requires_grad=False)
            training = training.index
            for epoch in range(50):
                for index in range(0, len(training), batch_size):
                    y = Y[index : index + batch_size]
                    batch_X = []
                    for patient in training[index : index + batch_size]:
                        p_data = pickle.load(open("MutSigPPI/DATASET/TheMostImportantWindows/" + patient + "_ppi.pickle", "rb"))
                        batch_X.append(p_data)
                    X = torch.FloatTensor(batch_X)
                    optimizer.zero_grad()
                    Y_hat = model(X)
                    loss = criterion(Y_hat, y)
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 0.5)
                    optimizer.step()
        test = specific_cancer_patients.iloc[indices_specific[1]]
        test_other = other_cancer_patients.iloc[indices_other[1]]
        test = shuffle(test.append(test_other))
        Y_test = test[['cancer_class']].values
        Y_test = Variable(torch.LongTensor(Y_test.flatten()), requires_grad=False)
        test = test.index
        test_list = []
        for patient in test:
            p_data = pickle.load(open("MutSigPPI/DATASET/TheMostImportantWindows/" + patient + "_ppi.pickle", "rb"))
            test_list.append(p_data)
        test_list = torch.FloatTensor(test_list)
        test_batch_Y_hat = model.forward(test_list)
        dummy, preds_test = torch.max (test_batch_Y_hat, dim = 1)
        print(preds_test)
        accuracy_test = (preds_test == Y_test).long().sum().float() /  preds_test.size()[0]
        print("ACC: " + str(accuracy_test))
        folds_accuracy.append(accuracy_test)
        indices_specific = next(parts_specific, None)
        indices_other = next(parts_other, None)
    pickle.dump(folds_accuracy, open("MutSigPPI/RESULTS/CancerClass_Accuracies_TheMostImportant/" + cancer_class + "_folds_accuracy.pickle","wb"))



#Predict Metastasis (Stage IV) or not (Stages I, II, and III)
#tcga_clinical_dataframe[tcga_clinical_dataframe['stage'] == 'Stage IVA']
which_clinicals = ['stage']
tcga_clinical_dataframe = tcga_clinical_dataframe[which_clinicals]
replace_statement = {}
metastasis_list = ['Stage IV','Stage IVA','Stage IVB','Stage IVC']
other_list = ['Stage I','Stage IA','Stage IB','Stage II','Stage IIA','Stage IIB','Stage IIC','Stage III','Stage IIIA','Stage IIIB','Stage IIIC']
#Metastasis Stage
for m in metastasis_list:
    replace_statement[m] = 1
#Non-metastasis Stage
for o in other_list:
    replace_statement[o] = 0

metastasis_patients = tcga_clinical_dataframe[tcga_clinical_dataframe["stage"].isin(metastasis_list)]
metastasis_patients = metastasis_patients.replace({'stage': replace_statement})
other_patients = tcga_clinical_dataframe[tcga_clinical_dataframe["stage"].isin(other_list)]
other_patients = other_patients.replace({'stage': replace_statement})

start_and_end_for_other = [0,793,1586,2379,3172,3965,4758,5554]
for i in range(7):
    print("PART: " + str(i))
    selected_other_patients = other_patients[start_and_end_for_other[i]:start_and_end_for_other[i+1]]
    folds_accuracy = []
    K = 10 #Kfold (number of parts = K)
    kf_other = KFold(n_splits = K, shuffle = True)
    kf_metastasis = KFold(n_splits = K, shuffle = True)
    parts_metastasis = kf_metastasis.split(metastasis_patients)
    parts_other = kf_other.split(selected_other_patients)
    indices_metastasis = next(parts_metastasis, None)
    indices_other = next(parts_other, None)
    fold_number = 1
    while(indices_metastasis):
        lr = 0.5
        #Define the model
        model = WinBinVec()
        criterion = nn.CrossEntropyLoss()
        optimizer = torch.optim.SGD(model.parameters(), lr=lr)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, 1.0, gamma=0.95)
        batch_size = 20
        for shuffled_epoch in range(10):
            print("Shuffled Epoch= " + str(shuffled_epoch))
            training = metastasis_patients.iloc[indices_metastasis[0]]
            training_other = selected_other_patients.iloc[indices_other[0]]
            training = shuffle(training.append(training_other))
            Y = training[['stage']].values
            Y = Variable(torch.LongTensor(Y.flatten()), requires_grad=False)
            training = training.index
            for epoch in range(50):
                for index in range(0, len(training), batch_size):
                    y = Y[index : index + batch_size]
                    batch_X = []
                    for patient in training[index : index + batch_size]:
                        p_data = pickle.load(open("MutSigPPI/DATASET/TheMostImportantWindows/" + patient + "_ppi.pickle", "rb"))
                        batch_X.append(p_data)
                    X = torch.FloatTensor(batch_X)
                    optimizer.zero_grad()
                    Y_hat = model(X)
                    loss = criterion(Y_hat, y)
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 0.5)
                    optimizer.step()
        test = metastasis_patients.iloc[indices_metastasis[1]]
        test_other = selected_other_patients.iloc[indices_other[1]]
        test = shuffle(test.append(test_other))
        Y_test = test[['stage']].values
        Y_test = Variable(torch.LongTensor(Y_test.flatten()), requires_grad=False)
        test = test.index
        test_list = []
        for patient in test:
            p_data = pickle.load(open("MutSigPPI/DATASET/TheMostImportantWindows/" + patient + "_ppi.pickle", "rb"))
            test_list.append(p_data)
        test_list = torch.FloatTensor(test_list)
        test_batch_Y_hat = model.forward(test_list)
        dummy, preds_test = torch.max (test_batch_Y_hat, dim = 1)
        accuracy_test = (preds_test == Y_test).long().sum().float() /  preds_test.size()[0]
        print("Fold: " + str(fold_number) + " ACC: " + str(accuracy_test))
        folds_accuracy.append(accuracy_test)
        indices_metastasis = next(parts_metastasis, None)
        indices_other = next(parts_other, None)
    pickle.dump(folds_accuracy, open("MutSigPPI/RESULTS/StagePrediction/Part" + str(i) + "_folds_accuracy.pickle","wb"))
    print("-----------------------------")

