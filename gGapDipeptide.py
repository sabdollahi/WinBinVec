import pickle
import pandas as pd
from sklearn.utils import shuffle
from sklearn.model_selection import KFold
import numpy as np
import torch.nn.functional as F
from torch.autograd import Variable
import torch
import torch.nn as nn

class gGapConv(nn.Sequential):
    def __init__(self):
        super(gGapConv, self).__init__()
        nn.Conv2d(1, 1, 3, stride=3) 
        #3, 1, 430, 430 -> 3, 1, 144, 144 (Conv2d)
        #3, 1, 144, 144 -> 3, 1, 72, 72 (MaxPool2d)
        #n by n * f by f -> floor[((n + 2p - f)/s) + 1] by floor[((n + 2p - f)/s) + 1]
        self.layer1 = nn.Sequential(
            nn.Conv2d(1, 1, kernel_size=3, stride=3, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2))
        #3, 1, 72, 72 -> 3, 1, 12, 12 (Conv2d and MaxPool2d)
        self.layer2 = nn.Sequential(
            nn.Conv2d(1, 1, kernel_size=3, stride=3, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2))
        self.fc1 = nn.Linear(12 * 12, 64)
        self.fc2 = nn.Linear(64, 16)
        self.fc3 = nn.Linear(16, 2)
        #Batch Normalization and Dropout Layers
        self.bn1 = nn.BatchNorm1d(num_features=64)
        self.drop1 = torch.nn.Dropout(0.2)
        self.bn2 = nn.BatchNorm1d(num_features=16)
        self.drop2 = torch.nn.Dropout(0.2)
        
    def forward(self,x):
        x = self.layer1(x)
        x = self.layer2(x)
        x = x.view(-1, 12 * 12)
        x = self.drop1(F.relu(self.bn1(self.fc1(x))))
        x = self.drop2(F.relu(self.bn2(self.fc2(x))))
        x = self.fc3(x)
        return x


main_folder = "MutSigPPI/"
tcga_clinical_dataframe = pickle.load(open(main_folder + "TCGA_clinical_dataframe.pickle","rb"))
classes = {"adrenal gland":0, "bladder":1, "breast":2, "GYN":3, "bile duct":4, "CRC":5, "bone marrow":6, "esophagus":7, "brain":8, "head and neck":9, "Kidney":10, "liver":11, "Lung":12, "pleura":13, "pancreatic":14, "male reproductive system":15, "other":16, "Melanoma":17, "stomach":18, "thyroid":19, "thymus":20}
which_clinicals = ['cancer_class']
tcga_clinical_dataframe = tcga_clinical_dataframe[which_clinicals]
batch_size = 10
dic_cancer_accuracies = {}
for cancer_class in classes:
    print("----------" + cancer_class + "----------")
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
        model = gGapConv()
        criterion = nn.CrossEntropyLoss()
        optimizer = torch.optim.SGD(model.parameters(), lr=lr)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, 1.0, gamma=0.95)
        for shuffled_epoch in range(20):
            print(" Shuffled Epoch=" + str(shuffled_epoch))
            training = specific_cancer_patients.iloc[indices_specific[0]]
            training_other = other_cancer_patients.iloc[indices_other[0]]
            training = shuffle(training.append(training_other))
            Y = training[['cancer_class']].values
            Y = Variable(torch.LongTensor(Y.flatten()), requires_grad=False)
            training = training.index
            for epoch in range(50):
                for index in range(0, len(training), batch_size):
                    y = Y[index : index + batch_size]
                    if(len(y) <= 1):
                        break
                    batch_X = []
                    for patient in training[index : index + batch_size]:
                        patient_df = pd.read_csv(main_folder + "gGap-PatientsDATA/g-Gap-Features/" + patient + ".tsv", sep='\t')
                        batch_X.append(patient_df.to_numpy()[:,1:])
                    batch_X = np.array(batch_X)
                    batch_X = batch_X.astype(np.float32)
                    batch_X = torch.FloatTensor(batch_X)
                    #Expand the input feature from 430*400 to 430*430 by adding zeros!
                    batch_X = F.pad(input=batch_X, pad=(15, 15, 0, 0), mode='constant', value=0)
                    batch_X = batch_X.view(len(y), 1, 430, 430)
                    batch_X = torch.FloatTensor(batch_X)
                    optimizer.zero_grad()
                    Y_hat = model(batch_X)
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
        ii = 0
        acc_avg = 0
        for index in range(0, len(test), batch_size):
            y = Y_test[index : index + batch_size]
            if(len(y) <= 1):
                break
            test_list = []
            for patient in test[index : index + batch_size]:
                patient_df = pd.read_csv(main_folder + "gGap-PatientsDATA/g-Gap-Features/" + patient + ".tsv", sep='\t')
                test_list.append(patient_df.to_numpy()[:,1:])
            test_list = np.array(test_list)
            test_list = test_list.astype(np.float32)
            test_list = torch.FloatTensor(test_list)
            #Expand the input feature from 430*400 to 430*430 by add zeros!
            test_list = F.pad(input=test_list, pad=(15, 15, 0, 0), mode='constant', value=0)
            test_list = test_list.view(len(y), 1, 430, 430)
            test_batch_Y_hat = model.forward(test_list)
            dummy, preds_test = torch.max (test_batch_Y_hat, dim = 1)
            print(preds_test)
            accuracy_test = (preds_test == y).long().sum().float() /  preds_test.size()[0]
            ii += 1
            acc_avg += accuracy_test
        acc_avg = acc_avg / ii
        print("Accuracy: " + str(acc_avg))
        folds_accuracy.append(acc_avg)
        indices_specific = next(parts_specific, None)
        indices_other = next(parts_other, None)
    dic_cancer_accuracies[cancer_class] = folds_accuracy
pickle.dump(dic_cancer_accuracies, open(main_folder + "RESULTS/gGapConv-Accuracies.pickle","wb"))
