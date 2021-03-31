#WinBinVec with Fully-Connected layers instead of Convolutional NN
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
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score

#WinBinVec-FC uses Fully-connected layers instead on 1D Convolutional layers
class WinBinVecFC(nn.Module):
    def __init__(self):
        super(WinBinVecFC, self).__init__()
        self.linears = torch.nn.ModuleList([torch.nn.Linear(25, 4) for i in range(430)])
        self.fc3 = torch.nn.Linear(1720, 1024)
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
        outs = []
        #Consider each Protein-Protein Interaction
        for i in range(430):
            out = self.linears[i](X[:,i])
            outs.append(out)    
        #Concatenate all obtained PPI vectors
        out = outs[0]
        for i in range(1,430):
            out = torch.cat((out, outs[i]), 1)
        out = self.drop3(self.bn3(self.fc3(out)))
        out = self.drop4(self.bn4(self.fc4(out)))
        out = self.drop5(self.bn5(self.fc5(out)))
        out = self.fc6(out)
        out = self.fc7(out)
        return self.fc8(out)

#Each Cancer Organ Prediction
tcga_clinical_dataframe = pickle.load(open("DATASET/TCGA_clinical_dataframe.pickle","rb"))
classes = {"adrenal gland":0, "bladder":1, "breast":2, "GYN":3, "bile duct":4, "CRC":5, "bone marrow":6, "esophagus":7, "brain":8, "head and neck":9, "Kidney":10, "liver":11, "Lung":12, "pleura":13, "pancreatic":14, "male reproductive system":15, "other":16, "Melanoma":17, "stomach":18, "thyroid":19, "thymus":20}
which_clinicals = ['cancer_class']
tcga_clinical_dataframe = tcga_clinical_dataframe[which_clinicals]
for cancer_class in classes:
    print(">>>>>>" + cancer_class)
    folds_accuracy = []
    folds_roc_auc = []
    folds_PR_auc = []
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
    fold = 1
    while(indices_specific):
        lr = 0.5
        #Define the model
        model = WinBinVecFC()
        criterion = nn.CrossEntropyLoss()
        optimizer = torch.optim.SGD(model.parameters(), lr=lr)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, 1.0, gamma=0.95)
        batch_size = 20
        print("Shuffled Epoch (20): ", end="")
        for shuffled_epoch in range(20):
            if(shuffled_epoch == 19):
                print((shuffled_epoch+1))
            else:
                print((shuffled_epoch+1), end=", ")
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
                        p_data = pickle.load(open("DATASET/WinBinVecInput/" + patient + "_ppi.pickle", "rb"))
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
            p_data = pickle.load(open("DATASET/WinBinVecInput/" + patient + "_ppi.pickle", "rb"))
            test_list.append(p_data)
        test_list = torch.FloatTensor(test_list)
        test_batch_Y_hat = model.forward(test_list)
        dummy, preds_test = torch.max (test_batch_Y_hat, dim = 1)
        accuracy_test = (preds_test == Y_test).long().sum().float() /  preds_test.size()[0]
        Y_prediction = torch.softmax(test_batch_Y_hat, dim=1)
        Y_prediction = np.array(Y_prediction.tolist())
        Y_real = np.array([[1,0] if y == 0 else [0,1] for y in Y_test])
        fpr = dict()
        tpr = dict()
        precision = dict()
        recall = dict()
        roc_auc = dict()
        PR_auc = dict()
        for i in range(2):
            fpr[i], tpr[i], _ = roc_curve(Y_real[:, i], Y_prediction[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])
            precision[i], recall[i], _ = precision_recall_curve(Y_real[:, i], Y_prediction[:, i])
            PR_auc[i] = auc(recall[i], precision[i])
        print("Fold " + str(fold) + " Accuracy: " + str(accuracy_test))
        print("Fold " + str(fold) + " ROC AUC: " + str(roc_auc[1]))
        print("Fold " + str(fold) + " PR AUC: " + str(PR_auc[1]))
        fold += 1
        folds_accuracy.append(accuracy_test)
        folds_roc_auc.append(roc_auc[1])
        folds_PR_auc.append(PR_auc[1])
        indices_specific = next(parts_specific, None)
        indices_other = next(parts_other, None)
    if not os.path.exists('RESULTS/WinBinVecFCResults'):
        os.makedirs('RESULTS/WinBinVecFCResults')          
    pickle.dump(folds_accuracy, open("RESULTS/WinBinVecFCResults/" + cancer_class + "_Accuracy.pickle","wb"))
    pickle.dump(folds_roc_auc, open("RESULTS/WinBinVecFCResults/" + cancer_class + "_ROC_AUC.pickle","wb"))
    pickle.dump(folds_PR_auc, open("RESULTS/WinBinVecFCResults/" + cancer_class + "_PR_AUC.pickle","wb"))
