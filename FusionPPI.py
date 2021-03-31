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

#Fusion of Window-based binary vectors, Gene Expression, Proteins' Pathogenicity, and Binding Affinity of Protein-Protein Interactions
class WinBinVec(nn.Module):
    def __init__(self):
        super(WinBinVec, self).__init__()
        self.fc1 = torch.nn.Linear(430, 256)
        self.fc2 = torch.nn.Linear(256, 430)
        self.conv = torch.nn.Conv1d(430, 430, kernel_size=5, stride=5, padding=0)
        self.fc3 = torch.nn.Linear(2150, 1024)
        self.fc4 = torch.nn.Linear(1024, 512)
        self.fc5 = torch.nn.Linear(512, 256)
        self.fc6 = torch.nn.Linear(256, 64)
        self.fc7 = torch.nn.Linear(64, 32)
        self.bn3 = nn.BatchNorm1d(num_features=1024)
        self.drop3 = torch.nn.Dropout(0.5)
        self.bn4 = nn.BatchNorm1d(num_features=512)
        self.drop4 = torch.nn.Dropout(0.5)
        self.bn5 = nn.BatchNorm1d(num_features=256)
        self.drop5 = torch.nn.Dropout(0.5)

    def forward(self, X):
        X = self.conv(X)
        mX = torch.mean(X, dim=2)
        mX = F.relu(self.fc1(mX))
        mX = torch.sigmoid(self.fc2(mX))
        mX = mX.view(X.size(0), -1, 1).expand_as(X)
        X = X * mX
        X = X.view(X.size(0),-1)
        X = self.drop3(self.bn3(self.fc3(X)))
        X = self.drop4(self.bn4(self.fc4(X)))
        X = self.drop5(self.bn5(self.fc5(X)))
        X = self.fc6(X)
        return self.fc7(X)

class FullyConnected(nn.Module):
    def __init__(self, x_dim):
        super(FullyConnected, self).__init__()
        self.fc1 = nn.Linear(x_dim, 256)
        self.bn1 = nn.BatchNorm1d(num_features = 256)
        self.relu1 =  nn.ReLU()
        self.d1 = nn.Dropout(0.5)
        self.fc2 = nn.Linear(256, 128)
        self.bn2 = nn.BatchNorm1d(num_features = 128)
        self.relu2 =  nn.ReLU()
        self.d2 = nn.Dropout(0.5)
        self.fc3 = nn.Linear(128, 64)
        self.bn3 = nn.BatchNorm1d(num_features = 64)
        self.relu3 = nn.ReLU()
        self.d3 = nn.Dropout(0.5)
        self.fc4 = nn.Linear(64, 32)
        
    def forward(self,x):
        x = self.d1(self.relu1(self.bn1(self.fc1(x))))
        x = self.d2(self.relu2(self.bn2(self.fc2(x))))
        x = self.d3(self.relu3(self.bn3(self.fc3(x))))
        x = self.fc4(x)
        return x

class FusionPPI(nn.Module):
    def __init__(self):
        super(FusionPPI, self).__init__() 
        self.winbinvec = WinBinVec()
        self.pathogenicity = FullyConnected(195)
        self.expression = FullyConnected(572)
        self.bindingaffinity = FullyConnected(357)
        self.fc1 = nn.Linear(128, 32)
        self.bn1 = nn.BatchNorm1d(num_features = 32)
        self.relu1 =  nn.ReLU()
        self.d1 = nn.Dropout(0.5)
        self.fc2 = nn.Linear(32, 8)
        self.bn2 = nn.BatchNorm1d(num_features = 8)
        self.relu2 =  nn.ReLU()
        self.d2 = nn.Dropout(0.5)
        self.fc3 = nn.Linear(8, 2)
        
    def forward(self, win_binary_vectors, smfm_pathogenicity, expression_values, binding_affinity):
        wb = self.winbinvec(win_binary_vectors)
        smfm = self.pathogenicity(smfm_pathogenicity)
        exp = self.expression(expression_values)
        ba = self.bindingaffinity(binding_affinity)
        fusion = torch.stack((wb,smfm, exp, ba), dim=1)
        fusion = fusion.view(wb.size(0), 128)
        output = self.d1(self.relu1(self.bn1(self.fc1(fusion))))
        output = self.d2(self.relu2(self.bn2(self.fc2(output))))
        output = self.fc3(output)
        return output
    
    
#Each Cancer Organ Prediction
tcga_clinical_dataframe = pickle.load(open("DATASET/TCGA_clinical_dataframe.pickle","rb"))
tcga_SMFM_dataframe = pickle.load(open("DATASET/OncomineSMFM.pickle","rb"))
tcga_clinical_dataframe = tcga_clinical_dataframe[['cancer_type']]
cancer_organ_dataframe = pd.read_excel("DATASET/cancer_type.xlsx")
cancer_organ_dataframe = cancer_organ_dataframe[['cancer_type','organ']]
organs_dataframe = tcga_clinical_dataframe.join(cancer_organ_dataframe.set_index('cancer_type'), on='cancer_type')
SMFM_cancerType_df = tcga_SMFM_dataframe.join(organs_dataframe)
tcga_clinical_dataframe = pickle.load(open("DATASET/TCGA_clinical_dataframe.pickle","rb"))
tcga_clinical_dataframe = tcga_clinical_dataframe[['cancer_type']]
cancer_organ_dataframe = pd.read_excel("DATASET/cancer_type.xlsx")
cancer_organ_dataframe = cancer_organ_dataframe[['cancer_type','organ']]
binding_affinity_dataframe = pickle.load(open("DATASET/BindingAffinilityDataframe.pickle","rb"))
organs_dataframe = tcga_clinical_dataframe.join(cancer_organ_dataframe.set_index('cancer_type'), on='cancer_type')
bindaff_cancerType_df = binding_affinity_dataframe.join(organs_dataframe)
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
    specific_cancer_patients = bindaff_cancerType_df[bindaff_cancerType_df["organ"] == cancer_class]
    specific_cancer_patients = specific_cancer_patients.replace({'organ': replace_statement})
    other_cancer_patients = bindaff_cancerType_df[bindaff_cancerType_df["organ"] != cancer_class]
    other_cancer_patients = shuffle(other_cancer_patients).sample(n = len(specific_cancer_patients))
    other_cancer_patients = other_cancer_patients.replace({'organ': replace_statement})
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
        model = FusionPPI()
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
            Y = training[['organ']].values
            Y = Variable(torch.LongTensor(Y.flatten()), requires_grad=False)
            training = training.index
            for epoch in range(50):
                for index in range(0, len(training), batch_size):
                    y = Y[index : index + batch_size]
                    batch_WBV = []
                    batch_Expr = []
                    for patient in training[index : index + batch_size]:
                        wbv_data = pickle.load(open("DATASET/WinBinVecInput/" + patient + "_ppi.pickle", "rb"))
                        expr_data = pickle.load(open("DATASET/ExpressionInputs/" + patient + "_expressions.pickle", "rb"))
                        batch_WBV.append(wbv_data)
                        batch_Expr.append(expr_data)
                    WBV = torch.FloatTensor(batch_WBV)
                    batch_SMFM = SMFM_cancerType_df.loc[training[index : index + batch_size]]
                    batch_SMFM = batch_SMFM.to_numpy()
                    SMFM = batch_SMFM[:,0:195].astype(float)
                    SMFM = torch.FloatTensor(SMFM)
                    EXPR = np.asarray(batch_Expr)
                    EXPR = EXPR.astype(np.float32)
                    EXPR = torch.FloatTensor(EXPR)
                    batch_BA = bindaff_cancerType_df.loc[training[index : index + batch_size]]
                    batch_BA = batch_BA.to_numpy()
                    BA = batch_BA[:,0:357].astype(float)
                    BA = torch.FloatTensor(BA)
                    optimizer.zero_grad()
                    Y_hat = model(WBV, SMFM, EXPR, BA)
                    loss = criterion(Y_hat, y)
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 0.5)
                    optimizer.step()
        test = specific_cancer_patients.iloc[indices_specific[1]]
        test_other = other_cancer_patients.iloc[indices_other[1]]
        test = shuffle(test.append(test_other))
        Y_test = test[['organ']].values
        Y_test = Variable(torch.LongTensor(Y_test.flatten()), requires_grad=False)
        test = test.index
        test_WBV = []
        test_Expr = []
        for patient in test:
            wbv_data = pickle.load(open("DATASET/WinBinVecInput/" + patient + "_ppi.pickle", "rb"))
            expr_data = pickle.load(open("DATASET/ExpressionInputs/" + patient + "_expressions.pickle", "rb"))
            test_WBV.append(wbv_data)
            test_Expr.append(expr_data)
        WBV = torch.FloatTensor(test_WBV)
        test_SMFM = SMFM_cancerType_df.loc[test]
        test_SMFM = test_SMFM.to_numpy()
        SMFM = test_SMFM[:,0:195].astype(float)
        SMFM = torch.FloatTensor(SMFM)
        EXPR = np.asarray(test_Expr)
        EXPR = EXPR.astype(np.float32)
        EXPR = torch.FloatTensor(EXPR)
        test_BA = bindaff_cancerType_df.loc[test]
        test_BA = test_BA.to_numpy()
        BA = test_BA[:,0:357].astype(float)
        BA = torch.FloatTensor(BA)
        test_batch_Y_hat = model.forward(WBV, SMFM, EXPR, BA)
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
    if not os.path.exists('RESULTS/FusionPPIResults'):
        os.makedirs('RESULTS/FusionPPIResults')          
    pickle.dump(folds_accuracy, open("RESULTS/FusionPPIResults/" + cancer_class + "_Accuracy.pickle","wb"))
    pickle.dump(folds_roc_auc, open("RESULTS/FusionPPIResults/" + cancer_class + "_ROC_AUC.pickle","wb"))
    pickle.dump(folds_PR_auc, open("RESULTS/FusionPPIResults/" + cancer_class + "_PR_AUC.pickle","wb"))


#Predict Metastasis (Stage IV) or not (Stages I, II, and III)
#tcga_clinical_dataframe[tcga_clinical_dataframe['stage'] == 'Stage IVA']
tcga_clinical_dataframe = pickle.load(open("DATASET/TCGA_clinical_dataframe.pickle","rb"))
tcga_SMFM_dataframe = pickle.load(open("DATASET/OncomineSMFM.pickle","rb"))
tcga_clinical_dataframe = tcga_clinical_dataframe[['cancer_type']]
cancer_organ_dataframe = pd.read_excel("DATASET/cancer_type.xlsx")
cancer_organ_dataframe = cancer_organ_dataframe[['cancer_type','organ']]
organs_dataframe = tcga_clinical_dataframe.join(cancer_organ_dataframe.set_index('cancer_type'), on='cancer_type')
SMFM_cancerType_df = tcga_SMFM_dataframe.join(organs_dataframe)
tcga_clinical_dataframe = pickle.load(open("DATASET/TCGA_clinical_dataframe.pickle","rb"))
binding_affinity_dataframe = pickle.load(open("DATASET/BindingAffinilityDataframe.pickle","rb"))
tcga_clinical_dataframe = tcga_clinical_dataframe[['stage']]
bindaff_cancerType_df = binding_affinity_dataframe.join(tcga_clinical_dataframe)
tcga_clinical_dataframe = pickle.load(open("DATASET/TCGA_clinical_dataframe.pickle","rb"))    
#tcga_clinical_dataframe = pickle.load(open("DATASET/TCGA_clinical_dataframe.pickle","rb"))
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
metastasis_patients = bindaff_cancerType_df[bindaff_cancerType_df["stage"].isin(metastasis_list)]
metastasis_patients = metastasis_patients.replace({'stage': replace_statement})
other_patients = bindaff_cancerType_df[bindaff_cancerType_df["stage"].isin(other_list)]
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
        model = FusionPPI()
        criterion = torch.nn.CrossEntropyLoss()
        optimizer = torch.optim.SGD(model.parameters(), lr=0.001)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, 1.0, gamma=0.95)
        batch_size = 20
        print("Shuffled Epoch (20): ", end="")
        for shuffled_epoch in range(20):
            if(shuffled_epoch == 19):
                print((shuffled_epoch+1))
            else:
                print((shuffled_epoch+1), end=", ")
            training = metastasis_patients.iloc[indices_metastasis[0]]
            training_other = selected_other_patients.iloc[indices_other[0]]
            training = shuffle(training.append(training_other))
            Y = training[['stage']].values
            Y = Variable(torch.LongTensor(Y.flatten()), requires_grad=False)
            training = training.index
            for epoch in range(50):
                for index in range(0, len(training), batch_size):
                    y = Y[index : index + batch_size]
                    batch_WBV = []
                    batch_Expr = []
                    for patient in training[index : index + batch_size]:
                        wbv_data = pickle.load(open("DATASET/WinBinVecInput/" + patient + "_ppi.pickle", "rb"))
                        expr_data = pickle.load(open("DATASET/ExpressionInputs/" + patient + "_expressions.pickle", "rb"))
                        batch_WBV.append(wbv_data)
                        batch_Expr.append(expr_data)
                    WBV = torch.FloatTensor(batch_WBV)
                    batch_SMFM = SMFM_cancerType_df.loc[training[index : index + batch_size]]
                    batch_SMFM = batch_SMFM.to_numpy()
                    SMFM = batch_SMFM[:,0:195].astype(float)
                    SMFM = torch.FloatTensor(SMFM)
                    EXPR = np.asarray(batch_Expr)
                    EXPR = EXPR.astype(np.float32)
                    EXPR = torch.FloatTensor(EXPR)
                    batch_BA = bindaff_cancerType_df.loc[training[index : index + batch_size]]
                    batch_BA = batch_BA.to_numpy()
                    BA = batch_BA[:,0:357].astype(float)
                    BA = torch.FloatTensor(BA)
                    optimizer.zero_grad()
                    Y_hat = model(WBV, SMFM, EXPR, BA)
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
        test_WBV = []
        test_Expr = []
        for patient in test:
            wbv_data = pickle.load(open("DATASET/WinBinVecInput/" + patient + "_ppi.pickle", "rb"))
            expr_data = pickle.load(open("DATASET/ExpressionInputs/" + patient + "_expressions.pickle", "rb"))
            test_WBV.append(wbv_data)
            test_Expr.append(expr_data)
        WBV = torch.FloatTensor(test_WBV)
        test_SMFM = SMFM_cancerType_df.loc[test]
        test_SMFM = test_SMFM.to_numpy()
        SMFM = test_SMFM[:,0:195].astype(float)
        SMFM = torch.FloatTensor(SMFM)
        EXPR = np.asarray(test_Expr)
        EXPR = EXPR.astype(np.float32)
        EXPR = torch.FloatTensor(EXPR)
        test_BA = bindaff_cancerType_df.loc[test]
        test_BA = test_BA.to_numpy()
        BA = test_BA[:,0:357].astype(float)
        BA = torch.FloatTensor(BA)
        test_batch_Y_hat = model.forward(WBV, SMFM, EXPR, BA)
        dummy, preds_test = torch.max (test_batch_Y_hat, dim = 1)
        accuracy_test = (preds_test == Y_test).long().sum().float() /  preds_test.size()[0]
        print("Fold: " + str(fold_number) + " ACC: " + str(accuracy_test))
        fold_number += 1
        folds_accuracy.append(accuracy_test)
        indices_metastasis = next(parts_metastasis, None)
        indices_other = next(parts_other, None)
    if not os.path.exists('RESULTS/FusionPPI-StagePrediction'):
        os.makedirs('RESULTS/FusionPPI-StagePrediction')        
    pickle.dump(folds_accuracy, open("RESULTS/FusionPPI-StagePrediction/Part" + str(i) + "_folds_accuracy.pickle","wb"))     
