#Feedforward Neural Network model and Obtain the accuracy of cancer type prediction  (Energy-based PPI)
import pickle
import pandas as pd
from sklearn.utils import shuffle
from sklearn.model_selection import KFold
import torch.nn.functional as F
from torch.autograd import Variable
import torch
import torch.nn as nn
import numpy as np
import seaborn as sns
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
import os
import seaborn as sns
sns.set(style="whitegrid")

#Binding Affinity-based Feedforward Neural Network Architecture
class BindAffFFNN(nn.Sequential):
    def __init__(self, x_dim):
        super(BindAffFFNN, self).__init__()
        self.linear1 = nn.Linear(x_dim, 512)
        self.bn1 = nn.BatchNorm1d(num_features = 512)
        self.relu1 =  nn.ReLU()
        self.d1 = nn.Dropout(0.05)
        self.linear2 = nn.Linear(512, 256)
        self.bn2 = nn.BatchNorm1d(num_features = 256)
        self.relu2 =  nn.ReLU()
        self.d2 = nn.Dropout(0.05)
        self.linear3 = nn.Linear(256, 128)
        self.bn3 = nn.BatchNorm1d(num_features = 128)
        self.relu3 = nn.ReLU()
        self.d3 = nn.Dropout(0.05)
        self.linear4 = nn.Linear(128, 64)
        self.bn4 = nn.BatchNorm1d(num_features = 64)
        self.relu4 = nn.ReLU()
        self.d4 = nn.Dropout(0.05)
        self.linear5 = nn.Linear(64, 16)
        self.bn5 = nn.BatchNorm1d(num_features = 16)
        self.relu5 = nn.ReLU()
        self.d5 = nn.Dropout(0.05)
        self.linear6 = nn.Linear(16, 2)
         
    def forward(self,x):
        x = self.d1(self.relu1(self.bn1(self.linear1(x))))
        x = self.d2(self.relu2(self.bn2(self.linear2(x))))
        x = self.d3(self.relu3(self.bn3(self.linear3(x))))
        x = self.d4(self.relu4(self.bn4(self.linear4(x))))
        x = self.d5(self.relu5(self.bn5(self.linear5(x))))
        x = self.linear6(x)
        return x


tcga_clinical_dataframe = pickle.load(open("DATASET/TCGA_clinical_dataframe.pickle","rb"))
binding_affinity_dataframe = pickle.load(open("DATASET/BindingAffinilityDataframe.pickle","rb"))
tcga_clinical_dataframe = tcga_clinical_dataframe[['cancer_type']]
cancer_organ_dataframe = pd.read_excel("DATASET/cancer_type.xlsx")
cancer_organ_dataframe = cancer_organ_dataframe[['cancer_type','organ']]
organs_dataframe = tcga_clinical_dataframe.join(cancer_organ_dataframe.set_index('cancer_type'), on='cancer_type')
bindaff_cancerType_df = binding_affinity_dataframe.join(organs_dataframe)
x_dim = 357

classes = {"adrenal gland":0, "bladder":1, "breast":2, "GYN":3, "bile duct":4, "CRC":5, "bone marrow":6, "esophagus":7, "brain":8, "head and neck":9, "Kidney":10, "liver":11, "Lung":12, "pleura":13, "pancreatic":14, "male reproductive system":15, "other":16, "Melanoma":17, "stomach":18, "thyroid":19, "thymus":20}
which_clinicals = ['cancer_class']
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
        #Define the model
        model = BindAffFFNN(x_dim)
        criterion = torch.nn.CrossEntropyLoss()
        optimizer = torch.optim.SGD(model.parameters(), lr=0.001)
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
                    batch = bindaff_cancerType_df.loc[training[index : index + batch_size]]
                    batch = batch.to_numpy()
                    X = batch[:,0:357].astype(float)
                    X = torch.FloatTensor(X)
                    optimizer.zero_grad()
                    Y_hat = model(X)
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
        avg_acc = 0
        ii = 0
        isFirstTime = True
        output_predicted = ""
        for index in range(0, len(test), batch_size):
            y = Y_test[index : index + batch_size]
            test_list = bindaff_cancerType_df.loc[test[index : index + batch_size]]
            test_list = test_list.to_numpy()
            test_list = test_list[:,0:357].astype(float)
            test_list = torch.FloatTensor(test_list)
            if(len(test_list) <= 1):
                break
            test_batch_Y_hat = model.forward(test_list)
            if(isFirstTime):
                output_predicted = test_batch_Y_hat
                isFirstTime = False
            else:
                output_predicted = torch.cat((output_predicted, test_batch_Y_hat), 0)
            dummy, preds_test = torch.max (test_batch_Y_hat, dim = 1)
            accuracy_test = (preds_test == y).long().sum().float() /  preds_test.size()[0]
            avg_acc += accuracy_test
            ii += 1
        avg_acc = avg_acc / ii
        Y_prediction = torch.softmax(output_predicted, dim=1)
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
        print("Fold " + str(fold) + " Accuracy: " + str(avg_acc))
        print("Fold " + str(fold) + " ROC AUC: " + str(roc_auc[1]))
        print("Fold " + str(fold) + " PR AUC: " + str(PR_auc[1]))
        fold += 1
        folds_accuracy.append(avg_acc)
        folds_roc_auc.append(roc_auc[1])
        folds_PR_auc.append(PR_auc[1])
        indices_specific = next(parts_specific, None)
        indices_other = next(parts_other, None)
    if not os.path.exists('RESULTS/BindAffFFNNResults'):
        os.makedirs('RESULTS/BindAffFFNNResults')        
    pickle.dump(folds_accuracy, open("RESULTS/BindAffFFNNResults/" + cancer_class + "_Accuracy.pickle","wb"))
    pickle.dump(folds_roc_auc, open("RESULTS/BindAffFFNNResults/" + cancer_class + "_ROC_AUC.pickle","wb"))
    pickle.dump(folds_PR_auc, open("RESULTS/BindAffFFNNResults/" + cancer_class + "_PR_AUC.pickle","wb"))



#Predict Metastasis (Stage IV) or not (Stages I, II, and III)
#tcga_clinical_dataframe[tcga_clinical_dataframe['stage'] == 'Stage IVA']    
x_dim = 357 
tcga_clinical_dataframe = pickle.load(open("DATASET/TCGA_clinical_dataframe.pickle","rb"))
binding_affinity_dataframe = pickle.load(open("DATASET/BindingAffinilityDataframe.pickle","rb"))
tcga_clinical_dataframe = tcga_clinical_dataframe[['stage']]
bindaff_cancerType_df = binding_affinity_dataframe.join(tcga_clinical_dataframe)
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
        model = BindAffFFNN(x_dim)
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
                    batch = bindaff_cancerType_df.loc[training[index : index + batch_size]]
                    batch = batch.to_numpy()
                    X = batch[:,0:357].astype(float)
                    X = torch.FloatTensor(X)
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
        avg_acc = 0
        ii = 0
        for index in range(0, len(training), batch_size):
            y = Y_test[index : index + batch_size]
            test_list = bindaff_cancerType_df.loc[test[index : index + batch_size]]
            test_list = test_list.to_numpy()
            test_list = test_list[:,0:357].astype(float)
            test_list = torch.FloatTensor(test_list)
            if(len(test_list) <= 1):
                break
            test_batch_Y_hat = model.forward(test_list)
            dummy, preds_test = torch.max (test_batch_Y_hat, dim = 1)
            accuracy_test = (preds_test == y).long().sum().float() /  preds_test.size()[0]
            avg_acc += accuracy_test
            ii += 1
        avg_acc = avg_acc / ii
        print("Fold: " + str(fold_number) + " ACC: " + str(avg_acc))
        fold_number += 1
        folds_accuracy.append(avg_acc)
        indices_metastasis = next(parts_metastasis, None)
        indices_other = next(parts_other, None)
    if not os.path.exists('RESULTS/BindAffFFNN-StagePrediction'):
        os.makedirs('RESULTS/BindAffFFNN-StagePrediction')          
    pickle.dump(folds_accuracy, open("RESULTS/BindAffFFNN-StagePrediction/Part" + str(i) + "_folds_accuracy.pickle","wb"))     
