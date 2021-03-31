import pickle
import pandas as pd
from sklearn.utils import shuffle
from sklearn.model_selection import KFold
import numpy as np
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


class DeepCAT(nn.Module):
    def __init__(self):
        super(DeepCAT, self).__init__()
        self.conv1 = nn.Sequential(
            nn.Conv2d(1, 8, kernel_size=(15,2)),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=(1,2), stride=1))
        self.conv2 = nn.Sequential(
            nn.Conv2d(8, 16, kernel_size=(1,2)),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=(1,2), stride=1))
        self.fc1 = nn.Linear(16 * 21, 64)
        self.fc2 = nn.Linear(64, 10)
        self.fc3 = nn.Linear(10, 2)
        #Batch Normalization and Dropout Layers
        self.drop1 = torch.nn.Dropout(0.4)
        self.drop2 = torch.nn.Dropout(0.4)

    def forward(self, X):
        #X is [1, 1, 15, 25]
        #Pass the input through the first convolutional module
        #[1, 8, 1, 23]
        X = self.conv1(X)
        #Pass through the second convolutional module
        #[1, 16, 1, 21]
        X = self.conv2(X)
        #[1, 336]
        X = X.view(-1, 16 * 21)
        X = self.drop1(F.relu(self.fc1(X)))
        X = self.drop2(F.relu(self.fc2(X)))
        X = self.fc3(X)
        return X
    
class PhyChemPPI(nn.Module):
    def __init__(self):
        super(PhyChemPPI, self).__init__()
        self.number_of_sequences = 430
        self.deepCAT = torch.nn.ModuleList([DeepCAT() for i in range(self.number_of_sequences)])
        self.fc1 = nn.Linear(860, 256)
        self.fc2 = nn.Linear(256, 8)
        self.fc3 = nn.Linear(8, 2)
        self.bn1 = nn.BatchNorm1d(num_features=256)
        self.drop1 = torch.nn.Dropout(0.4)
        self.bn2 = nn.BatchNorm1d(num_features=8)
        self.drop2 = torch.nn.Dropout(0.4)
        
    def forward(self, X):
        outs = []
        #Consider each Protein-Protein Interaction (430 Protein-Protein Interactions)
        for i in range(self.number_of_sequences):
            #Number of Amino Acids in current Protein-Protein Interaction
            #seq_number_of_aminoacids = len(X_array[0][i])
            seq_number_of_aminoacids = 25   
            #Obtain DeepCAT features
            deepCAT_input = X[:,i]
            deepCAT_input = deepCAT_input.view(deepCAT_input.size(0), 1, 15, 25)
            out = self.deepCAT[i](deepCAT_input)
            outs.append(out) 
        concat_dc = outs[0]
        for i in range(1,self.number_of_sequences):
            concat_dc = torch.cat((concat_dc, outs[i]), 1)
        output = self.drop1(F.relu(self.bn1(self.fc1(concat_dc))))
        output = self.drop2(F.relu(self.bn2(self.fc2(output))))
        output = self.fc3(output)
        return output
    

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
        #Define the model
        model = PhyChemPPI()
        # Mean Squared Error
        criterion = torch.nn.CrossEntropyLoss()
        # Stochastic Gradient Descent
        optimizer = torch.optim.SGD(model.parameters(), lr=0.001)
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
                        p_data = pickle.load(open("DATASET/PhysicochemicalPropsInputs/" + patient + "_ppi.pickle", "rb"))
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
        avg_acc = 0
        ii = 0
        isFirstTime = True
        output_predicted = ""        
        for index in range(0, len(training), batch_size):
            y = Y_test[index : index + batch_size]
            test_list = []
            for patient in test[index : index + batch_size]:
                p_data = pickle.load(open("DATASET/PhysicochemicalPropsInputs/" + patient + "_ppi.pickle", "rb"))
                test_list.append(p_data)
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
    if not os.path.exists('RESULTS/PHYCHEMResults'):
        os.makedirs('RESULTS/PHYCHEMResults')         
    pickle.dump(folds_accuracy, open("RESULTS/PHYCHEMResults/" + cancer_class + "_Accuracy.pickle","wb"))
    pickle.dump(folds_roc_auc, open("RESULTS/PHYCHEMResults/" + cancer_class + "_ROC_AUC.pickle","wb"))
    pickle.dump(folds_PR_auc, open("RESULTS/PHYCHEMResults/" + cancer_class + "_PR_AUC.pickle","wb"))



#Predict Metastasis (Stage IV) or not (Stages I, II, and III)
#tcga_clinical_dataframe[tcga_clinical_dataframe['stage'] == 'Stage IVA']
tcga_clinical_dataframe = pickle.load(open("DATASET/TCGA_clinical_dataframe.pickle","rb"))
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
        model = PhyChemPPI()
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
                    batch_X = []
                    for patient in training[index : index + batch_size]:
                        p_data = pickle.load(open("DATASET/PhysicochemicalPropsInputs/" + patient + "_ppi.pickle", "rb"))
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
        avg_acc = 0
        ii = 0
        for index in range(0, len(training), batch_size):
            y = Y_test[index : index + batch_size]
            test_list = []
            for patient in test[index : index + batch_size]:
                p_data = pickle.load(open("DATASET/PhysicochemicalPropsInputs/" + patient + "_ppi.pickle", "rb"))
                test_list.append(p_data)
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
    if not os.path.exists('RESULTS/PHYCHEM-StagePrediction'):
        os.makedirs('RESULTS/PHYCHEM-StagePrediction')        
    pickle.dump(folds_accuracy, open("RESULTS/PHYCHEM-StagePrediction/Part" + str(i) + "_folds_accuracy.pickle","wb")) 
