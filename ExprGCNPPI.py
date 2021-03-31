from __future__ import division
from __future__ import print_function
import time
import argparse
import numpy as np
import torch
import torch.nn.functional as F
import torch.optim as optim
import math
import torch.nn as nn
from torch.nn.parameter import Parameter
from torch.nn.modules.module import Module
import numpy as np
import seaborn as sns
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
import os
import pickle
from sklearn.model_selection import KFold
from torch.autograd import Variable
from sklearn.utils import shuffle

#"GraphConvolution" and "GCN" classes are obtained from the following GitHub repository
#https://github.com/tkipf/pygcn
class GraphConvolution(Module):
    def __init__(self, in_features, out_features, bias=True):
        super(GraphConvolution, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.weight = Parameter(torch.FloatTensor(in_features, out_features))
        if bias:
            self.bias = Parameter(torch.FloatTensor(out_features))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()

    def reset_parameters(self):
        stdv = 1. / math.sqrt(self.weight.size(1))
        self.weight.data.uniform_(-stdv, stdv)
        if self.bias is not None:
            self.bias.data.uniform_(-stdv, stdv)

    def forward(self, input, adj):
        support = torch.mm(input, self.weight)
        output = torch.spmm(adj, support)
        if self.bias is not None:
            return output + self.bias
        else:
            return output

    def __repr__(self):
        return self.__class__.__name__ + ' (' \
               + str(self.in_features) + ' -> ' \
               + str(self.out_features) + ')'

class GCN(nn.Module):
    def __init__(self, nfeat, nhid, nclass, dropout):
        super(GCN, self).__init__()
        self.gc1 = GraphConvolution(nfeat, nhid)
        self.gc2 = GraphConvolution(nhid, nclass)
        self.dropout = dropout

    def forward(self, x, adj):
        x = F.relu(self.gc1(x, adj))
        x = F.dropout(x, self.dropout, training=self.training)
        x = self.gc2(x, adj)
        return F.log_softmax(x, dim=1)


#ExprGCNPPI implemented by Sina Abdollahi for WinBinVec paper
class ExprGCNPPI(nn.Module):
    def __init__(self):
        super(ExprGCNPPI, self).__init__()
        self.gcn = GCN(nfeat=1, nhid=8, nclass=1, dropout=0.5)
        #626: The number of partner proteins involve in the PPIs
        self.fc1 = nn.Linear(572, 256)
        self.fc2 = nn.Linear(256, 8)
        self.fc3 = nn.Linear(8, 2)
        self.bn1 = nn.BatchNorm1d(num_features=256)
        self.drop1 = torch.nn.Dropout(0.4)
        self.bn2 = nn.BatchNorm1d(num_features=8)
        self.drop2 = torch.nn.Dropout(0.4)
        
    def forward(self, adj, expr, batch_size):
        outs = []
        for i in range(batch_size):
            outs.append(self.gcn(expr[i,:], adj[i,:]).view(1,adj.size(1)))
        concat_gcn = outs[0]
        for i in range(1,batch_size):
            concat_gcn = torch.cat([concat_gcn, outs[i]], dim=0)
        output = self.drop1(F.relu(self.bn1(self.fc1(concat_gcn))))
        output = self.drop2(F.relu(self.bn2(self.fc2(output))))
        output = self.fc3(output)
        return output
    
    
ppi_adj_matrix = pickle.load(open("DATASET/Adj.pickle", "rb"))
#adj -> Adjacency matrix (N by N ---> N is the number of nodes)
#expr -> Features (Expression) matrix (N by 1 ---> 1 is for gene expression value for each protein)
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
        model = ExprGCNPPI()
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
                    kk = 0
                    for patient in training[index : index + batch_size]:
                        kk += 1
                        p_data = pickle.load(open("DATASET/ExpressionInputs/" + patient + "_expressions.pickle", "rb"))
                        batch_X.append(p_data)
                    X = np.asarray(batch_X)
                    adj = np.array([ppi_adj_matrix]*batch_size)
                    adj = adj.astype(np.float32)
                    adj = torch.FloatTensor(adj)
                    X = X.astype(np.float32)
                    X = torch.FloatTensor(X)
                    X = X.view(X.size(0), X.size(1), 1)
                    optimizer.zero_grad()
                    Y_hat = model(adj, X, kk)
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
        isFirstTime = True
        output_predicted = ""
        ii = 0
        for index in range(0, len(training), batch_size):
            y = Y_test[index : index + batch_size]
            test_list = []
            kk = 0
            for patient in test[index : index + batch_size]:
                kk += 1
                p_data = pickle.load(open("DATASET/ExpressionInputs/" + patient + "_expressions.pickle", "rb"))
                test_list.append(p_data)
            #test_list = torch.FloatTensor(test_list)
            if(len(test_list) <= 1):
                break
            X_test = np.asarray(test_list)
            adj = np.array([ppi_adj_matrix]*batch_size)
            adj = adj.astype(np.float32)
            adj = torch.FloatTensor(adj)
            X_test = X_test.astype(np.float32)
            X_test = torch.FloatTensor(X_test)
            X_test = X_test.view(X_test.size(0), X_test.size(1), 1)  
            test_batch_Y_hat = model.forward(adj, X_test, kk)
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
    if not os.path.exists('RESULTS/ExprGCNPPIResults'):
        os.makedirs('RESULTS/ExprGCNPPIResults')        
    pickle.dump(folds_accuracy, open("RESULTS/ExprGCNPPIResults/" + cancer_class + "_Accuracy.pickle","wb"))
    pickle.dump(folds_roc_auc, open("RESULTS/ExprGCNPPIResults/" + cancer_class + "_ROC_AUC.pickle","wb"))
    pickle.dump(folds_PR_auc, open("RESULTS/ExprGCNPPIResults/" + cancer_class + "_PR_AUC.pickle","wb"))



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
        model = ExprGCNPPI()
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
                    kk = 0
                    for patient in training[index : index + batch_size]:
                        kk += 1
                        p_data = pickle.load(open("DATASET/ExpressionInputs/" + patient + "_expressions.pickle", "rb"))
                        batch_X.append(p_data)
                    X = np.asarray(batch_X)
                    adj = np.array([ppi_adj_matrix]*batch_size)
                    adj = adj.astype(np.float32)
                    adj = torch.FloatTensor(adj)
                    X = X.astype(np.float32)
                    X = torch.FloatTensor(X)
                    X = X.view(X.size(0), X.size(1), 1)
                    optimizer.zero_grad()
                    Y_hat = model(adj, X, kk)
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
            kk = 0
            for patient in test[index : index + batch_size]:
                kk += 1
                p_data = pickle.load(open("DATASET/ExpressionInputs/" + patient + "_expressions.pickle", "rb"))
                test_list.append(p_data)
            if(len(test_list) <= 1):
                break
            X_test = np.asarray(test_list)
            adj = np.array([ppi_adj_matrix]*batch_size)
            adj = adj.astype(np.float32)
            adj = torch.FloatTensor(adj)
            X_test = X_test.astype(np.float32)
            X_test = torch.FloatTensor(X_test)
            X_test = X_test.view(X_test.size(0), X_test.size(1), 1)  
            test_batch_Y_hat = model.forward(adj, X_test, kk)
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
    if not os.path.exists('RESULTS/ExprGCNPPI-StagePrediction'):
        os.makedirs('RESULTS/ExprGCNPPI-StagePrediction')         
    pickle.dump(folds_accuracy, open("RESULTS/ExprGCNPPI-StagePrediction/Part" + str(i) + "_folds_accuracy.pickle","wb")) 
