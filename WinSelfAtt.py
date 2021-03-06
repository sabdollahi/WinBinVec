import math
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import TransformerEncoder, TransformerEncoderLayer
from sklearn.model_selection import KFold
from sklearn.utils import shuffle
from torch.autograd import Variable
import pickle
import numpy as np
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
import os

#nhead: Number of Heads in Transformer
#nlayers: Number of Layers in Feedforward Neural Network in Transformer
#embedding_size: Input Embedding Size for each amino acid (In this case is 20 - The number of amino acids in human body)
#sequence_vector_size: Embedding Size obtain for each sequence (for each PPI)

#Do care about the position of each amino acid (Using sine and cosine)
class PositionalEncoding(torch.nn.Module):
    def __init__(self, d_model, dropout=0.1, max_len=5000):
        super(PositionalEncoding, self).__init__()
        self.dropout = torch.nn.Dropout(p=dropout)
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0).transpose(0, 1)
        self.register_buffer('pe', pe)

    def forward(self, x):
        x = x + self.pe[:x.size(0), :]
        return self.dropout(x)
    
#Transformer based on the Positional Encoding (The input passes through the PositionalEncoding)    
class PosTransformerModel(nn.Module):
    def __init__(self, embedding_size, nhead, ffn_dim, nlayers, dropout=0.5):
        super(PosTransformerModel, self).__init__()
        self.model_type = 'Transformer'
        self.pos_encoder = PositionalEncoding(embedding_size, dropout)
        encoder_layers = TransformerEncoderLayer(embedding_size, nhead, ffn_dim, dropout)
        self.transformer_encoder = TransformerEncoder(encoder_layers, nlayers)
        self.embedding_size = embedding_size

    def forward(self, src):
        #Pass the input through Positional Encoding (Because we do care about the position of each amino acid)
        src = self.pos_encoder(src)
        #Then, pass the result through the Transformer encoder
        output = self.transformer_encoder(src)
        return output
    
#In the case that we DO NOT CARE about the position
#After we will extract the vector of each Protein-Protein interaction sequence, 
#We DO NOT CARE about the position of the obtained vectors when we want to make a prediction using Self-Attention
class TransformerModel(nn.Module):
    def __init__(self, embedding_size, nhead, ffn_dim, nlayers, dropout=0.5):
        super(TransformerModel, self).__init__()
        self.model_type = 'Transformer'
        encoder_layers = TransformerEncoderLayer(embedding_size, nhead, ffn_dim, dropout)
        self.transformer_encoder = TransformerEncoder(encoder_layers, nlayers)
        self.embedding_size = embedding_size

    def forward(self, src):
        output = self.transformer_encoder(src)
        return output
    
#Protein-Protein Interaction Self-Attention (The final model)
class WinSelfAtt(nn.Module):
    def __init__(self, number_of_sequences, sequence_vector_size, embedding_size, nhead, ffn_dim, nlayers, dropout):
        super(WinSelfAtt, self).__init__()
        self.non_position_transformer = TransformerModel(embedding_size=sequence_vector_size, nhead=nhead, ffn_dim=ffn_dim, nlayers=nlayers, dropout=dropout)
        # 2 over here means that the classification is a binary classification
        self.patients_linear = torch.nn.Linear(number_of_sequences*sequence_vector_size, 8)
        self.last_linear = torch.nn.Linear(8, 2)
        self.embedding_size = embedding_size
        self.sequence_vector_size = sequence_vector_size
        self.position_transformers = torch.nn.ModuleList([PosTransformerModel(embedding_size=embedding_size, nhead=nhead, ffn_dim=ffn_dim, nlayers=nlayers, dropout=dropout) for i in range(number_of_sequences)])
        self.linears = torch.nn.ModuleList([torch.nn.Linear(25*embedding_size, sequence_vector_size) for i in range(number_of_sequences)])
            
    def forward(self, X_array):
        #Number of Protein-Protein Interactions = 430
        number_of_sequences = 430
        outs = []
        #Consider each Protein-Protein Interaction
        for i in range(number_of_sequences):
            #Number of Amino Acids in current Protein-Protein Interaction
            #seq_number_of_aminoacids = len(X_array[0][i])
            seq_number_of_aminoacids = 25   
            #Use Postion-based Transformer and extract vectors for all amino acids
            out = self.position_transformers[i](X_array[:,i])
            #Concatenate all obtained amino acid vectors
            out = out.view(-1, seq_number_of_aminoacids*self.embedding_size)
            #Use a dense layer with output size of 'sequence_vector_size'
            #To obtain the vector representation of current PPI
            out = self.linears[i](out)
            outs.append(out)    
        #Concatenate all obtained PPI vectors
        second_stage_tensor = outs[0]
        for i in range(1,number_of_sequences):
            second_stage_tensor = torch.cat((second_stage_tensor, outs[i]), 1)
        second_stage_tensor = second_stage_tensor.view(-1, number_of_sequences, self.sequence_vector_size)
        #Pass the Concatenated PPI vectors through a Non-Positional Transformer 
        out = self.non_position_transformer(second_stage_tensor)
        out = out.view(-1, number_of_sequences*self.sequence_vector_size)
        #Then, pass the result through a dense layer
        out = self.patients_linear(out)
        out = self.last_linear(out)
        #No need to use Softmax for Multi-class classification in PyTorch when we use CrossEntropy Loss function
        #out = torch.nn.functional.softmax(out, dim=1)
        return out



#Each Cancer Organ Prediction in an NLP-based Neural Network Architecture
#Number of PPIs
number_of_sequences = 430
#Input Embedding Size
em_size = 20
#Embedding Size obtain for each sequence
sequence_vector_size = 10
#Number of heads (in Multihead Attention)
nhead = 5
#Feedforward Neural Network dimension
ffn_dim = 10
#Number of FFNN layers
nlayers = 4
#Dropout amount
dropout = 0.5
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
        model = WinSelfAtt(number_of_sequences, sequence_vector_size, em_size, nhead, ffn_dim, nlayers, dropout)
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
                        p_data = pickle.load(open("DATASET/WinSelfAttInput/" + patient + "_ppi.pickle", "rb"))
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
        isFirstTime = True
        output_predicted = ""
        ii = 0
        for index in range(0, len(training), batch_size):
            y = Y_test[index : index + batch_size]
            test_list = []
            for patient in test[index : index + batch_size]:
                p_data = pickle.load(open("DATASET/WinSelfAttInput/" + patient + "_ppi.pickle", "rb"))
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
    if not os.path.exists('RESULTS/WinSelfAttResults'):
        os.makedirs('RESULTS/WinSelfAttResults')
    pickle.dump(folds_accuracy, open("RESULTS/WinSelfAttResults/" + cancer_class + "_Accuracy.pickle","wb"))
    pickle.dump(folds_roc_auc, open("RESULTS/WinSelfAttResults/" + cancer_class + "_ROC_AUC.pickle","wb"))
    pickle.dump(folds_PR_auc, open("RESULTS/WinSelfAttResults/" + cancer_class + "_PR_AUC.pickle","wb"))



#Predict Metastasis (Stage IV) or not (Stages I, II, and III)
#tcga_clinical_dataframe[tcga_clinical_dataframe['stage'] == 'Stage IVA']
number_of_sequences = 430
em_size = 20
sequence_vector_size = 10
nhead = 5
ffn_dim = 10
nlayers = 4
dropout = 0.5
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
        model = WinSelfAtt(number_of_sequences, sequence_vector_size, em_size, nhead, ffn_dim, nlayers, dropout)
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
                        p_data = pickle.load(open("DATASET/WinSelfAttInput/" + patient + "_ppi.pickle", "rb"))
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
                p_data = pickle.load(open("DATASET/WinSelfAttInput/" + patient + "_ppi.pickle", "rb"))
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
    if not os.path.exists('RESULTS/WinSelfAtt-StagePrediction'):
        os.makedirs('RESULTS/WinSelfAtt-StagePrediction') 
    pickle.dump(folds_accuracy, open("RESULTS/WinSelfAtt-StagePrediction/Part" + str(i) + "_folds_accuracy.pickle","wb"))   
