#Extract One-hot-Vector
import os
import pickle
import pandas as pd
from sklearn.utils import shuffle
from sklearn.model_selection import KFold
from torch.autograd import Variable
import torch

main_folder = "MutSigPPI/"
aa_oneHotVec = pickle.load(open(main_folder + "AminoAcids_OneHotVectors.pickle", "rb"))
only_mutated = main_folder + "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
patients_properties = main_folder + "PatientsSequencesInEachRefSeqChain/"
classes = ["adrenal gland", "bladder", "breast", "GYN", "bile duct", "CRC", "bone marrow", "esophagus", "brain", "head and neck", "Kidney", "liver", "Lung", "pleura", "pancreatic", "male reproductive system", "other", "Melanoma", "stomach", "thyroid", "thymus"]
cancer_names = ["acc","blca","brca","cesc","chol","coad","dlbc","esca","gbm","hnsc","kich","kirc","kirp","laml","lgg","lihc","luad","lusc","meso","ov","paad","pcpg","prad","read","sarc","skcm","stad","tgct","thca","thym","ucec","ucs","uvm"]
main_folder = "MutSigPPI/"
which_clinicals = ['vital_status', 'cancer_type']
vital_status = ["Alive", "Dead"]
tcga_clinical_dataframe = pickle.load(open(main_folder + "TCGA_clinical_dataframe.pickle","rb"))
tcga_clinical_dataframe = tcga_clinical_dataframe[which_clinicals]
#Make the dataset balance for outputs (In this case: 2018 patients for Dead dataset and 2032 patients for Alive dataset)
dead_df = tcga_clinical_dataframe[tcga_clinical_dataframe['vital_status'] == 'Dead']
alive_dataframe = tcga_clinical_dataframe[tcga_clinical_dataframe['vital_status'] == 'Alive']
alive_df = shuffle(alive_dataframe[alive_dataframe['cancer_type'] == cancer_names[0]]).head(60)
for i in range(1, len(cancer_names)):
    alive_df = alive_df.append(shuffle(alive_dataframe[alive_dataframe['cancer_type'] == cancer_names[i]]).head(69))
#Put 4050 Alive/Dead balanced dataset only!
tcga_clinical_dataframe = alive_df.append(dead_df)

K = 10 #Kfold (number of parts = K)
kf = KFold(n_splits = K, shuffle = True)
kfold_parts = []
status_dataframes = []
for status in vital_status:
    tmp_df = shuffle(tcga_clinical_dataframe[tcga_clinical_dataframe['vital_status'] == status])
    status_dataframes.append(tmp_df)
    tmp_parts = kf.split(tmp_df)
    kfold_parts.append(tmp_parts)
    
#Indices for different cancer types
indices = []
for i in range(len(kfold_parts)):
    indices.append(next(kfold_parts[i], None)) 
    
fold_num = 1
while(indices[0]):
    training_set = status_dataframes[0].iloc[indices[0][0]]
    test_set = status_dataframes[0].iloc[indices[0][1]]
    for i in range(1,len(indices)):
        train = status_dataframes[i].iloc[indices[i][0]]
        training_set = training_set.append(train)
        test =  status_dataframes[i].iloc[indices[i][1]]
        test_set = test_set.append(test)
    #Training set outputs (Alive = 1 and Dead = 0)
    training_set = shuffle(training_set)
    Y = training_set[['vital_status']].replace({'vital_status': {'Dead': 0, 'Alive': 1}}).values
    Y = Variable(torch.LongTensor(Y.flatten()), requires_grad=False)
    training_set = training_set.index
    #Test set outputs
    Y_test = test_set[['vital_status']].replace({'vital_status': {'Dead': 0, 'Alive': 1}}).values
    Y_test = Variable(torch.LongTensor(Y_test.flatten()), requires_grad=False)
    test_set = test_set.index
    
    #STORE test set 
    test_list = []
    dic_test_set = {}
    dic_test_set["X"] = ""
    dic_test_set["Y"] = Y_test
#------------------------------------------TEST One-hot VECTORS-------------------------------------    
    for patient in test_set:
        all_patients_current_ppi_oneHotVectors = []
        with open(main_folder + "PPIs-WithRefSeqs.txt","rb") as f:
            for interaction in f.readlines():
                interaction = str(interaction,'utf-8')
                interaction = interaction.rstrip()
                splittedLine = interaction.split()
                #First Chain Protein name and Chain name
                prot_1 = splittedLine[0]
                chain_1 = splittedLine[1]
                #Second Chain Protein name and Chain name
                prot_2 = splittedLine[3]
                chain_2 = splittedLine[4]
                current_ppi_oneHotVectors = []
                address1 = ""
                if(os.path.exists(patients_properties + patient + "/" + prot_1 + "/" + chain_1)):
                    address1 = patients_properties + patient + "/" + prot_1 + "/" + chain_1 + "/mutated_sequence.seq"
                else:
                    address1 = only_mutated + prot_1 + "/" + chain_1 + "/mutatedTrimedRefSeq.seq"
                with open(address1, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        l = l.rstrip()
                        for aa in l:
                            current_ppi_oneHotVectors.append(aa_oneHotVec[aa])
                address2 = ""
                if(os.path.exists(patients_properties + patient + "/" + prot_2 + "/" + chain_2)):
                    address2 = patients_properties + patient + "/" + prot_2 + "/" + chain_2 + "/mutated_sequence.seq" 
                else:
                    address2 = only_mutated + prot_2 + "/" + chain_2 + "/mutatedTrimedRefSeq.seq"
                with open(address2, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        l = l.rstrip()
                        for aa in l:
                            current_ppi_oneHotVectors.append(aa_oneHotVec[aa])
                current_len = len(current_ppi_oneHotVectors)
                if(current_len <= 149):
                    for k in range(current_len, 149):
                        current_ppi_oneHotVectors.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                    all_patients_current_ppi_oneHotVectors.append(current_ppi_oneHotVectors)
            test_list.append(all_patients_current_ppi_oneHotVectors)
    dic_test_set["X"] = test_list
    pickle.dump(dic_test_set, open("MutSigPPI/DATASET/test_set_fold" + str(fold_num) + ".pickle", "wb"))
    print("TEST: fold=" + str(fold_num) + " has already finished!")
#------------------------------------------TRAINING One-hot VECTORS------------------------------------- 
    batch_size = 5
    batch_number = 0
    for index in range(0, len(training_set), batch_size):
        batch_number += 1
        #STORE input 'X' ready to feed our neural network model
        input_list = []
        #Dictionary to save training/test sets
        dic_batch_set = {}
        dic_batch_set['X'] = ""
        dic_batch_set['Y'] = Y[index : index + batch_size]
        batch_set = training_set[index : index + batch_size]
        for patient in batch_set:
            all_patients_current_ppi_oneHotVectors = []
            with open(main_folder + "PPIs-WithRefSeqs.txt","rb") as f:
                for interaction in f.readlines():
                    interaction = str(interaction,'utf-8')
                    interaction = interaction.rstrip()
                    splittedLine = interaction.split()
                    #First Chain Protein name and Chain name
                    prot_1 = splittedLine[0]
                    chain_1 = splittedLine[1]
                    #Second Chain Protein name and Chain name
                    prot_2 = splittedLine[3]
                    chain_2 = splittedLine[4]
                    current_ppi_oneHotVectors = []
                    address1 = ""
                    if(os.path.exists(patients_properties + patient + "/" + prot_1 + "/" + chain_1)):
                        address1 = patients_properties + patient + "/" + prot_1 + "/" + chain_1 + "/mutated_sequence.seq"
                    else:
                        address1 = only_mutated + prot_1 + "/" + chain_1 + "/mutatedTrimedRefSeq.seq"
                    with open(address1, "rb") as seqf:
                        for l in seqf.readlines():
                            l = str(l,'utf-8')
                            l = l.rstrip()
                            for aa in l:
                                current_ppi_oneHotVectors.append(aa_oneHotVec[aa])
                    address2 = ""
                    if(os.path.exists(patients_properties + patient + "/" + prot_2 + "/" + chain_2)):
                        address2 = patients_properties + patient + "/" + prot_2 + "/" + chain_2 + "/mutated_sequence.seq" 
                    else:
                        address2 = only_mutated + prot_2 + "/" + chain_2 + "/mutatedTrimedRefSeq.seq"
                    with open(address2, "rb") as seqf:
                        for l in seqf.readlines():
                            l = str(l,'utf-8')
                            l = l.rstrip()
                            for aa in l:
                                current_ppi_oneHotVectors.append(aa_oneHotVec[aa])
                    current_len = len(current_ppi_oneHotVectors)
                    if(current_len <= 149):
                        for k in range(current_len, 149):
                            current_ppi_oneHotVectors.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                        all_patients_current_ppi_oneHotVectors.append(current_ppi_oneHotVectors)
                input_list.append(all_patients_current_ppi_oneHotVectors)
        dic_batch_set['X'] = input_list
        pickle.dump(dic_batch_set, open("MutSigPPI/DATASET/training_set_fold" + str(fold_num) + "_batch" + str(batch_number) + ".pickle", "wb"))        
        print("batch=" + str(batch_number))
    print("TRAINING: fold=" + str(fold_num) + " has already finished!")
    fold_num += 1
    indices = []
    for i in range(len(kfold_parts)):
        indices.append(next(kfold_parts[i], None))
