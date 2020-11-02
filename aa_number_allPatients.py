import pickle
import os 

patients_properties = "PatientsSequencesInEachRefSeqChain/"
only_mutated = "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
amino_acids_num = {'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'K': 8, 'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17, 'V': 18, 'W': 19, 'Y': 10}
tcga_clinical_dataframe = pickle.load(open("TCGA_clinical_dataframe.pickle","rb"))
patient_AAnumbers = {}
for patient in tcga_clinical_dataframe.index:
    patient_AAnumbers[patient] = []
    with open("PPIs-WithRefSeqsGreaterEqualLength5.txt","rb") as f:
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
            first_chain_mutated = False
            second_chain_mutated = False
            interaction_aa = {}
            if(os.path.exists(patients_properties + patient + "/" + prot_1 + "/" + chain_1)):
                first_chain_mutated = True
                current_ppi_AAnumbers = []
                address1 = patients_properties + patient + "/" + prot_1 + "/" + chain_1 + "/mutated_sequence.seq"
                with open(address1, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        l = l.rstrip()
                        current_ppi_AAnumbers = [amino_acids_num[c] for c in l]
                interaction_aa['first5aa'] = []
                for idx in range(5, len(current_ppi_AAnumbers), 5):
                    interaction_aa['first5aa'].append(current_ppi_AAnumbers[idx - 5: idx])
                    
            if(os.path.exists(patients_properties + patient + "/" + prot_2 + "/" + chain_2)):
                second_chain_mutated = True
                current_ppi_AAnumbers = []
                address2 = patients_properties + patient + "/" + prot_2 + "/" + chain_2 + "/mutated_sequence.seq" 
                with open(address2, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        l = l.rstrip()
                        current_ppi_AAnumbers = [amino_acids_num[c] for c in l]
                interaction_aa['second5aa'] = []
                for idx in range(5, len(current_ppi_AAnumbers), 5):
                    interaction_aa["second5aa"].append(current_ppi_AAnumbers[idx - 5: idx])
                     
            if(first_chain_mutated and not second_chain_mutated):
                current_ppi_AAnumbers = []
                address2 = only_mutated + prot_2 + "/" + chain_2 + "/mutatedTrimedRefSeq.seq"
                with open(address2, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        l = l.rstrip()
                        current_ppi_AAnumbers = [amino_acids_num[c] for c in l]
                interaction_aa['second5aa'] = []
                for idx in range(5, len(current_ppi_AAnumbers), 5):
                    interaction_aa["second5aa"].append(current_ppi_AAnumbers[idx - 5: idx])  
            
            if(second_chain_mutated and not first_chain_mutated):
                current_ppi_AAnumbers = []
                address1 = only_mutated + prot_1 + "/" + chain_1 + "/mutatedTrimedRefSeq.seq"
                with open(address1, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        l = l.rstrip()
                        current_ppi_AAnumbers = [amino_acids_num[c] for c in l]
                interaction_aa['first5aa'] = []
                for idx in range(5, len(current_ppi_AAnumbers), 5):
                    interaction_aa["first5aa"].append(current_ppi_AAnumbers[idx - 5: idx])  
            
            if(first_chain_mutated or second_chain_mutated):
                patient_AAnumbers[patient].append(interaction_aa)
pickle.dump(patient_AAnumbers, open("DATASET/patientsAAnumbers.pickle", "wb"))                                                                                
