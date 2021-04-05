import os
import pickle
import requests
import json
import numpy as np
from selenium import webdriver
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt



class TCGADataPreprocessing:
    def __init__(self, dataset_adr):
        self.dataset_adr = dataset_adr

    #Find the proteins sequences involved in the Protein-Protein Interactions
    def extract_partner_proteins(self):
        all_proteins = []
        with open(self.dataset_adr + "ppis.txt", "rb") as file:
            for linep in file.readlines():
                linep = str(linep,'utf-8')
                prots = linep.split(",")
                prot1 = prots[0].split()[0]
                prot2 = prots[1].split()[0]
                if(prot1 not in all_proteins):
                    all_proteins.append(prot1)
                if(prot2 not in all_proteins):
                    all_proteins.append(prot2)
        if not os.path.exists('Proteins_sequences/'):
            os.makedirs('Proteins_sequences/')
        prot_seqs_adr = "Proteins_sequences/"
        start_to_write = False
        with open(self.dataset_adr + "uniprot-proteomes.fasta", "rb") as file:
            for line in file.readlines():
                line = str(line,'utf-8')
                if(line.find(">sp") != -1):
                    if(line.find("Homo sapiens") != -1 and line.find("GN=") != -1):
                        infos = line.split("GN=")
                        prot_name = infos[1].split()[0]
                        if(prot_name in all_proteins):
                            start_to_write = True
                            all_proteins.remove(prot_name)
                            prot_seq_file = open(prot_seqs_adr + prot_name + ".txt", "w+")
                        else:
                            start_to_write = False
                    else:
                        start_to_write = False
                else:
                    if(line.find(">") == -1):
                        if(start_to_write):
                            prot_seq_file.write(line)

    #Obtain the proteins dictionary with their corresponding interaction ids 
    #The interactions must be have Identities >= 90
    def save_biophysical_PPIs(self, rcsb_filename):
        proteins_interactions_id = {}
        current_protein = ""
        current_interaction = ""
        with open(self.dataset_adr+rcsb_filename, "rb") as file:
            for line in file.readlines():
                line = str(line,'utf-8')
                if(line.find("GN=") != -1):
                    a = line.split("GN=")[1]
                    current_protein = a.split()[0]
                    proteins_interactions_id[current_protein] = []
                if(line.find(">") != -1 and line.find("Chain") != -1):
                    current_interaction = line.split(">")[1].split()[0]
                if(line.find("Identities =") != -1):
                    percent = float(line.split("(")[1].split("%")[0])
                    if(percent > 89.0):
                        proteins_interactions_id[current_protein].append(current_interaction)
        pickle.dump(proteins_interactions_id, open("protein_physical_interactions_90percent.pickle","wb"))

    def get_set(self, l1,l2):
        a1 = []
        for i in l1:
            a1.append(i.split("_")[0])
        a2 = []
        for i in l2:
            a2.append(i.split("_")[0])
        a1 = set(a1)
        a2 = set(a2)
        any_common = a1.intersection(a2)
        commons = []
        for i in l1:
            pdb = i.split("_")[0]
            if(pdb in any_common):
                commons.append(i)
        for i in l2:
            pdb = i.split("_")[0]
            if(pdb in any_common):
                commons.append(i)
        return commons
    
    #We also need to remove those interactions with bit_score < 50 and e_value > 0.000001
    def keep_PPIs_based_on_bit_score_and_evalue(self, oncomine_filename, rcsb_filename):
        proteins_interactions_id = pickle.load(open("protein_physical_interactions_90percent.pickle","rb"))
        oncomine_prots = pickle.load(open(self.dataset_adr+oncomine_filename,"rb"))
        start_to_find_interactions = False
        current_protein = ""
        with open(self.dataset_adr+rcsb_filename, "rb") as file:
            for line in file.readlines():
                line = str(line,'utf-8')
                if(line.find("GN=") != -1):
                    a = line.split("GN=")[1]
                    is_proposed_protein = False
                    start_to_find_interactions = False
                    p = a.split()[0]
                    if(p in proteins_interactions_id.keys()):
                        current_protein = p
                        is_proposed_protein = True
                if(line.find("Sequences producing significant alignments:") != -1 and is_proposed_protein):
                    start_to_find_interactions = True
                if(line.find(">") != -1 and line.find("Chain") != -1):
                    start_to_find_interactions = False
                if(line.find("_") != -1 and start_to_find_interactions):
                    i_informations = line.split()
                    interaction_id = i_informations[0].split("_")[0]
                    bit_score = float(i_informations[len(i_informations)-2])
                    e_value = float(i_informations[len(i_informations)-1])
                    if(interaction_id in proteins_interactions_id[current_protein]):
                        if(bit_score < 50 or e_value > 0.000001):
                            proteins_interactions_id[current_protein].remove(interaction_id)
        interaction_file = open("BLAST_protein_interactions_90.txt","a+")
        num_of_physical_interactions = 0
        proteins = list(proteins_interactions_id.keys())
        for prot1 in range(len(proteins)):
            for prot2 in range(prot1 + 1, len(proteins)):
                any_common = self.get_set(proteins_interactions_id[proteins[prot1]],proteins_interactions_id[proteins[prot2]])
                #if( (proteins[prot1] in oncomine_prots) or (proteins[prot1] in oncomine_prots) ):
                if(len(any_common) > 0):
                    num_of_physical_interactions += 1
                    interaction_file.write(proteins[prot1] + "\t" + proteins[prot2] + "\t")
                    for pdbid in any_common:
                        interaction_file.write(pdbid + "\t")
                    interaction_file.write("\n")    
        interaction_file.close()

    #Create ATOM files and Obtain the Chain sequences for Protein-Protein Interactions 
    #This method sends multiple requests to the RCSB website to obtain the ATOM files
    def create_ATOM_files(self):
        with open(self.dataset_adr+"ppis.txt", "rb") as file:
            for linep in file.readlines():
                linep = str(linep,'utf-8')
                prots = linep.split(",")
                pdb_full_name1 = prots[0].split()[1]
                pdb_full_name2 = prots[1].split()[1]
                pdb_info1 = pdb_full_name1.split("_")
                pdb_id = pdb_info1[0]
                chain1 = pdb_info1[1]
                chain2 = pdb_full_name2.split("_")[1]
                print("ATOM file with PDB_id=" + pdb_id + " has already started!")
                target_url = "http://files.rcsb.org/view/" + pdb_id + ".pdb"
                save_to = "ATOM_files/"
                if not os.path.exists('PDB_sequences'):
                    os.makedirs('PDB_sequences')
                pdb_seq_adr = "PDB_sequences/"
                work_on_atm1 = True
                work_on_atm2 = True
                if(not os.path.exists(save_to + pdb_id)):
                    os.mkdir(save_to + pdb_id)
                if(not os.path.exists(save_to + pdb_id + "/" + chain1)):
                    os.mkdir(save_to + pdb_id + "/" + chain1)
                    atm_file1 = open(save_to + pdb_id + "/" + chain1 + "/" + pdb_id + ".atm", "w+")
                    chain1_seq_file = open(pdb_seq_adr + pdb_id + "_" + chain1 + ".txt", "w+")
                else:
                    work_on_atm1 = False
                if(not os.path.exists(save_to + pdb_id + "/" + chain2)):
                    os.mkdir(save_to + pdb_id + "/" +  chain2)
                    atm_file2 = open(save_to + pdb_id + "/" + chain2 + "/" + pdb_id + ".atm", "w+")
                    chain2_seq_file = open(pdb_seq_adr + pdb_id + "_" + chain2 + ".txt", "w+")
                else:
                    work_on_atm2 = False
                if(work_on_atm1 or work_on_atm2):
                    response = requests.get(target_url)
                    data = response.text
                    atom_number1 = -1
                    atom_number2 = -2
                    continue_until_ter1 = True
                    continue_until_ter2 = True
                    for line in data.split("\n"):
                        infos = line.split()
                        if(len(infos) > 0):
                            if(infos[0] == "SEQRES"):
                                if(infos[2] == chain1 and work_on_atm1):
                                    lineseq = ""
                                    for infos_index in range(4,len(infos)):
                                        lineseq += infos[infos_index] + " "
                                    chain1_seq_file.write(lineseq + "\n")
                                if(infos[2] == chain2 and work_on_atm2):
                                    lineseq = ""
                                    for infos_index in range(4,len(infos)):
                                        lineseq += infos[infos_index] + " "
                                    chain2_seq_file.write(lineseq + "\n")
                            if(infos[0] == "ATOM"):
                                infos4 = infos[4]
                                is_sticked = False
                                if(len(infos[4]) > 1):
                                    infos4 = infos[4][0]
                                    is_sticked = True
                                if(infos4 == chain1 and work_on_atm1 and continue_until_ter1):
                                    atom_number1 = int(infos[1])
                                    if(is_sticked):
                                        wr = infos[0]+"\t"+infos[1]+"\t"+infos[2]+"\t"+infos[3]+"\t"+infos[4][0]+"\t"+infos[4][1:]
                                        for jj in range(5,len(infos)):
                                            wr = wr + "\t" + infos[jj]
                                        atm_file1.write(wr + "\n")
                                    else:
                                        atm_file1.write(line + "\n")
                                if(infos4 == chain2 and work_on_atm2 and continue_until_ter2):
                                    atom_number2 = int(infos[1])
                                    if(is_sticked):
                                        wr = infos[0]+"\t"+infos[1]+"\t"+infos[2]+"\t"+infos[3]+"\t"+infos[4][0]+"\t"+infos[4][1:]
                                        for jj in range(5,len(infos)):
                                            wr = wr + "\t" + infos[jj]
                                        atm_file2.write(wr + "\n")                                
                                    else:
                                        atm_file2.write(line + "\n")
                            if(infos[0] == "TER"):
                                if(int(infos[1]) == atom_number1 + 1 and work_on_atm1 and continue_until_ter1):
                                    if(len(infos[3]) > 1):
                                            atm_file1.write(infos[0]+"\t"+infos[1]+"\t"+infos[2]+"\t"+infos[3][0]+"\t"+infos[3][1:]+"\n")
                                    else:
                                        atm_file1.write(line + "\n")
                                    continue_until_ter1 = False
                                if(int(infos[1]) == atom_number2 + 1 and work_on_atm2 and continue_until_ter2):
                                    if(len(infos[3]) > 1):
                                            atm_file2.write(infos[0]+"\t"+infos[1]+"\t"+infos[2]+"\t"+infos[3][0]+"\t"+infos[3][1:]+"\n")
                                    else:
                                        atm_file2.write(line + "\n")
                                    continue_until_ter2 = False
                    if(work_on_atm1):
                        atm_file1.close()
                        chain1_seq_file.close()
                    if(work_on_atm2):
                        atm_file2.close()
                        chain2_seq_file.close()

    #Find the longest common substring between two strings 
    def lcs(self, S,T):
        m = len(S)
        n = len(T)
        counter = [[0]*(n+1) for x in range(m+1)]
        longest = 0
        lcs_set = set()
        for i in range(m):
            for j in range(n):
                if S[i] == T[j]:
                    c = counter[i][j] + 1
                    counter[i+1][j+1] = c
                    if c > longest:
                        lcs_set = set()
                        longest = c
                        lcs_set.add(S[i-c+1:i+1])
                    elif c == longest:
                        lcs_set.add(S[i-c+1:i+1])
        return lcs_set

    #Extract the start and end positions for those proteins which their corresponding PDB and chain sequence has only one part
    #Something like this: ----XXXXXXX-
    def extract_start_end_positions_one_part(self):
        proteins_refseqs_ids = pickle.load(open(self.dataset_adr + "proteins_refseqs_ids.pickle","rb"))
        o = 0
        indicesfile = open("indices.txt","a+")
        with open(self.dataset_adr + "ppis.txt", "rb") as file:
            for linep in file.readlines():
                linep = str(linep,'utf-8')
                prots = linep.split(",")
                pdb_full_name1 = prots[0].split()[1]
                pdb_full_name2 = prots[1].split()[1]
                pdb_info1 = pdb_full_name1.split("_")
                pdb_id = pdb_info1[0]
                chains = []
                chain1 = pdb_info1[1]
                chain2 = pdb_full_name2.split("_")[1]
                chains.append(chain1)
                chains.append(chain2)
                prots = []
                first_prot = linep.split()[0]
                second_prot = linep.split()[2]
                prots.append(first_prot)
                prots.append(second_prot)
                #Obtain the refseq ids of the proteins
                two_refseq_ids = []
                first_refseq_ids = proteins_refseqs_ids[first_prot]
                second_refseq_ids = proteins_refseqs_ids[second_prot]
                two_refseq_ids.append(first_refseq_ids)
                two_refseq_ids.append(second_refseq_ids)
                for rs in range(0,2):
                    prot_seqs = {}
                    if(len(two_refseq_ids[rs]) > 0):
                        for refseq_id in two_refseq_ids[rs]:
                            prot_seq = ""
                            with open(self.dataset_adr + "RefSeqSequences/" + refseq_id + ".seq", "rb") as seqfile:
                                for lineseq in seqfile.readlines():
                                    lineseq = str(lineseq,'utf-8')
                                    prot_seq = prot_seq + lineseq.rstrip()
                            prot_seqs[refseq_id] = prot_seq
                    else:
                        prot_seq = ""
                        with open(self.dataset_adr + "RefSeqSequences/" + prots[rs] + ".seq", "rb") as seqfile:
                            for lineseq in seqfile.readlines():
                                lineseq = str(lineseq,'utf-8')
                                prot_seq = prot_seq + lineseq.rstrip()
                        prot_seqs[prots[rs]] = prot_seq
                    chain_seq = ""
                    with open(self.dataset_adr + "Available_AminoAcids_in_seq/" + pdb_id + "_" + chains[rs] + ".seq", "rb") as availfile:
                        for lineseq in availfile.readlines():
                            lineseq = str(lineseq,'utf-8')
                            chain_seq += lineseq.rstrip()
                    chain_seqs = chain_seq.split("-")
                    chain_seqs = [elem for elem in chain_seqs if len(elem) > 1]            
                    total_common_seq = 0
                    if(len(chain_seqs) == 1):
                        for refseq_id in prot_seqs:
                            a = prot_seqs[refseq_id]
                            b = chain_seqs[0]
                            lcs1 = list(self.lcs(a,b))[0]
                            #Real index for starting_index=1 instead of 0
                            start_index = a.find(lcs1) + 1
                            end_index = start_index + len(lcs1) - 1
                            indicesfile.write(refseq_id + "\t" + pdb_id + "_" + chains[rs] + "\t" + str(start_index) + "\t" + str(end_index) + "\n")
        indicesfile.close()

    #Extract the start and end positions (more than one <start,end>) for those proteins which their corresponding PDB and chain sequence has more than one part
    #Something like this: ----XXXXXXX---XXX-XXXX---   
    def extract_start_end_positions_more_parts(self):
        o = 0
        indicesfile = open("indices_multipart.txt","a+")
        with open(self.dataset_adr + "ppis.txt", "rb") as file:
            for linep in file.readlines():
                linep = str(linep,'utf-8')
                prots = linep.split(",")
                pdb_full_name1 = prots[0].split()[1]
                pdb_full_name2 = prots[1].split()[1]
                pdb_info1 = pdb_full_name1.split("_")
                pdb_id = pdb_info1[0]
                chains = []
                chain1 = pdb_info1[1]
                chain2 = pdb_full_name2.split("_")[1]
                chains.append(chain1)
                chains.append(chain2)
                prots = []
                first_prot = linep.split()[0]
                second_prot = linep.split()[2]
                prots.append(first_prot)
                prots.append(second_prot)
                #Obtain the refseq ids of the proteins
                two_refseq_ids = []
                first_refseq_ids = proteins_refseqs_ids[first_prot]
                second_refseq_ids = proteins_refseqs_ids[second_prot]
                two_refseq_ids.append(first_refseq_ids)
                two_refseq_ids.append(second_refseq_ids)
                for rs in range(0,2):
                    prot_seqs = {}
                    if(len(two_refseq_ids[rs]) > 0):
                        for refseq_id in two_refseq_ids[rs]:
                            prot_seq = ""
                            with open(self.dataset_adr + "RefSeqSequences/" + refseq_id + ".seq", "rb") as seqfile:
                                for lineseq in seqfile.readlines():
                                    lineseq = str(lineseq,'utf-8')
                                    prot_seq = prot_seq + lineseq.rstrip()
                            prot_seqs[refseq_id] = prot_seq
                    else:
                        prot_seq = ""
                        with open(self.dataset_adr + "RefSeqSequences/" + prots[rs] + ".seq", "rb") as seqfile:
                            for lineseq in seqfile.readlines():
                                lineseq = str(lineseq,'utf-8')
                                prot_seq = prot_seq + lineseq.rstrip()
                        prot_seqs[prots[rs]] = prot_seq
                    chain_seq = ""
                    with open(self.dataset_adr + "Available_AminoAcids_in_seq/" + pdb_id + "_" + chains[rs] + ".seq", "rb") as availfile:
                        for lineseq in availfile.readlines():
                            lineseq = str(lineseq,'utf-8')
                            chain_seq += lineseq.rstrip()
                    chain_seqs = chain_seq.split("-")
                    chain_seqs = [elem for elem in chain_seqs if len(elem) > 1]            
                    total_common_seq = 0
                    if(len(chain_seqs) > 1):
                        for refseq_id in prot_seqs:
                            a = prot_seqs[refseq_id]
                            indicesfile.write(refseq_id + "\t" + pdb_id + "_" + chains[rs])
                            for chain_seq in chain_seqs:
                                b = chain_seq
                                if(len(list(lcs(a,b))) > 0):
                                    lcs1 = list(lcs(a,b))[0]
                                    #Real index for starting_index=1 instead of 0
                                    start_index = a.find(lcs1) + 1
                                    end_index = start_index + len(lcs1) - 1
                                    a = a[end_index : ]
                                    indicesfile.write("\t" + str(start_index) + "\t" + str(end_index))
                            indicesfile.write("\n")
        indicesfile.close() 

    # Organize the indices (Show the exact position of start and ends)
    def organize_indices(self):
        indicesfile = open("indices_multipart_organized.txt","a+")
        with open("indices_multipart.txt", "rb") as file:
            for linep in file.readlines():
                linep = str(linep,'utf-8')
                linep = linep.rstrip()
                infos = linep.split()
                indicesfile.write(infos[0] + "\t" + infos[1] + "\t" + infos[2] + "\t" + infos[3])
                add_me = int(infos[3])
                for i in range(4,len(infos), 2):
                    f = int(infos[i]) + add_me
                    add_me = int(infos[i+1]) + add_me
                    indicesfile.write("\t" + str(f) + "\t" + str(add_me))
                indicesfile.write("\n")
        indicesfile.close()

    #Pickle the RefSeq positions
    def save_RefSeq_positions(self):
        protein_positions = {}
        with open("indices.txt", "rb") as file:
            for linep in file.readlines():
                linep = str(linep,'utf-8')
                linep = linep.rstrip()
                if(linep.find("<") == -1):
                    infos = linep.split()
                    protein_name = infos[0]
                    PDBchain_id = infos[1]
                    if(protein_name not in protein_positions):
                        protein_positions[protein_name] = {}
                    protein_positions[protein_name][PDBchain_id] = []
                    for i in range(2, len(infos)):
                        protein_positions[protein_name][PDBchain_id].append(infos[i])
        pickle.dump(protein_positions, open("RefSeqs_positions.pickle","wb"))

    #Create pdb2rosetta and rosetta2pdb json files (for flex ddG)
    def create_rosetta_json_files(self):
        with open(self.dataset_adr + "ppis.txt", "rb") as fm:
            lines = fm.readlines()
            for linem in lines:
                linem = linem.rstrip()
                linem = str(linem,'utf-8')
                pinfos = linem.split()
                chain1 = pinfos[1].split(",")[0].split("_")[1]
                pdb_chain = pinfos[3].split("_")
                pdb_name = pdb_chain[0]
                chain2 = pdb_chain[1]
                pdb_chains_adr = "ATOM_files/" + pdb_name+"_"+chain1+chain2 + "/"
                if(not os.path.exists(pdb_chains_adr + pdb_name+"_"+chain1+chain2+".pdb")):
                    print("ERROR: The following ATOM file does not exist:")
                    print(pdb_name+"_"+chain1+chain2)
                    print("Therefore, we ignore it!")
                else:
                    if(not os.path.exists(pdb_chains_adr + "pdb2rosetta.resmap.json")):
                        print(pdb_name+"_"+chain1+chain2)
                        with open(pdb_chains_adr + pdb_name+"_"+chain1+chain2+".pdb","rb") as fp:
                            pdb2rosetta_file = open(pdb_chains_adr + "pdb2rosetta.resmap.json","a+")
                            rosetta2pdb_file = open(pdb_chains_adr + "rosetta2pdb.resmap.json","a+")
                            alllines = fp.readlines()
                            previous_num = "DifferentThing"
                            number_to_assign = 1
                            pdb2rosetta_content = "{\n"
                            rosetta2pdb_content = "{\n"
                            for chainsline in alllines:
                                chainsline = chainsline.rstrip()
                                chainsline = str(chainsline,'utf-8')
                                current_num = chainsline[22:27].strip()
                                if(current_num != previous_num):
                                    pdb2rosetta_content += "\"" + chainsline[21:27] + "\": " + str(number_to_assign) + ",\n"
                                    rosetta2pdb_content += "\"" + str(number_to_assign) + "\": " + "\"" + chainsline[21:27] + "\",\n"
                                    number_to_assign += 1
                                previous_num = current_num

                            pdb2rosetta_file.write(pdb2rosetta_content[:len(pdb2rosetta_content)-2] + "\n}")
                            pdb2rosetta_file.close()
                            rosetta2pdb_file.write(rosetta2pdb_content[:len(rosetta2pdb_content)-2] + "\n}")
                            rosetta2pdb_file.close()

    #Prepare 'resfile' and 'mutfile' files in folder 'PatientsPPImutationsData'
    #It is for extracting the mutations and the exactly position of the mutations in two different standard positions
    #1. Rosetta 2. pdb ATOM position
    #Notice: Before using this function you must prepare the output annotated files using ANNOVAR
    def prepare_rosetta_files(self, tcga_annovar_folder):
        if not os.path.exists('PatientsPPImutationsData/'):
            os.makedirs('PatientsPPImutationsData/')
        target_folder = "PatientsPPImutationsData/"
        chains_start_numbers_rosetta = pickle.load(open(self.dataset_adr + "chains_start_numbers_rosetta.pickle", "rb"))
        exact_positions_protein_refseq = pickle.load(open(self.dataset_adr + "exact_positions_protein_refseq.pickle","rb"))
        dirr = tcga_annovar_folder
        for tcga in os.listdir(dirr):
            os.mkdir(target_folder + tcga)
            folder = dirr + tcga + "/annovar_annotates/"
            for ann in os.listdir(folder):
                if(ann.find(".exonic") != -1):
                    patient_id = ann.split(".exonic")[0]
                    patient_resfile = {}
                    patient_mutfile = {}
                    with open(folder + ann, "rb") as fex:
                        for lineex in fex.readlines():
                            lineex = str(lineex,'utf-8')
                            lineex = lineex.rstrip()
                            mut_line = lineex.split()[3]
                            if(mut_line.find("p.") != -1 and lineex.split()[2] == "SNV"):
                                findings = mut_line.split(":")
                                prot = findings[0]
                                if(prot in exact_positions_protein_refseq.keys()):
                                    dd = mut_line.split(prot + ":")
                                    dd.remove('')
                                    for aa_mut in dd:
                                        if(aa_mut.find("p.") != -1):
                                            allchanges = aa_mut.split(":")
                                            current_refseq = allchanges[0]
                                            if('NOavailableRefSeq' not in exact_positions_protein_refseq[prot]):
                                                if(current_refseq in exact_positions_protein_refseq[prot]):
                                                    change = allchanges[len(allchanges)-1].split(",")[0][2:]
                                                    from_aa = change[0]
                                                    to_aa = change[len(change)-1]
                                                    pos = int(change[1:len(change)-1])
                                                    if(from_aa != to_aa):
                                                        pdbchain_infos = exact_positions_protein_refseq[prot][current_refseq]
                                                        for PDBchain in pdbchain_infos:
                                                            pdb_infos = PDBchain.split("_")
                                                            pdb_id = pdb_infos[0]
                                                            chain = pdb_infos[1]

                                                            NM_seq = ""
                                                            with open(self.dataset_adr + "RefSeqSequences/" + current_refseq + ".seq", "rb") as refseqf:
                                                                liners = refseqf.readlines()
                                                                for aline in liners:
                                                                    aline = aline.rstrip()
                                                                    aline = str(aline,'utf-8')
                                                                    NM_seq += aline
                                                            chain_part_seq = ""
                                                            with open(self.dataset_adr + "Available_AminoAcids_in_seq/" + PDBchain + ".seq", "rb") as pcf:
                                                                seq = ""
                                                                liners = pcf.readlines()
                                                                for aline in liners:
                                                                    aline = aline.rstrip()
                                                                    aline = str(aline,'utf-8')
                                                                    seq += aline
                                                                chain_seqs = seq.split("-")
                                                                chain_seqs = [elem for elem in chain_seqs if len(elem) > 1]
                                                            previous_distance = 0
                                                            base_index = 0
                                                            for start_end in range(1, len(pdbchain_infos[PDBchain]), 2):
                                                                start = int(pdbchain_infos[PDBchain][start_end - 1])
                                                                end = int(pdbchain_infos[PDBchain][start_end])
                                                                NM_seq = NM_seq[start-1:end-1]
                                                                chain_part_seq = chain_seqs[int(start_end/2)]
                                                                base_index += chain_part_seq.find(NM_seq)
                                                                if(chain_part_seq.find(NM_seq) == -1):
                                                                    print("ERR: Could not find the sequence!")
                                                                if(pos >= start and pos <= end):
                                                                    #Index which starts from 1 (NOT 0)
                                                                    exact_index_current_part = pos - start + 1
                                                                    exact_index = exact_index_current_part + previous_distance

                                                                    for chains in chains_start_numbers_rosetta[pdb_id]:
                                                                        if(chain in chains):
                                                                            #Chain Variable =>  chain
                                                                            #Original AMINO ACID Variable => aminoacid
                                                                            #ATOM id in the PDB file => atom_id
                                                                            mutfile_position = chains_start_numbers_rosetta[pdb_id][chains][chain] + exact_index + base_index - 1
                                                                            #Original AMINO ACID
                                                                            aminoacid = ""
                                                                            atom_id = ""
                                                                            with open('ATOM_files/' + pdb_id + "_" + chains + '/pdb2rosetta_complete.json') as json_file:
                                                                                rosetta2pdb = json.load(json_file)
                                                                                atom_id_aminoacid = rosetta2pdb[str(mutfile_position)][1:].split()
                                                                                atom_id = atom_id_aminoacid[0]
                                                                                aminoacid = atom_id_aminoacid[1]
                                                                            #ResultsAminoAcid (after mutation)
                                                                            #resfile structure: AtomId chainId PIKAA ResultAminoAcid
                                                                            #mutfile structure: OriginalAminoAcid rosettaId ResultAminoAcid
                                                                            if(aminoacid != to_aa):
                                                                                if(pdb_id + "_" + chains not in patient_mutfile):
                                                                                    patient_mutfile[pdb_id + "_" + chains] = []
                                                                                writein = False
                                                                                if(aminoacid + " " + str(mutfile_position) + " " + to_aa not in patient_mutfile[pdb_id + "_" + chains]):
                                                                                    patient_mutfile[pdb_id + "_" + chains].append(aminoacid + " " + str(mutfile_position) + " " + to_aa)
                                                                                    writein = True
                                                                                if(pdb_id + "_" + chains not in patient_resfile):
                                                                                    patient_resfile[pdb_id + "_" + chains] = []
                                                                                if(writein):
                                                                                    patient_resfile[pdb_id + "_" + chains].append(atom_id + " " + chain + " PIKAA " + to_aa)
                                                                previous_distance += end - start + 1
                                            else:
                                                change = allchanges[len(allchanges)-1].split(",")[0][2:]
                                                from_aa = change[0]
                                                to_aa = change[len(change)-1]
                                                pos = int(change[1:len(change)-1])
                                                if(from_aa != to_aa):
                                                    pdbchain_infos = exact_positions_protein_refseq[prot]['NOavailableRefSeq']
                                                    for PDBchain in pdbchain_infos:
                                                        pdb_infos = PDBchain.split("_")
                                                        pdb_id = pdb_infos[0]
                                                        chain = pdb_infos[1]
                                                        NM_seq = ""
                                                        with open(self.dataset_adr + "RefSeqSequences/" + prot + ".seq", "rb") as refseqf:
                                                            liners = refseqf.readlines()
                                                            for aline in liners:
                                                                aline = aline.rstrip()
                                                                aline = str(aline,'utf-8')
                                                                NM_seq += aline

                                                        chain_part_seq = ""
                                                        with open(self.dataset_adr + "Available_AminoAcids_in_seq/" + PDBchain + ".seq", "rb") as pcf:
                                                            seq = ""
                                                            liners = pcf.readlines()
                                                            for aline in liners:
                                                                aline = aline.rstrip()
                                                                aline = str(aline,'utf-8')
                                                                seq += aline
                                                            chain_seqs = seq.split("-")
                                                            chain_seqs = [elem for elem in chain_seqs if len(elem) > 1]
                                                        previous_distance = 0
                                                        base_index = 0
                                                        for start_end in range(1, len(pdbchain_infos[PDBchain]), 2):
                                                            start = int(pdbchain_infos[PDBchain][start_end - 1])
                                                            end = int(pdbchain_infos[PDBchain][start_end])
                                                            NM_seq = NM_seq[start-1:end-1]
                                                            chain_part_seq = chain_seqs[int(start_end/2)]
                                                            base_index += chain_part_seq.find(NM_seq)
                                                            if(chain_part_seq.find(NM_seq) == -1):
                                                                print("ERR: Could not find the sequence!")
                                                            if(pos >= start and pos <= end):
                                                                #Index which starts from 1 (NOT 0)
                                                                exact_index_current_part = pos - start + 1
                                                                exact_index = exact_index_current_part + previous_distance
                                                                for chains in chains_start_numbers_rosetta[pdb_id]:
                                                                    if(chain in chains):
                                                                        #Chain Variable =>  chain
                                                                        #Original AMINO ACID Variable => aminoacid
                                                                        #ATOM id in the PDB file => atom_id
                                                                        mutfile_position = chains_start_numbers_rosetta[pdb_id][chains][chain] + exact_index + base_index - 1
                                                                        #Original AMINO ACID
                                                                        aminoacid = ""
                                                                        atom_id = ""
                                                                        with open('ATOM_files/' + pdb_id + "_" + chains + '/pdb2rosetta_complete.json') as json_file:
                                                                            rosetta2pdb = json.load(json_file)
                                                                            atom_id_aminoacid = rosetta2pdb[str(mutfile_position)][1:].split()
                                                                            atom_id = atom_id_aminoacid[0]
                                                                            aminoacid = atom_id_aminoacid[1]
                                                                        #ResultsAminoAcid (after mutation)
                                                                        #resfile structure: AtomId chainId PIKAA ResultAminoAcid
                                                                        #mutfile structure: OriginalAminoAcid rosettaId ResultAminoAcid
                                                                        if(aminoacid != to_aa):
                                                                            if(pdb_id + "_" + chains not in patient_mutfile):
                                                                                patient_mutfile[pdb_id + "_" + chains] = []
                                                                            writein = False
                                                                            if(aminoacid + " " + str(mutfile_position) + " " + to_aa not in patient_mutfile[pdb_id + "_" + chains]):
                                                                                patient_mutfile[pdb_id + "_" + chains].append(aminoacid + " " + str(mutfile_position) + " " + to_aa)
                                                                                writein = True
                                                                            if(pdb_id + "_" + chains not in patient_resfile):
                                                                                patient_resfile[pdb_id + "_" + chains] = []
                                                                            if(writein):
                                                                                patient_resfile[pdb_id + "_" + chains].append(atom_id + " " + chain + " PIKAA " + to_aa)

                                                            previous_distance += end - start + 1
                    pickle.dump(patient_mutfile, open(target_folder + tcga + "/" + patient_id + "_mutfile.pickle","wb"))
                    pickle.dump(patient_resfile, open(target_folder + tcga + "/" + patient_id + "_resfile.pickle","wb"))

    #Each chain in an PDB complex has a starting Rosetta number
    def PDB_chain_starting_Rosetta_Num(self):
        chains_start_numbers_rosetta = {}
        main_folder = "ATOM_files/"
        for folder in os.listdir(main_folder):
            pdb = folder.split("_")[0]
            chains = folder.split("_")[1]
            if(pdb not in chains_start_numbers_rosetta):
                chains_start_numbers_rosetta[pdb] = {}
            chains_start_numbers_rosetta[pdb][chains] = {}
            with open(main_folder + folder + "/pdb2rosetta.resmap.json") as p2r_file:
                prev_character = "qq"
                for liner in p2r_file.readlines():
                    liner = liner.rstrip()
                    if(liner.find(":") != -1):
                        number = liner.split(":")[1].split(",")[0].strip()
                        if(prev_character != liner[1]):
                            chains_start_numbers_rosetta[pdb][chains][liner[1]] = int(number)
                        prev_character = liner[1]
        #chains_start_numbers_rosetta
        pickle.dump(chains_start_numbers_rosetta, open("chains_start_numbers_rosetta.pickle","wb"))

    def extract_exact_positions(self):
        avail_prots = []
        with open(self.dataset_adr + "ppis.txt", "rb") as file:
            for linep in file.readlines():
                linep = str(linep,'utf-8')
                first_prot = linep.split()[0]
                second_prot = linep.split()[2]
                if(first_prot not in avail_prots):
                    avail_prots.append(first_prot)
                if(second_prot not in avail_prots):
                    avail_prots.append(second_prot)
        bb = pickle.load(open(self.dataset_adr + "proteins_refseqs_ids.pickle","rb"))
        protein_refseq_positions = pickle.load(open("RefSeqs_positions.pickle","rb"))
        all_protein_refseqs_positions = {}
        for prot in bb:
            if(prot in avail_prots):
                all_protein_refseqs_positions[prot] = {}
                if(len(bb[prot]) == 0):
                    all_protein_refseqs_positions[prot]["NOavailableRefSeq"] = protein_refseq_positions[prot]
                else:
                    for refseq in bb[prot]:
                        if(refseq in protein_refseq_positions):
                            all_protein_refseqs_positions[prot][refseq] = protein_refseq_positions[refseq]
        pickle.dump(all_protein_refseqs_positions,open("exact_positions_protein_refseq.pickle","wb"))

    #Extract RefSeqGene ID protein (amino acid) sequences from NCBI website
    def extract_RefSeqID_sequences_from_NCBI(self):
        refseq_adr = "RefSeqSequences/"
        all_without_seq = "all_refseqs_without_seq.txt"
        fall = open(all_without_seq,"w+")
        with open(self.dataset_adr + "ProteinsFullInformation.txt","rb") as allrefseqfile:
            refseq_lines = allrefseqfile.readlines()
            for refseq_line in refseq_lines:
                refseq_line = str(refseq_line,'utf-8')
                refseq_line = refseq_line.rstrip()
                allpinfos = refseq_line.split()
                if(len(allpinfos) == 3):
                    refseq_id = allpinfos[2]
                    driver = webdriver.Firefox()
                    driver.get("https://www.ncbi.nlm.nih.gov/nuccore/" + refseq_id)
                    element = driver.find_elements_by_class_name('feature')
                    is_not_inside = True
                    seq = ""
                    for e in element:
                        a = str(e.text)
                        if(a.find("/translation=\"") != -1):
                            seq = a.split("/translation=\"")[1].split("\"")[0]
                            seq = seq.replace(" ", "")
                            seq = seq.replace("\n", "")
                            is_not_inside = False
                    if(is_not_inside):
                        print(refseq_id)
                        fall.write(refseq_id + "\n")
                    else:
                        refseqfile = open(refseq_adr + refseq_id + ".seq","w+")
                        refseqfile.write(seq)
                        refseqfile.close()
                    driver.close()
                elif(len(allpinfos) > 3):
                    for refseq_index in range(2,len(allpinfos)):
                        refseq_id = allpinfos[refseq_index].split(",")[0]
                        driver = webdriver.Firefox()
                        driver.get("https://www.ncbi.nlm.nih.gov/nuccore/" + refseq_id)
                        element = driver.find_elements_by_class_name('feature')
                        is_not_inside = True
                        seq = ""
                        for e in element:
                            a = str(e.text)
                            if(a.find("/translation=\"") != -1):
                                seq = a.split("/translation=\"")[1].split("\"")[0]
                                seq = seq.replace(" ", "")
                                seq = seq.replace("\n", "")
                                is_not_inside = False
                        if(is_not_inside):
                            print(refseq_id)
                            fall.write(refseq_id + "\n")
                        else:
                            refseqfile = open(refseq_adr + refseq_id + ".seq","w+")
                            refseqfile.write(seq)
                            refseqfile.close()
                        driver.close()
        fall.close()

    #To report the mutations occurred in each PDB complex
    #Create files for mutations.resfile, mutations.mutfile, nataa_mutations.resfile, and chains_to_move.txt for each PDB complex
    def create_mutfile_resfile_Rosetta(self):
        target_folder = "PatientsPPImutationsData/"
        for tcga in os.listdir(target_folder):
            addr = target_folder + tcga + "/"
            for mutfile in os.listdir(addr):
                if(mutfile.find("mutfile") != -1):
                    resfile = mutfile.replace("mutfile", "resfile")
                    mutfile_dictionary = pickle.load(open(addr + mutfile,"rb"))
                    resfile_dictionary = pickle.load(open(addr + resfile,"rb"))
                    #Consider all mutated interactions
                    for interaction in resfile_dictionary:
                        chains_included = []
                        chains_to_move = ""
                        counter = 0
                        for eachline in mutfile_dictionary[interaction]:
                            counter += 1
                            mutations_mutfile += eachline + "\n"
                        mutations_mutfile = "total " + str(counter) + "\n" + str(counter) + "\n" + mutations_mutfile
                        nataa_resfile = "NATAA\n"
                        nataa_resfile += "start\n"
                        mutations_resfile = "NATRO\n"
                        mutations_resfile += "start\n"
                        for eachline in resfile_dictionary[interaction]:
                            mutations_resfile += eachline + "\n"
                            nataa_resfile += eachline + "\n"
                            chain = eachline.split()[1]
                            if(chain not in chains_included):
                                chains_included.append(chain)
                        if(len(chains_included) > 1):
                            chains_to_move = chains_included[0] + "," + chains_included[1]
                        else:
                            chains_to_move = chains_included[0]
                        #For one interaction in a patient    
                        mutations_resfile = mutations_resfile[:-1]
                        nataa_resfile = nataa_resfile[:-1]
                        mutations_mutfile = mutations_mutfile[:-1]

    #Extract the positions of all amino acids which are mutated in all patients  (These position are the positions in the available amino acids in PDB file)
    #The positions in "mut_positions_RefSeq" are sorted!
    #We will keep only these position for our analysis process
    #For example: original_AA_sequence => ABCDEFGHIJKLMN and mut_positions_RefSeq => [1,3,6,8,9]
    #The result will be => ACFHI
    #Notice: Before Using this method, you must obtain annotated files using ANNOVAR 
    def extract_position_mutated_aminoacids(self, tcga_annovar_folder):
        exact_positions_protein_refseq = pickle.load(open(self.dataset_adr + "exact_positions_protein_refseq.pickle","rb"))
        mut_positions_RefSeq = {}
        dirr = tcga_annovar_folder
        for tcga in os.listdir(dirr):
            folder = dirr + tcga + "/annovar_annotates/"
            for ann in os.listdir(folder):
                if(ann.find(".exonic") != -1):
                    patient_id = ann.split(".exonic")[0]
                    with open(folder + ann, "rb") as fex:
                        for lineex in fex.readlines():
                            lineex = str(lineex,'utf-8')
                            lineex = lineex.rstrip()
                            mut_line = lineex.split()[3]
                            if(mut_line.find("p.") != -1 and lineex.split()[2] == "SNV"):
                                findings = mut_line.split(":")
                                prot = findings[0]
                                if(prot in exact_positions_protein_refseq.keys()):
                                    dd = mut_line.split(prot + ":")
                                    dd.remove('')
                                    for aa_mut in dd:
                                        if(aa_mut.find("p.") != -1):
                                            allchanges = aa_mut.split(":")
                                            current_refseq = allchanges[0]
                                            if('NOavailableRefSeq' not in exact_positions_protein_refseq[prot]):
                                                if(current_refseq in exact_positions_protein_refseq[prot]):
                                                    change = allchanges[len(allchanges)-1].split(",")[0][2:]
                                                    from_aa = change[0]
                                                    to_aa = change[len(change)-1]
                                                    pos = int(change[1:len(change)-1])
                                                    if(from_aa != to_aa):
                                                        pdbchain_infos = exact_positions_protein_refseq[prot][current_refseq]
                                                        for PDBchain in pdbchain_infos:
                                                            pdb_infos = PDBchain.split("_")
                                                            pdb_id = pdb_infos[0]
                                                            chain = pdb_infos[1]
                                                            previous_distance = 0
                                                            base_index = 0
                                                            for start_end in range(1, len(pdbchain_infos[PDBchain]), 2):
                                                                start = int(pdbchain_infos[PDBchain][start_end - 1])
                                                                end = int(pdbchain_infos[PDBchain][start_end])
                                                                if(pos >= start and pos <= end):
                                                                    #Index which starts from 1 (NOT 0)
                                                                    exact_index_current_part = pos - start + 1
                                                                    exact_index = exact_index_current_part + previous_distance

                                                                    if(current_refseq not in mut_positions_RefSeq):
                                                                        mut_positions_RefSeq[current_refseq] = {}
                                                                    if(PDBchain not in mut_positions_RefSeq[current_refseq]):
                                                                        mut_positions_RefSeq[current_refseq][PDBchain] = []
                                                                    if(exact_index not in mut_positions_RefSeq[current_refseq][PDBchain]):
                                                                        mut_positions_RefSeq[current_refseq][PDBchain].append(exact_index)
                                                                previous_distance += end - start + 1                
                                            else:
                                                change = allchanges[len(allchanges)-1].split(",")[0][2:]
                                                from_aa = change[0]
                                                to_aa = change[len(change)-1]
                                                pos = int(change[1:len(change)-1])
                                                if(from_aa != to_aa):
                                                    pdbchain_infos = exact_positions_protein_refseq[prot]['NOavailableRefSeq']
                                                    for PDBchain in pdbchain_infos:
                                                        pdb_infos = PDBchain.split("_")
                                                        pdb_id = pdb_infos[0]
                                                        chain = pdb_infos[1]
                                                        previous_distance = 0
                                                        base_index = 0
                                                        for start_end in range(1, len(pdbchain_infos[PDBchain]), 2):
                                                            start = int(pdbchain_infos[PDBchain][start_end - 1])
                                                            end = int(pdbchain_infos[PDBchain][start_end])
                                                            if(pos >= start and pos <= end):
                                                                #Index which starts from 1 (NOT 0)
                                                                exact_index_current_part = pos - start + 1
                                                                exact_index = exact_index_current_part + previous_distance
                                                                if(prot not in mut_positions_RefSeq):
                                                                    mut_positions_RefSeq[prot] = {}
                                                                if(PDBchain not in mut_positions_RefSeq[prot]):
                                                                    mut_positions_RefSeq[prot][PDBchain] = []
                                                                if(exact_index not in mut_positions_RefSeq[prot][PDBchain]):
                                                                    mut_positions_RefSeq[prot][PDBchain].append(exact_index)
                                                            previous_distance += end - start + 1                                    
        for refSeq in mut_positions_RefSeq:
            for PDBchain in mut_positions_RefSeq[refSeq]:
                mut_positions_RefSeq[refSeq][PDBchain].sort()            
        pickle.dump(mut_positions_RefSeq, open("Mutation_Positions_In_Trimed_RefSeqs.pickle","wb"))

    #Trim the sequences (Keep those Amino Acids which their ATOMs has reported in PDB)
    def keep_aa_reported_in_ATOM_files(self):
        exact_positions_protein_refseq = pickle.load(open(self.dataset_adr + "exact_positions_protein_refseq.pickle","rb"))
        destination = "RefSeq-TheSeqsContributeInPPI/"
        for protein in exact_positions_protein_refseq:
            the_prot_allRefSeqs = exact_positions_protein_refseq[protein]
            for current_refseq in the_prot_allRefSeqs:
                NM_seq = ""
                adr_folder_refseq = ""
                if(current_refseq != "NOavailableRefSeq"):
                    adr_folder_refseq = destination + current_refseq
                    if(not os.path.exists(adr_folder_refseq)):
                        os.mkdir(adr_folder_refseq)
                    with open(self.dataset_adr + "RefSeqSequences/" + current_refseq + ".seq", "rb") as refseqf:
                        liners = refseqf.readlines()
                        for aline in liners:
                            aline = aline.rstrip()
                            aline = str(aline,'utf-8')
                            NM_seq += aline
                else:
                    adr_folder_refseq = destination + protein
                    if(not os.path.exists(adr_folder_refseq)):
                        os.mkdir(adr_folder_refseq)
                    with open(self.dataset_adr + "RefSeqSequences/" + protein + ".seq", "rb") as refseqf:
                        liners = refseqf.readlines()
                        for aline in liners:
                            aline = aline.rstrip()
                            aline = str(aline,'utf-8')
                            NM_seq += aline
                the_refSeq_allChains = the_prot_allRefSeqs[current_refseq]        
                for pdb_chain in the_refSeq_allChains:
                    final_adr = adr_folder_refseq + "/" + pdb_chain + "/"
                    if(not os.path.exists(final_adr)):
                        os.mkdir(final_adr)     
                    file_trim_seq = open(final_adr + "trimedRefSeq_seq.seq","w")
                    pdbchain_ranges = the_refSeq_allChains[pdb_chain]
                    prot_sequence = ""
                    for start_end in range(1, len(pdbchain_ranges), 2):
                        start = int(pdbchain_ranges[start_end - 1])
                        end = int(pdbchain_ranges[start_end])
                        prot_sequence += NM_seq[start-1:end]
                    file_trim_seq.write(prot_sequence)
                    file_trim_seq.close()

    #We extract amino acids in the specific positions reported in "mut_positions_RefSeq"
    #"mut_positions_RefSeq" only keeps the positions of the amino acids which mutated in at least one patient
    def extract_aa_mutated(self):
        refseqs_seq_folder = self.dataset_adr + "RefSeq-TheSeqsContributeInPPI/"
        refseqs_only_mutated_aa_folder = "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
        mut_positions_RefSeq = pickle.load(open(self.dataset_adr + "Mutation_Positions_In_Trimed_RefSeqs.pickle","rb"))
        for refseq in mut_positions_RefSeq:
            for PDBchain in mut_positions_RefSeq[refseq]:
                if(not os.path.exists(refseqs_only_mutated_aa_folder + refseq)):
                    os.mkdir(refseqs_only_mutated_aa_folder + refseq)
                if(not os.path.exists(refseqs_only_mutated_aa_folder + refseq + "/" + PDBchain)):
                    os.mkdir(refseqs_only_mutated_aa_folder + refseq + "/" + PDBchain)
                mutated_seq_file_adr = refseqs_only_mutated_aa_folder + refseq + "/" + PDBchain + "/mutatedTrimedRefSeq.seq" 
                if(not os.path.exists(mutated_seq_file_adr)):
                    NM_seq = ""
                    with open(refseqs_seq_folder + refseq + "/" + PDBchain + "/trimedRefSeq_seq.seq", "rb") as refseqf:
                        liners = refseqf.readlines()
                        for aline in liners:
                            aline = aline.rstrip()
                            aline = str(aline,'utf-8')
                            NM_seq += aline       
                    original_mutatedAminoAcids = ""
                    for aa_pos in mut_positions_RefSeq[refseq][PDBchain]:
                        original_mutatedAminoAcids += NM_seq[aa_pos-1]
                    file = open(mutated_seq_file_adr,"w+")
                    file.write(str(original_mutatedAminoAcids))
                    file.close()

    #Obtain Each Patient final sequence after mutations (just mutated sequences) to feed as a sequence to our model
    #If there is no mutation occurred in a PPI, then we do not save its sequence in "PatientsSequencesInEachRefSeqChain" folder
    #In the case that there is no mutation in a PPI, you must use its sequence in "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients" folder
    def obtain_patients_sequences(self):
        target_folder = "PatientsSequencesInEachRefSeqChain/"
        refseqs_only_mutated_aa_folder = "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
        forEachPatient_mut_positions_RefSeq = pickle.load(open(self.dataset_adr + "Mutation_Positions_In_Trimed_RefSeqs_ForEachPatient.pickle","rb"))
        mut_positions_RefSeq = pickle.load(open(self.dataset_adr + "Mutation_Positions_In_Trimed_RefSeqs.pickle","rb"))
        for patient_id in forEachPatient_mut_positions_RefSeq:
            if(not os.path.exists(target_folder + patient_id)):
                    os.mkdir(target_folder + patient_id)
            for refseq in forEachPatient_mut_positions_RefSeq[patient_id]:
                if(not os.path.exists(target_folder + patient_id + "/" + refseq)):
                    os.mkdir(target_folder + patient_id + "/" + refseq)
                for PDBchain in forEachPatient_mut_positions_RefSeq[patient_id][refseq]:
                    if(not os.path.exists(target_folder + patient_id + "/" + refseq + "/" + PDBchain)):
                        os.mkdir(target_folder + patient_id + "/" + refseq + "/" + PDBchain)
                    NM_seq = ""
                    with open(refseqs_only_mutated_aa_folder + refseq + "/" + PDBchain + "/mutatedTrimedRefSeq.seq", "rb") as refseqf:
                        liners = refseqf.readlines()
                        for aline in liners:
                            aline = aline.rstrip()
                            aline = str(aline,'utf-8')
                            NM_seq += aline
                    aminoacid_mutations = forEachPatient_mut_positions_RefSeq[patient_id][refseq][PDBchain]
                    for aa_mutation in aminoacid_mutations:
                        to_aa = aa_mutation[len(aa_mutation)-1]        
                        original_index = int(aa_mutation[1:len(aa_mutation)-1])
                        final_seq_index = mut_positions_RefSeq[refseq][PDBchain].index(original_index)
                        NM_seq = NM_seq[:final_seq_index] + to_aa + NM_seq[final_seq_index + 1:]
                    fname = target_folder + patient_id + "/" + refseq + "/" + PDBchain + "/mutated_sequence.seq"
                    if(not os.path.exists(fname)):
                        file = open(fname,"w+")
                        file.write(str(NM_seq))
                        file.close()

    #Extract Mutation Windows with high Standard Deviation of the frequency of mutations occur in different cancer types
    def extract_important_windows(self):
        only_mutated = "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
        patients_properties = "PatientsSequencesInEachRefSeqChain/"
        tcga_clinical_dataframe = pickle.load(open(self.dataset_adr + "TCGA_clinical_dataframe.pickle","rb"))
        classes = {"adrenal gland":0, "bladder":1, "breast":2, "GYN":3, "bile duct":4, "CRC":5, "bone marrow":6, "esophagus":7, "brain":8, "head and neck":9, "Kidney":10, "liver":11, "Lung":12, "pleura":13, "pancreatic":14, "male reproductive system":15, "other":16, "Melanoma":17, "stomach":18, "thyroid":19, "thymus":20}
        classes_map = {"adrenal gland": "AdrGl", "bladder":"Blad", "breast":"brst", "GYN":"GYN", "bile duct":"BileD", "CRC":"CRC", "bone marrow":"BoneM", "esophagus":"Esoph", "brain":"Brain", "head and neck":"H&N", "Kidney":"Kidney", "liver":"Liver", "Lung":"Lung", "pleura":"Pleu", "pancreatic":"Panc", "male reproductive system":"MaleRep", "other":"Other", "Melanoma":"Melan", "stomach":"Stomc", "thyroid":"Thyr", "thymus":"Thym"}
        patients = {}
        for cancer_class in classes:
            patients[cancer_class] = tcga_clinical_dataframe[tcga_clinical_dataframe['cancer_class'] == cancer_class].index
        with open(self.dataset_adr + "PPIs-WithRefSeqs.txt","rb") as f:
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
            #CHAIN-1    
                address = only_mutated + prot_1 + "/" + chain_1 + "/mutatedTrimedRefSeq.seq"
                original_seq = ""
                with open(address, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        original_seq = l.rstrip()
            #CHAIN-2
                address = only_mutated + prot_2 + "/" + chain_2 + "/mutatedTrimedRefSeq.seq"
                with open(address, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        original_seq += l.rstrip()
                ppi_dic = {}
                ppi_dic['Window'] = []
                ppi_dic['Mutation Accumulation'] = []
                ppi_dic['Cancer Class'] = []
                final_counts = []
                for cancer_class in classes:
                    counts = []
                    j = 0
                    for index in range(0, len(original_seq), window_size):
                        counts.append(0)
                        ppi_dic['Cancer Class'].append(cancer_class)
                        ppi_dic['Window'].append(j)
                        j += 1
                    for patient in patients[cancer_class]:
                        patient_seq = ""
                        if(os.path.exists(patients_properties + patient + "/" + prot_1 + "/" + chain_1)):
                            address = patients_properties + patient + "/" + prot_1 + "/" + chain_1 + "/mutated_sequence.seq"
                            with open(address, "rb") as seqf:
                                for l in seqf.readlines():
                                    l = str(l,'utf-8')
                                    patient_seq = l.rstrip()
                        else:
                            address = only_mutated + prot_1 + "/" + chain_1 + "/mutatedTrimedRefSeq.seq"
                            with open(address, "rb") as seqf:
                                for l in seqf.readlines():
                                    l = str(l,'utf-8')
                                    patient_seq = l.rstrip()

                        if(os.path.exists(patients_properties + patient + "/" + prot_2 + "/" + chain_2)):
                            address = patients_properties + patient + "/" + prot_2 + "/" + chain_2 + "/mutated_sequence.seq"
                            with open(address, "rb") as seqf:
                                for l in seqf.readlines():
                                    l = str(l,'utf-8')
                                    patient_seq += l.rstrip()
                        else:
                            address = only_mutated + prot_2 + "/" + chain_2 + "/mutatedTrimedRefSeq.seq"
                            with open(address, "rb") as seqf:
                                for l in seqf.readlines():
                                    l = str(l,'utf-8')
                                    patient_seq += l.rstrip()
                        #Count the number of mutations in each windows
                        j = 0
                        for index in range(0, len(original_seq), window_size):
                            original_piece = original_seq[index : index + window_size]
                            patient_piece = patient_seq[index : index + window_size]
                            for i in range(len(original_piece)):
                                if(original_piece[i] != patient_piece[i]):
                                    counts[j] += 1
                                    break
                            j += 1
                    counts = np.array(counts)
                    counts = counts / len(patients[cancer_class])
                    final_counts.append(counts)
                final_counts = np.array(final_counts)
                ppi_dic['Mutation Accumulation'] = final_counts.flatten()
                df = pd.DataFrame(data=ppi_dic)
                #plt.figure(figsize=(20,12))
                #ax = sns.lineplot(x="Window", y="Mutation Accumulation", hue="Cancer Class", data=df, palette=["#7B241C","#E74C3C","#C39BD3","#6C3483","#2980B9","#AED6F1","#48C9B0","#117A65","#52BE80","#ABEBC6","#F7DC6F","#F5B041","#E67E22","#BA4A00","#B3B6B7","#797D7F","#839192","#515A5A","#5D6D7E","#17202A","#D7BDE2"])
                #lgd = plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), ncol=1)
                #fig = ax.get_figure()
                #fig.savefig('MutationAccumulations/' + prot_1 + '_' + chain_1 + '_and_' + prot_2 + '_' + chain_2 + '.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
                std_df = df.groupby("Window").std().reset_index().rename({'Mutation Accumulation': 'Standard Deviation'}, axis='columns')
                #ax = sns.lineplot(x="Window", y="Standard Deviation", data=std_df)
                #fig = ax.get_figure()
                #fig.savefig('MutationAccumulations/' + prot_1 + '_' + chain_1 + '_and_' + prot_2 + '_' + chain_2 + '_STD.png')
                #selected_df = std_df[std_df['Standard Deviation'] > std_df['Standard Deviation'].mean(axis = 0)]
                selected_df = std_df[std_df['Standard Deviation'] > std_df['Standard Deviation'].mean(axis = 0) + std_df['Standard Deviation'].std()]
                if(selected_df.empty):
                    selected_df = std_df[std_df['Standard Deviation'] > std_df['Standard Deviation'].mean(axis = 0)]
                selected_windows = selected_df['Window'].values

                pickle.dump(selected_windows, open("MutationAccumulations/TheMostImportantWindows/" + prot_1 + '_' + chain_1 + '_and_' + prot_2 + '_' + chain_2 + ".pickle" ,"wb"))
                #pickle.dump(selected_windows, open("MutationAccumulations/ImportantWindows/" + prot_1 + '_' + chain_1 + '_and_' + prot_2 + '_' + chain_2 + ".pickle" ,"wb"))

    #Use Mutation Windows (with high Standard Deviation of the frequency of mutations occur in different cancer types) to construct each patients' input features
    def extract_important_windows_binary_vectors(self, window_size):
        only_mutated = "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
        patients_properties = "PatientsSequencesInEachRefSeqChain/"
        where_windows = "TheMostImportantWindows/"
        if(where_windows == "TheMostImportantWindows/"):
            max_num_of_windows = 5
        else:
            max_num_of_windows = 7
        for patient in os.listdir("PatientsSequencesInEachRefSeqChain"):
            patient_vecs = []
            with open(self.dataset_adr + "PPIs-WithRefSeqs.txt","rb") as f:
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
                    original_seq = ""
                    #CHAIN-1    
                    address = only_mutated + prot_1 + "/" + chain_1 + "/mutatedTrimedRefSeq.seq"
                    original_seq = ""
                    with open(address, "rb") as seqf:
                        for l in seqf.readlines():
                            l = str(l,'utf-8')
                            original_seq = l.rstrip()
                    #CHAIN-2
                    address = only_mutated + prot_2 + "/" + chain_2 + "/mutatedTrimedRefSeq.seq"
                    with open(address, "rb") as seqf:
                        for l in seqf.readlines():
                            l = str(l,'utf-8')
                            original_seq += l.rstrip()
                    #Important windows with high std in the mutations frequency
                    windows_list = pickle.load(open("MutationAccumulations/" + where_windows + prot_1 + '_' + chain_1 + '_and_' + prot_2 + '_' + chain_2 + ".pickle" ,"rb"))        
                    patient_seq = ""
                    if(os.path.exists(patients_properties + patient + "/" + prot_1 + "/" + chain_1)):
                        address = patients_properties + patient + "/" + prot_1 + "/" + chain_1 + "/mutated_sequence.seq"
                        with open(address, "rb") as seqf:
                            for l in seqf.readlines():
                                l = str(l,'utf-8')
                                patient_seq = l.rstrip()
                    else:
                        address = only_mutated + prot_1 + "/" + chain_1 + "/mutatedTrimedRefSeq.seq"
                        with open(address, "rb") as seqf:
                            for l in seqf.readlines():
                                l = str(l,'utf-8')
                                patient_seq = l.rstrip()
                    if(os.path.exists(patients_properties + patient + "/" + prot_2 + "/" + chain_2)):
                        address = patients_properties + patient + "/" + prot_2 + "/" + chain_2 + "/mutated_sequence.seq"
                        with open(address, "rb") as seqf:
                            for l in seqf.readlines():
                                l = str(l,'utf-8')
                                patient_seq += l.rstrip()
                    else:
                        address = only_mutated + prot_2 + "/" + chain_2 + "/mutatedTrimedRefSeq.seq"
                        with open(address, "rb") as seqf:
                            for l in seqf.readlines():
                                l = str(l,'utf-8')
                                patient_seq += l.rstrip()
                    ppi_binary_vec = []
                    j = 0
                    for index in range(0, len(original_seq), window_size):
                        if(j in windows_list):
                            patient_ppi_seq = patient_seq[index : index + window_size]
                            wildtype_ppi_seq = original_seq[index : index + window_size]
                            for i in range(len(wildtype_ppi_seq)):
                                if(wildtype_ppi_seq[i] == patient_ppi_seq[i]):
                                    ppi_binary_vec.append(0)
                                else:
                                    ppi_binary_vec.append(1)
                            for i in range(len(wildtype_ppi_seq), window_size):
                                ppi_binary_vec.append(0)
                        j += 1
                    this_num_of_windows = len(ppi_binary_vec)
                    for i in range(this_num_of_windows, max_num_of_windows*window_size):
                        ppi_binary_vec.append(0)
                    patient_vecs.append(ppi_binary_vec)
            pickle.dump(patient_vecs, open(where_windows + patient + "_ppi.pickle","wb"))

    #Use Mutation Windows (with high Standard Deviation of the frequency of mutations occur in different cancer types) to construct each patients' input features
    #This version is for NLP-based Data Extraction
    def extract_important_windows_selfatt(self, window_size)
        zero_vec = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        only_mutated = "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
        patients_properties = "PatientsSequencesInEachRefSeqChain/"
        aa_w2v = pickle.load(open(self.dataset_adr + "AminoAcids_OneHotVectors.pickle", "rb"))
        where_windows = "TheMostImportantWindows/"
        max_num_of_windows = 5
        for patient in os.listdir("PatientsSequencesInEachRefSeqChain"):
            patient_vecs = []
            with open(self.dataset_adr + "PPIs-WithRefSeqs.txt","rb") as f:
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
                    original_seq = ""
                    #CHAIN-1    
                    address = only_mutated + prot_1 + "/" + chain_1 + "/mutatedTrimedRefSeq.seq"
                    original_seq = ""
                    with open(address, "rb") as seqf:
                        for l in seqf.readlines():
                            l = str(l,'utf-8')
                            original_seq = l.rstrip()
                    #CHAIN-2
                    address = only_mutated + prot_2 + "/" + chain_2 + "/mutatedTrimedRefSeq.seq"
                    with open(address, "rb") as seqf:
                        for l in seqf.readlines():
                            l = str(l,'utf-8')
                            original_seq += l.rstrip()
                    #Important windows with high std in the mutations frequency
                    windows_list = pickle.load(open("MutationAccumulations/" + where_windows + prot_1 + '_' + chain_1 + '_and_' + prot_2 + '_' + chain_2 + ".pickle" ,"rb"))        
                    patient_seq = ""
                    if(os.path.exists(patients_properties + patient + "/" + prot_1 + "/" + chain_1)):
                        address = patients_properties + patient + "/" + prot_1 + "/" + chain_1 + "/mutated_sequence.seq"
                        with open(address, "rb") as seqf:
                            for l in seqf.readlines():
                                l = str(l,'utf-8')
                                patient_seq = l.rstrip()
                    else:
                        address = only_mutated + prot_1 + "/" + chain_1 + "/mutatedTrimedRefSeq.seq"
                        with open(address, "rb") as seqf:
                            for l in seqf.readlines():
                                l = str(l,'utf-8')
                                patient_seq = l.rstrip()
                    if(os.path.exists(patients_properties + patient + "/" + prot_2 + "/" + chain_2)):
                        address = patients_properties + patient + "/" + prot_2 + "/" + chain_2 + "/mutated_sequence.seq"
                        with open(address, "rb") as seqf:
                            for l in seqf.readlines():
                                l = str(l,'utf-8')
                                patient_seq += l.rstrip()
                    else:
                        address = only_mutated + prot_2 + "/" + chain_2 + "/mutatedTrimedRefSeq.seq"
                        with open(address, "rb") as seqf:
                            for l in seqf.readlines():
                                l = str(l,'utf-8')
                                patient_seq += l.rstrip()
                    ppi_binary_vec = []
                    j = 0
                    for index in range(0, len(original_seq), window_size):
                        if(j in windows_list):
                            patient_ppi_seq = patient_seq[index : index + window_size]
                            wildtype_ppi_seq = original_seq[index : index + window_size]
                            for i in range(len(wildtype_ppi_seq)):
                                ppi_binary_vec.append(aa_w2v[patient_ppi_seq[i]])
                            for i in range(len(wildtype_ppi_seq), window_size):
                                ppi_binary_vec.append(zero_vec)
                        j += 1
                    this_num_of_windows = len(ppi_binary_vec)
                    for i in range(this_num_of_windows, max_num_of_windows*window_size):
                        ppi_binary_vec.append(zero_vec)
                    patient_vecs.append(ppi_binary_vec)
            pickle.dump(patient_vecs, open("NLPbased/" + where_windows + patient + "_ppi.pickle","wb"))



#p = TCGADataPreprocessing("DATASET/")
#p.extract_partner_proteins() 
#p.save_biophysical_PPIs("rcsb_p90.txt")
#p.keep_PPIs_based_on_bit_score_and_evalue("OncomineProteins.pickle","rcsb_p90.txt") 
#p.create_ATOM_files()    
#p.extract_start_end_positions_one_part()     
#p.extract_start_end_positions_more_parts() 
#p.organize_indices()
#p.save_RefSeq_positions()    
#p.create_rosetta_json_files()
#p.prepare_rosetta_files("TCGA_ANNOVAR/")  
#p.PDB_chain_starting_Rosetta_Num()
#p.extract_exact_positions()  
#p.extract_RefSeqID_sequences_from_NCBI()  
#p.create_mutfile_resfile_Rosetta() 
#p.extract_position_mutated_aminoacids("TCGA_ANNOVAR/")
#p.keep_aa_reported_in_ATOM_files()
#p.extract_aa_mutated()
#p.obtain_patients_sequences()
#p.extract_important_windows()
#p.extract_important_windows_binary_vectors(5)
#p.extract_important_windows_selfatt(5)
