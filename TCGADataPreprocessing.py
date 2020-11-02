#Find the proteins sequences involved in the Protein-Protein Interactions
all_proteins = []
with open("MutSigPPI/ppis.txt", "rb") as file:
    for linep in file.readlines():
        linep = str(linep,'utf-8')
        prots = linep.split(",")
        prot1 = prots[0].split()[0]
        prot2 = prots[1].split()[0]
        if(prot1 not in all_proteins):
            all_proteins.append(prot1)
        if(prot2 not in all_proteins):
            all_proteins.append(prot2)
prot_seqs_adr = "MutSigPPI/Proteins_sequences/"
start_to_write = False
with open("MutSigPPI/uniprot-proteomes.fasta", "rb") as file:
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


#------------------------------------------------------------------------------------------------------------------------------------------
#Obtain the proteins dictionary with their corresponding interaction ids 
#The interactions must be have Identities >= 90
import pickle

proteins_interactions_id = {}
current_protein = ""
current_interaction = ""
with open("MutSigPPI/rcsb_p70.txt", "rb") as file:
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
pickle.dump(proteins_interactions_id, open("MutSigPPI/protein_physical_interactions_70percent.pickle","wb"))


#------------------------------------------------------------------------------------------------------------------------------------------
#Moreover, we need to remove those interactions with bit_score < 50 and e_value > 0.000001

import pickle

def get_set(l1,l2):
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
    
proteins_interactions_id = pickle.load(open("MutSigPPI/protein_physical_interactions_90percent.pickle","rb"))
oncomine_prots = pickle.load(open("MutSigPPI/OncomineProteins.pickle","rb"))
start_to_find_interactions = False
current_protein = ""
with open("MutSigPPI/rcsb_p90.txt", "rb") as file:
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
interaction_file = open("MutSigPPI/BLAST_protein_interactions_90.txt","a+")
num_of_physical_interactions = 0
proteins = list(proteins_interactions_id.keys())
for prot1 in range(len(proteins)):
    for prot2 in range(prot1 + 1, len(proteins)):
        any_common = get_set(proteins_interactions_id[proteins[prot1]],proteins_interactions_id[proteins[prot2]])
        #if( (proteins[prot1] in oncomine_prots) or (proteins[prot1] in oncomine_prots) ):
        if(len(any_common) > 0):
            num_of_physical_interactions += 1
            interaction_file.write(proteins[prot1] + "\t" + proteins[prot2] + "\t")
            for pdbid in any_common:
                interaction_file.write(pdbid + "\t")
            interaction_file.write("\n")    
interaction_file.close()
print(num_of_physical_interactions)


#------------------------------------------------------------------------------------------------------------------------------------------
#Create ATOM files and Obtain the Chain sequences for Protein-Protein Interactions
import requests
import os

with open("MutSigPPI/ppis.txt", "rb") as file:
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
        save_to = "MutSigPPI/ATOM_files/"
        pdb_seq_adr = "MutSigPPI/PDB_sequences/"
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


#------------------------------------------------------------------------------------------------------------------------------------------
#Create pickle of proteins and their corresponding RefSeq IDs
import pickle

main_folder = "MutSigPPI/"
proteins_refseqs_ids = {}
with open(main_folder + "ProteinsFullInformation.txt", "rb") as pfi:
    for linepinf in pfi.readlines():
        linepinf = str(linepinf,'utf-8')
        linepinf = linepinf.rstrip()
        allinfosp = linepinf.split()
        proteins_refseqs_ids[allinfosp[0]] = []
        if(len(allinfosp) == 3):
            proteins_refseqs_ids[allinfosp[0]].append(allinfosp[2])
        elif(len(allinfosp) > 3):
            for refseq_index in range(2, len(allinfosp)):
                refseq_id = allinfosp[refseq_index].split(",")[0]
                proteins_refseqs_ids[allinfosp[0]].append(refseq_id)
pickle.dump(proteins_refseqs_ids, open(main_folder + "proteins_refseqs_ids.pickle","wb"))


#------------------------------------------------------------------------------------------------------------------------------------------
#Find the longest common substring between two strings 
def lcs(S,T):
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


#------------------------------------------------------------------------------------------------------------------------------------------
#Extract the start and end positions for those proteins which their corresponding PDB and chain sequence has only one part
#Something like this: ----XXXXXXX-
import os
import pickle

main_folder = "MutSigPPI/"
proteins_refseqs_ids = pickle.load(open(main_folder + "proteins_refseqs_ids.pickle","rb"))
o = 0
indicesfile = open(main_folder + "indices.txt","a+")
with open(main_folder + "ppis.txt", "rb") as file:
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
                    with open(main_folder + "RefSeqSequences/" + refseq_id + ".seq", "rb") as seqfile:
                        for lineseq in seqfile.readlines():
                            lineseq = str(lineseq,'utf-8')
                            prot_seq = prot_seq + lineseq.rstrip()
                    prot_seqs[refseq_id] = prot_seq
            else:
                prot_seq = ""
                with open(main_folder + "RefSeqSequences/" + prots[rs] + ".seq", "rb") as seqfile:
                    for lineseq in seqfile.readlines():
                        lineseq = str(lineseq,'utf-8')
                        prot_seq = prot_seq + lineseq.rstrip()
                prot_seqs[prots[rs]] = prot_seq
            chain_seq = ""
            with open(main_folder + "Available_AminoAcids_in_seq/" + pdb_id + "_" + chains[rs] + ".seq", "rb") as availfile:
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
                    lcs1 = list(lcs(a,b))[0]
                    #Real index for starting_index=1 instead of 0
                    start_index = a.find(lcs1) + 1
                    end_index = start_index + len(lcs1) - 1
                    indicesfile.write(refseq_id + "\t" + pdb_id + "_" + chains[rs] + "\t" + str(start_index) + "\t" + str(end_index) + "\n")
indicesfile.close()                


#------------------------------------------------------------------------------------------------------------------------------------------
#Extract the start and end positions (more than one <start,end>) for those proteins which their corresponding PDB and chain sequence has more than one part
#Something like this: ----XXXXXXX---XXX-XXXX---   
import os
main_folder = "MutSigPPI/"
o = 0
indicesfile = open(main_folder + "indices_multipart.txt","a+")
with open(main_folder + "ppis.txt", "rb") as file:
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
                    with open(main_folder + "RefSeqSequences/" + refseq_id + ".seq", "rb") as seqfile:
                        for lineseq in seqfile.readlines():
                            lineseq = str(lineseq,'utf-8')
                            prot_seq = prot_seq + lineseq.rstrip()
                    prot_seqs[refseq_id] = prot_seq
            else:
                prot_seq = ""
                with open(main_folder + "RefSeqSequences/" + prots[rs] + ".seq", "rb") as seqfile:
                    for lineseq in seqfile.readlines():
                        lineseq = str(lineseq,'utf-8')
                        prot_seq = prot_seq + lineseq.rstrip()
                prot_seqs[prots[rs]] = prot_seq
            chain_seq = ""
            with open(main_folder + "Available_AminoAcids_in_seq/" + pdb_id + "_" + chains[rs] + ".seq", "rb") as availfile:
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


#------------------------------------------------------------------------------------------------------------------------------------------
# Organize the indices (Show the exact position of start and ends)
main_folder = "MutSigPPI/"
indicesfile = open(main_folder + "indices_multipart_organized.txt","a+")
with open(main_folder + "indices_multipart.txt", "rb") as file:
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


#------------------------------------------------------------------------------------------------------------------------------------------
#Pickle the RefSeq positions
import pickle

main_folder = "MutSigPPI/"
protein_positions = {}
with open(main_folder + "indices.txt", "rb") as file:
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
pickle.dump(protein_positions, open(main_folder + "RefSeqs_positions.pickle","wb"))


#------------------------------------------------------------------------------------------------------------------------------------------
#Create pdb2rosetta and rosetta2pdb json files (for flex ddG)
import os

main_folder = "MutSigPPI/"
with open(main_folder + "ppis.txt", "rb") as fm:
    lines = fm.readlines()
    for linem in lines:
        linem = linem.rstrip()
        linem = str(linem,'utf-8')
        pinfos = linem.split()
        chain1 = pinfos[1].split(",")[0].split("_")[1]
        pdb_chain = pinfos[3].split("_")
        pdb_name = pdb_chain[0]
        chain2 = pdb_chain[1]
        pdb_chains_adr = main_folder + "ATOM_files/" + pdb_name+"_"+chain1+chain2 + "/"
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
                        #chain_id = chainsline[21]
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


#------------------------------------------------------------------------------------------------------------------------------------------
#Create pdb2rosetta_complete json files (NOT for flex ddG, just for myself) It gives me the original amino acids
import os

amino_acids_mapper = {   'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                         'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
main_folder = "MutSigPPI/"
with open(main_folder + "ppis.txt", "rb") as fm:
    lines = fm.readlines()
    for linem in lines:
        linem = linem.rstrip()
        linem = str(linem,'utf-8')
        pinfos = linem.split()
        chain1 = pinfos[1].split(",")[0].split("_")[1]
        pdb_chain = pinfos[3].split("_")
        pdb_name = pdb_chain[0]
        chain2 = pdb_chain[1]
        pdb_chains_adr = main_folder + "ATOM_files/" + pdb_name+"_"+chain1+chain2 + "/"
        if(not os.path.exists(pdb_chains_adr + pdb_name+"_"+chain1+chain2+".pdb")):
            print("ERROR: The following ATOM file does not exist:")
            print(pdb_name+"_"+chain1+chain2)
            print("Therefore, we ignore it!")
        else:
            if(not os.path.exists(pdb_chains_adr + "pdb2rosetta_complete.json")):
                print(pdb_name+"_"+chain1+chain2)
                with open(pdb_chains_adr + pdb_name+"_"+chain1+chain2+".pdb","rb") as fp:
                    rosetta2pdb_file = open(pdb_chains_adr + "pdb2rosetta_complete.json","a+")
                    alllines = fp.readlines()
                    previous_num = "DifferentThing"
                    number_to_assign = 1
                    rosetta2pdb_content = "{\n"
                    for chainsline in alllines:
                        chainsline = chainsline.rstrip()
                        chainsline = str(chainsline,'utf-8')
                        aminoacid = chainsline[17:20]
                        #chain_id = chainsline[21]
                        current_num = chainsline[22:27].strip()
                        if(current_num != previous_num):
                            rosetta2pdb_content += "\"" + str(number_to_assign) + "\": " + "\"" + chainsline[21:27] + " " + amino_acids_mapper[aminoacid] + "\",\n"
                            number_to_assign += 1
                        previous_num = current_num

                    rosetta2pdb_file.write(rosetta2pdb_content[:len(rosetta2pdb_content)-2] + "\n}")
                    rosetta2pdb_file.close()


#------------------------------------------------------------------------------------------------------------------------------------------
#Prepare 'resfile' and 'mutfile' files in folder 'PatientsPPImutationsData'
#It is for extracting the mutations and the exactly position of the mutations in two different standard positions
#1. Rosetta 2. pdb ATOM position
import os
import pickle
import json

main_folder = "MutSigPPI/"
target_folder = main_folder + "PatientsPPImutationsData/"
chains_start_numbers_rosetta = pickle.load(open(main_folder + "chains_start_numbers_rosetta.pickle", "rb"))
exact_positions_protein_refseq = pickle.load(open(main_folder + "exact_positions_protein_refseq.pickle","rb"))
dirr = "TCGA_ANNOVAR/"
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
                                                    with open(main_folder + "RefSeqSequences/" + current_refseq + ".seq", "rb") as refseqf:
                                                        liners = refseqf.readlines()
                                                        for aline in liners:
                                                            aline = aline.rstrip()
                                                            aline = str(aline,'utf-8')
                                                            NM_seq += aline
                                                    chain_part_seq = ""
                                                    with open(main_folder + "Available_AminoAcids_in_seq/" + PDBchain + ".seq", "rb") as pcf:
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
                                                                    main_folder = "MutSigPPI-01112019/"
                                                                    with open(main_folder + 'ATOM_files/' + pdb_id + "_" + chains + '/pdb2rosetta_complete.json') as json_file:
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
                                                with open(main_folder + "RefSeqSequences/" + prot + ".seq", "rb") as refseqf:
                                                    liners = refseqf.readlines()
                                                    for aline in liners:
                                                        aline = aline.rstrip()
                                                        aline = str(aline,'utf-8')
                                                        NM_seq += aline

                                                chain_part_seq = ""
                                                with open(main_folder + "Available_AminoAcids_in_seq/" + PDBchain + ".seq", "rb") as pcf:
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
                                                                main_folder = "MutSigPPI-01112019/"
                                                                with open(main_folder + 'ATOM_files/' + pdb_id + "_" + chains + '/pdb2rosetta_complete.json') as json_file:
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


#------------------------------------------------------------------------------------------------------------------------------------------
#Each chain in an PDB complex has a starting Rosetta number
import os
import pickle

chains_start_numbers_rosetta = {}
main_folder = "MutSigPPI/ATOM_files/"
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
pickle.dump(chains_start_numbers_rosetta, open("MutSigPPI-01112019/chains_start_numbers_rosetta.pickle","wb"))


#------------------------------------------------------------------------------------------------------------------------------------------
import pickle

main_folder = "MutSigPPI/"
protein_refseq_positions = pickle.load(open(main_folder + "RefSeqs_positions.pickle","rb"))
refseqs_of_prots = pickle.load(open(main_folder + "proteins_refseqs_ids.pickle","rb"))
protein_refseq_positions['NM_000612']


#------------------------------------------------------------------------------------------------------------------------------------------
#Protein_name -> RefSeq_ID -> Positions
import pickle

main_folder = "MutSigPPI/"
avail_prots = []
with open(main_folder + "ppis.txt", "rb") as file:
    for linep in file.readlines():
        linep = str(linep,'utf-8')
        first_prot = linep.split()[0]
        second_prot = linep.split()[2]
        if(first_prot not in avail_prots):
            avail_prots.append(first_prot)
        if(second_prot not in avail_prots):
            avail_prots.append(second_prot)
bb = pickle.load(open(main_folder + "proteins_refseqs_ids.pickle","rb"))
protein_refseq_positions = pickle.load(open(main_folder + "RefSeqs_positions.pickle","rb"))
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
pickle.dump(all_protein_refseqs_positions,open(main_folder + "exact_positions_protein_refseq.pickle","wb"))


#------------------------------------------------------------------------------------------------------------------------------------------
#Extract RefSeqGene ID protein (amino acid) sequences from NCBI website
from selenium import webdriver

main_folder = "MutSigPPI/"
refseq_adr = main_folder + "RefSeqSequences/"
all_without_seq = main_folder + "all_refseqs_without_seq.txt"
fall = open(all_without_seq,"w+")
with open(main_folder + "ProteinsFullInformation.txt","rb") as allrefseqfile:
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


#------------------------------------------------------------------------------------------------------------------------------------------
#The proteins which we have no RefSeqGene ID for them
with open(main_folder + "ProteinsFullInformation.txt","rb") as allrefseqfile:
    refseq_lines = allrefseqfile.readlines()
    for refseq_line in refseq_lines:
        refseq_line = str(refseq_line,'utf-8')
        refseq_line = refseq_line.rstrip()
        allpinfos = refseq_line.split()
        if(len(allpinfos) < 3):
            print(allpinfos[0])


#------------------------------------------------------------------------------------------------------------------------------------------
#To report the mutations occurred in each PDB complex
#Create files for mutations.resfile, mutations.mutfile, nataa_mutations.resfile, and chains_to_move.txt for each PDB complex
import os
import pickle
import json

main_folder = "MutSigPPI/"
target_folder = main_folder + "PatientsPPImutationsData/"
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


#------------------------------------------------------------------------------------------------------------------------------------------
#Extract the positions of all amino acids which are mutated in all patients  (These position are the positions in the available amino acids in PDB file)
#The positions in "mut_positions_RefSeq" are sorted!
#We will keep only these position for our analysis process
#For example: original_AA_sequence => ABCDEFGHIJKLMN and mut_positions_RefSeq => [1,3,6,8,9]
#The result will be => ACFHI 
import os
import pickle
import json

main_folder = "MutSigPPI/"
exact_positions_protein_refseq = pickle.load(open(main_folder + "exact_positions_protein_refseq.pickle","rb"))
mut_positions_RefSeq = {}
dirr = "TCGA_ANNOVAR/"
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


#------------------------------------------------------------------------------------------------------------------------------------------
#Trim the sequences (Keep those Amino Acids which their ATOMs has reported in PDB)
import pickle
import os

main_folder = "MutSigPPI-01112019/"
exact_positions_protein_refseq = pickle.load(open(main_folder + "exact_positions_protein_refseq.pickle","rb"))
destination = "RefSeq-TheSeqsContributeInPPI/"
for protein in exact_positions_protein_refseq:
    the_prot_allRefSeqs = exact_positions_protein_refseq[protein]
    for current_refseq in the_prot_allRefSeqs:
        NM_seq = ""
        adr_folder_refseq = ""
        if(current_refseq != "NOavailableRefSeq"):
            adr_folder_refseq = main_folder + destination + current_refseq
            if(not os.path.exists(adr_folder_refseq)):
                os.mkdir(adr_folder_refseq)
            with open(main_folder + "RefSeqSequences/" + current_refseq + ".seq", "rb") as refseqf:
                liners = refseqf.readlines()
                for aline in liners:
                    aline = aline.rstrip()
                    aline = str(aline,'utf-8')
                    NM_seq += aline
        else:
            adr_folder_refseq = main_folder + destination + protein
            if(not os.path.exists(adr_folder_refseq)):
                os.mkdir(adr_folder_refseq)
            with open(main_folder + "RefSeqSequences/" + protein + ".seq", "rb") as refseqf:
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


#------------------------------------------------------------------------------------------------------------------------------------------
#We extract amino acids in the specific positions reported in "mut_positions_RefSeq"
#"mut_positions_RefSeq" only keeps the positions of the amino acids which mutated in at least one patient
import pickle
import os

main_folder = "MutSigPPI-01112019/"
refseqs_seq_folder = main_folder + "RefSeq-TheSeqsContributeInPPI/"
refseqs_only_mutated_aa_folder = main_folder + "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
mut_positions_RefSeq = pickle.load(open(main_folder + "Mutation_Positions_In_Trimed_RefSeqs.pickle","rb"))
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

#------------------------------------------------------------------------------------------------------------------------------------------
#Obtain Each Patient final sequence after mutations (just mutated sequences) to feed as a sequence to our model
#If there is no mutation occurred in an PPI, then we do not save its sequence in "PatientsSequencesInEachRefSeqChain" folder
#In the case that there is no mutation in an PPI, you must use its sequence in "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients" folder
import pickle
import os

main_folder = "MutSigPPI-01112019/"
target_folder = main_folder + "PatientsSequencesInEachRefSeqChain/"
refseqs_only_mutated_aa_folder = main_folder + "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
forEachPatient_mut_positions_RefSeq = pickle.load(open(main_folder + "Mutation_Positions_In_Trimed_RefSeqs_ForEachPatient.pickle","rb"))
mut_positions_RefSeq = pickle.load(open(main_folder + "Mutation_Positions_In_Trimed_RefSeqs.pickle","rb"))
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

#------------------------------------------------------------------------------------------------------------------------------------------
#Extract the inputs/outputs of the Self-Attention 
import torch
import pickle
aa_numbers = pickle.load(open("MutSigPPI/DATASET/patientsAAnumbers.pickle","rb"))
for patient in aa_numbers:
    inputs = []
    outputs = []
    for ppi in range(len(aa_numbers[patient])):
        for an_input in aa_numbers[patient][ppi]["first5aa"]:
            for an_output in aa_numbers[patient][ppi]["second5aa"]:
                inputs.append(an_input)
                outputs.append(an_output)
                inputs.append(an_output)
                outputs.append(an_input)
    inputs = torch.tensor(inputs)
    outputs = torch.tensor(outputs)
    outputs = outputs.flatten()
    pickle.dump(inputs, open("MutSigPPI/DATASET/patientsAAfull/" + patient + "-inputs.pickle", "wb"))
    pickle.dump(outputs, open("MutSigPPI/DATASET/patientsAAfull/" + patient + "-outputs.pickle", "wb"))


#------------------------------------------------------------------------------------------------------------------------------------------
#Extract Mutational Binary Vectors for each patient
import os
import numpy as np
import pickle

main_folder = "MutSigPPI/"
patients_properties = main_folder + "PatientsSequencesInEachRefSeqChain/"
only_mutated = main_folder + "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
for patient in os.listdir(patients_properties):
    final_binary_vec = np.array([])
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
        #CHAIN-1    
            address = only_mutated + prot_1 + "/" + chain_1 + "/mutatedTrimedRefSeq.seq"
            original_seq = ""
            with open(address, "rb") as seqf:
                for l in seqf.readlines():
                    l = str(l,'utf-8')
                    original_seq =l.rstrip()
            mutation_binary_vec1 = np.zeros(len(original_seq))
            if(os.path.exists(patients_properties + patient + "/" + prot_1 + "/" + chain_1)):
                address = patients_properties + patient + "/" + prot_1 + "/" + chain_1 + "/mutated_sequence.seq"
                mutated_seq = ""
                with open(address, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        mutated_seq = l.rstrip()
                for i in range(len(original_seq)):
                    if(mutated_seq[i] != original_seq[i]):
                        mutation_binary_vec1[i] = 1
        #CHAIN-2
            address = only_mutated + prot_2 + "/" + chain_2 + "/mutatedTrimedRefSeq.seq"
            original_seq = ""
            with open(address, "rb") as seqf:
                for l in seqf.readlines():
                    l = str(l,'utf-8')
                    original_seq =l.rstrip()
            mutation_binary_vec2 = np.zeros(len(original_seq))
            if(os.path.exists(patients_properties + patient + "/" + prot_2 + "/" + chain_2)):
                address = patients_properties + patient + "/" + prot_2 + "/" + chain_2 + "/mutated_sequence.seq"
                mutated_seq = ""
                with open(address, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        mutated_seq = l.rstrip()
                for i in range(len(original_seq)):
                    if(mutated_seq[i] != original_seq[i]):
                        mutation_binary_vec2[i] = 1
            final_binary_vec = np.concatenate([final_binary_vec, mutation_binary_vec1, mutation_binary_vec2])
    pickle.dump(final_binary_vec, open(main_folder + "PatientsMutBinaryVectors/" + patient + ".pickle","wb"))


#------------------------------------------------------------------------------------------------------------------------------------------
#Extract Mutation Windows with high Standard Deviation of the frequency of mutations occur in different cancer types
import numpy as np
import pickle
import os
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt

window_size = 5
only_mutated = main_folder + "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
patients_properties = main_folder + "PatientsSequencesInEachRefSeqChain/"
main_folder = "MutSigPPI/"
tcga_clinical_dataframe = pickle.load(open(main_folder + "TCGA_clinical_dataframe.pickle","rb"))
classes = {"adrenal gland":0, "bladder":1, "breast":2, "GYN":3, "bile duct":4, "CRC":5, "bone marrow":6, "esophagus":7, "brain":8, "head and neck":9, "Kidney":10, "liver":11, "Lung":12, "pleura":13, "pancreatic":14, "male reproductive system":15, "other":16, "Melanoma":17, "stomach":18, "thyroid":19, "thymus":20}
classes_map = {"adrenal gland": "AdrGl", "bladder":"Blad", "breast":"brst", "GYN":"GYN", "bile duct":"BileD", "CRC":"CRC", "bone marrow":"BoneM", "esophagus":"Esoph", "brain":"Brain", "head and neck":"H&N", "Kidney":"Kidney", "liver":"Liver", "Lung":"Lung", "pleura":"Pleu", "pancreatic":"Panc", "male reproductive system":"MaleRep", "other":"Other", "Melanoma":"Melan", "stomach":"Stomc", "thyroid":"Thyr", "thymus":"Thym"}
patients = {}
for cancer_class in classes:
    patients[cancer_class] = tcga_clinical_dataframe[tcga_clinical_dataframe['cancer_class'] == cancer_class].index
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


#------------------------------------------------------------------------------------------------------------------------------------------
#Use Mutation Windows (with high Standard Deviation of the frequency of mutations occur in different cancer types) to construct each patients' input features
import pickle
import os

window_size = 5
only_mutated = main_folder + "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
patients_properties = main_folder + "PatientsSequencesInEachRefSeqChain/"
main_folder = "MutSigPPI/"
where_windows = "TheMostImportantWindows/"
if(where_windows == "TheMostImportantWindows/"):
    max_num_of_windows = 5
else:
    max_num_of_windows = 7
for patient in os.listdir(main_folder + "PatientsSequencesInEachRefSeqChain"):
    patient_vecs = []
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
    pickle.dump(patient_vecs, open(main_folder + "DATASET/" + where_windows + patient + "_ppi.pickle","wb"))


#------------------------------------------------------------------------------------------------------------------------------------------
#Use Mutation Windows (with high Standard Deviation of the frequency of mutations occur in different cancer types) to construct each patients' input features
#This version is for NLP-based Data Extraction
import pickle
import os

main_folder = "MutSigPPI/"
window_size = 5
zero_vec = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
only_mutated = main_folder + "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
patients_properties = main_folder + "PatientsSequencesInEachRefSeqChain/"
aa_w2v = pickle.load(open(main_folder + "AminoAcids_OneHotVectors.pickle", "rb"))
where_windows = "TheMostImportantWindows/"
max_num_of_windows = 5
for patient in os.listdir(main_folder + "PatientsSequencesInEachRefSeqChain"):
    patient_vecs = []
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
            windows_list = pickle.load(open(main_folder + "MutationAccumulations/" + where_windows + prot_1 + '_' + chain_1 + '_and_' + prot_2 + '_' + chain_2 + ".pickle" ,"rb"))        
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
    pickle.dump(patient_vecs, open(main_folder + "DATASET/NLPbased/" + where_windows + patient + "_ppi.pickle","wb"))


#------------------------------------------------------------------------------------------------------------------------------------------
#To evaluate our model PPI predictions, we used PPI Essentiality extracted from MEDICI paper
import pickle
import pandas as pd
import math
import seaborn as sns
sns.set(style="whitegrid")

column_headers = ["Cancer Organ", "nDCG"]
values = {}
values["Cancer Organ"] = []
values["nDCG"] = []
column_headers2 = ["Cancer Organ", "quantity", "evaluation"]
values2 = {}
values2["Cancer Organ"] = []
values2["quantity"] = []
values2["evaluation"] = []
a = pd.read_excel('MutSigPPI/ESSENTIALITY-PPI/Essentiality.xlsx')
aa = pickle.load(open('MutSigPPI/available_ppi_in_essentiality.pickle','rb'))
how_many_cols = {'male reproductive system': 3, 'liver': 1, 'esophagus': 8, 'GYN': 23, 'Melanoma': 7, 'pancreatic': 12, 'Lung': 10, 'Kidney': 2, 'stomach': 4, 'brain': 16, 'breast': 11, 'CRC': 17, 'bladder': 1, 'bone marrow': 4}
for cancer_type in how_many_cols:
    cancer_essentialities = pickle.load(open("MutSigPPI/ESSENTIALITY-PPI/" + "ESSENTIALITY-" + cancer_type + ".pickle", "rb"))
    which_cols = []
    which_cols.append(cancer_type)
    for k in range(1,how_many_cols[cancer_type]):
        which_cols.append(cancer_type + "_" + str(k))
    i = 1
    DCG = 0
    with open("MutSigPPI/RESULTS/RelatedPPIs/" + cancer_type + "-PPI-Contributions.txt", "rb") as f:
        for line in f.readlines():
            line = str(line,'utf-8')
            line = line.rstrip()
            line = line.split()[0]
            pr1 = line.split("-")[0]
            pr2 = line.split("-")[1]
            if((pr1 + "-" + pr2) in cancer_essentialities):
                avg = cancer_essentialities[pr1 + "-" + pr2]
                DCG += avg/math.log2(i+1)
                i += 1
            elif((pr2 + "-" + pr1) in cancer_essentialities):
                avg = cancer_essentialities[pr2 + "-" + pr1]
                DCG += avg/math.log2(i+1)
                i += 1
    sorted_l = sorted(cancer_essentialities, key=cancer_essentialities.get, reverse=True)
    i = 1
    IDCG = 0
    for ppi in sorted_l:
        avg = cancer_essentialities[ppi]
        IDCG += avg/math.log2(i+1)
        i += 1
    nDCG = DCG/IDCG
    print("Cancer Type= " + cancer_type + " DCG= " + str(DCG) + " IDCG= " + str(IDCG) + " NDCG= " + str(nDCG))
    values["Cancer Organ"].append(cancer_type)
    values["nDCG"].append(nDCG)
    values2["Cancer Organ"].append(cancer_type)
    values2["evaluation"].append('DCG')
    values2["quantity"].append(DCG)
    values2["Cancer Organ"].append(cancer_type)
    values2["evaluation"].append('IDCG')
    values2["quantity"].append(IDCG)
df = pd.DataFrame(values, columns = column_headers)
ax = sns.barplot(x="Cancer Organ", y="nDCG", data=df)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_ylim([0,1])
figure = ax.get_figure()    
figure.savefig('MutSigPPI/RESULTS/nDCG.png', dpi=350, bbox_inches='tight')
df2 = pd.DataFrame(values2, columns = column_headers2)
ax = sns.barplot(x="Cancer Organ", y="quantity", hue='evaluation', palette="Set2", data=df2)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_ylim([0,7])
figure = ax.get_figure()    
figure.savefig('MutSigPPI/RESULTS/DCG-IDCG.png', dpi=350, bbox_inches='tight')

#------------------------------------------------------------------------------------------------------------------------------------------
#Draw Survival Plot (Kaplan-Meier)
import pandas as pd
import pickle
import os
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt

patients_of_cancerTypes = pickle.load(open("MutSigPPI/CancerTypePatients.pickle","rb"))
proteins_refseq_ids = pickle.load(open("MutSigPPI/proteins_refseqs_ids.pickle", "rb"))
clinical_info = pickle.load(open("MutSigPPI/TCGA_clinical_dataframe.pickle","rb"))
ppis_in_order = pickle.load(open("MutSigPPI/PPI_Chains_In_Order.pickle","rb"))
clinical_info = clinical_info[['days_to_last_followup','vital_status','days_to_death','lost_follow_up']]
clinical_info = clinical_info.replace({'vital_status': 'Dead'}, 1)
clinical_info = clinical_info.replace({'vital_status': 'Alive'}, 0)
cancer_class = "GYN"
for ppi_chain in ppis_in_order:
    if(ppi_chain.find(")-") == -1):
        splitted_ppi = ppi_chain.split(";")
        pr1 = splitted_ppi[0]
        pr2 = splitted_ppi[1].split("(")[0]
        refseqs1 = proteins_refseq_ids[pr1]
        refseqs2 = proteins_refseq_ids[pr2]
        mutated_patients = []
        non_mutated_pateints = []
        ppi_refseqs = []
        if(len(refseqs1) == 0):
            ppi_refseqs.append(pr1)
        else:
            for refseq in refseqs1:
                ppi_refseqs.append(refseq)
        if(len(refseqs2) == 0):
            ppi_refseqs.append(pr2)
        else:
            for refseq in refseqs2:
                ppi_refseqs.append(refseq)
        for patient in patients_of_cancerTypes[cancer_class]:
            at_least_one_mutation = False
            for refseq in ppi_refseqs:
                if(os.path.exists("MutSigPPI/PatientsSequencesInEachRefSeqChain/" + patient + "/" + refseq)):
                    at_least_one_mutation = True
                    break
            if(at_least_one_mutation):
                mutated_patients.append(patient)
            else:
                non_mutated_pateints.append(patient)
        mutated_df = clinical_info.loc[mutated_patients]
        mutated_df = mutated_df[(mutated_df['days_to_last_followup'] != '-') | (mutated_df['days_to_death'] != '-')]
        n_df1 = mutated_df[(mutated_df['days_to_last_followup'] == '-') & (mutated_df['days_to_death'] != '-')][['vital_status','days_to_death']].rename(columns={"days_to_death": "Time"})
        n_df2 = mutated_df[mutated_df['days_to_last_followup'] != '-'][['vital_status','days_to_last_followup']].rename(columns={"days_to_last_followup": "Time"})
        mutated_df = n_df1.append(n_df2)
        mutated_df = mutated_df[mutated_df['vital_status'] != '-']

        nonmutated_df = clinical_info.loc[non_mutated_pateints]
        nonmutated_df = nonmutated_df[(nonmutated_df['days_to_last_followup'] != '-') | (nonmutated_df['days_to_death'] != '-')]
        n_df1 = nonmutated_df[(nonmutated_df['days_to_last_followup'] == '-') & (nonmutated_df['days_to_death'] != '-')][['vital_status','days_to_death']].rename(columns={"days_to_death": "Time"})
        n_df2 = nonmutated_df[nonmutated_df['days_to_last_followup'] != '-'][['vital_status','days_to_last_followup']].rename(columns={"days_to_last_followup": "Time"})
        nonmutated_df = n_df1.append(n_df2)
        nonmutated_df = nonmutated_df[nonmutated_df['vital_status'] != '-']
        if((not mutated_df.empty) and (not nonmutated_df.empty)):
            if((len(mutated_df) > 10) and (len(nonmutated_df) > 10)):
                kmf = KaplanMeierFitter()
                kmf.fit(mutated_df["Time"], event_observed=mutated_df["vital_status"], label='PPI Lost')
                ax = kmf.plot(show_censors=True, ci_show=False)
                kmf.fit(nonmutated_df["Time"], event_observed=nonmutated_df["vital_status"], label='PPI Intact')
                kmf.plot(ax=ax,show_censors=True, ci_show=False)
                ax.get_figure().savefig("MutSigPPI/RESULTS/Survival/" + cancer_class + "_" + pr1 + "-" + pr2 + ".png", dpi=500)
                plt.close()


#------------------------------------------------------------------------------------------------------------------------------------------
#Obtain the Cox Proportional Hazard Model results
#Coef represents the Hazard Ratio
import pandas as pd
import pickle
import os
from lifelines.datasets import load_rossi
from lifelines import CoxPHFitter
import matplotlib.pyplot as plt
import numpy as np

main_folder = "MutSigPPI/"
cancer_class = "bladder"
patients_of_cancerTypes = pickle.load(open("MutSigPPI/CancerTypePatients.pickle","rb"))
proteins_refseq_ids = pickle.load(open("MutSigPPI/proteins_refseqs_ids.pickle", "rb"))
clinical_info = pickle.load(open("MutSigPPI/TCGA_clinical_dataframe.pickle","rb"))
clinical_info = clinical_info[['days_to_last_followup','vital_status','days_to_death','lost_follow_up']]
clinical_info = clinical_info.replace({'vital_status': 'Dead'}, 1)
clinical_info = clinical_info.replace({'vital_status': 'Alive'}, 0)
final_sruvival_df = ""
is_first_time = True
ll = []
for sp in os.listdir(main_folder + "RESULTS/Survival/" + cancer_class):
    sp = sp.split(".")[0].split("_")[1]
    if(sp not in ll):
        sps = sp.split("-")
        pr1 = ""
        pr2 = ""
        if(len(sps) == 3):
            pr1 = sps[0] + "-" + sps[1]
            pr2 = sps[2]
        elif(len(sps) == 4):
            pr1 = sps[0] + "-" + sps[1]
            pr2 = sps[2] + "-" + sps[3]
        else:
            pr1 = sps[0]
            pr2 = sps[1]
        refseqs1 = proteins_refseq_ids[pr1]
        refseqs2 = proteins_refseq_ids[pr2]
        mutated_patients = []
        non_mutated_pateints = []
        ppi_refseqs = []
        if(len(refseqs1) == 0):
            ppi_refseqs.append(pr1)
        else:
            for refseq in refseqs1:
                ppi_refseqs.append(refseq)
        if(len(refseqs2) == 0):
            ppi_refseqs.append(pr2)
        else:
            for refseq in refseqs2:
                ppi_refseqs.append(refseq)
        for patient in patients_of_cancerTypes[cancer_class]:
            at_least_one_mutation = False
            for refseq in ppi_refseqs:
                if(os.path.exists("MutSigPPI/PatientsSequencesInEachRefSeqChain/" + patient + "/" + refseq)):
                    at_least_one_mutation = True
                    break
            if(at_least_one_mutation):
                mutated_patients.append(patient)
            else:
                non_mutated_pateints.append(patient)
        mutated_df = clinical_info.loc[mutated_patients]
        mutated_df = mutated_df[(mutated_df['days_to_last_followup'] != '-') | (mutated_df['days_to_death'] != '-')]
        n_df1 = mutated_df[(mutated_df['days_to_last_followup'] == '-') & (mutated_df['days_to_death'] != '-')][['vital_status','days_to_death']].rename(columns={"days_to_death": "Time"})
        n_df2 = mutated_df[mutated_df['days_to_last_followup'] != '-'][['vital_status','days_to_last_followup']].rename(columns={"days_to_last_followup": "Time"})
        mutated_df = n_df1.append(n_df2)
        mutated_df = mutated_df[mutated_df['vital_status'] != '-']
        mutated_df[pr1 + "-" + pr2] = np.zeros(len(mutated_df))
        if(is_first_time):
            is_first_time = False
            final_sruvival_df = mutated_df
        else:
            final_sruvival_df = pd.concat([final_sruvival_df, mutated_df], axis=0, sort=False).fillna(0)
        nonmutated_df = clinical_info.loc[non_mutated_pateints]
        nonmutated_df = nonmutated_df[(nonmutated_df['days_to_last_followup'] != '-') | (nonmutated_df['days_to_death'] != '-')]
        n_df1 = nonmutated_df[(nonmutated_df['days_to_last_followup'] == '-') & (nonmutated_df['days_to_death'] != '-')][['vital_status','days_to_death']].rename(columns={"days_to_death": "Time"})
        n_df2 = nonmutated_df[nonmutated_df['days_to_last_followup'] != '-'][['vital_status','days_to_last_followup']].rename(columns={"days_to_last_followup": "Time"})
        nonmutated_df = n_df1.append(n_df2)
        nonmutated_df = nonmutated_df[nonmutated_df['vital_status'] != '-']
        nonmutated_df[pr1 + "-" + pr2] = np.ones(len(nonmutated_df))
        final_sruvival_df = pd.concat([final_sruvival_df, nonmutated_df], axis=0, sort=False).fillna(0)
cph = CoxPHFitter()
cph.fit(final_sruvival_df, 'Time', event_col='vital_status')
cph.print_summary()
pickle.dump(cph.summary, open(main_folder + "RESULTS/CoxProportionalHazard/" + cancer_class + "/Multivariate.pickle","wb"))
for ppi in cph.summary.index:
    cph2 = CoxPHFitter()
    cph2.fit(final_sruvival_df[['Time', 'vital_status', ppi]], 'Time', event_col='vital_status')
    pickle.dump(cph2.summary, open(main_folder + "RESULTS/CoxProportionalHazard/" + cancer_class + "/" + ppi + ".pickle", "wb"))


#------------------------------------------------------------------------------------------------------------------------------------------
#Extract g-Gap method input files 
import os
import numpy as np
import pickle

main_folder = "MutSigPPI/"
patients_properties = main_folder + "PatientsSequencesInEachRefSeqChain/"
only_mutated = main_folder + "RefSeq_KeepOnlyMutatedAminoAcidsInAllPatients/"
for patient in os.listdir(patients_properties):
    pfile = open(main_folder + "gGap-PatientsDATA/raw_data/" + patient + ".txt","w")
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
            PDB = chain_1.split("_")
            PDB_id = PDB[0]
            PDB_chain1 = PDB[1]
            PDB_chain2 = chain_2.split("_")[1]
            PDB = PDB_id + "_" + PDB_chain1 + PDB_chain2
            PPI_sequence = ""
        #CHAIN-1    
            chain1_seq = ""
            if(os.path.exists(patients_properties + patient + "/" + prot_1 + "/" + chain_1)):
                address = patients_properties + patient + "/" + prot_1 + "/" + chain_1 + "/mutated_sequence.seq"
                with open(address, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        chain1_seq = l.rstrip()
            else:
                address = only_mutated + prot_1 + "/" + chain_1 + "/mutatedTrimedRefSeq.seq"
                with open(address, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        chain1_seq = l.rstrip()
        #CHAIN-2
            chain2_seq = ""
            if(os.path.exists(patients_properties + patient + "/" + prot_2 + "/" + chain_2)):
                address = patients_properties + patient + "/" + prot_2 + "/" + chain_2 + "/mutated_sequence.seq"
                with open(address, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        chain2_seq = l.rstrip()
            else:
                address = only_mutated + prot_2 + "/" + chain_2 + "/mutatedTrimedRefSeq.seq"
                with open(address, "rb") as seqf:
                    for l in seqf.readlines():
                        l = str(l,'utf-8')
                        chain2_seq = l.rstrip()
            PPI_sequence = chain1_seq + chain2_seq
            pfile.write(">" + PDB + "\n")
            pfile.write(PPI_sequence + "\n")
        pfile.close()

#------------------------------------------------------------------------------------------------------------------------------------------
#Patients' Ratio Mutations (Number of Patients with mutations / Total Number of Patients ; in each cancer class)
import matplotlib.pyplot as plt
import matplotlib
import os
import pickle
import pandas as pd
import seaborn as sns
sns.set(style="whitegrid")

main_folder = "MutSigPPI/"
proteins_refseq_ids = pickle.load(open("MutSigPPI/proteins_refseqs_ids.pickle", "rb"))
patients_of_cancerTypes = pickle.load(open("MutSigPPI/CancerTypePatients.pickle","rb"))
patients_properties = main_folder + "PatientsSequencesInEachRefSeqChain/"
pkmn_type_colors = ['#201923', '#fcff5d', '#7dfc00', '#228c68', '#8ad8e8', '#235b54', '#3998f5', '#3750db', '#f22020', '#991919', '#ffcba5', '#e68f66', '#c56133', '#b732cc', '#f07cab', '#c3a5b4', '#29bdab', '#772b9d', '#0ec434', '#F39C12', '#B2BABB']
column_headers = ["Cancer Class", "Number of patients / Total number of patients", "Protein Protein Interactions"]
values = {}
values["Cancer Class"] = []
values["Number of patients / Total number of patients"] = []
values["Protein Protein Interactions"] = []
ppi_already_passed = {}
for cancer_class in patients_of_cancerTypes:
    print(cancer_class)
    ppi_already_passed = {}
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
            avg_of_mutated_patients = 0.0
            for patient in patients_of_cancerTypes[cancer_class]:
                at_least_one_mutation = False
                if(os.path.exists(patients_properties + patient + "/" + prot_1 + "/" + chain_1)):
                    at_least_one_mutation = True
                if(os.path.exists(patients_properties + patient + "/" + prot_2 + "/" + chain_2)):
                    at_least_one_mutation = True
                if(at_least_one_mutation):
                    avg_of_mutated_patients += 1.0
            avg_of_mutated_patients /= len(patients_of_cancerTypes[cancer_class])
            protein1 = ""
            protein2 = ""
            for pr in proteins_refseq_ids:
                if(prot_1 in proteins_refseq_ids[pr]):
                    protein1 = pr
                if(prot_2 in proteins_refseq_ids[pr]):
                    protein2 = pr
            if(protein1 == ""):
                protein1 = prot_1
            if(protein2 == ""):
                protein2 = prot_2
            if((protein1 + "-" + protein2) in ppi_already_passed):
                ppi_already_passed[protein1 + "-" + protein2] = (avg_of_mutated_patients + ppi_already_passed[protein1 + "-" + protein2])/2
            else:
                ppi_already_passed[protein1 + "-" + protein2] = avg_of_mutated_patients
    for ppi in ppi_already_passed:
        values["Cancer Class"].append(cancer_class)
        values["Number of patients / Total number of patients"].append(ppi_already_passed[ppi])
        values["Protein Protein Interactions"].append(ppi)
df = pd.DataFrame(values, columns = column_headers)
i = 1
j = 0
collection = []
for ppi in ppi_already_passed: 
    collection.append(ppi)
    if(i%37 == 0 and j != 10):
        new_df = df.loc[df["Protein Protein Interactions"].isin(collection)]
        ax = sns.swarmplot(x='Protein Protein Interactions', y='Number of patients / Total number of patients', data=new_df, hue='Cancer Class', size=1.5, palette=pkmn_type_colors)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 6})
        plt.xlabel('Protein-Protein Interactions', fontsize=8)
        plt.ylabel('Number of patients / Total number of patients', fontsize=8)
        plt.xticks(fontsize=4, rotation=90)
        plt.yticks(fontsize=6)
        plt.tight_layout()
        figure = ax.get_figure()    
        figure.savefig('MutSigPPI/RESULTS/PatientsRatioMutations-InEachCancerClass/Mutated-Patients-Ratio-' + str(j) +'.png', dpi=500)
        plt.close()
        collection = [] 
        j += 1
    i += 1
new_df = df.loc[df["Protein Protein Interactions"].isin(collection)]
ax = sns.swarmplot(x='Protein Protein Interactions', y='Number of patients / Total number of patients', data=new_df, hue='Cancer Class', size=1.5, palette=pkmn_type_colors)
#ax = sns.boxplot(x="Type", y="Expression", data=expression_df, width=0.1, palette="Pastel1",linewidth=0.5)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 6})
plt.xlabel('Protein-Protein Interactions', fontsize=8)
plt.ylabel('Number of patients / Total number of patients', fontsize=8)
plt.xticks(fontsize=4, rotation=90)
plt.yticks(fontsize=6)
plt.tight_layout()
figure = ax.get_figure()    
figure.savefig('MutSigPPI/RESULTS/PatientsRatioMutations-InEachCancerClass/Mutated-Patients-Ratio-' + str(j) +'.png', dpi=500)
plt.close()


#------------------------------------------------------------------------------------------------------------------------------------------
#Draw ranges of Hazard Ratio  (FOREST PLOTS)
%matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np
import pickle

cancer_class = "thyroid"
main_folder = "MutSigPPI/"
ppi_hazard = pickle.load(open(main_folder + "RESULTS/CoxProportionalHazard/"+cancer_class+"/Multivariate.pickle", "rb"))
ppi_hazard = ppi_hazard[(ppi_hazard['exp(coef)'] > 2.15) | (ppi_hazard['exp(coef)'] < 1.73)]
#ppi_hazard = ppi_hazard[ppi_hazard['exp(coef)'] > 2.15]
y = list(ppi_hazard.index)
x = list(ppi_hazard['exp(coef)'])
lower_error = list(list(ppi_hazard['exp(coef)']) - ppi_hazard['exp(coef) lower 95%'])
upper_error = list(ppi_hazard['exp(coef) upper 95%'] - list(ppi_hazard['exp(coef)']))
dx = [lower_error,upper_error]
ax = plt.errorbar(x, y, xerr=dx, fmt='.', color='#3d5c5c',ecolor='#94b8b8', elinewidth=2, capsize=0)
plt.axvline(1, ls='--')
plt.tight_layout()
plt.yticks(fontsize=6)
plt.xticks(fontsize=8)
plt.savefig(main_folder + 'RESULTS/CoxProportionalHazard/'+cancer_class+'/ForestPlot.png', dpi=350)


#------------------------------------------------------------------------------------------------------------------------------------------
#Draw ranges of Hazard Ratio 
%matplotlib inline
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np
import pickle

cancer_class = "thyroid"
main_folder = "MutSigPPI/"
ppi_hazard = pickle.load(open(main_folder + "RESULTS/CoxProportionalHazard/"+cancer_class+"/Multivariate.pickle", "rb"))
y = list(ppi_hazard.index)
x = list(ppi_hazard['exp(coef)'])
lower_error = list(list(ppi_hazard['exp(coef)']) - ppi_hazard['exp(coef) lower 95%'])
upper_error = list(ppi_hazard['exp(coef) upper 95%'] - list(ppi_hazard['exp(coef)']))
dx = [lower_error,upper_error]
ax = plt.errorbar(x, y, xerr=dx, fmt='o', color='#5F6A6A',ecolor='#BFC9CA', elinewidth=3, capsize=0)
plt.axvline(1, ls='--')
plt.tight_layout()
plt.yticks(fontsize=6)
plt.xticks(fontsize=8)
plt.savefig(main_folder + 'RESULTS/CoxProportionalHazard/'+cancer_class+'/ForestPlot.png', dpi=350)
