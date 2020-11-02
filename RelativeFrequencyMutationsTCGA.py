#------------------------------------------------------------------------------------------------------------------------------------------
#Relative Frequency of 96 Mutations of TCGA dataset for different cancer classes
import pandas as pd 
import matplotlib.pyplot as plt
import pickle

main_folder = "MutSigPPI/"
cancer_classes = ['thyroid', 'liver', 'male reproductive system', 'adrenal gland', 'GYN', 'head and neck', 'Lung', 'other', 'Melanoma', 'thymus', 'CRC', 'bone marrow', 'pancreatic', 'stomach', 'Kidney', 'bladder', 'breast', 'brain', 'pleura', 'esophagus', 'bile duct']
mutations_df = pd.read_csv(main_folder + "MutationalSignature/Mutation-Counts-In-Each-Patient.txt", delimiter="\t") 
tcga_clinical_dataframe = pickle.load(open(main_folder + "TCGA_clinical_dataframe.pickle","rb"))
tcga_clinical_dataframe = tcga_clinical_dataframe[['cancer_class']]
mutations_df = mutations_df.set_index('Mutation Types')
mutations_df = mutations_df.T
first_round = True
final_df = ""
for cancer_class in cancer_classes:
    cancer_mut_info = mutations_df.loc[list(tcga_clinical_dataframe[tcga_clinical_dataframe['cancer_class'] == cancer_class].index)].sum().to_frame().rename(columns={0: cancer_class}).T
    if(first_round):
        first_round = False
        final_df = cancer_mut_info
    else:
        final_df = pd.concat([final_df, cancer_mut_info], axis=0, sort=False)
in_order_columns = ['A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T',
 'A[C>G]A', 'A[C>G]C', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T', 'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T', 'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T', 'A[C>G]G', 'A[C>G]T',
 'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T',
 'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T',
 'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T', 'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T',
 'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T']        
final_df = final_df[in_order_columns]                                            
stacked_data = final_df.apply(lambda x: x*100/sum(x), axis=1)
all_colors = ['#ffe6f0','#ffcce0','#ffb3d1','#ff99c2','#ff80b3','#ff66a3','#ff4d94','#ff3385','#ff1a75','#ff0066','#e6005c','#cc0052','#b30047','#99003d','#800033','#660029',
              '#ccffdd','#b3ffcc','#99ffbb','#80ffaa','#66ff99','#4dff88','#33ff77','#1aff66','#00ff55','#00e64d','#00cc44','#00b33c','#009933','#00802b','#006622','#004d1a', 
              '#cce0ff','#b3d1ff','#99c2ff','#80b3ff','#66a3ff','#4d94ff','#3385ff','#1a75ff','#0066ff','#005ce6','#0052cc','#0047b3','#003d99','#003380','#002966','#001f4d',
              '#ffffe6','#ffffcc','#ffffb3','#ffff99','#ffff80','#ffff66','#ffff4d','#ffff33','#ffff1a','#ffff00','#e6e600','#cccc00','#b3b300','#999900','#808000','#666600',
              '#fff2e6','#ffe6cc','#ffd9b3','#ffcc99','#ffbf80','#ffb366','#ffa64d','#ff9933','#ff8c1a','#ff8000','#e67300','#cc6600','#b35900','#994d00','#804000','#663300', 
              '#f5e6ff','#ebccff','#e0b3ff','#d699ff','#cc80ff','#c266ff','#b84dff','#ad33ff','#a31aff','#9900ff','#8a00e6','#7a00cc','#6b00b3','#5c0099','#4c0080','#3d0066']
stacked_data.plot(kind="bar", stacked=True, color=all_colors)
plt.title("Relative Frequency of base substitutions")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel("Cancer Class")
plt.ylabel("Relative Frequency (%)")
plt.savefig(main_folder + 'RESULTS/Relative_Frequency.png', dpi = 500, bbox_inches = "tight")

