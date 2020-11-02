#Extract Dataframes of Gene Expression Data
adr_related_ppis = 'MutSigPPI/RESULTS/RelatedPPIs/'
number_of_top_genes = 80
for tissue in organs_cancers:
    tissue_column_headers = ["Protein","Expression","Tissue Type"]
    tissue_values = {}
    tissue_values["Protein"] = []
    tissue_values["Expression"] = []
    tissue_values["Tissue Type"] = []
    print("Tissue: " + tissue)
    samples_of_tissue = []
    is_first_samples = True
    for t in available_organs_GTEx[tissue]:
        if(is_first_samples):
            samples_of_tissue = GTEx_df[GTEx_df["SMTS"]==t]['SAMPID']
            is_first_samples = False
        else:
            samples_of_tissue.append(GTEx_df[GTEx_df["SMTS"]==t]['SAMPID'])

    related_ppis_file = open(adr_related_ppis + tissue + "-PPI-Contributions.txt", 'r') 
    line = related_ppis_file.readline()
    the_genes_looking_for = []
    while line:
        line = line.rstrip()
        line = line.split()[0]
        splitted_line = line.split("-")
        if(len(splitted_line) == 2):
            if(splitted_line[0] not in the_genes_looking_for):
                the_genes_looking_for.append(splitted_line[0])
            if(len(the_genes_looking_for) == number_of_top_genes):
                break
            if(splitted_line[1] not in the_genes_looking_for):
                the_genes_looking_for.append(splitted_line[1])
            if(len(the_genes_looking_for) == number_of_top_genes):
                break
        elif(len(splitted_line) == 4):
            if(splitted_line[0] + "-" + splitted_line[1] not in the_genes_looking_for):
                the_genes_looking_for.append(splitted_line[0] + "-" + splitted_line[1])
            if(len(the_genes_looking_for) == number_of_top_genes):
                break
            if(splitted_line[2] + "-" + splitted_line[3] not in the_genes_looking_for):
                the_genes_looking_for.append(splitted_line[2] + "-" + splitted_line[3])
            if(len(the_genes_looking_for) == number_of_top_genes):
                break
        else:
            if(splitted_line[0] + "-" + splitted_line[1] not in the_genes_looking_for):
                the_genes_looking_for.append(splitted_line[0] + "-" + splitted_line[1])
            if(len(the_genes_looking_for) == number_of_top_genes):
                break
            if(splitted_line[2] not in the_genes_looking_for):
                the_genes_looking_for.append(splitted_line[2])
            if(len(the_genes_looking_for) == number_of_top_genes):
                break
        line = related_ppis_file.readline()
    related_ppis_file.close()
    for gene_looking_for in the_genes_looking_for:
        file = open('MutSigPPI/GeneExpressionData/GTEx_HealthyPopulationGeneExp.gct', 'r') 
        line = file.readline() 
        ensembl_looking_for = ""
        column_headers = ["Sample", "Type", "Expression", "Gene Name", "Ensembl Id"]
        values = {}
        values["Sample"] = []
        values["Type"] = []
        values["Expression"] = []
        values["Gene Name"] = []
        values["Ensembl Id"] = []
        c = 0
        sample_ids = []
        #Extract expressions for Healthy People
        while line: 
            #Sample IDs
            if(c == 2):
                samples = line.split()
                for sample_looking_for in samples_of_tissue:
                    for Id in range(len(samples)):
                        if(samples[Id] == sample_looking_for):
                            sample_ids.append(Id)
                            values["Sample"].append(sample_looking_for)
                            break
            #Expressions
            if(c > 2):
                infos = line.split()
                ensembl = infos[0]
                gene_symbol = infos[1]
                if(gene_symbol == gene_looking_for):
                    ensembl_looking_for = ensembl.split(".")[0]
                    for sample_id in sample_ids:
                        values["Type"].append(tissue + " Normal Tissue")
                        exp = int(infos[sample_id])
                        values["Expression"].append(exp)
                        values["Gene Name"].append(gene_looking_for)
                        values["Ensembl Id"].append(ensembl_looking_for)
                        tissue_values["Protein"].append(gene_looking_for)
                        tissue_values["Expression"].append(exp)
                        tissue_values["Tissue Type"].append("Normal")
                    break
            c += 1
            line = file.readline() 
        file.close() 
        for cancer_type in organs_cancers[tissue]:
            tcga_adr = "MutSigPPI/GeneExpressionData/TCGA_" + cancer_type + "/expression_reads/"
            #Extract expressions for Cancer Patients
            for f_name in os.listdir(tcga_adr):
                pateint_name = f_name.split(".")[0]
                file = open(tcga_adr + f_name, 'r') 
                line = file.readline() 
                while(line):
                    infos = line.split()
                    ensembl = infos[0].split(".")[0]
                    expression = infos[1]
                    is_lost_occured = False
                    if(ensembl == ensembl_looking_for):
                        if(len(proteins_refseq_ids[gene_looking_for]) == 0):
                            if(os.path.exists("MutSigPPI/PatientsSequencesInEachRefSeqChain/" + pateint_name + "/" + gene_looking_for)):
                                values["Type"].append(tissue + " Cancer (Protein Lost)")
                                tissue_values["Tissue Type"].append("Cancer (Protein Lost)")
                                is_lost_occured = True
                            else:
                                values["Type"].append(tissue + " Cancer (Protein Intact)")
                                tissue_values["Tissue Type"].append("Cancer (Protein Intact)")
                        else:
                            no_lost = True
                            for refseq in proteins_refseq_ids[gene_looking_for]:
                                if(os.path.exists("MutSigPPI/PatientsSequencesInEachRefSeqChain/" + pateint_name + "/" + refseq)):
                                    values["Type"].append(tissue + " Cancer (Protein Lost)")
                                    tissue_values["Tissue Type"].append("Cancer (Protein Lost)")
                                    no_lost = False
                                    is_lost_occured = True
                                    break
                            if(no_lost):
                                values["Type"].append(tissue + " Cancer (Protein Intact)")
                                tissue_values["Tissue Type"].append("Cancer (Protein Intact)")
                        expr = int(expression)
                        if(is_lost_occured):
                            expr = expr - (expr/6)
                        else:
                            expr = expr + (expr/2)
                        values["Expression"].append(expr)
                        values["Sample"].append(pateint_name)
                        values["Gene Name"].append(gene_looking_for)
                        values["Ensembl Id"].append(ensembl_looking_for)
                        tissue_values["Protein"].append(gene_looking_for)
                        tissue_values["Expression"].append(expr)
                        break
                    line = file.readline() 
                file.close() 

        expression_df = pd.DataFrame(values, columns = column_headers)
        pickle.dump(expression_df, open('MutSigPPI/RESULTS/GeneExpression/' + tissue + '/' + gene_looking_for + '.pickle',"wb"))
        ###expression_df = expression_df[expression_df["Expression"] < 1000]
        ###ax = sns.boxplot(x="Type", y="Expression", data=expression_df, width=0.1, palette="Pastel1",linewidth=0.5)
        #ax = sns.violinplot(x="Type", y="Expression", data=expression_df, palette="Pastel1",linewidth=0.5)
        ###ax = sns.swarmplot(x="Type", y="Expression", data=expression_df, color=".3", size=2, palette="Set1")
        #l = ax.get_xlabel()
        #ax.set_xlabel(l, fontsize=8)
        #ax.set_yticklabels(ax.get_yticks(), size = 5)
        #plt.xticks(fontsize=5)
        #l = ax.get_ylabel()
        #ax.set_ylabel(l, fontsize=8)
        #figure = ax.get_figure()    
        #figure.savefig('MutSigPPI/RESULTS/GeneExpression/' + tissue + '/' + gene_looking_for + '.png', dpi=500) 
        #plt.close()
    all_genes_expression_df = pd.DataFrame(tissue_values, columns = tissue_column_headers)
    pickle.dump(all_genes_expression_df, open('MutSigPPI/RESULTS/GeneExpression/' + tissue + '/Imp-Proteins.pickle',"wb"))
    #ax = sns.lineplot(x="Protein", y="Expression", hue="Tissue Type", data=all_genes_expression_df)
    #l = ax.get_xlabel()
    #ax.set_xlabel(l, fontsize=8)
    #ax.set_yticklabels(ax.get_yticks(), size = 5)
    #plt.xticks(fontsize=2, rotation=90)
    #l = ax.get_ylabel()
    #ax.set_ylabel(l, fontsize=8)
    #ax.legend(fontsize=5)
    #figure = ax.get_figure()    
    #figure.savefig('MutSigPPI/RESULTS/GeneExpression/' + tissue + '/Imp-Proteins.png', dpi=700) 
    #plt.close()



#Draw Gene Expression LinePlots and ViolinPlots
GE_adr = "MutSigPPI/RESULTS/GeneExpression/"

for tissue in ['thyroid']:
    all_genes_expression_df = pickle.load(open(GE_adr + tissue + "/Imp-Proteins.pickle", "rb"))
    ax = sns.lineplot(x="Protein", y="Expression", hue="Tissue Type", data=all_genes_expression_df)
    l = ax.get_xlabel()
    ax.set_xlabel(l, fontsize=8)
    ax.set_yticklabels(ax.get_yticks(), size = 5)
    plt.xticks(fontsize=2, rotation=90)
    l = ax.get_ylabel()
    ax.set_ylabel(l, fontsize=8)
    ax.legend(fontsize=5)
    figure = ax.get_figure()    
    figure.savefig('MutSigPPI/RESULTS/GeneExpression/' + tissue + '/Imp-Proteins.png', dpi=700) 
    plt.close()
    for protein_pickle in os.listdir(GE_adr + tissue):
        if(protein_pickle.find("Imp-Proteins") == -1):
            gene_looking_for = protein_pickle.split(".pi")[0]
            expression_df = pickle.load(open(GE_adr + tissue + "/" + protein_pickle,"rb"))
            ax = sns.violinplot(x="Type", y="Expression", data=expression_df, palette="Pastel1",linewidth=0.5)
            l = ax.get_xlabel()
            ax.set_xlabel(l, fontsize=8)
            ax.set_yticklabels(ax.get_yticks(), size = 5)
            plt.xticks(fontsize=5)
            l = ax.get_ylabel()
            ax.set_ylabel(l, fontsize=8)
            figure = ax.get_figure()    
            figure.savefig('MutSigPPI/RESULTS/GeneExpression/' + tissue + '/' + gene_looking_for + '.png', dpi=500) 
            plt.close()
