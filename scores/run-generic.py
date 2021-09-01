import seaborn as sns
import numpy as np
import pandas as pd

def gen_scores(data):
  # find size (length) of generic scoring vector
  rws = data.shape[0]

  # initialize empty vectors to score each gene (generic scoring)
  counts_score = np.zeros(rws) #counts
  tau_score = np.zeros(rws) #np.zeros((rws,16)) #tau
  gini_score = np.zeros(rws) #gini
  HS_score = np.zeros(rws) #shannon
  simpson_score = np.zeros(rws) #simpson
  roku_score = np.zeros(rws) #roku
  spmdpm_score = np.zeros(rws) #SPM DPM
  jssdpm_score = np.zeros(rws) #JSS DDM
  #kurtosis_score = np.zeros(rws) #kurtosis

  for i in range(rws):
    counts_score[i] = counts(data[i,:])
    tau_score[i] = tau(data[i,:])
    gini_score[i] = gini(data[i,:])
    HS_score[i] = HS(data[i,:])
    simpson_score[i] = simpson(data[i,:])
    roku_score[i] = ROKU_spec(data[i,:])
    spmdpm_score[i] = SPM_DPM(data[i,:])
    jssdpm_score[i] = JSS_DPM(data[i,:])
    #kurtosis_score[i] = kurto(data[i,:])

  counts_score = counts_score.reshape((rws,1))
  tau_score = tau_score.reshape((rws,1))
  gini_score = gini_score.reshape((rws,1))
  HS_score = HS_score.reshape((rws,1))
  simpson_score = simpson_score.reshape((rws,1))
  roku_score = roku_score.reshape((rws,1))
  spmdpm_score = spmdpm_score.reshape((rws,1))
  jssdpm_score = jssdpm_score.reshape((rws,1))
  #kurtosis_score = kurtosis_score.reshape((rws,1))

  return [counts_score, tau_score, gini_score, HS_score, simpson_score, roku_score, spmdpm_score, jssdpm_score]

def filter_tissue_before(df1, df2, df3, df4, df5, df6):
    list_of_cols = ['amygdala', 'basal_ganglia_basal_ganglion', 'cerebellum', 'cortex', 'hippocampal', 'midbrain',
                    'pituitary', 'eye', 'spinal_cord', 'retina']

    a = [i for i in list_of_cols if i in df1.columns]
    df1.drop(a, axis=1, inplace=True)
    b = [i for i in list_of_cols if i in df2.columns]
    df2.drop(b, axis=1, inplace=True)
    c = [i for i in list_of_cols if i in df3.columns]
    df3.drop(c, axis=1, inplace=True)
    d = [i for i in list_of_cols if i in df4.columns]
    df4.drop(d, axis=1, inplace=True)
    e = [i for i in list_of_cols if i in df5.columns]
    df5.drop(e, axis=1, inplace=True)
    f = [i for i in list_of_cols if i in df6.columns]
    df6.drop(f, axis=1, inplace=True)

    return df1, df2, df3, df4, df5, df5


def reshaping(df):  # reshape data for computation purposes
    df_genes = df.iloc[:, 0:1].to_numpy()
    df = df.set_index('Gene')
    df_head = df.columns.values
    df_head = df_head.reshape((1, len(df_head)))
    return df_genes, df.to_numpy(), df_head

def workflow_generic(h1, h2, h3, d1, d2, d3, filtering='n'):
    #list1 = [h1, h2, h3, d1, d2, d3]
    ht_genes, ht_df, ht_head = reshaping(h1)
    hpc_genes, hpc_df, hpc_head = reshaping(h2)
    hcl_genes, hcl_df, hcl_head = reshaping(h3)

    dt_genes, dt_df, dt_head = reshaping(d1)
    dpc_genes, dpc_df, dpc_head = reshaping(d2)
    dcl_genes, dcl_df, dcl_head = reshaping(d3)

    list2 = [ht_df, hpc_df, hcl_df, dt_df, dpc_df, dcl_df]
    list4 = ['healthy_tissue', 'healthy_primary_cells', 'healthy_cell_line', 'disease_tissue', 'disease_primary_cells', 'disease_cell_line']
    list5 = []

    for df in list2:
        list5.append(gen_scores(df))

    list2 = [ht_df, hpc_df, hcl_df, dt_df, dpc_df, dcl_df]
    list5 = []
    for df in list2:
        list5.append(gen_scores(df))

    dpi = 600
    font_size = 18

    sns.set(rc={'figure.figsize': (21, 7)})
    sns.set_theme(style="whitegrid")
    sns.color_palette("vlag", as_cmap=True)

    for i in range(len(list4)):
        sns.distplot(list5[i][0], hist=False, label=f'{list4[i]}_Counts')
        sns.distplot(list5[i][1], hist=False, label=f'{list4[i]}_Tau')
        sns.distplot(list5[i][2], hist=False, label=f'{list4[i]}_Gini')
        sns.distplot(list5[i][3], hist=False, label=f'{list4[i]}_HS')
        sns.distplot(list5[i][4], hist=False, label=f"{list4[i]}_Simpson")
        sns.distplot(list5[i][5], hist=False, label=f'{list4[i]}_Roku')
        sns.distplot(list5[i][6], hist=False, label=f'{list4[i]}_SPM DPM')
        sns.distplot(list5[i][7], hist=False, label=f"{list4[i]}_JSS DPM")
        plt.legend()
        name = '/scores/distributions/generic/' + 'kde_generic' + f'_{list4[i]}.png'
        plt.ylabel('Density')
        title = 'Density of ' + ' generic scores ' + f' for {list4[i]} data'
        plt.title(title, fontdict={'fontsize': font_size})
        plt.savefig(name, dpi=dpi)
        plt.show()

    ### SAVE GENERIC SCORES AS DATAFRAMES ###
    for i in range(len(list4)):
        scores = np.hstack((ht_genes, list5[i][0], list5[i][1], list5[i][2], list5[i][3], list5[i][4], list5[i][5], list5[i][6], list5[i][7]))
        HPA_scores = pd.DataFrame(data=scores, columns=["Genes", "counts_score", "tau_score", "gini_score", "HS_score", "simpson_score", "roku_score", "spmdpm_score", "jssdpm_score"])
        HPA_scores = HPA_scores.set_index('Genes')
        if filtering == 'y':
            HPA_scores.to_csv("/scores/dataframes/generic/filtered_before/"+"df_tissue_filtered_before"+f"_generic_scores_{list4[i]}.csv", sep='\t', index=True, header=True)
        else:
            HPA_scores.to_csv("/scores/dataframes/generic/filtered_after/" + "df" + f"_generic_scores_{list4[i]}.csv", sep='\t',
                              index=True, header=True)

genes = pd.read_csv("/scores/inputdata/Filter_subcellularLocation_FL75_WLG.csv").to_numpy()
genes = np.insert(genes, 0, 'Genes')
genes = pd.DataFrame(genes.reshape((4025,1)))
new_header = genes.iloc[0]
genes = genes[1:]
genes.columns = new_header
genes.drop_duplicates(inplace=True)
list_of_genes = genes['Genes']

def create_gene_filtered_data(path): # take away the version  of the gene
    # set the path to your file location
    dataframes = []
    for csv in os.listdir(path):
        fullpath = os.path.join(path, csv)
        print(fullpath)
        if os.path.isfile(fullpath):
            # Read a csv into a dataframe and append to a list of dataframes
            dataframe = pd.read_csv(fullpath, sep='\t')
            # dataframe = dataframe.rename(columns={'Unnamed: 0': 'Genes'})
            list1=[]
            a=np.array(dataframe['Genes'])
            for i in range(len(a)):
                list1.append(a[i].split('.')[0])
            dataframe['Genes'] = list1

            dataframe.drop_duplicates(inplace=True)
            dataframe = dataframe[dataframe['Genes'].isin(list_of_genes)]
            df_name = fullpath.split("/")[4]
            print(df_name)
            dataframe.to_csv(os.path.join('/scores/dataframes/generic/filtered_after/', f'gene_filtered_{df_name}.csv'), sep= '\t', index=False)
    return

# Load datasets for TPM (mean of samples)
h1 = pd.read_csv('/scores/inputdata/healthy_tissue_folder/healthy_tissue_folder_mean')
h2 = pd.read_csv('/scores/inputdata/healthy_primary_cells_folder/healthy_primary_cells_folder_mean')
h3 = pd.read_csv('/scores/inputdata/healthy_cell_line_folder/healthy_cell_line_folder_mean')
d1 = pd.read_csv('/scores/inputdata/disease_tissue_folder/disease_tissue_folder_mean')
d2 = pd.read_csv('/scores/inputdata/disease_primary_cells_folder/disease_primary_cells_folder_mean')
d3 = pd.read_csv('/scores/inputdata/disease_cell_line_folder/disease_cell_line_folder_mean')

tissue_filtering_before = input('Do you want tissue filtering before calculating the score? (y/n)')
if tissue_filtering_before == 'y':
    h1, h2, h3, d1, d2, d3 = filter_tissue_before(h1, h2, h3, d1, d2, d3)
    workflow_generic(h1, h2, h3, d1, d2, d3, filtering='y')
    gene_filtering_after = input('Do you want subcellular location filtering after calculating the score? (y/n)')
    if gene_filtering_after == 'y':
        create_gene_filtered_data("/scores/dataframes/specific/filtered_before/")
        sys.exit()
else:
    workflow_generic(h1, h2, h3, d1, d2, d3, filtering='n')
    gene_filtering_after = input('Do you want subcellular location filtering after calculating the score? (y/n)')
        if gene_filtering_after == 'y':
            create_gene_filtered_data("/scores/dataframes/specific/filtered_after/")
            sys.exit()
