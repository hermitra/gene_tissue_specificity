import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
from scores.tissuespecific import TSI_spec, zscore, SPM, JSS


def spec_scores(data):
    # find size (length) of specific scoring vector
    rws = data.shape[0]
    cls = data.shape[1]

    # initialize empty arrays to score each gene per tissue (specific scoring)
    tsi_score = np.zeros((rws, cls))  # TSI
    z_score = np.zeros((rws, cls))  # Z-Score
    spm_score = np.zeros((rws, cls))  # SPM
    jss_score = np.zeros((rws, cls))  # JSS
    #specs_score = np.zeros((rws, cls))  # SPECS

    for i in range(rws):
        tsi_score[i, :] = TSI_spec(data[i, :])
        z_score[i] = zscore(data[i, :])
        spm_score[i] = SPM(data[i, :])
        jss_score[i] = JSS(data[i, :])
        #specs_score[i] = SPECS(HPA[i, :])

    return [tsi_score, z_score, spm_score, jss_score]


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


def workflow(h1, h2, h3, d1, d2, d3, filtering='n'):
    #list1 = [h1, h2, h3, d1, d2, d3]
    ht_genes, ht_df, ht_head = reshaping(h1)
    hpc_genes, hpc_df, hpc_head = reshaping(h2)
    hcl_genes, hcl_df, hcl_head = reshaping(h3)

    dt_genes, dt_df, dt_head = reshaping(d1)
    dpc_genes, dpc_df, dpc_head = reshaping(d2)
    dcl_genes, dcl_df, dcl_head = reshaping(d3)

    list2 = [ht_df, hpc_df, hcl_df, dt_df, dpc_df, dcl_df]
    list3 = []
    for df in list2:
        list3.append(spec_scores(df))

    ### SAVE SPECIFIC SCORES AS DATAFRAMES ###
    list4 = ['healthy_tissue', 'healthy_primary_cells', 'healthy_cell_line', 'disease_tissue', 'disease_primary_cells',
             'disease_cell_line']
    list6 = [ht_genes, hpc_genes, hcl_genes, dt_genes, dpc_genes, dcl_genes]
    list7 = [ht_head, hpc_head, hcl_head, dt_head, dpc_head, dcl_head]

    for i in range(len(list4)):
        ### ADD TISSUE NAMES ###
        tsi_scoore_no = np.vstack((list7[i], list3[i][0]))
        z_scoore_no = np.vstack((list7[i], list3[i][1]))
        spm_scoore_no = np.vstack((list7[i], list3[i][2]))
        jss_scoore_no = np.vstack((list7[i], list3[i][3]))

        ### RESHAPE GENES VECTOR ###
        HPA_genes_2 = np.insert(list6[i], 0, 'Genes')
        HPA_genes_3 = HPA_genes_2.reshape((63857, 1))

        ### CREATE ARRAY WITH GENE NAMES ###
        tsi_score_3_no = np.hstack((HPA_genes_3, tsi_scoore_no))
        z_score_3_no = np.hstack((HPA_genes_3, z_scoore_no))
        spm_score_3_no = np.hstack((HPA_genes_3, spm_scoore_no))
        jss_score_3_no = np.hstack((HPA_genes_3, jss_scoore_no))

        ### TRANSFORM ARRAY TO DATAFRAME
        tsi_score_4_no = pd.DataFrame(tsi_score_3_no)
        z_score_4_no = pd.DataFrame(z_score_3_no)
        spm_score_4_no = pd.DataFrame(spm_score_3_no)
        jss_score_4_no = pd.DataFrame(jss_score_3_no)

        ### MAKE FIRST ROW HEADER
        new_header_no = tsi_score_4_no.iloc[0]  # grab the first row for the header

        tsi_score_4_no = tsi_score_4_no[1:]  # take the data less the header row
        z_score_4_no = z_score_4_no[1:]
        spm_score_4_no = spm_score_4_no[1:]
        jss_score_4_no = jss_score_4_no[1:]

        tsi_score_4_no.columns = new_header_no  # set the header row as the df header
        z_score_4_no.columns = new_header_no
        spm_score_4_no.columns = new_header_no
        jss_score_4_no.columns = new_header_no

        tsi_score_4_no = tsi_score_4_no.fillna(0)
        z_score_4_no = z_score_4_no.fillna(0)
        spm_score_4_no = spm_score_4_no.fillna(0)
        jss_score_4_no = jss_score_4_no.fillna(0)

        tsi_score_5_no = tsi_score_4_no.set_index('Genes')
        z_score_5_no = z_score_4_no.set_index('Genes')
        spm_score_5_no = spm_score_4_no.set_index('Genes')
        jss_score_5_no = jss_score_4_no.set_index('Genes')

        ### KERNEL DENSITY ESTIMATES OF SPECIFIC SCORES ###

        #specs1 = pd.read_csv('/scores/dataframes/specs/specs_healthy_tissue.csv',sep='\t', skiprows=[0]).iloc[:, 1:].to_numpy()
        #specs2 = pd.read_csv('/scores/dataframes/specs/specs_healthy_primary_cells.csv', sep='\t',skiprows=[0]).iloc[:, 1:].to_numpy()
        #specs3 = pd.read_csv('/scores/dataframes/specs/specs_healthy_cell_line.csv',sep='\t', skiprows=[0]).iloc[:, 1:].to_numpy()
        #specs4 = pd.read_csv('/scores/dataframes/specs/specs_disease_tissue.csv',sep='\t', skiprows=[0]).iloc[:, 1:].to_numpy()
        #specs5 = pd.read_csv('/scores/dataframes/specs/specs_disease_primary_cells.csv', sep='\t',skiprows=[0]).iloc[:, 1:].to_numpy()
        #specs6 = pd.read_csv('/scores/dataframes/specs/specs_disease_cell_line.csv',sep='\t', skiprows=[0]).iloc[:, 1:].to_numpy()

        dpi = 600
        font_size = 18

        sns.set(rc={'figure.figsize': (21, 7)})
        sns.set_theme(style="whitegrid")
        sns.color_palette("vlag", as_cmap=True)

        list4 = ['healthy_tissue', 'healthy_primary_cells', 'healthy_cell_line', 'disease_tissue', 'disease_primary_cells', 'disease_cell_line']
        #list5 = [specs1, specs2, specs3, specs4, specs5, specs6]
        for i in range(len(list4)):
            sns.distplot(list3[i][0], hist=False, label=f'{list4[i]}_Tsi')
            sns.distplot(list3[i][1], hist=False, label=f'{list4[i]}_Z')
            sns.distplot(list3[i][2], hist=False, label=f'{list4[i]}_SPM')
            sns.distplot(list3[i][3], hist=False, label=f'{list4[i]}_JSS')
            #sns.distplot(list5[i], hist=False, label=f"{list4[i]}_SPECs")

            plt.legend()
            name = '/scores/distributions/specific/' + 'kde_tissue_specific' + f'_{list4[i]}.png'
            plt.ylabel('Density')
            title = 'Density of ' + ' specific scores ' + f' for {list4[i]} data'
            plt.title(title, fontdict={'fontsize': font_size})
            plt.savefig(name, dpi=dpi)
            plt.show()

        if filtering == 'y':
            pd.DataFrame(tsi_score_5_no).to_csv(os.path.join("/scores/dataframes/specific/filtered_before/" + "df_tissue_filtered_before" + f"" + f"_tsi_score_{list4[i]}.csv"),
                                            sep='\t', index=True, header=True)
            pd.DataFrame(z_score_5_no).to_csv(os.path.join("/scores/dataframes/specific/filtered_before/" + "df_tissue_filtered_before" + f"_z_score_{list4[i]}.csv"),
                                          sep='\t', index=True, header=True)
            pd.DataFrame(spm_score_5_no).to_csv(os.path.join("/scores/dataframes/specific/filtered_before/" + "df_tissue_filtered_before" + f"_spm_score_{list4[i]}.csv"),
                                            sep='\t', index=True, header=True)
            pd.DataFrame(jss_score_5_no).to_csv(os.path.join("/scores/dataframes/specific/filtered_before/" + "df_tissue_filtered_before" + f"_jss_score_{list4[i]}.csv"),
                                            sep='\t', index=True, header=True)
        else:
            pd.DataFrame(tsi_score_5_no).to_csv(os.path.join("/scores/dataframes/specific/filtered_after/" + "df" + f"" + f"_tsi_score_{list4[i]}.csv"),
                                            sep='\t', index=True, header=True)
            pd.DataFrame(z_score_5_no).to_csv(os.path.join("/scores/dataframes/specific/filtered_after/" + "df" + f"_z_score_{list4[i]}.csv"),
                                          sep='\t', index=True, header=True)
            pd.DataFrame(spm_score_5_no).to_csv(os.path.join("/scores/dataframes/specific/filtered_after/" + "df" + f"_spm_score_{list4[i]}.csv"),
                                            sep='\t', index=True, header=True)
            pd.DataFrame(jss_score_5_no).to_csv(os.path.join("/scores/dataframes/specific/filtered_after/" + "df" + f"_jss_score_{list4[i]}.csv"),
                                            sep='\t', index=True, header=True)

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
            dataframe.to_csv(os.path.join('/scores/dataframes/specific/filtered_after/', f'gene_filtered_{df_name}.csv'), sep= '\t', index=False)
    return

def create_tissue_filtered_data(path):
    # set the path to your file location
    dataframes = []
    list_of_cols = ['amygdala', 'basal_ganglia_basal_ganglion', 'cerebellum', 'cortex', 'hippocampal', 'midbrain',
                    'pituitary', 'eye', 'spinal_cord', 'retina']
    for csv in os.listdir(path):
        fullpath = os.path.join(path, csv)
        if os.path.isfile(fullpath):
            dataframe = pd.read_csv(fullpath, sep='\t')
            #if 'amygdala' in dataframe.columns:
                #dataframe.drop(['amygdala'], axis=1 inplace=True)
            a = [i for i in list_of_cols if i in dataframe.columns ]
            dataframe.drop(a, axis=1, inplace=True)
            # create a dataframe and save
            df_name = fullpath.split("/")[8]
            dataframe.to_csv(os.path.join('/scores/dataframes/specific/filtered_after/', f'tissue_filtered_after_{df_name}.csv'), sep='\t', index=False)
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
    workflow(h1, h2, h3, d1, d2, d3, filtering='y')
    gene_filtering_after = input('Do you want subcellular location filtering after calculating the score? (y/n)')
    if gene_filtering_after == 'y':
        create_gene_filtered_data("/scores/dataframes/specific/filtered_before/")
        sys.exit()
else:
    workflow(h1, h2, h3, d1, d2, d3, filtering='n')
    tissue_filtering_after = input('Do you want tissue filtering after calculating the score? (y/n)')
    if tissue_filtering_after == 'y':
        create_tissue_filtered_data("/scores/dataframes/specific/")
        gene_filtering_after = input('Do you want subcellular location filtering after calculating the score? (y/n)')
        if gene_filtering_after == 'y':
            create_gene_filtered_data("/scores/dataframes/specific/filtered_after/")
            sys.exit()
    else:
        gene_filtering_after = input('Do you want subcellular location filtering after calculating the score? (y/n)')
        if gene_filtering_after == 'y':
            create_gene_filtered_data("/scores/dataframes/specific/filtered_after/")
            sys.exit()