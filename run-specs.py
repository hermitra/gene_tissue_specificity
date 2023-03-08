import multiprocessing
from tissuespecific import SPECS
from joblib import Parallel, delayed

num_cores = multiprocessing.cpu_count()

def create_data_specs(path):
    # set the path to your file location
    dataframes = []
    for csv in os.listdir(path):
        fullpath = os.path.join(path, csv)
        if os.path.isfile(fullpath):  # csv.startswith('.') and
            # Read a csv into a dataframe and append to a list of dataframes
            dataframe = pd.read_csv(fullpath)
            tissue = fullpath.split("nb")[1]
            #print(tissue)
            new_col_names = [f'{tissue}' + f'_{x}' for x in dataframe.columns]
            new_col_names[0] = 'Genes'
            #print(new_col_names)
            dataframe = dataframe.set_axis(new_col_names, axis=1)
            dataframes.append(dataframe)
            #print('_______________________________________________________________________________')

    # Concatenate all created dataframes into one
    df = reduce(lambda left, right: pd.merge(left, right, on=['Genes'],
                                            how='outer'), dataframes)
    # create the dataframes and save
    df_name = path.split("/")[6]

    df_disease_tissue_bis = df.set_index('Genes')
    idx_disease_tissue = df_disease_tissue_bis.iloc[:, 0:0]
    idx_disease_tissue.to_csv(os.path.join('/scores/specs/', f'idx_{df_name}.txt'), sep='\t')
    df_disease_tissue_bis.to_csv(os.path.join('/scores/specs/', f'df_{df_name}.txt'), sep='\t', header=True, index=True)

    return ()

list1 = ['scores/inputdata/healthy_tissue_folder/healthy_tissue_folder_mean', 'scores/inputdata/healthy_primary_cells_folder/healthy_primary_cells_folder_mean',
         'scores/inputdata/healthy_cell_line_folder/healthy_cell_line_folder_mean', 'scores/inputdata/disease_tissue_folder/disease_tissue_folder_mean',
         'scores/inputdata/disease_primary_cells_folder/disease_primary_cells_folder_mean', 'scores/inputdata/disease_cell_line_folder/disease_cell_line_folder_mean']
list2 = ['scores/inputdata/healthy_tissue_folder_filtered/healthy_tissue_folder_mean', 'scores/inputdata/healthy_primary_cells_folder_filtered/healthy_primary_cells_folder_mean',
         'scores/inputdata/healthy_cell_line_folder_filtered/healthy_cell_line_folder_mean', 'scores/inputdata/disease_tissue_folder_filtered/disease_tissue_folder_mean',
         'scores/inputdata/disease_primary_cells_folder_filtered/disease_primary_cells_folder_mean', 'scores/inputdata/disease_cell_line_folder_filtered/disease_cell_line_folder_mean']

tissue_filtering_before = input('Do you want tissue filtering before calculating the score? (y/n)')
if tissue_filtering_before == 'y':
    Parallel(n_jobs=num_cores)(delayed(create_data_specs)(path)
                               for path in list2)
else:
    Parallel(n_jobs=num_cores)(delayed(create_data_specs)(path)
                               for path in list1)

# finally compute specs
h11= open('/scores/specs/df_healthy_tissue.txt', 'r')
h12= open('/scores/specs/idx_healthy_tissue.txt', 'r')
h21= open('/scores/specs/df_healthy_primary_cells.txt', 'r')
h22= open('/scores/specs/idx_healthy_primary_cells.txt', 'r')
h31= open('/scores/specs/df_healthy_cell_line.txt', 'r')
h32= open('/scores/specs/idx_healthy_cell_line.txt', 'r')

d11= open('/scores/specs/df_disease_tissue.txt', 'r')
d12= open('/scores/specs/idx_disease_tissue.txt', 'r')
d21= open('/scores/specs/df_disease_primary_cells.txt', 'r')
d22= open('/scores/specs/idx_disease_primary_cells.txt', 'r')
d31= open('/scores/specs/df_disease_cell_line.txt', 'r')
d32= open('/scores/specs/idx_disease_cell_line.txt', 'r')

list3=[(h11, h12), (h21, h22), (h31, h32), (d11, d12), (d21, d22), (d31, d32)]
Parallel(n_jobs=num_cores)(delayed(SPECS)((i,j)) for (i,j) in list3)
