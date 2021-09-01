from tissuespecific import SPECS
#joblib parallel

def parallelize(x):
    # get the paths of your dataframe
    df_paths = [scores/dataframes/ , , , , , ]

    df_path_chosen = df_path[x]

    df = pd.read_csv(df_path_chosen)

    df.tocsv(f'dsfgfsgh_{x}.csv')

num_cores = multiprocessing.cpu_count()
Parallel(n_jobs=num_cores)(delayed(parallelize)(x)
                           for x in range(len(6)))é

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
    #df.to_csv(os.path.join('/scores/specs/', f'{df_name}.csv'))
    return ()

d11= open('/scores/specs/df_healthy_tissue_filtered.txt', 'r')#disease tissue df
d12= open('/scores/specs/idx_healthy_tissue_filtered.txt', 'r')#disease tissue index

for (d11,d12) in
    list1.append(SPECS(d11,d12))