# Utility functions for use with scanpy

# Calculate the number of PCs that contain some proportion (default is 95%) of the variance
def npcs(
  ADATA,
  var_perc=0.95,
  reduction="pca"
):
    from numpy import sum, var
    get_var = lambda i: var(ADATA.obsm[reduction][:,i])

    if ADATA.obsm[reduction] is None:
        print(f"Reduction '{reduction}', not found!")
        return None
    else:
        var_tmp = [get_var(i) for i in list(range(0,ADATA.obsm[reduction].shape[1]))]
        var_cut = var_perc * sum(var_tmp)
        n_pcs = 0
        var_sum = 0
        while var_sum<var_cut and n_pcs<ADATA.obsm[reduction].shape[1]-1:
            var_sum = var_sum + var_tmp[n_pcs]
            n_pcs = n_pcs + 1

        return(n_pcs)


# Re-order a dimensions of a reduction by decreasing % variance
def reorder_reduction(
    ADATA,
    reduction="pca",
    verbose=False
):
    from numpy import var, argsort

    if reduction in ADATA.obsm:
        get_var = lambda i: var(ADATA.obsm[reduction][:,i])
        var_tmp = [get_var(i) for i in list(range(0,ADATA.obsm[reduction].shape[1]))]
        if verbose:
            print("Reduction variance by dimension:")
            print(var_tmp)

        pc_order = argsort(var_tmp)[::-1]
        ADATA.obsm[reduction] = ADATA.obsm[reduction][:,pc_order]
    else:
        print(f"The reduction '{reduction}' was not found...")


# Read in a list of gene lists from .csv (each column is a gene list)
## Useful for plotting
def read_csv_to_dict(filename, names2check=""):
    import csv

    # Open the CSV file
    with open(filename, 'r') as file:

        # Create a CSV reader object
        reader = csv.reader(file)

        # Read the first row as header
        header = next(reader)

        # Create an empty dictionary to store the columns
        dict_out = {col: [] for col in header}

        # Loop through each row in the CSV file
        for row in reader:

            # Loop through each column in the row
            for col, value in zip(header, row):

                # Add the value to the corresponding column in the dictionary
                if value: # skip empty strings
                    dict_out[col].append(value)

    # Filter out unwanted entries based on the list in `names2check`
    if len(names2check) > 1:
        for KEY in dict_out.keys():
            dict_out[KEY] = [k for k in dict_out[KEY] if k in names2check]

    # Return the dictionary
    return dict_out



# Function to export DGEA results to a .csv file
def export_dgea_to_csv(
    adata, # adata
    dgea_name, # 'rank_gene_groups_clusters_leiden_0.5', name in adata.uns[]
    n_features,
    csv_out,
    axis=0,     # how to write results for each group (1=horizontal, 0=vertical)
    wide=False
):    
    import pandas as pd
    import scanpy as sc

    result = adata.uns[dgea_name]
    groups = result['names'].dtype.names

    if wide:
        celltype_markers = pd.DataFrame(
            {group + '_' + key[:-1]: result[key][group]
            for group in groups for key in ['names', 'logfoldchanges','pvals']}).head(n_features)
        celltype_markers.to_csv(csv_out, index=False)
    else:
        marker_list = list()
        for group in adata.uns[dgea_name]['names'].dtype.names:
            markers = sc.get.rank_genes_groups_df(adata, key=dgea_name, group = group).head(n_features)
            markers['celltypes'] = group
            marker_list.append(markers)
        
        celltype_markers = pd.concat(
            marker_list, 
            axis=axis
        )
        celltype_markers.to_csv(
            csv_out,
            index=False
        )

# Function to convert feature names 
def convert_feature_names(
        adata: ad.AnnData, 
        gtf_info: pd.DataFrame, 
        from_col: str='GENEID',
        to_col: str='GeneSymbol',
        inplace: bool=True,
        verbose: bool=True
) -> ad.AnnData:
    if not inplace:
        adata = adata.copy()

    # Filter gtf_info to keep only the gene names found in the anndata object
    gtf_info_filtered = gtf_info[gtf_info[from_col].isin(adata.var_names)]
    
    if verbose:
        num_found = len(gtf_info_filtered)
        num_total = len(adata.var_names)
        # fraction_found = num_found / num_total
        print(f"Fraction of adata.var_names found in gtf_info[{from_col}]: {num_found} out of {num_total}")
    
    gene_name_mapping = dict(zip(
        gtf_info[from_col], 
        gtf_info[to_col]
        ))
    
    adata.var[from_col] = adata.var_names
    adata.var[to_col] = adata.var[from_col].map(gene_name_mapping)

    adata.var.dropna(subset=[to_col], inplace=True)
    adata.var.reset_index(drop=True, inplace=True)

    # mask = adata.var_names.isin(adata.var[to_col].values)
    # adata = adata[:, mask].copy()
    adata.var_names = adata.var[to_col]
    adata.var_names_make_unique()

    if not inplace:
        return adata