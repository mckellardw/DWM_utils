import os
import pandas as pd

def gene_count_dict(adata):
    genes = adata.var_names
    counts = adata.X.toarray().flatten()
    
    # Only include genes that have counts > 0
    count_dict = {gene: count for gene, count in zip(genes, counts) if count > 0}
    return count_dict

def counts_to_df(adata_cell):
    # get gene count dict
    count_dict = gene_count_dict(adata_cell)

    # get spatial coordinates
    x_coord = adata_cell.obsm['spatial'][0, 0]
    y_coord = adata_cell.obsm['spatial'][0, 1]
 
    # Get total number of molcs
    n_molcs = int(sum(count_dict.values()))

    # Build gene_list
    gene_list = [[key]*int(count_dict[key]) for i, key in enumerate(count_dict)]
    
    # build and return output df
    return pd.DataFrame({
        'x': [x_coord]*n_molcs,
        'y': [y_coord]*n_molcs,
        'gene': sum(gene_list,[]) # concatenate list of lists
    })

def baysor_export(adata, output_file, overwrite=False):
    if os.path.isfile(output_file) and overwrite:
        os.remove(output_file)
        print(f"Removed file `{output_file}`...")
    if os.path.isfile(output_file) and not overwrite:
        print(f"File `{output_file}` already found!")
        return

    #iterate through each cell in adata
    for j, cell in enumerate(adata):
        # get dataframe of counts
        df = counts_to_df(cell)
        
        # append counts to output file
        df.to_csv(output_file, mode='a', header=(j==0), index=False)
