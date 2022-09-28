# Utility functions for use with scanpy

# Calculate the number of PCs that contain some proportion (default is 95%) of the variance
def npcs(
  ADATA,
  var_perc=0.95,
  reduction="pca"
):
    import numpy as np
    get_var = lambda i: np.var(ADATA.obsm[reduction][:,i])

    if ADATA.obsm[reduction] is None:
        print(f"Reduction {reduction}, not found!")
        return None
    else:
        var_tmp = [get_var(i) for i in list(range(0,ADATA.obsm[reduction].shape[1]))]
        var_cut = var_perc * np.sum(var_tmp)
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
    get_var = lambda i: np.var(ADATA.obsm[reduction][:,i])
    var_tmp = [get_var(i) for i in list(range(0,ADATA.obsm[reduction].shape[1]))]
    print("Harmony PC variance:")
    print(var_tmp)

    pc_order = np.argsort(var_tmp)[::-1]
    ADATA.obsm[reduction] = ADATA.obsm[reduction][:,pc_order]
