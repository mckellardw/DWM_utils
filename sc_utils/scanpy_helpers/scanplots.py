# Plotting functions for use with scanpy

# Knee plot to quality check UMI counts for single-cell data
def knee_plot(
    ADATA,
    x_lim=[0, 20000],
    line_width=2,
    line_color="b",
    verbose=False
):
    import matplotlib.pyplot as plt
    expected_num_cells = 10000

    for i in range(0,meta.shape[0]):
        knee = np.sort((np.array(ADATA.X.sum(axis=1))).flatten())[::-1]

        fig, ax = plt.subplots(figsize=(10, 7))

        ax.loglog(
            knee,
            range(len(knee)),
            linewidth=line_width,
            color=line_color
        )
    #     ax.axvline(x=knee[expected_num_cells], linewidth=3, color="k")
    #     ax.axhline(y=expected_num_cells, linewidth=3, color="k")

        ax.set_xlabel("UMI Counts")
        ax.set_ylabel("Set of Barcodes")
        ax.set_title(ADATA.obs["sample"][0])

        plt.xlim(x_lim)
        plt.grid(True, which="both")
        plt.show()

# Faceted plot for any embedding
# scanpy github issue reference- https://github.com/scverse/scanpy/issues/955
def facet_embedding(
    adata, 
    clust_key, 
    basis, 
    size=60, 
    frameon=False, 
    legend_loc=None, 
    **kwargs
):
    # import scanpy as sc

    tmp = adata.copy()

    for i,clust in enumerate(adata.obs[clust_key].cat.categories):
        tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
        tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]

    sc.pl.embedding(
        tmp, 
        groups=tmp.obs[clust].cat.categories[1:].values, 
        color=adata.obs[clust_key].cat.categories.tolist(), 
        basis=basis,
        size=size, frameon=frameon, legend_loc=legend_loc, 
        **kwargs
    )