# `scanpy_helpers` tips and code snippets

## Preprocessing
#TODO

## Plotting
### DotPlots
I prefer to keep my genes organized as a dict, where keys refer to the corresponding cell types, biological process (i.e. GO terms), etc. When plotting with these, `scanpy` annotates the DotPlot with the key, which I think is very helpful in adding cell type annotations to unsupervised clustering results.

This snippet will subset your giant dict (`genes_dict` here) of all cell types/etc. so that your DotPlot is actually legible:
```
{key:genes_dict[key] for key in cell_types}
```
This can be passed directly to `var_names` in `sc.pl.dotplot`:
```
cell_types=['Fibroblasts','Macrophages']

sc.pl.dotplot(
    adata,
    var_names={key:genes_dict[key] for key in cell_types},
    groupby="leiden_clusters",
    color_map="plasma",
    swap_axes=True,
    standard_scale ='var',
    dendrogram=True
)
```
