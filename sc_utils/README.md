# seurat_helpers/
#### Wrapper/helper functions written to smooth out Seurat pipelines... See [this wonderful website](https://satijalab.org/seurat/index.html) from the Satija Lab to learn more about `Seurat`!

## **seutils.R** (seurat utils)
`Seurat` utility functions, to make your life easier.

#### `Features()`
Grab feature names from any `Assay` within a `Seurat` object - analogous to `Seurat::Cells()`

#### `grepGene()`
Quick function to look for genes in a `Seurat` object that match a certain pattern. Can also set a `filter.pattern` to exclude genes which contain the string. Built on `grep` and `Seurat`

#### `ens2gene()`
Convert a list of ensembl IDs to gene symbols. Compatible with outputs from `biomaRt`.

#### `npcs()`
Calculate how many principal components account for X% of the variance.

#### `collapseMultimappers()`
Quickly collapse multimapper genes in a `Seurat` assay (non-unique features, marked by a "." in the feature name)

#### `seuPreProcess()`
Generic single-cell pipeline, based on [this vignette](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) from the Satija Lab

#### `AddCellTypeIdents()`


## **seuplots.R** (seurat plots)
Additional plots not included in Seurat including

#### `visListPlot()`
Wrapper function for neatly plotting multiple Visium objects and multiple features into a grid. Input is a `list()` of `Seurat` objects. Includes automatic color scaling and lots of options for customization. Built on `patchwork`. Combine with `coord_fixed()` to maintain true spatial coordinates.
```
visListPlot <- function(
  seu.list,
  features=NULL,
  alt.titles=NULL, # alternative titles for genes/features being passed
  sample.titles=NULL, # Sample IDs (y-axis labels)
  assay='Spatial',
  reduction="space",
  slot="data",
  legend.position="bottom",
  pt.size=1,
  font.size=8,
  axis.title.angle.y=90, #y-axis title angle (parameter for ggplot::element_text)
  combine=TRUE,
  abs.heights=TRUE, # Use absolute heights to size each subplot
  nrow=NULL,
  ncol=NULL,
  colormap="viridis", # either a viridis option, a vector of colors, or a list of options corresponding to `features`
  colormap.direction=1,
  colormap.same.scale=F, #whether (T) or not (F) to set all features to the same colormap scale
  na.value=gray(0.69), # color for na.value (spot where gene is not detected)
  verbose=FALSE
)
```

#### `visCoMap()`
Plot co-expression of any two genes/features, where "co-expression" can be abstracted to an generic function (default is just `prod`, to see expression of Gene_A times expression of Gene_B). Built on top of visListPlot(), so it also takes a list of `Seurat` objects.
```
visCoMap <- function(
  seu.list,
  features=NULL, # pair of genes
  alt.titles=NULL, # alternative titles for genes/features being passed
  sample.titles=NULL, # Sample IDs (y-axis labels)
  assay='Spatial', # Either a single assay, or a pair of assays (useful for spliced vs. unspliced plotting)
  reduction="space",
  slot="data", # Either a single slot, or a pair of slots
  legend.position="bottom",
  pt.size=1,
  font.size=8,
  axis.title.angle.y=90, #y-axis title angle (parameter for ggplot::element_text)
  combine=TRUE,
  abs.heights=TRUE, # Use absolute heights to size each subplot
  nrow=NULL,
  ncol=NULL,
  colormap=NULL, # either a viridis option, a vector of colors, or a list of options corresponding to `features`
  colormap.direction=1,
  colormap.same.scale=F, #whether (T) or not (F) to set all features to the same colormap scale
  na.value=gray(0.69), # color for na.value (spot where gene is not detected)
  comap.fxn = prod,
  coex.name = NULL, # plot title for computed co-expression values
  verbose=FALSE
)
```
#TODO
- Volcano plots (ggvolcano_v2)
- Silhouette plot

## `DWM_DoubletFinder_v1.R`
Has not been used in a while, *user beware*.

Partially re-written `DoubletFinder` functions with added functionalities like naming the output metadata columns, `DoubletFinder` in integrated reduced dimensional space, etc.
