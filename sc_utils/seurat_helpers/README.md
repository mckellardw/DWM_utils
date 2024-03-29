# `seurat_helpers` tips and code snippets

## Preprocessing


## Plotting
### Spatial plotting with `visListPlot()`
I tend to work with all of my Seurat objects in a `list()`, especially when working with Visium data. This keeps things tidy and also means I don't have to remember (or come up with...) variable names for every object individually. Another habit I have is to add the spatial dimensions for each spot/bead barcode into the Seurat object as a `reduction`. This allows me to use standard Seurat plotting functions like `DimPlot` (for categorical data) or `FeaturePlot()` (for continuous data). `visListPlot()` is a wrapper function that takes a list of Seurat objects and features, and generates a tidy grid of plots. I have added in features that I needed like maintaining the same color scale across all samples, and/or all features. Here are some examples of how to use it:

#### Plotting gene expression data:
```
tmp.feat <- c( # Features you want to plot
  "mmu-miR-1a-3p",
  "mmu-miR-206-3p"
)

visListPlot(
  skm.list,                             # List of 4 skeletal muscle samples (from our STRS paper!)
  sample.titles = c(                    # Vector of sample titles
    "Uninjured", "2dpi", "5dpi","7dpi"
  ),
  reduction="space",                    # The reduction in which the spatial locations are saved
  assay="kallisto",                     # Name of the `Assay()` in which gene expression data is sotred
  slot = 'data',                        # Plot log-normalized ("data") or raw counts ("counts"), following Seurat norms
  pt.size=0.6,                          # Size of each spot/bead
  legend.position = "right",
  font.size = 8,
  axis.title.angle.y=90,
  nrow = length(tmp.feat),              #
  flip_row_col = T,                     #
  combine = T,                          # Return either a grid of plots (T), or a list of columns of plots (F)
  verbose=F,
  colormap = "viridis",                 # Viridis color maps can be used automatically (also "plasma", "magma", and "inferno")
  colormap.direction = 1,               # Flip the color scale with `colormap.direction = -1`
  colormap.same.scale = T,              # Force all features to use the same color scale
  features=tmp.feat,                    # The features you want to plot - either gene names or metadata columns
  alt.titles = stringr::str_remove(     # Alternative feature names you want to use - here I am removing the "mmu-" from these miRNA names
    tmp.feat,
    pattern="mmu-"
  )
)&theme(                                # Use the `&` operator to apply some function (`theme()` here) to every plot!
  legend.text = element_text(
    size=6
  )
)&coord_fixed(                          # Use the `coord_fixed()` function to keep the spatial dimensions from being distorted
  ratio=1
)
```
