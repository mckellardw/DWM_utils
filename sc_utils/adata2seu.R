#' @title adata2seu
#' 
#' @description Convert an r-anndata object (which can be preprocessed in python-anndata) to a SeuratObject
#' 
#' @param adata_path
#' @param counts
#' @param meta.data 
#' @param reductions
#' @param graphs
#' @param
#' @param
#' @param verbose 
adata2seu <- function(
  adata_path=NULL,
  counts = "X",
  data = NULL,
  project = "SeuratProject",
  meta.data = "obs",
  reductions = "obsm",
  graphs = "obsp",
  verbose=T
){
  require(anndata)
  require(SeuratObject)
  # require(reticulate)
  
  if(verbose){message(paste0("Loading ", adata_path, " into anndata..."))}
  adata <- read_h5ad(adata_path)
  if(verbose){message("Done!")}
  
  # Initialize seurat object
  if(verbose){message("Initializing Seurat object...")}
  if(counts=="X"){
    seu <- CreateSeuratObject(
      assay = counts,
      project = project,
      counts = t(adata[[counts]]),
      meta.data = adata[[meta.data]]
    )
  }else{
    if(counts %in% adata$layers)
    seu <- CreateSeuratObject(
      assay = counts,
      counts = as.sparse(t(adata$layers[[counts]])),
      meta.data = adata[[meta.data]]
    )
  }
  
  if(verbose){message("Done!")}
  
  #TODO adata$var
  
  
  # Add additional assays
  #TODO
  
  # Add reductions
  if(verbose){message("Adding reductions...")}
  for(reduction in names(adata$obsm)){
    if(verbose){message(paste0("adding ", reduction, "..."))}
    #PCA loadings- ad$varm$PCs
    
    #TODO adata$varm
    
    seu[[reduction]] = CreateDimReducObject(
      embeddings = adata$obsm[[reduction]],
      key = paste0(reduction,"_")
    )
  }
  if(verbose){message("Done!")}
  
  # Add graphs
  for(graph in names(adata$obsp)){
    #TODO adata$varp
    seu[[graph]] = as.Graph(
      adata$obsp[[graph]]
    )
  }
  if(verbose){message("Done!")}
  
  # Return
  return(seu)
}
