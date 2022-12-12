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
  counts = "X", # where the counts are stored. Can also be the name of a `layer`
  data = NULL,
  layers = NULL, #additional `layer`s, that will be added as `Assay`s 
  project = "SeuratProject",
  meta.data = "obs",
  reductions = "obsm",
  graphs = "obsp", # set to NULL to ignore graphs
  verbose=T
){
  require(anndata)
  require(SeuratObject)
  # require(reticulate)
  
  # param check
  #TODO
  
  if(verbose){message(paste0("Loading ", adata_path, " into anndata..."))}
  adata <- read_h5ad(adata_path)
  if(verbose){message("Done!")}
  
  # Initialize seurat object
  if(verbose){message("Initializing Seurat object...")}
  #Note- anndata stores matrices as a dgTMatrix, where R/Seurat need a dgCMatrix. Need to figure out the best/fastest conversion strategy...
  if(counts=="X"){
    counts.mat <- t(adata[[counts]])
    
    seu <- CreateSeuratObject(
      assay = counts,
      project = project,
      counts = counts.mat,
      meta.data = adata[[meta.data]]
    )
  }else if(counts %in% names(adata$layers)){
    counts.mat <- t(ad$layers[[counts]])
    print(dim(counts.mat))
    
    if(!is.null(data)){
      message("I haven't implemented the conversion of normalized counts yet... Only raw counts will be converted.")
      # data.mat  <- t(ad$layers[[data]])
      
      seu <- CreateSeuratObject(
        assay = counts,
        counts = as.sparse(counts.mat),
        # data = as.sparse(data.mat),
        meta.data = adata[[meta.data]]
      )
    }else{
      seu <- CreateSeuratObject(
        assay = counts,
        counts = as.sparse(counts.mat),
        meta.data = adata[[meta.data]]
      )
    }
  }else{
    message("Counts not found, please check argument")
    return(NULL)
  }
  
  if(verbose){message("Done!")}
  
  #TODO adata$var
  
  
  # Add additional assays
  #TODO
  
  # Add reductions
  if(verbose){message("Adding reductions...")}
  if(!is.null(reductions)){
    for(reduction in names(adata$obsm)){
      if(verbose){message(paste0("adding ", reduction, "..."))}
      #PCA loadings- ad$varm$PCs
      
      #TODO adata$varm
      embedding = adata$obsm[[reduction]]
      rownames(embedding) <- Cells(seu)
      
      seu[[reduction]] <- CreateDimReducObject(
        embeddings = embedding,
        key = paste0(reduction,"_"),
        assay=counts
      )
    }
  }
  if(verbose){message("Done!")}
  
  # Add graphs
  if(!is.null(graphs)){
    for(graph in names(adata$obsp)){
      #TODO adata$varp
      seu[[graph]] <- as.Graph(
        adata$obsp[[graph]]
      )
    }
  }
  if(verbose){message("Done!")}
  
  # Return
  return(seu)
}
