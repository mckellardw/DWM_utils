###############################
# seutils -  Seurat utility functions
# Written by: David McKellar
# version: 1.0
###############################

# Calculate the number of PCs that contain some proportion (default is 95%) of the variance
npcs <- function(
  SEU,
  var.total=0.95,
  reduction="pca"
){
  if(is.null(SEU@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }

  tmp.var <- (SEU@reductions[[reduction]]@stdev)^2
  var.cut <- var.total*sum(tmp.var)
  n.pcs=0
  var.sum = 0
  while(var.sum < var.cut){
    n.pcs = n.pcs + 1
    var.sum <- var.sum + tmp.var[n.pcs]
  }

  return(n.pcs)
}


# Borrowed/adapted from the Marioni Lab, DropletUtils package (https://rdrr.io/github/MarioniLab/DropletUtils/src/R/write10xCounts.R)
#   (Had R version issues getting it to work as a dependency)
#' @importFrom utils write.table
#' @importFrom Matrix writeMM
#' @importFrom R.utils gzip
write_sparse <- function(
  path, # name of new directory
  x, # matrix to write as sparse
  barcodes, # cell IDs, colnames
  features # gene IDs, rownames

  # gene.symbol,#not used
  # gene.type
  ){
  require(utils,quietly = T)
  require(Matrix,quietly = T)
  require(R.utils,quietly = T)

  if(!dir.exists(path)){
    dir.create(path, showWarnings=FALSE)
  }

  # gene.info <- data.frame(gene.id, gene.symbol, stringsAsFactors=FALSE)

  # gene.info$gene.type <- rep(gene.type, length.out=nrow(gene.info))
  mhandle <- file.path(path, "matrix.mtx")
  bhandle <- gzfile(file.path(path, "barcodes.tsv.gz"), open="wb")
  fhandle <- gzfile(file.path(path, "features.tsv.gz"), open="wb")
  on.exit({
    close(bhandle)
    close(fhandle)
  })

  writeMM(x, file=mhandle)
  write(barcodes, file=bhandle)
  write.table(features, file=fhandle, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

  # Annoyingly, writeMM doesn't take connection objects.
  gzip(mhandle)

  return(NULL)
}


# TODO- haven't used this in a while, make sure it works...
# calculate entropy across groups in a seurat object
#     output: returns data.frame ("group.by" rows by 1 col)
sc_entropy <- function(
  SEU,
  group.by="sample",
  entropy.on="factorIDs",
  out.name=NULL,
  weighted=T,
  norm2one=T,
  verbose=T
){
  group.levels <- unlist(unique(SEU[[group.by]]))

  sub.meta <- SEU@meta.data[,c(group.by,entropy.on)]

  entropy.out <- list()
  for(lev in group.levels){
    if(verbose){cat("calculating entropy for ", lev, "... ",sep = "")}
    tmp <- sub.meta[sub.meta[[group.by]]==lev,] # subset metadata by sample

    entro.levels <- unlist(unique(sub.meta[[entropy.on]]))

    perc <- table(tmp[[entropy.on]])/nrow(tmp) # find proportions for each entropy level
    if(weighted){
      w <- (table(sub.meta[[entropy.on]])/nrow(sub.meta))[perc!=0] # weights for present levels
    }else{
      w <- rep(1,length(table(sub.meta[[entropy.on]])))[perc!=0]
    }

    perc <- perc[perc!=0] # remove zeroes

    entropy.out[[lev]] <- sum(-1*perc*log2(perc)*w) #calculate entropy

    if(verbose){cat("Done!\n",sep = "")}
  }

  # re-format entropy output
  if(is.null(out.name)){
    out.name <- paste0("entropy_",group.by,".by.",entropy.on)
  }
  entropy.out <- t(as.data.frame(entropy.out))
  colnames(entropy.out) <- out.name

  #replace NaN's with zeroes
  # entropy.out[is.nan(entropy.out)] <- 0

  if(norm2one){
    entropy.out <- entropy.out/log2(length(unique(sub.meta[[entropy.on]])))
  }

  return(entropy.out)
}


# Calculate silhouette
seu_silhouette <- function(
  SEU,
  group.by,
  reduction,
  #TODO- add graph input as option
  meta.name.out=NULL,
  dims=1:10
){
  require(cluster)

  if(is.null(meta.name.out)){
    meta.name.out=paste0('sil.',reduction,'.',group.by)
  }
  if(is.null(group.by)){
    message("group.by is missing!")
    return(SEU)
  }
  if(is.null(reduction)){
    message("reduction is missing!")
    return(SEU)
  }

  # Calculate silhouette coefficient
  sil.out <- silhouette(
    x = as.numeric(x = as.factor(x = unlist(SEU[[group.by]]))),
    dist = dist(x = Embeddings(object = SEU, reduction=reduction)[,dims])
  )
  SEU[[meta.name.out]] <- sil.out[,3]

  return(SEU)
}


# Basic function to convert human to mouse gene names
#   From @leonfodoulian (https://github.com/satijalab/seurat/issues/462)
#   https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
#
#   x = list of genes to be converted
#   Usage: mouse.genes <- lapply(X = human.genes, ConvertHumanGeneListToMM)
#
ConvertHumanGeneListToMM <- function(x){
  require(biomaRt)

  # Load human ensembl attributes
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  # Load mouse ensembl attributes
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  # Link both datasets and retrieve mouse genes from the human genes
  genes.list = biomaRt::getLDS(attributes = c("hgnc_symbol"),
                               filters = "hgnc_symbol",
                               values = x ,
                               mart = human,
                               attributesL = c("mgi_symbol"),
                               martL = mouse,
                               uniqueRows = TRUE)

  # Get unique names of genes (in case gene names are duplicated)
  mouse.gene.list <- unique(genes.list[, 2])

  # # Print the first 6 genes found to the screen
  # print(head(mouse.gene.list))
  return(mouse.gene.list)
}


# Add a new ident seurat metadata filed) based on a list of cell types
#
#     object:     seurat object
#     old.idents: name of the idents metadata you will be assigning cell types to
#     new.idents: vector of cell types, in order of cluster number
#     newName:    string of the new idents name
#
AddCellTypeIdents <- function(seu=NULL, old.name, new.name=NULL, new.idents, verbose=FALSE){
  old.idents = as.vector(names(table(seu[[old.name]])))

  if(is.null(new.name)){
    cat("**Need a new.name for the new idents**\n")
  }else{
    seu[[new.name]] <- as.vector(seu[[old.name]])
    for(i in 1:length(old.idents)){
      if(verbose){cat("Adding ", new.idents[i],"...", sep = "")}
      seu[[new.name]][ seu[[new.name]]==old.idents[i] ] <- new.idents[i]
      if(verbose){cat("Done!\n", sep = "")}
    }
  }

   return(seu)

}
