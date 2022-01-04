###############################
# seutils -  Seurat utility functions
# Written by: David McKellar
# version: 1.0
###############################

########################################
## Helpers for grabbing/searching feature names
########################################

# Quickly return genes/feature names from a Seurat object
Features <- function(
  SEU,
  assay=NULL
){
  if(is.null(assay)){
    assay <- SEU@active.assay
  }

  return(rownames(GetAssayData(SEU,assay=assay)))
}

# Check gene names for a pattern, using grep
grepGenes <- function(
  SEU,
  pattern=NULL, # pattern to look for
  assay=NULL,
  filter.pattern=NULL, # pattern to remove
  sort.by = c("expression","abc"),
  verbose=T
){
  if(is.null(pattern)){
    if(verbose){message("Need a pattern to look for!")}
    return(NULL)
  }
  if(is.list(SEU)){
    if(verbose){message("Don't pass a list!")}
    return(NULL)
  }
  if(is.null(assay)){
    assay <- SEU@active.assay
  }

  genes = SEU@assays[[assay]]@counts@Dimnames[[1]]

  

  if(length(pattern)>1){
    if(verbose){
      message(paste0("Found ", length(genes), " features in the assay '",assay,"'..."))
      message(paste0("Looking for multiple patterns in these features..."))
    }
    out.genes <- lapply(
      pattern,
      FUN = function(PAT) genes[grep(pattern=PAT,x=genes)] #get genes
    )
    out.genes <- unlist(out.genes)
  }else{
    if(verbose){
      message(paste0("Found ", length(genes), " features in the assay '",assay,"'..."))
      message(paste0("Looking for '", pattern, "' in these features..."))
    }
    out.genes <- genes[grep(pattern=pattern,x=genes)] #get genes
  }
  
  if(length(out.genes)==0){
    message("Nothing found!\n")
    return(NULL)
  }

  if(!is.null(filter.pattern)){ # filter out filter.pattern
    if(verbose){message(paste0("Removing features containing '", filter.pattern,"' from output..."))}
    out.genes <- out.genes[!grepl(pattern=filter.pattern,x=out.genes)]
  }

  if(is.null(sort.by[1])){
    if(verbose){message("Output genes not sorted...")}
  }else if(length(out.genes)==1){
    # Do nothing
  }else if(sort.by[1]=="expression"){
    if(verbose){message("Output genes sorted by expression...")}

    out.genes <- GetAssayData(
      SEU,
      assay=assay,
      slot="counts"
    )[out.genes,] %>%
        rowSums() %>% # sum counts across all cells
        sort(decreasing=T) %>% # sort genes according to total expression
        names() # get gene names back
  }else if(sort.by[1]=="abc"){
    if(verbose){message("Output genes sorted alphabetically...")}

    out.genes <- sort(out.genes, decreasing=T)
  }

  # Return matching gene names!
  return(
    out.genes
  )
}

# Convert ensembl IDs to gene IDs using a biomaRt reference
ens2gene <- function(
  ens=NULL, # vector of ensembl IDs to convert
  biomart.info=NULL, # biomaRt database
  ens.colname="ensembl_gene_id",
  gene.colname="mgi_symbol",
  verbose=F
){
  if(is.null(ens)){
    message("Need ensembl IDs to convert!")
    return(NULL)
  }
  if(is.null(biomart.info)){
    message("Need biomaRt reference for conversion!")
    return(NULL)
  }
  if(!ens.colname %in% colnames(biomart.info)){
    message("ensembl ID column not found in biomaRt reference. Check input for 'ens.colname'")
    return(NULL)
  }
  if(!gene.colname %in% colnames(biomart.info)){
    message("Gene ID column not found in biomaRt reference. Check input for 'gene.colname'")
    return(NULL)
  }
  print("Haven't finished this yet...")
  # genes <- biomart.info[,gene.colname]
  return(ens)
}

########################################
## General Seurat workflow helpers
########################################
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


# Collapse cell/nuclei/spot counts for multimapped genes
#TODO: parallelize!
collapseMultimappers <- function(
  SEU,
  assay=NULL,
  new.assay.name=NULL,
  verbose=F
){

  if(is.null(new.assay.name)){
    new.assay.name = paste0(assay,"_collpased")
    message("Using ",new.assay.name, " as new.assay.name...")
  }
  if(is.null(new.assay.name)){
    message("Need new.assay.name!")
    return(SEU)
  }

  SEU@active.assay <- assay

  multi.feats <- grepGenes( #Find genes with a period in them
    SEU, 
    assay = assay, 
    pattern="\\.", 
    sort.by="abc",
    verbose=verbose
  ) 
  if(length(multi.feats)==0){
    message("No multimappers found!")
    return(SEU)
  }

  multi.patterns <- stringr::str_split( #extract actual gene names
    multi.feats, 
    pattern = "\\.",
    n = 2
  ) %>% 
    lapply(FUN=function(X) X[1]) %>%
    unlist() %>%
    unique()

  if(verbose){
    message(paste0("Found ", length(multi.patterns), " multimappers and ", length(multi.feats)," loci..."))
  }

  # Collapse counts for each gene
  mat.multi <- GetAssayData(SEU, assay=assay, slot="counts") # count mat

  collapsed.list <- lapply(
    multi.patterns,
    FUN=function(X){
      tmp.genes = rownames(mat.multi)[grep(rownames(mat.multi),pattern=X)]
      tmp.mat = mat.multi[tmp.genes,]

      if(length(tmp.genes)==1){
        return(tmp.mat)
      }else{
        return(colSums(tmp.mat))
      }
    }
  )
  collapsed.mat <- do.call(rbind, collapsed.list) %>% as.sparse()
  rownames(collapsed.mat) <- multi.patterns

  print(multi.patterns)

  # Add new assay with collapsed counts + the rest of the genes
  if(verbose){cat(paste0("Adding back ", nrow(collapsed.mat), " features...\n"))}

  solo.feats <- rownames(SEU)[!rownames(SEU)%in%multi.feats]

  out.mat <- rbind(
    GetAssayData(SEU,assay=assay, slot="counts")[solo.feats,],
    collapsed.mat
  )
  SEU[[new.assay.name]] <- CreateAssayObject(counts=out.mat)

  SEU@active.assay <- new.assay.name

  # Return Seurat object!
  return(SEU)
}


# Preprocessing wrapper function
#   (1) NormalizeData(SEU) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
#   (2) FindNeighbors %>% RunUMAP, FindClusters
seuPreProcess <- function(
  SEU=NULL,
  assay='RNA',
  n.pcs=50,
  res=0.8,
  verbose=F
){
  if(is.null(SEU)){
    cat("Need a Seurat object to preprocess!\n")
    return(NULL)
  }
  if(!assay %in% Assays(SEU)){
    cat(paste0(assay, " not found in the seurat object! Not preprocessed.\n"))
    return(SEU)
  }

  # NormalizeData(SEU) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  pca.name = paste0('pca_', assay)
  pca.key = paste0(pca.name,'_')
  umap.name = paste0('umap_', assay)

  SEU = NormalizeData(
    SEU
  ) %>% FindVariableFeatures(
    assay = assay,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = verbose
  ) %>% ScaleData(
    assay = assay
  ) %>% RunPCA(
    assay = assay,
    reduction.name = pca.name,
    reduction.key = pca.key,
    verbose = verbose,
    npcs = n.pcs
  )

  #find pcs to use
  n.pcs.use = npcs(SEU=SEU, var.total = 0.95, reduction = pca.name)

  # FindNeighbors %>% RunUMAP, FindClusters
  SEU <- FindNeighbors(
    SEU,
    reduction = pca.name,
    dims = 1:n.pcs.use,
    force.recalc = TRUE,
    verbose = verbose
  ) %>% RunUMAP(
    reduction = pca.name,
    dims = 1:n.pcs.use,
    verbose = verbose,
    reduction.name=umap.name
  )

  SEU@reductions[[umap.name]]@misc$n.pcs.used <- n.pcs.use

  SEU <- FindClusters(
    object = SEU,
    resolution = res,
    verbose = verbose
  )
  # SEU[[paste0('RNA_res.',res)]] <- as.numeric(SEU@active.ident)
  gc()
  return(
    tryCatch(
      SEU,
      error=function(e) NULL
    )
  )
}


# Add a new ident seurat metadata filed) based on a list of cell types
#
#     object:     seurat object
#     old.idents: name of the idents metadata you will be assigning cell types to
#     new.idents: vector of cell types, in order of cluster number
#     newName:    string of the new idents name
#
AddCellTypeIdents <- function(
  SEU=NULL, old.name, new.name=NULL, new.idents, verbose=FALSE
){
  old.idents = as.vector(names(table(SEU[[old.name]])))

  if(is.null(new.name)){
    cat("**Need a new.name for the new idents**\n")
  }else{
    SEU[[new.name]] <- as.vector(SEU[[old.name]])
    for(i in 1:length(old.idents)){
      if(verbose){cat("Adding ", new.idents[i],"...", sep = "")}
      SEU[[new.name]][ SEU[[new.name]]==old.idents[i] ] <- new.idents[i]
      if(verbose){cat("Done!\n", sep = "")}
    }
  }

   return(SEU)

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
seu_entropy <- function(
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

########################################
## biomaRt helper functions
########################################

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
  genes.list = biomaRt::getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = x ,
    mart = human,
    attributesL = c("mgi_symbol"),
    martL = mouse,
    uniqueRows = TRUE
  )

  # Get unique names of genes (in case gene names are duplicated)
  mouse.gene.list <- unique(genes.list[, 2])

  return(mouse.gene.list)
}

# Add gene biotype percentages to a seurat object, given a biomaRt object.
seu_biotypes <- function(
  SEU,
  biomart=NULL, # biomaRt or data.frame containing gene biotypes
  gene.colname,
  biotype.colname,
  add.as=c("metadata","assay"), # how percent features should be added
  assay="RNA",
  verbose=TRUE
){

  if(is.null(biomart)){
    message("Need a list of gene biotypes! Nothing done.")
    return(SEU)
  }

  if(add.as[1]=="assay"){
    message("add.as='assay' is not yet implemented")
    new.assay.name = paste0(assay,".biotypes")
    return(SEU)
  }

  if(verbose){message(paste0("Adding gene biotype percentage values as ", add.as, " ..."))}

  biotypes = unique(biomart[[biotype.colname]])
  for(biotype in biotypes){
    tmp.mart =  biomart[biomart[[biotype.colname]]==biotype,] #subset out current gene biotype
    tmp.feat = unique( #get unique gene names which are present in the SEU object
      tmp.mart[[gene.colname]][tmp.mart[[gene.colname]]%in%Features(SEU, assay=assay)]
    )

    if(length(tmp.feat)==0){
      if(verbose){message(paste0("  No ", biotype, " genes found..."))}
    }else{
      if(add.as[1]=="metadata"){
        SEU <- PercentageFeatureSet(
          SEU,
          col.name = paste0("pct.",biotype),
          assay = assay,
          features=tmp.feat
        )
      }
      if(add.as[1]=="assay"){

      }
    }
  }
  return(SEU)
}

########################################
## Other helpers...
########################################

# Run PHATE on reduction, with a Seurat objects
seuPHATE <- function(
  SEU=NULL,
  reduction="pca",
  ndims=50,
  reduction.name=NULL, # output reduction name
  reduction.key=NULL, # output reduction key
  # phate() defaults
  ndim = 2,
  knn = 5,
  decay = 40,
  n.landmark = 2000,
  gamma = 1,
  t = "auto",
  mds.solver = "sgd",
  knn.dist.method = "euclidean",
  knn.max = NULL,
  init = NULL,
  mds.method = "metric",
  mds.dist.method = "euclidean",
  t.max = 100,
  npca = 100,
  plot.optimal.t = FALSE,
  verbose = 1,
  n.jobs = 1,
  seed = NULL,
  potential.method = NULL,
  k = NULL,
  alpha = NULL,
  use.alpha = NULL
){
  require(Seurat)
  require(phateR)

  dims <- SEU@reductions[[reduction]]@stdev %>% order(decreasing = T)
  dims <- dims[1:ndims]

  # run PHATE
  tmp.phate <- phate(
    SEU@reductions[[reduction]]@cell.embeddings[,dims], #cells x reduc dims
    ndim = ndim,
    knn = knn,
    decay = decay,
    n.landmark = n.landmark,
    gamma = gamma,
    t = t,
    mds.solver = mds.solver,
    knn.dist.method = knn.dist.method,
    knn.max = knn.max,
    init = init,
    mds.method = mds.method,
    mds.dist.method = mds.dist.method,
    t.max = t.max,
    npca = npca,
    plot.optimal.t = plot.optimal.t,
    verbose = verbose,
    n.jobs = n.jobs,
    seed = seed,
    potential.method = potential.method,
    k = k,
    alpha = alpha,
    use.alpha = use.alpha
  )

  # Set reduction name & key for SEU
  if(is.null(reduction.name)){
    reduction.name <- paste0("phate_",reduction)
  }
  if(is.null(reduction.key)){
    reduction.key <- paste0("phate_",reduction,"_")
  }

  colnames(tmp.phate$embedding) <- paste0(reduction.key, 1:2)
  tmp.phate$params$data <- NULL

  SEU[[reduction.name]] <- CreateDimReducObject(
    embeddings = tmp.phate$embedding,
    key = reduction.key,
    assay = 'RNA',
    misc = tmp.phate$params
  )

  # Add std dev to reduction
  SEU@reductions[[reduction.name]]@stdev <-
    apply(SEU@reductions[[reduction.name]]@cell.embeddings, 2, sd) #find std dev for phate vals

  return(SEU)
}
