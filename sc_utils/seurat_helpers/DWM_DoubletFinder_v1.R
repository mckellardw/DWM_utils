# Custom functions for DoubletFinder, adapted from cmcginnis' GitHub repo
## Source: https://github.com/chris-mcginnis-ucsf/DoubletFinder


# Helper Functions  ####
#     EstimateDoubletRate ####
# Simple linear regression to spit out a list of estimated doublet rates, based on 10x's published rates
estimateDoubletRate.DWM <- function(
    seu.list,
    doublet.dist=NULL
){
  if(is.null(doublet.dist)){
    doublet.dist <- data.frame(
      cells.recovered=c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000),
      multi.rate=     c(0.4,  0.8,  1.6,  2.3,  3.1,  3.9,  4.6,  5.4,  6.1,  6.9,   7.6)
    )
  }

  fit <- lm(multi.rate~cells.recovered, data = doublet.dist)
  fit.func <- function(x){
    return(as.numeric(fit$coefficients['cells.recovered']*x + fit$coefficients['(Intercept)']))
  }

  ncells.list <- lapply(seu.list, ncol)
  out.list <- lapply(ncells.list, fit.func)
  return(unlist(out.list))

}

# Modified DoubletFinder functions  ####
#     parallel_paramSweep_V3 ####
parallel_paramSweep_v3_DWM <- function(
  n,
  n_real.cells,
  real.cells,
  pK,
  pN,
  data,
  orig.commands,
  PCs,
  assay='RNA',
  slot='counts',
  sct
){
  #TODO add param checks

    sweep.res.list = list()
  list.ind = 0

  ## Make merged real-artifical data
  cat(paste("Creating artificial doublets for pN = ", pN[n]*100,"%\n",sep=""))
  n_doublets <- round(n_real.cells/(1 - pN[n]) - n_real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)

  ## Initialize seurat object
  if(slot=='counts'){
    cat("Creating Seurat object with artificial doublets...\n")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
    seu_wdoublets <- NormalizeData(
      seu_wdoublets,
      normalization.method = orig.commands[[paste0("NormalizeData.",assay)]]@params$normalization.method,
      scale.factor = orig.commands[[paste0("NormalizeData.",assay)]]@params$scale.factor,
      margin = orig.commands[[paste0("NormalizeData.",assay)]]@params$margin
    )
  }else if(slot=='data'){
    cat("Creating Seurat object with artificial doublets...\n")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets) #don't renormalize if normalized data is provided
  }

  ## Pre-process Seurat object
  if(assay=='sct'){
    require(sctransform)
    if(slot=='counts'){
      cat('WARNING: running SCTransfrom on normalized data! \n')
    }

    cat("Running SCTransform...\n")
    seu_wdoublets <- SCTransform(seu_wdoublets)

    cat("Running PCA...\n")
    seu_wdoublets <- RunPCA(
      seu_wdoublets,
      npcs = length(PCs),
      reduction.name='DOUBLETFINDER_PCA'
    )
    pca.coord <- seu_wdoublets@reductions$DOUBLETFINDER_PCA@cell.embeddings[ , PCs]
    cell.names <- rownames(seu_wdoublets@meta.data)
    nCells <- length(cell.names)
    rm(seu_wdoublets); gc()
  }else if(assay=='RNA'){
    cat("     Piping, FindVariableGenes(), FindVariableFeatures(), ScaleData(), and RunPCA()...\n")
    seu_wdoublets <-  FindVariableFeatures(
      seu_wdoublets,
      selection.method = orig.commands[[paste0("FindVariableFeatures.",assay)]]$selection.method,
      loess.span = orig.commands[[paste0("FindVariableFeatures.",assay)]]$loess.span,
      clip.max = orig.commands[[paste0("FindVariableFeatures.",assay)]]$clip.max,
      mean.function = orig.commands[[paste0("FindVariableFeatures.",assay)]]$mean.function,
      dispersion.function = orig.commands[[paste0("FindVariableFeatures.",assay)]]$dispersion.function,
      num.bin = orig.commands[[paste0("FindVariableFeatures.",assay)]]$num.bin,
      binning.method = orig.commands[[paste0("FindVariableFeatures.",assay)]]$binning.method,
      nfeatures = orig.commands[[paste0("FindVariableFeatures.",assay)]]$nfeatures,
      mean.cutoff = orig.commands[[paste0("FindVariableFeatures.",assay)]]$mean.cutoff,
      dispersion.cutoff = orig.commands[[paste0("FindVariableFeatures.",assay)]]$dispersion.cutoff
    ) %>% ScaleData(
      features = orig.commands[[paste0("ScaleData.",assay)]]$features,
      model.use = orig.commands[[paste0("ScaleData.",assay)]]$model.use,
      do.scale = orig.commands[[paste0("ScaleData.",assay)]]$do.scale,
      do.center = orig.commands[[paste0("ScaleData.",assay)]]$do.center,
      scale.max = orig.commands[[paste0("ScaleData.",assay)]]$scale.max,
      block.size = orig.commands[[paste0("ScaleData.",assay)]]$block.size,
      min.cells.to.block = orig.commands[[paste0("ScaleData.",assay)]]$min.cells.to.block
    )%>% RunPCA(
      features = orig.commands[[paste0("RunPCA.",assay)]]$features,
      npcs = length(PCs),
      rev.pca =  orig.commands[[paste0("RunPCA.",assay)]]$rev.pca,
      weight.by.var = orig.commands[[paste0("RunPCA.",assay)]]$weight.by.var,
      reduction.name='DOUBLETFINDER_PCA',
      verbose=FALSE
    )

    pca.coord <- seu_wdoublets@reductions$DOUBLETFINDER_PCA@cell.embeddings[ , PCs]
    cell.names <- rownames(seu_wdoublets@meta.data)
    nCells <- length(cell.names)
    rm(seu_wdoublets); gc() # Free up memory

  }else{
    cat("     Piping FindVariableGenes(), FindVariableFeatures(), ScaleData(), and RunPCA()...\n")
    seu_wdoublets <- FindVariableFeatures(
      seu_wdoublets,
      selection.method = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$selection.method,
      loess.span = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$loess.span,
      clip.max = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$clip.max,
      mean.function = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$mean.function,
      dispersion.function = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$dispersion.function,
      num.bin = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$num.bin,
      binning.method = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$binning.method,
      nfeatures = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$nfeatures,
      mean.cutoff = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$mean.cutoff,
      dispersion.cutoff = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$dispersion.cutoff
    ) %>% ScaleData(
      features = orig.commands[[paste('ScaleData', assay, sep='.')]]$features,
      model.use = orig.commands[[paste('ScaleData', assay, sep='.')]]$model.use,
      do.scale = orig.commands[[paste('ScaleData', assay, sep='.')]]$do.scale,
      do.center = orig.commands[[paste('ScaleData', assay, sep='.')]]$do.center,
      scale.max = orig.commands[[paste('ScaleData', assay, sep='.')]]$scale.max,
      block.size = orig.commands[[paste('ScaleData', assay, sep='.')]]$block.size,
      min.cells.to.block = orig.commands[[paste('ScaleData', assay, sep='.')]]$min.cells.to.block
    )%>% RunPCA(
      features = orig.commands[[paste('ScaleData', assay, sep='.')]]$features,
      npcs = length(PCs),
      rev.pca =  orig.commands[[paste('RunPCA', assay, sep='.')]]$rev.pca,
      weight.by.var = orig.commands[[paste('RunPCA', assay, sep='.')]]$weight.by.var,
      reduction.name='DOUBLETFINDER_PCA',
      verbose=FALSE
    )

    pca.coord <- seu_wdoublets@reductions$DOUBLETFINDER_PCA@cell.embeddings[ , PCs]
    cell.names <- rownames(seu_wdoublets@meta.data)
    nCells <- length(cell.names)
    rm(seu_wdoublets); gc() # Free up memory
  }

  ## Compute PC distance matrix
  cat("Calculating PC distance matrix...\n")
  dist.mat <- fields::rdist(pca.coord)[,1:n_real.cells]

  ## Pre-order PC distance matrix prior to iterating across pK for pANN computations
  print("Defining neighborhoods...")
  for (i in 1:n_real.cells) {
    dist.mat[,i] <- order(dist.mat[,i])
  }

  ## Trim PC distance matrix for faster manipulations
  ind <- round(nCells * max(pK))+5
  dist.mat <- dist.mat[1:ind, ]

  ## Compute pANN across pK sweep
  cat("Computing pANN across all pK...\n")
  for (k in 1:length(pK)) {
    print(paste("pK = ", pK[k], "...", sep = ""))
    pk.temp <- round(nCells * pK[k])
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    colnames(pANN) <- "pANN"
    rownames(pANN) <- real.cells
    list.ind <- list.ind + 1

    for (i in 1:n_real.cells) {
      neighbors <- dist.mat[2:(pk.temp + 1),i]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/pk.temp
    }

    sweep.res.list[[list.ind]] <- pANN
  }

  return(sweep.res.list)
}

#     parallel_paramSweep_PCA ####
parallel_paramSweep_PCA <- function(
    n, n_real.cells, real.cells, pK, pN, pc_data, PCs
){

  cat("Creating artificial doublets for pN = ", pN[n]*100, "\n")

  ## Make merged real-artifical data
  n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (pc_data[, real.cells1] + pc_data[, real.cells2])/2 # Literally taking the average of two cells...
  colnames(doublets) <- paste0("SYNTHETIC_", 1:n_doublets)
  pc_data_wdoublets <- cbind(pc_data, doublets)
  cell.names <- colnames(pc_data_wdoublets)
  nCells <- length(cell.names)

  ## Compute PC distance matrix
  cat("Calculating PC distance matrix...\n")
  dist.mat <- fields::rdist(t(pc_data_wdoublets))[,1:n_real.cells] #PC mat passed with PCs in rows!!!

  ## Pre-order PC distance matrix prior to iterating across pK for pANN computations
  cat("Defining neighborhoods...\n")
  for(i in 1:n_real.cells){
    dist.mat[,i] <- order(dist.mat[,i])
  }

  ## Trim PC distance matrix for faster manipulations
  ind <- round(nCells * max(pK)) + 5
  dist.mat <- dist.mat[1:ind, ]

  ## Compute pANN across pK sweep
  cat("Computing pANN across ", length(pK), " pK values...\n")

  sweep.res.list = list()
  list.ind = 0

  for(k in 1:length(pK)){
    cat("pK = ", pK[k], "...")
    pk.temp <- round(nCells * pK[k])
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    colnames(pANN) <- "pANN"
    rownames(pANN) <- real.cells
    list.ind <- list.ind + 1

    for (i in 1:n_real.cells){
      neighbors <- dist.mat[2:(pk.temp + 1),i]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/pk.temp
    }

    sweep.res.list[[list.ind]] <- pANN
  }

  return(sweep.res.list)
}

#     paramSweep_V3_DWM ####
paramSweep_v3_DWM <- function(
    seu,
    PCs=1:10,
    assay='RNA',
    slot='counts',
    reduction='pca', #only needed if running on PC values instead of RNA
    num.cores=1,
    is.pca = F, #if performing paramSweep on PCA values, not gene expression
    sct = FALSE
) {
  require(Seurat); require(fields);

  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05, 0.3, by=0.05)

  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]
  cat("Using ", length(pK), " values for pK testing...\n")

  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands

  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if(nrow(seu@meta.data) > 10000){
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
  }else if (nrow(seu@meta.data) <= 10000){
    real.cells <- rownames(seu@meta.data)
  }else{
    message("ERROR IN SUBSETTING THE DATA")
  }

  if(is.pca){
    if(is.null(seu@reductions[[reduction]])){
      message("Reduction ", reduction, " is missing. Try again... \n")
      return(seu)
    }
    data <- t(seu@reductions[[reduction]]@cell.embeddings[real.cells,PCs])
  }else{
    data <- GetAssayData(seu, assay=assay, slot=slot)[ , real.cells]
  }
  n_real.cells <- ncol(data)

  cat("Using ", n_real.cells, " out of ", ncol(seu), " cells...\n")

  ## Iterate through pN, computing pANN vectors at varying pK
  if(num.cores>1){
    require(parallel)
    cat("     Running on ", num.cores, " cores...")
    cl <- makeCluster(num.cores)

    if(is.pca){ #workflow in PC space
      output2 <- mclapply(
        as.list(1:length(pN)),
        FUN = parallel_paramSweep_PCA,
        n_real.cells=n_real.cells,
        real.cells=real.cells,
        pK=pK,
        pN=pN,
        pc_data=data,
        PCs=PCs,
        mc.cores=num.cores
      )
    }else{ # workflow in feature space
      output2 <- mclapply(
        as.list(1:length(pN)),
        FUN = parallel_paramSweep_v3_DWM,
        n_real.cells,
        real.cells,
        pK,
        pN,
        data,
        orig.commands,
        PCs,
        assay,
        slot,
        sct,
        mc.cores=num.cores
      )
    }
    stopCluster(cl)
  }else{
    cat("     Running on ", num.cores, " core...")
    if(is.pca){ #workflow in PC space
      output2 <- lapply(
        as.list(1:length(pN)),
        FUN = parallel_paramSweep_PCA,
        n_real.cells=n_real.cells,
        real.cells=real.cells,
        pK=pK,
        pN=pN,
        pc_data=data,
        PCs=PCs
      )
    }else{ # workflow in feature space
      output2 <- lapply(
        as.list(1:length(pN)),
        FUN = parallel_paramSweep_v3_DWM,
        n_real.cells,
        real.cells,
        pK,
        pN,
        data,
        orig.commands,
        PCs,
        assay,
        slot,
        sct
      )
    }
  }

  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }

  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec

  cat("Done!\n")

  return(sweep.res.list)
}
#
#     doubletFinder_V3 ####
# Added ability to customize the meta data names
doubletFinder_V3.DWM <- function(
    seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE,
    classification.name=NULL, pANN.name=NULL
){
  require(Seurat); require(fields); require(KernSmooth)

  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }

  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- Cells(seu)
    data <- seu@assays[[assay]]@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)

    ## Store important pre-processing information
    orig.commands <- seu@commands

    ## Pre-process Seurat object
    if (sct == FALSE) {
      cat("Creating Seurat object...\n")
      seu_wdoublets <-CreateSeuratObject(counts = data_wdoublets)

      cat("     Normalizing Seurat object...\n")
      cat("     Finding variable genes...\n")
      cat("     Scaling data...\n")
      cat("     Running PCA...\n")
      seu_wdoublets <- NormalizeData(
        seu_wdoublets,
        normalization.method = orig.commands[[paste0("NormalizeData.",assay)]]@params$normalization.method,
        scale.factor = orig.commands[[paste0("NormalizeData.",assay)]]@params$scale.factor,
        margin = orig.commands[[paste0("NormalizeData.",assay)]]@params$margin
      ) %>% FindVariableFeatures(
        selection.method = orig.commands[[paste0("FindVariableFeatures.",assay)]]$selection.method,
        loess.span = orig.commands[[paste0("FindVariableFeatures.",assay)]]$loess.span,
        clip.max = orig.commands[[paste0("FindVariableFeatures.",assay)]]$clip.max,
        mean.function = orig.commands[[paste0("FindVariableFeatures.",assay)]]$mean.function,
        dispersion.function = orig.commands[[paste0("FindVariableFeatures.",assay)]]$dispersion.function,
        num.bin = orig.commands[[paste0("FindVariableFeatures.",assay)]]$num.bin,
        binning.method = orig.commands[[paste0("FindVariableFeatures.",assay)]]$binning.method,
        nfeatures = orig.commands[[paste0("FindVariableFeatures.",assay)]]$nfeatures,
        mean.cutoff = orig.commands[[paste0("FindVariableFeatures.",assay)]]$mean.cutoff,
        dispersion.cutoff = orig.commands[[paste0("FindVariableFeatures.",assay)]]$dispersion.cutoff
      ) %>% ScaleData(
        features = orig.commands[[paste0("ScaleData.",assay)]]$features,
        model.use = orig.commands[[paste0("ScaleData.",assay)]]$model.use,
        do.scale = orig.commands[[paste0("ScaleData.",assay)]]$do.scale,
        do.center = orig.commands[[paste0("ScaleData.",assay)]]$do.center,
        scale.max = orig.commands[[paste0("ScaleData.",assay)]]$scale.max,
        block.size = orig.commands[[paste0("ScaleData.",assay)]]$block.size,
        min.cells.to.block = orig.commands[[paste0("ScaleData.",assay)]]$min.cells.to.block
      )%>% RunPCA(
        features = orig.commands[[paste0("ScaleData.",assay)]]$features,
        npcs = length(PCs),
        rev.pca =  orig.commands[[paste0("RunPCA.",assay)]]$rev.pca,
        weight.by.var = orig.commands[[paste0("RunPCA.",assay)]]$weight.by.var,
        verbose=FALSE
      )

      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- Cells(seu_wdoublets)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    }

    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)

      cat("Running PCA...\n")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc()
    }

    ## Compute PC distance matrix
    cat("Calculating PC distance matrix...\n")
    dist.mat <- fields::rdist(pca.coord)

    ## Compute pANN
    cat("Computing pANN...\n")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    cat('   nCells = ', nCells,' \n')
    k <- round(nCells * pK)
    cat('   k = ', k,' \n')
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      neighbor.names <- rownames(dist.mat)[neighbors]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
    }

    #Doublet Classification
    cat("Classifying doublets...\n")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"

    # Add in metadata columns with classifications and pANN values
    if(is.null(classification.name)){
      seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    }else if(is.character(classification.name)){
      seu@meta.data[, classification.name] <- classifications
    }else{
      seu@meta.data[, "DF.classifications"] <- classifications
    }

    if(pANN.name){
      seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    }else if(is.character(pANN.name)){
      seu@meta.data[, pANN.name] <- classifications
    }else{
      seu@meta.data[, "DF.pANN"] <- pANN[rownames(seu@meta.data), 1]
    }

    return(seu)
  }
}

#     doubletFinder_V3_v2 ####
doubletFinder_V3.DWM_v2 <- function(
    seu, PCs, pN = 0.25, pK, nExp='auto',
    pANN.cutoff=NULL, # hard limit for distinguishing singlets/doublets
    assay='RNA', reduction='pca',slot='counts', #DWM
    classification.name=NULL, pANN.name=NULL, # meta data naming parameters
    reuse.pANN = FALSE, #changed in my implementation
    sct = FALSE
){
  require(Seurat); require(fields); require(KernSmooth)

  #TODO: check passed seurat object for assay, reduction, etc.
  #TODO: add parallelization with future to the Seurat section

  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if(is.character(reuse.pANN)){ #reuse pANN values, but rename classifications
    if(is.null(seu[[reuse.pANN]])){
      cat("pANN.name '", reuse.pANN, " does not exist... Please try again. \n")
      return(seu)
    }else{
      #Doublet Classification
      cat("Classifying doublets based on previous pANN values (", reuse.pANN, ")...\n",sep = '')

      pANN.old <- seu@meta.data[ , reuse.pANN]
      classifications <- rep("Singlet", length(pANN.old))

      if(is.null(pANN.cutoff)){
        classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
      }else{
        classifications[pANN.old>=pANN.cutoff] <- "Doublet"
      }

      # Add in metadata columns with classifications
      if(is.null(classification.name)){
        seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
      }else if(is.character(classification.name)){
        seu@meta.data[, classification.name] <- classifications
      }else{
        cat("Doublet classifications labeled as 'DF.classifications'...\n")
        seu@meta.data[, "DF.classifications"] <- classifications
      }

      #Don't need to add pANN values again

      return(seu)
    }
  }
  if(reuse.pANN){
    if(is.null(seu[[pANN.name]])){
      cat("pANN.name '", pANN.name, " does not exist... Please try again. \n")
      return(seu)
    }else{
      #Doublet Classification
      cat("Classifying doublets based on previous pANN values...\n")

      pANN.old <- seu[[pANN.name]]
      classifications <- rep("Singlet", length(pANN.old))

      if(is.null(pANN.cutoff)){
        classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
      }else{
        classifications[pANN.old>=pANN.cutoff] <- "Doublet"
      }

      # Add in metadata columns with classifications
      if(is.null(classification.name)){
        seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
      }else if(is.character(classification.name)){
        seu@meta.data[, classification.name] <- classifications
      }else{
        cat("Doublet classifications labeled as 'DF.classifications'...\n")
        seu@meta.data[, "DF.classifications"] <- classifications
      }

      #Don't need to add pANN values again

      return(seu)
    }
  }else{

    ## Make merged real-artifical data
    real.cells <- Cells(seu)
    data <- GetAssayData(seu, assay=assay, slot=slot)
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)

    cat("Creating ", n_doublets, " artificial doublets from ", n_real.cells, " cells...\n")

    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)

    doublets <- (data[, real.cells1] + data[, real.cells2])/2 # Literally taking the average of two cells...
    colnames(doublets) <- paste0("X", 1:n_doublets)
    data_wdoublets <- cbind(data, doublets)

    ## Store important pre-processing information
    orig.commands <- seu@commands

    ## Initialize Seurat object
    if(slot=='counts'){
      cat("Creating Seurat object with artificial doublets...\n")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      if(assay!='SCT'){
        seu_wdoublets <- NormalizeData(
          seu_wdoublets,
          normalization.method = orig.commands[[paste0("NormalizeData.",assay)]]@params$normalization.method,
          scale.factor = orig.commands[[paste0("NormalizeData.",assay)]]@params$scale.factor,
          margin = orig.commands[[paste0("NormalizeData.",assay)]]@params$margin
        )
      }
    }else if(slot=='data'){
      cat("Creating Seurat object with artificial doublets...\n")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets) #don't renormalize if normalized data is provided
    }

    ## Preprocess Seurat object
    if(assay=='SCT'){
      require(sctransform)

      if(slot=='data'){
        cat('WARNING: running SCTransform on normalized data!\n')
      }

      cat("Running SCTransform & PCA...\n")
      seu_wdoublets <- SCTransform(
        seu_wdoublets
      ) %>% RunPCA(
        npcs = length(PCs),
        reduction.name='DOUBLETFINDER_PCA'
      )

      pca.coord <- seu_wdoublets@reductions$DOUBLETFINDER_PCA@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)

      rm(seu_wdoublets); gc()

    }else if(assay=='RNA'){
      cat("     Piping FindVariableFeatures(), ScaleData(), and RunPCA()...\n")
      seu_wdoublets <- FindVariableFeatures(
        seu_wdoublets,
        selection.method = orig.commands[[paste0("FindVariableFeatures.",assay)]]$selection.method,
        loess.span = orig.commands[[paste0("FindVariableFeatures.",assay)]]$loess.span,
        clip.max = orig.commands[[paste0("FindVariableFeatures.",assay)]]$clip.max,
        mean.function = orig.commands[[paste0("FindVariableFeatures.",assay)]]$mean.function,
        dispersion.function = orig.commands[[paste0("FindVariableFeatures.",assay)]]$dispersion.function,
        num.bin = orig.commands[[paste0("FindVariableFeatures.",assay)]]$num.bin,
        binning.method = orig.commands[[paste0("FindVariableFeatures.",assay)]]$binning.method,
        nfeatures = orig.commands[[paste0("FindVariableFeatures.",assay)]]$nfeatures,
        mean.cutoff = orig.commands[[paste0("FindVariableFeatures.",assay)]]$mean.cutoff,
        dispersion.cutoff = orig.commands[[paste0("FindVariableFeatures.",assay)]]$dispersion.cutoff
      ) %>% ScaleData(
        features = orig.commands[[paste0("ScaleData.",assay)]]$features,
        model.use = orig.commands[[paste0("ScaleData.",assay)]]$model.use,
        do.scale = orig.commands[[paste0("ScaleData.",assay)]]$do.scale,
        do.center = orig.commands[[paste0("ScaleData.",assay)]]$do.center,
        scale.max = orig.commands[[paste0("ScaleData.",assay)]]$scale.max,
        block.size = orig.commands[[paste0("ScaleData.",assay)]]$block.size,
        min.cells.to.block = orig.commands[[paste0("ScaleData.",assay)]]$min.cells.to.block
      )%>% RunPCA(
        features = orig.commands[[paste0("ScaleData.",assay)]]$features,
        npcs = length(PCs),
        rev.pca =  orig.commands[[paste0("RunPCA.",assay)]]$rev.pca,
        weight.by.var = orig.commands[[paste0("RunPCA.",assay)]]$weight.by.var,
        reduction.name='DOUBLETFINDER_PCA',
        verbose=FALSE
      )

      pca.coord <- seu_wdoublets@reductions$DOUBLETFINDER_PCA@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    }else{
      cat("     Piping FindVariableFeatures(), ScaleData(), and RunPCA()...\n")
      seu_wdoublets <- FindVariableFeatures(
        seu_wdoublets,
        selection.method = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$selection.method,
        loess.span = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$loess.span,
        clip.max = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$clip.max,
        mean.function = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$mean.function,
        dispersion.function = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$dispersion.function,
        num.bin = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$num.bin,
        binning.method = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$binning.method,
        nfeatures = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$nfeatures,
        mean.cutoff = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$mean.cutoff,
        dispersion.cutoff = orig.commands[[paste('FindVariableFeatures', assay, sep='.')]]$dispersion.cutoff
      ) %>% ScaleData(
        features = orig.commands[[paste('ScaleData', assay, sep='.')]]$features,
        model.use = orig.commands[[paste('ScaleData', assay, sep='.')]]$model.use,
        do.scale = orig.commands[[paste('ScaleData', assay, sep='.')]]$do.scale,
        do.center = orig.commands[[paste('ScaleData', assay, sep='.')]]$do.center,
        scale.max = orig.commands[[paste('ScaleData', assay, sep='.')]]$scale.max,
        block.size = orig.commands[[paste('ScaleData', assay, sep='.')]]$block.size,
        min.cells.to.block = orig.commands[[paste('ScaleData', assay, sep='.')]]$min.cells.to.block
      )%>% RunPCA(
        features = orig.commands[[paste('ScaleData', assay, sep='.')]]$features,
        npcs = length(PCs),
        rev.pca =  orig.commands[[paste('RunPCA', assay, sep='.')]]$rev.pca,
        weight.by.var = orig.commands[[paste('RunPCA', assay, sep='.')]]$weight.by.var,
        reduction.name='DOUBLETFINDER_PCA',
        verbose=FALSE
      )

      pca.coord <- seu_wdoublets@reductions$DOUBLETFINDER_PCA@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)

      rm(seu_wdoublets); gc() # Free up memory
    }

    ## Compute PC distance matrix
    cat("Calculating PC distance matrix...\n")
    dist.mat <- fields::rdist(pca.coord)

    ## Compute pANN
    cat("Computing pANN...\n")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    cat('   nCells = ', nCells,' \n')
    k <- round(nCells * pK)
    cat('   k = ', k,' \n')
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      neighbor.names <- rownames(dist.mat)[neighbors]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
    }

    #Smooth pANN values, compute cutoff for doublet identification - DWM
    #TODO: find local minima in pANN values, and calculate the tiers of multiplets
    if(nExp=='auto'){
      #Smooth pANN

      #Find the two maxes and one local min

      # Set cutoff to the local min
    }



    #Doublet Classification
    cat("Classifying doublets...\n")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"

    # Add in metadata columns with classifications and pANN values
    if(is.null(classification.name)){
      seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    }else if(is.character(classification.name)){
      seu@meta.data[, classification.name] <- classifications
    }else{
      seu@meta.data[, "DF.classifications"] <- classifications
    }

    if(is.null(pANN.name)){
      seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    }else if(is.character(pANN.name)){
      seu@meta.data[, pANN.name] <- pANN[rownames(seu@meta.data), 1]
    }else{
      seu@meta.data[, "DF.pANN"] <- pANN[rownames(seu@meta.data), 1]
    }

    return(seu)
  }
}
#
#     doubletFinder_PCA ####
doubletFinder_PCA <- function(seu, PCs, pN = 0.25, pK, nExp='auto',
                              pANN.cutoff=NULL, # hard limit for distinguishing singlets/doublets
                              reduction='pca', #DWM
                              n.chunks = 1, # number of chunks to split pANN computaiton into
                              classification.name=NULL, pANN.name=NULL, # meta data naming parameters
                              reuse.pANN = FALSE #changed in my implementation to be a logical, whether or not to reuse
){
  require(Seurat); require(fields); require(KernSmooth)

  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if(reuse.pANN){
    if(is.null(seu[[pANN.name]])){
      cat("pANN.name '", pANN.name, " does not exist... Please try again. \n")
      return(seu)
    }else{
      #Doublet Classification
      cat("Classifying doublets based on previous pANN values...\n")

      pANN.old <- seu[[pANN.name]]
      classifications <- rep("Singlet", length(pANN.old))

      if(is.null(pANN.cutoff)){
        classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
      }else{
        classifications[pANN.old>=pANN.cutoff] <- "Doublet"
      }

      # Add in metadata columns with classifications
      if(is.null(classification.name)){
        seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
      }else if(is.character(classification.name)){
        seu@meta.data[, classification.name] <- classifications
      }else{
        cat("Doublet classifications labeled as 'DF.classifications'...\n")
        seu@meta.data[, "DF.classifications"] <- classifications
      }

      #Don't need to add pANN values again

      return(seu)
    }
  }else{
    if(is.null(PCs)){ #Use all PCs if ndims is not set
      PCs <-  1:ncol(tmp@reductions[[reduction]][])
    }
    cat("Using PCs ", PCs[1], " to " , PCs[length(PCs)], " for the reduction ", reduction, "\n")
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    pc_data <- seu@reductions[[reduction]]@cell.embeddings[,PCs] #cells x PCs matrix
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)

    cat("Creating ", n_doublets, " artificial doublets from ", n_real.cells, " cells...\n")
    cat("   (Using the reduction '", reduction, "')\n")

    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (pc_data[real.cells1, ] + pc_data[real.cells2, ])/2 # Literally taking the average of two cells...
    rownames(doublets) <- paste0("SYNTHETIC_", 1:n_doublets)

    pc_data_wdoublets <- rbind(pc_data, doublets)
    cell.names <- rownames(pc_data_wdoublets)
    nCells <- length(cell.names)

    ## Compute PC distance matrix
    if(n.chunks > 1){ # run in chunks to save memory usage
      cat("Calculating PC distance matrix in", n.chunks, "chunks ...\n")
      require(parallel)

      # find chunk indices
      inds <- list()
      tmp.inds <- seq(1, ncol(pc_data_wdoublets), n.chunks)
      inds[[1]] <- c(tmp.inds[1], tmp.inds[2])

      for(i in 2:(length(tmp.inds)-1)){
        inds[[i]] <- c(tmp.inds[i]+1, tmp.inds[i+1])
      }
      inds[[i+1]] <- c(tmp.inds[i], ncol(pc_data_wdoublets))

      print(inds)
      # compute pANN values, 1 chunk at a time
      pANN.list <- lapply(
        inds,
        FUN=function(inds, pc_data, k, n_real.cells){
          dist.mat <- fields::rdist(pc_data[inds,],pc_data)
          print(inds)
          print(dim(dist.mat))

          pANN <- as.data.frame(matrix(0L, nrow = abs(diff(inds)), ncol = 1))
          rownames(pANN) <- rownames(pc_data)[inds]
          colnames(pANN) <- "pANN"

          for(i in 1:nrow(dist.mat)){
            neighbors <- order(dist.mat[, i])
            neighbors <- neighbors[2:(k + 1)]
            neighbor.names <- rownames(dist.mat)[neighbors]

            pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
          }

        },
        pc_data = pc_data_wdoublets,
        k = round(nCells * pK),
        n_real.cells = n_real.cells

      )
      # collapse list
      pANN <- do.call(rbind, pANN.list)

      # remove synthetic doublets from pann list
      pANN <- pANN[real.cells]

    }else{
      ## Compute PC distance matrix
      cat("Calculating PC distance matrix...\n")
      print(dim(pc_data_wdoublets))
      dist.mat <- fields::rdist(pc_data_wdoublets)

      ## Compute pANN
      cat("Computing pANN...\n")
      pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
      rownames(pANN) <- real.cells
      colnames(pANN) <- "pANN"
      cat('   nCells = ', nCells,' \n')
      k <- round(nCells * pK)
      cat('   k = ', k,' \n')
      for(i in 1:n_real.cells){
        if((n_real.cells/20 %% i)==0){
          cat('=')
        }
        neighbors <- order(dist.mat[, i])
        neighbors <- neighbors[2:(k + 1)]
        neighbor.names <- rownames(dist.mat)[neighbors]
        pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      }
    }


    #Doublet Classification
    cat("Classifying doublets...\n")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"

    # Add in metadata columns with classifications and pANN values
    if(is.null(classification.name)){
      seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    }else if(is.character(classification.name)){
      seu@meta.data[, classification.name] <- classifications
    }else{
      seu@meta.data[, "DF.classifications"] <- classifications
    }

    if(is.null(pANN.name)){
      seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    }else if(is.character(pANN.name)){
      seu@meta.data[, pANN.name] <- pANN[rownames(seu@meta.data), 1]
    }else{
      seu@meta.data[, "DF.pANN"] <- pANN[rownames(seu@meta.data), 1]
    }

    return(seu)
  }
}

# end ####
