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
