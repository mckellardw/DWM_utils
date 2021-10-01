###############################
# seutils -  Seurat utility functions
# Written by: David McKellar
# version: 1.0
###############################

#############################################################################
## Utils
#############################################################################

# Split violin plots
.GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {

    data <-
      transform(
        data,
        xminv = x - violinwidth * (x - xmin),
        xmaxv = x + violinwidth * (xmax - x)
      )

    grp <- data[1, "group"]

    newdata <-
      plyr::arrange(transform(data, x = if (grp %% 2 == 1)
        xminv
        else
          xmaxv), if (grp %% 2 == 1)
            y
        else-y)

    newdata <- rbind(newdata[1,], newdata, newdata[nrow(newdata),], newdata[1,])

    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

    if (length(draw_quantiles) > 0 &
        !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                1))
      quantiles <-
        ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <-
        data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <-
        rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <-
        GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin",
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
      ggplot2:::ggname("geom_split_violin",
                       GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

.geom_split_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "identity",
           ...,
           draw_quantiles = NULL,
           trim = TRUE,
           scale = "area",
           na.rm = FALSE,
           show.legend = NA,
           inherit.aes = TRUE) {
    layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = .GeomSplitViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        trim = trim,
        scale = scale,
        draw_quantiles = draw_quantiles,
        na.rm = na.rm,
        ...
      )
    )
  }

#############################################################################
## Differential Gene Expression Plots
#############################################################################
# Build volcano plot after running FindMarkers()
#   *Note- requires Seurat v4 or greater (b/c colnames in the output file changed for some reason...)
ggVolcano_v2 <- function(
  markers=NULL,
  expression=NULL,
  seu=NULL,
                      logFC_filter = 1,
                      neg.log.pval.Thresh=50,
                      pct.thresh = 0.4,
                      plotTitle='Volcano Plot',
                      pseudocount=10^-300,
                      fill.cols = c("#000000","#cf4c59","#2c30a8"),
                      xlim=NULL, ylim=NULL,
                      genes=NULL,
                      gene.text.size=6,  repel=T, nudge_x=-0.1,
                      dot.scale=1,
                      line.width=0.5, segment.color="gray",
                      pt.size=1, pt.alpha=1
                    ){
  # TODO: add split.by - generalize?

  require(ggrepel)

  if(is.null(rownames(markers))){
    stop("Differential gene expression matrix needs gene names as rownames()...")
  }
  if(plotTitle=='Volcano Plot'){
    cat('Warning: no plot title given...\n')
  }

  # Build data.frame ####
  df <- data.frame(
    p_val_adj = -log10(markers$p_val_adj+pseudocount),
    avg_log2FC = markers$avg_log2FC,
    genes = rownames(markers),
    p_val_filter = -log10(markers$p_val_adj+pseudocount)>neg.log.pval.Thresh, # Values that pass p value threshold
    FC_filter.low = markers$avg_log2FC < -1*logFC_filter,
    FC_filter.high= markers$avg_log2FC > logFC_filter
  )

  if(is.null(seu)){
    df$pct =
      apply(markers,1,
            function(X){
              if(as.numeric(X["avg_log2FC"])>0){
                (return(as.numeric(X["pct.1"])))
              }else{
                return(as.numeric(X["pct.2"]))
              }
            }
      )

    df$expr.filter = df$pct > pct.thresh # which genes to label
  }else{
    #TODO: add in gene expression data
  }

  # Add colors ####
  df$cols <- rep(x='none', times=length(df$genes))
  df$cols[df$p_val_filter & df$FC_filter.high] <- 'high'
  df$cols[df$p_val_filter & df$FC_filter.low] <- 'low'
  df$cols <- factor(df$cols,levels=c("none", "high", "low"))

  cat(table(df$colors!="none" & df$expr.filter)[2], "genes drawn \n")
  # Build ggplot object ####
  out.gg <- ggplot(
    data = df,
    aes(x=avg_log2FC, y=p_val_adj)
  ) +
    geom_hline(# pval dotted line
      yintercept = neg.log.pval.Thresh,
      linetype='dashed', color='lightgray', size = line.width
    )  +
    geom_vline(# logFC dotted line
      xintercept = c(logFC_filter, -logFC_filter),
      linetype='dashed', color='lightgray', size = line.width
    ) +
    geom_vline(xintercept = 0,  color='black', size = line.width) +
    geom_point(
      pch=21, size=pt.size, alpha=pt.alpha,
      aes(
        fill = cols, col=cols
        # size = pct # scale dot size according to relative pct expression
      )
    )
  if(repel){
    out.gg <- out.gg + geom_text_repel(
      data=df[df$cols!="none" & df$expr.filter,],
      size=gene.text.size*(5/14), #convert to same scale as normal ggplot
      segment.size=0.25,
      segment.alpha = 0.8,
      segment.color = segment.color,
      point.padding = 0.4,
      aes(
        label=genes,
        col=cols
      )
    )
  }else{
    out.gg <- out.gg + geom_text(
      data=df[df$cols=="both" & df$expr.filter,],
      size=gene.text.size*(5/14), #convert to same scale as normal ggplot
      # position=position_dodge(width = 1),
      nudge_x = nudge_x,
      aes(
        label=genes,
        col=cols
      )
    )
  }
  # Finish plot ####
  out.gg <- out.gg +
    labs(
      x="log2_FC",
      y="-log10(adj._p_val)",
      size = "Pct. Expr."
    ) +
    ggtitle(plotTitle) +
    guides(color=FALSE, fill=FALSE) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),#element_line(colour = "black", size = line.width),
      axis.ticks=element_line(colour = "black", size = line.width),
      panel.border = element_rect(color='black', size=line.width, fill=NA),
      plot.title = element_text(face="bold", hjust=0.5)
    ) +
    scale_fill_manual(values=fill.cols,aesthetics = c("col","fill"))+
    scale_x_continuous(expand=c(0.01,0.01))+
    scale_y_continuous(expand=c(0.01,0.01))

  return(out.gg)
}


#############################################################################
## Spatial/Visium plot functions
#############################################################################

# Generate feature plots given a list of Visium/Seurat objects
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
  combine=TRUE,
  abs.heights=TRUE, # Use absolute heights to size each subplot
  nrow=NULL,
  ncol=NULL,
  option="viridis", #viridis option
  verbose=FALSE
){
  require(Seurat)
  require(ggplot2)
  require(viridis)

  if(verbose){cat(paste0("Plotting Visium data, using the assay ",assay,"!\n"))}

  if(is.null(alt.titles)){
    alt.titles=features
  }

  seu.list <- lapply(
    seu.list,
    FUN = function(SEU){
      SEU@active.assay <- assay
      return(SEU)
    }
  )

  # Get expression limits for each gene, across all datasets
  gene.lims <- lapply(
    features,
    FUN = function(FEAT){
      out.max <- lapply(
        seu.list,
        FUN = function(SEU){
          if(FEAT %in% rownames(SEU)){
            return(max(GetAssayData(SEU,assay=assay,slot=slot)[FEAT,]))
          }else if(FEAT %in% colnames(SEU@meta.data)){
            return(max(SEU@meta.data[,FEAT]))
          }else{
            message(FEAT, " not found!")
          }
        }
      ) %>% unlist() %>% max()
      return(c(10^-100,out.max))
    }
  )

# Get plot heights
if(abs.heights){
  heights <- lapply(
    seu.list,
    FUN=function(SEU) abs(diff(range(SEU@reductions[[reduction]]@cell.embeddings[,2])))
  ) %>% unlist()
  if(verbose){
    message(paste0("Using these plot heights:"))
    print(heights)
  }
}else{
  heights <- rep(1,length(features))
}

# Plot
  plot.list <- list()
  for(i in 1:length(features)){
    tmp <- lapply(
      seu.list,
      FUN = function(SEU)
        FeaturePlot(
          SEU,
          slot = slot,
          features = features[i],
          pt.size = pt.size,
          reduction=reduction
        ) +
        scale_color_viridis(
          limits=unlist(gene.lims[i]),
          option=option,
          na.value = gray(0.42)
        )+
        theme(
          plot.margin = unit(rep(0,4), "inches"),
          axis.ticks = element_blank(),
          axis.text=element_blank(),
          axis.title = element_blank(),
          axis.line=element_blank(),
          plot.title = element_blank(),
          legend.position=legend.position,
          legend.title = element_text(size=font.size,face="bold", hjust=0.5),
          legend.text = element_text(size=font.size,face="bold")
        )
    )
    tmp[[1]] <- tmp[[1]] +
      theme(
        plot.title = element_text(size=font.size,face="bold.italic",  vjust=1)
      ) +
      labs(title=alt.titles[i])
    plot.list[[i]] <- tmp
  }

  for(i in 1:length(plot.list[[1]]) ){
    plot.list[[1]][[i]] <- plot.list[[1]][[i]] +
      theme(axis.title.y = element_text(size=font.size, face="bold", color="black"))

      if(!is.null(sample.titles)){ # add sample titles
        plot.list[[1]][[i]] <- plot.list[[1]][[i]] +
          labs(y=sample.titles[i])
      }
  }

  plot.list <- lapply(
    plot.list,
    FUN = function(X){

      wrap_plots(
        X,
        ncol=1,
        heights=heights,
        guides="collect"
      )&theme(
        legend.position=legend.position,
        legend.margin = margin(0,0,0,0,"inches")
      )
    }
  )

  if(verbose){cat("Done plotting Visium data!\n")}

  if(combine){
    return(
      wrap_plots(
        plot.list,
        nrow=nrow,
        ncol=ncol
      )
    )
  }else{
    return(plot.list)
  }
}



# Generate split violin plots for Visium data to look at co-occurence of cell types & gene expression
#   Co-occurence is based on bayesprism theta values
visCoocVlnPlot<-function(
  SEU,
  assay.ge='Spatial',# gene expression assay
  assay.ct, #cell type assay
  slot.ge="data",
  features.lr, # ligand-receptor(in that order)  pair
  features.ct=NULL, # 2 celltypes to compare
  bp.thresh=0.1, # Lower threshold for bayesprism theta values to determine presence of a cell type
  scale.data=T,
  scale.min=-2.5,
  scale.max=2.5,
  width=0.9,
  #TODO
  legend.position="bottom",
  pt.size=0,
  font.size=8,
  verbose=T
){
  require(ggplot2)
  require(dplyr)

  #Check inputs
  if(is.null(features.lr)){
    message("Need genes to plot (features.lr)...")
    return(NULL)
  }
  if(is.null(features.ct)|length(features.ct)<2){
    message("Need 2 cell types to compare (features.ct)...")
    return(NULL)
  }

  # build df
  df <- data.frame(
    cell=Cells(SEU),
    GetAssayData(SEU, assay=assay.ge,slot=slot.ge)[features.lr[1:2],]%>%t(),
    GetAssayData(SEU, assay=assay.ct)[features.ct[1:2],]%>%t()
  )
  colnames(df)<-c("cell",features.lr[1:2],features.ct[1:2])

  df$cooc <- rep("L-/R-",nrow(df))

  df$cooc[df[features.ct[1]]>=bp.thresh & df[features.ct[2]]>=bp.thresh] <- "L+/R+" #both present
  df$cooc[df[features.ct[1]]>=bp.thresh & df[features.ct[2]]<bp.thresh] <- "L+/R-" #ligand cell type only
  df$cooc[df[features.ct[1]]<bp.thresh & df[features.ct[2]]>=bp.thresh] <- "L-/R+" # receptor cell type only

  if(verbose){
    print(table(df$cooc))
  }

  #TODO- scale gene expression data
  # if(scale.data){
  #   df[features.lr[1]] <- scale(x = df[features.lr[1]])
  #   df[features.lr[2]] <- scale(x = df[features.lr[2]])
  #
  #   df[features.lr[1]] <- MinMax(data = df[features.lr[1]], min = scale.min, max = scale.max)
  #   df[features.lr[2]] <- MinMax(data = df[features.lr[2]], min = scale.min, max = scale.max)
  # }

  # melt df for ggplot
  df <- reshape2::melt(
    df,
    measure.vars=features.ct,
    variable.name="cell_type",
    value.name = "theta"
  ) %>% reshape2::melt(
    measure.vars=features.lr[1:2],
    variable.name="gene",
    value.name = "log_norm_expr"
  )

  # plot!
  out.plot <-ggplot(
    df,
    aes(
      x=cooc,
      y=log_norm_expr,
      fill=gene
      # group=cooc
    )
  )+
    .geom_split_violin(
      width=width,
      show.legend = T
    )+
    labs(title = paste0(features.ct[1]," -> ", features.ct[2]))+
    vln.theme+
    theme(
      legend.position = "right",
      plot.title = element_text(color="black", face="bold")
    )

    return(out.plot)
}
