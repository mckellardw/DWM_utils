###############################
# seutils -  Seurat utility functions
# Written by: David McKellar
# version: 1.0
###############################


# Build volcano plot after running FindMarkers()
#   *Note- requires Seurat v4 or greater (b/c colnames in the output file changed for some reason...)
ggVolcano_v2 <- function(markers=NULL, expression=NULL,
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
                      pt.size=1, pt.alpha=1){
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
