################################################################################### 
# Custom plotting themes and color pallettes for single-cell analyses
# Compiled by David McKellar
###################################################################################

require(ggplot2)

scThemes <- list()

scThemes$umap <- theme(
  axis.line = element_line(color = "black", size = line.width),
  axis.title = element_text(face='bold',size = small.font, hjust = 0, vjust = 1), 
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  legend.text = element_text(size=small.font,color="black"),
  legend.title = element_text(size=big.font,color="black")
)

fig1bcd.theme <- theme_minimal()+
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_line(color='gray'),
    panel.grid.minor = element_line(color='light gray'),
    panel.border = element_rect(colour = "black", fill=NA, size=line.width),
    axis.title = element_text(face = 'bold',size = big.font, hjust = 0.5, vjust = 0.5), 
    axis.text = element_text(size = small.font, color = "black",hjust=0.5),
    axis.ticks = element_line(color="black"),
    legend.background = element_rect(color = "black", fill='white', size=0.5),
    legend.text = element_text(size = small.font, hjust = 0, vjust = 0.5),
    legend.title = element_text(face = 'bold',size = small.font, hjust = 0.5, vjust = 0.5),
    legend.position = c(0.2,0.6)
  ) 

scThemes$bar <- theme_minimal() +
  theme(
    title=element_text(face="bold.italic", hjust=0.5,vjust = 0.5, size=big.font),
    plot.title=element_text(face="bold.italic", hjust=0.5,vjust = 0.5, size=big.font),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=small.font, color="black"),
    legend.title = element_text(size=big.font, color="black"),
    axis.line = element_line(color='black',size = line.width),
    axis.title = element_text(face = 'bold',size = small.font, hjust = 0.5, vjust = 1), 
    axis.text.x=element_text(size = small.font, color = "black",angle = 0,hjust=0.5),
    axis.text.y=element_text(size = small.font, color = "black",hjust=0.5,angle = 90),
    axis.ticks = element_line(color="black")
  )+ NoLegend()

scThemes$pie <- 
  theme_minimal() + 
  theme(
    plot.background = element_blank(),
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    title=element_text(face="bold.italic", hjust=0.5,vjust = 0.5, size=big.font),
    plot.title=element_text(face="bold.italic", hjust=0.5,vjust = 0.5, size=big.font)
  ) + NoLegend()

scThemes$dot <- theme(
  axis.line = element_blank(),
  panel.border = element_rect(colour = "black", fill=NA, size=1),
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(),
  axis.text.x=element_text(angle = 90,hjust = 1,vjust= 0.5, size=small.font),
  axis.text.y=element_text(angle=0,hjust = 1,vjust= 0.5, size=small.font),
  legend.text = element_text(size=small.font,color="black"),
  legend.title = element_text(size=big.font,color="black"),
  panel.grid.major = element_line(colour = "gray", size = 0.5)
)

scThemes$vln <- theme(
  panel.background = element_blank(),
  axis.line=element_line(color="black"),
  legend.text = element_text(size=small.font, color="black"),
  panel.grid = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y=element_text(size=big.font, color="black"),
  axis.text = element_text(color="black",size=small.font),
  axis.ticks.x = element_blank()
)

scThemes$vln <- theme_minimal() +
  theme(
    axis.text=element_text(color="black",size=small.font),
    axis.title=element_text(color="black",size=big.font, face="bold"),
    panel.grid = element_blank()
  )