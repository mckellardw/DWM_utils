###########
# SeqSatUMI.R
#   R script which calculates sequencing saturation in a UMI-aware manner (designed for TXG datasets)
#   David W. McKellar
###########

# Usage: ----
# Rscript /path/to/SeqSatUMI.R BamPath tmpDir outDir N.Iter N.Cores
#  
#     FullBamPath:
#     DedupBamPath
#     outDir: path/name of file to write results to
#     N: number of iterations to check saturation at each sampling rate
#
# Rscript /home/dwm269/DWM_utils/seq_utils

# Outputs ---- 


# Libs & requirements ----
#TODO-add checks/download
require("data.table")
require("ggplot2")
require("parallel")
require("dplyr")

# Get args ----
args=(commandArgs(TRUE))
bampath <- args[1]
tmpdir <- args[2]
outdir <- args[3]
n.iter <- args[4] # default to 3
n.cores <- args[5] # default to 4

# TODO 
# s.range

# Arg checks ----

if(!exists("tmpdir")){
  if(dir.exists("/workdir/dwm269/tmp")){
    tmpdir <- "/workdir/dwm269/tmp" 
  }else{
    message("Temporary directory not found!")
  }
}

if(!exists("n.iter")){
  n.iter=3
}

if(!exists("s.range")){
  s.range <- seq(0.1,0.9,0.1) # sampling range
}

# Get read/UMI totals ----

# Total reads
cat("Counting total reads and UMIs...")
system(
  paste0(
    "samtools view ", bampath, 
    " | grep sS:Z: | sed 's/.*sS:Z:\\([ACGT]*\\).*/\\1/' ",
    " | sort | uniq -c > ", outdir,"reads_per_umi.txt"
  )
)
cat("Done.\n")

umitable <- data.table::fread(
  paste0(outdir,"reads_per_umi.txt"),
  nThread = n.cores,
  skip = 1
)

total.reads = sum(umitable$V1)
total.umis = nrow(umitable)

message(paste("Found", total.umis, "in", total.reads,"reads, now calculating saturation..."))

# Run saturation calculations ----
run.permuts <- expand.grid(s.range,1:n.iter)%>%
  apply(
    MARGIN = 1, 
    function(X) data.frame(s.perc=X[1], iter=X[2]),
    simplify = T
  )

# Check n.core request to make sure it's necessary; update if not
if(n.cores>length(run.permuts)){
  n.cores <- length(run.permuts)
}

saturation.outs <- mclapply(
  run.permuts[c(seq(1,20,9),seq(2,21,9))],
  FUN=function(X){
    tmpfile = paste0(tmpdir,"/sub_",X$s.perc,"_",X$iter,".txt") # file containing 
    
    # subsample bam file
    #TODO: not returning unique reads; need to incorporate sequence info
    tmp.cmd = paste(
      "samtools view -s", X$s.perc, bampath, 
      "| grep CR:Z: | grep CB:Z: ",
      "| sed 's/.*sS:Z:\\([ACGT]*\\).*/\\1/'",
      "| sort | uniq -c >",tmpfile
    )
    system(
      tmp.cmd
    )
    
    # count number of UMIs
    umitable <- data.table::fread(
      tmpfile,
      nThread = n.cores,
      skip = 1
    )
    unique.sub.reads = nrow(umitable)
    
    # compute saturation
    seqSat = 1 - (unique.sub.reads / total.reads)
    
    # remove tmpfile
    file.remove(tmpfile)
    
    # return formatted data.frame for plotting
    return(
      # umitable
      data.frame(
        sample.perc = X$s.perc,
        unique.sub.reads = unique.sub.reads,
        saturation = seqSat,
        iter = X$iter
      )
    )
  },
  mc.cores = n.cores
)

# Clean-up and plotting ----
## collapse saturation outputs into plot-able form ----
df = do.call(rbind, saturation.outs)

## write output for later ----
write.table(
  df,
  file = paste0(outdir,"/saturation_calculations.tsv"),
  delim="\t"
)

## plot it! ----
ggplot(
  df,
  aes(
    x = sample.perc,
    y = seqSat
  )
)+
  geom_point()+
  geom_line()+
  # geom_
  theme_minimal()+
  theme(
    axis.lines=element_line(color="black"),
    axis.ticks=element_line
  )

ggsave(
  file = paste0(outdir,"saturation_plot.pdf"),
  device="pdf"
)