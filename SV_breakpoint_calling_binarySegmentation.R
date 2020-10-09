require(changepoint)
require(zoo)
require(mcp)
require(data.table)
require(GenomicRanges)


# Define function that takes average depth of tiling windows of size w
tiling_average <- function(mydf=depth_BK364.03,w=1000){
  max_len <- nrow(mydf)
  myidx <- seq(1,max_len,w)
  
  # Array of START and END coordinates
  myranges <- (cbind(myidx[-length(myidx)],myidx[-1]))
  ranges_ir <- IRanges(start=(myranges[,1]),end=(myranges[,2]))
  
  # Setup binnedAverage function to find mean coverage across tiling windows
  depth_vector <- RleList(c(mydf$V3))
  
  # Calculate mean average in windows
  average_vector <- aggregate(depth_vector[[1]],ranges_ir,FUN=mean)
  
  # Bind start & end coordinates with mean average 
  final_out <- cbind(myranges,average_vector) 
  
  # Return function output
  return(final_out)
}


# Input read depth data (generated with "samtools depth").
# Update paths to where the output files to the directory where 'samtools depth' output is located.

depth_S016 <- fread("CNV_samples.tar/CNV_samples/S016_CTNND2_dup.tab",header=F,sep='\t')
depth_BK364.03 <- fread("CNV_samples.tar/CNV_samples/BK364-03_1p36.11_dup.tab",header=F,sep = '\t')
depth_BK144_03 <- fread("CNV_samples.tar/CNV_samples/BK144-03_22q13.3_del.tab",header=F,sep='\t')
depth_BK397_101 <- fread("CNV_samples.tar/CNV_samples/BK397-101_16p11.2_del.tab",header=F,sep='\t')
depth_BK430_103 <- fread("CNV_samples.tar/CNV_samples/BK430-103_16p11.2_dup.tab",header=F,sep='\t')
depth_BK506_03 <- fread("CNV_samples.tar/CNV_samples/BK506-03_5p15.33_del.tab",header=F,sep='\t')
depth_BK294_03 <- fread("CNV_samples.tar/CNV_samples/BK294-03_22q11.2_dup.tab",header=F,sep='\t')
depth_BK482_101 <- fread("CNV_samples.tar/CNV_samples/BK482-101_1q21.1_dup.tab",header=F,sep='\t')
depth_BK487_101 <- fread("CNV_samples.tar/CNV_samples/BK487-101_1q21_del.tab",header=F,sep='\t')
depth_BK180_03 <- fread("CNV_samples.tar/CNV_samples/BK180-03_15q11-13_dup.tab",header=F,sep='\t')
depth_S014 <- fread("CNV_samples.tar/CNV_samples/S014-control-CNV.tab",header=F,sep='\t')
depth_S020_1 <- fread("CNV_samples.tar/CNV_samples/S020-control-CNV-chr4-chr14_1.bed",header=F,sep='\t')
depth_S020_2 <- fread("CNV_samples.tar/CNV_samples/S020-control-CNV-chr4-chr14_2.bed",header=F,sep='\t')
depth_S020_3 <- fread("CNV_samples.tar/CNV_samples/S020-control-CNV-chr4-chr14_3.bed",header=F,sep='\t')
depth_S020_4 <- fread("CNV_samples.tar/CNV_samples/S020-control-CNV-chr4-chr14_4.bed",header=F,sep='\t')
depth_S021 <- fread("CNV_samples.tar/CNV_samples/S021-control-CNV-8p.tab",header=F,sep='\t')
depth_S022_1 <- fread("CNV_samples.tar/CNV_samples/S022-control-CNV-4q_1.bed",header=F,sep='\t')
depth_S022_2 <- fread("CNV_samples.tar/CNV_samples/S022-control-CNV-4q_2.bed",header=F,sep='\t')
depth_S023_1 <- fread("CNV_samples.tar/CNV_samples/S023-control-CNV-ring18_1.bed",header=F,sep='\t')
depth_S023_2 <- fread("CNV_samples.tar/CNV_samples/S023-control-CNV-ring18_2.bed",header=F,sep='\t')
depth_S035_1 <- fread("CNV_samples.tar/CNV_samples/S035-control-CNV-8q-16p_1.tab",header=F,sep='\t')
depth_S035_2 <- fread("CNV_samples.tar/CNV_samples/S035-control-CNV-8q-16p_2.tab",header=F,sep='\t')
depth_S036_1 <- fread("CNV_samples.tar/CNV_samples/S036-control-CNV-10q_1.tab",header=F,sep='\t')
depth_S036_2 <- fread("CNV_samples.tar/CNV_samples/S036-control-CNV-10q_2.tab",header=F,sep='\t')
depth_S036_3 <- fread("CNV_samples.tar/CNV_samples/S036-control-CNV-10q_3.tab",header=F,sep='\t')
depth_S036_4 <- fread("CNV_samples.tar/CNV_samples/S036-control-CNV-10q_4.tab",header=F,sep='\t') 



# Tiling averages with variable-sized windows:
depth_BK364.03_tile <- tiling_average(mydf = depth_BK364,w = 1000)
depth_S016_tile <- tiling_average(mydf = depth_S016,w = 1000)
depth_BK144.03_tile <- tiling_average(mydf = depth_BK144.03,w = 1000)
depth_BK290.03_tile <- tiling_average(mydf = depth_BK290.03,w = 1000)
depth_BK397.101_tile <- tiling_average(mydf = depth_BK397.101,w = 1000)
depth_BK430.103_tile <- tiling_average(mydf = depth_BK430.103,w = 1000)
depth_BK506.03_tile <- tiling_average(mydf = depth_BK506.03,w = 1000)
depth_BK180.03_tile <- tiling_average(mydf = depth_BK180.03,w = 1000)
depth_BK482.101_tile <- tiling_average(mydf = depth_BK482.101,w = 1000)
depth_BK487.101_tile <- tiling_average(mydf = depth_BK487.101,w = 1000)
depth_S016_tile <- tiling_average(mydf=depth_S016,w=1000)
depth_BK364_03_tile <- tiling_average(mydf=depth_BK364.03,w=1000)
depth_BK144_03_tile <- tiling_average(mydf=depth_BK144_03,w=1000)
depth_BK397_101_tile <- tiling_average(mydf=depth_BK397_101,w=1000)
depth_BK430_103_tile <- tiling_average(mydf=depth_BK430_103,w=1000)
depth_BK506_03_tile <- tiling_average(mydf=depth_BK506_03,w=1000)
depth_BK294_03_tile <- tiling_average(mydf=depth_BK294_03,w=1000)
depth_BK482_101_tile <- tiling_average(mydf=depth_BK482_101,w=1000)
depth_BK487_101_tile <- tiling_average(mydf=depth_BK487_101,w=1000)
depth_BK180_03_tile <- tiling_average(mydf=depth_BK180_03,w=1000)
depth_S014_tile <- tiling_average(mydf=depth_S014,w=1000)
depth_S020_1_tile <- tiling_average(mydf=depth_S020_1,w=1000)
depth_S020_2_tile <- tiling_average(mydf=depth_S020_2,w=1000)
depth_S020_3_tile <- tiling_average(mydf=depth_S020_3,w=1000)
depth_S020_4_tile <- tiling_average(mydf=depth_S020_4,w=1000)
depth_S021_tile <- tiling_average(mydf=depth_S021,w=1000)
depth_S022_1_tile <- tiling_average(mydf=depth_S022_1,w=1000)
depth_S022_2_tile <- tiling_average(mydf=depth_S022_2,w=1000)
depth_S023_1_tile <- tiling_average(mydf=depth_S023_1,w=1000)
depth_S023_2_tile <- tiling_average(mydf=depth_S023_2,w=1000)
depth_S035_1_tile <- tiling_average(mydf=depth_S035_1,w=1000)
depth_S035_2_tile <- tiling_average(mydf=depth_S035_2,w=1000)
depth_S036_1_tile <- tiling_average(mydf=depth_S036_1,w=1000)
depth_S036_2_tile <- tiling_average(mydf=depth_S036_2,w=1000)
depth_S036_3_tile <- tiling_average(mydf=depth_S036_3,w=1000)
depth_S036_4_tile <- tiling_average(mydf=depth_S036_4,w=1000)


###########################
####### CHANGEPOINT #######
###########################
# Refine CNV breakpoints

# List full sample list (26 CNV regions)

sample_vector <- c("depth_S016_tile", "depth_BK364_03_tile", "depth_BK144_03_tile", "depth_BK397_101_tile", "depth_BK430_103_tile", "depth_BK506_03_tile", "depth_BK294_03_tile", "depth_BK482_101_tile", "depth_BK487_101_tile", "depth_BK180_03_tile", "depth_S014_tile", "depth_S020_1_tile", "depth_S020_2_tile", "depth_S020_3_tile", "depth_S020_4_tile", "depth_S021_tile", "depth_S022_1_tile", "depth_S022_2_tile", "depth_S023_1_tile", "depth_S023_2_tile", "depth_S035_1_tile", "depth_S035_2_tile", "depth_S036_1_tile", "depth_S036_2_tile", "depth_S036_3_tile", "depth_S036_4_tile")


par(mfrow=c(4,7))

for(sample_index in 1:length(sample_vector)){
  
  # Conduct binary segmentation analysis using a user-defined Q-value (expected number of changepoints).
  binseg <- cpt.meanvar(eval(parse(text=sample_vector[sample_index]))[,3],test.stat='Normal',method='BinSeg',Q=10,penalty="BIC")
  
  # Print out index of each changepoint (i.e., CNV breakpoint)
  cpts(binseg)
  
  # Generate mean and variance of each region
  param.est(binseg)

  # Plot a read-depth profile plot with 
  plot(binseg,cpt.width=3,cpt.col='red',xlab="Position (x 1kb)",ylab="Read depth (fold)",main=toString(sample_vector[sample_index]))
}

