#Simple wrapper script for breakpoint R using run parameters Zeid optimized during his Sept 30 2019 experiment
#Set up for paired end reads-----just change to pairedEndReads=FALSE for SE reads
#############################################################################################################################################################

library(breakpointR)

#store the number of threads to use
args <- as.numeric(commandArgs(T))

#run breakpointr and save to PDF
pdf()
breakpointr(inputfolder="goodPredictions_JAN2021/", outputfolder="goodPredictions_JAN2021_BPR/", pairedEndReads=TRUE, numCPU=20,windowsize=175,binMethod="reads",peakTh=0.3875,min.mapq=7.75,trim=6.5,background=0.15)
dev.off()
