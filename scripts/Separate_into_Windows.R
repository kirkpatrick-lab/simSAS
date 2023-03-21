#SEPARATE CHROM BY WINDOW - 
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

file = args[1]

chrom = read.table(file, header = T)
windows = unique(chrom$Window)

for(win in windows){
  write.table(chrom[chrom$Window == win,], paste("data/chr_", chrom[chrom$Window == win,]$Chrom[1], "_win_", win, ".counts", sep = ""), quote = F, col.names = T, row.names = F)
}
