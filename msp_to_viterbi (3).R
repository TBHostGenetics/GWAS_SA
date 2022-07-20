#!/usr/bin/env Rscript
#Rscript msp_to_viterbi.R <msp file> <bim_file> </ouputdir/outputprefix>

list.of.packages <- c("magrittr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(magrittr)
args = commandArgs(trailingOnly=TRUE)

msp_file = args[1]
bim_file = args[2]
out_file = args[3]

#Read inputs
msp = read.csv(msp_file, 
               header = T, 
               skip = 1,
               sep = "\t")

bim = read.csv(bim_file, 
               header = F,
               sep = "\t")

#Get msp indexes for each snp 
lapply(bim[,4], function(e){
  which(msp$spos <= e & msp$epos > e)
}) ->
  snp_idx
  
snp_idx[[length(snp_idx)]] = nrow(msp)


#Get list of snps not in msp
sapply(snp_idx, function(e){
  length(e)
}) == 0 ->
  no_snp

#Write list of sites not in msp, but in bim
if(sum(no_snp) > 0){
    print(paste("Warning SNPs", paste(bim[which(no_snp),4], collapse = ", "), "missing from msp file."))
    write.table(cbind(bim[which(no_snp),1],bim[which(no_snp),4]), file = paste(out_file, ".notinmsp", sep=""), quote = F, sep = "\t", row.names = F, col.names = F)
}else{
    print("No SNPs missing from msp")
}

#Create viterbi file
if(sum(no_snp) > 0){
    vit = data.frame("pos" = bim[-which(no_snp),4], msp[unlist(snp_idx), 7:ncol(msp)] + 1)
}else{
    vit = data.frame("pos" = bim[,4], msp[unlist(snp_idx), 7:ncol(msp)] + 1)
}
colnames(vit) = gsub("\\.0|\\.1","", colnames(vit))
write.table(vit, file = paste(out_file, ".vit.tsv", sep=""), quote = F, sep = " ", row.names = F, col.names = T)
