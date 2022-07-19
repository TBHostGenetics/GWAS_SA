#!/usr/bin/env Rscript

list.of.packages <- c("magrittr", "rstatix", "purrr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(magrittr)
library(rstatix)
library(purrr)

args = commandArgs(trailingOnly=TRUE)

mspfile = args[1]
simfile = args[2]
output = args[3]

#Get Anc pops form msp
readLines(mspfile) %>%
  strsplit(., split = "\t") %>%
  sapply(., function(e){
    strsplit(e, split = " ")
  }) %>%
  unlist() %>%
  .[grepl("=", .)] %>%
  sapply(., function(e){
    strsplit(e, split = "=")[[1]][1]
  }) ->
  pops

msp = read.delim(mspfile, skip = 1, sep = "\t")

sim = readLines(simfile)
sim = do.call("rbind", lapply(sim, function(e) as.numeric(unlist(strsplit(e, split = "")))))
sim = t(sim)

sim %>% c %>% unique -> ancs

glob_anc = sapply(1:ncol(sim), function(e){
  sim[,e] %>%
    c(ancs) %>%
    freq_table %>%
    data.frame %>%
    .[,2]
})

correct = sapply(1:ncol(msp), function(e){
  sim[which(sim[,e] == msp[,e]),e] %>%
    c(ancs) %>%
    freq_table %>%
    data.frame %>%
    .[,2] %>%
    map(., function(f){f - 1}) %>%
    unlist / glob_anc[,e]
})

correct[which(glob_anc == 1)] = NA


write.table(correct, file = output, row.names = F, col.names = F, sep = "\t", quote = F)

