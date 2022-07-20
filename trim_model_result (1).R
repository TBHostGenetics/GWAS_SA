#!/usr/bin/env Rscript

list.of.packages <- c("magrittr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(magrittr)

args = commandArgs(trailingOnly=TRUE)

file = args[1]
p_val = args[2]
p_val = as.numeric(p_val)
output = args[3]

readLines(file, n=1) %>%
  strsplit(., split = "\t") %>%
  unlist() %>%
  grep("p_val", .) %>%
  paste(., collapse = ",") ->
  p_cols

system(paste("cut -f", p_cols, " ", file, sep = ""), intern = T) %>%
  strsplit(., "\t") %>%
  .[-1] %>%
  do.call("rbind", .) %>%
  apply(., 2, function(e){
    which(as.numeric(e) < p_val)
  }) %>% 
  unlist() %>%
  sort() %>%
  unique() %>%
  c(0, .) ->
  rows

readLines(file) %>%
  .[rows+1] %>%
  writeLines(., output)
