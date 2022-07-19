#!/usr/bin/env Rscript

list.of.packages <- c("doSNOW", "bigmemory", "magrittr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(doSNOW)
library(bigmemory)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)

out.file.name = args[1]
print(paste("Output file:", out.file.name, sep = " "))
print(paste("Reading phenotypes:", args[2], sep = " "))
pheno.frame = read.delim(args[2], stringsAsFactors = F, row.names = "sample_id")
print(paste("Reading positions:", args[3], sep = " "))
pos.frame = read.delim(args[3], stringsAsFactors = F, header=T)

print(paste("Reading allele frame:", args[4], sep = " "))
strsplit(args[4], split="/") %>% unlist() %>% head(.,-1) %>% paste(.,collapse="/") %>% paste(., "/", sep="") -> allele.path
strsplit(args[4], split="/") %>% unlist() %>% tail(., n=1) %>% strsplit(., split=".", fixed=T) %>% unlist() %>% head(.,-1) -> allele.base
if(!file.exists(paste(allele.path,allele.base,".bin", sep=""))){
  allele.frame = read.big.matrix(args[4],sep = "\t",header = T, type="char", backingfile = paste(allele.base,".bin", sep=""), backingpath=allele.path, descriptorfile = paste(allele.base,".desc", sep=""), shared = T)
}
print(paste("Reading pop frame:", args[5], sep = " "))
strsplit(args[5], split="/") %>% unlist() %>% head(.,-1) %>% paste(.,collapse="/") %>% paste(., "/", sep="") -> pop.path
strsplit(args[5], split="/") %>% unlist() %>% tail(., n=1) %>% strsplit(., split=".", fixed=T) %>% unlist() %>% head(.,-1) -> pop.base
if(!file.exists(paste(pop.path,pop.base,".bin", sep=""))){
  pop.frame = read.big.matrix(args[5],sep = "\t",header = T, type="char", backingfile = paste(pop.base,".bin", sep=""), backingpath=pop.path, descriptorfile = paste(pop.base,".desc", sep=""), shared = T)
}
print(paste("Reading allele pop frame:", args[6], sep = " "))
strsplit(args[6], split="/") %>% unlist() %>% head(.,-1) %>% paste(.,collapse="/") %>% paste(., "/", sep="") -> allele.pop.path
strsplit(args[6], split="/") %>% unlist() %>% tail(., n=1) %>% strsplit(., split=".", fixed=T) %>% unlist() %>% head(.,-1) -> allele.pop.base
if(!file.exists(paste(allele.pop.path,allele.pop.base,".bin", sep=""))){
  allele.pop.frame = read.big.matrix(args[6],sep = "\t",header = T, type="char", backingfile = paste(allele.pop.base,".bin", sep=""), backingpath=allele.pop.path, descriptorfile = paste(allele.pop.base,".desc", sep=""), shared = T)
}

print("Done reading!")

print(paste("Analysing components", args[7], sep = " "))
pop = args[7]
print(paste("Analysing trait", args[8], sep = " "))
trait = args[8]
print(paste("Using", args[9], "cores", sep = " "))
cores = args[9]
cores = strtoi(cores)

if(pop=="nama"){
  pheno.frame = pheno.frame[,c("age", "sex",trait, "NAMA", "GBR", "EP")]
  components = "NAMA + EP"
}else{
  pheno.frame = pheno.frame[,c("age", "sex",trait, "NAMA", "GBR", "MSL", "SAS", "EAS")]
  components = "NAMA + MSL + GBR + SAS"
}

###NULL
runGlobalModel <- function(model.frame, trait) {
  form = as.formula(paste(trait, "~", "age + sex + ", components, " + allele_dose"))
  return (lm(form, data=model.frame))
}

###GAO
runGlobalModelSummary <- function(model.frame, trait) {
  form = as.formula(paste(trait, "~", "age + sex + ", components, " + allele_dose"))
  return (summary(lm(form, data=model.frame)))
}

###LAO
runLocalModelSummary <- function(model.frame, trait) {
  form = as.formula(paste(trait, "~", "age + sex + ", components, " + pop_dose"))
  return (summary(lm(form, data=model.frame)))
}

###APA
runAlleleModelSummary <- function(model.frame, trait) {
  form = as.formula(paste(trait, "~", "age + sex + ", components, " + allele_dose + pop_dose"))
  return (summary(lm(form, data=model.frame)))
}

###LAAA
runLaaaModel <- function(model.frame, trait) {
  form = as.formula(paste(trait, "~", "age + sex + ", components, " + allele_dose + pop_dose + allele_pop_dose"))
  return (lm(form, data=model.frame))
}

lmp <- function (modelobject) {
  f <- modelobject$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

runModels = function(trait, position, model.frame){
  #Fit the model
  models = c()
  m.null <- runGlobalModel(model.frame, trait)
  m.laaa <- runLaaaModel(model.frame, trait)
  models[["null"]] <- summary(m.null)
  models[["laaa"]] <- summary(m.laaa)
  models[["global"]] <- runGlobalModelSummary(model.frame, trait)
  models[["local"]] <- runLocalModelSummary(model.frame, trait)
  models[["apa"]] <- runAlleleModelSummary(model.frame, trait)
  anova_p <- anova(m.null, m.laaa, test="Chisq")[2,5]
 
  lapply(1:length(models), function(m){
    model_name = names(models[m])
    m = models[[m]]
    coeff = m$coefficients
    cols = names(coeff[1,])
    rows = rownames(coeff)
    res_names = c(paste(model_name,"p_val", sep="_"),paste(model_name,sapply(rows,function(e){rep(e,length(cols))}), cols, sep="_"))
    
    p_val = lmp(m)
    
    res = c(p_val,t(coeff))
    names(res) = res_names
    res
  }) %>%
    do.call("c", .) ->
    model_results
  
  return(c(model_results, anova_p))
}

progress <- function (x, max) {
  percent <- x / max * 100
  cat(sprintf('\r[%-50s] %d%%',
              paste(rep('=', percent / 2), collapse = ''),
              floor(percent)))
  if (x == max)
    cat('\n')
}

cl = makeSOCKcluster(cores)
registerDoSNOW(cl)

n_positions = nrow(pos.frame)

pb <- txtProgressBar(max=n_positions, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

results = foreach(p=1:n_positions, .options.snow=opts, .packages=c("bigmemory", "magrittr"))%dopar%{

  allele.frame.line = attach.big.matrix(paste(allele.path, allele.base, ".desc", sep=""))[p,]
  pop.frame.line = attach.big.matrix(paste(pop.path, pop.base, ".desc", sep=""))[p,]
  allele.pop.frame.line = attach.big.matrix(paste(allele.pop.path, allele.pop.base, ".desc", sep=""))[p,]
  pos_alleles = pos.frame[p,]
    
  model.frame = data.frame(pheno.frame, 
                           allele_dose=as.numeric(allele.frame.line), 
                           pop_dose=as.numeric(pop.frame.line), 
                           allele_pop_dose=as.numeric(allele.pop.frame.line))
  
  #Estimate the alternate allele frequency
  frq <- sum(model.frame$allele_dose)/(dim(model.frame)[1]*2)
  
  model_results = runModels(trait, pos_alleles[1], model.frame)

  result = c(pos_alleles, frq, model_results)
  names(result) = c(names(pos_alleles), "frq", names(model_results))
  
  return(result)
}

stopCluster(cl)

#Write the output
write(paste(names(results[[1]]), collapse = "\t"), out.file.name, append = T)
sapply(results, function(e){
  write(paste(e, collapse = "\t"), out.file.name, append = T)
})

