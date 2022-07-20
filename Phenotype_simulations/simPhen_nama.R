library(PhenotypeSimulator)
library(magrittr)
library(bigmemory)

read_background = function(rdata, idx){
  load(rdata)
  genBg$shared = genBg$shared[idx,]
  genBg$independent = genBg$independent[idx,]
  genBg$evec_kinship = genBg$evec_kinship[idx,]
  return(genBg)
}

filter_genotypes = function(data, ids, vars=NULL){
  if(is.null(vars)){
    vars = c(1:ncol(data$genotypes))
  }
  genotypes <- as(data$genotypes[ids,vars], 'numeric')
  data$fam <-  data$fam[ids,]
  id_samples <- rownames(genotypes)
  id_snps <- colnames(genotypes)
  format_files <- list(plink_fam = data$fam, plink_map = data$map)
  colnames(genotypes) <- id_snps
  rownames(genotypes) <- id_samples
  genotypes <- list(genotypes=genotypes, id_snps=id_snps, id_samples=id_samples,
                    format_files = format_files)
  return(genotypes)
}

check_freqs = function(causalSNPs, sample_ids){
  sapply(sample_ids, function(e){
    filter_genotypes(data, e, colnames(causalSNPs)) %>%
      .$genotypes%>%
      apply(., 2, getAlleleFrequencies)
  })
}

#Read in data
"/home/gerald/Documents/PhD/papers/paper4/multi_sim_5000_all_maf05_biallelic_snps_ids_filtered_sampleids" %>%
  snpStats::read.plink() ->
  data

#Get population ids
ids = rownames(data$genotypes)
nama_idx = grep("Nama", ids)
sac_idx = grep("SAC", ids)
gbr_idx = grep("GBR", ids)
ep_idx = grep("EP", ids)
sas_idx = grep("SAS", ids)
eas_idx = grep("EAS", ids)
msl_idx = grep("MSL", ids)

# total SNP effect on phenotype: 0.01
genVar <- 0.6 #Proportion [double] of total genetic variance.
noiseVar <- 1 - genVar #Proportion [double] of total noise variance.
totalSNPeffect <- 0.12
totalLAeffect <- 0.06
totalInteraction <- 0.04
h2s <- totalSNPeffect/genVar #Proportion [double] of genetic variance of genetic variant effects.
lambda <- totalLAeffect/genVar #Proportion [double] of genetic variance of local ancestry effects.
psy <- totalInteraction/genVar #Proportion [double] of genetic variance of the interaction of allelic and local ancestry effects.
phi <- 0.4 #Proportion [double] of noise variance of observational noise effects
#rho <- 0.1 #Proportion [double] of noise variance of correlated effects
delta <- 0.6 #Proportion [double] of noise variance of non-genetic covariate effects
shared <- 0.8 #
independent <- 1 - shared #

# simulate infinitesimal genetic effects
genBg = read_background("/home/gerald/Documents/PhD/papers/paper4/bg_p10_3e5_ids.RData", nama_idx)
#genBg = read_background("/home/gveeden/PhD/data/gwas_sim/bg_p10_3e5_ids.RData", nama_idx)
noiseFixed <- noiseFixedEffects(N = 5000, P = 10, NrFixedEffects = 2, 
                                NrConfounders = c(1, 1),
                                pIndependentConfounders = c(0, 0),  
                                distConfounders = c("bin", "unif"),
                                mConfounders = 45,
                                sdConfounders = 25,
                                probConfounders = 0.4)

# simulate observational noise effects
noiseBg <- noiseBgEffects(N = 5000, P = 10)

#Rescale background phenotypic components
noiseBg_shared_scaled <- rescaleVariance(noiseBg$shared, shared * phi* noiseVar)
noiseBg_independent_scaled <- rescaleVariance(noiseBg$independent, 
                                              independent * phi* noiseVar)
noiseFixed_shared_scaled <- rescaleVariance(noiseFixed$shared, shared * delta * 
                                              noiseVar)
noiseFixed_independent_scaled <- rescaleVariance(noiseFixed$independent, 
                                                 independent * delta * noiseVar)

nama_genotypes = filter_genotypes(data, nama_idx)
#nama_genotypes = filter_genotypes(data, sac_idx)
nama_freqs = list(nama_idx, gbr_idx, ep_idx)


# simulate 10 genetic variant effects (from non-standardised SNP genotypes)
causalSNPs <- getCausalSNPs(N=5000, genotypes = nama_genotypes$genotypes, 
                            NrCausalSNPs = 10, verbose = FALSE)

#Check allele frequencies of causalSNPs
check_freqs(causalSNPs, nama_freqs)

genFixed <- geneticFixedEffects(N = 5000, 
                                P = 1, 
                                X_causal = causalSNPs)

#Write out copy of causalSNPs
writeLines(colnames(causalSNPs), paste("/home/gerald/Documents/PhD/papers/paper4/causalSNPs_3_nama.txt", sep = ""))

#anc = read.big.matrix("/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_ancs_nama_5000_bt.csv",sep = ",", type="char", backingfile = "nama_ancs_nama_5000_bt.bin", backingpath="/home/gerald/Documents/PhD/papers/paper4/true_anc/", descriptorfile = "nama_ancs_nama_5000_bt.desc", shared = T)

# rescale genetic phenotype components
#Allelic
phen_no = 1
# Total variance proportions have to add up to 1
total <- round(shared * h2s *genVar, digits = 3) +  round(independent * h2s * genVar, digits = 3) +
  round(shared * (1-h2s) * genVar, digits = 3) +   round(independent * (1-h2s) * genVar, digits = 3) +
  round(shared * phi* noiseVar, digits = 3) +  round(independent * phi* noiseVar, digits = 3) +
  #rho * noiseVar +
  round(shared * delta * noiseVar, digits = 3) +  round(independent * delta * noiseVar, digits = 3)

total == 1
#> [1] TRUE

genBg_shared_scaled <- rescaleVariance(genBg$shared, shared * (1-h2s) *genVar)
genBg_independent_scaled <- rescaleVariance(genBg$independent, 
                                            independent * (1-h2s) * genVar)
genFixed_shared_scaled <- rescaleVariance(genFixed$shared, shared * h2s *genVar)
genFixed_independent_scaled <- rescaleVariance(genFixed$independent, 
                                               independent * h2s *genVar)
ancFixed_shared_scaled <- 0
ancFixed_independent_scaled <- 0

# combine components into final phenotype
t_allele <- scale(noiseBg_shared_scaled$component[,phen_no] +
           noiseBg_independent_scaled$component[,phen_no] +
           noiseFixed_shared_scaled$component[,phen_no] +
           #noiseFixed_independent_scaled$component[,phen_no] + 
           ancFixed_shared_scaled +
           ancFixed_independent_scaled +
           genBg_shared_scaled$component[,phen_no] +
           genBg_independent_scaled$component[,phen_no] + 
           genFixed_shared_scaled$component + 
           genFixed_independent_scaled$component)

colnames(t_allele) = "t_allele"

#nama = attach.big.matrix("/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_ancs_nama_5000_bt.desc", type="double")
#gbr = attach.big.matrix("/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_ancs_gbr_5000_bt.desc", type="double")
#ep = attach.big.matrix("/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_ancs_ep_5000_bt.desc", type="double")
#ancFixed_nama = geneticFixedEffects(N = 5000, 
#                                    P = 1, 
#                                    X_causal = nama[-1,])
#gc()
#save(ancFixed_nama, file = "/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_ancs_nama_5000_fixed.RData")
#rm(ancFixed_nama)
#gc()
#ancFixed_gbr = geneticFixedEffects(N = 5000, 
#                                    P = 1, 
#                                    X_causal = gbr[-1,])
#gc()
#save(ancFixed_gbr, file = "/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_ancs_gbr_5000_fixed.RData")
#rm(ancFixed_gbr)
#gc()
#ancFixed_ep = geneticFixedEffects(N = 5000, 
#                                    P = 1, 
#                                    X_causal = ep[-1,])
#gc()
#save(ancFixed_ep, file = "/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_ancs_ep_5000_fixed.RData")
#rm(ancFixed_ep)
#gc()

#Ancestry
# Total variance proportions have to add up to 1
total <- round(shared * lambda *genVar, digits = 3) +  round(independent * lambda  * genVar, digits = 3) +
  round(shared * (1-lambda ) *genVar, digits = 3) +  round(independent * (1-lambda ) * genVar, digits = 3) +
  round(shared * phi* noiseVar, digits = 3) +  round(independent * phi* noiseVar, digits = 3) +
  #rho * noiseVar +
  round(shared * delta * noiseVar, digits = 3) +  round(independent * delta * noiseVar, digits = 3)

total == 1
#> [1] TRUE

genBg_shared_scaled <- 0
genBg_independent_scaled <- 0



causal_IDs = which(colnames(nama_genotypes$genotypes) %in% colnames(causalSNPs))

anc_files = c("/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_ancs_nama_5000_bt.csv", "/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_ancs_gbr_5000_bt.csv", "/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_ancs_ep_5000_bt.csv")

lapply(anc_files, function(a){
  lapply(causal_IDs, function(e){
    system(paste("sed -n ", e, "p ", a, sep=""), intern = T) %>%
      strsplit(., split = ",") %>%
      unlist() %>%
      .[-1] %>%
      as.numeric()
  }) %>%
    do.call("cbind", .) ->
    anc
  
  rownames(anc) = rownames(causalSNPs)
  geneticFixedEffects(N = 5000,
                      P = 1,
                      X_causal = anc)
}) ->
  anc

anc_names = c("t_anc_nama","t_anc_gbr","t_anc_ep")


loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

t_anc = c()
#Nama:4, MSL:6, GBR:7 EP:8, EAS:10, SAS:11

lapply(1:3, function(p){
  ancFixed = anc[[p]]
  ancFixed_shared_scaled <- rescaleVariance(ancFixed$shared, shared * (1-h2s) *genVar)
  ancFixed_independent_scaled <- rescaleVariance(ancFixed$independent, 
                                               independent * (1-h2s) *genVar)
  # combine components into final phenotype
  res = scale(noiseBg_shared_scaled$component[,phen_no] +
                      noiseBg_independent_scaled$component[,phen_no] +
                      noiseFixed_shared_scaled$component[,phen_no] +
                      #noiseFixed_independent_scaled$component[,phen_no] + 
                      ancFixed_shared_scaled$component[,phen_no] +
                      ancFixed_independent_scaled$component[,phen_no] +
                      genBg_shared_scaled +
                      genBg_independent_scaled + 
                      genFixed_shared_scaled$component + 
                      genFixed_independent_scaled$component)
  colnames(res) = anc_names[p]
  res
}) -> t_anc

#Ancestry + Allelic
# Total variance proportions have to add up to 1
total <- round(shared * h2s *genVar, digits = 3) +  round(independent * h2s * genVar, digits = 3) +
  round(shared * lambda *genVar, digits = 3) +  round(independent * lambda * genVar, digits = 3) +
  round(shared * (1-h2s-lambda) * genVar, digits = 3) +   round(independent * (1-h2s-lambda) * genVar, digits = 3) +
  round(shared * phi* noiseVar, digits = 3) +  round(independent * phi* noiseVar, digits = 3) +
  #rho * noiseVar +
  round(shared * delta * noiseVar, digits = 3) +  round(independent * delta * noiseVar, digits = 3)

total == 1
#> [1] TRUE

genBg_shared_scaled <- rescaleVariance(genBg$shared, shared * (1-h2s-lambda) *genVar)
genBg_independent_scaled <- rescaleVariance(genBg$independent, 
                                            independent * (1-h2s-lambda) * genVar)
                                            
anc_names = c("t_apa_nama","t_apa_gbr","t_apa_ep")                                            
t_allele_anc = c()
#Nama:4, MSL:6, GBR:7 EP:8, EAS:10, SAS:11
lapply(1:3, function(p){
  ancFixed = anc[[p]]
  
  ancFixed_shared_scaled <- rescaleVariance(ancFixed$shared, shared * lambda *genVar)
  ancFixed_independent_scaled <- rescaleVariance(ancFixed$independent,independent * lambda *genVar)
  
  # combine components into final phenotype
  res = scale(noiseBg_shared_scaled$component[,phen_no] +
              noiseBg_independent_scaled$component[,phen_no] +
              noiseFixed_shared_scaled$component[,phen_no] +
              #noiseFixed_independent_scaled$component[,phen_no] + 
              ancFixed_shared_scaled$component[,phen_no] +
              ancFixed_independent_scaled$component[,phen_no] +
              genBg_shared_scaled$component[,phen_no] +
              genBg_independent_scaled$component[,phen_no] + 
              genFixed_shared_scaled$component + 
              genFixed_independent_scaled$component)
  colnames(res) = anc_names[p]
  res
}) -> t_apa

#LAAA
# Total variance proportions have to add up to 1
total <- round(shared * h2s *genVar, digits = 3) +  round(independent * h2s * genVar, digits = 3) +
  round(shared * lambda *genVar, digits = 3) +  round(independent * lambda * genVar, digits = 3) +
  round(shared * psy *genVar, digits = 3) +  round(independent * psy * genVar, digits = 3) +
  round(shared * (1-h2s-lambda-psy) * genVar, digits = 3) +   round(independent * (1-h2s-lambda-psy) * genVar, digits = 3) +
  round(shared * phi* noiseVar, digits = 3) +  round(independent * phi* noiseVar, digits = 3) +
  #rho * noiseVar +
  round(shared * delta * noiseVar, digits = 3) +  round(independent * delta * noiseVar, digits = 3)

total == 1
#> [1] TRUE

genBg_shared_scaled <- rescaleVariance(genBg$shared, shared * (1-h2s-lambda-psy) *genVar)
genBg_independent_scaled <- rescaleVariance(genBg$independent, 
                                            independent * (1-h2s-lambda-psy) * genVar)
  
int_files = c("/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_laaa_nama_5000_bt.desc", "/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_laaa_gbr_5000_bt.desc", "/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_laaa_ep_5000_bt.desc")

lapply(int_files, function(a){
  anc = attach.big.matrix(a)
  
  selected = anc[,causal_IDs]
  
  rownames(selected) = rownames(causalSNPs)
  selected
  
  geneticFixedEffects(N = 5000,
                      P = 1,
                      X_causal = selected)
}) ->
  ints

anc_names = c("t_laaa_nama","t_laaa_gbr","t_laaa_ep")                                            
t_allele_anc = c()
#Nama:4, MSL:6, GBR:7 EP:8, EAS:10, SAS:11
lapply(1:3, function(p){
  ancFixed = anc[[p]]
  intFixed = ints[[p]]
  
  ancFixed_shared_scaled <- rescaleVariance(ancFixed$shared, shared * lambda *genVar)
  ancFixed_independent_scaled <- rescaleVariance(ancFixed$independent,independent * lambda *genVar)
  intFixed_shared_scaled <- rescaleVariance(intFixed$shared, shared * psy *genVar)
  intFixed_independent_scaled <- rescaleVariance(intFixed$independent,independent * psy *genVar)
  
  # combine components into final phenotype
  res = scale(noiseBg_shared_scaled$component[,phen_no] +
              noiseBg_independent_scaled$component[,phen_no] +
              noiseFixed_shared_scaled$component[,phen_no] +
              #noiseFixed_independent_scaled$component[,phen_no] + 
              ancFixed_shared_scaled$component[,phen_no] +
              ancFixed_independent_scaled$component[,phen_no] +
              intFixed_shared_scaled$component[,phen_no] +
              intFixed_independent_scaled$component[,phen_no] +
              genBg_shared_scaled$component[,phen_no] +
              genBg_independent_scaled$component[,phen_no] + 
              genFixed_shared_scaled$component + 
              genFixed_independent_scaled$component)
  colnames(res) = anc_names[p]
  res
}) -> t_laaa
gc()

                                                 
phenotypes = cbind(t_allele, t_anc[[1]], t_apa[[1]], t_laaa[[1]], t_anc[[2]], t_apa[[2]], t_laaa[[2]], t_anc[[3]], t_apa[[3]], t_laaa[[3]])

glob_ancs = read.csv("/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_global_ancs_5000_ALL.csv", header=T)

age_sex = noiseFixed$cov
age_sex[,2] = as.integer(age_sex[,2])
phen_age_sex = data.frame("age" = age_sex[,1], "sex" = age_sex[,2], phenotypes,  glob_ancs[,c("NAMA", "GBR", "EP")])
write.table(phen_age_sex, "/home/gerald/Documents/PhD/papers/paper4/true_anc/nama_5000_phenotypes_3.tsv", quote = F, sep = "\t", row.names = T, col.names = T)








