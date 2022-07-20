library(PhenotypeSimulator)

#Read in data
genotypefile = "./multi_sim_5000_all_maf05_biallelic_snps_ids_filtered"
data <- snpStats::read.plink(bed=genotypefile)

sds = sapply(1:ncol(data$genotypes),  function(e) {
  snp = as(data$genotypes[,e], 'numeric')
  af = getAlleleFrequencies(snp)
  var_geno = sqrt(2*af[1]*af[2])
  var_geno[var_geno == 0] <- 1
  geno_sd = (snp - 2*af[1]) / var_geno
})

rm(data)
gc()

kinship <- getKinship(N=15000, X=sds, verbose = T)
rm(sds)
gc()
save(kinship, file = "kinship_3e5.RData")
load(file = "../kinship_3e5.RData")
genBg <- geneticBgEffects(N=15000, kinship = kinship, P = 1, shared = F)
save(genBg, file = "bg_p10_3e5.RData")


rm(kinship)
gc()

