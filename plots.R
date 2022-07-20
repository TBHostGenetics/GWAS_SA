library(magrittr)
library(vcfR)
library(ldsep)
library(bigmemory)
library(doParallel)

bin_it = function(input, bin_size){
  block_start = 1
  blocks = c()
  blocks[[1]] = input[1]
  curr_block = 1
  for(h in 2:length(input)){
    if(input[block_start] + bin_size > input[h] &
       input[h] > input[block_start]){
      blocks[[curr_block]] = c(blocks[[curr_block]], input[h])
    } else{
      block_start = h
      curr_block = curr_block + 1
      blocks[[curr_block]] = input[block_start]
    }
  }
  return(blocks)
}

input = t_sac_c1
causal_set = causal_sac[[1]]
p_val = 1e-8
pops = s_pops

sapply(t_sac_c1, function(f){
  sapply(f, function(e){
    e[,"position"] %>% duplicated() %>% sum()
  })
})

chr_start = c(1)
input = t_nama_c1
for(i in 2:nrow(input$laaa$nama)){
  if(input$laaa$nama[i,"position"] < input$laaa$nama[i-1,"position"]){
    chr_start = c(chr_start, i)
  }
}

analyse = function(input, causal_set, c_ld, p_val = 1e-8, pops){
  
  traits = c("laaa", "apa", "anc", "allele")
  
  input$laaa$nama %>% names() %>% grep("p_val", .) -> p_cols
  input$laaa$nama[,p_cols] %>% names() %>% gsub("_p_val", "", .) -> mod_names
  
  tmp = c()
  
  for(trait in traits){
    for(pop in pops){
      #Get hits
      lapply(p_cols, function(p){
        input[[trait]][[pop]][which(input[[trait]][[pop]][,p] < p_val),1]
      }) -> 
        hits
      
      #Get f_pos by filtering hits for SNPs in LD
      lapply(hits, function(h){
        h[! h %in% unlist(c_ld)]
      }) ->
        hits_no_ld
      
      #Check causal
      sapply(causal_set$pos, function(e){
        sapply(p_cols, function(col){
          which(input[[trait]][[pop]][,col] < p_val) %>%
            input[[trait]][[pop]][.,"position"] %>%
            grepl(e, .) %>%
            sum()
        })
      }) ->
        tmp[[trait]][[pop]]
      
      tmp[[trait]][[pop]][which(tmp[[trait]][[pop]] != 0)] = 1
      
      rownames(tmp[[trait]][[pop]]) = mod_names
      colnames(tmp[[trait]][[pop]]) = paste(causal_set[,1],causal_set[,2], sep = "_")
      
      t_pos = rowSums(tmp[[trait]][[pop]])
      hits = sapply(hits, length)
      ld = hits - sapply(hits_no_ld, length)
      f_pos = hits- t_pos - ld
      f_neg = 10-t_pos
      
      tmp[[trait]][[pop]] = cbind(hits,
                                  t_pos,
                                  ld,
                                  f_pos,
                                  f_neg)
    }
  }
  
  return(tmp)
  
}

#Get all false positive hit postions
get_hits = function(input, p_val, causal_set, pops){
  hit_pos = c()
  traits = c("laaa", "apa", "anc", "allele")
  t_nama_c1$laaa$nama %>% names() %>% grep("p_val", .) -> p_cols
  t_nama_c1$laaa$nama[,p_cols] %>% names() %>% gsub("_p_val", "", .) -> mod_names
  hits = c()
  
  for(trait in traits){
    for(pop in pops){
      #Get hits
      lapply(p_cols, function(p){
        input[[trait]][[pop]][which(input[[trait]][[pop]][,p] < p_val),1]
      }) -> 
        hits[[trait]][[pop]]
      
      combined_hits = unlist(hits[[trait]][[pop]])
      
      #Get pos of causal hits
      #Exclude causal hits from hits
      #Get false positive hit postions
      combined_hits[-which(combined_hits %in% causal_set$pos) ] %>%
      sapply(., function(hit){
        which(pos == hit)
      }) %>% 
        unlist() %>%
        unique() %>%
        c(hit_pos, .) ->
        hit_pos
    }
  }
  return(list(hits, hit_pos))
}

path = "/home/gerald/Documents/PhD/papers/paper4/"

causal_nama = system(paste("ls ", path, "results/causalSNPs_*_nama.txt", sep = ""), intern = T)
lapply(causal_nama, function(e){
  read.delim(e, header = F, sep = " ", col.names = c("chr", "pos"))
}) ->
  causal_nama

causal_sac = system(paste("ls ", path, "results/causalSNPs_*_sac.txt", sep = ""), intern = T)
lapply(causal_sac, function(e){
  read.delim(e, header = F, sep = " ", col.names = c("chr", "pos"))
}) ->
  causal_sac

pos = read.delim("/home/gerald/Documents/PhD/papers/paper4/data/sim_alleles.txt", header = T)[,1]

nama.bm = attach.big.matrix("/home/gerald/Documents/PhD/papers/paper4/data/sim_allele_dose_nofrq.desc")

# t_nama_c1 = c()
# t_nama_c2 = c()
# t_nama_c3 = c()
# 
# i_nama_c1 = c()
# i_nama_c2 = c()
# i_nama_c3 = c()
# 
# t_sac_c1 = c()
# t_sac_c2 = c()
# t_sac_c3 = c()
# 
# i_sac_c1 = c()
# i_sac_c2 = c()
# i_sac_c3 = c()

traits = c("laaa", "apa", "anc", "allele")
n_pops = c("nama", "ep", "eur")
s_pops = c("nama", "msl", "gbr", "eas", "sas")


# for(trait in traits){
#   for(pop in n_pops){
#     t_nama_c1[[trait]][[pop]] = read.delim(paste(path, "results/true/nama/trim_c1_", trait, "_", pop, ".tsv", sep = ""))
#   }
# }
# 
# for(trait in traits){
#   for(pop in n_pops){
#     t_nama_c2[[trait]][[pop]] = read.delim(paste(path, "results/true/nama/trim_c2_", trait, "_", pop, ".tsv", sep = ""))
#   }
# }
# 
# for(trait in traits){
#   for(pop in n_pops){
#     t_nama_c3[[trait]][[pop]] = read.delim(paste(path, "results/true/nama/trim_c3_", trait, "_", pop, ".tsv", sep = ""))
#   }
# }
# 
# 
# for(trait in traits){
#   for(pop in n_pops){
#     i_nama_c1[[trait]][[pop]] = read.delim(paste(path, "results/inf/nama/trim_c1_", trait, "_", pop, ".txt", sep = ""))
#   }
# }
# 
# for(trait in traits){
#   for(pop in n_pops){
#     i_nama_c2[[trait]][[pop]] = read.delim(paste(path, "results/inf/nama/trim_c2_", trait, "_", pop, ".txt", sep = ""))
#   }
# }
# 
# for(trait in traits){
#   for(pop in n_pops){
#     i_nama_c3[[trait]][[pop]] = read.delim(paste(path, "results/inf/nama/trim_c3_", trait, "_", pop, ".txt", sep = ""))
#   }
# }
# 
# 
# for(trait in traits){
#   for(pop in s_pops){
#     t_sac_c1[[trait]][[pop]] = read.delim(paste(path, "results/true/sac/filtered_", trait, "_", pop, "_1.txt", sep = ""))
#   }
# }
# 
# for(trait in traits){
#   for(pop in s_pops){
#     t_sac_c2[[trait]][[pop]] = read.delim(paste(path, "results/true/sac/filtered_", trait, "_", pop, "_2.txt", sep = ""))
#   }
# }
# 
# for(trait in traits){
#   for(pop in s_pops){
#     t_sac_c3[[trait]][[pop]] = read.delim(paste(path, "results/true/sac/filtered_", trait, "_", pop, "_3.txt", sep = ""))
#   }
# }
# 
# 
# for(trait in traits){
#   for(pop in s_pops){
#     i_sac_c1[[trait]][[pop]] = read.delim(paste(path, "results/inf/sac/filtered_", trait, "_", pop, "_1.txt", sep = ""))
#   }
# }
# 
# for(trait in traits){
#   for(pop in s_pops){
#     i_sac_c2[[trait]][[pop]] = read.delim(paste(path, "results/inf/sac/filtered_", trait, "_", pop, "_2.txt", sep = ""))
#   }
# }
# 
# for(trait in traits){
#   for(pop in s_pops){
#     i_sac_c3[[trait]][[pop]] = read.delim(paste(path, "results/inf/sac/filtered_", trait, "_", pop, "_3.txt", sep = ""))
#   }
# }


# hits_t_nama_1 = get_hits(t_nama_c1, 1e-8, causal_nama[[1]], n_pops)
# hits_t_nama_2 = get_hits(t_nama_c2, 1e-8, causal_nama[[2]], n_pops)
# hits_t_nama_3 = get_hits(t_nama_c3, 1e-8, causal_nama[[3]], n_pops)
# 
# hits_i_nama_1 = get_hits(i_nama_c1, 1e-8, causal_nama[[1]], n_pops)
# hits_i_nama_2 = get_hits(i_nama_c2, 1e-8, causal_nama[[2]], n_pops)
# hits_i_nama_3 = get_hits(i_nama_c3, 1e-8, causal_nama[[3]], n_pops)
# 
# hits_t_sac_1 = get_hits(t_sac_c1, 1e-8, causal_sac[[1]], s_pops)
# hits_t_sac_2 = get_hits(t_sac_c2, 1e-8, causal_sac[[2]], s_pops)
# hits_t_sac_3 = get_hits(t_sac_c3, 1e-8, causal_sac[[3]], s_pops)
# 
# hits_i_sac_1 = get_hits(i_sac_c1, 1e-8, causal_sac[[1]], s_pops)
# hits_i_sac_2 = get_hits(i_sac_c2, 1e-8, causal_sac[[2]], s_pops)
# hits_i_sac_3 = get_hits(i_sac_c3, 1e-8, causal_sac[[3]], s_pops)
# 
# hits_pos_t_nama_1 = hits_t_nama_1[[2]]
# hits_pos_t_nama_2 = hits_t_nama_2[[2]]
# hits_pos_t_nama_3 = hits_t_nama_3[[2]]
# hits_pos_i_nama_1 = hits_i_nama_1[[2]]
# hits_pos_i_nama_2 = hits_i_nama_2[[2]]
# hits_pos_i_nama_3 = hits_i_nama_3[[2]]
# hits_pos_t_sac_1 = hits_t_sac_1[[2]]
# hits_pos_t_sac_2 = hits_t_sac_2[[2]]
# hits_pos_t_sac_3 = hits_t_sac_3[[2]]
# hits_pos_i_sac_1 = hits_i_sac_1[[2]]
# hits_pos_i_sac_2 = hits_i_sac_2[[2]]
# hits_pos_i_sac_3 = hits_i_sac_3[[2]]

# save(list = c("hits_pos_t_nama_1", "hits_pos_t_nama_2", "hits_pos_t_nama_3",
#               "hits_pos_i_nama_1", "hits_pos_i_nama_2", "hits_pos_i_nama_3",
#               "hits_pos_t_sac_1", "hits_pos_t_sac_2", "hits_pos_t_sac_3",
#               "hits_pos_i_sac_1", "hits_pos_i_sac_2", "hits_pos_i_sac_3",
#               "hits_t_nama_1", "hits_t_nama_2", "hits_t_nama_3",
#               "hits_i_nama_1", "hits_i_nama_2", "hits_i_nama_3",
#               "hits_t_sac_1", "hits_t_sac_2", "hits_t_sac_3",
#               "hits_i_sac_1", "hits_i_sac_2", "hits_i_sac_3",
#               "t_nama_c1","t_nama_c2","t_nama_c3",
#               "i_nama_c1","i_nama_c2","i_nama_c3",
#               "t_sac_c1","t_sac_c2","t_sac_c3",
#               "i_sac_c1","i_sac_c2","i_sac_c3"),
# file = "/home/gerald/Documents/PhD/papers/paper4/results/raw_and_hits.RData")

load(file = "/home/gerald/Documents/PhD/papers/paper4/results/raw_and_hits.RData")


c1_nama_f_pos_hits = unique(c(hits_pos_t_nama_1,hits_pos_i_nama_1))
c2_nama_f_pos_hits = unique(c(hits_pos_t_nama_2,hits_pos_i_nama_2))
c3_nama_f_pos_hits = unique(c(hits_pos_t_nama_3,hits_pos_i_nama_3))

c1_sac_f_pos_hits = unique(c(hits_pos_t_sac_1,hits_pos_i_sac_1))
c2_sac_f_pos_hits = unique(c(hits_pos_t_sac_2,hits_pos_i_sac_2))
c3_sac_f_pos_hits = unique(c(hits_pos_t_sac_3,hits_pos_i_sac_3))


#registerDoParallel(4)
#for each false pos hit, get LD scores for all causal SNPs
# LD_nama_f_pos_c1 = foreach(c_pos=which(pos %in% causal_nama[[1]]$pos)) %dopar% {
#   lapply(unique(c1_nama_f_pos_hits), function(h){
#     c(h, mldest(geno = nama.bm[c(h, c_pos),], K = 2, type = "comp")$r2)
#   }) %>%
#     do.call("rbind", .)
# }
# LD_nama_f_pos_c2 = foreach(c_pos=which(pos %in% causal_nama[[2]]$pos)) %dopar% {
#   lapply(unique(c2_nama_f_pos_hits), function(h){
#     c(h, mldest(geno = nama.bm[c(h, c_pos),], K = 2, type = "comp")$r2)
#   }) %>%
#     do.call("rbind", .)
# }
# LD_nama_f_pos_c3 = foreach(c_pos=which(pos %in% causal_nama[[3]]$pos)) %dopar% {
#   lapply(unique(c3_nama_f_pos_hits), function(h){
#     c(h, mldest(geno = nama.bm[c(h, c_pos),], K = 2, type = "comp")$r2)
#   }) %>%
#     do.call("rbind", .)
# }
# 
# LD_sac_f_pos_c1 = foreach(c_pos=which(pos %in% causal_sac[[1]]$pos)) %dopar% {
#   lapply(unique(c1_sac_f_pos_hits), function(h){
#     c(h, mldest(geno = nama.bm[c(h, c_pos),], K = 2, type = "comp")$r2)
#   }) %>%
#     do.call("rbind", .)
# }
# 
# LD_sac_f_pos_c2 = foreach(c_pos=which(pos %in% causal_sac[[2]]$pos)) %dopar% {
#   lapply(unique(c2_sac_f_pos_hits), function(h){
#     c(h, mldest(geno = nama.bm[c(h, c_pos),], K = 2, type = "comp")$r2)
#   }) %>%
#     do.call("rbind", .)
# }
# 
# LD_sac_f_pos_c3 = foreach(c_pos=which(pos %in% causal_sac[[3]]$pos)) %dopar% {
#   lapply(unique(c3_sac_f_pos_hits), function(h){
#     c(h, mldest(geno = nama.bm[c(h, c_pos),], K = 2, type = "comp")$r2)
#   }) %>%
#     do.call("rbind", .)
# }

# save(list = c("LD_nama_f_pos_c1", "LD_nama_f_pos_c2", "LD_nama_f_pos_c3", 
#        "LD_sac_f_pos_c1", "LD_sac_f_pos_c2", "LD_sac_f_pos_c3"), 
#      file = paste(path, "results/LD_f_pos_all.RData", sep=""))

load(file = paste(path, "results/LD_f_pos_all.RData", sep=""))
#load(file = paste(path, "results/LD_f_pos_sac.RData", sep=""))

#Get positions of sites for each causal SNP with LD > 0.2
causal_ld = c()
lapply(LD_nama_f_pos_c1, function(c_ld){
  pos[c1_nama_f_pos_hits[which(c_ld[,2] > 0.1)]]
}) -> causal_ld[["nama"]][["c1"]]
lapply(LD_nama_f_pos_c2, function(c_ld){
  pos[c2_nama_f_pos_hits[which(c_ld[,2] > 0.1)]]
}) -> causal_ld[["nama"]][["c2"]]
lapply(LD_nama_f_pos_c3, function(c_ld){
  pos[c3_nama_f_pos_hits[which(c_ld[,2] > 0.1)]]
}) -> causal_ld[["nama"]][["c3"]]
lapply(LD_sac_f_pos_c1, function(c_ld){
  pos[c1_sac_f_pos_hits[which(c_ld[,2] > 0.1)]]
}) -> causal_ld[["sac"]][["c1"]]
lapply(LD_sac_f_pos_c2, function(c_ld){
  pos[c2_sac_f_pos_hits[which(c_ld[,2] > 0.1)]]
}) -> causal_ld[["sac"]][["c2"]]
lapply(LD_sac_f_pos_c3, function(c_ld){
  pos[c3_sac_f_pos_hits[which(c_ld[,2] > 0.1)]]
}) -> causal_ld[["sac"]][["c3"]]


res = c()
res[["nama"]][["true"]][["c1"]] = analyse(t_nama_c1, causal_nama[[1]], causal_ld$nama$c1, 5e-8, n_pops)
res[["nama"]][["true"]][["c2"]] = analyse(t_nama_c2, causal_nama[[2]], causal_ld$nama$c2, 5e-8, n_pops)
res[["nama"]][["true"]][["c3"]] = analyse(t_nama_c3, causal_nama[[3]], causal_ld$nama$c3, 5e-8, n_pops)
#res[["nama"]][["inf"]][["c1"]]$allele$nama = res[["nama"]][["inf"]][["c1"]]$allele$ep
res[["nama"]][["inf"]][["c1"]] = analyse(i_nama_c1, causal_nama[[1]], causal_ld$nama$c1, 5e-8, n_pops)
res[["nama"]][["inf"]][["c2"]] = analyse(i_nama_c2, causal_nama[[2]], causal_ld$nama$c2, 5e-8, n_pops)
res[["nama"]][["inf"]][["c3"]] = analyse(i_nama_c3, causal_nama[[3]], causal_ld$nama$c3, 5e-8, n_pops)
res[["sac"]][["true"]][["c1"]] = analyse(t_sac_c1, causal_sac[[1]], causal_ld$sac$c1, 5e-8, s_pops)
res[["sac"]][["true"]][["c2"]] = analyse(t_sac_c2, causal_sac[[2]], causal_ld$sac$c2, 5e-8, s_pops)
res[["sac"]][["true"]][["c3"]] = analyse(t_sac_c3, causal_sac[[3]], causal_ld$sac$c3, 5e-8, s_pops)
res[["sac"]][["inf"]][["c1"]] = analyse(i_sac_c1, causal_sac[[1]], causal_ld$sac$c1, 5e-8, s_pops)
res[["sac"]][["inf"]][["c2"]] = analyse(i_sac_c2, causal_sac[[2]], causal_ld$sac$c2, 5e-8, s_pops)
res[["sac"]][["inf"]][["c3"]] = analyse(i_sac_c3, causal_sac[[3]], causal_ld$sac$c3, 5e-8, s_pops)



library(purrr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(viridis)
  
plot_stuff = function(input, target, plot_label, pop, excluder = NULL){
  
  if(length(target) > 1){
    names(data.frame(input)) %>% grep(target[1], .) ->
      df_nama_cols1
    data.frame(input)[,df_nama_cols1] %>%
      melt(.) ->
      plot_data1
    
    names(data.frame(input)) %>% grep(target[2], .) ->
      df_nama_cols2
    data.frame(input)[,df_nama_cols2] %>%
      melt(.) ->
      plot_data2
    
    plot_data = data.frame("variable" = plot_data1[,1], "value" = plot_data2[,2] / (plot_data2[,2]+plot_data1[,2]))
    
  } else{
    names(data.frame(input)) %>% grep(target, .) ->
      df_nama_cols
    data.frame(input)[,df_nama_cols] %>%
      melt(.) ->
      plot_data
  }
  
  if(!is.null(excluder)){
    plot_data = plot_data[grep(excluder, plot_data[,1], invert = T),]
  }
  
  if(pop == "nama"){
    reps = 4*(nrow(plot_data)/5)
  } else{
    reps = 4*(nrow(plot_data)/5)
  }
  
  model_names = rep(row.names(data.frame(input)[,df_nama_cols]), reps)
  model_names = gsub("null", "neutral", model_names)
  
  plot_data = cbind(plot_data, "models" = model_names)
  
  ggplot() + geom_col(aes(x = variable, y = value, fill = as.factor(models)), data = plot_data, size = 1, position = "dodge") +
    labs(title = plot_label, y = "Hits", x = "Phenotypes")
}

get_avs = function(input, var){
  result = c()
  for(trait in traits){
    for(pop in names(input[[1]][[1]])){
      result[[trait]][[pop]] = unlist(Map(function(x,y,z) mean(c(x,y,z)), input$c1[[trait]][[pop]][,var], input$c2[[trait]][[pop]][,var], input$c3[[trait]][[pop]][,var]))
    }
  }
  return(result)
}

combine_avs = function(input) {
  get_avs(input, "t_pos") %>% 
    data.frame() %>% 
    cbind(c("Standard", "LAAA", "GA", "LA", "APA"), .) %>% 
    melt() -> 
    t_pos
  
  colnames(t_pos) = c("model", "variable", "t_pos")
  
  get_avs(input, "f_pos") %>% 
    data.frame() %>% 
    melt() %>% 
    .$value ->
    f_pos
  
  exclude = which(!grepl(pattern = "anc", x = t_pos$variable))
  
  return(data.frame(t_pos, "f_pos" = f_pos)[exclude,])
}

combined = c()
combined[["nama"]][["inf"]] = combine_avs(res$nama$inf)
combined[["nama"]][["true"]] = combine_avs(res$nama$true)

combined[["sac"]][["inf"]] = combine_avs(res$sac$inf)
combined[["sac"]][["true"]] = combine_avs(res$sac$true)


plot_combined = function(data, pop, var, label = ""){
  if(pop == "nama"){
    x_lab = c("LAAA-Nama","LAAA-LWK","LAAA-GBR",
              "APA-Nama","APA-LWK","APA-GBR",
              #"LAO-Nama","LAO-EP","LAO-GBR",
              "AO-Nama","AO-LWK","AO-GBR")
  } else{
    x_lab = c("LAAA-Nama","LAAA-MSL","LAAA-GBR","LAAA-CHB","LAAA-GIH",
              "APA-Nama","APA-MSL","APA-GBR","APA-CHB","APA-GIH",
              #"LAO-Nama","LAO-MSL","LAO-GBR","LAO-CHB","LAO-GIH",
              "AO-Nama","AO-MSL","AO-GBR","AO-CHB","AO-GIH")
  }
  
  if(var == "t_pos"){
    breaks = 1:10
    maxlim = 10
  } else{
    breaks = seq(250,2500,250)
    maxlim = 2500 
  }
  
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
            "#D55E00", "#CC79A7")
  
  ggplot() + geom_col(aes(x = data[,var], y = variable, fill = as.factor(model)), data = data, size = 1, width = 0.8, position = position_dodge(width = 0.8)) +
    labs(title = label, x = "Hits", y = "Phenotypes") + 
    scale_fill_discrete(name = "GWAS model", type = cbp1) +
    coord_cartesian(xlim=c(0,maxlim)) + 
    scale_x_continuous(breaks = breaks) + 
    scale_y_discrete(labels = x_lab)+ 
    theme(axis.text = element_text(size=14),
          axis.text.x = element_text(angle = 0, margin = margin(10,0,0,0)),
          axis.title.x = element_text(size=20, margin = margin(10,0,0,0)),
          axis.title.y = element_text(size=20, margin = margin(0,10,0,0)),
          plot.margin = unit(c(1,1,1,1), units = "cm"),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          panel.grid.major.y = element_line(size = 0, linetype = 'solid',
                                            colour = "grey"),
          panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "grey"))
}


(plot_combined(combined$nama$inf, "nama", "t_pos", "") /
    plot_combined(combined$nama$inf, "nama", "f_pos", "")) +
  plot_layout(guides = 'collect')
    
(plot_combined(combined$nama$true, "nama", "t_pos", "") /
    plot_combined(combined$nama$true, "nama", "f_pos", "")) +
  plot_layout(guides = 'collect')


(plot_combined(combined$sac$inf, "sac", "t_pos", "") /
    plot_combined(combined$sac$inf, "sac", "f_pos", "")) +
  plot_layout(guides = 'collect')

(plot_combined(combined$sac$true, "sac", "t_pos", "") /
    plot_combined(combined$sac$true, "sac", "f_pos", "")) +
  plot_layout(guides = 'collect')


png(filename = "/home/gerald/Documents/PhD/papers/paper4/nama_inf.png", units = "cm", width = 46, height = 32, res = 300)

(plot_combined(combined$nama$inf, "nama", "t_pos", "") |
    plot_combined(combined$nama$inf, "nama", "f_pos", "")) +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

dev.off()

png(filename = "/home/gerald/Documents/PhD/papers/paper4/nama_true.png", units = "cm", width = 46, height = 32, res = 300)

  (plot_combined(combined$nama$true, "nama", "t_pos", "") |
       plot_combined(combined$nama$true, "nama", "f_pos", "")) +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

dev.off()


png(filename = "/home/gerald/Documents/PhD/papers/paper4/sac_true.png", units = "cm", width = 46, height = 32, res = 300)

(plot_combined(combined$sac$inf, "sac", "t_pos", "") |
    plot_combined(combined$sac$inf, "sac", "f_pos", "")) +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

dev.off()

png(filename = "/home/gerald/Documents/PhD/papers/paper4/sac_inf.png", units = "cm", width = 46, height = 32, res = 300)

(plot_combined(combined$sac$true, "sac", "t_pos", "") |
   plot_combined(combined$sac$true, "sac", "f_pos", "")) +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

dev.off()


n_inf = c()

n_inf[["t"]][[1]] = plot_stuff(res$nama$inf$c1, "t_pos", "", "nama") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
n_inf[["t"]][[2]] = plot_stuff(res$nama$inf$c2, "t_pos", "",  "nama") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
n_inf[["t"]][[3]] = plot_stuff(res$nama$inf$c3, "t_pos", "",  "nama") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))

n_inf[["f"]][[1]] = plot_stuff(res$nama$inf$c1, "f_pos", "", "nama") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
n_inf[["f"]][[2]] = plot_stuff(res$nama$inf$c2, "f_pos", "",  "nama") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
n_inf[["f"]][[3]] = plot_stuff(res$nama$inf$c3, "f_pos", "",  "nama") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))

# n_inf[["t/f"]][[1]] = plot_stuff(res$nama$inf$c1, c("t_pos","f_pos"), "", "nama") + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
# n_inf[["t/f"]][[2]] = plot_stuff(res$nama$inf$c2, c("t_pos","f_pos"), "",  "nama") + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
# n_inf[["t/f"]][[3]] = plot_stuff(res$nama$inf$c3, c("t_pos","f_pos"), "",  "nama") + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))


(n_inf[["t"]][[1]] | n_inf[["t"]][[2]] | n_inf[["t"]][[3]]) / 
  (n_inf[["f"]][[1]] | n_inf[["f"]][[2]] | n_inf[["f"]][[3]]) +
  plot_layout(guides = 'collect')

# n_inf[["t/f"]][[1]] / n_inf[["t/f"]][[2]] / n_inf[["t/f"]][[3]] +
#   plot_layout(guides = 'collect')

n_true = c()

n_true[["t"]][[1]] = plot_stuff(res$nama$true$c1, "t_pos", "", "nama") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
n_true[["t"]][[2]] = plot_stuff(res$nama$true$c2, "t_pos", "",  "nama") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
n_true[["t"]][[3]] = plot_stuff(res$nama$true$c3, "t_pos", "",  "nama") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))

n_true[["f"]][[1]] = plot_stuff(res$nama$true$c1, "f_pos", "", "nama") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
n_true[["f"]][[2]] = plot_stuff(res$nama$true$c2, "f_pos", "",  "nama") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
n_true[["f"]][[3]] = plot_stuff(res$nama$true$c3, "f_pos", "",  "nama") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))

# n_true[["t/f"]][[1]] = plot_stuff(res$nama$true$c1, c("t_pos","f_pos"), "", "nama") + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
# n_true[["t/f"]][[2]] = plot_stuff(res$nama$true$c2, c("t_pos","f_pos"), "",  "nama") + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))
# n_true[["t/f"]][[3]] = plot_stuff(res$nama$true$c3, c("t_pos","f_pos"), "",  "nama") + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))


(n_true[["t"]][[1]] | n_true[["t"]][[2]] | n_true[["t"]][[3]]) / 
  (n_true[["f"]][[1]] | n_true[["f"]][[2]] | n_true[["f"]][[3]]) +
  plot_layout(guides = 'collect')

# n_true[["t/f"]][[1]] / n_true[["t/f"]][[2]] / n_true[["t/f"]][[3]] +
#   plot_layout(guides = 'collect')
# 
# 
# n_true[["f"]][[1]] + scale_x_discrete(labels = c("laaa-nama","laaa-ep","laaa-gbr","apa-nama","apa-ep","apa-gbr","anc-nama","anc-ep","anc-gbr","allele-nama","allele-ep","allele-gbr"))



s_inf = c()

s_inf[["t"]][[1]] = plot_stuff(res$sac$inf$c1, "t_pos", "", "sac") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
s_inf[["t"]][[2]] = plot_stuff(res$sac$inf$c2, "t_pos", "",  "sac") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
s_inf[["t"]][[3]] = plot_stuff(res$sac$inf$c3, "t_pos", "",  "sac") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))

s_inf[["f"]][[1]] = plot_stuff(res$sac$inf$c1, "f_pos", "", "sac") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
s_inf[["f"]][[2]] = plot_stuff(res$sac$inf$c2, "f_pos", "",  "sac") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
s_inf[["f"]][[3]] = plot_stuff(res$sac$inf$c3, "f_pos", "",  "sac") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))

# s_inf[["t/f"]][[1]] = plot_stuff(res$sac$inf$c1, c("t_pos","f_pos"), "", "sac") + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
# s_inf[["t/f"]][[2]] = plot_stuff(res$sac$inf$c2, c("t_pos","f_pos"), "",  "sac") + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
# s_inf[["t/f"]][[3]] = plot_stuff(res$sac$inf$c3, c("t_pos","f_pos"), "",  "sac") + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))


(s_inf[["t"]][[1]] | s_inf[["t"]][[2]] | s_inf[["t"]][[3]]) / 
  (s_inf[["f"]][[1]] | s_inf[["f"]][[2]] | s_inf[["f"]][[3]]) +
  plot_layout(guides = 'collect')

# s_inf[["t/f"]][[1]] / s_inf[["t/f"]][[2]] / s_inf[["t/f"]][[3]] +
#   plot_layout(guides = 'collect')



s_true = c()

s_true[["t"]][[1]] = plot_stuff(res$sac$true$c1, "t_pos", "", "sac") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
s_true[["t"]][[2]] = plot_stuff(res$sac$true$c2, "t_pos", "",  "sac") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
s_true[["t"]][[3]] = plot_stuff(res$sac$true$c3, "t_pos", "",  "sac") + scale_y_continuous(breaks = 1:10, limits = c(0,10)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))

s_true[["f"]][[1]] = plot_stuff(res$sac$true$c1, "f_pos", "", "sac") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
s_true[["f"]][[2]] = plot_stuff(res$sac$true$c2, "f_pos", "",  "sac") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
s_true[["f"]][[3]] = plot_stuff(res$sac$true$c3, "f_pos", "",  "sac") + scale_y_continuous(limits = c(0,6000)) + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))

# s_true[["t/f"]][[1]] = plot_stuff(res$sac$true$c1, c("t_pos","f_pos"), "", "sac") + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
# s_true[["t/f"]][[2]] = plot_stuff(res$sac$true$c2, c("t_pos","f_pos"), "",  "sac") + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))
# s_true[["t/f"]][[3]] = plot_stuff(res$sac$true$c3, c("t_pos","f_pos"), "",  "sac") + scale_x_discrete(labels = c("laaa-nama","laaa-msl","laaa-gbr","laaa-eas","laaa-sas","apa-nama","apa-msl","apa-gbr","apa-eas","apa-sas","anc-nama","anc-msl","anc-gbr","anc-eas","anc-sas","allele-nama","allele-msl","allele-gbr","allele-eas","allele-sas"))


(s_true[["t"]][[1]]) / 
  (s_true[["f"]][[1]]) +
  plot_layout(guides = 'collect')

(s_true[["t"]][[2]]) / 
  (s_true[["f"]][[2]]) +
  plot_layout(guides = 'collect')

(s_true[["t"]][[3]]) / 
  (s_true[["f"]][[3]]) +
  plot_layout(guides = 'collect')

# s_true[["t/f"]][[1]] / s_true[["t/f"]][[2]] / s_true[["t/f"]][[3]] +
#   plot_layout(guides = 'collect')



# n_true[["t/f"]][[1]] + scale_y_continuous()



feature = "f_pos"
lim = 5000

trait = "laaa"
flat_c1 = map_df(res_c1[[trait]], ~as.data.frame(.x), .id="id")
flat_c2 = map_df(res_c2[[trait]], ~as.data.frame(.x), .id="id")
flat_c3 = map_df(res_c3[[trait]], ~as.data.frame(.x), .id="id")


plot_data = rbind(flat_c1[,c(1,12:15)], flat_c2[,c(1,12:15)], flat_c3[,c(1,12:15)])
plot_data[,"models"] = rep(c("null", "laaa", "global", "local", "apa"), 9)
plot_data[,"causal_set"] = c(rep(1, 15), rep(2, 15), rep(3, 15))

g_laaa = ggplot() + geom_col(aes(x = models, y = f_pos, fill = id, color = as.factor(causal_set)), data = plot_data, size = 1, position = "dodge") +
  scale_color_manual(values = c("black", "black", "black")) +
  labs(title = "LAAA", y = "True Positives") +
  scale_y_continuous(breaks = seq(1,lim,1), limits = c(0,lim))


trait = "apa"
flat_c1 = map_df(res_c1[[trait]], ~as.data.frame(.x), .id="id")
flat_c2 = map_df(res_c2[[trait]], ~as.data.frame(.x), .id="id")
flat_c3 = map_df(res_c3[[trait]], ~as.data.frame(.x), .id="id")


plot_data = rbind(flat_c1[,c(1,12:15)], flat_c2[,c(1,12:15)], flat_c3[,c(1,12:15)])
plot_data[,"models"] = rep(c("null", "laaa", "global", "local", "apa"), 9)
plot_data[,"causal_set"] = c(rep(1, 15), rep(2, 15), rep(3, 15))

g_apa = ggplot() + geom_col(aes(x = models, y = feature, fill = id, color = as.factor(causal_set)), data = plot_data, size = 1, position = "dodge") +
  scale_color_manual(values = c("black", "black", "black")) +
  labs(title = "APA", y = "True Positives") +
  scale_y_continuous(breaks = seq(1,10,1), limits = c(0,10))


trait = "anc"
flat_c1 = map_df(res_c1[[trait]], ~as.data.frame(.x), .id="id")
flat_c2 = map_df(res_c2[[trait]], ~as.data.frame(.x), .id="id")
flat_c3 = map_df(res_c3[[trait]], ~as.data.frame(.x), .id="id")


plot_data = rbind(flat_c1[,c(1,12:15)], flat_c2[,c(1,12:15)], flat_c3[,c(1,12:15)])
plot_data[,"models"] = rep(c("null", "laaa", "global", "local", "apa"), 9)
plot_data[,"causal_set"] = c(rep(1, 15), rep(2, 15), rep(3, 15))

g_anc = ggplot() + geom_col(aes(x = models, y = feature, fill = id, color = as.factor(causal_set)), data = plot_data, size = 1, position = "dodge") +
  scale_color_manual(values = c("black", "black", "black")) +
  labs(title = "LAO", y = "True Positives") +
  scale_y_continuous(breaks = seq(1,10,1), limits = c(0,10))


trait = "allele"
flat_c1 = map_df(res_c1[[trait]], ~as.data.frame(.x), .id="id")
flat_c2 = map_df(res_c2[[trait]], ~as.data.frame(.x), .id="id")
flat_c3 = map_df(res_c3[[trait]], ~as.data.frame(.x), .id="id")


plot_data = rbind(flat_c1[,c(1,12:15)], flat_c2[,c(1,12:15)], flat_c3[,c(1,12:15)])
plot_data[,"models"] = rep(c("null", "laaa", "global", "local", "apa"), 9)
plot_data[,"causal_set"] = c(rep(1, 15), rep(2, 15), rep(3, 15))

g_allele = ggplot() + geom_col(aes(x = models, y = feature, fill = id, color = as.factor(causal_set)), data = plot_data, size = 1, position = "dodge") +
  scale_color_manual(values = c("black", "black", "black")) +
  labs(title = "Allele", y = "True Positives") +
  scale_y_continuous(breaks = seq(1,10,1), limits = c(0,10))


(g_laaa | g_apa) / (g_anc | g_allele) +
  plot_layout(guides = 'collect')



ggplot() + geom_point(aes(x = models, y = t_pos, color = id, shape = as.factor(causal_set)), data = plot_data, size = 5, position = position_dodge(width = .5))


ggplot() + geom_col(aes(x = row.names(flat_c3_laaa), y = t_pos, fill = id), data = flat_c3_laaa)


laaa
apa
anc
allele


