library(magrittr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(rstatix)
library(reshape2)

###Combined maps####
read_accs = function(path, suffix, sep){
  accs = system(paste("ls ", path, suffix, sep = ""), intern = T)
  lapply(accs, function(e){
    readLines(accs[14]) %>%
      strsplit(., split = sep) %>%
      .[-1] %>%
      unlist(., recursive = F) %>%
      as.numeric()
  }) %>%
    do.call("cbind", .) -> 
    tabs
  
  sapply(accs, function(e){
    gsub(paste(path, "AA_sim_", sep = ""),"",e) %>%
      gsub(".acc.txt", "", .) %>%
      gsub(".merged.tsv.acc2", "", .) %>%
      gsub("_19_100", "", .) %>%
      gsub("_19_s1111", "", .)%>%
      gsub("_100", "", .)
  }) ->
    n
  
  colnames(tabs) = n
  
  return(tabs)
}

accs_n = system(paste("ls ", "/home/gerald/Documents/PhD/papers/paper4/accuracy/nama_ancs_1000_5000_*.acc.tsv", sep = ""), intern = T)

accs_s = system(paste("ls ", "/home/gerald/Documents/PhD/papers/paper4/accuracy/sac_ancs_1000_5000_*.acc.tsv", sep = ""), intern = T)

lapply(accs_n, function(chr){
  readLines(chr) %>%
    strsplit(., split = "\t") %>%
    .[-1] %>%
    lapply(., function(e){
      e[-c(1:3)] %>%
        as.numeric(e) ->
        res
      
      res[-1]/res[1]
    }) %>%
    do.call("rbind", .) %>%
    colMeans()
}) ->
  result_n

lapply(accs_s, function(chr){
  readLines(chr) %>%
    strsplit(., split = "\t") %>%
    .[-1] %>%
    lapply(., function(e){
      e[-c(1:3)] %>%
        as.numeric(e) ->
        res
      
      res[-1]/res[1]
    }) %>%
    do.call("rbind", .) %>%
    colMeans()
}) ->
  result_s


sapply(result_n, mean) %>%
  mean
sapply(result_s, mean) %>%
  mean


process_accs = function(raw_accs){
  melt(raw_accs)[,2:3] ->
    melted
  
  melted[,2] = as.numeric(melted[,2])
  melted = data.frame(melted)
  colnames(melted) = c("group", "acc")
  
  return(melted)
}

split_accs = function(raw_accs, haps, ancs){
  nacc = nrow(raw_accs)/haps
  
  lapply(1:nacc, function(e){
    raw_accs[(e*haps-haps+1):(e*haps),] %>%
      melt %>%
      .[-1] ->
      melted
    
    df = data.frame("group" = melted[,1], "acc" = melted[,2], "anc" = ancs[e])
    
    df
  })
}


lapply(accs_n, function(chr){
  readLines(chr) %>%
    strsplit(., split = "\t") %>%
    .[-1] %>%
    lapply(., function(e){
      e[-c(1:3)] %>%
        as.numeric(e) ->
        res
      
      res[-1]/res[1]
    })
}) ->
  raw_n














adm_2_raw = read_accs("/home/gerald/Documents/PhD/papers/paper3/lai_acc/2_way/", "*.acc.txt", " ")
adm_2_raw = adm_2_raw[,c(9,11,7,6,3,1,2)]
colnames(adm_2_raw) = c("LAI (100)", "LAI (2000)", "IBD (100)", "IBD (2000)", "LD (100)", "Recent", "Ancestral")
adm_3_raw = read_accs("/home/gerald/Documents/PhD/papers/paper3/lai_acc/3_way/", "*.acc.txt", " ")
adm_3_raw = adm_3_raw[,c(9,11,7,6,3,1,2)]
colnames(adm_3_raw) = c("LAI (100)", "LAI (2000)", "IBD (100)", "IBD (2000)", "LD (100)", "Recent", "Ancestral")
adm_2 = process_accs(adm_2_raw)
adm_2 = data.frame(adm_2, "pop" = "Two-way")
adm_3 = process_accs(adm_3_raw)
adm_3 = data.frame(adm_3, "pop" = "Three-way")

adm = rbind(adm_2, adm_3)


acc_2_raw = read_accs("/home/gerald/Documents/PhD/papers/paper3/lai_acc/2_way/", "*.merged.tsv.acc2", "\t")
acc_2_raw = acc_2_raw[,c(9,11,7,6,3,1,2)]
colnames(acc_2_raw) = c("LAI (100)", "LAI (2000)", "IBD (100)", "IBD (2000)", "LD (100)", "Recent", "Ancestral")
acc_3_raw = read_accs("/home/gerald/Documents/PhD/papers/paper3/lai_acc/3_way/", "*.merged.tsv.acc2", "\t")
acc_3_raw = acc_3_raw[,c(9,11,7,6,3,1,2)]
colnames(acc_3_raw) = c("LAI (100)", "LAI (2000)", "IBD (100)", "IBD (2000)", "LD (100)", "Recent", "Ancestral")

acc_2 = split_accs(acc_2_raw, 200, c("African", "European"))
acc_3 = split_accs(acc_3_raw, 200, c("African", "East Asian", "European"))

do.call("rbind", acc_2) ->
  anc_2

do.call("rbind", acc_3) ->
  anc_3

res.kruskal <- dat %>% kruskal_test(acc ~ group)
pwc <- dat %>% 
  dunn_test(acc ~ group, p.adjust.method = "bonferroni") 

nsig = sum(pwc$p.adj.signif != "ns")

get_test_label(kruskal_test(data = adm_2, formula = acc ~ group))

w2_w3 =  ggboxplot(adm, x = "group", y = "acc", color = "pop") +
  scale_y_reverse(breaks = rev(seq(1, 0.79, -0.1))) +
  labs(x = "Recombination map used for LAI", 
       y = "Accuracy", 
       caption = get_test_label(kruskal_test(data = adm_2, formula = acc ~ group)),
       color = "Admixture") +
  theme(axis.text = element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.caption = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(0.5,0,0,0,"cm")),
        axis.title.y = element_text(margin = margin(0,0.5,0,0, "cm")),
        plot.margin = margin(1,1,1,1,"cm"),
        legend.position = "right", 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16),
        legend.key.height = unit(3, "lines"),
        legend.key.width = unit(1.5, "lines")) + 
  coord_flip(ylim = c(1, 0.79))  +
  color_palette(rev(viridis(5)[c(1,4)]))

w2_all =  ggboxplot(anc_2, x = "group", y = "acc", color = "anc") +
  scale_y_reverse(breaks = rev(seq(1, 0, -0.1))) +
  labs(x = NULL, 
       #x = "Recombination map used", 
       y = NULL,
       color = "Ancestry") +
  theme(axis.text = element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.caption = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(0.5,0,0,0,"cm")),
        axis.title.y = element_text(margin = margin(0,0.5,0,0, "cm")),
        plot.margin = margin(0,0,0,0.5,"cm")) + 
  coord_flip(ylim = c(1, 0))  +
  color_palette(rev(inferno(5)[c(1,4)])) 

w3_all =  ggboxplot(anc_3, x = "group", y = "acc", color = "anc", size = 1.05) +
  scale_y_reverse(breaks = rev(seq(1, 0, -0.1))) +
  labs(#x = NULL, 
       x = "Recombination map used", 
       y = "Accuracy", 
       color = "Ancestry") +
  theme(axis.text = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold"),
        plot.caption = element_text(size = 14), 
        axis.title = element_text(size = 20, face = "bold"),
        axis.title.x = element_text(margin = margin(1,0,0,0,"cm")),
        axis.title.y = element_text(margin = margin(0,1,0,0, "cm"))) + 
  coord_flip(ylim = c(1, 0)) + 
  color_palette(rev(inferno(5)[c(1,3,4)])) + 
  theme(legend.position = "right", 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 20, face = "bold"),
        legend.key.height = unit(3, "lines"),
        legend.key.width = unit(1.5, "lines"),
        panel.grid.major.x = element_line(colour = rgb(0,0,0,0.1)))


png(filename = "/home/gerald/Documents/PhD/papers/paper3/writeup/acc_per_anc_3way_sw.png", units = "cm", width = 40, height = 30, res = 600)

w3_all +
  plot_layout(guides = "collect") & 
  theme(legend.position = ("right"), 
        plot.margin = margin(1,1,1,1,"cm"))

dev.off()


find_max_sig = function(dat, space){
  sapply(dat, function(e){
    pwc <- e %>% 
      dunn_test(acc ~ group, p.adjust.method = "bonferroni") 
    
    nsig = sum(pwc$p.adj.signif != "ns")
  }) %>%
    max * space - space + (1 + space)
}

plot_box <-  function(dat, title = NULL, remove_axis = Fm, ymin = 0.79, xticks = NULL, sig_space = 0.004, ymax = 1, vjust = 2, br_just = 0.0035, size = 3.88, combine = F, showmean = T, showsig = T){
  res.kruskal <- dat %>% kruskal_test(acc ~ group)
  pwc <- dat %>% 
    dunn_test(acc ~ group, p.adjust.method = "bonferroni") 
  
  nsig = sum(pwc$p.adj.signif != "ns")
  
  if(nsig < 1){
    br = -1
  } else{
    br = -seq((1+sig_space), 1+sig_space+br_just+(nsig*sig_space)-sig_space, sig_space) 
    l = br[1]
    sl = br[2]
    br[1] = sl
    br[2] = l
    }
  
  pwc <-  add_xy_position(pwc, x = "group")
  ggboxplot(dat, x = "group", y = "acc") +
    scale_y_reverse(breaks = rev(seq(1, ymin, -0.1))) +
    color_palette("black") + 
    #ylim(ymax, ymin) +
    stat_summary(fun=mean, geom="point", size=3, color="red", na.rm = T) +
    stat_pvalue_manual(pwc, 
                       hide.ns = TRUE,
                       y.position = if(showsig){br}else{1000000},
                       tip.length = -sig_space + br_just, 
                       vjust = vjust,
                       coord.flip = T, 
                       size = size) + 
    labs(title = title,
         caption = get_test_label(res.kruskal), 
         x = "Recombination map used", 
         #x = NULL, 
         y = if(!remove_axis){"Accuracy"}) +
    #legend(legend = c("LAI - LAI 100 individuals", "LAI 2k - LAI 2000 individuals", "IBD - IBD 100 individuals", "IBD 2k - IBD 2000 individuals"), x = "bottom", ) + 
    theme(axis.text = element_text(size = 14),
          axis.text.x = xticks,
          plot.title = element_text(size = 16),
          plot.caption = element_text(size = 14), 
          axis.title = element_text(size = 16, face = "bold"),
          axis.title.x = element_text(margin = margin(1,0,0,0,"cm")),
          axis.title.y = element_text(margin = margin(0,1,0,0, "cm")),
          plot.margin = margin(0,0,0,1,"cm"),
          #panel.background = element_rect(fill = rgb(0,0,0,0.25)),
          panel.grid.major.x = element_line(colour = rgb(0,0,0,0.1))) + 
    coord_flip(ylim = c(ymax, ymin)) 
    # coord_cartesian(ylim = c(ymax, ymin))
  
}

find_max_union = function(mat, arr){
  if(length(arr) == 0){
    find_max_union(mat[,-ncol(mat)], intersect(mat[,1], mat[,ncol(mat)]))
  } else if(length(dim(mat)) > 1){
    find_max_union(mat[,-ncol(mat)], intersect(mat[,ncol(mat)], arr))
  } else{
    return(intersect(mat, arr))
  }
}

# adm_2_raw = read_accs("/home/gerald/Documents/PhD/papers/paper3/lai_acc/2_way/", "*.acc.txt", " ")
# adm_2_raw = adm_2_raw[,c(9,11,7,6,3,1,2)]
# colnames(adm_2_raw) = c("LAI (100)", "LAI (2000)", "IBD (100)", "IBD (2000)", "FastEPRR", "Recent", "Ancestral")
# adm_3_raw = read_accs("/home/gerald/Documents/PhD/papers/paper3/lai_acc/3_way/", "*.acc.txt", " ")
# adm_3_raw = adm_3_raw[,c(9,11,7,6,3,1,2)]
# colnames(adm_3_raw) = c("LAI (100)", "LAI (2000)", "IBD (100)", "IBD (2000)", "FastEPRR", "Recent", "Ancestral")
# adm_2 = process_accs(adm_2_raw)
# adm_3 = process_accs(adm_3_raw)


#Compare groups
w2 = plot_box(dat = adm_2,
              vjust = 2.5,
              br_just = 0.004,
              size = 5,
              # title = "Two-way admixed", 
              remove_axis = T, ymin = 0.79, sig_space = 0.006,
              #xticks = element_text(angle = 45, vjust = 0.5), 
              ymax = find_max_sig(list(adm_2, adm_3), space = 0.006))
w3 = plot_box(dat = adm_3,
              vjust = 2.5,
              br_just = 0.004,
              size = 5, 
              # title = "Three-way admixed", 
              remove_axis = F, ymin = 0.79, sig_space = 0.006,
              #xticks = element_text(angle = 45, vjust = 0.5), 
              ymax = find_max_sig(list(adm_2, adm_3), space = 0.006))

png(filename = "/home/gerald/Documents/PhD/papers/paper3/writeup/acc_2_3way_sw.png", units = "cm", width = 24, height = 30, res = 600)
(w2  / w3) +  plot_annotation(tag_levels = "A",
                              theme = theme(plot.title = element_text(size = 24, 
                                                                      margin = margin(0,0,1,0, "cm")),
                                            plot.margin = margin(1,1,1,0,"cm"))) + 
  plot_layout(guides = "collect") & 
  theme(plot.tag = element_text(size = 20, face = "bold"))
dev.off()

acc_2_raw = read_accs("/home/gerald/Documents/PhD/papers/paper3/lai_acc/2_way/", "*.merged.tsv.acc2", "\t")
acc_2_raw = acc_2_raw[,c(9,11,7,6,3,1,2)]
colnames(acc_2_raw) = c("LAI (100)", "LAI (2000)", "IBD (100)", "IBD (2000)", "FastEPRR", "Recent", "Ancestral")
acc_3_raw = read_accs("/home/gerald/Documents/PhD/papers/paper3/lai_acc/3_way/", "*.merged.tsv.acc2", "\t")
acc_3_raw = acc_3_raw[,c(9,11,7,6,3,1,2)]
colnames(acc_3_raw) = c("LAI (100)", "LAI (2000)", "IBD (100)", "IBD (2000)", "FastEPRR", "Recent", "Ancestral")

acc_2 = split_accs(acc_2_raw, 200, c("African", "European", "East Asian"))
acc_3 = split_accs(acc_3_raw, 200, c("African", "East Asian", "European"))

# lapply(acc_2, function(e){
#   remove = c()
#   # remove = which(e[,2] == 0)
#   if(length(remove) > 0){
#     e[-remove,] 
#   } else{
#     e
#   }
# }) ->
#   acc_2
# 
# lapply(acc_3, function(e){
#   remove = c()
#   #remove = which(e[,2] == 0)
#   if(length(remove) > 0){
#     e[-remove,] 
#   } else{
#     e
#   }
# }) ->
#   acc_3

# ancs = c("African", "East Asian", "European")

do.call("rbind", acc_2) ->
  anc_2

do.call("rbind", acc_3) ->
  anc_3


w2_afr = plot_box(acc_2[[1]], 
                  # title = "African", 
                  remove_axis = T, ymin = 0, 
                  #xticks = element_text(angle = 45, vjust = 0.5),
                  sig_space = 0.01,
                  vjust = 2.1,
                  br_just = 0.0075,
                  size = 2.88,
                  ymax = find_max_sig(acc_2, space = 0.01)) +
  stat_summary(fun=mean, geom="point", size=3, color="red", na.rm = T)
w2_eur = plot_box(acc_2[[2]], 
                  # title = "European", 
                  remove_axis = F, ymin = 0, 
                  #xticks = element_text(angle = 45, vjust = 0.5),
                  sig_space = 0.01,
                  vjust = 2.1,
                  br_just = 0.0075,
                  size = 2.88,
                  ymax = find_max_sig(acc_2, space = 0.01)) +
  stat_summary(fun=mean, geom="point", size=3, color="red", na.rm = T)

w3_afr = plot_box(acc_3[[1]], 
                  # title = "African",
                  remove_axis = T, ymin = 0, 
                  #xticks = element_text(angle = 45, vjust = 0.5),
                  sig_space = 0.012,
                  vjust = 2.2,
                  br_just = 0.008,
                  size = 2.88,
                  ymax = find_max_sig(acc_3, space = 0.012)) +
  stat_summary(fun=mean, geom="point", size=3, color="red", na.rm = T)
w3_eur = plot_box(acc_3[[2]], 
                  # title = "European", 
                  remove_axis = T, ymin = 0, 
                  #xticks = element_text(angle = 45, vjust = 0.5),
                  sig_space = 0.012,
                  vjust = 2.2,
                  br_just = 0.008,
                  size = 2.88,
                  ymax = find_max_sig(acc_3, space = 0.012)) +
  stat_summary(fun=mean, geom="point", size=3, color="red", na.rm = T)
w3_eas = plot_box(acc_3[[3]], 
                  # title = "East Asian", 
                  remove_axis = F, ymin = 0, 
                  #xticks = element_text(angle = 45, vjust = 0.5),
                  sig_space = 0.012,
                  br_just = 0.008,
                  size = 2.88,
                  vjust = 2.2,
                  ymax = find_max_sig(acc_3, space = 0.012)) +
  stat_summary(fun=mean, geom="point", size=3, color="red", na.rm = T)

w2_all =  ggboxplot(anc_2, x = "group", y = "acc", color = "anc") +
  scale_y_reverse(breaks = rev(seq(1, 0, -0.1))) +
  labs(x = NULL, 
       #x = "Recombination map used", 
       y = NULL,
       color = "Ancestry") +
  theme(axis.text = element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.caption = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(0.5,0,0,0,"cm")),
        axis.title.y = element_text(margin = margin(0,0.5,0,0, "cm")),
        plot.margin = margin(0,0,0,0.5,"cm")) + 
  coord_flip(ylim = c(1, 0))  +
  color_palette(rev(inferno(5)[c(1,4)])) 
# +
#   theme(legend.position = "none", 
#         legend.key.size = unit(0, "mm"), 
#         legend.text = element_text(size = 0), 
#         legend.title = element_text(size = 0))

w3_all =  ggboxplot(anc_3, x = "group", y = "acc", color = "anc") +
  scale_y_reverse(breaks = rev(seq(1, 0, -0.1))) +
  labs(x = NULL, 
       #x = "Recombination map used", 
       y = "Accuracy", 
       color = "Ancestry") +
  theme(axis.text = element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.caption = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(0.5,0,0,0,"cm")),
        axis.title.y = element_text(margin = margin(0,0.5,0,0, "cm")),
        plot.margin = margin(0,0,0,0.5,"cm")) + 
  coord_flip(ylim = c(1, 0)) + 
  color_palette(rev(inferno(5)[c(1,3,4)])) + 
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16),
        legend.key.height = unit(3, "lines"),
        legend.key.width = unit(1.5, "lines"))

two = (w2_afr / w2_eur)
three = (w3_afr / w3_eur / w3_eas)
both = two | three



# png(filename = "/home/gerald/Documents/PhD/papers/paper3/writeup/acc_per_anc_2way_sw.png", units = "cm", width = 24, height = 40, res = 300)
# two + plot_annotation(
#   #title = "Comparison of LAI accuracy per ancestry", 
#   theme = theme(plot.title = element_text(size = 24, margin = margin(0,0,1,0, "cm")),
#                 plot.margin = margin(1,1,1,0,"cm"))
# )
# dev.off()
# 
# png(filename = "/home/gerald/Documents/PhD/papers/paper3/writeup/acc_per_anc_3way_sw.png", units = "cm", width = 24, height = 40, res = 300)
# three + plot_annotation(
#   #title = "Comparison of LAI accuracy per ancestry", 
#   theme = theme(plot.title = element_text(size = 24, margin = margin(0,0,1,0, "cm")),
#                 plot.margin = margin(1,1,1,0,"cm"))
# )
# dev.off()
# 
# png(filename = "/home/gerald/Documents/PhD/papers/paper3/writeup/acc_per_anc_2_3way_sw.png", units = "cm", width = 40, height = 30, res = 300)
# (two | three) +
#   plot_annotation(theme = theme(plot.margin = margin(0.5,0.5,0.5,0,"cm")), tag_levels = "A")
# dev.off()

# png(filename = "/home/gerald/Documents/PhD/papers/paper3/writeup/acc_per_anc_2_3way_sw2.png", units = "cm", width = 40, height = 30, res = 300)
# 
# ((w2 | w2_all & theme(legend.position = "none", 
#                       plot.margin = margin(0.5,0.5,0.5,0.5,"cm"))) / 
#     (w3 | w3_all & theme(legend.position = ("right"), 
#                          plot.margin = margin(0.5,0.5,0.5,0.5,"cm")))) + 
#   plot_annotation(theme = theme(plot.margin = margin(0.5,0.5,0.5,0.5,"cm")), tag_levels = "A") + 
#   plot_layout(guides = "collect") &
#   theme(plot.tag = element_text(size = 20, face = "bold"))
# 
# dev.off()




read_accs = function(path, suffix, sep){
  path = "/home/gerald/Documents/PhD/papers/paper3/lai_acc/2_way/"
  suffix = "*.merged.tsv.acc2"
  accs = system(paste("ls ", path, suffix, sep = ""), intern = T)
  lapply(accs, function(e){
    readLines(e) %>%
      strsplit(., split = "\t") %>%
      unlist  %>%
      as.numeric()
  }) %>%
    do.call("cbind", .) -> 
    tabs
  
  sapply(accs, function(e){
    gsub(paste(path, "AA_sim_", sep = ""),"",e) %>%
      gsub(".acc.txt", "", .) %>%
      gsub("merged.tsv.acc", "", .) %>%
      gsub("_19_100", "", .) %>%
      gsub("_19_s1111", "", .)%>%
      gsub("_100", "", .)
  }) ->
    n
  
  colnames(tabs) = n
  
  return(tabs)
}








acc_2[[2]] %>% na.omit() %>% .[order(.[,2]),] %>% .[1:20,]
acc_3[[2]] %>% na.omit() %>% .[order(.[,2]),] %>% .[1:20,]
acc_3[[3]] %>% na.omit() %>% .[order(.[,2]),] %>% .[1:20,]
















#Read truth files
simfile2 = "/home/gerald/Documents/PhD/papers/paper3/lai_acc/2_way/AA_sim_Afr_Eur_19_100.hanc"
sim2 = readLines(simfile2)
sim2 = do.call("rbind", lapply(sim2, function(e) as.numeric(unlist(strsplit(e, split = "")))))
sim2 = t(sim2)
simfile3 = "/home/gerald/Documents/PhD/papers/paper3/lai_acc/3_way/AA_sim_Afr_Eur_EA_19_100.hanc"
sim3 = readLines(simfile3)
sim3 = do.call("rbind", lapply(sim3, function(e) as.numeric(unlist(strsplit(e, split = "")))))
sim3 = t(sim3)

get_sim_summary = function(sim) {
  apply(sim, 2, function(e){
      x = 1
      l = c()
      init = 1
      for(i in 2:length(e)){
        if(e[i] != e[i-1]){
          l[[x]] = c(i-init, e[i-1])
          x = x+1
          init = i
        }
        if(i == length(e)){
          l[[x]] = c(i-init, e[i])
        }
      } 
      frame = do.call("rbind", l)
      data.frame("len" = frame[,1], "anc" = frame[,2])
  }) ->
    switches
}

v = get_sim_summary(sim3)

anc_sum = c()

sapply(v, function(e){
  unique(e[,2]) %>% 
    sort() ->
    ancs
  sapply(ancs, function(f){
    filter(e, anc == f) %>%
      unlist(.[,1]) %>%
      mean() %>%
      c(f,.)
  })
}) ->
  anc_sum

sapply(1:length(anc_sum), function(e){
  lapply(seq(0,1200,200), function(f){
    sapply(1:3, function(g){
      acc_3[[g]][e+f,2]
    }) %>%
      na.omit %>%
      rbind(anc_sum[[e]], .) %>%
      t
  }) %>%
    do.call("rbind", .)
}) ->
  w

do.call("rbind", w) %>%
  data.frame()->
  x

x[x[,1] == 0,1] = "African"
x[x[,1] == 1,1] = "East Asian"
x[x[,1] == 2,1] = "European"

library(viridis)
col1 <- colorRampPalette(head(viridis(10)))

png(filename = "/home/gerald/Documents/PhD/papers/paper3/writeup/acc_vs_len_per_anc_3way.png", units = "cm", width = 40, height = 24, res = 300)
ggplot() + 
  geom_point(mapping = aes(x = x[,2], y = x[,3], color = as.factor(x[,1])), size = 3, data = x) +
  color_palette(rev(inferno(5)[c(1,3,4)])) + 
  labs(x = "Segment length (bp)",
       y = "Accuracy",
       color = "Ancestry") + 
  scale_x_continuous(breaks = c(100000,200000,300000,400000),
                     labels = c("100 000","200 000","300 000","400 000")) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.title.x = element_text(margin = margin(0.5,0,0,0,"cm")),
        axis.title.y = element_text(margin = margin(0,0.5,0,0, "cm")),
        plot.margin = margin(1,1,1,1,"cm"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14), 
        panel.background =  element_rect(fill = "white", 
                                         colour = "gray50"), 
        panel.grid = element_line(colour = rgb(0,0,0,0.1)),
        # panel.grid.major = element_line(colour = rgb(0,0,0,0.1))
        )+
  guides(color = guide_legend(override.aes = list(size = 7, shape = "square") ) )
dev.off()


v[[164]]
v[[181]]
v[[183]]

# apply(acc_3_raw, 2, function(e){
#   data.frame("idx" = seq(1, length(e), 1), "val" = e) %>%
#     .[order(.[,2]),] %>%
#     .[1:100,]
# }) ->
#   mins3



acc_3[[2]] %>%
  data.frame("id" = rep(1:200,7), .) %>%
  dcast(., id ~ group, fun.aggregate = sum, value.var = "acc") %>%
  apply(., 2, function(e){
    data.frame("idx" = seq(1, length(e), 1), "val" = e) %>%
      .[order(.[,2]),] %>%
      .[.[,2] != 0,] %>%
     .[1:5,]
  }) ->
    mins3


sapply(mins3, function(e){e[,1]}) %>% 
  cbind() %>%
  find_max_union(., c()) ->
  lowest_shared3

# sapply(mins3, function(e){
#   e[e[,1] %in% lowest_shared3,2]
# }) %>%
#   process_accs %>%
#   plot_box
#   
# apply(adm_3_raw, 2, function(e){
#   quantile(e,0.1)
# })
# 
# 
# sapply(2:ncol(mins), function(e){
#   union(mins[,e-1],mins[,e])
# })
# 
# sapply(1:ncol(adm_2_raw), function(e){
#   c(mins) %>%
#     unique ->
#     umins
#   
#     which(umins %in% mins[,e])
# }) %>%
#   c %>%
#   dupli
  
#sim_2_low = sim2[,lowest_shared2]
#sim_3_low = sim3[,lowest_shared3]
sim_3_low = sim3[,c(164,167,181)]
#sim_2_high = sim2[,-lowest_shared2]
sim_3_high = sim3[,-lowest_shared3]


apply(sim_3_low[,1:ncol(sim_3_low)], 2, function(e){
  x = 1
  l = c()
  init = 1
  for(i in 2:length(e)){
    if(e[i] != e[i-1]){
      l[[x]] = i-init
      x = x+1
      init = i
    }
  }
  l
}) ->
  switches_h2

sapply(switches_h2, function(e){
  length(unlist(e))
}) ->
  num_switches_h2

res.kruskal2 <- adm_2 %>% kruskal_test(acc ~ group)
res.kruskal3 <- adm_3 %>% kruskal_test(acc ~ group)
pwc2 <- adm_2 %>% 
  dunn_test(acc ~ group, p.adjust.method = "bonferroni") 
pwc3 <- adm_3 %>% 
  dunn_test(acc ~ group, p.adjust.method = "bonferroni") 
pwc2
pwc3


pwc2[pwc2$p.adj.signif != "ns",]
pwc3[pwc3$p.adj.signif != "ns",]

pwc2[pwc2$p <= 0.05,]
pwc3[pwc3$p <= 0.05,]




msp = read.table("/home/gerald/Documents/PhD/papers/paper3/lai_acc/2_way/AA_sim_aa_19_100.merged.tsv", sep = "\t", header = T)

sapply(1:nrow(sim), function(e){
  sum(sim[e,] == msp[e,])
}) ->
  lat_sum


lat_sum_dat = cbind(1:nrow(sim), lat_sum)
colnames(lat_sum_dat) = c("pos", "sum")

lat_sum_dat = data.frame(lat_sum_dat)

ggplot(data = lat_sum_dat, mapping = aes(x = pos, y = sum))  + geom_line()





###LAI Accuracy
# read_accs = function(path){
#   accs = system(paste("ls ", path,"*.acc.tsv", sep = ""), intern = T)
#   lapply(accs, function(e){
#     read.table(e, header = T)
#   }) -> 
#     tabs
#   
#   sapply(accs, function(e){
#     gsub(paste(path, "AA_sim_", sep = ""),"",e) %>%
#       gsub(".acc.tsv", "", .)
#   }) ->
#     n
#   
#   names(tabs) = n
#   
#   return(tabs)
# }

# process_accs = function(raw_accs, filt){
#   if(length(filt) > 0){
#     lapply(raw_accs, function(f){
#       lapply(5:ncol(f), function(e){
#         mean(f[,e] / f[,4])
#       }) %>%
#         unlist() %>%
#         .[-filt]
#     }) ->
#       accs
#   } else{
#     lapply(raw_accs, function(f){
#       lapply(5:ncol(f), function(e){
#         mean(f[,e] / f[,4])
#       }) %>%
#         unlist()
#     }) ->
#       accs
#     }
#  
#   
#   lapply(1:length(accs), function(e){
#     cbind(names(raw_accs)[e],accs[[e]])
#   }) %>%
#     do.call("rbind", .) %>%
#     data.frame ->
#     melted
#   
#   melted[,2] = as.numeric(melted[,2])
#   melted = data.frame(melted)
#   colnames(melted) = c("group", "acc")
#   
#   return(melted)
# }



print_means = function(acc){
  unique(acc[,1]) %>%
    sapply(., function(e){
      mean(acc[acc[,1]==e, 2])
    })
}

# simfile = "/home/gerald/Documents/PhD/papers/paper3/lai_acc/2_way/AA_sim_Afr_Eur_19_100.hanc"
# sim = readLines(simfile)
# sim = do.call("rbind", lapply(sim, function(e) as.numeric(unlist(strsplit(e, split = "")))))
# sim = t(sim)
# 
# apply(sim[,1:ncol(sim)], 2, function(e){
#   x = 1
#   l = c()
#   init = 1
#   for(i in 2:length(e)){
#     if(e[i] != e[i-1]){
#       l[[x]] = i-init
#       x = x+1
#       init = i
#     }
#   }
#   l
# }) ->
#   switches
# 
# sapply(switches, function(e){
#   length(unlist(e))
# }) ->
#   num_switches
# 
# filt = which(num_switches == 0)

filt = c()

# ra1 = read_accs("/home/gerald/Documents/PhD/papers/paper3/lai_acc/2_way/1111/")
# ra2 = read_accs("/home/gerald/Documents/PhD/papers/paper3/lai_acc/2_way/22222/")
# ra3 = read_accs("/home/gerald/Documents/PhD/papers/paper3/lai_acc/2_way/333333/")
# a1 = process_accs(ra1, filt)
# a2 = process_accs(ra2, filt)
# a3 = process_accs(ra3, filt)
# #Means
# run_means = rbind(print_means(a1), print_means(a2), print_means(a3))
# colMeans(run_means)



library(tidyverse)

res.kruskal <- melted %>% kruskal_test(acc ~ group)
res.kruskal
pwc <- melted %>% 
  dunn_test(acc ~ group, p.adjust.method = "bonferroni") 
pwc
pwc <-  add_xy_position(pwc, x = "group")
ggboxplot(melted, x = "group", y = "acc") +
  stat_pvalue_manual(pwc, hide.ns = TRUE, y.position = seq(1.01, 1.04, 0.01)) +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )


## Box plots of lai accuracy
ggplot(melted, aes(x = melted[,1], y = as.numeric(melted[,2]))) + 
  geom_boxplot() +
  ylab("LAI Accuracy") + 
  xlab("RM used in LAI") + 
  labs(title = "LAI of 200 simulated AA haplotypes using different RM's") + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        title=element_text(size=24))

ggplot(melted2, aes(x = melted2[,1], y = as.numeric(melted2[,2]))) + 
  geom_boxplot() +
  ylab("LAI Accuracy") + 
  xlab("RM used in LAI") + 
  labs(title = "LAI of 200 simulated AA haplotypes using different RM's") + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18,face="bold"),
        title=element_text(size=24))

shapiro.test(as.numeric(melted[melted[,1]=="aapub",2]))

kruskal.test (melted[,1], as.numeric(melted[,2]))
# p-value = 0.0002869
melted[melted[,1]=="fast" | melted[,1]=="rfmix100" | melted[,1]=="rfmix2000",] %>%
  kruskal.test(g = .[,1], x = as.numeric(.[,2]))
#p-value = 0.02696
wilcox.test(as.numeric(melted[melted[,1]=="fast",2]), 
            as.numeric(melted[melted[,1]=="rfmix100",2]))
#p-value = 0.03124
wilcox.test(as.numeric(melted[melted[,1]=="fast",2]), 
            as.numeric(melted[melted[,1]=="rfmix2000",2]))
# p-value = 0.8571
wilcox.test(as.numeric(melted[melted[,1]=="rfmix100",2]), 
            as.numeric(melted[melted[,1]=="rfmix2000",2]))
#p-value = 0.01339

melted[melted[,1]=="aapub" | melted[,1]=="combpub" | melted[,1]=="rfmix100",] %>%
  kruskal.test(g = .[,1], x = as.numeric(.[,2]))
#p-value = 0.5507
wilcox.test(as.numeric(melted[melted[,1]=="aapub",2]), 
            as.numeric(melted[melted[,1]=="combpub",2]))
#p-value = 0.676



#Correlate number of switch points and accuracy
apply(sim[,1:ncol(sim)], 2, function(e){
  x = 1
  l = c()
  init = 1
  for(i in 2:length(e)){
    if(e[i] != e[i-1]){
      l[[x]] = i-init
      x = x+1
      init = i
    }
  }
  l
}) ->
  switches

sapply(switches, function(e){
  length(unlist(e))
}) ->
  num_switches

sapply(switches, function(e){
  mean(unlist(e))
}) ->
  mean_anc_seg_len
mean_anc_seg_len[is.na(mean_anc_seg_len)] = 1

acc_s_len = data.frame("acc" = accs[[3]], "s" = num_switches, "len" = mean_anc_seg_len)

ggplot(data = acc_s_len, mapping = aes(x = acc, y = s)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)+
  ylab("# Switchpoints") + 
  xlab("LAI Accuracy")

ggplot(data = acc_s_len, mapping = aes(x = acc, y = len)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)+
  ylab("Anc block size") + 
  xlab("LAI Accuracy")

cor(acc_s_len)


apply(msp_fast[,7:ncol(msp_fast)], 2, function(e){
  x = 1
  l = c()
  init = 1
  for(i in 2:length(e)){
    if(e[i] != e[i-1]){
      l[[x]] = (i-init)*sum(msp_fast[init:(i-1),4])
      x = x+1
      init = i
    }
  }
  l
}) ->
  msp_switches

sapply(msp_switches, function(e){
  length(unlist(e))
}) ->
  msp_num_switches

sapply(msp_switches, function(e){
  mean(unlist(e))
}) ->
  msp_mean_anc_seg_len

msp_acc_s_len = data.frame("acc" = accs[[1]], "s" = msp_num_switches, "len" = msp_mean_anc_seg_len)

apply(msp_combpub_50[,7:ncol(msp_combpub_50)], 2, function(e){
  x = 1
  l = c()
  init = 1
  for(i in 2:length(e)){
    if(e[i] != e[i-1]){
      l[[x]] = (i-init)*sum(msp_combpub_50[init:(i-1),4])
      x = x+1
      init = i
    }
  }
  l
}) ->
  msp_50_switches

sapply(msp_50_switches, function(e){
  length(unlist(e))
}) ->
  msp_50_num_switches

sapply(msp_50_switches, function(e){
  unlist(e) %>%
    mean
}) ->
  msp_50_mean_anc_seg_len
msp_50_mean_anc_seg_len[is.na(msp_50_mean_anc_seg_len)] = 1

msp_50_acc_s_len = data.frame("acc" = accs[[1]], "s" = msp_50_num_switches, "len" = msp_50_mean_anc_seg_len)

mean(acc_s_len$len)
mean(msp_acc_s_len$len)
mean(msp_50_acc_s_len$len)

unlist(switches) %>%
  .[.< 100000] %>%
  data.frame("win" = .) %>%
  ggplot(data = ., mapping = aes(x = win)) + 
  geom_histogram(binwidth = 10000) ->
  true_switches_l

unlist(msp_switches) %>%
  .[.< 100000] %>%
  data.frame("win" = .) %>%
  ggplot(data = ., mapping = aes(x = win)) + 
  geom_histogram(binwidth = 10000) ->
  comb_switches_l

unlist(msp_50_switches) %>%
  .[.< 100000] %>%
  data.frame("win" = .) %>%
  ggplot(data = ., mapping = aes(x = win)) + 
  geom_histogram(binwidth = 10000) ->
  comb_50_switches_l



unlist(switches) %>%
  .[.< 1000] %>%
  data.frame("win" = .) %>%
  ggplot(data = ., mapping = aes(x = win)) + 
  geom_histogram(binwidth = 10) ->
  true_switches_s

unlist(msp_switches) %>%
  .[.< 1000] %>%
  data.frame("win" = .) %>%
  ggplot(data = ., mapping = aes(x = win)) + 
  geom_histogram(binwidth = 10) ->
  comb_switches_s

unlist(msp_50_switches) %>%
  .[.< 1000] %>%
  data.frame("win" = .) %>%
  ggplot(data = ., mapping = aes(x = win)) + 
  geom_histogram(binwidth = 10) ->
  comb_50_switches_s

true_switches_l + comb_switches_l + comb_50_switches_l
true_switches_s + comb_switches_s + comb_50_switches_s
