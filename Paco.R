# Procrustes Application to Cophylogenetic Analysis (PACo)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tidytable)
library(ggpubr)

library(phyloseq)
library(ape)
library(paco)

library(future)
library(future.batchtools)
library(future.apply)
work_dir="/Users/caseymeili/Desktop/R"
threads = 6
my_seed = 8394

rescale_dist_mtx = function(m){
  m = m %>% as.matrix
  labs = m %>% colnames
  n_row = m %>% nrow
  n_col = m %>% ncol
  x = m %>% as.vector 
  x = scales::rescale(x) 
  m = matrix(x, nrow=n_row, ncol=n_col)
  colnames(m) = labs
  rownames(m) = labs
  m = m %>% as.dist
  return(m)
}
Host_Sp <-ape::read.tree(file="~/Desktop/R/Host.nwk")
Taxa <-ape::read.tree(file="~/Desktop/R/GeneraOG_mafft.tree")

library(readxl)
otu_mat <-read_excel("~/Desktop/R/Physeq.xlsx", sheet="OTU")

tax_mat<- read_excel("~/Desktop/R/Physeq.xlsx", sheet="taxon")

Meta <-read_excel("~/Desktop/R/Physeq.xlsx", sheet="Sample")

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("#OTU ID")
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("#OTU ID")
Meta <- Meta  %>% 
  tibble::column_to_rownames("Sample")
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(Meta)
Physeq <-phyloseq(OTU, TAX, samples, Taxa)
Physeq
library(HTSSIP)
# metadata
metadata = Physeq %>%
  phyloseq2df(sample_data) %>%
  dplyr::select(SampleNames, Animal) %>%
  as.data.frame

rownames(metadata) = metadata$SampleID
metadata[sample_names(Physeq),]
metadata

# number of subsampled trees
ntrees = 100

#' randomly selecting one per group
tree_subsample = function(L, df, tree){
  # get subsample (note: subsampling within each species)
  to_keep = df %>% 
    group_by(Animal) %>% 
    sample_n(1) %>%
    .$sample
  # subsampling tree
  to_rm = setdiff(tree$tip.label, to_keep)
  tree = drop.tip(tree, to_rm)
  return(tree)
}

# subsampling trees
df = metadata %>%
  mutate(sample = metadata$SampleNames) %>%
  dplyr::select(sample, Animal) 

doParallel::registerDoParallel(threads)
host_tree_l = plyr::llply(as.list(1:ntrees), 
                          function(x) tree_subsample(x, df, Host_Sp),
                          .parallel=TRUE)

# tree lengths
host_tree_l %>% length %>% print
lapply(host_tree_l, function(x) x$tip.label %>% length) %>% unlist %>% summary

make_paco_input = function(host_tree, physeq){
  # subsampled phyloseq object
  physeq_f = prune_samples(sample_names(physeq) %in% host_tree$tip.label, 
                           physeq) 
  
  # microbial distance matrix
  micro_D = physeq_f %>%
    phy_tree %>% cophenetic %>% 
    rescale_dist_mtx %>% as.matrix
  # host distance matrix
  host_D = host_tree %>% cophenetic %>% 
    rescale_dist_mtx %>% as.matrix
  
  # otu pres/abs matrix
  otu = physeq_f %>% 
    otu_table %>%
    as.data.frame %>%
    as.matrix %>%
    t %>%
    apply(2, function(x) ifelse(x > 0, 1, 0)) %>% as.matrix                       
  
  # checking overlap
  x = length(intersect(rownames(otu), rownames(host_D)))
  y = length(union(rownames(otu), rownames(host_D)))
  stopifnot(x == y)
  
  # preparing paco data
  D = prepare_paco_data(H=host_D, P=micro_D, HP=otu)
  D = add_pcoord(D, correction='cailliez')
  return(D)
}
doParallel::registerDoParallel(threads)
paco_l = plyr::llply(host_tree_l, make_paco_input, physeq=Physeq, .parallel=TRUE)
paco_l %>% length

paco_l[[1]] %>% .$H %>% dim
paco_l[[1]] %>% .$P %>% dim
paco_l[[1]] %>% .$HP %>% dim

paco_each = function(D, nperm=999, paco_seed=3843){
  PACo(D, nperm=nperm, seed=paco_seed, method='quasiswap', symmetric=TRUE)
}

paco_res_l = paco_l %>%
  future_lapply(paco_each, nperm=9, 
                future.packages=c('paco'), 
                future.seed=my_seed)
paco_res_l %>% length

gof = lapply(paco_res_l, function(x) as.vector(x$gof)) %>% 
  do.call(rbind, .) %>%
  as.data.frame

gof$p %>% as.Num %>% summary %>% print
gof$ss %>% as.Num %>% summary %>% print


D_links_l = future_lapply(paco_res_l, pacolinks_each, 
                          future.packages=c('paco'),
                          future.seed=my_seed)

# residuals
get_residuals = function(rep, links_l){
  # residuals
  res = residuals_paco(links_l[[rep]]$proc) %>% as.data.frame 
  colnames(res) = 'residuals'
  res = res %>%
    mutate(comparison = rownames(.),
           subsample_rep = rep) %>%
    separate(comparison, c('host', 'microbe'), sep='-') 
  
  # jackknife
  D_links_jk = do.call(rbind, as.list(links_l[[rep]]$jackknife)) %>%
    t %>% as.data.frame %>%
    mutate(comparison = rownames(.)) %>%
    separate(comparison, c('host', 'microbe'), sep='-') %>%
    inner_join(res, c('host'='host', 'microbe'='microbe'))
  
  # return
  return(D_links_jk)
}

## these below never worked so I had to improvise
D_links_l = lapply(as.list(1:length(D_links_l)), get_residuals, D_links_l=D_links_l)
D_links_l = do.call(rbind, D_links_l)
D_links_l

## how I improvised
## I extraced the $proc for each of the 100 trees and saved it as separate data frames
## You will then go and combined all these files into one file that you will save as Residuals.txt
res1 = residuals_paco(D_links_l[[1]]$proc) %>% as.data.frame
res1 = residuals_paco(D_links_l[[1]]$proc) %>% as.data.frame
res2 = residuals_paco(D_links_l[[2]]$proc) %>% as.data.frame
res3 = residuals_paco(D_links_l[[3]]$proc) %>% as.data.frame
res4 = residuals_paco(D_links_l[[4]]$proc) %>% as.data.frame
res5 = residuals_paco(D_links_l[[5]]$proc) %>% as.data.frame
res6 = residuals_paco(D_links_l[[6]]$proc) %>% as.data.frame
res7 = residuals_paco(D_links_l[[7]]$proc) %>% as.data.frame
res8 = residuals_paco(D_links_l[[8]]$proc) %>% as.data.frame
res9 = residuals_paco(D_links_l[[9]]$proc) %>% as.data.frame
res10 = residuals_paco(D_links_l[[10]]$proc) %>% as.data.frame
res11 = residuals_paco(D_links_l[[11]]$proc) %>% as.data.frame
res12 = residuals_paco(D_links_l[[12]]$proc) %>% as.data.frame
res13 = residuals_paco(D_links_l[[13]]$proc) %>% as.data.frame
res14 = residuals_paco(D_links_l[[14]]$proc) %>% as.data.frame
res15 = residuals_paco(D_links_l[[15]]$proc) %>% as.data.frame
res16 = residuals_paco(D_links_l[[16]]$proc) %>% as.data.frame
res17 = residuals_paco(D_links_l[[17]]$proc) %>% as.data.frame
res18 = residuals_paco(D_links_l[[18]]$proc) %>% as.data.frame
res19 = residuals_paco(D_links_l[[19]]$proc) %>% as.data.frame
res20 = residuals_paco(D_links_l[[20]]$proc) %>% as.data.frame
res21 = residuals_paco(D_links_l[[21]]$proc) %>% as.data.frame
res22 = residuals_paco(D_links_l[[22]]$proc) %>% as.data.frame
res23 = residuals_paco(D_links_l[[23]]$proc) %>% as.data.frame
res24 = residuals_paco(D_links_l[[24]]$proc) %>% as.data.frame
res25 = residuals_paco(D_links_l[[25]]$proc) %>% as.data.frame
res26 = residuals_paco(D_links_l[[26]]$proc) %>% as.data.frame
res27 = residuals_paco(D_links_l[[27]]$proc) %>% as.data.frame
res28 = residuals_paco(D_links_l[[28]]$proc) %>% as.data.frame
res29 = residuals_paco(D_links_l[[29]]$proc) %>% as.data.frame
res30 = residuals_paco(D_links_l[[30]]$proc) %>% as.data.frame
res31 = residuals_paco(D_links_l[[31]]$proc) %>% as.data.frame
res32 = residuals_paco(D_links_l[[32]]$proc) %>% as.data.frame
res33 = residuals_paco(D_links_l[[33]]$proc) %>% as.data.frame
res34 = residuals_paco(D_links_l[[34]]$proc) %>% as.data.frame
res35 = residuals_paco(D_links_l[[35]]$proc) %>% as.data.frame
res36 = residuals_paco(D_links_l[[36]]$proc) %>% as.data.frame
res37 = residuals_paco(D_links_l[[37]]$proc) %>% as.data.frame
res38 = residuals_paco(D_links_l[[38]]$proc) %>% as.data.frame
res39 = residuals_paco(D_links_l[[39]]$proc) %>% as.data.frame
res40 = residuals_paco(D_links_l[[40]]$proc) %>% as.data.frame
res41 = residuals_paco(D_links_l[[41]]$proc) %>% as.data.frame
res42 = residuals_paco(D_links_l[[42]]$proc) %>% as.data.frame
res43 = residuals_paco(D_links_l[[43]]$proc) %>% as.data.frame
res44 = residuals_paco(D_links_l[[44]]$proc) %>% as.data.frame
res45 = residuals_paco(D_links_l[[45]]$proc) %>% as.data.frame
res46 = residuals_paco(D_links_l[[46]]$proc) %>% as.data.frame
res47 = residuals_paco(D_links_l[[47]]$proc) %>% as.data.frame
res48 = residuals_paco(D_links_l[[48]]$proc) %>% as.data.frame
res49 = residuals_paco(D_links_l[[49]]$proc) %>% as.data.frame
res50 = residuals_paco(D_links_l[[50]]$proc) %>% as.data.frame
res51 = residuals_paco(D_links_l[[51]]$proc) %>% as.data.frame
res52 = residuals_paco(D_links_l[[52]]$proc) %>% as.data.frame
res53 = residuals_paco(D_links_l[[53]]$proc) %>% as.data.frame
res54 = residuals_paco(D_links_l[[54]]$proc) %>% as.data.frame
res55 = residuals_paco(D_links_l[[55]]$proc) %>% as.data.frame
res56 = residuals_paco(D_links_l[[56]]$proc) %>% as.data.frame
res57 = residuals_paco(D_links_l[[57]]$proc) %>% as.data.frame
res58 = residuals_paco(D_links_l[[58]]$proc) %>% as.data.frame
res59 = residuals_paco(D_links_l[[59]]$proc) %>% as.data.frame
res60 = residuals_paco(D_links_l[[60]]$proc) %>% as.data.frame
res61 = residuals_paco(D_links_l[[61]]$proc) %>% as.data.frame
res62 = residuals_paco(D_links_l[[62]]$proc) %>% as.data.frame
res63 = residuals_paco(D_links_l[[63]]$proc) %>% as.data.frame
res64 = residuals_paco(D_links_l[[64]]$proc) %>% as.data.frame
res65 = residuals_paco(D_links_l[[65]]$proc) %>% as.data.frame
res66 = residuals_paco(D_links_l[[66]]$proc) %>% as.data.frame
res67 = residuals_paco(D_links_l[[67]]$proc) %>% as.data.frame
res68 = residuals_paco(D_links_l[[68]]$proc) %>% as.data.frame
res69 = residuals_paco(D_links_l[[69]]$proc) %>% as.data.frame
res70 = residuals_paco(D_links_l[[70]]$proc) %>% as.data.frame
res71 = residuals_paco(D_links_l[[71]]$proc) %>% as.data.frame
res72 = residuals_paco(D_links_l[[72]]$proc) %>% as.data.frame
res73 = residuals_paco(D_links_l[[73]]$proc) %>% as.data.frame
res74 = residuals_paco(D_links_l[[74]]$proc) %>% as.data.frame
res75 = residuals_paco(D_links_l[[75]]$proc) %>% as.data.frame
res76 = residuals_paco(D_links_l[[76]]$proc) %>% as.data.frame
res77 = residuals_paco(D_links_l[[77]]$proc) %>% as.data.frame
res78 = residuals_paco(D_links_l[[78]]$proc) %>% as.data.frame
res79 = residuals_paco(D_links_l[[79]]$proc) %>% as.data.frame
res80 = residuals_paco(D_links_l[[80]]$proc) %>% as.data.frame
res81 = residuals_paco(D_links_l[[81]]$proc) %>% as.data.frame
res82 = residuals_paco(D_links_l[[82]]$proc) %>% as.data.frame
res83 = residuals_paco(D_links_l[[83]]$proc) %>% as.data.frame
res84 = residuals_paco(D_links_l[[84]]$proc) %>% as.data.frame
res85 = residuals_paco(D_links_l[[85]]$proc) %>% as.data.frame
res86 = residuals_paco(D_links_l[[86]]$proc) %>% as.data.frame
res87 = residuals_paco(D_links_l[[87]]$proc) %>% as.data.frame
res88 = residuals_paco(D_links_l[[88]]$proc) %>% as.data.frame
res89 = residuals_paco(D_links_l[[89]]$proc) %>% as.data.frame
res90 = residuals_paco(D_links_l[[90]]$proc) %>% as.data.frame
res91 = residuals_paco(D_links_l[[91]]$proc) %>% as.data.frame
res92 = residuals_paco(D_links_l[[92]]$proc) %>% as.data.frame
res93 = residuals_paco(D_links_l[[93]]$proc) %>% as.data.frame
res94 = residuals_paco(D_links_l[[94]]$proc) %>% as.data.frame
res95 = residuals_paco(D_links_l[[95]]$proc) %>% as.data.frame
res96 = residuals_paco(D_links_l[[96]]$proc) %>% as.data.frame
res97 = residuals_paco(D_links_l[[97]]$proc) %>% as.data.frame
res98 = residuals_paco(D_links_l[[98]]$proc) %>% as.data.frame
res99 = residuals_paco(D_links_l[[99]]$proc) %>% as.data.frame
res100 = residuals_paco(D_links_l[[100]]$proc) %>% as.data.frame
write.table(res1, file="res1")
write.table(res2, file="res2")
write.table(res3, file="res3")
write.table(res4, file="res4")
write.table(res5, file="res5")
write.table(res6, file="res6")
write.table(res7, file="res7")
write.table(res8, file="res8")
write.table(res9, file="res9")
write.table(res10, file="res10")
write.table(res11, file="res11")
write.table(res12, file="res12")
write.table(res13, file="res13")
write.table(res14, file="res14")
write.table(res15, file="res15")
write.table(res16, file="res16")
write.table(res17, file="res17")
write.table(res18, file="res18")
write.table(res19, file="res19")
write.table(res20, file="res20")
write.table(res21, file="res21")
write.table(res22, file="res22")
write.table(res23, file="res23")
write.table(res24, file="res24")
write.table(res25, file="res25")
write.table(res26, file="res26")
write.table(res27, file="res27")
write.table(res28, file="res28")
write.table(res29, file="res29")
write.table(res30, file="res30")
write.table(res31, file="res31")
write.table(res32, file="res32")
write.table(res33, file="res33")
write.table(res34, file="res34")
write.table(res35, file="res35")
write.table(res36, file="res36")
write.table(res37, file="res37")
write.table(res38, file="res38")
write.table(res39, file="res39")
write.table(res40, file="res40")
write.table(res41, file="res41")
write.table(res42, file="res42")
write.table(res43, file="res43")
write.table(res44, file="res44")
write.table(res45, file="res45")
write.table(res46, file="res46")
write.table(res47, file="res47")
write.table(res48, file="res48")
write.table(res49, file="res49")
write.table(res50, file="res50")
write.table(res51, file="res51")
write.table(res52, file="res52")
write.table(res53, file="res53")
write.table(res54, file="res54")
write.table(res55, file="res55")
write.table(res56, file="res56")
write.table(res57, file="res57")
write.table(res58, file="res58")
write.table(res59, file="res59")
write.table(res60, file="res60")
write.table(res61, file="res61")
write.table(res62, file="res62")
write.table(res63, file="res63")
write.table(res64, file="res64")
write.table(res65, file="res65")
write.table(res66, file="res66")
write.table(res67, file="res67")
write.table(res68, file="res68")
write.table(res69, file="res69")
write.table(res70, file="res70")
write.table(res71, file="res71")
write.table(res72, file="res72")
write.table(res73, file="res73")
write.table(res74, file="res74")
write.table(res75, file="res75")
write.table(res76, file="res76")
write.table(res77, file="res77")
write.table(res78, file="res78")
write.table(res79, file="res79")
write.table(res80, file="res80")
write.table(res81, file="res81")
write.table(res82, file="res82")
write.table(res83, file="res83")
write.table(res84, file="res84")
write.table(res85, file="res85")
write.table(res86, file="res86")
write.table(res87, file="res87")
write.table(res88, file="res88")
write.table(res89, file="res89")
write.table(res90, file="res90")
write.table(res91, file="res91")
write.table(res92, file="res92")
write.table(res93, file="res93")
write.table(res94, file="res94")
write.table(res95, file="res95")
write.table(res96, file="res96")
write.table(res97, file="res97")
write.table(res98, file="res98")
write.table(res99, file="res99")
write.table(res100, file="res100")
## again extract individual info per tree, write to a file, go and combine these files
## in one file called jackknifes.txt
D_links_jk1 = do.call(rbind, as.list(D_links_l[[1]]$jackknife)) %>% as.data.frame
D_links_jk2 = do.call(rbind, as.list(D_links_l[[2]]$jackknife)) %>% as.data.frame
D_links_jk3 = do.call(rbind, as.list(D_links_l[[3]]$jackknife)) %>% as.data.frame
D_links_jk4 = do.call(rbind, as.list(D_links_l[[4]]$jackknife)) %>% as.data.frame
D_links_jk5 = do.call(rbind, as.list(D_links_l[[5]]$jackknife)) %>% as.data.frame
D_links_jk6 = do.call(rbind, as.list(D_links_l[[6]]$jackknife)) %>% as.data.frame
D_links_jk7 = do.call(rbind, as.list(D_links_l[[7]]$jackknife)) %>% as.data.frame
D_links_jk8 = do.call(rbind, as.list(D_links_l[[8]]$jackknife)) %>% as.data.frame
D_links_jk9 = do.call(rbind, as.list(D_links_l[[9]]$jackknife)) %>% as.data.frame
D_links_jk10 = do.call(rbind, as.list(D_links_l[[10]]$jackknife)) %>% as.data.frame
D_links_jk11 = do.call(rbind, as.list(D_links_l[[11]]$jackknife)) %>% as.data.frame
D_links_jk12 = do.call(rbind, as.list(D_links_l[[12]]$jackknife)) %>% as.data.frame
D_links_jk13 = do.call(rbind, as.list(D_links_l[[13]]$jackknife)) %>% as.data.frame
D_links_jk14 = do.call(rbind, as.list(D_links_l[[14]]$jackknife)) %>% as.data.frame
D_links_jk15 = do.call(rbind, as.list(D_links_l[[15]]$jackknife)) %>% as.data.frame
D_links_jk16 = do.call(rbind, as.list(D_links_l[[16]]$jackknife)) %>% as.data.frame
D_links_jk17 = do.call(rbind, as.list(D_links_l[[17]]$jackknife)) %>% as.data.frame
D_links_jk18 = do.call(rbind, as.list(D_links_l[[18]]$jackknife)) %>% as.data.frame
D_links_jk19 = do.call(rbind, as.list(D_links_l[[19]]$jackknife)) %>% as.data.frame
D_links_jk20 = do.call(rbind, as.list(D_links_l[[20]]$jackknife)) %>% as.data.frame
D_links_jk21 = do.call(rbind, as.list(D_links_l[[21]]$jackknife)) %>% as.data.frame
D_links_jk22 = do.call(rbind, as.list(D_links_l[[22]]$jackknife)) %>% as.data.frame
D_links_jk23 = do.call(rbind, as.list(D_links_l[[23]]$jackknife)) %>% as.data.frame
D_links_jk24 = do.call(rbind, as.list(D_links_l[[24]]$jackknife)) %>% as.data.frame
D_links_jk25 = do.call(rbind, as.list(D_links_l[[25]]$jackknife)) %>% as.data.frame
D_links_jk26 = do.call(rbind, as.list(D_links_l[[26]]$jackknife)) %>% as.data.frame
D_links_jk27 = do.call(rbind, as.list(D_links_l[[27]]$jackknife)) %>% as.data.frame
D_links_jk28 = do.call(rbind, as.list(D_links_l[[28]]$jackknife)) %>% as.data.frame
D_links_jk29 = do.call(rbind, as.list(D_links_l[[29]]$jackknife)) %>% as.data.frame
D_links_jk30 = do.call(rbind, as.list(D_links_l[[30]]$jackknife)) %>% as.data.frame
D_links_jk31 = do.call(rbind, as.list(D_links_l[[31]]$jackknife)) %>% as.data.frame
D_links_jk32 = do.call(rbind, as.list(D_links_l[[32]]$jackknife)) %>% as.data.frame
D_links_jk33 = do.call(rbind, as.list(D_links_l[[33]]$jackknife)) %>% as.data.frame
D_links_jk34 = do.call(rbind, as.list(D_links_l[[34]]$jackknife)) %>% as.data.frame
D_links_jk35 = do.call(rbind, as.list(D_links_l[[35]]$jackknife)) %>% as.data.frame
D_links_jk36 = do.call(rbind, as.list(D_links_l[[36]]$jackknife)) %>% as.data.frame
D_links_jk37 = do.call(rbind, as.list(D_links_l[[37]]$jackknife)) %>% as.data.frame
D_links_jk38 = do.call(rbind, as.list(D_links_l[[38]]$jackknife)) %>% as.data.frame
D_links_jk39 = do.call(rbind, as.list(D_links_l[[39]]$jackknife)) %>% as.data.frame
D_links_jk40 = do.call(rbind, as.list(D_links_l[[40]]$jackknife)) %>% as.data.frame
D_links_jk41 = do.call(rbind, as.list(D_links_l[[41]]$jackknife)) %>% as.data.frame
D_links_jk42 = do.call(rbind, as.list(D_links_l[[42]]$jackknife)) %>% as.data.frame
D_links_jk43 = do.call(rbind, as.list(D_links_l[[43]]$jackknife)) %>% as.data.frame
D_links_jk44 = do.call(rbind, as.list(D_links_l[[44]]$jackknife)) %>% as.data.frame
D_links_jk45 = do.call(rbind, as.list(D_links_l[[45]]$jackknife)) %>% as.data.frame
D_links_jk46 = do.call(rbind, as.list(D_links_l[[46]]$jackknife)) %>% as.data.frame
D_links_jk47 = do.call(rbind, as.list(D_links_l[[47]]$jackknife)) %>% as.data.frame
D_links_jk48 = do.call(rbind, as.list(D_links_l[[48]]$jackknife)) %>% as.data.frame
D_links_jk49 = do.call(rbind, as.list(D_links_l[[49]]$jackknife)) %>% as.data.frame
D_links_jk50 = do.call(rbind, as.list(D_links_l[[50]]$jackknife)) %>% as.data.frame
D_links_jk51 = do.call(rbind, as.list(D_links_l[[51]]$jackknife)) %>% as.data.frame
D_links_jk52 = do.call(rbind, as.list(D_links_l[[52]]$jackknife)) %>% as.data.frame
D_links_jk53 = do.call(rbind, as.list(D_links_l[[53]]$jackknife)) %>% as.data.frame
D_links_jk54 = do.call(rbind, as.list(D_links_l[[54]]$jackknife)) %>% as.data.frame
D_links_jk55 = do.call(rbind, as.list(D_links_l[[55]]$jackknife)) %>% as.data.frame
D_links_jk56 = do.call(rbind, as.list(D_links_l[[56]]$jackknife)) %>% as.data.frame
D_links_jk57 = do.call(rbind, as.list(D_links_l[[57]]$jackknife)) %>% as.data.frame
D_links_jk58 = do.call(rbind, as.list(D_links_l[[58]]$jackknife)) %>% as.data.frame
D_links_jk59 = do.call(rbind, as.list(D_links_l[[59]]$jackknife)) %>% as.data.frame
D_links_jk60 = do.call(rbind, as.list(D_links_l[[60]]$jackknife)) %>% as.data.frame
D_links_jk61 = do.call(rbind, as.list(D_links_l[[61]]$jackknife)) %>% as.data.frame
D_links_jk62 = do.call(rbind, as.list(D_links_l[[62]]$jackknife)) %>% as.data.frame
D_links_jk63 = do.call(rbind, as.list(D_links_l[[63]]$jackknife)) %>% as.data.frame
D_links_jk64 = do.call(rbind, as.list(D_links_l[[64]]$jackknife)) %>% as.data.frame
D_links_jk65 = do.call(rbind, as.list(D_links_l[[65]]$jackknife)) %>% as.data.frame
D_links_jk66 = do.call(rbind, as.list(D_links_l[[66]]$jackknife)) %>% as.data.frame
D_links_jk67 = do.call(rbind, as.list(D_links_l[[67]]$jackknife)) %>% as.data.frame
D_links_jk68 = do.call(rbind, as.list(D_links_l[[68]]$jackknife)) %>% as.data.frame
D_links_jk69 = do.call(rbind, as.list(D_links_l[[69]]$jackknife)) %>% as.data.frame
D_links_jk70 = do.call(rbind, as.list(D_links_l[[70]]$jackknife)) %>% as.data.frame
D_links_jk71 = do.call(rbind, as.list(D_links_l[[71]]$jackknife)) %>% as.data.frame
D_links_jk72 = do.call(rbind, as.list(D_links_l[[72]]$jackknife)) %>% as.data.frame
D_links_jk73 = do.call(rbind, as.list(D_links_l[[73]]$jackknife)) %>% as.data.frame
D_links_jk74 = do.call(rbind, as.list(D_links_l[[74]]$jackknife)) %>% as.data.frame
D_links_jk75 = do.call(rbind, as.list(D_links_l[[75]]$jackknife)) %>% as.data.frame
D_links_jk76 = do.call(rbind, as.list(D_links_l[[76]]$jackknife)) %>% as.data.frame
D_links_jk77 = do.call(rbind, as.list(D_links_l[[77]]$jackknife)) %>% as.data.frame
D_links_jk78 = do.call(rbind, as.list(D_links_l[[78]]$jackknife)) %>% as.data.frame
D_links_jk79 = do.call(rbind, as.list(D_links_l[[79]]$jackknife)) %>% as.data.frame
D_links_jk80 = do.call(rbind, as.list(D_links_l[[80]]$jackknife)) %>% as.data.frame
D_links_jk81 = do.call(rbind, as.list(D_links_l[[81]]$jackknife)) %>% as.data.frame
D_links_jk82 = do.call(rbind, as.list(D_links_l[[82]]$jackknife)) %>% as.data.frame
D_links_jk83 = do.call(rbind, as.list(D_links_l[[83]]$jackknife)) %>% as.data.frame
D_links_jk84 = do.call(rbind, as.list(D_links_l[[84]]$jackknife)) %>% as.data.frame
D_links_jk85 = do.call(rbind, as.list(D_links_l[[85]]$jackknife)) %>% as.data.frame
D_links_jk86 = do.call(rbind, as.list(D_links_l[[86]]$jackknife)) %>% as.data.frame
D_links_jk87 = do.call(rbind, as.list(D_links_l[[87]]$jackknife)) %>% as.data.frame
D_links_jk88 = do.call(rbind, as.list(D_links_l[[88]]$jackknife)) %>% as.data.frame
D_links_jk89 = do.call(rbind, as.list(D_links_l[[89]]$jackknife)) %>% as.data.frame
D_links_jk90 = do.call(rbind, as.list(D_links_l[[90]]$jackknife)) %>% as.data.frame
D_links_jk91 = do.call(rbind, as.list(D_links_l[[91]]$jackknife)) %>% as.data.frame
D_links_jk92 = do.call(rbind, as.list(D_links_l[[92]]$jackknife)) %>% as.data.frame
D_links_jk93 = do.call(rbind, as.list(D_links_l[[93]]$jackknife)) %>% as.data.frame
D_links_jk94 = do.call(rbind, as.list(D_links_l[[94]]$jackknife)) %>% as.data.frame
D_links_jk95 = do.call(rbind, as.list(D_links_l[[95]]$jackknife)) %>% as.data.frame
D_links_jk96 = do.call(rbind, as.list(D_links_l[[96]]$jackknife)) %>% as.data.frame
D_links_jk97 = do.call(rbind, as.list(D_links_l[[97]]$jackknife)) %>% as.data.frame
D_links_jk98 = do.call(rbind, as.list(D_links_l[[98]]$jackknife)) %>% as.data.frame
D_links_jk99 = do.call(rbind, as.list(D_links_l[[99]]$jackknife)) %>% as.data.frame
D_links_jk100 = do.call(rbind, as.list(D_links_l[[100]]$jackknife)) %>% as.data.frame
write.table(D_links_jk1, file="D_links_jk1")
write.table(D_links_jk2, file="D_links_jk2")
write.table(D_links_jk3, file="D_links_jk3")
write.table(D_links_jk4, file="D_links_jk4")
write.table(D_links_jk5, file="D_links_jk5")
write.table(D_links_jk6, file="D_links_jk6")
write.table(D_links_jk7, file="D_links_jk7")
write.table(D_links_jk8, file="D_links_jk8")
write.table(D_links_jk9, file="D_links_jk9")
write.table(D_links_jk10, file="D_links_jk10")
write.table(D_links_jk11, file="D_links_jk11")
write.table(D_links_jk12, file="D_links_jk12")
write.table(D_links_jk13, file="D_links_jk13")
write.table(D_links_jk14, file="D_links_jk14")
write.table(D_links_jk15, file="D_links_jk15")
write.table(D_links_jk16, file="D_links_jk16")
write.table(D_links_jk17, file="D_links_jk17")
write.table(D_links_jk18, file="D_links_jk18")
write.table(D_links_jk19, file="D_links_jk19")
write.table(D_links_jk20, file="D_links_jk20")
write.table(D_links_jk21, file="D_links_jk21")
write.table(D_links_jk22, file="D_links_jk22")
write.table(D_links_jk23, file="D_links_jk23")
write.table(D_links_jk24, file="D_links_jk24")
write.table(D_links_jk25, file="D_links_jk25")
write.table(D_links_jk26, file="D_links_jk26")
write.table(D_links_jk27, file="D_links_jk27")
write.table(D_links_jk28, file="D_links_jk28")
write.table(D_links_jk29, file="D_links_jk29")
write.table(D_links_jk30, file="D_links_jk30")
write.table(D_links_jk31, file="D_links_jk31")
write.table(D_links_jk32, file="D_links_jk32")
write.table(D_links_jk33, file="D_links_jk33")
write.table(D_links_jk34, file="D_links_jk34")
write.table(D_links_jk35, file="D_links_jk35")
write.table(D_links_jk36, file="D_links_jk36")
write.table(D_links_jk37, file="D_links_jk37")
write.table(D_links_jk38, file="D_links_jk38")
write.table(D_links_jk39, file="D_links_jk39")
write.table(D_links_jk40, file="D_links_jk40")
write.table(D_links_jk41, file="D_links_jk41")
write.table(D_links_jk42, file="D_links_jk42")
write.table(D_links_jk43, file="D_links_jk43")
write.table(D_links_jk44, file="D_links_jk44")
write.table(D_links_jk45, file="D_links_jk45")
write.table(D_links_jk46, file="D_links_jk46")
write.table(D_links_jk47, file="D_links_jk47")
write.table(D_links_jk48, file="D_links_jk48")
write.table(D_links_jk49, file="D_links_jk49")
write.table(D_links_jk50, file="D_links_jk50")
write.table(D_links_jk51, file="D_links_jk51")
write.table(D_links_jk52, file="D_links_jk52")
write.table(D_links_jk53, file="D_links_jk53")
write.table(D_links_jk54, file="D_links_jk54")
write.table(D_links_jk55, file="D_links_jk55")
write.table(D_links_jk56, file="D_links_jk56")
write.table(D_links_jk57, file="D_links_jk57")
write.table(D_links_jk58, file="D_links_jk58")
write.table(D_links_jk59, file="D_links_jk59")
write.table(D_links_jk60, file="D_links_jk60")
write.table(D_links_jk61, file="D_links_jk61")
write.table(D_links_jk62, file="D_links_jk62")
write.table(D_links_jk63, file="D_links_jk63")
write.table(D_links_jk64, file="D_links_jk64")
write.table(D_links_jk65, file="D_links_jk65")
write.table(D_links_jk66, file="D_links_jk66")
write.table(D_links_jk67, file="D_links_jk67")
write.table(D_links_jk68, file="D_links_jk68")
write.table(D_links_jk69, file="D_links_jk69")
write.table(D_links_jk70, file="D_links_jk70")
write.table(D_links_jk71, file="D_links_jk71")
write.table(D_links_jk72, file="D_links_jk72")
write.table(D_links_jk73, file="D_links_jk73")
write.table(D_links_jk74, file="D_links_jk74")
write.table(D_links_jk75, file="D_links_jk75")
write.table(D_links_jk76, file="D_links_jk76")
write.table(D_links_jk77, file="D_links_jk77")
write.table(D_links_jk78, file="D_links_jk78")
write.table(D_links_jk79, file="D_links_jk79")
write.table(D_links_jk80, file="D_links_jk80")
write.table(D_links_jk81, file="D_links_jk81")
write.table(D_links_jk82, file="D_links_jk82")
write.table(D_links_jk83, file="D_links_jk83")
write.table(D_links_jk84, file="D_links_jk84")
write.table(D_links_jk85, file="D_links_jk85")
write.table(D_links_jk86, file="D_links_jk86")
write.table(D_links_jk87, file="D_links_jk87")
write.table(D_links_jk88, file="D_links_jk88")
write.table(D_links_jk89, file="D_links_jk89")
write.table(D_links_jk90, file="D_links_jk90")
write.table(D_links_jk91, file="D_links_jk91")
write.table(D_links_jk92, file="D_links_jk92")
write.table(D_links_jk93, file="D_links_jk93")
write.table(D_links_jk94, file="D_links_jk94")
write.table(D_links_jk95, file="D_links_jk95")
write.table(D_links_jk96, file="D_links_jk96")
write.table(D_links_jk97, file="D_links_jk97")
write.table(D_links_jk98, file="D_links_jk98")
write.table(D_links_jk99, file="D_links_jk99")
write.table(D_links_jk100, file="D_links_jk100")

## on your saved Jackknifes.txt add a header of Host, microbe, residuals
## read in your saved file
D_links_list=read.table("~/Desktop/GradStudents/Casey-Adrienne/FecesPaper/D_links_l.txt", header=TRUE)
D_links_list

D_links_list = D_links_list %>%
  group_by(Host) %>%
  summarize(mean_resid = mean(residuals),
            median_resid = median(residuals),
            sd_resid = sd(residuals),
            CV_resid = sd_resid / mean_resid * 100) %>%
  ungroup()
D_links_list

# host taxonomy
host_tax = Physeq %>% sample_data %>% 
  as.matrix %>% as.data.frame %>%
  mutate(sample = rownames(.)) %>%
  dplyr::select(sample, Animal, Family, Country, Lifestyle)
host_tax

D_links_list = D_links_list %>%
  inner_join(host_tax, c('Host'='sample'))
write.table(D_links_list, file="~/Desktop/GradStudents/Casey-Adrienne/FecesPaper/D_links_l_tax.txt")

## next we plot box plots of residuals
# by Animal
p = ggplot(D_links_list, aes(mean_resid, sd_resid, color=Animal)) +
  geom_point(alpha=0.2) +
  geom_point(alpha=0.5, shape='O') +
  theme_bw() 

p.dims(6,4)
plot(p)

# by host Family
p = ggplot(D_links_list, aes(mean_resid, sd_resid, color=Family)) +
  geom_point(alpha=0.2) +
  geom_point(alpha=0.5, shape='O') +
  theme_bw() 

p.dims(6,4)
plot(p)

# by Country
p = ggplot(D_links_list, aes(mean_resid, sd_resid, color=Country)) +
  geom_point(alpha=0.2) +
  geom_point(alpha=0.5, shape='O') +
  theme_bw() 

p.dims(7,4)
plot(p)

# by host Lifestyle
p = ggplot(D_links_list, aes(mean_resid, sd_resid, color=Lifestyle)) +
  geom_point(alpha=0.2) +
  geom_point(alpha=0.5, shape='O') +
  theme_bw() 

p.dims(7,4)
plot(p)


# We can attempt to add signficance to the plots as well
## ordering
tmp = D_links_list %>% 
  group_by(Animal) %>%
  mutate(median_resid = median(mean_resid)) %>%
  ungroup() %>%
  mutate(Animal = Animal %>% reorder(-median_resid))

p = ggplot(tmp, aes(Animal, mean_resid)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(x='Animal', y='Residuals') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1)
  )

p.dims(3.5,2.5)
plot(p)


# host Animal
tmp = D_links_list %>%
  group_by(Host, Animal) %>%
  summarize(mean_resid = mean(mean_resid)) %>%
  ungroup()

# stats
df.dims(20)
means=compare_means(mean_resid ~ Animal, data=tmp, method = 'wilcox.test', p.adjust.method = 'BH')
means
write.table(means, file="MeansAnimal.txt")

## plotting (change the list below to reflect the list of animals you have)
my_comparisons = combn(as.character(unique(tmp$Animal)), 2, simplify=FALSE) 
levs = c('Bison', 'Buffalo', 'Cow', 'Deer', 'Goat','Oryx', 'Sheep', 'Yak', 'Donkey', 'Elephant', 'Horse', 'Mule', 'Rhinoceros', 'Zebra', 'Alpaca', 'Camel')
p = tmp %>% 
  mutate(Animal = factor(Animal, levels=levs)) %>%
  ggboxplot(x = "Animal", y = "mean_resid", palette = "jco") + 
  stat_compare_means(comparisons = my_comparisons, 
                     label = 'p.signif', 
                     hide.ns = TRUE, vjust=0.5) +
  labs(x='Animal', y='Residuals') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1)
  )
p.dims(3.5,2.5)
plot(p)

# host Family
tmp = D_links_list %>%
  group_by(Host, Family) %>%
  summarize(mean_resid = mean(mean_resid)) %>%
  ungroup()

# stats
df.dims(20)
means=compare_means(mean_resid ~ Family, data=tmp, method = 'wilcox.test', p.adjust.method = 'BH')
means
write.table(means, file="MeansFam.txt")

## plotting (change the list below to reflect the list of families you have)
my_comparisons = combn(as.character(unique(tmp$Family)), 2, simplify=FALSE) 
levs = c('Bovidae', 'Cervidae', 'Giraffidae', 'Camelidae', 'Equidae', 'Caviidae', 'Elephantidae', 'Trichechidae', 'Rhinocerotidae')
p = tmp %>% 
  ggboxplot(x = "Family", y = "mean_resid", palette = "jco") + 
  stat_compare_means(comparisons = my_comparisons, 
                     label = 'p.signif', 
                     hide.ns = TRUE, vjust=0.5) +
  labs(x='Family', y='Residuals') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
p.dims(3.5,2.5)
plot(p)

# host Country
tmp = D_links_list %>%
  group_by(Host, Country) %>%
  summarize(mean_resid = mean(mean_resid)) %>%
  ungroup()

# stats
df.dims(20)
means=compare_means(mean_resid ~ Country, data=tmp, method = 'wilcox.test', p.adjust.method = 'BH')
means
write.table(means, file="MeansCountry.txt")

## plotting (change the list below to reflect the list of countries you have)
my_comparisons = combn(as.character(unique(tmp$Country)), 2, simplify=FALSE) 
levs = c('Foregut', 'Hindgut', 'Pseudoruminant')
p = tmp %>% 
  ggboxplot(x = "Country", y = "mean_resid", palette = "jco") + 
  stat_compare_means(comparisons = my_comparisons, 
                     label = 'p.signif', 
                     hide.ns = TRUE, vjust=0.5) +
  labs(x='Country', y='Residuals') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1)
  )
p.dims(3.5,2.5)
plot(p)


# host LS
tmp = D_links_list %>%
  group_by(Host, Lifestyle) %>%
  summarize(mean_resid = mean(mean_resid)) %>%
  ungroup()

# stats
df.dims(20)
means=compare_means(mean_resid ~ Lifestyle, data=tmp, method = 'wilcox.test', p.adjust.method = 'BH')
means
write.table(means, file="MeansLS.txt")

## plotting
my_comparisons = combn(as.character(unique(tmp$Lifestyle)), 2, simplify=FALSE) 
levs = c('Domesticated', 'Non-domesticated')
p = tmp %>% 
  ggboxplot(x = "Lifestyle", y = "mean_resid", palette = "jco") + 
  stat_compare_means(comparisons = my_comparisons, 
                     label = 'p.signif', 
                     hide.ns = TRUE, vjust=0.5) +
  labs(x='Lifestyle', y='Residuals') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1)
  )
p.dims(3.5,2.5)
plot(p)
