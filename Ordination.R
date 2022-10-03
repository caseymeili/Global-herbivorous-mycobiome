#  Ordination, anosim, adonis, and MRM (with and without subsampling)


library(readxl)
library(phyloseq)
library(ape)
library(plyr)
library(vegan)
library(rbiom)

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

dist_mtx_order = function(d, x){
  m = d %>% as.matrix
  d = as.dist(m[x,x])
  return(d)
}

## create phyloseq object with the genera tree
otu_mat <-read_excel("~/Desktop/GradStudents/Casey-Adrienne/FecesPaper/Physeq.xlsx", sheet="OTU")
tax_mat<- read_excel("~/Desktop/GradStudents/Casey-Adrienne/FecesPaper/Physeq.xlsx", sheet="taxon")
Meta <-read_excel("~/Desktop/GradStudents/Casey-Adrienne/FecesPaper/Physeq.xlsx", sheet="Sample")
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("#OTU ID")
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("#OTU ID")
Meta <- metadata  %>%
  tibble::column_to_rownames("Sample")
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(metadata)
samples
Tree <-ape::read.tree(file="~/Desktop/GradStudents/Casey-Adrienne/FecesPaper/Genera_rooted.nwk")

Physeq <-phyloseq(OTU, TAX, samples, Tree)
Physeq

## create ordination plots

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="Animal", shape="Family")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
BrayOrd <-ggplot(pdataframe, aes(Axis_1, Axis_2, color=Animal, shape=Family))
BrayOrd=BrayOrd+geom_point(size=2)
BrayOrd=BrayOrd+facet_wrap(~method, scales="free")
BrayOrd

dist="jaccard"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="Animal", shape="Family")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
JaccardOrd <-ggplot(pdataframe, aes(Axis_1, Axis_2, color=Animal, shape=Family))
JaccardOrd=JaccardOrd+geom_point(size=2)
JaccardOrd=JaccardOrd+facet_wrap(~method, scales="free")
JaccardOrd

dist="euclidean"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="Animal", shape="Family")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
EucOrd <-ggplot(pdataframe, aes(Axis_1, Axis_2, color=Animal, shape=Family))
EucOrd=EucOrd+geom_point(size=2)
EucOrd=EucOrd+facet_wrap(~method, scales="free")
EucOrd

dist="unifrac"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="Animal", shape="Family")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
UniOrd <-ggplot(pdataframe, aes(Axis_1, Axis_2, color=Animal, shape=Family))
UniOrd=UniOrd+geom_point(size=2)
UniOrd=UniOrd+facet_wrap(~method, scales="free")
UniOrd

dist="Wunifrac"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, Physeq, dist){
  ordi = ordinate(Physeq, method=i, distance=dist)
  plot_ordination(Physeq, ordi, "samples", color="Animal", shape="GutType")
}, Physeq, dist)
names(plist) <- ord_meths
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
WUniOrd <-ggplot(pdataframe, aes(Axis_1, Axis_2, color=Animal, shape=GutType))
WUniOrd=WUniOrd+geom_point(size=2)
WUniOrd=WUniOrd+facet_wrap(~method, scales="free")
WUniOrd

## reading in host tree, converting it to dist
Host_Sp <-ape::read.tree(file="~/Desktop/GradStudents/Casey-Adrienne/FecesPaper/Host.nwk")
host_tree_d = Host_Sp %>% cophenetic %>% as.dist %>% rescale_dist_mtx

X = labels(host_tree_d)

## creating beta diversity matrices then ordering them
## DM is the abundance file
## factors file contains the variables
## Biom is the abundance file but flipped with genera as rows and samples as columns
DM <-read.table("~/Desktop/GradStudents/Casey-Adrienne/WorkingUnclass/DM.txt", header=TRUE, row.names=1)
Factors <-read.table("~/Desktop/GradStudents/Casey-Adrienne/WorkingUnclass/Factors.txt", header=TRUE)
Biom <-read.table("~/Desktop/GradStudents/Casey-Adrienne/FecesPaper/biom.txt", header=TRUE, row.names=1)
Biom <-as.matrix(Biom)

Bray <-vegdist(DM, method="bray")
Bray_o <-dist_mtx_order(Bray, X)
Euc <-vegdist(DM, method="euclidean")
Euc_o <-dist_mtx_order(Euc, X)
Jac <-vegdist(DM, method="jaccard")
Jac_o <-dist_mtx_order(Jac, X)
Unifrac_uw <-rbiom::unifrac(Biom, weighted=FALSE, tree=Tree)
Uni_uw_o <-dist_mtx_order(Unifrac_uw, X)
Unifrac_w <-unifrac(Biom, weighted=TRUE, tree=Tree)
Uni_w_o <-dist_mtx_order(Unifrac_w, X)

## Testing for significance of factors
## Anosim

Anosim1 <-anosim(Bray_o, Factors$Animal)
Anosim2 <-anosim(Euc_o, Factors$Animal)
Anosim3 <-anosim(Jac_o, Factors$Animal)
Anosim4 <-anosim(Unifrac_uw_o, Factors$Animal)
Anosim5 <-anosim(Unifrac_w_o, Factors$Animal)

Anosim6 <-anosim(Bray_o, Factors$Family)
Anosim7 <-anosim(Euc_o, Factors$Family)
Anosim8 <-anosim(Jac_o, Factors$Family)
Anosim9 <-anosim(Unifrac_uw_o, Factors$Family)
Anosim10 <-anosim(Unifrac_w_o, Factors$Family)

Anosim11 <-anosim(Bray_o, Factors$Lifestyle)
Anosim12 <-anosim(Euc_o, Factors$Lifestyle)
Anosim13 <-anosim(Jac_o, Factors$Lifestyle)
Anosim14 <-anosim(Unifrac_uw_o, Factors$Lifestyle)
Anosim15 <-anosim(Unifrac_w_o, Factors$Lifestyle)

Anosim16 <-anosim(Bray_o, Factors$Country)
Anosim17 <-anosim(Euc_o, Factors$Country)
Anosim18 <-anosim(Jac_o, Factors$Country)
Anosim19 <-anosim(Unifrac_uw_o, Factors$Country)
Anosim20 <-anosim(Unifrac_w_o, Factors$Country)

Anosim1
Anosim2
Anosim3
Anosim4
Anosim5
Anosim6
Anosim7
Anosim8
Anosim9
Anosim10
Anosim11
Anosim12
Anosim13
Anosim14
Anosim15
Anosim16
Anosim17
Anosim18
Anosim19
Anosim20

##adonis
Adonis1 <-adonis(Bray_o ~Animal, Factors)
Adonis2 <-adonis(Euc_o ~Animal, Factors)
Adonis3 <-adonis(Jac_o ~Animal, Factors)
Adonis4 <-adonis(Unifrac_uw_o ~Animal, Factors)
Adonis5 <-adonis(Unifrac_w_o ~Animal, Factors)

Adonis6 <-adonis(Bray_o ~Family, Factors)
Adonis7 <-adonis(Euc_o ~Family, Factors)
Adonis8 <-adonis(Jac_o ~Family, Factors)
Adonis9 <-adonis(Unifrac_uw_o ~Family, Factors)
Adonis10 <-adonis(Unifrac_w_o ~Family, Factors)

Adonis11 <-adonis(Bray_o ~Lifestyle, Factors)
Adonis12 <-adonis(Euc_o ~Lifestyle, Factors)
Adonis13 <-adonis(Jac_o ~Lifestyle, Factors)
Adonis14 <-adonis(Unifrac_uw_o ~Lifestyle, Factors)
Adonis15 <-adonis(Unifrac_w_o ~Lifestyle, Factors)

Adonis16 <-adonis(Bray_o ~Country, Factors)
Adonis17 <-adonis(Euc_o ~Country, Factors)
Adonis18 <-adonis(Jac_o ~Country, Factors)
Adonis19 <-adonis(Unifrac_uw_o ~Country, Factors)
Adonis20 <-adonis(Unifrac_w_o ~Country, Factors)

Adonis1$aov
Adonis2$aov
Adonis3$aov
Adonis4$aov
Adonis5$aov
Adonis6$aov
Adonis7$aov
Adonis8$aov
Adonis9$aov
Adonis10$aov
Adonis11$aov
Adonis12$aov
Adonis13$aov
Adonis14$aov
Adonis15$aov
Adonis16$aov
Adonis17$aov
Adonis18$aov
Adonis19$aov
Adonis20$aov

## MRM and mantel
## creating factors distance matrices
Lifestyle_d <-vegan::vegdist(Factors$Lifestyle, metric="gower")
Lifestyle_d_o <-dist_mtx_order(Lifestyle_d, X)
Country_d <-vegan::vegdist(Factors$Country, metric="gower")
Country_d_o <-dist_mtx_order(Gut_d, X)
Family_d <-vegan::vegdist(Factors$Family, metric="gower")
Family_d_o <-dist_mtx_order(Family_d, X)
Lifestyle_d_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Country_d_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Family_d_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Uni_uw_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Uni_w_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Bray_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Euc_o %>% lapply(function(x) x %>% as.matrix %>% dim)
Jac_o %>% lapply(function(x) x %>% as.matrix %>% dim)

MRM_Euc <-ecodist::MRM(Euc_o ~ host_tree_d+Family_d_o+Country_d_o+Lifestyle_d_o)
MRM_Euc
MRM_Jac <-ecodist::MRM(Jac_o ~ host_tree_d+Family_d_o+Country_d_o+Lifestyle_d_o)
MRM_Jac
MRM_Bray <-ecodist::MRM(Bray_o ~ host_tree_d+Family_d_o+Country_d_o+Lifestyle_d_o)
MRM_Bray
MRM_Uni_uw <-ecodist::MRM(Uni_uw_o ~ host_tree_d+Family_d_o+Country_d_o+Lifestyle_d_o)
MRM_Uni_uw
MRM_Uni_w <-ecodist::MRM(Uni_w_o ~ host_tree_d+Family_d_o+Country_d_o+Lifestyle_d_o)
MRM_Uni_w
Mantel_Euc <-ecodist::mantel(Euc_o ~ host_tree_d)
Mantel_Euc
Mantel_Euc <-ecodist::mantel(Euc_o ~ Family_d_o)
Mantel_Euc
Mantel_Euc <-ecodist::mantel(Euc_o ~ Country_d_o)
Mantel_Euc
Mantel_Euc <-ecodist::mantel(Euc_o ~ Lifestyle_d_o)
Mantel_Euc
Mantel_Jac <-ecodist::mantel(Jac_o ~ host_tree_d)
Mantel_Jac
Mantel_Jac <-ecodist::mantel(Jac_o ~ Family_d_o)
Mantel_Jac
Mantel_Jac <-ecodist::mantel(Jac_o ~ Country_d_o)
Mantel_Jac
Mantel_Jac <-ecodist::mantel(Jac_o ~ Lifestyle_d_o)
Mantel_Jac
Mantel_Bray <-ecodist::mantel(Bray_o ~ host_tree_d)
Mantel_Bray
Mantel_Bray <-ecodist::mantel(Bray_o ~ Family_d_o)
Mantel_Bray
Mantel_Bray <-ecodist::mantel(Bray_o ~ Country_d_o)
Mantel_Bray
Mantel_Bray <-ecodist::mantel(Bray_o ~ Lifestyle_d_o)
Mantel_Bray
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ host_tree_d)
Mantel_Uni_uw
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ Family_d_o)
Mantel_Uni_uw
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ Country_d_o)
Mantel_Uni_uw
Mantel_Uni_uw <-ecodist::mantel(Uni_uw_o ~ Lifestyle_d_o)
Mantel_Uni_uw
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ host_tree_d)
Mantel_Uni_w
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ Family_d_o)
Mantel_Uni_w
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ Country_d_o)
Mantel_Uni_w
Mantel_Uni_w <-ecodist::mantel(Uni_w_o ~ Lifestyle_d_o)
Mantel_Uni_w

## subsampling for MRM
## the headers SampleNames, and Animal need to match the headers in your samples sheet of the phyloseq excel
## need the header of the column that has the sample name, and the header of the column that has the animal genus
## change dplyr::select(SampleNames, Animal) below to match these two headers
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
## need to change group_by(Animal) to the name of the header that contains the animal genus info

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
# metadata$SampleNames can be changed to the header that has the names of the samples
df = metadata %>%
  mutate(sample = metadata$SampleNames) %>%
  dplyr::select(sample, Animal) 
## this following will create a file with 100 trees each with just one from each type of animal
doParallel::registerDoParallel(threads)
host_tree_l = plyr::llply(as.list(1:ntrees), 
                          function(x) tree_subsample(x, df, Host_Sp),
                          .parallel=TRUE)

# tree lengths (should be 100)
# the output of the second line below should be equal to the number of different animal types you have
host_tree_l %>% length %>% print
lapply(host_tree_l, function(x) x$tip.label %>% length) %>% unlist %>% summary

# number of permutated SpecD datasetst
nperm_datasets = 100
# number of permutations per MRM analysis
nperm = 1000
#' randomly selecting one per group
#' L : list of distance matrixes used for MRM
#' df_grps : data.frame (sample, group)
one_per_group = function(L, df_grps, ...){
  # get subsample
  colnames(df_grps) = c('sample', 'Animal')
  df_grps = df_grps %>%
    group_by(Animal) %>%
    sample_n(1)
  # subsetting all matrices
  lapply(L, function(x) dist_mtx_order(x, df_grps$sample))
}

#' MRM on one subsample rep
#' i : rep number
#' L : list of list of distance matrices generated by `one_per_group()`
# nperm : nperm function for MRM
# f : MRM fomulat
mrm_each = function(i, L, f, nperm=99){
  m = L[[i]]
  f = as.formula(f)
  x = ecodist::MRM(f, nperm=nperm, mrank=TRUE)
  # coefficients
  df = x$coef %>% as.data.frame
  colnames(df) = c('coef', 'coef_pval')
  df$variable = rownames(df)
  df$R2 = x$r.squared[1]
  df$pval = x$r.squared[2]
  df$F = x$F.test[1]
  df$rep = i
  return(df)
}

# creating subsample permutations of the distance matrices (starting with Jaccard)
L = list(beta = Jac_o,
         host_phy = host_tree_d,
         Countrty = Country_d_o,
         Lifestyle = Lifestyle_d_o, 
         Family = Family_d_o)

m_perm = lapply(as.list(1:nperm_datasets), function(x) one_per_group(L, df, x))

# MRM on each permutation (in parallel)
doParallel::registerDoParallel(threads)
x = as.list(1:length(m_perm))
f='m$beta ~m$host_phy+ m$Gut+ m$Lifestyle+ m$Family'
mrm_res = plyr::llply(x, mrm_each, L=m_perm, f=f, nperm=nperm, .parallel=TRUE)
mrm_res = do.call(rbind, mrm_res)
mrm_res
mrm_res %>%
  filter(variable != 'Int') %>% 
  distinct(coef, coef_pval, R2, pval, rep) %>%
  summary

mrm_res %>%
  filter(variable != 'Int') %>%
  mutate(pval = pval %>% as.Num) %>%
  group_by(variable) %>%
  summarize(overall_pval = 1 - (sum(coef_pval < 0.05) / length(coef_pval))) %>%
  ungroup()


# formatting results
mrm_res_s = mrm_res %>% 
  filter(variable != 'Int') %>%
  gather(category, value, -variable, -R2, -pval, -F, -rep) %>%
  mutate(category = ifelse(category == 'coef', 'Coef.', 'Adj. P-value'))

mrm_res_s = mrm_res_s %>%
  inner_join(rename_df, c('variable'='old_name')) %>%
  dplyr::select(-variable) %>%
  rename('variable' = new_name)           

# plotting
X = data.frame(yint=c(0.05, NA), category=c('Adj. P-value', 'Coef.'))
p = ggplot(mrm_res_s, aes(variable, value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept=yint), linetype='dashed', alpha=0.3, data=X) +
  labs(x='Variable', y='Intra-species variance') +
  facet_grid(category ~ ., scales='free_y') +
  theme_bw()

p.dims(5,3)
plot(p)
write.table(mrm_res, file="mrm_res_Jac")

## below is the same set of commands but repeated by beta diversity measure
# creating subsample permutations of the distance matrices
L = list(beta = Bray_o,
         host_phy = host_tree_d,
         Country = Country_d_o,
         Lifestyle = Lifestyle_d_o, 
         Family = Family_d_o)

m_perm = lapply(as.list(1:nperm_datasets), function(x) one_per_group(L, df, x))

# MRM on each permutation (in parallel)
doParallel::registerDoParallel(threads)
x = as.list(1:length(m_perm))
f='m$beta ~m$host_phy+ m$Country+ m$Lifestyle+ m$Family'
mrm_res = plyr::llply(x, mrm_each, L=m_perm, f=f, nperm=nperm, .parallel=TRUE)
mrm_res = do.call(rbind, mrm_res)
mrm_res
mrm_res %>%
  filter(variable != 'Int') %>% 
  distinct(coef, coef_pval, R2, pval, rep) %>%
  summary
# significant?
df.dims(10)
mrm_res %>%
  filter(variable != 'Int') %>%
  mutate(pval = pval %>% as.Num) %>%
  group_by(variable) %>%
  summarize(overall_pval = 1 - (sum(coef_pval < 0.05) / length(coef_pval))) %>%
  ungroup()
df.dims()

# formatting results
mrm_res_s = mrm_res %>% 
  filter(variable != 'Int') %>%
  gather(category, value, -variable, -R2, -pval, -F, -rep) %>%
  mutate(category = ifelse(category == 'coef', 'Coef.', 'Adj. P-value'))

mrm_res_s = mrm_res_s %>%
  inner_join(rename_df, c('variable'='old_name')) %>%
  dplyr::select(-variable) %>%
  rename('variable' = new_name)           

# plotting
X = data.frame(yint=c(0.05, NA), category=c('Adj. P-value', 'Coef.'))
p = ggplot(mrm_res_s, aes(variable, value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept=yint), linetype='dashed', alpha=0.3, data=X) +
  labs(x='Variable', y='Intra-species variance') +
  facet_grid(category ~ ., scales='free_y') +
  theme_bw()

p.dims(5,3)
plot(p)
write.table(mrm_res, file="mrm_res_Bray")

# creating subsample permutations of the distance matrices
L = list(beta = Euc_o,
         host_phy = host_tree_d,
         Country = Country_d_o,
         Lifestyle = Lifestyle_d_o, 
         Family = Family_d_o)

m_perm = lapply(as.list(1:nperm_datasets), function(x) one_per_group(L, df, x))

# MRM on each permutation (in parallel)
doParallel::registerDoParallel(threads)
x = as.list(1:length(m_perm))
f='m$beta ~m$host_phy+ m$Country+ m$Lifestyle+ m$Family'
mrm_res = plyr::llply(x, mrm_each, L=m_perm, f=f, nperm=nperm, .parallel=TRUE)
mrm_res = do.call(rbind, mrm_res)
mrm_res
mrm_res %>%
  filter(variable != 'Int') %>% 
  distinct(coef, coef_pval, R2, pval, rep) %>%
  summary
# significant?
df.dims(10)
mrm_res %>%
  filter(variable != 'Int') %>%
  mutate(pval = pval %>% as.Num) %>%
  group_by(variable) %>%
  summarize(overall_pval = 1 - (sum(coef_pval < 0.05) / length(coef_pval))) %>%
  ungroup()
df.dims()

# formatting results
mrm_res_s = mrm_res %>% 
  filter(variable != 'Int') %>%
  gather(category, value, -variable, -R2, -pval, -F, -rep) %>%
  mutate(category = ifelse(category == 'coef', 'Coef.', 'Adj. P-value'))

mrm_res_s = mrm_res_s %>%
  inner_join(rename_df, c('variable'='old_name')) %>%
  dplyr::select(-variable) %>%
  rename('variable' = new_name)           

# plotting
X = data.frame(yint=c(0.05, NA), category=c('Adj. P-value', 'Coef.'))
p = ggplot(mrm_res_s, aes(variable, value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept=yint), linetype='dashed', alpha=0.3, data=X) +
  labs(x='Variable', y='Intra-species variance') +
  facet_grid(category ~ ., scales='free_y') +
  theme_bw()

p.dims(5,3)
plot(p)
write.table(mrm_res, file="mrm_res_Euc")

# creating subsample permutations of the distance matrices
L = list(beta = Uni_uw_o,
         host_phy = host_tree_d,
         Country = Country_d_o,
         Lifestyle = Lifestyle_d_o, 
         Family = Family_d_o)

m_perm = lapply(as.list(1:nperm_datasets), function(x) one_per_group(L, df, x))

# MRM on each permutation (in parallel)
doParallel::registerDoParallel(threads)
x = as.list(1:length(m_perm))
f='m$beta ~m$host_phy+ m$Country+ m$Lifestyle+ m$Family'
mrm_res = plyr::llply(x, mrm_each, L=m_perm, f=f, nperm=nperm, .parallel=TRUE)
mrm_res = do.call(rbind, mrm_res)
mrm_res
mrm_res %>%
  filter(variable != 'Int') %>% 
  distinct(coef, coef_pval, R2, pval, rep) %>%
  summary
# significant?
df.dims(10)
mrm_res %>%
  filter(variable != 'Int') %>%
  mutate(pval = pval %>% as.Num) %>%
  group_by(variable) %>%
  summarize(overall_pval = 1 - (sum(coef_pval < 0.05) / length(coef_pval))) %>%
  ungroup()
df.dims()

# formatting results
mrm_res_s = mrm_res %>% 
  filter(variable != 'Int') %>%
  gather(category, value, -variable, -R2, -pval, -F, -rep) %>%
  mutate(category = ifelse(category == 'coef', 'Coef.', 'Adj. P-value'))

mrm_res_s = mrm_res_s %>%
  inner_join(rename_df, c('variable'='old_name')) %>%
  dplyr::select(-variable) %>%
  rename('variable' = new_name)           

# plotting
X = data.frame(yint=c(0.05, NA), category=c('Adj. P-value', 'Coef.'))
p = ggplot(mrm_res_s, aes(variable, value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept=yint), linetype='dashed', alpha=0.3, data=X) +
  labs(x='Variable', y='Intra-species variance') +
  facet_grid(category ~ ., scales='free_y') +
  theme_bw()

p.dims(5,3)
plot(p)
write.table(mrm_res, file="mrm_res_UniUW")

# creating subsample permutations of the distance matrices
L = list(beta = Uni_w_o,
         host_phy = host_tree_d,
         Country = Country_d_o,
         Lifestyle = Lifestyle_d_o, 
         Family = Family_d_o)

m_perm = lapply(as.list(1:nperm_datasets), function(x) one_per_group(L, df, x))

# MRM on each permutation (in parallel)
doParallel::registerDoParallel(threads)
x = as.list(1:length(m_perm))
f='m$beta ~m$host_phy+ m$Country+ m$Lifestyle+ m$Family'
mrm_res = plyr::llply(x, mrm_each, L=m_perm, f=f, nperm=nperm, .parallel=TRUE)
mrm_res = do.call(rbind, mrm_res)
mrm_res
mrm_res %>%
  filter(variable != 'Int') %>% 
  distinct(coef, coef_pval, R2, pval, rep) %>%
  summary
# significant?
df.dims(10)
mrm_res %>%
  filter(variable != 'Int') %>%
  mutate(pval = pval %>% as.Num) %>%
  group_by(variable) %>%
  summarize(overall_pval = 1 - (sum(coef_pval < 0.05) / length(coef_pval))) %>%
  ungroup()
df.dims()

# formatting results
mrm_res_s = mrm_res %>% 
  filter(variable != 'Int') %>%
  gather(category, value, -variable, -R2, -pval, -F, -rep) %>%
  mutate(category = ifelse(category == 'coef', 'Coef.', 'Adj. P-value'))

mrm_res_s = mrm_res_s %>%
  inner_join(rename_df, c('variable'='old_name')) %>%
  dplyr::select(-variable) %>%
  rename('variable' = new_name)           
mrm_res_s
# plotting
X = data.frame(yint=c(0.05, NA), category=c('Adj. P-value', 'Coef.'))
p = ggplot(mrm_res_s, aes(variable, value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept=yint), linetype='dashed', alpha=0.3, data=X) +
  labs(x='Variable', y='Intra-species variance') +
  facet_grid(category ~ ., scales='free_y') +
  theme_bw()

p.dims(5,3)
plot(p)
write.table(mrm_res, file="mrm_res_UniW")

