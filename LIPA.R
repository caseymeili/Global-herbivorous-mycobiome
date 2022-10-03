# Local Indicator of Phylogenetic Association (LIPA) 

library(phylosignal)
library(ape)
library(phylobase)
library(phytools)

## read the percentage abundance file. This will have samples as rows, genera as columns 
## and percent abundance as cell values
PercAbund <-read.table("~/Desktop/R/GeneraTable.txt", header=TRUE)

## read host tree
Host_Sp <-ape::read.tree(file="~/Desktop/R/Host.nwk")

## read the genera tree
Taxa <-ape::read.tree(file="~/Desktop/R/Genera_rooted.nwk")

## create a phylo4d object
p4dAll <-phylo4d(Host_sp, data.frame(PecrAbund))

## calculate the Lipa index for each sample-genus pair
LipaAll <-lipaMoran(p4dAll)

## read the output from Lipa and copy it to excel
LipaALL$lipa
LipaAll$p.value

## calculate other indices of phylosignal
PS <-phylosignal::phyloSignal(p4dAll, methods="all")
PS$stat
PS$pvalue

## create two trees facing one another with association in between
## the association file is two columns with animal in first column and genus in second
## Decisions on which associations to include depend on the lipa value. I used >1 as strong
## read association file
Assoc <-read.table("~/Desktop/R/StrongAssoc.txt", header=TRUE)
cophylo(Host_Sp, Taxa, assoc=Assoc)
