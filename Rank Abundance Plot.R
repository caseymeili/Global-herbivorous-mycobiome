DM <-read.table("/Users/caseymeili/Desktop/R/GeneraTable.txt", header=TRUE)

library(BiodiversityR)
Rank_abundance <-rankabundance(DM)
rankabunplot(Rank_abundance, scale='abundance')
