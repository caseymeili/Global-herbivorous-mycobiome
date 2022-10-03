DM <-read.table("/Users/nohayoussef/Desktop/GradStudents/Casey-Adrienne/FecesPaper/GeneraTable.txt", header=TRUE)

library(BiodiversityR)
Rank_abundance <-rankabundance(DM)
rankabunplot(Rank_abundance, scale='abundance')
