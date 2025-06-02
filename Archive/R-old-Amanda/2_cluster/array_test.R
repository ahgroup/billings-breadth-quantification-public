
#Install/Load packages
#package is installed in my directory
library(base)
library(Racmacs)

#Load the data

args = commandArgs(TRUE)
file = args[1]

g <- readRDS(file)

map_g <- acmap(titer_table = g[,-1], ag_names = as.character(g[,1]), sr_names = colnames(g)[-1])

#Dimension Testing
d <- dimensionTestMap(map = map_g)

saveRDS(object = d, 
        file = paste("dimen_analysis_",
                      substr(file,
                             9,
                             nchar(file)),
                     sep = ""))
