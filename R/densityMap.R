#!/usr/bin/env Rscript
library(hexbin)
library(RColorBrewer)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pcoord1D2D = args[1]
pcoord3D = args[2]
prout = args[3]


pcoord1D2D = "/home/borrela2/ChemMaps/data_analysis/DSSTox/map_0.9-90/coord1D2D.csv"
pcoord3D = "/home/borrela2/ChemMaps/data_analysis/DSSTox/map_0.9-90/coord3D.csv"
prout = "/home/borrela2/ChemMaps/data_analysis/DSSTox/"


dcoord1D2D = read.csv(pcoord1D2D, sep = ",", header = TRUE)
dcoord1D2D = dcoord1D2D[,c(1,2, 3)]
rownames(dcoord1D2D) = dcoord1D2D[,1]
dcoord1D2D = dcoord1D2D[,-1]

dcoord3D = read.csv(pcoord1D2D, sep = ",", header = TRUE)
dcoord1D2D = dcoord1D2D[,c(1,2, 3)]
rownames(dcoord1D2D) = dcoord1D2D[,1]
dcoord1D2D = dcoord1D2D[,-1]


rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(100)
hexbinplot(DIM2~DIM1,data=dcoord1D2D,xlab="DIM1",
           ylab="DIM2",colramp=rf, xbins = 500)


dcoordCut = dcoord1D2D[which(dcoord1D2D[,1] <= 20 & dcoord1D2D[,2] <= 20),]

hexbinplot(DIM2~DIM1,data=dcoordCut,xlab="DIM1",
           ylab="DIM2",colramp=rf, xbins = 200)
