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
pmaptomap = "/home/borrela2/ChemMaps/data_analysis/PFAS/forDB/db.csv" 
prout = "/home/borrela2/ChemMaps/data_analysis/DSSTox/"


dcoord1D2D = read.csv(pcoord1D2D, sep = ",", header = TRUE)
dcoord1D2D = dcoord1D2D[,c(1,2, 3)]
rownames(dcoord1D2D) = dcoord1D2D[,1]
dcoord1D2D = dcoord1D2D[,-1]

dcoord3D = read.csv(pcoord3D, sep = ",", header = TRUE)
dcoord3D = dcoord3D[,c(1,2, 3)]
rownames(dcoord3D) = dcoord3D[,1]
dcoord3D = dcoord3D[,-1]

dcoord3D = dcoord3D[rownames(dcoord1D2D),]

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(100)
png(paste(prout, "global_xy_DSSTOX.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM2~DIM1,data=dcoord1D2D,xlab="DIM1",
           ylab="DIM2",colramp=rf, xbins = 500, cex.labels = 0.6, cex.title=0.6)
dev.off()


dcoordCut = dcoord1D2D[which(dcoord1D2D[,1] <= 30 & dcoord1D2D[,2] <= 30),]
png(paste(prout, "cut_xy_DSSTOX.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM2~DIM1,data=dcoordCut,xlab="DIM1",
           ylab="DIM2",colramp=rf, xbins = 500, cex.labels = 0.6, cex.title=0.6)
dev.off()


dcoord_xzD = data.frame(DIM1 = dcoord1D2D$DIM1, DIM3 = dcoord3D$DIM3.1)
png(paste(prout, "global_xz_DSSTOX.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM3~DIM1,data=dcoord_xzD,xlab="DIM1",
           ylab="DIM3",colramp=rf, xbins = 500, cex.labels = 0.6, cex.title=0.6)
dev.off()

dcoordcut = dcoord[which(dcoord[,1] <= 30 & dcoord[,2] <= 30),]
png(paste(prout, "cut_xz_DSSTOX.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM3~DIM1,data=dcoordcut,xlab="DIM1",
           ylab="DIM3",colramp=rf, xbins = 500, cex.labels = 0.6, cex.title=0.6)
dev.off()



# map chemical on the map
dmap = read.csv(pmaptomap, sep = "\t", header = TRUE)
