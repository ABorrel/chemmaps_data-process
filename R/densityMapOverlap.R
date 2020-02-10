#!/usr/bin/env Rscript
library(hexbin)
library(RColorBrewer)



################
#     MAIN     #
################

args <- commandArgs(TRUE)
pcoordMain1D2D = args[1]
pcoordMain3D = args[2]
plistChem = args[3]
name_map = args[4]
prout = args[5]


# linux 


# window
pcoordMain1D2D = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/DSSTox/map_0.9-90/coord1D2D_light.csv"
pcoordMain3D = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/DSSTox/map_0.9-90/coord3D_light.csv"
plistChem = "/Users/Aborrel/research/NIEHS/ChemMaps/data/NPAHs.csv"
name_map = "NPAH_test"
prout = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/DSSTox/density_map/"




dcoord1D2D = read.csv(pcoordMain1D2D, sep = ",", header = TRUE)
dcoord1D2D = dcoord1D2D[,c(1,2, 3)]
rownames(dcoord1D2D) = dcoord1D2D[,1]
dcoord1D2D = dcoord1D2D[,-1]

dcoord3D = read.csv(pcoordMain3D, sep = ",", header = TRUE)
dcoord3D = dcoord3D[,c(1,2, 3)]
rownames(dcoord3D) = dcoord3D[,1]
dcoord3D = dcoord3D[,-1]

dcoord3D = dcoord3D[rownames(dcoord1D2D),]


# remove outlier values =limit to a box 80*80
dcoord_xy_cut =  dcoord1D2D[which(dcoord1D2D[,1] <= 80 & dcoord1D2D[,2] <= 80),]

# for cut data
dcoord_xzD_cut = dcoord_xzD[which(dcoord_xzD[,1] <= 80 & dcoord_xzD[,2] <= 30 & dcoord_xzD[,2] >= -50),]


# color
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(100)










png(paste(prout, "global_xy_DSSTOX.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM2~DIM1,data=dcoord1D2D,xlab="DIM1",
           ylab="DIM2",colramp=rf, xbins = 500, cex.labels = 0.6, cex.title=0.6)
dev.off()


dcoordCut = dcoord1D2D[which(dcoord1D2D[,1] <= 80 & dcoord1D2D[,2] <= 80),]
png(paste(prout, "DSSTOX_xycut.png", sep = ""), res = 300, width = 2000, height = 2000)
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


############################
# PORJECTION ON THE DSSTOX #
############################
#pmaptomap = "/home/borrela2/ChemMaps/data/NPAHs.csv" #==> NPAHs
#pmaptomap = "/home/borrela2/ChemMaps/data/PAH-14Jan2020.csv"
pmaptomap = "/home/borrela2/ChemMaps/data_analysis/TOX21/forDB/db.csv"# ==> Tox21
prout = "/home/borrela2/ChemMaps/data_analysis/density_map/"
name_map = "Tox21"

# map chemical on the map
dmap = read.csv(pmaptomap, sep = "\t", header = TRUE)
dmap = na.omit(dmap)
rownames(dmap) = dmap[,1]
dmap = dmap[,-1]


dmap_plot = dmap
#dmap_cancer = dmap[dmap$cancer == "y", ]
#dmap_cancer = dmap[dmap$cancer == "n", ]


dcoorMap = dcoordCut[as.vector(dmap_plot$inchikey),]
dcoorMap = rbind(dcoorMap, c(min(dcoordCut$DIM1), min(dcoordCut$DIM2)))
dcoorMap = rbind(dcoorMap, c(max(dcoordCut$DIM1), max(dcoordCut$DIM2)))

png(paste(prout, paste(name_map, "VSDSSTOX_D1D2.png"), sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM2~DIM1,data=dcoorMap,xlab="DIM1",
           ylab="DIM2",colramp=rf, xbins = 50, cex.labels = 0.6, cex.title=0.6 ) +
  geom_segment(data=pts, aes(x, y, xend=xend, yend=yend))
dev.off()





