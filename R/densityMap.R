#!/usr/bin/env Rscript
library(hexbin)
library(RColorBrewer)
library(rgl)

################
#     MAIN     #
################

args <- commandArgs(TRUE)
pcoord1D2D = args[1]
pcoord3D = args[2]
name_map = args[3]
prout = args[4]


# linux 
pcoord1D2D = "/home/borrela2/ChemMaps/data_analysis/DSSTox/map_0.9-90/coord1D2D.csv"
pcoord3D = "/home/borrela2/ChemMaps/data_analysis/DSSTox/map_0.9-90/coord3D.csv"
prout = "/home/borrela2/ChemMaps/data_analysis/NPAHs/"


### window +> for test
# DSSTOX
pcoord1D2D = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/DSSTox/map_0.9-90/coord1D2D_light.csv"
pcoord3D = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/DSSTox/map_0.9-90/coord3D_light.csv"
name_map = "DSSTox_test"
prout = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/DSSTox/density_map/"

# Tox21
pcoord1D2D = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/TOX21/map_0.9-90/coord1D2D.csv"
pcoord3D = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/TOX21/map_0.9-90/coord3D.csv"
name_map = "Tox21"
prout = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/TOX21/density_map/"



dcoord1D2D = read.csv(pcoord1D2D, sep = ",", header = TRUE)
dcoord1D2D = dcoord1D2D[,c(1,2, 3)]
rownames(dcoord1D2D) = dcoord1D2D[,1]
dcoord1D2D = dcoord1D2D[,-1]

dcoord3D = read.csv(pcoord3D, sep = ",", header = TRUE)
dcoord3D = dcoord3D[,c(1,2, 3)]
rownames(dcoord3D) = dcoord3D[,1]
dcoord3D = dcoord3D[,-1]

lchem = intersect(rownames(dcoord1D2D), rownames(dcoord3D))
dcoord3D = dcoord3D[lchem,]
dcoord1D2D = dcoord1D2D[lchem,]

# remove outlier values =limit to a box 80*80
dcoord1D2D_cut =  dcoord1D2D[which(dcoord1D2D[,1] <= 80 & dcoord1D2D[,2] <= 80),]

# plot 3D +> html
#plot_ly(x=temp, y=pressure, z=dtime, type="scatter3d", mode="markers", color=temp)
plot3d(x=dcoord1D2D$DIM1, y=dcoord1D2D$DIM2, z=dcoord3D$DIM3.1)
writeWebGL(filename=paste(prout, name_map, "_3Dplot.html") ,width=2000, height=2000)


## make a projection on 2D ##
# define color
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(100)

# make the projection
# for all data => x, y
png(paste(prout, name_map, "_xy_all.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM2~DIM1,data=dcoord1D2D,xlab="DIM1",
           ylab="DIM2",colramp=rf, xbins = 500, cex.labels = 0.6, cex.title=0.6)
dev.off()

# for cut data
png(paste(prout, name_map, "_xy_cut.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM2~DIM1,data=dcoord1D2D_cut, xlab="DIM1",
           ylab="DIM2",colramp=rf, xbins = 500, cex.labels = 0.6, cex.title=0.6)
dev.off()


# for all data => x, z
dcoord_xzD = data.frame(DIM1 = dcoord1D2D$DIM1, DIM3 = dcoord3D$DIM3.1)
png(paste(prout, name_map, "_xz_all.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM3~DIM1,data=dcoord_xzD,xlab="DIM1",
           ylab="DIM3",colramp=rf, xbins = 500, cex.labels = 0.6, cex.title=0.6)
dev.off()

# for cut data
dcoord_xzD_cut = dcoord_xzD[which(dcoord_xzD[,1] <= 80 & dcoord_xzD[,2] <= 30 & dcoord_xzD[,2] >= -50),]

png(paste(prout, name_map, "_xz_cut.png", sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM3~DIM1,data=dcoord_xzD_cut,xlab="DIM1",
           ylab="DIM3",colramp=rf, xbins = 500, cex.labels = 0.6, cex.title=0.6)
dev.off()

