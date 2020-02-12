#!/usr/bin/env Rscript
library(hexbin)
library(RColorBrewer)
library(ggplot2)


# redine function for hexbio
library (reshape2)


bindata <- data.frame(x=rnorm(100), y=rnorm(100))
fac_probs <- dnorm(seq(-3, 3, length.out=26))
fac_probs <- fac_probs/sum(fac_probs)
bindata$factor <- sample(letters, 100, replace=TRUE, prob=fac_probs)


bindata$factor <- as.factor (bindata$factor)
h <- hexbin (bindata, xbins = 5, IDs = TRUE, 
             xbnds = range (bindata$x), 
             ybnds = range (bindata$y))

counts <- hexTapply (h, bindata$factor, table)
counts <- t (simplify2array (counts))
counts <- melt (counts)
colnames (counts)  <- c ("ID", "factor", "counts")


hexdf <- data.frame (hcell2xy (h),  ID = h@cell)
hexdf <- merge (counts, hexdf)

hexdf$counts [hexdf$counts == 0] <- NA

ggplot(hexdf, aes(x=x, y=y, fill = counts)) +
  geom_hex(stat="identity") +
  facet_wrap(~factor) +
  coord_equal () +
  scale_fill_continuous (low = "grey80", high = "#000040", na.value = "#00000000")


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
pcoordSet1_1D2D = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/TOX21/map_0.9-90/coord1D2D.csv"
pcoordSet1_3D = "/Users/Aborrel/research/NIEHS/ChemMaps/data_analysis/TOX21/map_0.9-90/coord3D.csv"

pcoordSet2_1D2D = "/Users/Aborrel/Desktop/Chemmaps_forWarren/skin_chem_dim1-2.csv"
pcoordSet2_3D = "/Users/Aborrel/Desktop/Chemmaps_forWarren/skin_chem_dim3.csv"
name_map = "Skin"
prout = "/Users/Aborrel/Desktop/Chemmaps_forWarren/"



# load set 1
dcoord1D2D_set1 = read.csv(pcoordSet1_1D2D, sep = ",", header = TRUE)
dcoord1D2D_set1 = dcoord1D2D_set1[,c(1,2, 3)]
rownames(dcoord1D2D_set1) = dcoord1D2D_set1[,1]
dcoord1D2D_set1 = dcoord1D2D_set1[,-1]

dcoord3D_set1 = read.csv(pcoordSet1_3D, sep = ",", header = TRUE)
dcoord3D_set1 = dcoord3D_set1[,c(1,2, 3)]
rownames(dcoord3D_set1) = dcoord3D_set1[,1]
dcoord3D_set1 = dcoord3D_set1[,-1]

lchem = intersect(rownames(dcoord1D2D_set1), rownames(dcoord3D_set1))
dcoord3D_set1 = dcoord3D_set1[lchem,]
dcoord1D2D_set1 = dcoord1D2D_set1[lchem,]


# load set2
dcoord1D2D_set2 = read.csv(pcoordSet2_1D2D, sep = ",", header = TRUE)
dcoord1D2D_set2 = dcoord1D2D_set2[,c(1,2, 3)]
rownames(dcoord1D2D_set2) = dcoord1D2D_set2[,1]
dcoord1D2D_set2 = dcoord1D2D_set2[,-1]

dcoord3D_set2 = read.csv(pcoordSet2_3D, sep = ",", header = TRUE)
dcoord3D_set2 = dcoord3D_set2[,c(1,2, 3)]
rownames(dcoord3D_set2) = dcoord3D_set2[,1]
dcoord3D_set2 = dcoord3D_set2[,-1]

lchem = intersect(rownames(dcoord1D2D_set2), rownames(dcoord3D_set2))
dcoord3D_set2 = dcoord3D_set2[lchem,]
dcoord1D2D_set2 = dcoord1D2D_set2[lchem,]


# color1
rf1 <- colorRampPalette(rev(brewer.pal(9,'Blues')))
r1 <- rf1(100)

rf2 <- colorRampPalette(rev(brewer.pal(9,'Oranges')))
r2 <- rf2(100)


# plot => DIM1 VS DIM2 #
########################
# calibrate the map
maxdim1 = max(rbind(dcoord1D2D_set1$DIM1, dcoord1D2D_set1$DIM2))
maxdim2 = max(rbind(dcoord1D2D_set1$DIM2, dcoord1D2D_set1$DIM2))

mindim1 = min(rbind(dcoord1D2D_set1$DIM1, dcoord1D2D_set1$DIM2))
mindim2 = min(rbind(dcoord1D2D_set1$DIM2, dcoord1D2D_set1$DIM2))

dcoorMap_set1 = rbind(dcoord1D2D_set1, c(mindim1, mindim2))
dcoorMap_set1 = rbind(dcoorMap_set1, c(maxdim1, maxdim2))

dcoorMap_set2 = rbind(dcoord1D2D_set2, c(mindim1, mindim2))
dcoorMap_set2 = rbind(dcoorMap_set2, c(maxdim1, maxdim2))

# set1
png(paste(prout, paste(name_map, "_set1.png", sep = ""), sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM2~DIM1,data=dcoorMap_set1,xlab="DIM1",
           ylab="DIM2",colramp=rf1, xbins = 50, cex.labels = 0.6, cex.title=0.6, aspect = 2 )
dev.off()


# set2
png(paste(prout, paste(name_map, "_set2.png", sep = ""), sep = ""), res = 300, width = 2000, height = 2000)
hexbinplot(DIM2~DIM1,data=dcoorMap_set2,xlab="DIM1",
           ylab="DIM2",colramp=rf2, xbins = 50, cex.labels = 0.6, cex.title=0.6, aspect = 2)
dev.off()
