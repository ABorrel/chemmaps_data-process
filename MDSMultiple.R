#!/usr/bin/env Rscript
source ("tool.R")
require(lattice)
library(scatterplot3d)




MDSMulti = function (d1D, d2D, d3D, pout, type_dist){
  
  if (type_dist == "corr"){
    MC1 = cor(t(scale(d1D)))
    dist1 = as.matrix(abs(1-MC1))
    
    MC2 = cor(t(scale(d2D)))
    dist2 = as.matrix(abs(1-MC2))
    
    MC3 = cor(t(scale(d3D)))
    dist3 = as.matrix(abs(1-MC3))
  }

  if (type_dist == "euc"){
    dist1 = dist(scale(d1D), method = "euclidean")
    dist2 = dist(scale(d2D), method = "euclidean")
    dist3 = dist(scale(d3D), method = "euclidean")
  }

  fit1 = cmdscale(as.dist(dist1),eig=TRUE, k=1)
  print("1")
  fit2 = cmdscale(as.dist(dist2),eig=TRUE, k=1)
  print("2")
  fit3 = cmdscale(as.dist(dist3),eig=TRUE, k=1)
  print("3")
  
  coords = cbind(fit1$points[,1], fit2$points[,1])
  coords  = cbind(coords, fit3$points[,1])
  
  print (coords)
  
  gifGeneration(pout, coords)
  
}


################
#     MAIN     #
################

args <- commandArgs(TRUE)
p1D = args[1]
p2D = args[2]
p3D = args[3]
prout = args[4]


p1D = "/home/aborrel/ChemMap/results/Desc1D.csv"
p2D = "/home/aborrel/ChemMap/results/Desc2D.csv"
p3D = "/home/aborrel/ChemMap/results/Desc3D.csv"
prout = "/home/aborrel/ChemMap/results/analysis/MCSs/"


valcor = 0

d1D = openData(p1D, valcor, prout, c(1,2))
d2D = openData(p2D, valcor, prout, c(1,2))


d1D_data = d1D[[1]]
rownames(d1D_data) = d1D_data[,1]
d1D_data = d1D_data[,-1]
d1D_data = d1D_data[,-1] # SMILES

d2D_data = d2D[[1]]
rownames(d2D_data) = d2D_data[,1]
d2D_data = d2D_data[,-1]
d2D_data = d2D_data[,-1] # SMILES


#3D need to remove empty cols and rows
d3 = read.csv(p3D, sep = "\t", header = TRUE)
# remove compound not computed
d3 = delete.na(d3,(dim(d3)[1] - 3))
# remove empty descriptor
d3 = delete.na(t(d3),(dim(d3)[2] - 20))
d3clean = delete.na(t(d3),(dim(d3)[1] - 3))

# rewrite descriptors3D
p3Dclean = paste(strsplit(p3D, "[.]")[[1]][1], "_clean.csv", sep = "")
write.table(d3clean, file = p3Dclean, row.names = FALSE, sep = "\t")
d3D = openData(p3Dclean, valcor, prout, c(1,2))
d3D_data = d3D[[1]]
rownames(d3D_data) = d3D_data[,1]
d3D_data = d3D_data[,-1]
d3D_data = d3D_data[,-1] # SMILES

d1D_data = delnohomogeniousdistribution(d1D_data, 75)
d2D_data = delnohomogeniousdistribution(d2D_data, 75)
d3D_data = delnohomogeniousdistribution(d3D_data, 75)



#fusion dataset
vcompound = rownames(d1D_data)
vcompound = intersect(vcompound,rownames(d2D_data))
vcompound = intersect(vcompound,rownames(d3D_data))
dglobal = cbind(d2D_data[vcompound,], d1D_data[vcompound,])
dglobal = cbind(dglobal, d3D_data[vcompound,])


MDSMulti(d1D_data[vcompound,], d2D_data[vcompound,], d3D_data[vcompound,], paste(prout, "euc", sep = ""), "euc")
MDSMulti(d1D_data[vcompound,], d2D_data[vcompound,], d3D_data[vcompound,], paste(prout, "corr", sep = ""), "corr")

#print(dim(dglobal))
#print (colnames(dglobal))
#print (dglobal[1:3,])