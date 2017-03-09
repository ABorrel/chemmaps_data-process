#!/usr/bin/env Rscript
source ("tool.R")
require(visNetwork, quietly = TRUE)



networkhtml = function(din, cutoff){
  
  dist1 = as.matrix(dist(scale(din), method = "euclidean", upper = TRUE))
  #print (dist1)
  print(dim(dist1))
  #fit1 = cmdscale(as.dist(dist1),eig=TRUE, k=2)
  
  #print(fit1)
  nodes = data.frame(id=1:dim(din)[1])
  nb_node = dim(dist1)[1]
  
  edges = NULL
  i = 1
  while (i < nb_node){
    j = i + 1
    while (j <= nb_node){
      dist_temp = dist1[i,j]
      if(dist_temp <= cutoff){
        edges = rbind(edges, c(i,j))
      }
      j = j + 1
    }
    i = i + 1
  }
  
  print (edges)
  colnames(edges) = c("from", "to")
  edges = as.data.frame(edges)
  
  network = visNetwork(nodes, edges, width = "100%")
  visSave(network, file = "network.html")
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


valcor = 0.9

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



#fusion dataset
vcompound = rownames(d1D_data)
vcompound = intersect(vcompound,rownames(d2D_data))
vcompound = intersect(vcompound,rownames(d3D_data))
dglobal = cbind(d2D_data[vcompound[1:100],], d1D_data[vcompound[1:100],])
dglobal = cbind(dglobal, d3D_data[vcompound[1:100],])

networkhtml (dglobal, 15)
