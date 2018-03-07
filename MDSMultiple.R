#!/usr/bin/env Rscript
source ("tool.R")
require(lattice)
library(scatterplot3d)


MDS3DGlobal = function(dglobal, prout, disttype){
  
  
  if(disttype == "corr"){
    MC = cor(t(scale(dglobal)))
    distG = as.matrix(abs(1-MC))
  }else{
    distG = dist(scale(dglobal), method = disttype)
  }
  dim(distG)
    fit1 = cmdscale(as.dist(distG),eig=TRUE, k=3)
 
    coords = cbind(fit1$points[,1], fit1$points[,2])
    coords = cbind(coords, fit1$points[,3])
 
    gifGeneration(paste(prout, "MDS3D_", disttype, sep = ""), coords)
}



MDSMulti = function (d1D2D, pout, disttype){
  
  if (type_dist == "corr"){
    MC1 = cor(t(scale(d1D)))
    dist1 = as.matrix(abs(1-MC1))
    
    MC3 = cor(t(scale(d3D)))
    dist3 = as.matrix(abs(1-MC3))
  }else{
    dist1 = dist(scale(d1D), method = disttype)
    dist3 = dist(scale(d3D), method = disttype)
  }

  fit1 = cmdscale(as.dist(dist1),eig=TRUE, k=2)
  print("1")
  fit3 = cmdscale(as.dist(dist3),eig=TRUE, k=1)
  print("3")
  
  coords = cbind(fit1$points[,1], fit1$points[,2])
  coords  = cbind(coords, fit3$points[,1])
  
  #print (coords)
  
  gifGeneration(paste(prout, "MDSCombined_", disttype, sep = ""), coords)
  
}
