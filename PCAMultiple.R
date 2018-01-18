#!/usr/bin/env Rscript
source ("tool.R")
source("cardMatrix.R")
library(fastICA)
library (vrmlgen)

generatePCAcoords = function(din){
  
  dinScale = scale(din)
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp

  return(list(data_plot, var_cap))
  
  
}


generateIPAcoords = function(din, path_result){
  
  a = fastICA(din, 3, alg.typ = "parallel", fun = "logcosh", alpha = 1, method = "R", row.norm = FALSE, maxit = 200, tol = 0.0001, verbose = TRUE)
  print (a)
  
  gifGeneration(paste(path_result, "ICA3D", sep = ""), a$S)
}


PCAcombined3plans = function(d1D, d2D, d3D, pfilout){
  
  lcoord1D = generatePCAcoords(d1D)
  lcoord2D = generatePCAcoords(d2D)
  lcoord3D = generatePCAcoords(d3D)
  
  coordSpace = cbind(lcoord1D[[1]][,1], lcoord2D[[1]][,1])
  coordSpace = cbind(coordSpace, lcoord3D[[1]][,1])
  
  print(paste(lcoord1D[[2]], lcoord2D[[2]], lcoord3D[[2]], sep = "%  "))
  
  gifGeneration(paste(pfilout, "PCA3D", sep = ""), coordSpace)
  
  
}



PCAcombined2plans = function(d1D, d2D, pfilout){
  
  lcoord1D = generatePCAcoords(d1D)
  lcoord3D = generatePCAcoords(d2D)
  
  coordSpace = cbind(lcoord1D[[1]][,1], lcoord1D[[1]][,2])
  coordSpace = cbind(coordSpace, lcoord3D[[1]][,1])
  
  print(paste(lcoord1D[[2]][1], lcoord1D[[2]][2], lcoord3D[[2]][1], sep = "%  "))
  
  gifGeneration(paste(pfilout, "PCA3D2plan", sep = ""), coordSpace)
  
  
}




PCA3D = function(din, path_result){
  
  dinScale = scale(din)
  
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp
  
  #col.desc = colorDesc(colnames(din))
  
  col.desc = "black"
  
  gifGeneration(paste(path_result, "PCA3D", sep = ""), data_plot)
  
}



PCAplot = function (din, path_result){
  
  dinScale = scale(din)
  
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp
  
  #col.desc = colorDesc(colnames(din))
  
  col.desc = "black"
  
  print(col.desc)
  
  png (paste (path_result, "_text.png", sep = ""), 1700, 1500)
  factor = factorACP (data_plot, cp)
  
  color_arrow = col.desc[rownames(cp)]
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  text (data_plot[,1],data_plot[,2], label = rownames (din), cex = 1.2)
  abline(h=0,v=0)
  warnings ()
  dev.off()
  
  colpoint <- colorRampPalette(c("white", "black"))
  
  png (paste (path_result, "_color.png", sep = ""), 1700, 1500)
  factor = factorACP (data_plot, cp)
  color_arrow =col.desc
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], pch=20, col = colpoint(dim(data_plot)[1]), xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4)
  abline(h=0,v=0)
  warnings ()
  dev.off()
  
  
  png (paste (path_result, "_descriptor.png", sep = ""), 1700, 1500)
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), pch=20, cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
  #points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)
  abline(h=0,v=0)
  arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = color_arrow, lwd = 4 )
  text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = color_arrow, cex = 2.5)
  dev.off()
  
  
  svg (file = paste (path_result, "_descriptor.svg", sep = ""), 25, 25)
  par(mar=c(8,8,8,8))
  plot(data_plot[,1],data_plot[,2], main = "", xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 6, cex.main = 4, cex.axis = 1.75, cex = 6, type = "n")
  #points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)
  abline(h=0,v=0)
  arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = color_arrow, lwd = 3 )
  text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = color_arrow, cex = 3.5)
  dev.off()
  
}



model3D = function(d1D, d2D, dcol, pfilout){
  
  lcoord1D = generatePCAcoords(d1D)
  lcoord3D = generatePCAcoords(d2D)
  
  coordSpace = cbind(lcoord1D[[1]][,1], lcoord1D[[1]][,2])
  coordSpace = cbind(coordSpace, lcoord3D[[1]][,1])
  
  coordApproved = coordSpace[which(dcol[,1] == "approved"),]
  coordInDev = coordSpace[which(dcol[,2] == "grey"),]
  coordwhidraW = coordSpace[which(dcol[,2] == "red"),]
  
  write.csv(coordApproved, file = paste(pfilout, "approvedcorrd3D.csv"))
  write.csv(coordInDev, file = paste(pfilout, "indevcorrd3D.csv"))
  write.csv(coordwhidraW, file = paste(pfilout, "withDrawcorrd3D.csv"))
  #print (coordSpace[seq(1,10),])
  #coordSpace = 100*coordSpace
  #print (coordSpace[seq(1,10),])
  
  vrml.open(file = "3Dmodel.wrl", scale = 1,5)
  points3d(coordApproved, col = "green", scale = 2.6, pointstyle = "b" )
  points3d(coordInDev, col = "white", scale = 1, transparency = 0.5,  pointstyle = "b")
  points3d(coordwhidraW, col = "red", scale = 2.2, transparency = 0.2,  pointstyle = "b")
  #text3d(coordApproved[1,1], coordApproved[1,2], coordApproved[1,3] + 8, text = rownames(coordApproved)[1], col = "green", scale = 5)
  lines3d(x = c(0,10), y = c(0,0), z = c(0,0), lwd = 3, col = "blue")
  lines3d(x = c(0,0), y = c(0,10), z = c(0,0), lwd = 3, col = "blue")
  lines3d(x = c(0,0), y = c(0,0), z = c(0,10), lwd = 3, col = "blue")
  vrml.close()
  #print(paste(lcoord1D[[2]][1], lcoord1D[[2]][2], lcoord3D[[2]][1], sep = "%  "))
  
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
pprop = "/home/aborrel/ChemMap/structures.csv"
prout = "/home/aborrel/ChemMap/results/analysis/PCAs/"

# manually define
outlier = c("DB00793", "DB00516", "DB06690", "DB09157", "DB03627", "DB01751", "DB03853" ,"DB04711")

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


# prop for color
dprop = read.csv(pprop, sep = "\t", stringsAsFactors = F, header = TRUE)
dcol = as.matrix(dprop[,c("DRUGBANK_ID","DRUG_GROUPS")])
dcol = cbind(dcol, rep("grey", dim(dcol)[1]))
dcol[which(dcol[,2] == "approved"),3] = "green"
dcol[which(dcol[,2] == "withdrawn"),3] = "red"
dcol[which(dcol[,2] == "illicit"),3] = "orange"
rownames(dcol) = dcol[,1]
dcol = dcol[,-1]

# selection using distribution
#print (dim(d1D_data))
#print (dim(d2D_data))
#print (dim(d3D_data))


d1D_data = delnohomogeniousdistribution(d1D_data, 80)
d2D_data = delnohomogeniousdistribution(d2D_data, 80)
d3D_data = delnohomogeniousdistribution(d3D_data, 80)

#print (dim(d1D_data))
#print (dim(d2D_data))
#print (dim(d3D_data))

#fusion dataset
vcompound = rownames(d1D_data)
vcompound = intersect(vcompound,rownames(d2D_data))
vcompound = intersect(vcompound,rownames(d3D_data))
vcompound = vcompound[!vcompound %in% outlier]
dglobal = cbind(d2D_data[vcompound,], d1D_data[vcompound,])
dglobal = cbind(dglobal, d3D_data[vcompound,])
dcol = dcol[vcompound,]

#cardMatrixCor(cor(cbind(d1D_data[vcompound,], d2D_data[vcompound,])), paste(prout, "matrixCor1D2D", sep = ""), 6)
#cardMatrixCor(cor(cbind(d2D_data[vcompound,], d3D_data[vcompound,])), paste(prout, "matrixCor2D3D", sep = ""), 6)
#cardMatrixCor(cor(cbind(d1D_data[vcompound,], d3D_data[vcompound,])), paste(prout, "matrixCor1D3D", sep = ""), 6)

cardMatrix(scale(dglobal), paste(prout, "data.png", sep = ""), 6)

#plot histogram
#histDataOne(data1 = d1D_data[vcompound,], paste(prout, "homodishist1D.pdf"))
#histDataOne(data1 = d2D_data[vcompound,], paste(prout, "homodishist2D.pdf"))
#histDataOne(data1 = d3D_data[vcompound,], paste(prout, "homodishist3D.pdf"))

#generateIPAcoords(dglobal, prout)
#print (colnames(dglobal))
#print (dglobal[1:3,])

#PCAplot(dglobal, paste(prout, "global", sep = ""))
#PCA3D(dglobal, paste(prout, "global", sep = ""))
#PCAplot(d1D_data[vcompound,], paste(prout, "1DDesc", sep = ""))
#PCAplot(d2D_data[vcompound,], paste(prout, "2DDesc", sep = ""))
#PCAplot(d3D_data[vcompound,], paste(prout, "3DDesc", sep = ""))


#PCA multiple
#PCAcombined3plans(d1D_data[vcompound,], d2D_data[vcompound,], d3D_data[vcompound,], prout)
#PCAcombined2plans(cbind(d1D_data[vcompound,], d2D_data[vcompound,]), d3D_data[vcompound,], prout)

#model3D(cbind(d1D_data[vcompound,], d2D_data[vcompound,]), d3D_data[vcompound,], dcol, prout)