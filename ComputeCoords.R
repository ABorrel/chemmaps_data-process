#!/usr/bin/env Rscript
source ("tool.R")
source("cardMatrix.R")
source("MDSMultiple.R")
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


generateICAcoords = function(din, path_result){
  
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


generateCoordCombinedPCA = function(d1D2D, d3D, prout){

  lcoord1D2D = generatePCAcoords(d1D2D)
  lcoord3D = generatePCAcoords(d3D)
  
  coordSpace = cbind(lcoord1D2D[[1]][,1], lcoord1D2D[[1]][,2])
  coordSpace = cbind(coordSpace, lcoord3D[[1]][,1])
  
  
  colnames(coordSpace) = c("1D", "2D", "3D")
  rownames(coordSpace) = rownames(d1D2D)
  
  write.csv(coordSpace, file = paste(prout, "coordPCAcombined.csv", sep = ""), row.names = TRUE, col.names = "TRUE")
  
  variancePlane = c(lcoord1D2D[[2]][1], lcoord1D2D[[2]][2], lcoord3D[[2]][1])
  names(variancePlane) = c("1D", "2D", "3D")
  write.csv(variancePlane, file = paste(prout, "variancePlanPCAcombined.csv", sep = ""), row.names = TRUE, col.names = "TRUE")
  
  vrml.open(file = "SpacePCAcombined.wrl", scale = 1,5)
  points3d(coordSpace, col = "white", scale = 1, transparency = 0.5,  pointstyle = "b")
  lines3d(x = c(0,10), y = c(0,0), z = c(0,0), lwd = 3, col = "blue")
  lines3d(x = c(0,0), y = c(0,10), z = c(0,0), lwd = 3, col = "blue")
  lines3d(x = c(0,0), y = c(0,0), z = c(0,10), lwd = 3, col = "blue")
  vrml.close()
  print(paste(lcoord1D2D[[2]][1], lcoord1D2D[[2]][2], lcoord3D[[2]][1], sep = "%  "))
  
    
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
p1D2D = args[2]
p3D = args[3]
prout = args[4]
valcor = as.double(args[5])
maxquantile = as.integer(args[6])


p1D2D = "c://Users/Aborrel/chemmaps/drugbankDesc/1D2D.csv"
p3D = "c://Users/Aborrel/chemmaps/drugbankDesc/3D.csv"
#pprop = "/home/aborrel/ChemMap/structures.csv"
prout = "c://Users/Aborrel/chemmaps/drugbankDesc/PCAs/"
valcor = 0.9
maxquantile = 80

# manually define
# outlier = c("DB00793", "DB00516", "DB06690", "DB09157", "DB03627", "DB01751", "DB03853" ,"DB04711")

#### 1D2D matrix ####
#####################
d1D2D = openData(p1D2D, valcor, prout, c(1,2))
d1D2D_data = d1D2D[[1]]
rownames(d1D2D_data) = d1D2D_data[,1]
d1D2D_data = d1D2D_data[,-1] # remove name
d1D2D_data = d1D2D_data[,-1] # remove SMILES

# remove not well distributed descriptors #
###########################################
d1D2D_data = delnohomogeniousdistribution(d1D2D_data, 80)

p1D2Dclean = paste(prout, "1D2D_clean.csv", sep = "")
write.table(d1D2D_data, file = p1D2Dclean, row.names = TRUE, sep = "\t")




##### 3D matrix #####
#####################
d3D = openData(p3D, valcor, prout, c(1))
d3D_data = d3D[[1]]
rownames(d3D_data) = d3D_data[,1]
d3D_data = d3D_data[,-1] # remove name
# remove not well distributed descriptors #
###########################################
d3D_data = delnohomogeniousdistribution(d3D_data, 80)

# write selected descriptors
p3Dclean = paste(prout, "3D_clean.csv", sep = "")
write.table(d3D_data, file = p3Dclean, row.names = TRUE, sep = "\t")



#################
# merge dataset #
#################

vcompound = rownames(d1D2D_data)
vcompound = intersect(vcompound,rownames(d3D_data))

## case of outlier ##
#####################
#vcompound = vcompound[!vcompound %in% outlier]

dglobal = cbind(d1D2D_data[vcompound,], d3D_data[vcompound,])



#####################################################
#  Analyse descriptor correlation  and distribution #
#####################################################
#cardMatrixCor(cor(cbind(d1D2D_data[vcompound,], d3D_data[vcompound,])), paste(prout, "cor1D2DVS3D", sep = ""), 6)

#plot histogram #
#histDataOne(data1 = d1D2D_data[vcompound,], paste(prout, "homodishist1D2D.pdf"))
#histDataOne(data1 = d3D_data[vcompound,], paste(prout, "homodishist3D.pdf"))


#######################
# analyse projection  #
#######################
# ICA
#generateICAcoords(dglobal, prout)

# PCA 2D
#PCAplot(dglobal, paste(prout, "PCA_DescAll2D", sep = ""))

# PCA 3D
#PCA3D(dglobal, paste(prout, "PCA_DescAll3D", sep = ""))

# PCA combined
#PCAcombined2plans(d1D2D_data[vcompound,], d3D_data[vcompound,], paste(prout, "combined-1D2D_3D"))
generateCoordCombinedPCA(d1D2D_data[vcompound,], d3D_data[vcompound,], prout)

# MDSglobal
#MDS3DGlobal(dglobal, prout, "corr")
#MDS3DGlobal(dglobal, prout, "euclidean")
#MDS3DGlobal(dglobal, prout, "manhattan")

#MDS combined
#MDSMulti(d1D2D, d3D, prout, "corr")
#MDSMulti(d1D2D, d3D, prout, "euclidean")
#MDSMulti(d1D2D, d3D, prout, "manhattan")

# model for AR
#model3D(cbind(d1D2D_data[vcompound,] d3D_data[vcompound,], dcol, prout)




####################### old
### color not used ####
# prop for color
#dprop = read.csv(pprop, sep = "\t", stringsAsFactors = F, header = TRUE)
#dcol = as.matrix(dprop[,c("DRUGBANK_ID","DRUG_GROUPS")])
#dcol = cbind(dcol, rep("grey", dim(dcol)[1]))
#dcol[which(dcol[,2] == "approved"),3] = "green"
#dcol[which(dcol[,2] == "withdrawn"),3] = "red"
#dcol[which(dcol[,2] == "illicit"),3] = "orange"
#rownames(dcol) = dcol[,1]
#dcol = dcol[,-1]


