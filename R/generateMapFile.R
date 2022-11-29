source("tool.R")



scaling = function(din){
  
  dinScale = scale(din)
  lscale = attr(dinScale, "scaled:scale")
  lcenter = attr(dinScale, "scaled:center")
  
  dforscaling = rbind(lscale, lcenter)
  rownames(dforscaling) = c("scale", "center")
  
  return(list(dinScale, dforscaling))
  
}


generatePCAcoords = function(dinScale){
  
  data.cor=cor(dinScale)
  data.eigen=eigen(data.cor)
  lambda = data.eigen$values
  var_cap = lambda/sum(lambda)*100
  cp = data.eigen$vectors
  rownames (cp) = colnames (dinScale)
  colnames (cp) = colnames (dinScale)
  data_plot = as.matrix(dinScale)%*%cp
  rownames(data_plot) = rownames(dinScale)
  
  return(list(data_plot, var_cap, cp))
}



prepMatrixDesc = function(pin, valcor, maxquantile, vexclude, homo){
  
  #### 1D2D matrix ####
  #####################
  din = openData(pin, valcor, "", vexclude)
  din_data = din[[1]]
  rownames(din_data) = din_data[,1]
  din_data = din_data[,-1] # remove name
  
  # remove not well distributed descriptors #
  ###########################################
  if(homo == 1){
    din_data = delnohomogeniousdistribution(din_data, maxquantile)
  }else{
    din_temp = apply(din_data,2,as.double)
    rownames(din_temp) = rownames(din_data)
    din_data = din_temp
  }
  din_data = na.omit(din_data)
  
  return(din_data)
  
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p1D2D = args[1]
p3D = args[2]
prout = args[3]
valcor = as.double(args[4])
maxquantile = as.integer(args[5])

# dsstoxmap
#p1D2D = "c://Users/aborr/research/sandbox/chemmaps_data-process/results/updateDSSTOX/coords/1D2D.csv"
#p3D = "c://Users/aborr/research/sandbox/chemmaps_data-process/results/updateDSSTOX/coords/3D.csv"
#prout = "c://Users/aborr/research/sandbox/chemmaps_data-process/results/updateDSSTOX/coords/proj_0.9-90/"
#valcor = 0.9
#maxquantile = 90


# drugmap
#p3D = "/home/borrela2/ChemMaps/data_analysis/drugBankAnalysis/Desc/3D.csv"
#p1D2D = "/home/borrela2/ChemMaps/data_analysis/drugBankAnalysis/Desc/1D2D.csv"
#valcor = 0.9
#maxquantile = 95
#prout = "/home/borrela2/ChemMaps/data_analysis/drugBankAnalysis/Desc/map/"

d1D2D = prepMatrixDesc(p1D2D, valcor, maxquantile, c(1), 1)
d3D = prepMatrixDesc(p3D, valcor, maxquantile, c(1), 1)

#################
# merge dataset #
#################

vcompound = rownames(d1D2D)
vcompound = intersect(vcompound,rownames(d3D))

d1D2D = d1D2D[vcompound,]
d3D = d3D[vcompound,]

d1D2D = delSDnull(d1D2D)
d3D = delSDnull(d3D)
#########################
#  generate coordinates #
#########################

# scaling
ld1D2Dscale = scaling(d1D2D)
ld3Dscale = scaling(d3D)

write.csv(ld1D2Dscale[[2]], file = paste(prout, "1D2Dscaling.csv", sep = ""))
write.csv(ld3Dscale[[2]], file = paste(prout, "3Dscaling.csv", sep = ""))

dscale1D2D = delSDnull(ld1D2Dscale[[1]])

# generate PCA -1D2D
lcoord = generatePCAcoords(dscale1D2D)
# write cp
write.csv(lcoord[[3]], file=paste(prout, "CP1D2D.csv", sep = ""))
# write coord
dcoord = lcoord[[1]]
colnames(dcoord) = paste("DIM", seq(1,dim(dcoord)[2]), sep = "")
write.csv(dcoord[,1:10], file=paste(prout, "coord1D2D.csv", sep = ""))#plot only 10 dimensions
write.csv(lcoord[[2]], file=paste(prout, "VarPlan1D2D.csv", sep = ""))

dscale3D = delSDnull(ld3Dscale[[1]])

# generate PCA -3D
lcoord = generatePCAcoords(dscale3D)
# write cp
write.csv(lcoord[[3]], file=paste(prout, "CP3D.csv", sep = ""))
# write coord
dcoord = lcoord[[1]]
colnames(dcoord) = paste("DIM3-", seq(1,dim(dcoord)[2]), sep = "")
write.csv(dcoord[,1:10], file=paste(prout, "coord3D.csv", sep = ""))
write.csv(lcoord[[2]], file=paste(prout, "VarPlan3D.csv", sep = ""))
