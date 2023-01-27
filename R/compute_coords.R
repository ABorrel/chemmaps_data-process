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
p_desc = args[1]
type_desc = args[2]
prout = args[3]
valcor = as.double(args[4])
maxquantile = as.integer(args[5])


#############
# Prep desc #
#############

df_desc = prepMatrixDesc(p_desc, valcor, maxquantile, c(1), 1)
vcompound = rownames(df_desc)
df_desc = delSDnull(df_desc) # double check a SD NULL

#########################
#  generate coordinates #
#########################

# scaling
ldescscale = scaling(df_desc)

write.csv(ldescscale[[2]], file = paste(prout, type_desc, "scaling.csv", sep = ""))

dscale = delSDnull(ldescscale[[1]])

# generate PCA -1D2D
lcoord = generatePCAcoords(dscale)
# write cp
write.csv(lcoord[[3]], file=paste(prout, "CP", type_desc, ".csv", sep = ""))
# write coord
dcoord = lcoord[[1]]
colnames(dcoord) = paste("DIM", seq(1,dim(dcoord)[2]), sep = "")
write.csv(dcoord[,1:10], file=paste(prout, "coord", type_desc, ".csv", sep = "")) #take only 10 dimensions
write.csv(lcoord[[2]], file=paste(prout, "VarPlan", type_desc, ".csv", sep = ""))