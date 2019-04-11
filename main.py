import pathFolder
import loadDB
import computeDB
import runExternalSoft
import toolbox
import createJS

from os import path, remove, listdir
import math

from formatDatabase import drugbank

import prepChemLib# to replace main

#########################################
# Tox train dataset from Kamel Monsouri #
#########################################

#psdf = "/home/aborrel/ChemMap/generateCords/ToxTrainingSet_3D.sdf"
#pranalysis = "/home/aborrel/ChemMap/generateCords/ToxAnalysis/"
#kname = "CASRN"

#main(psdf, pranalysis, kname, Desc1D2D=1, generation3D=0, Desc3D=1)

#import descriptors3D
#print descriptors3D.get3Ddesc("/home/aborrel/ChemMap/generateCords/ToxAnalysisGlobal/Desc/SDF/361442-04-8.sdf")


##########################################
# Tox global dataset from Kamel Monsouri #
##########################################

#psdf = "/home/borrela2/ChemMaps/data_analysis/ToxTrainTest_3D.sdf"
#pranalysis = "/home/borrela2/ChemMaps/data_analysis/ToxAnalysisGlobal/"
#kname = "CASRN"
#lkinfo = ["CASRN", "LD50_mgkg", "GHS_category", "EPA_category", "Weight", "LogP"]

#corval = 0.9
#maxquantile = 90

#prepChemLib.Run(psdf, pranalysis, kname, corval=corval, maxquantile=maxquantile, lkinfo=lkinfo, Desc1D2D=0,
#                generation3D=0, Desc3D=0, projection=1, JS=1)


###################################
#drugbank https://www.drugbank.ca/#
###################################

import prepChemLib

psdf = "/home/borrela2/ChemMaps/data/drugbank-20-12-2018.sdf"
pranalysis = "/home/borrela2/ChemMaps/data_analysis/drugBankAnalysisV2018/"
kname = "DATABASE_ID"

corval = 0.9
maxquantile = 90
prepChemLib.Run(psdf, pranalysis, kname, corval=corval, maxquantile=maxquantile, Desc1D2D=0,
                generation3D=0, Desc3D=0, projection=0, map=1)
jjj

##############################
# by list from EPA dashboard #
##############################

#prListChemical = "/home/borrela2/ChemMaps/data/toxEPA-lists/"
#llistChem = listdir(prListChemical)

#kname = "<CASRN>"
#ksmile = "SMILES"
#corval = 0.9
#maxquantile = 90


#for listChem in llistChem:
#    plist = prListChemical + listChem
#    print plist
    # main(psdf, pranalysis, kname, corval=corval, maxquantile=maxquantile, lkinfo=lkinfo, Desc1D2D=0, generation3D=0, Desc3D=0, projection=1, JS=1)


    # have to had generation 3D
    # acess





#################
# prep DSSTOX   #
#################
import DSSTOXlib

pDSSTOX = "/home/borrela2/ChemMaps/data/DSSTox_QSAR-r_1-15.csv"
prDSSTOX = "/home/borrela2/ChemMaps/data/DSSTox/"
prDSSTOXPred = "/home/borrela2/ChemMaps/data/DSSTOX_pred/"
pknownSDF = "/home/borrela2/ChemMaps/data/ToxTrainTest_3D.sdf"
pLD50 = "/home/borrela2/ChemMaps/data/LD50_data.csv"
corval = 0.9
maxquantile = 90

db = DSSTOXlib.DSSTOX(pDSSTOX, 1, 0, prDSSTOX)
db.computeDesc2D3D(compute=0)
db.writeDescMatrix("1D2D")
db.writeDescMatrix("3D")
db.generateFileMap(corval, maxquantile)
db.splitMap(64)
db.generateTableProp(prDSSTOXPred, pknownSDF, pLD50)
db.generateNeighborMatrix(20)


#db = DSSTOXlib.DSSTOX(pDSSTOX, 50000, 100000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 100000, 150000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 150000, 200000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 200000, 250000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 250000, 300000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 300000, 350000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 350000, 400000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 400000, 450000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 450000, 500000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 500000, 550000, prDSSTOX)
#b.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 550000, 600000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 650000, 700000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 700000, 750000, prDSSTOX)
#db.computeDesc2D3D()

#db = DSSTOXlib.DSSTOX(pDSSTOX, 750000, 0, prDSSTOX)
#db.computeDesc2D3D()



# NEW DSSTOX
#pDSSTOX = "/home/borrela2/ChemMaps/data/DSSTox_NewChems_20181129_QSAR-r_1-2.csv"
#db = DSSTOXlib.DSSTOX(pDSSTOX, 1, 0, prDSSTOX)
#db.computeDesc2D3D()