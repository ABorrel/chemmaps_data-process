###################################
#drugbank https://www.drugbank.ca/#
###################################

import DrugBank

psdf = "/home/borrela2/ChemMaps/data/drugbank-20-12-2018.sdf"
psdf = "C:\\Users\\borrela2\\development\\trash\\drugbank-02-07-2019.sdf"

pranalysis = "/home/borrela2/ChemMaps/data_analysis/drugBankAnalysisV2019/"
pranalysis = "C:\\Users\\borrela2\\development\\trash\\drugBankAnalysisV2019\\"

kname = "DATABASE_ID"
prDESC = "/home/borrela2/ChemMaps/data_analysis/DESC/"
prDESC = "C:\\Users\\borrela2\\development\\trash\\DESC\\"

corval = 0.9
maxquantile = 90

cDrugBank = DrugBank.DrugBank(psdf, prDESC, pranalysis)
cDrugBank.parseSDFDB()
cDrugBank.computeDesc()





##################################################
eee

prepChemLib.SDFDrugBanktoDB(psdf, prDESC, pranalysis)
dd
prepChemLib.Run(psdf, pranalysis, kname, corval=corval, maxquantile=maxquantile, Desc1D2D=0,
                generation3D=0, Desc3D=0, projection=0, map=0)



eee














import pathFolder
import loadDB
import computeDB
import runExternalSoft
import toolbox
import createJS

from os import path, remove, listdir
import math
from random import shuffle

from formatDatabase import drugbank

import prepChemLib# to replace main



# function for list in the EPA Comptox

import DSSTOXlib

def generateCoordFromEPAlist(plist, prout, computeDesc, computePNG, corval=0.9, maxquantile=90, splitMap=1, istart=0, iend=0):

    prDSSTOXPred = "/home/borrela2/ChemMaps/data/DSSTOX_pred/"
    pknownSDF = "/home/borrela2/ChemMaps/data/ToxTrainTest_3D.sdf"
    pLD50 = "/home/borrela2/ChemMaps/data/LD50_data.csv"
    pDSTOXIDmap = "/home/borrela2/ChemMaps/data/DSSTox_Identifiers_Map.csv"
    prDESC = "/home/borrela2/ChemMaps/data_analysis/DESC/"

    db = DSSTOXlib.DSSTOX(plist, istart, iend, prDESC, prout)
    db.computeDesc2D3D(compute=computeDesc)
    db.writeDescMatrix("1D2D")
    db.writeDescMatrix("3D")
    if computePNG == 1:
        db.computepng()
    # only if compute all map
    if iend == 0:
        db.generateFileMap(corval, maxquantile)
        db.projection(corval, maxquantile)
        db.splitMap(splitMap, 1)
        db.generateMapSplitFile(1)
        db.generateTableProp(prDSSTOXPred, pknownSDF, pLD50, pDSTOXIDmap)
        db.generateNeighborMatrix(20)
        db.splitMap(splitMap, 2)
        db.generateMapSplitFile(2)
        db.generateTableProp(prDSSTOXPred, pknownSDF, pLD50, pDSTOXIDmap)
        db.generateNeighborMatrix(20)
        db.splitMap(splitMap, 3)
        db.generateMapSplitFile(3)
        db.generateTableProp(prDSSTOXPred, pknownSDF, pLD50, pDSTOXIDmap)
        db.generateNeighborMatrix(20)




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



#############
# prep PFAS #
#############
pPFAS = "/home/borrela2/ChemMaps/data/PFAS/All_PFAS_Compounds_From_EPA.csv"
prPFAS = "/home/borrela2/ChemMaps/data_analysis/PFAS/"
#generateCoordFromEPAlist(pPFAS, prPFAS, computeDesc=0, computePNG=0, corval=0.9, maxquantile=90, splitMap=1, istart=0, iend=0)


################
# Prep Tox21   #
################
pTOX21 = "/home/borrela2/ChemMaps/data/TOX21/list_chemicals-2019-05-13-10-11-38.csv"
prTOX21 = "/home/borrela2/ChemMaps/data_analysis/TOX21/"

#generateCoordFromEPAlist(pTOX21, prTOX21, computeDesc=0, computePNG=0, corval=0.9, maxquantile=90, splitMap=1, istart=0, iend=0)
#generateCoordFromEPAlist(pTOX21, prTOX21, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=1000, iend=2000)
#generateCoordFromEPAlist(pTOX21, prTOX21, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=2000, iend=3000)
#generateCoordFromEPAlist(pTOX21, prTOX21, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=3000, iend=4000)
#generateCoordFromEPAlist(pTOX21, prTOX21, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=4000, iend=5000)
#generateCoordFromEPAlist(pTOX21, prTOX21, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=5000, iend=6000)
#generateCoordFromEPAlist(pTOX21, prTOX21, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=6000, iend=7000)
#generateCoordFromEPAlist(pTOX21, prTOX21, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=7000, iend=9000)

#################
# prep DSSTOX   #
#################
pDSSTOX = "/home/borrela2/ChemMaps/data/DSSTox_QSAR-r_1-15.csv"
prDSSTOX = "/home/borrela2/ChemMaps/data_analysis/DSSTox/"

generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=0, computePNG=0, corval=0.9, maxquantile=90, splitMap=100, istart=0, iend=0)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=650000, iend=700000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=100000, iend=150000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=150000, iend=200000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=200000, iend=250000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=250000, iend=300000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=300000, iend=350000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=350000, iend=400000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=400000, iend=450000)

#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=450000, iend=500000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=500000, iend=550000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=550000, iend=600000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=600000, iend=9000000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=650000, iend=700000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=61100, iend=900000)


#runPNG("/home/borrela2/ChemMaps/data_analysis/DESC/SMI/", "/home/borrela2/ChemMaps/data_analysis/DESC/PNG/")
