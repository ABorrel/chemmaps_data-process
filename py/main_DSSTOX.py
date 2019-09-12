
#############################
# DSSTOX for different list #
#############################

import DSSTOXlib
def generateCoordFromEPAlist(plist, prout, computeDesc, computePNG, corval=0.9, maxquantile=90, splitMap=1, istart=0, iend=0):

    prDSSTOXPred = "/home/borrela2/ChemMaps/data/DSSTOX_pred/"
    pknownSDF = "/home/borrela2/ChemMaps/data/ToxTrainTest_3D.sdf"
    pLD50 = "/home/borrela2/ChemMaps/data/LD50_data.csv"
    pDSTOXIDmap = "/home/borrela2/ChemMaps/data/DSSTox_Identifiers_Map.csv"
    prDESC = "/home/borrela2/ChemMaps/data_analysis/DESC/"

    db = DSSTOXlib.DSSTOX(plist, istart, iend, prDESC, prout)
    db.loadlistChem()
    db.computeDesc(insertDB =0, w=1)
    db.runRprojection(corval, maxquantile)
    db.computeCoords(corval, maxquantile, insertDB=0)


    #db.writeDescMatrix("1D2D")
    #db.writeDescMatrix("3D")
    #if computePNG == 1:
    #    db.computepng()
    # only if compute all map
    #if iend == 0:
    #    db.generateFileMap(corval, maxquantile)
    #    db.projection(corval, maxquantile)
    #    db.splitMap(splitMap, 1)
    #    db.generateMapSplitFile(1)
    #    db.generateTableProp(prDSSTOXPred, pknownSDF, pLD50, pDSTOXIDmap)
    #    db.generateNeighborMatrix(20)
    #    db.splitMap(splitMap, 2)
    #    db.generateMapSplitFile(2)
    #    db.generateTableProp(prDSSTOXPred, pknownSDF, pLD50, pDSTOXIDmap)
    #    db.generateNeighborMatrix(20)
    #    db.splitMap(splitMap, 3)
    #    db.generateMapSplitFile(3)
    #    db.generateTableProp(prDSSTOXPred, pknownSDF, pLD50, pDSTOXIDmap)
    #    db.generateNeighborMatrix(20)




#############
# prep PFAS #
#############
#pPFAS = "/home/borrela2/ChemMaps/data/PFAS/list_chemicals-2019-09-12-10-05-50.csv"
#prPFAS = "/home/borrela2/ChemMaps/data_analysis/PFAS/"
#generateCoordFromEPAlist(pPFAS, prPFAS, computeDesc=1, computePNG=0, corval=0.9, maxquantile=90, splitMap=1, istart=0, iend=0)



################
# Prep Tox21   #
################
pTOX21 = "/home/borrela2/ChemMaps/data/TOX21/list_chemicals-2019-05-13-10-11-38.csv"
prTOX21 = "/home/borrela2/ChemMaps/data_analysis/TOX21/"

#generateCoordFromEPAlist(pTOX21, prTOX21, computeDesc=1, computePNG=0, corval=0.9, maxquantile=90, splitMap=1, istart=0, iend=0)



#################
# prep DSSTOX   #
#################
pDSSTOX = "/home/borrela2/ChemMaps/data/DSSTox_QSAR-r_1-15.csv"
prDSSTOX = "/home/borrela2/ChemMaps/data_analysis/DSSTox/"

#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=0, computePNG=0, corval=0.9, maxquantile=90, splitMap=100, istart=0, iend=0)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=0, iend=20000)
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



# do witj
#runPNG("/home/borrela2/ChemMaps/data_analysis/DESC/SMI/", "/home/borrela2/ChemMaps/data_analysis/DESC/PNG/")
