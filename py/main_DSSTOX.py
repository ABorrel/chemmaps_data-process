
#############################
# DSSTOX for different list #
#############################

import DSSTOXlib
def generateCoordFromEPAlist(plist, prout, nameMap, computeDesc, computePNG, corval=0.9, maxquantile=90, splitMap=1, istart=0, iend=0, project = 0):

    prDSSTOXPred = "/home/borrela2/ChemMaps/data/DSSTOX_pred/"
    pknownSDF = "/home/borrela2/ChemMaps/data/ToxTrainTest_3D.sdf"
    pLD50 = "/home/borrela2/ChemMaps/data/LD50_data.csv"
    
    #prDSSTOXPred = "C:\\Users\\borrela2\\development\\trash\\DSSTOX_pred\\"
    #pknownSDF = "C:\\Users\\borrela2\\development\\trash\\ToxTrainTest_3D.sdf"
    #pLD50 = "C:\\Users\\borrela2\\development\\trash\\LD50_data.csv"
    #pDSTOXIDmap = "C:\\Users\\borrela2\\development\\trash\\DSSTox_Identifiers_Map.csv"

    prDESC = "/home/borrela2/ChemMaps/data_analysis/DESC/"

    db = DSSTOXlib.DSSTOX(plist, nameMap, istart, iend, prDESC, prout)
    db.loadlistChem()
    #db.pushChemInDB()
    #db.generateTablePropAllDSSTOX(prDSSTOXPred, pknownSDF, pLD50, plist, insertDB=0)
    #db.pushTablePropAllInDB()
    #db.computepng()
    #db.computeDesc(insertDB=1, w=0)
    #db.pushDssToxNamePropInDB()
    if project == 1:
        if nameMap != "dsstox":
            db.runRprojection(corval, maxquantile)
    if nameMap != "dsstox":
        #db.pushPropInDB()
        #db.computeCoords(corval, maxquantile, insertDB=1)
        #db.generateNeighborMatrix(20, [2,1])
        #db.generateNeighborMatrix(20, [])
        db.pushNeighbors()
        #db.updateTableProp(nameMap)
    else:
        db.computeCoords(corval, maxquantile, insertDB=0)
        #db.splitMap(splitMap, 1, insertDB=1)
        #db.splitMap(splitMap, 2, insertDB=1)
        #db.splitMap(splitMap, 3, insertDB=1)
        #db.generateCentroidFile()
        # add here function to put map in the the DB
        #db.generateNeighborMatrix(20, [2,1])
        #db.pushDSSTOXNeighbors(prout + "Neighbors/")
        #print("*")






#############
# prep PFAS #
#############
pPFAS = "/home/borrela2/ChemMaps/data/PFAS/list_chemicals-2019-09-12-10-05-50.csv"
prPFAS = "/home/borrela2/ChemMaps/data_analysis/PFAS/"
generateCoordFromEPAlist(pPFAS, prPFAS, "pfas", computeDesc=0, computePNG=0, corval=0.9, maxquantile=90, splitMap=1, istart=0, iend=0)
ddd

################
# Prep Tox21   #
################
pTOX21 = "/home/borrela2/ChemMaps/data/TOX21/list_chemicals-2019-05-13-10-11-38.csv"
prTOX21 = "/home/borrela2/ChemMaps/data_analysis/TOX21/"
#generateCoordFromEPAlist(pTOX21, prTOX21, "tox21", computeDesc=0, computePNG=0, corval=0.9, maxquantile=90, splitMap=1, istart=0, iend=0)
#################
# prep DSSTOX   #
#################
pDSSTOX = "/home/borrela2/ChemMaps/data/DSSTox_QSAR-r_1-15.csv"
prDSSTOX = "/home/borrela2/ChemMaps/data_analysis/DSSTox/"
#pDSSTOX = "C:/Users/borrela2/development/trash/dsstox/DSSTox_QSAR-r_1-15.csv"
#prDSSTOX = "C:/Users/borrela2/development/trash/dsstox/"

#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=0, computePNG=0, corval=0.9, maxquantile=90, splitMap=100, istart=0, iend=0)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=600000, iend=0)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=0, iend=100000)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=100000, iend=150000)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=150000, iend=200000)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=200000, iend=250000)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=250000, iend=300000)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=300000, iend=350000)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=350000, iend=400000)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=400000, iend=450000)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=450000, iend=500000)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=500000, iend=550000)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=550000, iend=600000)
generateCoordFromEPAlist(pDSSTOX, prDSSTOX, "dsstox", computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=600000, iend=9000000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=650000, iend=700000)
#generateCoordFromEPAlist(pDSSTOX, prDSSTOX, computeDesc=1, computePNG=1, corval=0.9, maxquantile=90, splitMap=1, istart=61100, iend=900000)


# do with
#runPNG("/home/borrela2/ChemMaps/data_analysis/DESC/SMI/", "/home/borrela2/ChemMaps/data_analysis/DESC/PNG/")


# add desc OPERA
import DBrequest
import toolbox

def loadOPERADesc(pOPERA, table):

    ddesc = toolbox.loadMatrixToDict(pOPERA, sep = ",")
    LPROP = list(ddesc[list(ddesc.keys())[0]].keys())
    cDB = DBrequest.DBrequest()
    cDB.verbose = 1
    i = 1
    for PROP in LPROP[1:]:
        cDB.addElement(table, ["id", "name"], [i, PROP])
        i = i + 1

#loadOPERADesc("/home/borrela2/sandbox/VM/ChemMap2Site/temp/949289/OPERA/ZWRUINPWMLAQRD-UHFFFAOYSA-N.csv", "desc_opera_name")