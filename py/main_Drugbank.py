###################################
#drugbank https://www.drugbank.ca/#
###################################

import DrugBank

psdf = "/home/borrela2/ChemMaps/data/drugbank-02-07-2019.sdf"
#psdf = "C:\\Users\\borrela2\\development\\trash\\drugbank-02-07-2019.sdf"

pranalysis = "/home/borrela2/ChemMaps/data_analysis/drugBankAnalysisV2019/"
#pranalysis = "C:\\Users\\borrela2\\development\\trash\\drugBankAnalysisV2019\\"

kname = "DATABASE_ID"
prDESC = "/home/borrela2/ChemMaps/data_analysis/DESC/"
#prDESC = "C:\\Users\\borrela2\\development\\trash\\DESC\\"

corval = 0.9
maxquantile = 90


cDrugBank = DrugBank.DrugBank(psdf, prDESC, pranalysis)
#cDrugBank.parseSDFDB()
cDrugBank.computePNG()
#cDrugBank.pushChemInDB()
#cDrugBank.pushDrugBankNamePropInDB()
#cDrugBank.pushPropInDB()
#cDrugBank.computeDesc(insertDB=1)
#cDrugBank.pushDescNameInDB("1D2D")
#cDrugBank.pushDescNameInDB("3D")
#cDrugBank.runRprojection(corval, maxquantile)
#cDrugBank.computeCoords(corval, maxquantile, insertDB=1)


#nbNeighbor = 20
#lnDim = [2,1] # dimension in the space with 1D2D desc and
#lnDim = []
#cDrugBank.neighbormatrix(nbNeighbor, lnDim, insertDB=1)
#cDrugBank.pushNeighbors()