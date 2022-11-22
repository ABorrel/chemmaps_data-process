
# script use to rebuild the data for chemmaps in the lastest version
#####

from os import path
import pathFolder
import DrugBank

# format the input for modeling
PDIR_SCRIPT = path.realpath(path.dirname(__file__)) + "\\"
PDIR_ROOT = path.abspath(PDIR_SCRIPT + "..\\..\\") + "\\"
PDIR_INPUT = PDIR_ROOT + "input\\"
PDIR_OUTPUT = PDIR_ROOT + "output\\"

# create output
pathFolder.createFolder(PDIR_OUTPUT)
PR_DESC = pathFolder.createFolder(PDIR_OUTPUT + "DESC/")



# start with the drugbank
p_sdf_all = PDIR_INPUT + "drugbank_structures_11-2022.sdf"
p_dir_drugbank = pathFolder.createFolder(PDIR_OUTPUT + "DRUGBANK/")

cDrugBank = DrugBank.DrugBank(p_sdf_all, PR_DESC, p_dir_drugbank)
cDrugBank.parse_SDF()
cDrugBank.pushDB_prop_name()
cDrugBank.pushDB_chem()
cDrugBank.pushDB_desc_name() # that work be for every maps
cDrugBank.compute_pushDB_Desc()
here



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