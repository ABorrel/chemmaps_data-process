
# script use to rebuild the data for chemmaps in the lastest version
#####

from os import path
import pathFolder
import DrugBank
import DSSTOX
import PFAS

# format the input for modeling
PDIR_SCRIPT = path.realpath(path.dirname(__file__)) + "\\"
PDIR_ROOT = path.abspath(PDIR_SCRIPT + "..\\..\\") + "\\"
PDIR_INPUT = PDIR_ROOT + "input\\"
PDIR_OUTPUT = PDIR_ROOT + "output\\"

# create output
pathFolder.createFolder(PDIR_OUTPUT)
PR_DESC = pathFolder.createFolder(PDIR_OUTPUT + "DESC/")



# DrugBank #
# ########## 
# p_sdf_all = PDIR_INPUT + "drugbank/drugbank_structures_11-2022.sdf"
# p_dir_drugbank = pathFolder.createFolder(PDIR_OUTPUT + "DRUGBANK/")
# cDrugBank = DrugBank.DrugBank(p_sdf_all, PR_DESC, p_dir_drugbank)
# cDrugBank.pushDB_all()



# DSSTOX #
##########
# p_smiles_dsstox = PDIR_INPUT + "dsstox/DSSTox_082021_Structures.csv"
# p_mapping = PDIR_INPUT + "dsstox/DSSTox_082021_IDs.csv"
# p_dir_dsstox = pathFolder.createFolder(PDIR_OUTPUT + "DSSTOX/")
# cDsstox = DSSTOX.DSSTOX(p_smiles_dsstox, p_mapping, "dsstox", 0, 0, PR_DESC, p_dir_dsstox)
# cDsstox.pushDB_all()



# PFAS #
########
p_smiles_pfas = PDIR_INPUT + "pfas/PFASSTRUCT-2023-03-01.csv"
p_dir_PFAS = pathFolder.createFolder(PDIR_OUTPUT + "PFAS/")
cPFAS = PFAS.PFAS(p_smiles_pfas, "pfas", 0, 0, PR_DESC, p_dir_PFAS)
# cPFAS.pushDB_chem()
cPFAS.compute_pushDB_Desc()
# cPFAS.compute_coords(0.9, 90) # coord are compute in R
# cPFAS.pushDB_coords()
# cPFAS.draw_map()
# cPFAS.compute_onDB_neighbors() 