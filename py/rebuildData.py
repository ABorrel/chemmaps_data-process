
# script use to rebuild the data for chemmaps in the lastest version
#####

from os import path
import pathFolder
import DrugBank
import DSSTOX
import PFAS
import Tox21
import CATMOS

# format the input for modeling
PDIR_SCRIPT = path.realpath(path.dirname(__file__)) + "\\"
PDIR_ROOT = path.abspath(PDIR_SCRIPT + "..\\..\\") + "\\"
PDIR_INPUT = PDIR_ROOT + "input\\"
PDIR_OUTPUT = PDIR_ROOT + "output\\"

# create output
pathFolder.createFolder(PDIR_OUTPUT)
PR_DESC = pathFolder.createFolder(PDIR_OUTPUT + "DESC/")


# push CATMOS #
###############
## function to push CATMOS acute tox
p_CATMOS = PDIR_INPUT + "CATMOS\\CATMOS.xlsx"
c_CATMOS = CATMOS.CATMOS(p_CATMOS)
c_CATMOS.load_data()
# c_CATMOS.pushDB_name()
c_CATMOS.pushDB_Tox()
sss


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
# cDsstox.pushDB_casn()

# PFAS #
########
p_smiles_pfas = PDIR_INPUT + "pfas/PFASSTRUCT-2023-03-01.csv"
p_dir_PFAS = pathFolder.createFolder(PDIR_OUTPUT + "PFAS/")
# cPFAS = PFAS.PFAS(p_smiles_pfas, "pfas", 0, 0, PR_DESC, p_dir_PFAS)
# cPFAS.pushDB_chem()
# cPFAS.compute_pushDB_Desc()
# cPFAS.compute_coords(0.9, 90) # coord are compute in R
# cPFAS.pushDB_coords()
# cPFAS.draw_map()
# cPFAS.compute_onDB_neighbors() 


# Tox21 #
#########

p_tsv = PDIR_INPUT + "tox21/tox21_10k_library.tsv"
p_dir_Tox21 = pathFolder.createFolder(PDIR_OUTPUT + "Tox21/")
p_annotation = PDIR_INPUT + "tox21/CCD-Batch-Search_2023-03-02_07_29_30.csv"
cTox21 = Tox21.Tox21(p_tsv, p_annotation, "tox21", 0, 0, PR_DESC, p_dir_Tox21)
# cTox21.pushDB_chem()
cTox21.compute_pushDB_Desc()
