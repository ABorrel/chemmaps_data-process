from os import path
import updateChemDB

# Run functions to update the chemicals database from OPERA files 
PR_DATA = path.abspath("../../data/") + "/"
PR_OUTPUT = path.abspath("../../results/") + "/"

#################
# UPDATE DSSTOX

name_update = "2020-01"
pr_OPERA_pred = PR_DATA + "dsstox_01-2020/"

cUpdate = updateChemDB.updateChemDB(name_update, pr_OPERA_pred, PR_OUTPUT)
#cUpdate.formatChemForToolChem(50000)
#cUpdate.updateChemicalsTableFromOPERAFile("chemicals")


# update DB
#cUpdate.extractOnlyNewChem("chemicals", "dsstox_id")
cUpdate.computeDescNewChem()