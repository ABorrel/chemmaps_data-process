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

#########
##  Format for toolchem
#
#cUpdate.formatChemForToolChem(50000)
#cUpdate.formatPrepChemForToolChem()
#cUpdate.formatDesc2DForToolChem()
#cUpdate.formatDesc3DForToolChem()
#cUpdate.formatOPERAForToolChem()

#########
## update OPERA name table
#
#cUpdate.pushOPERANameTable("chem_descriptor_opera_name_new")

#############
## update chemicals in the DB
#
#cUpdate.updateMissingDTXSID("chemicals")# check missing dtxid with smiles
#cUpdate.updateSMILES("chemicals")# check different original smiles from previous version
#cUpdate.updateNameAndCAS("chemicals")# update name and casrn
#cUpdate.updateDescOPERA()# update name and casrn


####
## update new chemicals to push in the DB
#
#cUpdate.extractOnlyNewChem("chemicals", "dsstox_id")
#cUpdate.pushNewChemInDB()
cUpdate.extractChemicalWithDescNoComputed()
cUpdate.computeDescNewChem()
ddd
#############
## update desc in DB
#
