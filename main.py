import pathFolder
import formatDatabase
import liganddescriptors

from os import path, remove




###############
##   MAIN    ##
###############

# drugbank https://www.drugbank.ca/
psdf = pathFolder.PR_REF + "structures.sdf"
#pmacrocycle = "/home/aborrel/macrocycle_phyophyo/len_12_0.sdf"

# formate database
DB = formatDatabase.drugbank(psdf)
DB.parseAll()

print len(DB.lc)

pdesc = pathFolder.PR_RESULT + "drugbank.desc"

#clean
if path.exists(pdesc):
    remove(pdesc)

for compound in DB.lc:

    desc = liganddescriptors.Descriptors(compound)
    #desc.get_descriptorOD1D()
    #desc.get_descriptor2D()
    desc.get_descriptor3D()

    desc.writeTablesDesc(pdesc)







