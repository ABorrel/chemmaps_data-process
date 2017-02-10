import pathFolder
import formatDatabase
import liganddescriptors






###############
##   MAIN    ##
###############

# drugbank https://www.drugbank.ca/
psdf = pathFolder.PR_REF + "structures.sdf"

# formate database
DB = formatDatabase.drugbank(psdf)
DB.parseAll()

print len(DB.lc)

for compound in DB.lc:

    desc = liganddescriptors.Descriptors(compound)
    desc.get_descriptorOD1D()
    desc.get_descriptor2D()
    desc.get_descriptor3D()
    desc.writeTablesDesc(pathFolder.PR_RESULT)







