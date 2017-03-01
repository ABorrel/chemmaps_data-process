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

pdesc = pathFolder.PR_RESULT + "Desc"
plog = pathFolder.PR_RESULT + "log.txt"
log = open(plog, "w")

#clean
#if path.exists(pdesc):
#    remove(pdesc)

i = 0
for compound in DB.lc:
    print i
    desc = liganddescriptors.Descriptors(compound, log)
    if desc.log == "ERROR":
        continue
    desc.get_descriptorOD1D()
    desc.get_descriptor2D()
    #desc.get_descriptor3D(log)

    desc.writeTablesDesc(pdesc)
    i += 1


log.close()




