import pathFolder
import formatDatabase
import liganddescriptors
import analysis

from os import path, remove



def computeDesc(psdf):

    dout = {}
    dout["1D"] = ""
    dout["2D"] = ""
    dout["3D"] = ""

    pdesc = pathFolder.PR_RESULT + "Desc"

    if path.exists(pdesc + "3D.csv"):
        dout["3D"] = pdesc + "3D.csv"
    if path.exists(pdesc + "2D.csv"):
        dout["2D"] = pdesc + "2D.csv"
    if path.exists(pdesc + "1D.csv"):
        dout["2D"] = pdesc + "1D.csv"
    return dout

    # formate database
    DB = formatDatabase.drugbank(psdf)
    DB.parseAll()

    print len(DB.lc)


    plog = pathFolder.PR_RESULT + "log.txt"
    log = open(plog, "w")

    # clean
    # if path.exists(pdesc):
    #    remove(pdesc)

    i = 0
    for compound in DB.lc:
        print i

        desc = liganddescriptors.Descriptors(compound, log)
        if desc.log == "ERROR":
            continue

        desc.get_descriptorOD1D()
        desc.get_descriptor2D()
        desc.get_descriptor3D(log)

        desc.writeTablesDesc(pdesc)
        i += 1


###############
##   MAIN    ##
###############

#drugbank https://www.drugbank.ca/
psdf = pathFolder.PR_REF + "structures.sdf"
#pmacrocycle = "/home/aborrel/macrocycle_phyophyo/len_12_0.sdf"

#formate database and descriptor computation
dpfiledesc = computeDesc(psdf)
analysis.PCAplot(dpfiledesc["1D"], dpfiledesc["2D"], dpfiledesc["3D"])

