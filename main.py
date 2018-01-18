import pathFolder
import loadDB
#import liganddescriptors
import analysis

from os import path, remove
import math


def computeDesc(psdf):

    dout = {}
    dout["1D"] = ""
    dout["2D"] = ""
    dout["3D"] = ""

    pdesc = pathFolder.PR_RESULT + "Desc"

    if path.exists(pdesc + "3D.csv") and path.exists(pdesc + "2D.csv") and path.exists(pdesc + "1D.csv"):
        dout["3D"] = pdesc + "3D.csv"
        dout["2D"] = pdesc + "2D.csv"
        dout["2D"] = pdesc + "1D.csv"
        return dout

    # formate database
    DB = loadDB.sdfDB(psdf)
    DB.parseAll()
    DB.writeTable()

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


def extractCloseCompounds(pfilin, nneighbor, pfilout):

    filin = open(pfilin, "r")
    lcords = filin.readlines()
    filin.close()

    dcor = {}

    for cord in lcords[1:]:
        lelem = cord.strip().split(",")
        ID = lelem[0]
        x = lelem[1]
        y = lelem[2]
        z = lelem[3]
        dcor[ID] = [float(x),float(y),float(z)]

    ddist = {}
    for ID in dcor.keys():
        ddist[ID] = {}
        for ID2 in dcor.keys():
            if ID != ID2:
                ddist[ID][ID2] = math.sqrt(sum([(xi-yi)**2 for xi,yi in zip(dcor[ID], dcor[ID2])]))

        lID = [i[0] for i in sorted(ddist[ID].items(), key=lambda x:x[1])][:nneighbor]
        ddist[ID] = lID


    filout = open(pfilout, "w")
    filout.write("var lneighbor = {")

    lwrite = []
    for ID in ddist.keys():
        w = str(ID) + ":[" + ",".join(ddist[ID]) + "]"
        lwrite.append(w)

    filout.write(",".join(lwrite))
    filout.write("}")
    filout.close()





###############
##   MAIN    ##
###############

#drugbank https://www.drugbank.ca/
psdf = "/home/aborrel/chemmaps/drugData/structures.sdf"
pranalysis = "/home/aborrel/chemmaps/drugData/dbanalyse/"
pathFolder.createFolder(pranalysis)
db = loadDB.sdfDB(psdf, pranalysis)
#db.writeTable()
db.writeTableSpecific(["DRUG_GROUPS", "GENERIC_NAME", "FORMULA", "MOLECULAR_WEIGHT", "ALOGPS_LOGP", "SYNONYMS"], "tableChemmaps")
#db.writeNameforJS(pranalysis + "list.js")
#db.drawMolecules()


#pcord3D = pranalysis + "cor3D.csv"

#extractCloseCompounds(pcord3D, 20, pranalysis + "lneighbor.js")



#formate database and descriptor computation
#dpfiledesc = computeDesc(psdf)
#analysis.PCAplot(dpfiledesc["1D"], dpfiledesc["2D"], dpfiledesc["3D"])



