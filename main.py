import pathFolder
import loadDB
import liganddescriptors
import analysis

from os import path, remove
import math


def computeDesc(psdf, prdesc):

    dout = {}
    dout["1D2D"] = ""
    dout["3D"] = ""


    if path.exists(prdesc + "3D.csv") and path.exists(prdesc + "1D2D.csv"):
        dout["3D"] = prdesc + "3D.csv"
        dout["1D2D"] = prdesc + "1D2D.csv"
        return dout

    # formate database
    DB = loadDB.sdfDB(psdf)
    DB.parseAll()
    DB.writeTable()

    print len(DB.lc)


    plog = prdesc + "log.txt"
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

        desc.get_descriptor1D2D()
        desc.get_descriptor3D(log)

        desc.writeTablesDesc(prdesc)
        i += 1
    log.close()

    return dout


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

#drugbank https://www.drugbank.ca/#
###################################
psdf = "/home/aborrel/chemmaps/drugData/structures.sdf"
pranalysis = "/home/aborrel/chemmaps/drugData/dbanalyse/"
prforjsupdate = "/home/aborrel/chemmaps/drugData/dbanalyse/JS/"
pathFolder.createFolder(prforjsupdate)


###################################
# analysis and compute descriptor #
###################################

pathFolder.createFolder(pranalysis)
db = loadDB.sdfDB(psdf, pranalysis)
db.writeTable()
db.writeTableSpecific(["DRUG_GROUPS", "GENERIC_NAME", "FORMULA", "MOLECULAR_WEIGHT", "ALOGPS_LOGP", "SYNONYMS"], "tableChemmaps")
db.writeNameforJS(prforjsupdate + "listSearch.js")

# draw compound
prpng = "/home/aborrel/chemmaps/drugData/dbanalyse/cpdpng/"
pathFolder.createFolder(prpng)
db.drawMolecules()



# descriptor computation #
##########################
prcoords = "/home/aborrel/chemmaps/drugData/dbanalyse/desc/"
pathFolder.createFolder()

dpfiledesc = computeDesc(psdf)
analysis.PCAplot(dpfiledesc["1D2D"], dpfiledesc["3D"])



#####################
# update webserver  #
#####################





#pcord3D = pranalysis + "cor3D.csv"

#extractCloseCompounds(pcord3D, 20, pranalysis + "lneighbor.js")



#formate database and descriptor computation




