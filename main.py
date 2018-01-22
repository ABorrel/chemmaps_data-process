import pathFolder
import loadDB
import liganddescriptors
import analysis

from os import path, remove
import math


def computeDesc(psdf, prdesc, Desc1D2D=0, Desc3D=1, control=0):

    dout = {}
    dout["1D2D"] = ""
    dout["3D"] = ""


    if path.exists(prdesc + "3D.csv") and path.exists(prdesc + "1D2D.csv" and control == 1):
        dout["3D"] = prdesc + "3D.csv"
        dout["1D2D"] = prdesc + "1D2D.csv"
        return dout


    # formate database
    DB = loadDB.sdfDB(psdf, prdesc)
    DB.parseAll()
    DB.writeTable()

    print len(DB.lc)


    plog = prdesc + "log.txt"
    log = open(plog, "w")

    # clean
    # if path.exists(pdesc):
    #    remove(pdesc)

    # Descriptor 1D and 2D
    if Desc1D2D == 1:
        i = 0
        for compound in DB.lc:
            print i
            desc = liganddescriptors.Descriptors(compound, log, prdesc)
            if desc.log == "ERROR":
                continue
            desc.get_descriptor1D2D()

    # descriptor 3D #
    #################

    if Desc3D == 1:
        i = 0
        for compound in DB.lc:
            print i

            desc = liganddescriptors.Descriptors(compound, log, prdesc)
            prSDF3D = desc.generate3DFromSMILES(log)
            i += 1

        pdesc3D = prdesc + "3D.csv"
        pdesc3D = liganddescriptors.get_descriptor3D(prSDF3D, pdesc3D)

        dout["3D"] = pdesc3D

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
psdf = "/home/aborrel/ChemMap/updateDrugBank/structures-20-12-2017.sdf"
pranalysis = "/home/aborrel/ChemMap/updateDrugBank/dbanalyse/"
pathFolder.createFolder(pranalysis)

prforjsupdate = "/home/aborrel/ChemMap/updateDrugBank/dbanalyse/JS/"
pathFolder.createFolder(prforjsupdate)


###################################
# analysis and compute descriptor #
###################################

db = loadDB.sdfDB(psdf, pranalysis)
db.writeTable()
db.writeTableSpecific(["DRUG_GROUPS", "GENERIC_NAME", "FORMULA", "MOLECULAR_WEIGHT", "ALOGPS_LOGP", "SYNONYMS"], "tableChemmaps")
db.writeNameforJS(prforjsupdate + "listSearch.js")


# draw compounds #
##################
#prpng = "/home/aborrel/ChemMap/updateDrugBank/dbanalyse/cpdpng/"
#pathFolder.createFolder(prpng)
#db.drawMolecules()



# descriptor computation #
##########################
prDesc = "/home/aborrel/ChemMap/updateDrugBank/dbanalyse/Desc/"
pathFolder.createFolder(prDesc, clean = 1)

dpfiledesc = computeDesc(psdf, prDesc)
#analysis.PCAplot(dpfiledesc["1D2D"], dpfiledesc["3D"])




#####################
# update webserver  #
#####################





#pcord3D = pranalysis + "cor3D.csv"

#extractCloseCompounds(pcord3D, 20, pranalysis + "lneighbor.js")



#formate database and descriptor computation




