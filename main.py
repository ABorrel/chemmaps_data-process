import pathFolder
import loadDB
import liganddescriptors
import analysis

from os import path, remove
import math


def computeDesc(psdf, prdesc, Desc1D2D=1, generation3D = 0, Desc3D=1, control=2, namek="DATABASE_ID"):

    dout = {}
    dout["1D2D"] = prdesc + "1D2D.csv"
    dout["3D"] = prdesc + "3D.csv"


    # shortcut
    if path.exists(prdesc + "3D.csv") or path.exists(prdesc + "1D2D.csv") and control == 1:
        return dout

    if path.exists(prdesc + "3D.csv") or path.exists(prdesc + "1D2D.csv") and control == 2:
        try:remove(dout["1D2D"])
        except:pass

        try:remove(dout["3D"])
        except:pass

        return dout


    # formate database
    DB = loadDB.sdfDB(psdf, namek, prdesc)
    DB.parseAll()
    DB.renameHeader() # change first line of sdf to wrtie the good name

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
            desc = liganddescriptors.Descriptors(compound, log, prdesc, dout, 1, namek)
            if desc.log == "ERROR":
                continue
            desc.get_descriptor1D2D()
            desc.writeTablesDesc("1D2D")


    # descriptor 3D #
    #################

    if Desc3D == 1:
        i = 0
        for compound in DB.lc:
            print i, "generation 3D or split file"
            if generation3D == 1:
                desc = liganddescriptors.Descriptors(compound, log, prdesc, dout, 1, namek)
                prSDF3D = desc.generate3DFromSMILES(log)
            else:
                # write check == 1 to generate sdf
                desc = liganddescriptors.Descriptors(compound, log, prdesc, dout, 1, namek)
                prSDF3D = prdesc + "SDF/"
            i += 1

        liganddescriptors.get_descriptor3D(prSDF3D, dout["3D"])

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


###################################
# Tox dataset from Kamel Monsouri #
###################################
psdf = "/home/aborrel/ChemMap/generateCords/ToxTrainingSet_3D.sdf"
pranalysis = "/home/aborrel/ChemMap/generateCords/ToxAnalysis/"
pathFolder.createFolder(pranalysis)

prforjsupdate = "/home/aborrel/ChemMap/generateCords/ToxAnalysis/JS/"
pathFolder.createFolder(prforjsupdate)


###################################
# analysis and compute descriptor #
###################################

db = loadDB.sdfDB(psdf, "CASRN", pranalysis)
db.writeTable("ToxTrain.csv")
#db.writeTableSpecific(["DRUG_GROUPS", "GENERIC_NAME", "FORMULA", "MOLECULAR_WEIGHT", "ALOGPS_LOGP", "SYNONYMS"], "tableChemmaps")
db.writeNameforJS(prforjsupdate + "listSearch.js", ["CASRN", "DTXSID", "Name"])



# descriptor computation #
##########################
prDesc = "/home/aborrel/ChemMap/generateCords/ToxAnalysis/Desc/"
pathFolder.createFolder(prDesc, clean = 1)

dpfiledesc = computeDesc(psdf, prDesc, namek="CASRN")
#analysis.PCAplot(dpfiledesc["1D2D"], dpfiledesc["3D"])





dddd

###### HAVE TO UPDATE


#############################################################################################################################
#############################################################################################################################

###################################
#drugbank https://www.drugbank.ca/#
###################################
psdf = "/home/aborrel/ChemMap/generateCords/drugbank-20-12-2017.sdf"
pranalysis = "/home/aborrel/ChemMap/generateCords/drugBankAnalysis/"
pathFolder.createFolder(pranalysis)

prforjsupdate = "/home/aborrel/ChemMap/generateCords/drugBankAnalysis/JS/"
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




