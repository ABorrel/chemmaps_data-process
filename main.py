import pathFolder
import loadDB
import liganddescriptors
import runExternalSoft
import toolbox

from os import path, remove, listdir
import math


def computeDesc(psdf, prdesc, Desc1D2D=1, generation3D = 0, Desc3D=1, control=0, namek="DATABASE_ID"):

    dout = {}
    dout["1D2D"] = prdesc + "1D2D.csv"
    dout["3D"] = prdesc + "3D.csv"


    # shortcut
    if path.exists(prdesc + "3D.csv") and path.exists(prdesc + "1D2D.csv") and Desc3D == 0 and Desc1D2D == 0:
        return dout

    if control == 2:
        try:remove(dout["1D2D"])
        except:pass

        try:remove(dout["3D"])
        except:pass


    if control == 0:
        if Desc1D2D == 1:
            try: remove(dout["1D2D"])
            except: pass
        if Desc3D == 1:
            try: remove(dout["3D"])
            except: pass


    # formate database
    DB = loadDB.sdfDB(psdf, namek, prdesc)
    DB.parseAll()
    DB.renameHeader() # change first line of sdf to write the good name

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
            desc = liganddescriptors.Descriptors(compound, log, prdesc, dout, namek)
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
                desc = liganddescriptors.Descriptors(compound, log, prdesc, dout, namek)
                prSDF3D = desc.generate3DFromSMILES(log)
            else:
                # write check == 1 to generate sdf
                desc = liganddescriptors.Descriptors(compound, log, prdesc, dout, namek)
                prSDF3D = prdesc + "SDF/"
            i += 1
        # rename header
        if generation3D == 1:
            lsdf3D = listdir(prSDF3D)
            for sdf3D in lsdf3D:
                toolbox.renameHeaderSDF(prSDF3D + str(sdf3D))
        liganddescriptors.get_descriptor3DPadel(prSDF3D, dout["3D"])
    log.close()
    return dout







###############
##   MAIN    ##
###############


def main(psdf, pranalysis, kname, control=1, Desc1D2D=1, generation3D = 1, Desc3D=1):

    pathFolder.createFolder(pranalysis)

    ###################################
    # analysis and compute descriptor #
    ###################################

    db = loadDB.sdfDB(psdf, kname, pranalysis)
    db.writeTable("TaleProp.csv")
    #db.writeTableSpecific(["DRUG_GROUPS", "GENERIC_NAME", "FORMULA", "MOLECULAR_WEIGHT", "ALOGPS_LOGP", "SYNONYMS"], "tableChemmaps")
    #db.writeNameforJS(prforjsupdate + "listSearch.js", ["CASRN", "DTXSID", "Name"])


    # descriptor computation #
    ##########################
    prDesc = pranalysis + "Desc/"
    pathFolder.createFolder(prDesc, clean=0)
    dpfiledesc = computeDesc(psdf, prDesc, control=control, Desc1D2D=Desc1D2D, generation3D = generation3D, Desc3D=Desc3D, namek=kname)

    # analyse projection  and compute coordinate #
    ##############################################
    prproject = pathFolder.createFolder(pranalysis + "projection/")
    corval = 0.9
    maxQuantile = 80
    runExternalSoft.RComputeCor(dpfiledesc["1D2D"], dpfiledesc["3D"], prproject, corval, maxQuantile)

    ###################
    # for the website #
    ###################
    prforjsupdate = pranalysis + "JS/"
    pathFolder.createFolder(prforjsupdate)


    # 1. compute png #
    ##################
    prpng = pranalysis + "cpdpng/"
    pathFolder.createFolder(prpng)

    # draw from desc descriptors
    prSDF = prDesc + "SDF/"
    if path.exists(prSDF):
        lsdfs = listdir(prSDF)
        # control if nSDF = nPNG
        if len(lsdfs) != len(listdir(prpng)):
            for sdfile in lsdfs:
                runExternalSoft.molconvert(prSDF + sdfile, prpng + sdfile.split(".")[0] + ".png")
    else:
        db.drawMolecules(prpng)


    # 2. update JS coords #
    #######################

    # 3. update JS properties #
    ###########################

    # 4. update neighborhood #
    ##########################



#########################################
# Tox train dataset from Kamel Monsouri #
#########################################

psdf = "/home/aborrel/ChemMap/generateCords/ToxTrainingSet_3D.sdf"
pranalysis = "/home/aborrel/ChemMap/generateCords/ToxAnalysis/"
kname = "CASRN"

#main(psdf, pranalysis, kname, Desc1D2D=1, generation3D=0, Desc3D=1)





##########################################
# Tox global dataset from Kamel Monsouri #
##########################################

psdf = "/home/aborrel/ChemMap/generateCords/ToxGlobal_3D.sdf"
pranalysis = "/home/aborrel/ChemMap/generateCords/ToxAnalysisGlobal/"
kname = "CASRN"

#main(psdf, pranalysis, kname, control =1, Desc1D2D=0, generation3D=0, Desc3D=1)







# test 3D desciptor
import descriptors3D

d = descriptors3D.get3Ddesc("/home/aborrel/ChemMap/generateCords/ToxAnalysisGlobal/Desc/SDF/28610-84-6.sdf")
print d
d = descriptors3D.get3Ddesc("/home/aborrel/ChemMap/generateCords/ToxAnalysisGlobal/Desc/SDF/16120-70-0.sdf")
print d
d = descriptors3D.get3Ddesc("/home/aborrel/ChemMap/generateCords/ToxAnalysisGlobal/Desc/SDF/74050-98-9.sdf")
print d








###################################
#drugbank https://www.drugbank.ca/#
###################################

psdf = "/home/aborrel/ChemMap/generateCords/drugbank-20-12-2017.sdf"
pranalysis = "/home/aborrel/ChemMap/generateCords/drugBankAnalysis/"
kname = "DATABASE_ID"

#main(psdf, pranalysis, kname, Desc1D2D=0, generation3D=1, Desc3D=1)


###################################
# analysis and compute descriptor #
###################################

#db = loadDB.sdfDB(psdf, pranalysis)
#db.writeTable()
#db.writeTableSpecific(["DRUG_GROUPS", "GENERIC_NAME", "FORMULA", "MOLECULAR_WEIGHT", "ALOGPS_LOGP", "SYNONYMS"], "tableChemmaps")
#db.writeNameforJS(prforjsupdate + "listSearch.js")




# descriptor computation #
##########################
#prDesc = "/home/aborrel/ChemMap/updateDrugBank/dbanalyse/Desc/"
#pathFolder.createFolder(prDesc, clean = 1)

#dpfiledesc = computeDesc(psdf, prDesc)
#analysis.PCAplot(dpfiledesc["1D2D"], dpfiledesc["3D"])




#####################
# update webserver  #
#####################





#pcord3D = pranalysis + "cor3D.csv"

#extractCloseCompounds(pcord3D, 20, pranalysis + "lneighbor.js")



#formate database and descriptor computation




