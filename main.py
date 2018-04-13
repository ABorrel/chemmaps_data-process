import pathFolder
import loadDB
import computeDB
import runExternalSoft
import toolbox
import createJS

from os import path, remove, listdir
import math

from formatDatabase import drugbank






###############
##   MAIN    ##
###############


def main(psdf, pranalysis, kname, lkinfo=[], corval=0.8, maxquantile=80, control=1, Desc1D2D=1, generation3D = 1, Desc3D=1, projection=1, JS=1, drawPNG=1):

    pathFolder.createFolder(pranalysis)

    ###################################
    # analysis and compute descriptor #
    ###################################

    db = loadDB.sdfDB(psdf, kname, pranalysis)
    db.writeTable("TaleProp.csv")


    # descriptor computation #
    ##########################
    prDesc = pranalysis + "Desc/"
    pathFolder.createFolder(prDesc, clean=0)
    dpfiledesc = computeDB.computeDesc(psdf, prDesc, control=control, Desc1D2D=Desc1D2D, generation3D=generation3D, Desc3D=Desc3D, namek=kname)
    # analyse projection  and compute coordinate #
    ##############################################
    prproject = pathFolder.createFolder(pranalysis + "projection" + str(corval) + "-" + str(maxquantile) + "/")
    if projection == 1:
        print dpfiledesc["1D2D"], dpfiledesc["3D"], prproject, corval, maxquantile
        runExternalSoft.RComputeCor(dpfiledesc["1D2D"], dpfiledesc["3D"], prproject, corval, maxquantile)

    ###################
    # for the website #
    ###################


    # 1. compute png #
    ##################
    if drawPNG == 1:
        prpng = pranalysis + "cpdpng/"
        pathFolder.createFolder(prpng)

        # draw from desc descriptors
        prSMI = prDesc + "SMI/"
        if path.exists(prSMI):
            lsmi = listdir(prSMI)
            # control if nSDF = nPNG
            if len(lsmi) != len(listdir(prpng)):
                for smifile in lsmi:
                    runExternalSoft.molconvert(prSMI + smifile, prpng + smifile.split(".")[0] + ".png")
        else:
            # case where considering only original map
            db.drawMolecules(prpng)


    ##################
    # JS file update #
    ##################
    if JS == 1:
        prforjsupdate = pranalysis + "JS" + str(corval) + "-" + str(maxquantile) + "/"
        pathFolder.createFolder(prforjsupdate)
        pfileDataJS = prforjsupdate + "data.js"

        if path.exists(pfileDataJS):
            remove(pfileDataJS)

        # 2. update JS coords #
        #######################
        pcoordsCombine = prproject + "coordPCAcombined.csv"
        if path.exists(pcoordsCombine):
            createJS.formatCoordinates(pcoordsCombine, pfileDataJS)


        # 3. update JS properties #
        ###########################
        createJS.formatInfo(db, dpfiledesc["1D2D"], lkinfo, pfileDataJS)


        # 4. update neighborhood #
        ##########################
        createJS.extractCloseCompounds(pcoordsCombine, 20, pfileDataJS)


#########################################
# Tox train dataset from Kamel Monsouri #
#########################################

#psdf = "/home/aborrel/ChemMap/generateCords/ToxTrainingSet_3D.sdf"
#pranalysis = "/home/aborrel/ChemMap/generateCords/ToxAnalysis/"
#kname = "CASRN"

#main(psdf, pranalysis, kname, Desc1D2D=1, generation3D=0, Desc3D=1)

#import descriptors3D
#print descriptors3D.get3Ddesc("/home/aborrel/ChemMap/generateCords/ToxAnalysisGlobal/Desc/SDF/361442-04-8.sdf")


##########################################
# Tox global dataset from Kamel Monsouri #
##########################################

#psdf = "/home/borrela2/ChemMaps/data_analysis/ToxTrainTest_3D.sdf"
#pranalysis = "/home/borrela2/ChemMaps/data_analysis/ToxAnalysisGlobal/"
#kname = "CASRN"
#lkinfo = ["CASRN", "LD50_mgkg", "GHS_category", "EPA_category", "Weight", "LogP"]

#corval = 0.9
#maxquantile = 90
#main(psdf, pranalysis, kname, lkinfo=lkinfo, corval=corval, maxquantile=maxquantile, control=1, Desc1D2D=0, generation3D=0, Desc3D=0, drawPNG=0, projection=1, JS=0)




###################################
#drugbank https://www.drugbank.ca/#
###################################

psdf = "/home/borrela2/ChemMaps/data_analysis/drugbank-20-12-2017.sdf"
pranalysis = "/home/borrela2/ChemMaps/data_analysis/drugBankAnalysis/"
kname = "DATABASE_ID"
lkinfo = ["DRUG_GROUPS", "GENERIC_NAME", "FORMULA", "MOLECULAR_WEIGHT", "ALOGPS_LOGP", "SYNONYMS"]

corval = 0.9
maxquantile = 90
main(psdf, pranalysis, kname, corval=corval, maxquantile=maxquantile, lkinfo=lkinfo, Desc1D2D=0, generation3D=0,
     Desc3D=1, projection=1, JS=1)




##############################
# by list from EPA dashboard #
##############################

prListChemical = "/home/borrela2/ChemMaps/data/toxEPA-lists/"
llistChem = listdir(prListChemical)

kname = "<CASRN>"
ksmile = "SMILES"
corval = 0.9
maxquantile = 90


for listChem in llistChem:
    plist = prListChemical + listChem
    print plist
    # main(psdf, pranalysis, kname, corval=corval, maxquantile=maxquantile, lkinfo=lkinfo, Desc1D2D=0, generation3D=0, Desc3D=0, projection=1, JS=1)


    # have to had generation 3D
    # acess





