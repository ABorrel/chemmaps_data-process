import pathFolder
import loadDB
import computeDB
import runExternalSoft
import toolbox
import createJS

from os import path, remove, listdir
import math

from formatDatabase import drugbank
from rdkit import Chem



def Run (psdf, pranalysis, kname, corval=0.8, maxquantile=80, control=1, Desc1D2D=0, generation3D = 1, Desc3D=0, projection=1, map=1, drawPNG=1):

    pathFolder.createFolder(pranalysis)

    ###################################
    # analysis and compute descriptor #
    ###################################

    # load db from drugbank
    db = loadDB.sdfDB(psdf, kname, pranalysis)
    db.parseAll()


    # descriptor computation #
    ##########################
    prDesc = pranalysis + "Desc/"
    pathFolder.createFolder(prDesc, clean=0)
    dpfiledesc = computeDB.computeDesc(psdf, prDesc, control=control, Desc1D2D=Desc1D2D, generation3D=generation3D, Desc3D=Desc3D, namek=kname)


    # do not think it usefull
    #writeDescMatrix("1D2D", prDesc + "1D2DbyChem/", prDesc + "3DbyChem/", prDesc + "SMIclean/", prDesc)
    #writeDescMatrix("3D", prDesc + "3DbyChem/", prDesc + "3DbyChem/", prDesc + "SMIclean/", prDesc)

    # analyse projection  and compute coordinate #
    ##############################################
    prproject = pathFolder.createFolder(pranalysis + "projection" + str(corval) + "-" + str(maxquantile) + "/")
    if projection == 1:
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


    ###################
    # map file update #
    ###################
    if map == 1:

        # 1. create map file #
        ######################

        prmap = pathFolder.createFolder(pranalysis + "map_" + str(corval) + "-" + str(maxquantile) + "/")
        runExternalSoft.RComputeMapFiles(dpfiledesc["1D2D"], dpfiledesc["3D"], prmap, corval, maxquantile)

        pcoordDim1Dim2 = prmap + "coord1D2D.csv"
        pcoordDim3Dim4 = prmap + "coord3D.csv"

        if not path.exists(pcoordDim1Dim2) or not path.exists(pcoordDim3Dim4):
            print "ERROR file map"
            return

        # 2. write properties #
        #######################
        ldesc = ["SMILES", "inchikey", 'JCHEM_ROTATABLE_BOND_COUNT', 'JCHEM_POLAR_SURFACE_AREA', 'MOLECULAR_WEIGHT',
                 'JCHEM_PHYSIOLOGICAL_CHARGE', 'SYNONYMS', 'JCHEM_RULE_OF_FIVE', 'JCHEM_VEBER_RULE', 'FORMULA',
                 'JCHEM_GHOSE_FILTER', 'GENERIC_NAME', 'JCHEM_TRADITIONAL_IUPAC', 'ALOGPS_SOLUBILITY',
                 'JCHEM_MDDR_LIKE_RULE', 'SECONDARY_ACCESSION_NUMBERS', 'DRUG_GROUPS', 'PRODUCTS',
                 'JCHEM_IUPAC', 'ALOGPS_LOGP', 'ALOGPS_LOGS', 'JCHEM_PKA_STRONGEST_BASIC',
                 'JCHEM_NUMBER_OF_RINGS', 'JCHEM_PKA', 'JCHEM_ACCEPTOR_COUNT',
                 'JCHEM_PKA_STRONGEST_ACIDIC', 'ORIGINAL_SOURCE', 'EXACT_MASS', 'JCHEM_DONOR_COUNT',
                 'INTERNATIONAL_BRANDS', 'JCHEM_AVERAGE_POLARIZABILITY', 'JCHEM_BIOAVAILABILITY',
                 'JCHEM_REFRACTIVITY', 'JCHEM_LOGP', 'JCHEM_FORMAL_CHARGE', 'SALTS']


        dinfo = {}
        ddesc = toolbox.loadMatrixToDict(dpfiledesc["1D2D"])
        for chemID in ddesc.keys():
            dinfo[chemID] = {}
            for desc in ldesc:
                dinfo[chemID][desc] = "NA"

            dinfo[chemID]["SMILES"] = ddesc[chemID]["SMILES"]

            # generate inchkey
            chemMol = Chem.MolFromSmiles(dinfo[chemID]["SMILES"])
            inchi = Chem.inchi.MolToInchi(chemMol)
            inchikey = Chem.inchi.InchiToInchiKey(inchi)
            dinfo[chemID]["inchikey"] = inchikey



            i = 0
            imax = len(db.lc)
            while i < imax:
                if db.lc[i]["DATABASE_ID"] == chemID:
                    for k in db.lc[i].keys():
                        if k in ldesc and k != "SMILES":
                            dinfo[chemID][k] = db.lc[i][k]

                    del db.lc[i]
                    break
                i = i + 1

        pfilout = prmap + "TableProp.csv"
        filout = open(pfilout, "w")
        filout.write("ID\t%s\n" % ("\t".join(ldesc)))
        for chemID in dinfo.keys():
            filout.write("%s\t%s\n" % (chemID, "\t".join([dinfo[chemID][k] for k in ldesc])))
        filout.close()



        # 3. update JS properties #
        ###########################
        #createJS.formatInfo(db, dpfiledesc["1D2D"], lkinfo, pfileDataJS, pranalysis)


        # 4. update neighborhood #
        ##########################
        #createJS.extractCloseCompounds(pcoordsCombine, 20, pfileDataJS, pranalysis)







def writeDescMatrix(typeDesc, pr2D, pr3D, prSMI, prout):

    # write 2D
    if typeDesc == "1D2D" or typeDesc == "1D" or typeDesc == "2D":
        prin = pr2D
        pfilout = prout + "1D2D.csv"
    elif typeDesc == "3D":
        prin = pr3D
        pfilout = prout + "3D.csv"

    if path.exists(pfilout):
        return pfilout

    filout = open(pfilout, "w")

    ldesc = []
    lfSMI = listdir(prSMI)
    for fSMI in lfSMI:
        DBID = fSMI[0:-4]
        pSMI = prSMI + fSMI
        fSMILES = open(pSMI, "r")
        SMILES = fSMILES.readlines()[0].strip()

        pdesc = prin + fSMI[0:-3] + "txt"
        if path.exists(pdesc):
            dink = toolbox.loadMatrixToDict(pdesc, sep="\t")

            if ldesc == []:
                ldesc = dink[DBID].keys()
                del ldesc[ldesc.index("ID")]
                filout.write("ID\tSMILES\t" + "\t".join(ldesc) + "\n")

            # put back DSSTOX
            dink[DBID]["ID"] = DBID
            dink[DBID]["SMILES"] = SMILES
            filout.write("%s\t%s"%(dink[DBID]["ID"],  dink[DBID]["SMILES"]))
            for desc in ldesc:
                try: filout.write("\t%s"%(str(dink[DBID][desc])))
                except: filout.write("\tNA")
            filout.write("\n")
    filout.close()