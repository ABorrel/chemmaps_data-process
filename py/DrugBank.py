import pathFolder
import parseSDF
import toolbox
import DBrequest

#=> to import ToxCast librairy
import sys
#sys.path.insert(0, "/home/borrela2/development/descriptor/")
sys.path.insert(0, "C:\\Users\\borrela2\\development\\molecular-descriptors\\")
import Chemical

from os import path, remove
from random import shuffle


LPROP = ['JCHEM_ROTATABLE_BOND_COUNT', 'JCHEM_POLAR_SURFACE_AREA', 'MOLECULAR_WEIGHT',
                 'JCHEM_PHYSIOLOGICAL_CHARGE', 'SYNONYMS', 'JCHEM_RULE_OF_FIVE', 'JCHEM_VEBER_RULE', 'FORMULA',
                 'JCHEM_GHOSE_FILTER', 'GENERIC_NAME', 'JCHEM_TRADITIONAL_IUPAC', 'ALOGPS_SOLUBILITY',
                 'JCHEM_MDDR_LIKE_RULE', 'SECONDARY_ACCESSION_NUMBERS', 'DRUG_GROUPS', 'PRODUCTS',
                 'JCHEM_IUPAC', 'ALOGPS_LOGP', 'ALOGPS_LOGS', 'JCHEM_PKA_STRONGEST_BASIC',
                 'JCHEM_NUMBER_OF_RINGS', 'JCHEM_PKA', 'JCHEM_ACCEPTOR_COUNT',
                 'JCHEM_PKA_STRONGEST_ACIDIC', 'ORIGINAL_SOURCE', 'EXACT_MASS', 'JCHEM_DONOR_COUNT',
                 'INTERNATIONAL_BRANDS', 'JCHEM_AVERAGE_POLARIZABILITY', 'JCHEM_BIOAVAILABILITY',
                 'JCHEM_REFRACTIVITY', 'JCHEM_LOGP', 'JCHEM_FORMAL_CHARGE', 'SALTS']


class DrugBank:
    def __init__(self, pDBsdf, prDesc, prAnalysis):

        self.pSDF = pDBsdf
        self.prDesc = prDesc
        self.prout = prAnalysis


    def parseSDFDB (self):

        prForDB = pathFolder.createFolder(self.prout + "forDB/")

        cSDF = parseSDF.parseSDF(self.pSDF, "DATABASE_ID", self.prout)
        cSDF.parseAll()

        pfilout = prForDB + "db.csv"
        try:remove(pfilout)
        except:pass
        print(pfilout)
        if path.exists(pfilout):
            dchem = toolbox.loadMatrixToDict(pfilout)
        else:


            dchem = {}
            for chemSDF in cSDF.lc[:100]:
                drugbank_id = chemSDF["DATABASE_ID"]
                SMILES_origin = chemSDF["SMILES"]

                dchem[drugbank_id] = {}
                dchem[drugbank_id]["drugbank_id"] = drugbank_id
                dchem[drugbank_id]["smiles_origin"] = SMILES_origin

                # prepare ligand
                cchem = Chemical.Chemical(SMILES_origin, self.prDesc)
                cchem.prepChem()
                if cchem.err == 1:
                    qsar_ready = 0
                    cleanSMILES = "NA"
                    inchikey = "NA"
                else:
                    qsar_ready = 1
                    cleanSMILES = cchem.smi
                    inchikey = cchem.generateInchiKey()
                    cchem.writeSMIClean()

                dchem[drugbank_id]["smiles_clean"] = cleanSMILES
                dchem[drugbank_id]["inchikey"] = inchikey
                dchem[drugbank_id]["qsar_ready"] = qsar_ready
                db_prop = []
                for PROP in LPROP:
                    if PROP in list(chemSDF.keys()):
                        db_prop.append(chemSDF[PROP])
                    else:
                        db_prop.append("NA")
                dchem[drugbank_id]["DB_prop"] = db_prop

        # write table for control -> after open and put in the DB
        filout = open(pfilout, "w", encoding="utf8")
        filout.write("drugbank_id\tsmiles_origin\tsmiles_clean\tinchikey\tqsar_ready\tDB_prop\n")
        for chem in dchem.keys():
            filout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chem, dchem[chem]["smiles_origin"], dchem[chem]["smiles_clean"], dchem[chem]["inchikey"], dchem[chem]["qsar_ready"], "___".join(dchem[chem]["DB_prop"])))
        filout.close()
        self.dchem = dchem


    def pushChemInDB(self):

        if not "dchem" in self.__dict__:
            self.parseSDFDB()

        cDB = DBrequest.DBrequest()
        cDB.connOpen()
        for chem in self.dchem.keys():
            lprop = self.dchem[chem]["DB_prop"]
            wprop = "{" + ",".join([ "\"%s\"" % (prop.replace("\'", "").replace("\"", "")) for prop in lprop]) + "}"
            cDB.connOpen()
            cDB.addElement("drugbank", ["drugbank_id", "smiles_origin", "smiles_clean", "inchikey", "qsar_ready", "DB_prop"],
                                 [self.dchem[chem]["drugbank_id"], self.dchem[chem]["smiles_origin"], self.dchem[chem]["smiles_clean"],
                                  self.dchem[chem]["inchikey"], self.dchem[chem]["qsar_ready"], wprop])
            cDB.connClose()



    def computeDesc(self):

        if not "dchem" in self.__dict__:
            self.parseSDFDB()

        lchemID = list(self.dchem.keys())
        shuffle(lchemID)

        for chemID in lchemID:
            SMILESClean = self.dchem[chemID]["smiles_clean"]
            cChem = Chemical.Chemical(SMILESClean, self.prDesc)
            cChem.prepChem()
            if cChem.err == 0:
                cChem.computeAll2D(update=0)
                cChem.writeMatrix("2D")
                cChem.set3DChemical()
                cChem.computeAll3D(update=0)
                cChem.writeMatrix("3D")


    def computeCoords(self):

        return


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
        prSMI = prDesc + "SMIclean/"
        #print prSMI
        #ddd
        if path.exists(prSMI):
            lsmi = listdir(prSMI)
            shuffle(lsmi)
            print ("ddd")
            # control if nSDF = nPNG
            if len(lsmi) != len(listdir(prpng)):
                for smifile in lsmi:
                    fSMI = open(prSMI + smifile, "r")
                    SMI = fSMI.readlines()[0].strip()
                    fSMI.close()
                    inchikey = toolbox.convertSMILEStoINCHIKEY(SMI)

                    runExternalSoft.molconvert(prSMI + smifile, prpng + inchikey + ".png")
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
            print ("ERROR file map")
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
                ldesc = list(dink[DBID].keys())
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