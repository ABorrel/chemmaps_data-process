import pathFolder
import parseSDF
import toolbox
import DBrequest
import runExternalSoft

#=> to import ToxCast librairy
import sys
sys.path.insert(0, "/home/borrela2/development/descriptor/")
#sys.path.insert(0, "C:\\Users\\borrela2\\development\\molecular-descriptors\\")
import Chemical

from os import path, remove
from random import shuffle
from math import sqrt


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


    def pushDrugBankNamePropInDB(self):

        cDB = DBrequest.DBrequest()
        i = 1
        for PROP in LPROP:
            cDB.addElement("drugbank_prop", ["id", "property"], [i, PROP])
            i = i + 1



    def parseSDFDB (self):

        prForDB = pathFolder.createFolder(self.prout + "forDB/")

        cSDF = parseSDF.parseSDF(self.pSDF, "DATABASE_ID", self.prout)
        cSDF.parseAll()
        cSDF.writeTableSpecific(LPROP, "Table_prop")

        pfilout = prForDB + "db.csv"
        #try:remove(pfilout)
        #except:pass
        #print(pfilout)
        if path.exists(pfilout):
            dchem = toolbox.loadMatrixToDict(pfilout, sep = "\t")
            for chem in dchem.keys():
                dchem[chem]["DB_prop"] = dchem[chem]["DB_prop"].split("___")
        else:
            dchem = {}
            for chemSDF in cSDF.lc:
                if not "SMILES" in list(chemSDF.keys()):
                    continue
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
                        db_prop.append(chemSDF[PROP].replace("\t", ""))
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
        for chem in self.dchem.keys():
            cDB.addElement("drugbank_chem", ["db_id", "smiles_origin", "smiles_clean", "inchikey", "qsar_ready"],
                                 [self.dchem[chem]["drugbank_id"], self.dchem[chem]["smiles_origin"], self.dchem[chem]["smiles_clean"],
                                  self.dchem[chem]["inchikey"], self.dchem[chem]["qsar_ready"]])



    def pushPropInDB(self):

        if not "dchem" in self.__dict__:
            self.parseSDFDB()

        cDB = DBrequest.DBrequest()
        for chem in self.dchem.keys():
            lprop = self.dchem[chem]["DB_prop"]
            wprop = "{" + ",".join(["\"%s\"" % (prop.replace("\'", "").replace("\"", "")) for prop in lprop]) + "}"
            cDB.addElement("drugbank_prop", ["drugbank_id", "prop_value"], [self.dchem[chem]["drugbank_id"], wprop])

        return

    def pushPropNameInDB(self, typeDesc):
        ldesc = Chemical.getLdesc(typeDesc)
        if typeDesc == "1D2D" or typeDesc == "1D" or typeDesc == "2D":
            DBname = "desc_1d2d_prop"
        elif typeDesc == "3D":
            DBname = "desc_3d_prop"

        cDB = DBrequest.DBrequest()
        i = 1
        for desc in ldesc:
            cDB.addElement(DBname, ["id", "descriptor"], [i, desc])
            i = i + 1



    def computeDesc(self, insertDB=1, update=0):

        pfilout1D2D = self.prout + "1D2D.csv"
        pfilout3D = self.prout + "3D.csv"

        if update == 0 and path.exists(pfilout3D) and path.exists(pfilout1D2D) and insertDB == 0:
            self.p1D2D = pfilout1D2D
            self.p3D = pfilout3D
            return


        if not "dchem" in self.__dict__:
            self.parseSDFDB()

        lchemID = list(self.dchem.keys())
        shuffle(lchemID)

        ldesc1D2D = Chemical.getLdesc("1D2D")
        ldesc3D = Chemical.getLdesc("3D")

        filout1D2D = open(pfilout1D2D, "w")
        filout1D2D.write("inchikey\t" + "\t".join(ldesc1D2D) + "\n")

        filout3D = open(pfilout3D, "w")
        filout3D.write("inchikey\t" + "\t".join(ldesc3D) + "\n")

        if insertDB == 1:
            cDB = DBrequest.DBrequest()

        for chemID in lchemID:
            SMILESClean = self.dchem[chemID]["smiles_clean"]
            if SMILESClean == "NA":
                continue
            cChem = Chemical.Chemical(SMILESClean, self.prDesc)
            cChem.prepChem()

            #print(SMILESClean)
            if cChem.err == 0 :
                # 2D descriptors
                cChem.computeAll2D(update=0)
                cChem.writeMatrix("2D")
                filout1D2D.write(
                    "%s\t%s\n" % (cChem.inchikey, "\t".join([str(cChem.all2D[desc]) for desc in ldesc1D2D])))

                # insert in DB
                if insertDB == 1:
                    cDB.verbose = 0
                    out1D2D = cDB.getRow("desc_1d2d", "inchikey='%s'" % (cChem.inchikey))
                    if out1D2D == []:
                        valDesc = [cChem.all2D[desc1D2D] for desc1D2D in ldesc1D2D]
                        valDesc = ['-9999' if desc == "NA" else desc for desc in valDesc]

                        w1D2D = "{" + ",".join(["\"%s\"" % (desc) for desc in valDesc]) + "}"
                        cDB.addElement("desc_1d2d", ["inchikey", "desc_value"], [cChem.inchikey, w1D2D])


                # 3D descriptors
                cChem.set3DChemical()
                # control if 3D generated
                if cChem.err == 0:
                    cChem.computeAll3D(update=0)
                    cChem.writeMatrix("3D")

                    # write master table
                    filout3D.write("%s\t%s\n" % (cChem.inchikey, "\t".join([str(cChem.all3D[desc]) for desc in ldesc3D])))

                    # put in table descriptors
                    if insertDB == 1:
                        out3D = cDB.getRow("desc_3d", "inchikey='%s'" % (cChem.inchikey))
                        if out3D == []:
                            valDesc = [cChem.all3D[desc3D] for desc3D in ldesc3D]
                            valDesc = ['-9999' if desc == "NA" else desc for desc in valDesc]

                            w3D = "{" + ",".join(["\"%s\"" % (desc) for desc in valDesc]) + "}"
                            cDB.addElement("desc_3d", ["inchikey", "desc_value"], [cChem.inchikey, w3D])

        filout1D2D.close()
        filout3D.close()

        self.p1D2D = pfilout1D2D
        self.p3D = pfilout3D



    def computeCoords(self, corVal, distributionVal, insertDB=1):


        if not "p1D2D" in self.__dict__ and not "p3D" in self.__dict__:
            self.computeDesc(insertDB=0)

        err = 0
        # create coords
        prmap = pathFolder.createFolder(self.prout + "map_" + str(corVal) + "-" + str(distributionVal) + "/")
        pcoordDim1Dim2 = prmap + "coord1D2D.csv"
        pcoordDim3D = prmap + "coord3D.csv"
        if path.exists(pcoordDim1Dim2) and path.exists(pcoordDim3D):
            self.pcoords1D2D = pcoordDim1Dim2
            self.pcoords3D = pcoordDim3D

        elif not path.exists(pcoordDim1Dim2) or not path.exists(pcoordDim3D):
            runExternalSoft.RComputeMapFiles(self.p1D2D, self.p3D, prmap, corVal, distributionVal)
        elif not path.exists(pcoordDim1Dim2) or not path.exists(pcoordDim3D):
            print("ERROR file map")
            err = 1

        self.pcoords1D2D = pcoordDim1Dim2
        self.pcoords3D = pcoordDim3D


        if insertDB == 1 and err ==0:
            dcoord1D2D = toolbox.loadMatrixToDict(pcoordDim1Dim2, sep = ",")
            dcoord3D = toolbox.loadMatrixToDict(pcoordDim3D, sep = ",")
            cDB = DBrequest.DBrequest()
            cDB.verbose = 0
            for chem in dcoord1D2D.keys():
                #print(chem)
                #out1D2D = cDB.getRow("drugbank_coords", "inchikey='%s'" % (chem))
                out1D2D = []
                if out1D2D == []:
                    nbdim1d2d = len(dcoord1D2D[chem].keys()) - 1
                    nbdim3d = len(dcoord3D[chem].keys()) - 1

                    w1D2D = "{" + ",".join(["\"%s\"" % (dcoord1D2D[chem]["DIM" + str(i)]) for i in range(1, nbdim1d2d + 1)]) + "}"
                    w3D = "{" + ",".join(["\"%s\"" % (dcoord3D[chem]["DIM3-" + str(i)]) for i in range(1, nbdim3d + 1)]) + "}"
                    cDB.addElement("drugbank_coords", ["inchikey", "dim1d2d", "dim3d", "indrugbank"], [chem, w1D2D, w3D, "True"])


    def runRprojection(self, corVal, distributionVal):

        if not "p1D2D" in self.__dict__ and not "p3D" in self.__dict__:
            self.computeDesc(insertDB=0)

        # create coords
        prproj = pathFolder.createFolder(self.prout + "proj_" + str(corVal) + "-" + str(distributionVal) + "/")
        runExternalSoft.RComputeCor(self.p1D2D, self.p3D, prproj, corVal, distributionVal)


    def neighbormatrix(self, nbNeighbor, lnDim):

        if not "pcoords1D2D" in self.__dict__:
            print("Compute Coord first")
            return 1
        else:
            prNeighbor = pathFolder.createFolder(self.prout + "Neighbors/")
            pfilout = prNeighbor + "Table_DIM1D2D-" + str(lnDim[0]) + "_" + str(lnDim[1]) + ".csv"
            if path.exists(pfilout):
                ddist = toolbox.loadMatrixToDictp(pfilout)
                for chem in ddist.keys():
                    ddist[chem] = ddist[chem].split(" ")
            else:
                dDim1D2D = toolbox.loadMatrixToDict(self.pcoords1D2D, sep = ",")
                dDim3D = toolbox.loadMatrixToDict(self.pcoords3D, sep=",")

                # compute all dimension
                if lnDim ==[]:
                    chem1 = list(dDim1D2D.keys())[0]
                    n1D2D = len(list(dDim1D2D[chem1].keys())) - 1
                    n3D = len(list(dDim3D[chem1].keys())) - 1
                    lnDim = [n1D2D, n3D]


                dcor = {}
                # from 1D2D coord
                for inch in dDim1D2D.keys():
                    dcor[inch] = []

                    i = 1
                    while i <= lnDim[0]:
                        dcor[inch].append(float(dDim1D2D[inch]["DIM" + str(i)]))
                        i = i + 1

                    i = 1
                    while i <= lnDim[1]:
                        dcor[inch].append(float(dDim3D[inch]["DIM3-" + str(i)]))
                        i = i + 1

                ddist = {}
                for ID in dcor.keys():
                    ddist[ID] = {}
                    for ID2 in dcor.keys():
                        if ID != ID2:
                            ddist[ID][ID2] = sqrt(sum([(xi - yi) ** 2 for xi, yi in zip(dcor[ID], dcor[ID2])]))

                    lID = [i[0] for i in sorted(ddist[ID].items(), key=lambda x: x[1])][:nbNeighbor]
                    ddist[ID] = lID

                # write in table
                ftable = open(pfilout, "w")
                ftable.write("ID\tNeighbors\n")
                for ID in ddist.keys():
                    ftable.write("%s\t%s\n" % (ID, " ".join(ddist[ID])))
                ftable.close()

    def pushNeighbors(self):
        prneighbor = pathFolder.createFolder(self.prout + "Neighbors/")
        ptable3Dim = prneighbor + "Table_DIM1D2D-2_1.csv"
        ptableNDim = prneighbor + "Table_DIM1D2D-120_164.csv"
        if path.exists(ptable3Dim) and path.exists(ptableNDim):
            ddist3D = toolbox.loadMatrixToDict(ptable3Dim)
            for chem in ddist3D.keys():
                ddist3D[chem] = ddist3D[chem]["Neighbors"].split(" ")
            ddistND = toolbox.loadMatrixToDict(ptableNDim)
            for chem in ddistND.keys():
                ddistND[chem] = ddistND[chem]["Neighbors"].split(" ")

            cDB = DBrequest.DBrequest()
            cDB.verbose = 0
            for chem in ddist3D.keys():
                # print(chem)
                out1D2D = cDB.getRow("drugbank_neighbors", "inchikey='%s'" % (chem))
                if out1D2D == []:
                    w3D = "{" + ",".join(["\"%s\"" % (neighbor) for neighbor in ddist3D[chem]]) + "}"
                    wND = "{" + ",".join(["\"%s\"" % (neighbor) for neighbor in ddistND[chem]]) + "}"
                    cDB.addElement("drugbank_neighbors", ["inchikey", "neighbors_dim3", "neighbors_dimn"], [chem, w3D, wND])


