import pathFolder
import parseSDF
import toolbox
import DBrequest
import runExternalSoft

#=> to import ToxCast librairy
import sys
sys.path.insert(0, "/home/borrela2/development/descriptor/")
#sys.path.insert(0, "C:\\Users\\borrela2\\development\\molecular-descriptors\\")
import CompDesc

from os import path, remove
from random import shuffle
from math import sqrt


class DrugBank:
    def __init__(self, pDBsdf, prDesc, prAnalysis):

        self.pSDF = pDBsdf
        self.prDesc = prDesc
        self.prout = prAnalysis
        self.c_DB = DBrequest.DBrequest()

    def pushDB_prop_name(self, name_table = "chem_prop_drugbank_name"):

        # check if not empty
        size_table = self.c_DB.execCMD("select count(*) from %s"%(name_table))
        size_table = int(size_table[0][0])
   
        if size_table == 0:
            for prop in self.l_prop:
                self.c_DB.addElement(name_table, ["name"], [prop])


        l_prop = self.c_DB.execCMD("select * from %s"%(name_table))
        l_prop = [ prop[1] for prop in l_prop]
        print(l_prop)
        
        self.l_prop = l_prop

    def parse_SDF (self):

        cSDF = parseSDF.parseSDF(self.pSDF, "DATABASE_ID", self.prout)
        cSDF.parseAll()
        cSDF.write_table_prop("Table_prop")
        
        self.l_prop = cSDF.l_prop
        self.l_chem_sdf = cSDF.lc

    def prep_chem(self):
        
        # can load the file
        p_filout = self.prout + "chem.csv"
        
        for d_chem in self.l_chem_sdf:
            if not "SMILES" in list(d_chem.keys()):
                continue
            else:
                SMILES_origin = d_chem["SMILES"]

                cchem = CompDesc.CompDesc(SMILES_origin, self.prDesc)
                cchem.prepChem()
                if cchem.err == 1:
                    cleanSMILES = ""
                    inchikey = ""
                else:
                    cleanSMILES = cchem.smi
                    cchem.generateInchiKey()
                    if cchem.err == 0:
                        inchikey = cchem.inchikey

                d_chem["smiles_clean"] = cleanSMILES
                d_chem["inchikey"] = inchikey          
            

        # write table for control -> after open and put in the DB
        filout = open(p_filout, "w", encoding="utf8")
        filout.write("drugbank_id\tsmiles_origin\tsmiles_clean\tinchikey\tname\n")
        for d_chem in self.l_chem_sdf:
            filout.write("%s\t%s\t%s\t%s\t%s\n" % (d_chem["DRUGBANK_ID"], d_chem["SMILES"], d_chem["smiles_clean"], d_chem["inchikey"], d_chem["GENERIC_NAME"]))
        filout.close()

    def pushDB_chem(self, name_table="chemicals"):

        # first check if chemicals are already in the DB
        # check if not empty
        size_table = self.c_DB.execCMD("SELECT count(*) FROM %s WHERE drugbank_id IS NOT NULL"%(name_table))
        size_table = int(size_table[0][0])
        print(size_table)

        # case empty --> build
        if size_table == 0:
            self.prep_chem()
            self.c_DB.connOpen()
            for d_chem in self.l_chem_sdf:
                self.c_DB.addElement(name_table, ["smiles_origin", "smiles_clean", "inchikey", "drugbank_id", "name"],
                                    [d_chem["SMILES"], d_chem["smiles_clean"],
                                    d_chem["inchikey"], d_chem["DRUGBANK_ID"], d_chem["GENERIC_NAME"].replace("'", "\\'").replace("\u03b1", "alpha").replace("'", "''").replace("\u03b2", "beta").replace("\u03ba", "kappa").replace("\u03bb", "lamda").replace("\u2032", "`"),])
            self.c_DB.connClose


        l_chem_out = []
        l_chem = self.c_DB.execCMD("select * from %s WHERE drugbank_id IS NOT NULL"%(name_table))
        for tup_chem in l_chem:
            d_chem = {"id":tup_chem[0], "smiles_clean": tup_chem[2]}
            l_chem_out.append(d_chem)
        
        self.l_chem_prep = l_chem_out


    def pushDB_prop(self):

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


