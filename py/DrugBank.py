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
from re import search
from rdkit import Chem

class DrugBank:
    def __init__(self, pDBsdf, prDesc, prAnalysis):

        self.pSDF = pDBsdf
        self.prDesc = prDesc
        self.prout = prAnalysis
        self.c_DB = DBrequest.DBrequest()
        self.c_CompDesc = CompDesc.CompDesc("", self.prDesc)

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
                if d_chem["smiles_clean"] == "":
                    self.c_DB.addElement(name_table, ["smiles_origin", "drugbank_id", "name"],
                                    [d_chem["SMILES"], d_chem["DRUGBANK_ID"], d_chem["GENERIC_NAME"].replace("'", "\\'").replace("\u03b1", "alpha").replace("'", "''").replace("\u03b2", "beta").replace("\u03ba", "kappa").replace("\u03bb", "lamda").replace("\u2032", "`"),])

                else:
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

    def pushDB_desc_name(self):

        # check if not empty
        DBname = "chem_descriptor_1d2d_name"
        size_table = self.c_DB.execCMD("select count(*) from %s"%(DBname))
        size_table = int(size_table[0][0])
        #1D and 2D descriptor
        if size_table == 0:
            self.c_DB.connOpen()
            l_desc = self.c_CompDesc.getLdesc("2D")
            for desc in l_desc:
                self.c_DB.addElement(DBname, ["name"], [desc])
            self.c_DB.connClose()

        #for 3D desc
        DBname = "chem_descriptor_3d_name"
        size_table = self.c_DB.execCMD("select count(*) from %s"%(DBname))
        size_table = int(size_table[0][0])

        if size_table == 0:
            self.c_DB.connOpen()
            l_desc = self.c_CompDesc.getLdesc("3D")
            for desc in l_desc:
                self.c_DB.addElement(DBname, ["name"], [desc])
            self.c_DB.connClose()

    def compute_pushDB_Desc(self, table_chemicals="chemicals"):
        """
        !! value NA will be -9999 because of the database format
        """

        # first extract SMILES clean with ID from the chemicals database
        cmd_SQL = "SELECT id, smiles_clean, inchikey FROM chemicals WHERE smiles_clean IS NOT NULL AND drugbank_id IS NOT NULL"
        l_chem = self.c_DB.execCMD(cmd_SQL)
        
        l_chem_desc = "SELECT inchikey FROM chemical_description WHERE map_name = drugbank"
        l_chem = self.c_DB.execCMD(cmd_SQL)

        l_desc_2D = self.c_CompDesc.getLdesc("1D2D")
        l_desc_3D = self.c_CompDesc.getLdesc("3D")

        self.c_DB.connOpen()
        for chem in l_chem[99:]:
            SMILES = chem[1]
            inchikey = chem[2]

            # the SMILES is already prepared
            self.c_CompDesc = CompDesc.CompDesc(input=SMILES, prdesc=self.prDesc)
            self.c_CompDesc.smi = SMILES
            self.c_CompDesc.inchikey = inchikey
            self.c_CompDesc.mol = Chem.MolFromSmiles(SMILES)
            self.c_CompDesc.computeAll2D()
            l_compute_desc = [self.c_CompDesc.all2D[desc_2D] if self.c_CompDesc.all2D[desc_2D] != "NA" else "-9999" for desc_2D in l_desc_2D]
            w1D2D = "{" + ",".join(["\"%s\"" % (compute_desc) for compute_desc in l_compute_desc]) + "}"

            self.c_CompDesc.set3DChemical(psdf3D = "")
            if self.c_CompDesc.err == 0:
                self.c_CompDesc.computeAll3D()
                l_compute_desc = [self.c_CompDesc.all3D[desc_3D] if self.c_CompDesc.all3D[desc_3D] != "NA" else "-9999" for desc_3D in l_desc_3D]
                w3D = "{" + ",".join(["\"%s\"" % (compute_desc) for compute_desc in l_compute_desc]) + "}"
                self.c_DB.addElement("chemical_description", ["inchikey", "desc_1d2d", "desc_3D", "map_name"], [inchikey, w1D2D, w3D, "drugbank"])

            else:
                self.c_DB.addElement("chemical_description", ["inchikey", "desc_1d2d", "map_name"], [inchikey, w1D2D, "drugbank"])

        self.c_DB.connClose()
        

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


