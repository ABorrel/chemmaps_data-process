import pathFolder
import toolbox
import runExternalSoft
import DBrequest
import calculate
import parseSDF

#=> to import ToxCast librairy
import sys
sys.path.insert(0, "/home/borrela2/development/descriptor/")
#sys.path.insert(0, "C:\\Users\\borrela2\\development\\molecular-descriptors\\")
import CompDesc



from os import path, listdir, remove
from copy import deepcopy
from math import sqrt
from rdkit import Chem
from re import search
from random import shuffle
import time


LPROP = ["inchikey", "SMILES", "preferred_name", "GHS_category", "EPA_category", "consensus_LD50", "LD50_mgkg", "MolWeight", "LogOH_pred", "CATMoS_VT_pred", "CATMoS_NT_pred", "CATMoS_EPA_pred",
                 "CATMoS_GHS_pred", "CATMoS_LD50_pred", "CERAPP_Ago_pred", "CERAPP_Anta_pred", "CERAPP_Bind_pred", "Clint_pred", "CoMPARA_Ago_pred", "CoMPARA_Anta_pred",
                 "CoMPARA_Bind_pred", "FUB_pred", "LogHL_pred", "LogKM_pred", "LogKOA_pred", "LogKoc_pred", "LogBCF_pred", "LogD55_pred", "LogP_pred", "MP_pred", "pKa_a_pred",
                 "pKa_b_pred", "ReadyBiodeg_pred", "RT_pred", "LogVP_pred", "LogWS_pred", "BioDeg_LogHalfLife_pred", "BP_pred", "nbLipinskiFailures"]



class DSSTOX:

    def __init__(self, p_listChem, p_mapping, name_map, istart, iend, p_dir_Desc, p_dir_out):

        self.p_listChem = p_listChem
        self.p_mapping = p_mapping
        self.name_map = name_map
        self.istart = istart
        if iend == 0:
            self.iend = -1
        else:
            self.iend = iend
        self.plog = p_dir_out + "log.txt"
        self.p_dir_Desc = p_dir_Desc
        self.p_dir_out = p_dir_out
        self.c_DB = DBrequest.DBrequest()
        self.c_CompDesc = CompDesc.CompDesc("", self.p_dir_Desc)

    def load_annotation(self):

        #load 
        l_d_chemID = toolbox.loadMatrixToList(self.p_mapping, sep = ",")
        d_structure = toolbox.loadMatrixToDict(self.p_listChem, sep = ",")

        # load ID
        for d_chemID in l_d_chemID:
            dtxcid = d_chemID["DSSTOX_COMPOUND_ID"]
            d_chemID["SMILES"] = d_structure[dtxcid]["Original_SMILES"]
        self.l_chem = l_d_chemID[self.istart:self.iend]
        
    def prep_chem(self):
        
        if not "l_chem" in self.__dict__:
            self.load_annotation()
        
        p_filout = self.p_dir_out + "chemical_table.csv"

        for d_chem in self.l_chem:
            if not "SMILES" in list(d_chem.keys()):
                continue
            else:
                SMILES_origin = d_chem["SMILES"]

                cchem = CompDesc.CompDesc(SMILES_origin, self.p_dir_Desc)
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
        filout.write("dsstox_id\tsmiles_origin\tsmiles_clean\tinchikey\tname\n")
        for d_chem in self.l_chem:
            filout.write("%s\t%s\t%s\t%s\t%s\n" % (d_chem["DSSTOX_SUBSTANCE_ID"], d_chem["SMILES"], d_chem["smiles_clean"], d_chem["inchikey"], d_chem["PREFERRED_NAME"]))
        filout.close() 

    def pushDB_chem(self, name_table="chemicals"):

        # check if not empty
        size_table = self.c_DB.execCMD("SELECT count(*) FROM %s WHERE dsstox_id IS NOT NULL"%(name_table))
        print(size_table)
        size_table = int(size_table[0][0])
        print(size_table)
        if size_table == 0:
            # load all chemicals already in the database -> just case where not dsstox included
            l_smiles = self.c_DB.execCMD("SELECT smiles_origin FROM %s WHERE dsstox_id IS NULL"%(name_table))
            l_smiles = [smiles[0] for smiles in l_smiles]
            l_smiles = list(set(l_smiles))
            l_smiles.sort()

            l_dsstox = self.c_DB.execCMD("SELECT dsstox_id FROM %s WHERE dsstox_id IS NOT NULL"%(name_table))
            l_dsstox = [dsstox[0] for dsstox in l_dsstox]
            l_dsstox = list(set(l_dsstox))
            l_dsstox.sort()

            # case empty --> build all
            self.prep_chem()
            i = self.istart
            for d_chem in self.l_chem:
                i = i + 1
                if i%10000 == 0:
                    print(i)
                # need to check is the chem is not already in the DB
                if toolbox.binary_search(l_smiles, d_chem["SMILES"]) != -1:
                    # need to update database
                    cmd_update = "UPDATE %s SET dsstox_id='%s' where smiles_origin='%s'"%(name_table, d_chem["DSSTOX_SUBSTANCE_ID"], d_chem["SMILES"])
                    self.c_DB.execCMD(cmd_update)
                    continue

                # case id already in
                if toolbox.binary_search(l_dsstox, d_chem["DSSTOX_SUBSTANCE_ID"]) != -1:
                    continue

                if d_chem["smiles_clean"] == "":
                    self.c_DB.addElement(name_table, ["smiles_origin", "dsstox_id", "name"],
                                    [d_chem["SMILES"], d_chem["DSSTOX_SUBSTANCE_ID"], toolbox.update_str4DB(d_chem["PREFERRED_NAME"])])

                else:
                    self.c_DB.addElement(name_table, ["smiles_origin", "smiles_clean", "inchikey", "dsstox_id", "name"],
                                    [d_chem["SMILES"], d_chem["smiles_clean"],
                                    d_chem["inchikey"], d_chem["DSSTOX_SUBSTANCE_ID"], toolbox.update_str4DB(d_chem["PREFERRED_NAME"])])


        # load chemical in the database
        l_chem_out = []
        l_chem = self.c_DB.execCMD("select * from %s WHERE dsstox_id IS NOT NULL"%(name_table))
        for tup_chem in l_chem:
            d_chem = {"id":tup_chem[0], "smiles_clean": tup_chem[2]}
            l_chem_out.append(d_chem)
        
        self.l_chem_prep = l_chem_out

    def compute_pushDB_Desc(self):
        """
        !! value NA will be -9999 because of the database format
        """

        # first extract SMILES clean with ID from the chemicals database
        # extract all
        #cmd_SQL = "SELECT id, smiles_clean, inchikey FROM chemicals WHERE smiles_clean IS NOT NULL AND drugbank_id IS NOT NULL"
        
        #extract only missing value ++++
        # extract with exclusion duplicate
        # cmd = "select distinct c1.smiles_clean, c1.inchikey from chemicals c1 where c1.smiles_clean is not null except (select c.smiles_clean, cd.inchikey from chemical_description cd inner join chemicals c on c.inchikey = cd.inchikey where cd.map_name='dsstox')"
        cmd_SQL = "select smiles_clean, inchikey from chemicals where dsstox_id is not null and smiles_clean is not null and inchikey is not null"
        l_chem = self.c_DB.execCMD(cmd_SQL)

        # sample list
        l_chem = l_chem[self.istart:self.iend]
        shuffle(l_chem)
        
        # need to check if inchkey is already in database
        cmd_inch = "select inchikey from chemical_description where map_name = 'dsstox'"
        l_inch = self.c_DB.execCMD(cmd_inch)

        l_inch = [inch[0] for inch in l_inch]
        l_inch = list(set(l_inch))
        l_inch.sort()

        if len(l_chem) == 0:
            return

        l_desc_2D = self.c_CompDesc.getLdesc("1D2D")
        l_desc_3D = self.c_CompDesc.getLdesc("3D")

        l_chem_to_compute = [chem for chem in l_chem if toolbox.binary_search(l_inch, chem[1]) == -1]
        print("To compute:", len(l_chem_to_compute))

        i = 0
        for chem in l_chem_to_compute:
            i = i + 1 
            SMILES = chem[0]
            inchikey = chem[1]

            if i%1000 == 0:
                print(i, len(l_chem_to_compute))

            # check in bfore - only 23s for 1M chemicals
            # check is in db already -- save time
            # print(toolbox.binary_search(l_inch, inchikey))
            # if toolbox.binary_search(l_inch, inchikey) != -1:
            #     continue
            

            # the SMILES is already prepared
            self.c_CompDesc = CompDesc.CompDesc(input=SMILES, prdesc=self.p_dir_Desc)
            self.c_CompDesc.smi = SMILES
            self.c_CompDesc.inchikey = inchikey
            self.c_CompDesc.mol = Chem.MolFromSmiles(SMILES)
            try: self.c_CompDesc.computeAll2D() # case of the name has a *
            except: continue
            l_compute_desc = [self.c_CompDesc.all2D[desc_2D] if self.c_CompDesc.all2D[desc_2D] != "NA" else "-9999" for desc_2D in l_desc_2D]
            w1D2D = "{" + ",".join(["\"%s\"" % (compute_desc) for compute_desc in l_compute_desc]) + "}"

            self.c_CompDesc.set3DChemical(psdf3D = "")
            if self.c_CompDesc.err == 0:
                self.c_CompDesc.computeAll3D()
                l_compute_desc = [self.c_CompDesc.all3D[desc_3D] if self.c_CompDesc.all3D[desc_3D] != "NA" else "-9999" for desc_3D in l_desc_3D]
                w3D = "{" + ",".join(["\"%s\"" % (compute_desc) for compute_desc in l_compute_desc]) + "}"
                self.c_DB.addElement("chemical_description", ["inchikey", "desc_1d2d", "desc_3D", "map_name"], [inchikey, w1D2D, w3D, "dsstox"])

            else:
                self.c_DB.addElement("chemical_description", ["inchikey", "desc_1d2d", "map_name"], [inchikey, w1D2D, "dsstox"])
            
    def compute_coords(self, corVal, distributionVal, nb_load = 20000):
        """
        Function to compute coordinates in R --> do only the coordinates
        """
        p_desc_1D2D = self.p_dir_out + "1D2D.csv"
        if not path.exists(p_desc_1D2D):
            l_desc_1D2D = self.c_DB.execCMD("SELECT * FROM chem_descriptor_1d2d_name")
            l_desc_1D2D = [desc1D2D[1] for desc1D2D in l_desc_1D2D]

            nb_row = self.c_DB.execCMD("SELECT count(*)FROM chemical_description cd WHERE map_name ='dsstox'" )[0][0]
            
            # open file
            f_1D2D = open(p_desc_1D2D, "w")
            f_1D2D.write("inchikey\t%s\n"%("\t".join(l_desc_1D2D)))
            
            i = 1
            imax = nb_row
            while i < imax:
                l_chem = self.c_DB.execCMD("SELECT inchikey, desc_1d2d FROM chemical_description cd where map_name = 'dsstox' ORDER BY inchikey OFFSET %s ROWS FETCH NEXT %s ROWS ONLY"%(i, nb_load))
                i = i + nb_load
                
                for chem in l_chem:
                    f_1D2D.write("%s\t%s\n"%(chem[0], "\t".join(str(desc) for desc in chem[1])))
            f_1D2D.close()

        
            
        p_desc_3D = self.p_dir_out + "3D.csv"
        if not path.exists(p_desc_3D):
            l_desc_3D = self.c_DB.execCMD("SELECT * FROM chem_descriptor_3d_name")
            l_desc_3D = [desc3D[1] for desc3D in l_desc_3D]

            nb_row = self.c_DB.execCMD("SELECT count(*)FROM chemical_description cd WHERE map_name ='dsstox'" )[0][0]
            
            # open 3D
            f_3D = open(p_desc_3D, "w")
            f_3D.write("inchikey\t%s\n"%("\t".join(l_desc_3D)))

            i = 1
            imax = nb_row
            while i < imax:
                l_chem = self.c_DB.execCMD("SELECT inchikey, desc_3d FROM chemical_description cd where map_name = 'dsstox' AND desc_3d is not null ORDER BY inchikey OFFSET %s ROWS FETCH NEXT %s ROWS ONLY"%(i, nb_load))
                i = i + nb_load
                
                for chem in l_chem:
                    f_3D.write("%s\t%s\n"%(chem[0], "\t".join(str(desc) for desc in chem[1])))
            f_3D.close()            
        
        # create coords
        prmap = pathFolder.createFolder(self.p_dir_out + "map_" + str(corVal) + "-" + str(distributionVal) + "/")
        self.p_dir_map = prmap

        pcoordDim1Dim2 = prmap + "coord1D2D.csv"
        pcoordDim3D = prmap + "coord3D.csv"
        if not path.exists(pcoordDim1Dim2):
            runExternalSoft.RComputeMapFiles(p_desc_1D2D, "1D2D", prmap, corVal, distributionVal)

        self.p_desc_1D2D = p_desc_1D2D
        self.p_coords_1D2D = pcoordDim1Dim2


        if not path.exists(pcoordDim3D):
            runExternalSoft.RComputeMapFiles(p_desc_3D, "3D", prmap, corVal, distributionVal)
       
        self.p_desc_3D = p_desc_3D
        self.p_coords_3D = pcoordDim3D

    def draw_map(self):
        if not "p_desc_1D2D" in self.__dict__ or not "p_desc_3D" in self.__dict__:
            print("Compute coords first")
            return 
        
        runExternalSoft.RDrawHexaView(self.p_coords_1D2D, self.p_coords_3D, "DsstoxMap", "50", self.p_dir_map)

    def pushDB_coords(self):
        """
        will push in DB the coordinates
        - push 10 coords for 1D2D and 10 coords for 3D
        """
        if not path.exists(self.p_coords_1D2D) or not path.exists(self.p_coords_3D):
            print("ERROR: compute coords first")
            return 

        name_map = "dsstox"
        # load coords
        d_coords_1D2D = toolbox.loadMatrixToDict(self.p_coords_1D2D, sep = ",")
        d_coord_3D = toolbox.loadMatrixToDict(self.p_coords_3D, sep = ",")

        # only load l_chem with coord are not included
        cmd_extract = "SELECT inchikey from chemical_description cd  WHERE map_name='%s' AND d3_cube is NULL"%(name_map)
        l_inchikey = self.c_DB.execCMD(cmd_extract)
        self.c_DB.connOpen()
        
        i = 0
        imax = len(l_inchikey)
        while i < imax: 

            try:
                w_1D2D = "{" + ",".join(["\"%s\"" % (d_coords_1D2D[l_inchikey[i][0]]["DIM%s"%(k)]) for k in range(1, 11)]) + "}"
            except:
                w_1D2D = ""

            try:
                w_3D = "{" + ",".join(["\"%s\"" % (d_coord_3D[l_inchikey[i][0]]["DIM%s"%(k)]) for k in range(1, 11)]) + "}"
            except:
                w_3D = ""
            
            if w_3D != "" and w_1D2D != "":
                w_cube = "{\"%s\", \"%s\", \"%s\"}"%(d_coords_1D2D[l_inchikey[i][0]]["DIM1"], d_coords_1D2D[l_inchikey[i][0]]["DIM2"], d_coord_3D[l_inchikey[i][0]]["DIM1"])
                cmd_sql = "UPDATE chemical_description SET dim1d2d='%s', dim3d='%s', d3_cube='%s' WHERE inchikey='%s' AND map_name='%s'"%(w_1D2D, w_3D, w_cube, l_inchikey[i][0], name_map)
                self.c_DB.updateTable_run(cmd_sql)
            i = i + 1
        self.c_DB.connClose()


    def compute_onDB_neighbors(self):

        # select inchikey where neighbor is empty
        cmd_extract = "select inchikey from chemical_description where d3_cube is not null and map_name = 'dsstox' and neighbors_dim3 is null"
        l_inch = self.c_DB.execCMD(cmd_extract)

        for inch in l_inch:
            inchikey = inch[0]
            cmd_extract_neighbor = "SELECT inchikey FROM chemical_description "\
                "WHERE map_name = 'dsstox' ORDER BY cube(d3_cube) <->  (select cube (d3_cube) "\
                "FROM chemical_description where inchikey='%s' AND map_name = 'dsstox' limit (1))  limit (21)" %(inchikey)
            l_neighbor = self.c_DB.execCMD(cmd_extract_neighbor)
            l_neighbor = [neighbor[0] for neighbor in l_neighbor[1:]]
            

            # update database
            w_neighbors = "{" + ",".join(["\"%s\"" % (str(neighbor)) for neighbor in l_neighbor]) + "}" # remove the inchikey in the list
            cmd_update = "UPDATE chemical_description SET neighbors_dim3 = '%s' WHERE inchikey='%s' AND map_name='dsstox';" %(w_neighbors, inchikey)
            # print(cmd_update)
            self.c_DB.updateTable(cmd_update)  

        return 





#### propbably to remove

    def computeDesc (self, insertDB =0, w=0):

        if not "dchem" in self.__dict__:
            self.loadlistChem()

        pfilout1D2D = self.prout + "1D2D.csv"
        pfilout3D = self.prout + "3D.csv"

        if path.exists(pfilout1D2D) and path.exists(pfilout3D) and w == 1 and insertDB == 0:
            self.p1D2D = pfilout1D2D
            self.p3D = pfilout3D
            return


        lchemID = list(self.dchem.keys()) # can be shuffle
        if len(lchemID) < 50000 and self.iend == 0:
            shuffle(lchemID)

        imax = len(lchemID)
        if self.iend == 0 or self.iend > imax:
            iend = imax
        else:
            iend = self.iend

        lchemID = lchemID[self.istart:iend]
        shuffle(lchemID)
        ldesc1D2D = Chemical.getLdesc("1D2D")
        ldesc3D = Chemical.getLdesc("3D")

        if insertDB == 1:
            cDB = DBrequest.DBrequest()

        if w == 1:
            filout1D2D = open(pfilout1D2D, "w")
            filout1D2D.write("inchikey\t" + "\t".join(ldesc1D2D) + "\n")

            filout3D = open(pfilout3D, "w")
            filout3D.write("inchikey\t" + "\t".join(ldesc3D) + "\n")

        i = 0
        imax = len(lchemID)
        while i < imax:
            if i%1000 == 0:
                print (i)

            SMILESClean = self.dchem[lchemID[i]]["smiles_clean"]
            if SMILESClean == "NA":
                i = i + 1
                continue
            cChem = Chemical.Chemical(SMILESClean, self.prDesc)
            cChem.prepChem()

            # print(SMILESClean)
            if cChem.err == 0:
                # 2D descriptors
                cChem.computeAll2D(update=0)
                if cChem.err == 1:
                    i = i + 1
                    continue
                cChem.writeMatrix("2D")
                if w == 1:
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

                    if cChem.err == 1:
                        i = i + 1
                        continue
                    # write master table
                    if w == 1:
                        filout3D.write(
                            "%s\t%s\n" % (cChem.inchikey, "\t".join([str(cChem.all3D[desc]) for desc in ldesc3D])))

                    # put in table descriptors
                    if insertDB == 1:
                        out3D = cDB.getRow("desc_3d", "inchikey='%s'" % (cChem.inchikey))
                        if out3D == []:
                            valDesc = [cChem.all3D[desc3D] for desc3D in ldesc3D]
                            valDesc = ['-9999' if desc == "NA" else desc for desc in valDesc]

                            w3D = "{" + ",".join(["\"%s\"" % (desc) for desc in valDesc]) + "}"
                            cDB.addElement("desc_3d", ["inchikey", "desc_value"], [cChem.inchikey, w3D])
            i = i + 1

        if w == 1:
            filout1D2D.close()
            filout3D.close()

        self.p1D2D = pfilout1D2D
        self.p3D = pfilout3D

    def pushPropInDB(self):
        
        tableDB = "dsstox_chem"
        print("Push prop in DB: ")

        if not "dchem" in self.__dict__:
            self.loadlistChem()
        cDB = DBrequest.DBrequest()
        cDB.verbose = 0
        for chem in self.dchem.keys():
            cmdSQL = "SELECT count(*) FROM %s WHERE inchikey = '%s';"%(tableDB, self.dchem[chem]["inchikey"])
            findInch = cDB.execCMD(cmdSQL)[0][0]
            if findInch != 0:
                cmdSQL = "UPDATE %s SET %s = TRUE WHERE inchikey = '%s'"%(tableDB, self.nameMap, self.dchem[chem]["inchikey"])
                cDB.updateTable(cmdSQL)
        print("Finish push prop in DB <-")
        return 

    def computeCoords(self, corVal, distributionVal, insertDB=1):


        if not "p1D2D" in self.__dict__ and not "p3D" in self.__dict__:
            self.computeDesc(insertDB=0, w=1)

        # create coords
        prmap = pathFolder.createFolder(self.prout + "map_" + str(corVal) + "-" + str(distributionVal) + "/")
        self.prmap = prmap

        pcoordDim1Dim2 = prmap + "coord1D2D.csv"
        pcoordDim3D = prmap + "coord3D.csv"
        if not path.exists(pcoordDim1Dim2) or not path.exists(pcoordDim3D):
            runExternalSoft.RComputeMapFiles(self.p1D2D, self.p3D, prmap, corVal, distributionVal)

        if not path.exists(pcoordDim1Dim2) or not path.exists(pcoordDim3D):
            print("ERROR file map")
            return
        else:
            self.pcoords1D2D = pcoordDim1Dim2
            self.pcoords3D = pcoordDim3D

        if insertDB == 1:
            if self.nameMap == "dsstox":
                dcoord1D2D = toolbox.loadMatrixCoords(pcoordDim1Dim2, 10)
                dcoord3D = toolbox.loadMatrixCoords(pcoordDim3D, 10)
            else:
                dcoord1D2D = toolbox.loadMatrixToDict(pcoordDim1Dim2, ",")
                dcoord3D = toolbox.loadMatrixToDict(pcoordDim3D, ",")

            cDB = DBrequest.DBrequest()
            cDB.verbose = 0
            lchem = list(dcoord1D2D.keys())
            i = 0
            imax = len(lchem)
            while i < imax: 
                #out1D2D = cDB.getRow("%s_coords"%(self.nameMap), "inchikey='%s'" % (chem))
                #if out1D2D == []:
                    if self.nameMap == "dsstox":

                        w1D2D = "{" + ",".join(["\"%s\"" % (str(coord)) for coord in dcoord1D2D[lchem[i]]]) + "}"
                        w3D = "{" + ",".join(["\"%s\"" % (str(coord)) for coord in dcoord3D[lchem[i]]]) + "}"
                        cDB.addElement("%s_coords"%(self.nameMap), ["inchikey", "dim1d2d", "dim3d", "in_db"], [lchem[i], w1D2D, w3D, "1"])
                        
                        del dcoord1D2D[lchem[i]]
                        del dcoord3D[lchem[i]]
                        del lchem[i]
                        imax = imax - 1

                    else: 
                        nbdim1d2d = len(dcoord1D2D[lchem[i]].keys()) - 1
                        nbdim3d = len(dcoord3D[lchem[i]].keys()) - 1

                        w1D2D = "{" + ",".join(["\"%s\"" % (dcoord1D2D[lchem[i]]["DIM" + str(i)]) for i in range(1, nbdim1d2d + 1)]) + "}"
                        w3D = "{" + ",".join(["\"%s\"" % (dcoord3D[lchem[i]]["DIM3-" + str(i)]) for i in range(1, nbdim3d + 1)]) + "}"
                        cDB.addElement("%s_coords"%(self.nameMap), ["inchikey", "dim1d2d", "dim3d", "in_db"], [lchem[i], w1D2D, w3D, "1"])
                        
                        del dcoord1D2D[lchem[i]]
                        del dcoord3D[lchem[i]]
                        del lchem[i]
                        imax = imax - 1

    def runRprojection(self, corVal, distributionVal):

        if not "p1D2D" in self.__dict__ and not "p3D" in self.__dict__:
            self.computeDesc(insertDB=0)

        # create coords
        prproj = pathFolder.createFolder(self.prout + "proj_" + str(corVal) + "-" + str(distributionVal) + "/")
        runExternalSoft.RDrawProjection(self.p1D2D, self.p3D, prproj, corVal, distributionVal)

    def splitMap(self, nbsplit, dim, insertDB = 0):

        if not "prmap" in self.__dict__:
            print ("Generate the map files first")
            return

        else:
            prout = pathFolder.createFolder(self.prmap + "split_" + str(nbsplit) + "/")
            self.prmaps = prout

            if not "psplitMap" in self.__dict__:
                self.psplitMap = {}

            # generate only one file with chem and map
            if dim == 1:
                pfilout = prout + "mapx_split.csv"
                
            elif dim == 2: 
                pfilout = prout + "mapy_split.csv"

            else:
                pfilout = prout + "mapz_split.csv"

            self.psplitMap[dim] = pfilout
            if path.exists(pfilout) and insertDB == 0:
                return
            elif not path.exists(pfilout):
                coord1D2D = self.prmap + "coord1D2D.csv"
                coord3D = self.prmap + "coord3D.csv"

                if dim == 1 or dim == 2:
                    din = toolbox.loadMatrixCoords(coord1D2D, 2)
                else:
                    din = toolbox.loadMatrixCoords(coord3D, 2)

                # max and min 1D2D
                maxDim = 0.0
                minDim = 0.0

                nbchem = len(list(din.keys()))
                nbchembymap = int(nbchem/nbsplit)

                # calibrate max and min
                print("== Initiate calibration ==")
                for chem in din.keys():

                    if dim == 1 or dim == 3:
                        dimVal = din[chem][0]
                    elif dim == 2:
                        dimVal = din[chem][1]

                    if dimVal > maxDim:
                        maxDim = dimVal
                    if dimVal < minDim:
                        minDim = dimVal
                print("== End calibration ==")

                dmap = {}
                imap = 1
                dmap[imap] = []

                dimVal = minDim
                while dimVal < maxDim:
                    dimVal = dimVal + 0.10
                    if len(dmap[imap]) > nbchembymap:
                        imap = imap + 1
                        dmap[imap] = []
                    ichem = 0
                    lchem = list(din.keys())
                    nbchem = len(lchem)
                    while ichem < nbchem:

                        if dim == 1 or dim == 3:
                            valtemp = din[lchem[ichem]][0]
                        elif dim == 2:
                            valtemp = din[lchem[ichem]][1]


                        if valtemp < dimVal:

                            dmap[imap].append(deepcopy(lchem[ichem]))
                            del din[lchem[ichem]]
                            del lchem[ichem]
                            nbchem = nbchem - 1
                            continue
                        else:
                            ichem = ichem + 1

                print("==== Write output ====")
                filout = open(pfilout, "w")
                filout.write("inchikey\tmap\n")        
                for d in dmap.keys():
                    for chem in dmap[d]:
                        filout.write("%s\t%s\n"%(chem, d))
                filout.close()

        if insertDB == 1:
            cDB = DBrequest.DBrequest()
            #cDB.verbose = 1

            dmap = toolbox.loadMatrixToDict(pfilout)
            tableIn = "dsstox_coords"
            if dim == 1:
                mapIn = "mapx"
            elif dim == 2: 
                mapIn = "mapy"
            else:
                mapIn = "mapz"


            for chem in dmap.keys():
                inch = chem.replace("\"", "")
                
                cmdSQL = "UPDATE %s SET %s=%s WHERE inchikey='%s';" %(tableIn, mapIn, dmap[chem]["map"], inch)
                cDB.updateTable(cmdSQL)

    # have to be optimize
    def generateCentroidFile(self):

        if not "prmaps" in self.__dict__:
            print ("Generate Maps first")
            return

        if not "psplitMap" in self.__dict__:
            print("Generate the split map first")
            return

        else:
            lpfmap = list(self.psplitMap.values())
            #print(lpfmap)

            pfilout = self.prmaps + "centroids.csv"
            #if path.exists(pfilout):
            #    return 
            
            coords1D2D = toolbox.loadMatrixCoords(self.pcoords1D2D, 2)
            coords3D = toolbox.loadMatrixCoords(self.pcoords3D, 2)
            
            dout = {}
            for pmap in lpfmap:
                
                print(pmap)
                nameMap = pmap.split("/")[-1].split("_")[0]
                dmap = toolbox.loadMatrixToDict(pmap)
                print(nameMap)

                i = 1
                while 1:
                    lcoords = []
                    for chem in dmap.keys():
                        if int(dmap[chem]["map"]) == i :
                            lcoords.append([coords1D2D[chem][0], coords1D2D[chem][1], coords3D[chem][0]])
                    if lcoords == []:
                        break
                    else:
                        print(len(lcoords))
                        print(lcoords[0])
                        coordCentroid = calculate.centroid(lcoords)
                        dout[nameMap + "_" + str(i)] = coordCentroid
                    i = i + 1


            #print(dout)
            filout = open(pfilout, "w")
            filout.write("map\tx\ty\tz\n")
            for map in dout.keys():
                filout.write("%s\t%s\t%s\t%s\n"%(map, dout[map][0], dout[map][1], dout[map][2]))
            filout.close()

    def updateTableProp(self, nameMap):

        cDB = DBrequest.DBrequest()
        self.verbose = 0

        if not "dchem" in self.__dict__:
            self.loadlistChem()
        
        for chem in self.dchem:
            dsstox = chem.replace("\"", "")
            cmdSQL = "UPDATE dsstox_prop SET %s=true WHERE db_id='%s';" %(self.nameMap, dsstox)
            cDB.updateTable(cmdSQL)

    def pushDssToxNamePropInDB(self):

        cDB = DBrequest.DBrequest()
        cDB.verbose = 1
        i = 1
        for PROP in LPROP:
            cDB.addElement("dsstox_name_prop", ["id", "name"], [i, PROP])
            i = i + 1

    def generateTablePropAllDSSTOX(self, prDSSTOXPred, pknownSDF, pLD50, pDSSTOXMapOnCID, insertDB=0):

        
        pTableinfo = self.prout + "tablePropForDB.csv"
        if path.exists(pTableinfo) and insertDB == 0:
            self.pTableInAll = pTableinfo
            return

        #print ("LOAD INFO FROM DCHEM")
        #if not "dchem" in self.__dict__:
        #    self.loadlistChem()

        # intialisation ful dictionnary
        dDSSTOX = {}
        dmapCIDtoSID = {}
        print ("LOAD INFO MAP SID to CID")

        filMap = open(pDSSTOXMapOnCID, "r", encoding="utf8", errors="ignore")
        llines = filMap.readlines()
        filMap.close()
        
        lhead = llines[0].replace("\"", "")
        lhead = lhead.strip().split(",")
        #print(lhead)
        iDSSSID = lhead.index("dsstox_substance_id")
        iDSSCID = lhead.index("DSSTox_Structure_Id")
        iname = lhead.index("preferred_name")
        i = 1
        imax = len(llines)
        while i < imax:#####################################
            lineClean = toolbox.formatLine(llines[i])
            lelem = lineClean.strip().split(",")
            try:
                dDSSTOX[lelem[iDSSSID]] = {}
                dDSSTOX[lelem[iDSSSID]]["preferred_name"] = lelem[iname]
                dDSSTOX[lelem[iDSSSID]]["SMILES"] = self.dchem[lelem[iDSSSID]]["smiles_clean"]
                dDSSTOX[lelem[iDSSSID]]["inchikey"] = self.dchem[lelem[iDSSSID]]["inchikey"]
                dmapCIDtoSID[lelem[iDSSCID]] = lelem[iDSSSID]
            except:
                pass
            i = i + 1
        filMap.close()


        print("INIT DICTIONNARY")
        # put in dict out -> initialization to NA
        for chem in dDSSTOX.keys():
            for PROP in LPROP[3:]:
                try: dDSSTOX[chem][PROP] = "NA"
                except: break


        print("LOAD PRED")
        # load prediction and update table
        lppred = listdir(prDSSTOXPred)
        for ppred in lppred:##########################################
            if ppred[-3:] == "csv":
                print (ppred, "Load file")
                dtemp = toolbox.loadMatrixToDict(prDSSTOXPred + ppred, sep=",")
                k1 = list(dtemp.keys())[0]
                #print(dtemp[k1])
                #dddd
                for chemIDtemp in dtemp.keys():
                    DTXCID = dtemp[chemIDtemp]["MoleculeID"]
                    try: DTXSID = dmapCIDtoSID[DTXCID]
                    except: continue
                    for k in dtemp[chemIDtemp].keys():
                        if k in LPROP[3:]:
                            dDSSTOX[DTXSID][k] = dtemp[chemIDtemp][k]

        print ("PRED LOAD")

        print("LOAD SDF AND LD50")
        #load sdf
        dsdf = parseSDF.parseSDF(pknownSDF, "InChI Key_QSARr", self.prout)
        dsdf.parseAll()


        #load LD50 file
        dLD50 = toolbox.loadMatrixToDict(pLD50)
        print ("SDF and table LD50 loaded")

        
        for chem in dDSSTOX.keys():
            tempinchKey = dDSSTOX[chem]["inchikey"]
            # look sdf -> map on the sdf
            for dchemIDsdf in dsdf.lc:
                if dchemIDsdf["InChI Key_QSARr"] == tempinchKey:
                    for ksdf in dchemIDsdf.keys():
                        if ksdf in LPROP[3:]:
                            dDSSTOX[chem][ksdf] = dchemIDsdf[ksdf]

            # look in LD50 file -> map on the LD50
            for chemIDLD50 in dLD50.keys():
                if dLD50[chemIDLD50]["InChI Key_QSARr"] == tempinchKey:
                    for kLD50 in dLD50[chemIDLD50].keys():
                        if kLD50 in LPROP[3:]:
                            dDSSTOX[chem][kLD50] = dLD50[chemIDLD50][kLD50]



        print("WRITE TABLE")
        # load MAP

        filout = open(pTableinfo, "w")
        filout.write("ID\t%s\n"%("\t".join(LPROP)))
        for chem in dDSSTOX.keys():
            filout.write("%s\t%s\n"% (chem, "\t".join([str(dDSSTOX[chem][prop]) for prop in LPROP])))
        filout.close()

        self.pTableInAll = pTableinfo

    def pushTablePropAllInDB(self):

        if not "pTableInAll" in self.__dict__:
            print ("GENERATE TABLE FIRST")

        dtopush = toolbox.loadMatrixToDict(self.pTableInAll)
        cDB = DBrequest.DBrequest()
        cDB.verbose = 0
        i = 0
        i = 500000
        lchem = list(dtopush.keys())
        imax = len(lchem)
        while i < imax:
        #for chem in dtopush.keys():
            chem = lchem[i]
            outChem = cDB.getRow("dsstox_prop", "db_id='%s'" % (chem))
            if outChem == []:
                wprop = "{" + ",".join(["\"%s\"" % (dtopush[chem][PROP].replace("'", "")) for PROP in LPROP]) + "}"
                cDB.addElement("dsstox_prop", ["db_id", "prop_value"], [chem, wprop])
            i = i + 1
        return 

    def generateNeighborMatrix(self, nbNeighbor, lnDim):

        if not "pcoords1D2D" in self.__dict__:
            print("Compute Coord first")
            return 1
        else:
            if self.nameMap == "dsstox":
                # no N dimension because to slow
                prNeighbor = pathFolder.createFolder(self.prout + "Neighbors/")
                #pfilout = prNeighbor + "Table_DIM1D2D-2_1.csv"
                #if path.exists(pfilout):
                #    return
                

                    
                dDim1D2D = toolbox.loadMatrixCoords(self.pcoords1D2D, 2)
                dDim3D = toolbox.loadMatrixCoords(self.pcoords3D, 2)
                lpfmap = self.psplitMap
                lmap = []
                for imap in lpfmap.keys():
                    lmap.append(toolbox.loadMatrixToDict(lpfmap[imap]))
                #print(lmap)

                # from 1D2D coord
                lchem = list(dDim1D2D.keys())
                shuffle(lchem)
                i = 0
                imax = len(lchem)
                while i < imax:
                    inch = lchem[i]
                    
                    pfilout = prNeighbor + inch
                    if path.exists(pfilout):
                        i = i + 1
                        continue
                    filout = open(pfilout, "w")
                    filout.write("ID\tNeighbors\n")

                    # define map where we inspect
                    linmap = []
                    for dmap in lmap:
                        mapin = int(dmap[inch]["map"])
                        for chem in dmap.keys():
                            if int(dmap[chem]["map"]) == mapin or int(dmap[chem]["map"]) == (mapin + 1) or int(dmap[chem]["map"]) == (mapin - 1):
                                if not chem in lmap:
                                    linmap.append(chem)
                        
                    print(len(linmap))
                    ddist = {}
                    ddist[inch] = {}
                    for ID in linmap:
                        if ID != inch:
                            ddist[inch][ID] = sqrt(sum([(xi - yi) ** 2 for xi, yi in zip([dDim1D2D[ID][0], dDim1D2D[ID][1], dDim3D[ID][0]], [dDim1D2D[inch][0], dDim1D2D[inch][1], dDim3D[inch][0]],)]))

                    lID = [i[0] for i in sorted(ddist[inch].items(), key=lambda x: x[1])][:nbNeighbor]
                    filout.write("%s\t%s\n" % (inch, " ".join(lID)))
                    filout.close()
                    i = i + 1
            else:
                # compute all dimension openning withou restriction
                if lnDim == []:
                    dDim1D2D = toolbox.loadMatrixToDict(self.pcoords1D2D, sep = ",")
                    dDim3D = toolbox.loadMatrixToDict(self.pcoords3D, sep=",")
                    chem1 = list(dDim1D2D.keys())[0]
                    n1D2D = len(list(dDim1D2D[chem1].keys())) - 1
                    n3D = len(list(dDim3D[chem1].keys())) - 1
                    lnDim = [n1D2D, n3D]
                
                prNeighbor = pathFolder.createFolder(self.prout + "Neighbors/")
                pfilout = prNeighbor + "Table_DIM1D2D-" + str(lnDim[0]) + "_" + str(lnDim[1]) + ".csv"
                if path.exists(pfilout):
                    return
                else:
                    dDim1D2D = toolbox.loadMatrixToDict(self.pcoords1D2D, sep = ",")
                    dDim3D = toolbox.loadMatrixToDict(self.pcoords3D, sep=",")

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

    def pushDSSTOXNeighbors(self, prin):

        cDB = DBrequest.DBrequest()
        cDB.verbose = 0
        lfile = listdir(prin)
        for fileNeighbor in lfile:
            try:dneighbor = toolbox.loadMatrixToDict(prin + fileNeighbor)
            except: 
                remove(prin + fileNeighbor)
                continue
            inchkey = list(dneighbor.keys())[0]
            dneighbor[inchkey]["Neighbors"] = dneighbor[inchkey]["Neighbors"].split(" ")
            w3D = "{" + ",".join(["\"%s\"" % (neighbor) for neighbor in dneighbor[inchkey]["Neighbors"]]) + "}"
            cDB.addElement("dsstox_neighbors", ["inchikey", "neighbors_dim3"], [inchkey, w3D])
        return 

    def pushNeighbors(self):
        prneighbor = pathFolder.createFolder(self.prout + "Neighbors/")
        ptable3Dim = prneighbor + "Table_DIM1D2D-2_1.csv"
        ptableNDim = prneighbor + "Table_DIM1D2D-170_207.csv"
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
                #out1D2D = cDB.getRow("%s_neighbors"%(self.nameMap), "inchikey='%s'" % (chem))
                out1D2D = []
                if out1D2D == []:
                    w3D = "{" + ",".join(["\"%s\"" % (neighbor) for neighbor in ddist3D[chem]]) + "}"
                    wND = "{" + ",".join(["\"%s\"" % (neighbor) for neighbor in ddistND[chem]]) + "}"
                    cDB.addElement("%s_neighbors"%(self.nameMap), ["inchikey", "neighbors_dim3", "neighbors_dimn"], [chem, w3D, wND])










