import pathFolder
import toolbox
import runExternalSoft
import DBrequest
import calculate
import parseSDF

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
        shuffle(l_inch)

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

    def pushDB_all(self):
        self.pushDB_chem()
        self.compute_pushDB_Desc()
        self.compute_coords(0.9, 90) # coord are compute in R
        self.pushDB_coords()
        self.draw_map()
        self.compute_onDB_neighbors() 




