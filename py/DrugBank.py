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
        
        self.l_prop = l_prop

    def pushDB_prop_value(self):

        
        # select chemical to update
        cmd_extract_chem = "select c.drugbank_id from chemicals c except select cpdv.drugbank_id  from chem_prop_drugbank_value cpdv"
        l_drugbank_id = self.c_DB.execCMD(cmd_extract_chem)
        l_drugbank_id = list(set([drugbank_id[0] for drugbank_id in l_drugbank_id]))
        l_drugbank_id.sort()

        if not "l_chem_sdf" in self.__dict__:
            self.parse_SDF()

        if not "l_prop" in self.__dict__:
            self.pushDB_prop_name()

        for d_chem_sdf in self.l_chem_sdf:
            drugbank_id = d_chem_sdf["DATABASE_ID"]
            if toolbox.binary_search(l_drugbank_id, drugbank_id) != -1:
                l_val = []
                for prop in self.l_prop:
                    try: val = d_chem_sdf[prop]
                    except: val = "NA"
                    if val == "":
                        val = "NA"
                    l_val.append(val)
                wprop = "{%s}"%(",".join(["\"%s\""%(toolbox.update_str4DB(val)) for val in l_val]))
                self.c_DB.addElement("chem_prop_drugbank_value", ["drugbank_id", "prop_value"], [drugbank_id, wprop])

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
            for d_chem in self.l_chem_sdf:
                if d_chem["smiles_clean"] == "":
                    self.c_DB.addElement(name_table, ["smiles_origin", "drugbank_id", "name"],
                                    [d_chem["SMILES"], d_chem["DRUGBANK_ID"], toolbox.update_str4DB(d_chem["GENERIC_NAME"])])

                else:
                    self.c_DB.addElement(name_table, ["smiles_origin", "smiles_clean", "inchikey", "drugbank_id", "name"],
                                    [d_chem["SMILES"], d_chem["smiles_clean"],
                                    d_chem["inchikey"], d_chem["DRUGBANK_ID"], toolbox.update_str4DB(d_chem["GENERIC_NAME"])])


        l_chem_out = []
        l_chem = self.c_DB.execCMD("select * from %s WHERE drugbank_id IS NOT NULL"%(name_table))
        for tup_chem in l_chem:
            d_chem = {"id":tup_chem[0], "smiles_clean": tup_chem[2]}
            l_chem_out.append(d_chem)
        
        self.l_chem_prep = l_chem_out

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

    def compute_pushDB_Desc(self):
        """
        !! value NA will be -9999 because of the database format
        """

        # first extract SMILES clean with ID from the chemicals database
        # extract all
        #cmd_SQL = "SELECT id, smiles_clean, inchikey FROM chemicals WHERE smiles_clean IS NOT NULL AND drugbank_id IS NOT NULL"
        
        #extract only missing value ++++
        cmd_SQL = "select c1.id, c1.smiles_clean, c1.inchikey from chemicals c1 where c1.smiles_clean is not null except select c.id, c.smiles_clean, c.inchikey from chemicals c join (select * from chemical_description where chemical_description.map_name = 'drugbank' ) cd on c.inchikey = cd.inchikey"
        l_chem = self.c_DB.execCMD(cmd_SQL)
        
        if len(l_chem) == 0:
            return

        l_desc_2D = self.c_CompDesc.getLdesc("1D2D")
        l_desc_3D = self.c_CompDesc.getLdesc("3D")

        self.c_DB.connOpen()
        for chem in l_chem:
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
        
    def compute_pushDB_Coords(self, corVal, distributionVal, n_dim=10):

        #create folder of coords
        p_dir_coords = pathFolder.createFolder(self.prout + "coords/")

        # load descriptors 1D and 2D
        # -> only take if desc computed
        cmd_SQL_desc = "select inchikey, desc_1d2d, desc_3d from chemical_description cd where desc_1d2d is not null and desc_3d is not null and d3_cube is null and map_name='drugbank'"
        l_chem_desc = self.c_DB.execCMD(cmd_SQL_desc)
        if len(l_chem_desc) == 0:
            return  

        # -> take desc name 1D2D
        cmd_sql_desc_name = "select id, name from chem_descriptor_1d2d_name"
        l_desc1D2D_name = self.c_DB.execCMD(cmd_sql_desc_name)
        l_desc1D2D_name = [chem_name[1] for chem_name in l_desc1D2D_name]

        cmd_sql_desc_name = "select id, name from chem_descriptor_3d_name"
        l_desc3D_name = self.c_DB.execCMD(cmd_sql_desc_name)
        l_desc3D_name = [chem_name[1] for chem_name in l_desc3D_name]

        # wrtie matrix
        # - 1D2D
        p_desc1D2D = p_dir_coords + "desc1D2D.csv"
        f_desc1D2D = open(p_desc1D2D, "w")
        f_desc1D2D.write("ID\t%s\n"%("\t".join(l_desc1D2D_name)))

        # -3D
        p_desc3D = p_dir_coords + "desc3D.csv"
        f_desc3D = open(p_desc3D, "w")
        f_desc3D.write("ID\t%s\n"%("\t".join(l_desc3D_name)))

        for chem_desc in l_chem_desc:
            chemid = chem_desc[0]
            l_desc_1D2D = [str(x) if x != -9999 else "NA" for x in chem_desc[1]]
            l_desc_3D = [str(x) if x != -9999 else "NA" for x in chem_desc[2]]
            f_desc1D2D.write("%s\t%s\n"%(chemid, "\t".join(l_desc_1D2D)))
            f_desc3D.write("%s\t%s\n"%(chemid, "\t".join(l_desc_3D)))

        f_desc1D2D.close()
        f_desc3D.close()

        # create coords
        p_dir_map = pathFolder.createFolder(p_dir_coords + "map_" + str(corVal) + "-" + str(distributionVal) + "/")
        p_coordDim1Dim2 = p_dir_map + "coord1D2D.csv"
        p_coordDim3D = p_dir_map + "coord3D.csv"
        
        if path.exists(p_coordDim1Dim2) and path.exists(p_coordDim3D):
            self.pcoords1D2D = p_coordDim1Dim2
            self.pcoords3D = p_coordDim3D
        elif not path.exists(p_coordDim1Dim2) or not path.exists(p_coordDim3D):
            runExternalSoft.RComputeMapFiles(p_desc1D2D, "1D2D", p_dir_map, corVal, distributionVal)
            runExternalSoft.RComputeMapFiles(p_desc1D2D, "3D", p_dir_map, corVal, distributionVal)
        elif not path.exists(p_coordDim1Dim2) or not path.exists(p_coordDim1Dim2):
            print("ERROR file map")
            err = 1

        self.pcoords1D2D = p_coordDim1Dim2
        self.pcoords3D = p_coordDim3D

        # push in database
        dcoord1D2D = toolbox.loadMatrixToDict(self.pcoords1D2D, sep = ",")
        dcoord3D = toolbox.loadMatrixToDict(self.pcoords3D, sep = ",")

        for chem in dcoord1D2D.keys():
            w1D2D = "{" + ",".join(["\"%s\"" % (dcoord1D2D[chem]["DIM" + str(i)]) for i in range(1, n_dim + 1)]) + "}"
            w3D = "{" + ",".join(["\"%s\"" % (dcoord3D[chem]["DIM3-" + str(i)]) for i in range(1, n_dim + 1)]) + "}"
            wcube = "{" + ",".join([dcoord1D2D[chem]["DIM1"], dcoord1D2D[chem]["DIM2"], dcoord3D[chem]["DIM3-1"]]) + "}"

            self.c_DB.updateTable("UPDATE chemical_description SET dim1d2d='%s', dim3d='%s', d3_cube='%s' WHERE inchikey='%s' AND map_name='drugbank';" %(w1D2D, w3D, wcube, chem))

    def draw_map(self, corVal, distributionVal):

        #create folder of coords
        p_dir_coords = pathFolder.createFolder(self.prout + "coords/")
        p_desc1D2D = p_dir_coords + "desc1D2D.csv"
        p_desc3D = p_dir_coords + "desc3D.csv"
        
        if path.exists(p_desc1D2D) and path.exists(p_desc3D):

            # create coords
            p_dir_proj = pathFolder.createFolder(p_dir_coords + "proj_" + str(corVal) + "-" + str(distributionVal) + "/")
            runExternalSoft.RComputeCor(p_desc1D2D, p_desc3D, p_dir_proj, corVal, distributionVal)
        else:
            print("ERROR - no descriptors files")

    def compute_onDB_neighbors(self):

        # select inchikey where neighbor is empty
        cmd_extract = "select inchikey from chemical_description where d3_cube is not null and map_name = 'drugbank' and neighbors_dim3 is null"
        l_inch = self.c_DB.execCMD(cmd_extract)

        for inch in l_inch:
            inchikey = inch[0]
            cmd_extract_neighbor = "SELECT inchikey FROM chemical_description "\
                "WHERE map_name = 'drugbank' ORDER BY cube(d3_cube) <->  (select cube (d3_cube) "\
                "FROM chemical_description where inchikey='%s' AND map_name = 'drugbank' limit (1))  limit (21)" %(inchikey)
            l_neighbor = self.c_DB.execCMD(cmd_extract_neighbor)
            l_neighbor = [neighbor[0] for neighbor in l_neighbor[1:]]
            

            # update database
            w_neighbors = "{" + ",".join(["\"%s\"" % (str(neighbor)) for neighbor in l_neighbor]) + "}" # remove the inchikey in the list
            cmd_update = "UPDATE chemical_description SET neighbors_dim3 = '%s' WHERE inchikey='%s';" %(w_neighbors, inchikey)
            self.c_DB.updateTable(cmd_update)  

        return 
    
    def pushDB_all(self):
        #self.parse_SDF()
        #self.pushDB_prop_name()
        #self.pushDB_prop_value()
        #self.pushDB_chem()
        #self.pushDB_desc_name() # that work be for every maps
        #self.compute_pushDB_Desc()
        #self.compute_pushDB_Coords(0.9, 95)
        self.draw_map(0.85, 90)
        self.draw_map(0.80, 95)
        self.draw_map(0.9, 95)
        self.draw_map(0.95, 90)
        #self.compute_onDB_neighbors()