import DBrequest
import CompDesc
import toolbox
import pathFolder
import runExternalSoft

from os import path
import random
from rdkit import Chem

class Tox21:

    def __init__(self, p_listChem, p_annotation, name_map, istart, iend, p_dir_Desc, p_dir_out):

        self.p_listChem = p_listChem
        self.p_annotation = p_annotation
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
    
    def prep_chem(self):
        
        l_chem = toolbox.loadMatrixToList(self.p_listChem, sep = "\t")
        p_filout = self.p_dir_out + "chemical_table.csv"

        if path.exists(p_filout):
            self.l_chem = toolbox.loadMatrixToList(p_filout)
            return

        for d_chem in l_chem:
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
        
        # Add the dsstox in the mapping
        d_annotation = toolbox.loadMatrixToDict(self.p_annotation, sep = ",")
        i = 0
        imax = len(l_chem)
        while i < imax:
            casrn = l_chem[i]["CAS"]
            DTXSID = d_annotation[casrn]["DTXSID"]
            if DTXSID == "N/A":
                del l_chem[i] 
                imax = imax - 1
                continue    
            l_chem[i]["DTXSID"] = DTXSID
            i = i + 1


        # write table for control -> after open and put in the DB
        filout = open(p_filout, "w", encoding="utf8")
        filout.write("dsstox_id\tsmiles_origin\tsmiles_clean\tinchikey\tname\n")
        for d_chem in l_chem:
            filout.write("%s\t%s\t%s\t%s\t%s\n" % (d_chem["DTXSID"], d_chem["SMILES"], d_chem["smiles_clean"], d_chem["inchikey"], d_chem["SAMPLE_NAME"]))
        filout.close() 
        self.l_chem = toolbox.loadMatrixToList(p_filout)
    
    def pushDB_chem(self, name_table="chemicals"):

        # check is chem is in it
        if not "l_chem" in self.__dict__:
            self.prep_chem()

        for d_chem in self.l_chem:
            cmd = "SELECT COUNT(*) FROM %s WHERE dsstox_id='%s'"%(name_table, d_chem["dsstox_id"])
            count = self.c_DB.execCMD(cmd)  
            if count[0][0] == 0:
                print(d_chem)
                if d_chem["smiles_clean"] == "":
                    self.c_DB.addElement(name_table, ["smiles_origin", "dsstox_id", "name"],
                                    [d_chem["smiles_origin"], d_chem["dsstox_id"], toolbox.update_str4DB(d_chem["name"])])

                else:
                    self.c_DB.addElement(name_table, ["smiles_origin", "smiles_clean", "inchikey", "dsstox_id", "name"],
                                    [d_chem["smiles_origin"], d_chem["smiles_clean"],
                                    d_chem["inchikey"], d_chem["dsstox_id"], toolbox.update_str4DB(d_chem["name"])])

    def compute_pushDB_Desc(self):
        """
        !! value NA will be -9999 because of the database format
        """
        # check if already in 
        if not "l_chem" in self.__dict__:
            self.prep_chem()        
        
        # need to check if inchkey is already in database
        cmd_inch = "select inchikey from chemical_description where map_name = 'tox21'"
        l_inch = self.c_DB.execCMD(cmd_inch)

        l_inch = [inch[0] for inch in l_inch]
        l_inch = list(set(l_inch))
        l_inch.sort()

        if len(self.l_chem) == 0:
            return

        l_desc_2D = self.c_CompDesc.getLdesc("1D2D")
        l_desc_3D = self.c_CompDesc.getLdesc("3D")

        l_chem_to_compute = [d_chem for d_chem in self.l_chem if toolbox.binary_search(l_inch, d_chem["inchikey"]) == -1]
        print("To compute:", len(l_chem_to_compute))

        # shuffle list
        random.shuffle(l_chem_to_compute)
        
        i = 0
        l_inch_added = []
        for d_chem_to_compute in l_chem_to_compute:
            i = i + 1 
            SMILES = d_chem_to_compute["smiles_clean"]
            inchikey = d_chem_to_compute["inchikey"]
            if inchikey in l_inch_added:
                continue
            else:
                l_inch_added.append(inchikey)

            if i%1000 == 0:
                print(i, len(l_chem_to_compute))


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
                self.c_DB.addElement("chemical_description", ["inchikey", "desc_1d2d", "desc_3D", "map_name"], [inchikey, w1D2D, w3D, "tox21"])

            else:
                self.c_DB.addElement("chemical_description", ["inchikey", "desc_1d2d", "map_name"], [inchikey, w1D2D, "tox21"])

    def compute_coords(self, corVal, distributionVal, nb_load = 20000):
        """
        Function to compute coordinates in R --> do only the coordinates
        """
        p_desc_1D2D = self.p_dir_out + "1D2D.csv"
        if not path.exists(p_desc_1D2D):
            l_desc_1D2D = self.c_DB.execCMD("SELECT * FROM chem_descriptor_1d2d_name")
            l_desc_1D2D = [desc1D2D[1] for desc1D2D in l_desc_1D2D]

            nb_row = self.c_DB.execCMD("SELECT count(*)FROM chemical_description cd WHERE map_name ='tox21'" )[0][0]
            
            # open file
            f_1D2D = open(p_desc_1D2D, "w")
            f_1D2D.write("inchikey\t%s\n"%("\t".join(l_desc_1D2D)))
            
            i = 1
            imax = nb_row
            while i < imax:
                l_chem = self.c_DB.execCMD("SELECT inchikey, desc_1d2d FROM chemical_description cd where map_name = 'tox21' ORDER BY inchikey OFFSET %s ROWS FETCH NEXT %s ROWS ONLY"%(i, nb_load))
                i = i + nb_load
                
                for chem in l_chem:
                    f_1D2D.write("%s\t%s\n"%(chem[0], "\t".join(str(desc) for desc in chem[1])))
            f_1D2D.close()

        
            
        p_desc_3D = self.p_dir_out + "3D.csv"
        if not path.exists(p_desc_3D):
            l_desc_3D = self.c_DB.execCMD("SELECT * FROM chem_descriptor_3d_name")
            l_desc_3D = [desc3D[1] for desc3D in l_desc_3D]

            nb_row = self.c_DB.execCMD("SELECT count(*)FROM chemical_description cd WHERE map_name ='tox21'" )[0][0]
            
            # open 3D
            f_3D = open(p_desc_3D, "w")
            f_3D.write("inchikey\t%s\n"%("\t".join(l_desc_3D)))

            i = 1
            imax = nb_row
            while i < imax:
                l_chem = self.c_DB.execCMD("SELECT inchikey, desc_3d FROM chemical_description cd where map_name = 'tox21' AND desc_3d is not null ORDER BY inchikey OFFSET %s ROWS FETCH NEXT %s ROWS ONLY"%(i, nb_load))
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
        
        return 
    
    def draw_map(self):
        if not "p_desc_1D2D" in self.__dict__ or not "p_desc_3D" in self.__dict__:
            print("Compute coords first")
            return 
        
        runExternalSoft.RDrawHexaView(self.p_coords_1D2D, self.p_coords_3D, "Tox21Map", "50", self.p_dir_map)
    
    def pushDB_coords(self):
        """
        will push in DB the coordinates
        - push 10 coords for 1D2D and 10 coords for 3D
        """
        if not path.exists(self.p_coords_1D2D) or not path.exists(self.p_coords_3D):
            print("ERROR: compute coords first")
            return 

        name_map = "tox21"
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

        return 
    
    def compute_onDB_neighbors(self):
        # select inchikey where neighbor is empty
        cmd_extract = "select inchikey from chemical_description where d3_cube is not null and map_name = 'tox21' and neighbors_dim3 is null"
        l_inch = self.c_DB.execCMD(cmd_extract)
        random.shuffle(l_inch)

        for inch in l_inch:
            inchikey = inch[0]
            cmd_extract_neighbor = "SELECT inchikey FROM chemical_description "\
                "WHERE map_name = 'tox21' ORDER BY cube(d3_cube) <->  (select cube (d3_cube) "\
                "FROM chemical_description where inchikey='%s' AND map_name = 'tox21' limit (1))  limit (21)" %(inchikey)
            l_neighbor = self.c_DB.execCMD(cmd_extract_neighbor)
            l_neighbor = [neighbor[0] for neighbor in l_neighbor[1:]]
            

            # update database
            w_neighbors = "{" + ",".join(["\"%s\"" % (str(neighbor)) for neighbor in l_neighbor]) + "}" # remove the inchikey in the list
            cmd_update = "UPDATE chemical_description SET neighbors_dim3 = '%s' WHERE inchikey='%s' AND map_name='tox21';" %(w_neighbors, inchikey)
            # print(cmd_update)
            self.c_DB.updateTable(cmd_update) 
