import DBrequest
import CompDesc
import toolbox

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
