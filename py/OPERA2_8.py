from re import search
from os import path, listdir

import DBrequest
import toolbox
from random import shuffle
from copy import deepcopy


L_STRP = ["MolWeight", "nbAtoms", "nbHeavyAtoms", "nbC", "nbO", "nbN", "nbAromAtom", "nbRing", "nbHeteroRing", "Sp3Sp2HybRatio", "nbRotBd", "nbHBdAcc", "ndHBdDon", "nbLipinskiFailures", "TopoPolSurfAir", "MolarRefract", "CombDipolPolariz"]


class OPERA2_8:
    def __init__(self, pr_pred, p_dir_out):
        self.pr_pred = pr_pred
        self.p_dir_out = p_dir_out
        self.version = "OPERA2.8"

        # load list of prop
        self.c_DB = DBrequest.DBrequest()
        cmd = "SELECT * FROM chem_descriptor_opera_name"
        l_prop = self.c_DB.execCMD(cmd)
        l_prop = [prop[1] for prop in l_prop]
        self.l_prop_OPERA = l_prop

    def load_prop(self):

        if self.prop == "":
            print("Need to load a valid prop")
            return 
        
        # add split for model with several prediction in it
        if self.prop in L_STRP:
            p_prop = self.pr_pred + "OPERA_Pred_StrP.csv"
        elif self.prop == "LogD55_pred" or self.prop == "LogD74_pred":
            p_prop = self.pr_pred + "OPERA_Pred_LogD.csv"
        elif search("CERAPP", self.prop):
            p_prop = self.pr_pred + "OPERA_Pred_CERAPP.csv"
        elif search("CoMPARA", self.prop):
            p_prop = self.pr_pred + "OPERA_Pred_CoMPARA.csv"
        elif search("CATMoS", self.prop):
            p_prop = self.pr_pred + "OPERA_Pred_CATMoS.csv"    
        elif search("pKa", self.prop) or search("ionization", self.prop):
            p_prop = self.pr_pred + "OPERA_Pred_pKa.csv"   
        elif search("BioDeg_LogHalfLife_pred", self.prop):
            p_prop = self.pr_pred + "OPERA_Pred_LogBioDeg.csv"  
        elif search("LogOH_pred", self.prop):
            p_prop = self.pr_pred + "OPERA_Pred_LogAOH.csv"
        else:
            p_prop = self.pr_pred + "OPERA_Pred_" + self.prop.split("_")[0] + ".csv"

        # check if file exist
        if path.exists(p_prop):
            d_prop = toolbox.loadMatrixToDict(p_prop, sep = ",")
            
        print("Loaded: ", self.prop, len(list(d_prop.keys())))    
        self.d_prop = d_prop

    def load_annotation(self):

        pr_anotation_opera28 = self.pr_pred + "DSSTox_082021_IDs_Structures/"
        #load 
        l_d_chemID = toolbox.loadMatrixToList(pr_anotation_opera28 + "DSSTox_082021_IDs.csv", sep = ",")
        d_structure = toolbox.loadMatrixToDict(pr_anotation_opera28 + "DSSTox_082021_Structures.csv", sep = ",")

        d_out = {}
        d_map_dtxsid_dtxcid = {}
        d_map_inch_dtxcid = {}

        # load ID
        for d_chemID in l_d_chemID:
            dtxcid = d_chemID["DSSTOX_COMPOUND_ID"]
            try: d_out[dtxcid]
            except:
                d_out[dtxcid] = {}
                d_out[dtxcid]["DSSTOX_COMPOUND_ID"] = []
                d_out[dtxcid]["DSSTOX_SUBSTANCE_ID"] = []
                d_out[dtxcid]["CASRN"] = []
                d_out[dtxcid]["PREFERRED_NAME"] = []
                d_out[dtxcid]["QC_LEVEL"] = []
                d_out[dtxcid]["Original_SMILES"] = d_structure[dtxcid]["Original_SMILES"]
                d_out[dtxcid]["Number of connected components"] = d_structure[dtxcid]["Number of connected components"]
                d_out[dtxcid]["Canonical_QSARr"] = d_structure[dtxcid]["Canonical_QSARr"]
                d_out[dtxcid]["Salt_Solvent"] = d_structure[dtxcid]["Salt_Solvent"]
                d_out[dtxcid]["InChI_Code_QSARr"] = d_structure[dtxcid]["InChI_Code_QSARr"]
                d_out[dtxcid]["InChI Key_QSARr"] = d_structure[dtxcid]["InChI Key_QSARr"]
                d_out[dtxcid]["Salt_Solvent_ID"] = d_structure[dtxcid]["Salt_Solvent_ID"]
            d_out[dtxcid]["DSSTOX_COMPOUND_ID"].append(dtxcid)
            d_out[dtxcid]["DSSTOX_SUBSTANCE_ID"].append(d_chemID["DSSTOX_SUBSTANCE_ID"])
            d_out[dtxcid]["CASRN"].append(d_chemID["CASRN"])
            d_out[dtxcid]["PREFERRED_NAME"].append(d_chemID["PREFERRED_NAME"])
            d_out[dtxcid]["QC_LEVEL"].append(d_chemID["QC_LEVEL"])
            
            d_map_dtxsid_dtxcid[d_chemID["DSSTOX_SUBSTANCE_ID"]] = dtxcid
            d_map_inch_dtxcid[d_out[dtxcid]["InChI Key_QSARr"]] = d_chemID["DSSTOX_SUBSTANCE_ID"]
        
        self.d_annotation = d_out
        self.d_map_dtxsid_dtxcid = d_map_dtxsid_dtxcid
        self.d_map_inch_dtxcid = d_map_inch_dtxcid
    
    def load_chem_OPERA_missing(self, table = "chemical_description"):

        if not "d_annotation" in self.__dict__:
            self.load_annotation()

        cmd_extract = "SELECT inchikey FROM %s WHERE desc_opera IS NULL"%(table)
        l_inch = self.c_DB.execCMD(cmd_extract)
        l_inch = [inch [0] for inch in l_inch]
        shuffle(l_inch)
        print(len(l_inch))

        # reduce l_inch 
        i = 0
        imax = len(l_inch)
        while i < imax:
            try:
                self.d_map_inch_dtxcid[l_inch[i]]
                i = i + 1
            except:
                del l_inch[i]
                imax = imax - 1
                
        for prop in self.l_prop_OPERA:
            self.prop = prop
            p_filout_prop = self.p_dir_out + "opera_%s_to_push.csv"%(self.prop)
            if path.exists(p_filout_prop):
                continue
            
            self.load_prop()
            
            filout = open(p_filout_prop, "w")
            filout.write("inchikey,%s\n"%(self.prop))
            
            for inchkey in l_inch:
                if prop == "LogOH_pred":
                    prop_table = "LogAOH_pred"
                elif prop == "BioDeg_LogHalfLife_pred":
                    prop_table = "LogBioDeg_pred"
                elif prop == "ionization":
                    prop_table = "ionization_pred"
                elif prop == "LogKoc_pred":
                    prop_table = "LogKOC_pred"
                else:
                    prop_table = prop
                try: val = self.d_prop[self.d_map_inch_dtxcid[inchkey]][prop_table]
                except: val = "-9999"

                if val == "" or val == "?":
                    val = "-9999"

                filout.write("%s,%s\n"%(inchkey, val))
            filout.close()
    
    def pushDB_OPERA_missing(self, table = "chemical_description"):

        cmd_extract = "SELECT inchikey FROM %s WHERE desc_opera IS NULL"%(table)
        l_inch_to_update = self.c_DB.execCMD(cmd_extract)
        l_inch_to_update = [inch [0] for inch in l_inch_to_update]
        
        l_files = listdir(self.p_dir_out)
        file_0 = l_files[0]
        d_prop = toolbox.loadMatrixToDict(self.p_dir_out + file_0, sep = ",")
        l_inch_computed = list(d_prop.keys())
        
        # take only the intersect to not reproduce too much
        l_inch = list(set(l_inch_to_update) & set(l_inch_computed))
        
        print("Chem to update: %s"%(len(l_inch)))
        d_topush = {}
        for inch in l_inch:
            d_topush[inch] = []

        for prop in self.l_prop_OPERA:
            p_prop = self.p_dir_out + "opera_%s_to_push.csv"%(prop)
            d_prop = toolbox.loadMatrixToDict(p_prop, sep=",")
            for inch in l_inch:
                d_topush[inch].append(str(d_prop[inch][prop]))
        
        
        for inch in l_inch:
            l_pred = d_topush[inch]
            l_pred = list(map(lambda x: x.replace('?', '-9999'), l_pred))
            l_pred = [pred if pred != "" else "-9999" for pred in l_pred ]
            d_topush[inch] = l_pred
            w_opera = "{" + ",".join(["\"%s\"" % (d_topush[inch][k]) for k in range(0, 49)]) + "}"
            cmd_sql = "UPDATE chemical_description SET desc_opera='%s' WHERE inchikey='%s'"%(w_opera, inch)
            self.c_DB.connOpen()
            self.c_DB.updateTable_run(cmd_sql)
            self.c_DB.connClose()

    def load_all_chem_info_opera_missing(self):

        if not "d_annotation" in self.__dict__:
            self.load_annotation()

        p_filout = self.p_dir_out + "chem_to_compute.csv"

        cmd_extract = "SELECT c.smiles_origin, c.smiles_clean, c.dsstox_id, c.casn FROM chemical_description cd inner join chemicals c ON c.inchikey = cd.inchikey WHERE desc_opera IS NULL"
        l_chem_to_extract = self.c_DB.execCMD(cmd_extract)
        filout = open(p_filout, "w")
        filout.write("smiles_origin,smiles_clean,dsstox_id,casn\n")
        for chem_to_extract in l_chem_to_extract:
            filout.write("%s,%s,%s,%s\n"%(chem_to_extract[0], chem_to_extract[1], chem_to_extract[2], chem_to_extract[3]))
        filout.close()

        filout_smi = open(self.p_dir_out + "chem_to_compute.smi", "w")
        for chem_to_extract in l_chem_to_extract:
            filout_smi.write("%s\t%s\n"%(chem_to_extract[1], chem_to_extract[2]))
        filout_smi.close()