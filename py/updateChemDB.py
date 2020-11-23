from re import search
from os import listdir, path
from copy import deepcopy
from random import shuffle
import toolbox
import pathFolder
import DBrequest
import CompDesc

class updateChemDB:

    def __init__(self, name_update, pr_OPERA_preproc, pr_out):
        self.name_update = name_update
        self.pr_OPERA_preproc = pr_OPERA_preproc
        self.pr_out = pr_out

        self.loadOPERAFileInClass()    
        self.cDB = DBrequest.DBrequest()


    def loadOPERAFileInClass(self):

        l_fopera = listdir(self.pr_OPERA_preproc)
        for foprea in l_fopera:
            if search("dsstox_names_casrns", foprea):
                self.p_chem_name = self.pr_OPERA_preproc + foprea
            elif search("InChIKey_QSAR", foprea):
                self.p_chem_SMILES = self.pr_OPERA_preproc + foprea
    
    def formatChemForToolChem(self, split_nbChem):

        pr_out = pathFolder.createFolder(self.pr_out + "ForToolchem/")
        d_dsstox_name = toolbox.loadMatrixToDict(self.p_chem_name)
        d_dsstox_SMILES = toolbox.loadMatrixToDict(self.p_chem_SMILES, sep = ",")

        i = 0
        ifile = 0
        l_dsstoxid = list(d_dsstox_name.keys())
        l_dsstoxid = shuffle(l_dsstoxid)
        imax = len(l_dsstoxid)
        i_file = 0
        while i < imax:
            if i_file == split_nbChem or i == 0:
                i_file = 0
                if i == 0:
                    f_out = open("%schemlist_1.csv"%(pr_out), "w")
                else:
                    print(i/split_nbChem + 1)
                    f_out.close()
                    f_out = open("%schemlist_%i.csv"%(pr_out, i/split_nbChem + 1), "w")
                f_out.write("smiles_origin\tdsstox_id\tdrugbank_id\tname\tcasn\n")
            
            dsstox_id = l_dsstoxid[i]
            name = d_dsstox_name[dsstox_id]["preferred_name"]
            casrn = d_dsstox_name[dsstox_id]["casrn"]
            try:
                smiles_origin = d_dsstox_SMILES[dsstox_id]["Original_SMILES"]
            except:
                smiles_origin = ""
            f_out.write("%s\t%s\tNA\t%s\t%s\n"%(smiles_origin, dsstox_id, name, casrn))
            
            i = i + 1
            i_file = i_file + 1
        f_out.close()
    
    def updateChemicalsTableFromOPERAFile(self, name_table):
        """Function use to update the chemical table => error in the prefered name"""

        self.cDB.connOpen()

        d_dsstox_name = toolbox.loadMatrixToDict(self.p_chem_name)
        d_dsstox_SMILES = toolbox.loadMatrixToDict(self.p_chem_SMILES, sep = ",")

        i = 0
        l_dsstoxid = list(d_dsstox_name.keys())
        imax = len(l_dsstoxid)
        l_cmd_sql = []
        while i < imax:
            if i % 100 == 0:
                print(i)
            dsstox_id = l_dsstoxid[i]
            name = d_dsstox_name[dsstox_id]["preferred_name"]
            casrn = d_dsstox_name[dsstox_id]["casrn"]
            try:
                smiles_origin = d_dsstox_SMILES[dsstox_id]["Original_SMILES"]
            except:
                smiles_origin = ""

            #cmdCount = "SELECT COUNT(*) FROM chemicals WHERE dsstox_id='%s'"%(dsstox_id)
            #inDB = self.cDB.execCMDrun(cmdCount)[0][0]
            
            #if inDB == 1:
            cmd_sql = "UPDATE %s SET casn = '%s', name = '%s' WHERE dsstox_id='%s';"%(name_table, casrn, name.replace("'", "''"), dsstox_id)
            self.cDB.execCMDrun(cmd_sql)

            i = i + 1

        #pr_out = pathFolder.createFolder(self.pr_out + "SQL/")
        #filout = open(pr_out + "updateChem.sql", "w")
        #filout.write("UPDATE chemicals " + "THEN ".join(l_cmd_sql))
        #filout.close()

        self.cDB.connClose()


    def extractOnlyNewChem(self, name_table, field_comparison):

        pr_out = pathFolder.createFolder(self.pr_out + "updateDSSTOX/")
        p_filout = pr_out + "chem_list.txt"
        
        if path.exists(p_filout):
            filout = open(p_filout, "r")
            self.l_chem_toadd = filout.read().split("\n")
            filout.close()
            return

        filout = open(p_filout, "w")
        d_dsstox_name = toolbox.loadMatrixToDict(self.p_chem_name)
        d_dsstox_SMILES = toolbox.loadMatrixToDict(self.p_chem_SMILES, sep = ",")
        # extract list of chemicals in the DB
        l_chem_DB = self.cDB.execCMD("SELECT %s FROM %s"%(field_comparison, name_table))
        for chem_DB in l_chem_DB:
            chem = chem_DB[0]
            if chem == None:
                continue
            try:
                del d_dsstox_name[chem]
            except:
                pass

        for chem in d_dsstox_name.keys():
            try: 
                smi = d_dsstox_SMILES[chem]["Original_SMILES"]
                filout.write(chem + "\n")
            except:
                pass

            
        filout.close()

        self.l_chem_toadd = list(d_dsstox_name.keys())



    def computeDescNewChem(self):
        
        if not "l_chem_toadd" in self.__dict__:
            self.extractOnlyNewChem("chemicals", "dsstox_id")

        self.pr_desc = pathFolder.createFolder(self.pr_out + "DESC/")
        d_dsstox_SMILES = toolbox.loadMatrixToDict(self.p_chem_SMILES, sep = ",")


        l_chem_add = self.l_chem_toadd
        shuffle(l_chem_add)

        i = 0
        imax = len(self.l_chem_toadd)
        print(imax)
        while i < imax:
            if i % 1000 == 0:
                print(i)
            chem = l_chem_add[i]

            try:smiles = d_dsstox_SMILES[chem]["Original_SMILES"]
            except:
                print(i, ": ERROR in SMILES - ", chem)
                i = i + 1
                continue

            cChem = CompDesc.CompDesc(smiles, self.pr_desc)
            cChem.prepChem()

            if cChem.err == 0:
                cChem.smi
                cChem.generateInchiKey()
                if cChem.err == 1:
                    print("Error inch: %s"%(l_chem_add[i]))
                    i = i + 1
                    continue

                # 2D desc
                cChem.computeAll2D()

                #3D desc
                cChem.set3DChemical()
                if cChem.err == 0:
                    cChem.computeAll3D()
                else:
                    print("Error 3D generation: %s -- %s"%(l_chem_add[i], i))
            else:
                print("Error prep: %s -- %s"%(l_chem_add[i], i))
            i = i + 1



