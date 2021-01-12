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
        self.cDB = DBrequest.DBrequest(verbose=0)

    def loadOPERAFileInClass(self):

        l_fopera = listdir(self.pr_OPERA_preproc)
        l_pr_OPERA_pred = []
        for foprea in l_fopera:
            if search("dsstox_names_casrns", foprea):
                self.p_chem_name = self.pr_OPERA_preproc + foprea
            elif search("InChIKey_QSAR", foprea):
                self.p_chem_SMILES = self.pr_OPERA_preproc + foprea
            elif search("dtxsid_dtxcid_map", foprea):
                self.p_dtxsid_dtxcid_map = self.pr_OPERA_preproc + foprea
            elif search("OPERA", foprea):
                pr_OPERA_pred = self.pr_OPERA_preproc + foprea
                if path.isdir(pr_OPERA_pred):
                    l_pr_OPERA_pred.append(pr_OPERA_pred + "/")
        self.l_prOPERA_pred = l_pr_OPERA_pred
    
    def formatChemForToolChem(self, split_nbChem):

        pr_out = pathFolder.createFolder(self.pr_out + "ForToolchem/")
        self.pr_forToolChem = pr_out
        l_filin = listdir(pr_out)
        if len(l_filin) > 0:# case files are already computed 
            return 

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
    
    def formatPrepChemForToolChem(self):

        pr_desc = self.pr_out + "DESC/"

        l_file_chem = listdir(self.pr_forToolChem)
        for file_chem in l_file_chem:
            if file_chem != "chemicals_listNew.csv":##############################################################
                continue##########################################################################################
            p_filout = self.pr_forToolChem + file_chem[:-4] + "_chemPrep.csv"
            filout = open(p_filout, "w")
            filout.write("dsstox_id\tsmiles_origin\tsmiles_cleaned\tinchikey\tdrugbank_id\tcasn\tname\n")
            l_chem = toolbox.loadMatrixToList(self.pr_forToolChem + file_chem)
            i = 0
            imax = len(l_chem)
            while i < imax:
                d_chem = l_chem[i]
                c_chem = CompDesc.CompDesc(d_chem["smiles_origin"], pr_desc)
                c_chem.prepChem()
                if c_chem.err == 0:
                    c_chem.computeAll2D()
                    if c_chem.err == 0:
                            filout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(d_chem["dsstox_id"], d_chem["smiles_origin"], c_chem.smi, c_chem.inchikey, "NA", d_chem["casn"], d_chem["name"]))
                i = i + 1
            filout.close()

    def formatDesc2DForToolChem(self):

        c_chem = CompDesc.CompDesc("", "")
        l_desc2D = c_chem.getLdesc("2D")

        pr_desc = self.pr_out + "DESC/"

        l_file_chem = listdir(self.pr_forToolChem)
        for file_chem in l_file_chem:
            if file_chem != "chemicals_listNew.csv":##############################################################
                continue##########################################################################################
            p_filout = self.pr_forToolChem + file_chem[:-4] + "_desc2D.csv"
            filout = open(p_filout, "w")
            filout.write("inchikey\t" + "\t".join(l_desc2D))
            l_chem = toolbox.loadMatrixToList(self.pr_forToolChem + file_chem)
            i = 0
            imax = len(l_chem)
            while i < imax:
                d_chem = l_chem[i]
                c_chem = CompDesc.CompDesc(d_chem["smiles_origin"], pr_desc)
                c_chem.prepChem()
                if c_chem.err == 0:
                    c_chem.computeAll2D()
                    if c_chem.err == 0:
                        filout.write("%s\t%s\n"%(c_chem.inchikey, "\t".join([str(c_chem.all2D[desc]) for desc in l_desc2D])))
                i = i + 1
            filout.close()

    def formatDesc3DForToolChem(self):
        c_chem = CompDesc.CompDesc("", "")
        l_desc3D = c_chem.getLdesc("3D")

        pr_desc = self.pr_out + "DESC/"

        l_file_chem = listdir(self.pr_forToolChem)
        for file_chem in l_file_chem:
            if file_chem != "chemicals_listNew.csv":##############################################################
                continue##########################################################################################
            p_filout = self.pr_forToolChem + file_chem[:-4] + "_desc3D.csv"
            filout = open(p_filout, "w")
            filout.write("inchikey\t" + "\t".join(l_desc3D))
            l_chem = toolbox.loadMatrixToList(self.pr_forToolChem + file_chem)
            i = 0
            imax = len(l_chem)
            while i < imax:
                d_chem = l_chem[i]
                c_chem = CompDesc.CompDesc(d_chem["smiles_origin"], pr_desc)
                c_chem.prepChem()
                if c_chem.err == 0:
                    c_chem.set3DChemical()
                    if c_chem.err == 0:
                        c_chem.computeAll3D()
                        if c_chem.err == 0:
                            filout.write("%s\t%s\n"%(c_chem.inchikey, "\t".join([str(c_chem.all3D[desc]) for desc in l_desc3D])))
                i = i + 1
            filout.close()


        return

    def formatOPERAForToolChem(self):

        # to change
        c_chem = CompDesc.CompDesc("", "")
        l_desc2D = c_chem.getLdesc("2D")

        pr_desc = self.pr_out + "DESC/"

        l_file_chem = listdir(self.pr_forToolChem)
        for file_chem in l_file_chem:
            if file_chem != "chemicals_listNew.csv":##############################################################
                continue##########################################################################################
            p_filout = self.pr_forToolChem + file_chem[:-4] + "_desc2D.csv"
            filout = open(p_filout, "w")
            filout.write("inchikey\t" + "\t".join(l_desc2D))
            l_chem = toolbox.loadMatrixToList(self.pr_forToolChem + file_chem)
            i = 0
            imax = len(l_chem)
            while i < imax:
                d_chem = l_chem[i]
                c_chem = CompDesc.CompDesc(d_chem["smiles_origin"], pr_desc)
                c_chem.prepChem()
                if c_chem.err == 0:
                    c_chem.computeAll2D()
                    if c_chem.err == 0:
                        filout.write("%s\t%s\n"%(c_chem.inchikey, "\t".join([str(c_chem.all2D[desc]) for desc in l_desc2D])))
                i = i + 1
            filout.close()

    def updateNameAndCAS(self, name_table):
        """Function use to update the chemical table => error in the prefered name"""

        cmd_SQL = "SELECT id, dsstox_id, casn, name FROM %s "%(name_table)
        l_chem_DB = self.cDB.execCMD(cmd_SQL)

        d_chem_DB = {}
        for chem_DB in l_chem_DB:
            d_chem_DB[chem_DB[1]] = [chem_DB[0], chem_DB[2], chem_DB[3]]


        d_dsstox_name = toolbox.loadMatrixToDict(self.p_chem_name)
        d_dsstox_SMILES = toolbox.loadMatrixToDict(self.p_chem_SMILES, sep = ",")

        i = 0
        l_dsstoxid = list(d_dsstox_name.keys())
        imax = len(l_dsstoxid)
        j=0
        while i < imax:
            if i % 50000 == 0:
                print(i)
            dsstox_id = l_dsstoxid[i]
            name = d_dsstox_name[dsstox_id]["preferred_name"].replace("'", "''")
            casrn = d_dsstox_name[dsstox_id]["casrn"]
            
            try:
                id_db = d_chem_DB[dsstox_id][0]
                name_db = d_chem_DB[dsstox_id][2]
                cas_db = d_chem_DB[dsstox_id][1]
            except:
                i = i + 1
                continue
            
            if name_db != name or casrn != cas_db:
                cmd_sql = "UPDATE %s SET casn = '%s', name = '%s' WHERE id='%s';"%(name_table, casrn, name, id_db)
                j = j + 1
                self.cDB.updateTable(cmd_sql)

            i = i + 1

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
            filout.write(chem + "\n")
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
                cChem.generateInchiKey()
                if cChem.err == 1:
                    print("Error inch: %s"%(l_chem_add[i]))
                    i = i + 1
                    continue

                # 2D desc
                cChem.computeAll2D()
                cChem.writeMatrix("2D")

                #3D desc
                cChem.set3DChemical()
                if cChem.err == 0:
                    cChem.computeAll3D()
                    if cChem.err == 1:
                        print("Error 3D desc: %s -- %s"%(l_chem_add[i], i))
                    else:
                        cChem.writeMatrix("3D")
                else:
                    print("Error 3D generation: %s -- %s"%(l_chem_add[i], i))
            else:
                print("Error prep: %s -- %s"%(l_chem_add[i], i))
            i = i + 1

    def pushOPERANameTable(self, name_table):
        """Extract only descriptors from OPERA pred"""

        l_desc_opera = []
        for pr_OPERA_pred in self.l_prOPERA_pred:
            l_p_fopera = listdir(pr_OPERA_pred)
            for p_fopera in l_p_fopera:
                fopera = open(pr_OPERA_pred + p_fopera, "r", encoding="utf8", errors='ignore' )
                line_first = fopera.readline()
                line_first = line_first.replace("\"", "")
                fopera.close()
                l_desc_opera = l_desc_opera + line_first.split(",")

        
        l_desc_opera = list(set(l_desc_opera))
        l_toadd = []
        for desc_opera in l_desc_opera:
            if desc_opera.endswith("_pred"):
                l_toadd.append(desc_opera)
        l_desc_prop =  ["MolWeight", "nbAtoms", "nbHeavyAtoms", "nbC", "nbO", "nbN" ,"nbAromAtom","nbRing","nbHeteroRing","Sp3Sp2HybRatio","nbRotBd","nbHBdAcc","ndHBdDon","nbLipinskiFailures","TopoPolSurfAir","MolarRefract","CombDipolPolariz","ionization"]
        l_toadd =  l_desc_prop + l_toadd 

        cmd_add = "INSERT INTO %s (id, name) VALUES"%(name_table)
        l_val_add = []
        i = 1
        for desc_toadd in l_toadd:
            l_val_add.append("(%s, '%s')"%(i, desc_toadd))
            i = i + 1

        cmd_add = cmd_add + ",".join(l_val_add)
        
        print(cmd_add)

    def updateMissingDTXSID(self, name_table):
        """Check if we can populate chemicals with no DTXSID with new update"""


        d_dsstox_SMILES = toolbox.loadMatrixToDict(self.p_chem_SMILES, sep = ",")
        d_dsstox_name = toolbox.loadMatrixToDict(self.p_chem_name)

        #extract chemical without DTXSID
        # see if chem included
        cmd_SQL = "SELECT id, smiles_origin FROM %s WHERE dsstox_id is null"%(name_table)
        l_chem_DB = self.cDB.execCMD(cmd_SQL)

        d_chem_DB = {}
        for chem_DB in l_chem_DB:
            print(chem_DB)
            d_chem_DB[chem_DB[1]] = chem_DB[0]
        
        for chem in d_dsstox_SMILES.keys():
            smiles = d_dsstox_SMILES[chem]["Original_SMILES"]
            try:
                id_chem = d_chem_DB[smiles]
                # update chemical
                cmd_sql = "UPDATE %s SET casn = '%s', name = '%s', dsstox_id='%s' WHERE id='%s';"%(name_table,  d_dsstox_SMILES[chem]["casrn"], d_dsstox_name[d_dsstox_SMILES[chem]["dsstox_substance_id"]]["preferred_name"].replace("'", "''"), d_dsstox_SMILES[chem]["dsstox_substance_id"], id_chem)
                print(cmd_sql)
                self.cDB.updateTable(cmd_sql)
            except:
                continue

    def updateSMILES(self, name_table = "chemicals"):
        """Function use to update the chemical table => check if smiles origin change"""

        d_dsstox_SMILES = toolbox.loadMatrixToDict(self.p_chem_SMILES, sep = ",")
        d_dsstox_name = toolbox.loadMatrixToDict(self.p_chem_name)
        self.pr_desc = pathFolder.createFolder(self.pr_out + "DESC/")

        #extract chemical without DTXSID
        # see if chem included
        cmd_SQL = "SELECT id, dsstox_id, smiles_origin, inchikey, smiles_clean FROM %s "%(name_table)
        l_chem_DB = self.cDB.execCMD(cmd_SQL)

        d_chem_DB = {}
        for chem_DB in l_chem_DB:
            d_chem_DB[chem_DB[1]] = [chem_DB[0], chem_DB[2], chem_DB[3], chem_DB[4]]
        
        i = 0
        for chem in d_dsstox_SMILES.keys():
            dsstox_id = d_dsstox_SMILES[chem]["dsstox_substance_id"]
            inchkey = d_dsstox_SMILES[chem]["InChI Key_QSARr"]
            smiles = d_dsstox_SMILES[chem]["Original_SMILES"]
            smiles_cleaned = d_dsstox_SMILES[chem]["Canonical_QSARr"]
            try:smiles_indb = d_chem_DB[dsstox_id][1] # case of chemical is not in the DB
            except:continue
            inchkey_db = d_chem_DB[dsstox_id][2]
            smiles_cleaned_db = d_chem_DB[dsstox_id][3]
            smiles_db = d_chem_DB[dsstox_id][1]
            if smiles != smiles_db:
                # recompute cleaned SMILES
                c_chem = CompDesc.CompDesc(smiles, self.pr_desc)
                c_chem.prepChem()
                if c_chem.err == 0:
                    c_chem.generateInchiKey()
                else:
                     c_chem.smi = None
                if c_chem.err == 0:
                    inchikey = c_chem.inchikey
                else:
                    inchikey = None

                if d_chem_DB[dsstox_id][2] != inchikey:
                    cmd_sql = "UPDATE %s SET smiles_origin = '%s', smiles_clean = '%s', inchikey='%s' WHERE id='%s';"%(name_table, smiles, c_chem.smi, inchikey, d_chem_DB[dsstox_id][0])

                else:
                    continue#cmd_sql = "UPDATE %s SET smiles_origin = '%s' WHERE id='%s';"%(name_table, smiles, d_chem_DB[dsstox_id][0])

                #print(smiles_cleaned,smiles_indb, dsstox_id)
                print(i)
                i = i + 1
                self.cDB.updateTable(cmd_sql)

        return 

    def pushNewChemInDB(self, name_table = "chemicals"):

        if not "l_chem_toadd" in self.__dict__:
            self.extractOnlyNewChem(name_table, field_comparison)

        d_dsstox_name = toolbox.loadMatrixToDict(self.p_chem_name)
        d_dsstox_SMILES = toolbox.loadMatrixToDict(self.p_chem_SMILES, sep = ",")
        
        self.pr_desc = pathFolder.createFolder(self.pr_out + "DESC/")
        id_chem = self.cDB.execCMD("SELECT MAX(id) FROM %s"%(name_table))[0][0]

        i = 0
        imax = len(self.l_chem_toadd)
        #imax = 100
        print(imax)
        cmd_all = ""

        while i < imax:
            if i % 1000 == 0:
                print(i)
            chem = self.l_chem_toadd[i]

            # case of 50K chemicals without SMILES
            try:
                smiles = d_dsstox_SMILES[chem]["Original_SMILES"]
            except:
                smiles = ""
            
            name = d_dsstox_name[chem]["preferred_name"]
            casrn = d_dsstox_name[chem]["casrn"]

            cChem = CompDesc.CompDesc(smiles, self.pr_desc)
            if smiles != "":
                cChem.prepChem()
            else:
                cChem.err = 1

            if cChem.err == 0:
                smi_cleaned = cChem.smi
                if cChem.err == 1:
                    inch = ""
                else:
                    inch = cChem.generateInchiKey()
            else:
                smi_cleaned = ""
                inch = ""
                
            id_chem = id_chem + 1
            if smiles == "":
                cmd_sql = "INSERT INTO %s (id, dsstox_id, casn, name) VALUES (%s, '%s', '%s', '%s');\n"%(name_table, id_chem, chem, casrn, name.replace("'", "''"))
            elif smi_cleaned == "":
                cmd_sql = "INSERT INTO %s (id, smiles_origin, dsstox_id, casn, name) VALUES (%s, '%s', '%s', '%s', '%s');\n"%(name_table, id_chem, smiles, chem, casrn, name.replace("'", "''"))
            elif inch == "":
                cmd_sql = "INSERT INTO %s (id, smiles_origin, smiles_clean, dsstox_id, casn, name) VALUES (%s, '%s', '%s', '%s', '%s', '%s');\n"%(name_table, id_chem, smiles, smi_cleaned, chem, casrn, name.replace("'", "''"))
            else:
                cmd_sql = "INSERT INTO %s (id, smiles_origin, smiles_clean, inchikey, dsstox_id, casn, name) VALUES (%s, '%s', '%s', '%s', '%s', '%s', '%s');\n"%(name_table, id_chem, smiles, smi_cleaned, inch, chem, casrn, name.replace("'", "''"))
            i = i + 1

            cmd_all = cmd_all + cmd_sql


        self.cDB.runCMDaddElement(cmd_all)
    
    def extractChemicalWithDescNoComputed(self, name_table = "chemical_description", map_name = "dsstox"):
        
        # selection list inchikey and dsstox from chemicals table 
        
        pr_out = pathFolder.createFolder(self.pr_out + "updateDSSTOX/")
        p_filout = pr_out + "chem_list_desc.txt"
        
        # interesting to speed up test but none in multi run
        #if path.exists(p_filout):
        #    filout = open(p_filout, "r")
        #    l_chem_toadd = filout.read().split("\n")
        #    filout.close()
        #    self.l_chem_toadd = l_chem_toadd
            
        #    return

        filout = open(p_filout, "w")
        # extract list of chemicals in the DB
        l_chem_chemicalsDB = self.cDB.execCMD("SELECT dsstox_id, inchikey FROM chemicals WHERE inchikey is not null AND dsstox_id is not NULL") # extract dsstoxID when inchikey is not null
        l_chem_chemicalsDescDB = self.cDB.execCMD("SELECT inchikey FROM %s WHERE inchikey is not null AND map_name='%s'"%(name_table, map_name)) # extract dsstoxID when inchikey is not null
        d_chem_chemicalsDB = {}
        for chem_DB in l_chem_chemicalsDB:
            d_chem_chemicalsDB[chem_DB[1]] = chem_DB[0]

        # l_inchikey computed
        l_inchikeyInChemDesc = []
        for inch in l_chem_chemicalsDescDB:
            l_inchikeyInChemDesc.append(inch[0])
        l_inchikeyInChemDesc.sort()

        l_out = []
        for inchInChemicals in d_chem_chemicalsDB.keys():
            if toolbox.binary_search(l_inchikeyInChemDesc, inchInChemicals) == -1:
                # check if is in DESC folder
                #print(self.pr_out + "DESC/2D/" + inchInChemicals + ".csv")
                if path.exists(self.pr_out + "DESC/2D/" + inchInChemicals + ".csv") and path.getsize(self.pr_out + "DESC/2D/" + inchInChemicals + ".csv") > 0:
                    continue
                else:
                    filout.write( d_chem_chemicalsDB[inchInChemicals] + "\n")
                    l_out.append(d_chem_chemicalsDB[inchInChemicals])

        filout.close()
        self.l_chem_toadd = l_out

    def updateDescOPERA(self):


        # 1. load opera descriptor used for the DB
        cmd_SQL = "SELECT id, name FROM chem_descriptor_opera_name_new "
        l_opera_desc = self.cDB.execCMD(cmd_SQL)

        d_opera_desc = {}
        for opera_desc in l_opera_desc:
            id_opera = int(opera_desc[0])
            desc_name = opera_desc[1]
            d_opera_desc[id_opera] = desc_name


        # 2. load precomputed OPERA desc
        d_dsstox2dsctox = toolbox.loadMatrixToDict(self.p_dtxsid_dtxcid_map, sep = "\t")
        #print(list(d_dsstox2dsctox.keys())[1])
        #print(d_dsstox2dsctox[list(d_dsstox2dsctox.keys())[1]])

        # 3. load inch and dtx id from chemicals table
        cmd_SQL = "SELECT dsstox_id, inchikey FROM chemicals WHERE inchikey is not null"
        l_chemicalsDB = self.cDB.execCMD(cmd_SQL)

        # 3.2 load inch without desc opera 
        cmd_SQL = "SELECT inchikey FROM chemical_description WHERE desc_opera is null"
        l_inchDB = self.cDB.execCMD(cmd_SQL)

        l_inch_toupdate = []
        for inch in l_inchDB:
            l_inch_toupdate.append(inch[0])
        l_inch_toupdate.sort()


        d_chemicalsDB = {}
        for chemicalsDB in l_chemicalsDB:
            inchikey = chemicalsDB[1]
            if toolbox.binary_search(l_inch_toupdate, inchikey) == -1:
                continue
            dsstoxid = chemicalsDB[0]
            if dsstoxid == None:
                continue

            try:
                d_chemicalsDB[inchikey].append(dsstoxid)
            except:
                d_chemicalsDB[inchikey] = [dsstoxid]


        d_update = {}
        for inchikey in d_chemicalsDB.keys():
            for dsstox_sid in d_chemicalsDB[inchikey]:
                try:dsstox_cid = d_dsstox2dsctox[dsstox_sid]["dsstox_compound_id"]
                except:continue
                d_update[dsstox_cid] = {}
                d_update[dsstox_cid]["inchikey"] = inchikey
                d_update[dsstox_cid]["opera"] = {}
                for i_desc_opera in d_opera_desc.keys():
                    d_update[dsstox_cid]["opera"][d_opera_desc[i_desc_opera]] = "-9999"
                break


        # 4. load OPERA desccriptors
        for pr_OPERA_pred in self.l_prOPERA_pred:
            l_p_fopera = listdir(pr_OPERA_pred)
            for p_fopera in l_p_fopera:
                print(pr_OPERA_pred + p_fopera)
                d_temp = toolbox.loadMatrixToDict(pr_OPERA_pred + p_fopera, sep = ",")

                for dtx_cid in d_update.keys():
                    for desc in d_update[dtx_cid]["opera"].keys():
                        try:
                            d_update[dtx_cid]["opera"][desc] = float(d_temp[dtx_cid][desc])
                        except:
                            pass


        # 5. update DB
        j = 0
        l_dtx_cid = list(d_update.keys())
        jmax = len(l_dtx_cid)
        shuffle(l_dtx_cid)
        while j < jmax:
            dtx_cid = l_dtx_cid[j]
            inchikey = d_update[dtx_cid]["inchikey"]

            l_toadd = []
            i = 1
            imax = len(list(d_opera_desc.keys())) + 1
            while i < imax:
                l_toadd.append(d_update[dtx_cid]["opera"][d_opera_desc[i]])
                i = i + 1

            wOPERA = "{" + ",".join(["\"%s\"" % (str(desc_val)) for desc_val in l_toadd]) + "}"
            cmd_sql = "UPDATE chemical_description SET desc_opera = '%s' WHERE inchikey='%s';"%(wOPERA,inchikey)
            self.cDB.verbose = 0
            self.cDB.updateTable(cmd_sql)

            if j % 1000 == 0:
                print(j)

            j = j + 1

    def pushChemDescInDB(self, map_name, table_descriptor):

        self.pr_desc = pathFolder.createFolder(self.pr_out + "DESC/")

        #1. extract list of descriptors
        #1.1 OPERA
        # NONE need to run updateDescOPERA after loading

        #1.2 2D descriptors
        d_namedesc_2D = self.pullDescName("chem_descriptor_1d2d_name")

        #1.3 3D descriptors
        d_namedesc_3D = self.pullDescName("chem_descriptor_3d_name")

        #2. Extract chemicals to add
        pr_out = pathFolder.createFolder(self.pr_out + "updateDSSTOX/")
        p_filout = pr_out + "chem_list_to_update.txt"
        filout = open(p_filout, "w")
        # extract list of chemicals in the DB
        l_chem_chemicalsDB = self.cDB.execCMD("SELECT dsstox_id, inchikey FROM chemicals WHERE inchikey is not null AND dsstox_id is not NULL") # extract dsstoxID when inchikey is not null
        l_chem_chemicalsDescDB = self.cDB.execCMD("SELECT inchikey FROM %s WHERE inchikey is not null AND map_name='%s'"%(table_descriptor, map_name)) # extract dsstoxID when inchikey is not null
        d_chem_chemicalsDB = {}
        for chem_DB in l_chem_chemicalsDB:
            d_chem_chemicalsDB[chem_DB[1]] = chem_DB[0]

        # l_inchikey computed
        l_inchikeyInChemDesc = []
        for inch in l_chem_chemicalsDescDB:
            l_inchikeyInChemDesc.append(inch[0])
        l_inchikeyInChemDesc.sort()

        l_inch_out = []
        for inchInChemicals in d_chem_chemicalsDB.keys():
            if toolbox.binary_search(l_inchikeyInChemDesc, inchInChemicals) == -1:
                filout.write( d_chem_chemicalsDB[inchInChemicals] + "\t" + inchInChemicals + "\n")
                l_inch_out.append(inchInChemicals)
        filout.close()

        shuffle(l_inch_out)# case of multirun
        self.l_chem_push_desc = l_inch_out


        # 3. push in database and check if descriptor are computed
        i = 0
        imax = len(l_inch_out)
        while i < imax:
            p_desc1d2 = "%s2D/%s.csv"%(self.pr_desc, l_inch_out[i])
            p_desc3d = "%s3D/%s.csv"%(self.pr_desc, l_inch_out[i])

            w2D = self.formDescToPushInDB(p_desc1d2, d_namedesc_2D)
            w3D = self.formDescToPushInDB(p_desc3d, d_namedesc_3D)  

            if w2D != "" and w3D != "":
                self.cDB.addElement(table_descriptor, ["inchikey", "map_name", "desc_1d2d", "desc_3d"], [l_inch_out[i], map_name, w2D, w3D])

            i = i + 1


        return 


    # toolbox to extract name of descriptor for description table
    def pullDescName(self, table):

        cmd_SQL = "SELECT id, name FROM %s"%(table)
        l_desc = self.cDB.execCMD(cmd_SQL)

        d_desc = {}
        for desc in l_desc:
            id_desc = int(desc[0])
            desc_name = desc[1]
            d_desc[id_desc] = desc_name
        
        return d_desc

    
    def formDescToPushInDB(self, p_desc, d_desc_name):

        if not path.exists(p_desc) or path.getsize(p_desc) == 0:
            return ""
        l_d_desc = toolbox.loadMatrixToList(p_desc)
        d_desc = l_d_desc[0]
        for desc in d_desc.keys():
            if d_desc == "NA":
                d_desc[desc] = -999

        l_id_desc = list(d_desc_name.keys())
        min_id = min(l_id_desc)
        max_id = max(l_id_desc)
        wSQL = "{" + ",".join(["\"%s\"" % ('-9999' if str(d_desc[d_desc_name[i + min_id]]) == "NA" else str(d_desc[d_desc_name[i + min_id]]) ) for i in range(0,len(list(d_desc_name.keys()))) ]) + "}"
        return wSQL

