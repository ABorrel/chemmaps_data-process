from os import listdir
import toolbox
import DBrequest

from re import search

class tripod:

    def __init__(self, pr_tripodFiles):

        self.pr_tripodFiles = pr_tripodFiles
        self.DB = DBrequest.DBrequest(verbose=0)

    
    def pushInDB(self, table_name):
        
        if not "l_col" in self.get_colnames():
            self.get_colnames()

        self.DB.connOpen()
        lp_tripodFile = listdir(self.pr_tripodFiles)
        for p_tripodFile in lp_tripodFile:
            if not search("aggregrated", p_tripodFile):
                continue
            l_samples = toolbox.loadMatrixToList(self.pr_tripodFiles + "/" + p_tripodFile, sep = "\t")


            print(self.pr_tripodFiles + "/" + p_tripodFile)
            for d_sample in l_samples:
                l_val = []
                for col in self.l_col:
                    if col == "ac50" and d_sample[col.upper()] == "":
                        l_val.append(0.0)
                    elif col == "curve_rank" and d_sample[col.upper()] == "":
                        l_val.append(0.0)
                    elif d_sample[col.upper()] == "":
                        l_val.append("NA")
                    else:
                        l_val.append(d_sample[col.upper()].replace("'", "''"))
                self.DB.addMultipleElem(table_name, self.l_col, l_val)
        self.DB.connClose()


    def get_colnames(self):
        
        lp_tripodFile = listdir(self.pr_tripodFiles)
        p_tripodFile = lp_tripodFile[0] 
        print(p_tripodFile)
        d_tripod = toolbox.loadMatrixToDict(self.pr_tripodFiles + "/" + p_tripodFile)


        l_col = [col.lower() for col in list(d_tripod[list(d_tripod.keys())[0]].keys())]
        self.l_col = l_col

        return l_col
    
    def insertAssaysTox21(self):

        cmd = "INSERT INTO tox21_assays (protocol_name, assay_target,target_category, cell_line, cell_target) VALUES('tox21-ar-mda-kb2-luc-antagonist-p1', 'AR-MDA agonist', 'NR', 'MDA-MB-453', 'Breast Cancer');"
        cmd = "INSERT INTO tox21_assays (protocol_name, assay_target,target_category, cell_line, cell_target) VALUES('tox21-sbe-bla-agonist-p1', 'SBE-BLA (TGF-beta) agonist', 'Developmental Toxicity', 'HEK293', 'Kidney');"