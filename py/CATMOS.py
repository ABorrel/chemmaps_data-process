import DBrequest
import toolbox

import pandas as pd
import math

class CATMOS:

    def __init__(self, p_CATMOS):

        self.c_DB = DBrequest.DBrequest()
        self.p_CATMOS = p_CATMOS

    def pushDB_name(self, table = "chem_toxexp_name"):

        d_name = {"1":"very_toxic", "2":"nontoxic","3":"LD50_mgkg", "4": "EPA_category", "5":"GHS_category"}

        for id in d_name.keys():
            self.c_DB.addElement(table, ["id", "name"], [id, d_name[id]])


    def load_data(self):

        df_training = pd.read_excel(self.p_CATMOS, sheet_name="Training")
        df_test = pd.read_excel(self.p_CATMOS, sheet_name="Test")
        df_all = pd.concat([df_training, df_test])
        print(df_all.shape)

        self.df_all = df_all
    

    def pushDB_Tox(self, table = "chem_toxexp_value"):
        
        if not "df_all" in self.__dict__:
            self.load_data()

        for index, row in self.df_all.iterrows():
            dtxsid = row['DTXSID']
            try:
                float(dtxsid)
                continue
            except: pass
            l_val = [row["very_toxic"], row["nontoxic"], row["LD50_mgkg"], row["EPA_category"], row["GHS_category"]]
            l_w = []
            for val in l_val:
                if math.isnan(val) == True:
                    l_w.append("NA")
                elif val == True:
                    l_w.append("1")
                elif val == False:
                    l_w.append("0")
                else:
                    l_w.append(val)
            if l_w[2] != "NA":
                l_w[2] = str(l_w[2])
            if l_w[3] != "NA":
                l_w[3] = str(int(l_w[3]))
            if l_w[4] != "NA":
                l_w[4] = str(int(l_w[4]))
            
            w_tox = "{%s}"%(",".join(["\"%s\""%(toolbox.update_str4DB(val)) for val in l_w]))
            self.c_DB.addElement(table, ["dsstox_id", "prop_value"], [dtxsid, w_tox])
