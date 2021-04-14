from re import search
import parseSDF
import toolbox
import pathFolder
import DBrequest

class toxExp:

    def __init__(self, p_sdf, p_LD50, pr_out):

        self.p_sdf = p_sdf
        self.p_LD50 = p_LD50
        self.pr_out = pathFolder.createFolder(pr_out + "TOX_DB/")

        self.cDB = DBrequest.DBrequest()

    def update(self):
        
        self.loadsdf()
        self.loadLC50()
        #self.pushName("chem_toxexp_name")
        self.pushInDB("chem_toxexp_value")



    def loadsdf(self):

        self.dout = {}
        csdf = parseSDF.parseSDF(self.p_sdf, "Original_SMILES", self.pr_out)
        csdf.parseAll()

        for chem in csdf.lc:
            if chem["DTXSID"] != "":
                self.dout[chem["DTXSID"]] = {}
                try:self.dout[chem["DTXSID"]]["very_toxic"] = chem["very_toxic"]
                except:self.dout[chem["DTXSID"]]["very_toxic"] = "NA"
                try:self.dout[chem["DTXSID"]]["nontoxic"] = chem["nontoxic"]
                except:self.dout[chem["DTXSID"]]["nontoxic"] = "NA"
                try:self.dout[chem["DTXSID"]]["LD50_mgkg"] = chem["LD50_mgkg"]
                except:self.dout[chem["DTXSID"]]["LD50_mgkg"] = "NA"
                try:self.dout[chem["DTXSID"]]["EPA_category"] = chem["EPA_category"]
                except:self.dout[chem["DTXSID"]]["EPA_category"] = "NA"
                try:self.dout[chem["DTXSID"]]["GHS_category"] = chem["GHS_category"]
                except:self.dout[chem["DTXSID"]]["GHS_category"] = "NA"

    def loadLC50(self):

        d_LC50 = toolbox.loadMatrixToDict(self.p_LD50)
        for chem in d_LC50.keys():
            dtxsid = d_LC50[chem]["DTXSID"]
            if not dtxsid in list(self.dout.keys()):
                self.dout[dtxsid] = {}
                self.dout[dtxsid]["very_toxic"] = "NA"
                self.dout[dtxsid]["nontoxic"] = "NA"
                self.dout[dtxsid]["LD50_mgkg"] = "NA"
                self.dout[dtxsid]["EPA_category"] = "NA"
                self.dout[dtxsid]["GHS_category"] = "NA"
            
            self.dout[dtxsid]["LD50_mgkg_Literature"] = d_LC50[chem]["LD50_mgkg_Literature"]
            self.dout[dtxsid]["log(LD50_Literature)"] = d_LC50[chem]["log(LD50_Literature)"]
            self.dout[dtxsid]["consensus_LD50"] = d_LC50[chem]["consensus_LD50"]
            self.dout[dtxsid]["concordance_LD50"] = d_LC50[chem]["concordance_LD50"]

    def pushName(self, table_name):

        l_name = ["very_toxic", "nontoxic", "LD50_mgkg", "EPA_category", "GHS_category", "LD50_mgkg_Literature", "log(LD50_Literature)", "consensus_LD50", "concordance_LD50"]
        
        i = 1
        for name in l_name:
            self.cDB.addElement(table_name, ["id", "name"], [i, name])
            i = i + 1

    
    def pushInDB(self, table_value):

        l_name = ["very_toxic", "nontoxic", "LD50_mgkg", "EPA_category", "GHS_category", "LD50_mgkg_Literature", "log(LD50_Literature)", "consensus_LD50", "concordance_LD50"]
        
        for dtxsid in self.dout.keys():
            for name in l_name:
                if not name in list(self.dout[dtxsid].keys()):
                    self.dout[dtxsid][name] = "NA"
                elif self.dout[dtxsid][name] == "":
                    self.dout[dtxsid][name] = "NA"

            wadd = "{" + ",".join(["\"%s\"" % (self.dout[dtxsid][name]) for name in l_name]) + "}"
            self.cDB.addElement(table_value, ["dsstox_id", "prop_value"], [dtxsid, wadd])




class drugbank:

    def __init__(self, psdf):
        self.psdf = psdf


    def parseAll (self):

        lout = []

        filin = open(self.psdf, "r")
        handle_read = filin.read()
        filin.close()
        #print handle_read

        l_compound = handle_read.split("$$$$\n")

        for compound in l_compound:
            dcompound = {}
            dcompound["sdf"] = compound + "$$$$"
            llines = compound.split("\n")
            i = 0
            nblines = len(llines)
            while i < nblines:
                if search("> <", llines[i]):
                    kin = llines[i].split("> <")[1]
                    kin = kin.split(">")[0]
                    valuek = llines[i+1].strip()
                    dcompound[kin] = valuek
                i += 1
            lout.append(dcompound)
        self.lc=lout
        print (len(lout))





