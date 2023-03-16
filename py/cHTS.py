import DBrequest
import pandas as pd
import toolbox

class cHTS:
    def __init__(self, p_dir_out):
        self.c_DB = DBrequest.DBrequest()
        self.p_dir_out = p_dir_out

    def pushDB_chemicalQC(self, p_chemQC, table):

        df_chemicals = pd.read_excel(p_chemQC, sheet_name="cHTS Chemical QC DATA")
        for i_chem, row in df_chemicals.iterrows():
            self.c_DB.addElement(table, ["ncgc_id", "tox21_id", "casn", "dsstox_id", "chem_type", "tox21_qc_t0", "tox21_qc_t4", "NICEATM_qc_t0", "NICEATM_qc_t4", "NICEATM_qc_summary_call"], 
                                [row["ncgc_id"], row["tox21_id"], row["casrn"], row["dtxsid"], row["chem_type"], row["tox21_qc_t0"], row["tox21_qc_t4"], row["NICEATM_qc_t0"], row["NICEATM_qc_t4"], row["NICEATM_qc_summary_call"]])

    def pushDB_assay(self, p_assay, table):
        
        df_assay = pd.read_excel(p_assay, sheet_name="invitrodb34_AssayAnnotation")
        df_mechnistic = pd.read_excel(p_assay, sheet_name="AllMTMOA")
        for i_chem, row in df_assay.iterrows():
            assay = row["Assay"]
            mechanistic_target = df_mechnistic.loc[df_mechnistic['new_AssayEndpointName'] == assay, "MechanisticTarget"].to_list()
            mechanistic_target = toolbox.update_str4DB(",".join(list(set(mechanistic_target))))

            self.c_DB.addElement(table, ["assay", "assay_source", "species", "tissue", "invitro_assay_format", "gene", "entrez_gene_id", "mechanistic_target"], 
                                [row["Assay"], row["Assay Source"], row["Species"], row["Tissue"], row["Invitro Assay Format"], row["Gene"], row["Entrez Gene ID"], mechanistic_target])

    
    def pushDB_invitroDB(self, p_cHTS_invitro_DB, table):
        '''
        To low not use prefer an import in block
        '''
        l_row_invitroDB = toolbox.loadMatrixToList(p_cHTS_invitro_DB, sep = "\t")
        i=0
        imax = 10000
        self.c_DB.connOpen()
        for d_row in l_row_invitroDB:
            if i%imax == 0:
                self.c_DB.connClose()
                self.c_DB.connOpen()
            self.c_DB.addElement(table, ["recordid", "casn", "dsstox_id", "assay", "curve_flag_description", "chemical_qc_description", "test_range", "test_range_unit", "endpoint", "response", "response_unit", "reference"], 
                                [d_row["RecordID"], d_row["CASRN"], d_row["DTXSID"], d_row["Assay"], d_row["Curve Flag Description"], d_row["Chemical QC Description"], d_row["TestRange"], d_row["TestRange.Unit"], d_row["Endpoint"], d_row["Response"], 
                                d_row["ResponseUnit"], d_row["Reference"]], openDB=0)


        self.c_DB.connClose()
