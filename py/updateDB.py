import DBrequest
import toolbox
import parseSDF
from random import randint
from os import listdir
from copy import deepcopy


def UpdateDBChemPropVal(prPred, pInde, psdf, LPROP, prout):

    # load SDF
    cSDF = parseSDF.parseSDF(psdf, "CASRN", prout)
    cSDF.parseAll()

    dSDF = {}
    for chem in cSDF.lc:
        if not chem["Original_SMILES"] in list(dSDF.keys()):
            dSDF[chem["Original_SMILES"]] = deepcopy(chem) 
     

    # pIdentifier
    dSMILES = toolbox.loadMatrixToDict(pInde, sep = ",")
    dSMILES_out = {}
    for chem in dSMILES.keys():
        try:dSMILES_out[dSMILES[chem]["DTXSID"]] = deepcopy(dSMILES[chem])
        except: continue

    # load chem prediction
    lfilePred = listdir(prPred)
    dpred = {}
    for filePred in lfilePred[:1]:
        print(filePred)
        dtemp = toolbox.loadMatrixToDict(prPred + filePred, sep=",")

        # primary key will be DTXSID
        for DTXCID in dtemp.keys():
            DTXSID = dtemp[DTXCID]["dsstox_substance_id"]
            if not DTXSID in list(dpred.keys()):
                dpred[DTXSID] = dtemp[DTXCID]


    
        
    # write for DB Update
    pfiloutDesc = prout + "OPERA_desc_" + "update.csv"
    filoutDesc = open(pfiloutDesc, "w")
    filoutDesc.write("DTXSID\t%s\n"%"\t".join(LPROP))

    for DTXSID in dpred.keys():

        try:SMILES = dSMILES_out[DTXSID]["SMILES"]
        except: SMILES = "ERROR"
        iprop = 0
        imax = len(LPROP)
        lval = []
        err = 0
        while iprop < imax:
            PROP = LPROP[iprop]
            #print(PROP)
            try:
                val = str(dpred[DTXSID][PROP])
            except:
                try:
                    val = str(dSMILES_out[DTXSID][PROP])
                except:
                    try:
                        val = str(dSDF[SMILES][PROP])
                    except:
                        val = "NA"
                        print(PROP)
                        #err = 1
                        #    break
            
            if val == "NaN":
                val = "NA"
            lval.append(val)
            iprop = iprop + 1
        if err == 0:
            filoutDesc.write("%s\t%s\n"%(DTXSID, "\t".join(lval)))
    filoutDesc.close()



def UpdateDBName():

    return


# define file where DB update are included

ptrainTestOPERA = "/home/borrela2/data/updateDB_01-2020/ToxTrainTest_3D.sdf" #file with endpoint
#ptrainTestOPERA = "C:\\Users\\borrela2\\development\\trash\\ToxTrainTest_3D.sdf" #file with endpoint

pIdentifier = "/home/borrela2/data/updateDB_01-2020/DSSToxMS-Ready_Rprep_corrected.csv" #identification with name
#pIdentifier = "/home/borrela2/data/updateDB_01-2020/DSSTox_Identifiers_and_CASRN_04-2019.csv" #identification with name
#pIdentifier = "C:\\Users\\borrela2\\development\\trash\\DSSTox_Identifiers_and_CASRN_04-2019_corrected.csv"


prOPERAPred = "/home/borrela2/data/updateDB_01-2020/DSSTox_opera_01-2020/" #folder included prediction for OPERA
#prOPERAPred = "C:\\Users\\borrela2\\development\\trash\\DSSTox_opera_01-2020\\"

prout = "/home/borrela2/ChemMaps/update_server/" #out directory
#prout = "C:\\Users\\borrela2\\development\\trash\\"

LPROPCHEMAPS = ["inchikey", "SMILES", "preferred_name", "GHS_category", "EPA_category", "consensus_LD50", "LD50_mgkg", "MolWeight", "LogOH_pred", "CATMoS_VT_pred", "CATMoS_NT_pred", "CATMoS_EPA_pred",
                 "CATMoS_GHS_pred", "CATMoS_LD50_pred", "CERAPP_Ago_pred", "CERAPP_Anta_pred", "CERAPP_Bind_pred", "Clint_pred", "CoMPARA_Ago_pred", "CoMPARA_Anta_pred",
                 "CoMPARA_Bind_pred", "FUB_pred", "LogHL_pred", "LogKM_pred", "LogKOA_pred", "LogKoc_pred", "LogBCF_pred", "LogD55_pred", "LogP_pred", "MP_pred", "pKa_a_pred",
                 "pKa_b_pred", "ReadyBiodeg_pred", "RT_pred", "LogVP_pred", "LogWS_pred", "BioDeg_LogHalfLife_pred", "BP_pred", "nbLipinskiFailures"]

LPROPOPERA = ["MolWeight", "nbAtoms", "nbHeavyAtoms", "nbC", "nbO", "nbN", "nbAromAtom", "nbRing", "nbHeteroRing", "Sp3Sp2HybRatio", "nbRotBd", "nbHBdAcc", "ndHBdDon","nbLipinskiFailures", "TopoPolSurfAir", "MolarRefract", "CombDipolPolariz", "LogP_pred", "AD_LogP", "AD_index_LogP", "Conf_index_LogP", "MP_pred", "AD_MP", "AD_index_MP","Conf_index_MP","BP_pred","AD_BP","AD_index_BP","Conf_index_BP","LogVP_pred","AD_VP","AD_index_VP","Conf_index_VP","LogWS_pred","AD_WS","AD_index_WS","Conf_index_WS","LogHL_pred","AD_HL","AD_index_HL","Conf_index_HL","RT_pred","AD_RT","AD_index_RT","Conf_index_RT","LogKOA_pred","AD_KOA","AD_index_KOA","Conf_index_KOA","ionization","pKa_a_pred","pKa_b_pred","AD_pKa","AD_index_pKa","Conf_index_pKa","LogD55_pred","LogD74_pred","AD_LogD","AD_index_LogD","Conf_index_LogD","LogOH_pred","AD_AOH","AD_index_AOH","Conf_index_AOH","LogBCF_pred","AD_BCF","AD_index_BCF","Conf_index_BCF","BioDeg_LogHalfLife_pred","AD_BioDeg","AD_index_BioDeg","Conf_index_BioDeg","ReadyBiodeg_pred","AD_ReadyBiodeg","AD_index_ReadyBiodeg","Conf_index_ReadyBiodeg","LogKM_pred","AD_KM","AD_index_KM","Conf_index_KM","LogKoc_pred","AD_LogKoc","AD_index_LogKoc","Conf_index_LogKoc"]


UpdateDBChemPropVal(prOPERAPred, pIdentifier, ptrainTestOPERA, LPROPCHEMAPS, prout)

