import DBrequest
import toolbox
import parseSDF
from random import randint
from os import listdir


def UpdateDBVal(prPred, pIdentifier, pSDF, prout):

    #dIdentifier = toolbox.loadMatrixToDict(pIdentifier)
    #print(len(list(dIdentifier.keys())))
    csdf = parseSDF.parseSDF(pSDF, "InChI Key_QSARr", prout)
    csdf.parseAll()

    #print(len(csdf.lc))
    #print(csdf.lc[10])
    #dddd


    # load chem prediction
    lfilePred = listdir(prPred)
    dpred = {}
    for filePred in lfilePred:
        print(filePred)
        dtemp = toolbox.loadMatrixToDict(prPred + filePred, sep=",")

        for DTXCID in dtemp.keys():
            DTXSID = dtemp[DTXCID]["dsstox_substance_id"]
            if not DTXSID in list(dpred.keys()):
                dpred[DTXSID] = dtemp[DTXCID]

    k = list(dpred.keys())[0]

    return



def UpdateDBName():

    return


# define file where DB update are included

ptrainTestOPERA = "/home/borrela2/data/updateDB_01-2020/ToxTrainTest_3D.sdf" #file with endpoint
ptrainTestOPERA = "C:\\Users\\borrela2\\development\\trash\\ToxTrainTest_3D.sdf" #file with endpoint

pIdentifier = "/home/borrela2/data/updateDB_01-2020/DSSTox_Identifiers_and_CASRN_04-2019_corrected.csv" #identification with name
pIdentifier = "/home/borrela2/data/updateDB_01-2020/DSSTox_Identifiers_and_CASRN_04-2019.csv" #identification with name
pIdentifier = "C:\\Users\\borrela2\\development\\trash\\DSSTox_Identifiers_and_CASRN_04-2019_corrected.csv"

prOPERAPred = "/home/borrela2/data/updateDB_01-2020/DSSTox_opera_01-2020/" #folder included prediction for OPERA
prOPERAPred = "C:\\Users\\borrela2\\development\\trash\\DSSTox_opera_01-2020\\"

prout = "/home/borrela2/ChemMaps/update_server/" #out directory
prout = "C:\\Users\\borrela2\\development\\trash\\"

LPROPCHEMAPS = ["inchikey", "SMILES", "preferred_name", "GHS_category", "EPA_category", "consensus_LD50", "LD50_mgkg", "MolWeight", "LogOH_pred", "CATMoS_VT_pred", "CATMoS_NT_pred", "CATMoS_EPA_pred",
                 "CATMoS_GHS_pred", "CATMoS_LD50_pred", "CERAPP_Ago_pred", "CERAPP_Anta_pred", "CERAPP_Bind_pred", "Clint_pred", "CoMPARA_Ago_pred", "CoMPARA_Anta_pred",
                 "CoMPARA_Bind_pred", "FUB_pred", "LogHL_pred", "LogKM_pred", "LogKOA_pred", "LogKoc_pred", "LogBCF_pred", "LogD55_pred", "LogP_pred", "MP_pred", "pKa_a_pred",
                 "pKa_b_pred", "ReadyBiodeg_pred", "RT_pred", "LogVP_pred", "LogWS_pred", "BioDeg_LogHalfLife_pred", "BP_pred", "nbLipinskiFailures"]

UpdateDBVal(prOPERAPred, pIdentifier, ptrainTestOPERA, prout)

