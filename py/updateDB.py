import DBrequest
import toolbox
import parseSDF
from random import randint
from os import listdir


def UpdateDBVal(prPred, pIdentifier, pSDF, prout):

    dIdentifier = toolbox.loadMatrixToDict(pIdentifier)
    sss
    csdf = parseSDF.parseSDF(pSDF, "InChI Key_QSARr", prout)
    csdf.parseAll()

    # load chem prediction
    lfilePred = listdir(prPred)
    dpred = {}
    for filePred in lfilePred:
        
        print(filePred)


    return



def UpdateDBName():
    
    



    
    return


# define file where DB update are included

ptrainTestOPERA = "/home/borrela2/data/updateDB_01-2020/ToxTrainTest_3D.sdf" #file with endpoint
pIdentifier = "/home/borrela2/data/updateDB_01-2020/DSSTox_Identifiers_and_CASRN_04-2019_corrected.csv" #identification with name
pIdentifier = "/home/borrela2/data/updateDB_01-2020/DSSTox_Identifiers_and_CASRN_04-2019.csv" #identification with name

prOPERAPred = "/home/borrela2/data/updateDB_01-2020/DSSTox_opera_01-2020/" #folder included prediction for OPERA
prout = "/home/borrela2/ChemMaps/update_server/" #out directory

LPROP = ["inchikey", "SMILES", "preferred_name", "GHS_category", "EPA_category", "consensus_LD50", "LD50_mgkg", "MolWeight", "LogOH_pred", "CATMoS_VT_pred", "CATMoS_NT_pred", "CATMoS_EPA_pred",
                 "CATMoS_GHS_pred", "CATMoS_LD50_pred", "CERAPP_Ago_pred", "CERAPP_Anta_pred", "CERAPP_Bind_pred", "Clint_pred", "CoMPARA_Ago_pred", "CoMPARA_Anta_pred",
                 "CoMPARA_Bind_pred", "FUB_pred", "LogHL_pred", "LogKM_pred", "LogKOA_pred", "LogKoc_pred", "LogBCF_pred", "LogD55_pred", "LogP_pred", "MP_pred", "pKa_a_pred",
                 "pKa_b_pred", "ReadyBiodeg_pred", "RT_pred", "LogVP_pred", "LogWS_pred", "BioDeg_LogHalfLife_pred", "BP_pred", "nbLipinskiFailures"]



UpdateDBVal(prOPERAPred, pIdentifier, ptrainTestOPERA, prout)

