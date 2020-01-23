import DBrequest
import toolbox
import parseSDF
from random import randint
from os import listdir


def UpdateDBVal(prPred, pSMILES, pSDF, prout):

    csdf = parseSDF.parseSDF(pSDF, "InChI Key_QSARr", prout)
    csdf.parseAll()

    dSDF = {}
    for chem in csdf.lc:
        SMILES = chem["Original_SMILES"]
        dSDF[SMILES] = chem


    dSMILES = toolbox.loadMatrixToDict(pSMILES, ",")
    dSMILES_out = {}
    for ID in dSMILES.keys():
        dSMILES_out[dSMILES[ID]["DTXSID"]] = dSMILES[ID]

    # load chem prediction
    lfilePred = listdir(prPred)
    dpred = {}
    for filePred in lfilePred[:1]:
        print(filePred)
        dtemp = toolbox.loadMatrixToDict(prPred + filePred, sep=",")

        for DTXCID in dtemp.keys():
            DTXSID = dtemp[DTXCID]["dsstox_substance_id"]
            if not DTXSID in list(dpred.keys()):
                dpred[DTXSID] = dtemp[DTXCID]

    k = list(dpred.keys())[0]
    print(len(dpred.keys()))

    # write for DB Update
    pfiloutDesc = prout + "DescDescription.csv"
    filoutDesc = open(pfiloutDesc, "w")
    filoutDesc.write("ID\t" + "\t".join(LPROPCHEMAPS) + "\n")

    for DTXSID in dpred:
        print(DTXSID)
        try:
            inch = dSMILES_out[DTXSID]["INCHIKEY"]
            SMILES = dSMILES_out[DTXSID]["SMILES"]
            name = dSMILES_out[DTXSID]["name"]
        except:
            continue

        if SMILES in list(dSDF.keys()):
            GHS_category = dSDF[SMILES]["GHS_category"]
            if GHS_category == "":
                GHS_category = "NA"
            EPA_category = dSDF[SMILES]["EPA_category"]
            if EPA_category == "":
                EPA_category = "NA"
            consensus_LD50 = dSDF[SMILES]["consensus_LD50"]
            if consensus_LD50 == "":
                consensus_LD50 = "NA" 
            LD50_mgkg = dSDF[SMILES]["LD50_mgkg"]
            if LD50_mgkg == "":
                LD50_mgkg = "NA"
        else:
            GHS_category = "NA"
            LD50_mgkg = "NA"
            consensus_LD50 = "NA"
            EPA_category = "NA"

        lw = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(DTXSID, inch, SMILES, name, GHS_category, EPA_category, consensus_LD50, LD50_mgkg, "\t".join([str(dpred[DTXSID][PROP]) for PROP in LPROPCHEMAPS[7:]]))
        filoutDesc.write(lw)

    filoutDesc.close()




def UpdateDBName():

    return


# define file where DB update are included

ptrainTestOPERA = "/home/borrela2/data/updateDB_01-2020/ToxTrainTest_3D.sdf" #file with endpoint
#ptrainTestOPERA = "C:\\Users\\borrela2\\development\\trash\\ToxTrainTest_3D.sdf" #file with endpoint

pIdentifier = "/home/borrela2/data/updateDB_01-2020/DSSToxMS-Ready_corrected.csv" #identification with name
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

LPROPOPERA = []

UpdateDBVal(prOPERAPred, pIdentifier, ptrainTestOPERA, prout)

