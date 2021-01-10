from os import system, path
from shutil import copy
from copy import deepcopy


def loadMatrixToList(pmatrixIn, sep = "\t"):

    filin = open(pmatrixIn, "r", encoding="utf8", errors='ignore')
    llinesMat = filin.readlines()
    filin.close()

    l_out = []
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1)-1) and lval1[-1] != "":
        lheaders.append("ID")


    i = 1
    imax = len(llinesMat)
    while i < imax:
        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        j = 0
        if len(lvalues) != len(lheaders):
            print("Check different size - line: ", i)
            #print(lvalues)
            #print(lheaders)
        jmax = len(lheaders)
        dtemp = {}
        while j < jmax:
            try:dtemp[lheaders[j]] = lvalues[j]
            except:pass
            j += 1
        l_out.append(dtemp)
        i += 1

    return l_out

def loadMatrixCoords(pccord, nbcoord):

    dout = {}
    filin = open(pccord, "r")
    line = filin.readline()# header
    line = filin.readline()
    while line:
        lelem = line.split(",")
        i = 1
        dout[lelem[0].replace("\"", "")] = []
        while i <= nbcoord:
            dout[lelem[0].replace("\"", "")].append(float(lelem[i])) 
            i = i + 1
        line = filin.readline()
    return dout

def loadMatrixTolistFromDB(pfilin, sep):
    lout = []
    filin = open(pfilin, "r")
    llines = filin.readlines()
    filin.close()

    for l in llines:
        ltemp = l.split(",")
        val1 = ltemp[0]
        val2 = ltemp[1]
        val3 = l.split(",", 2)[-1].strip()
        llout = [val1, val2, val3]
        lout.append(llout)
    
    return lout

def loadMatrixToDict(pmatrixIn, sep ="\t"):

    filin = open(pmatrixIn, "r", encoding="utf-8", errors="ignore")
    llinesMat = filin.readlines()
    filin.close()

    dout = {}
    line0 = formatLine(llinesMat[0])
    line1 = formatLine(llinesMat[1])
    lheaders = line0.split(sep)
    lval1 = line1.split(sep)

    # case where R written
    if len(lheaders) == (len(lval1)-1) and lval1[-1] != "":
        lheaders.append("val")

    i = 1
    imax = len(llinesMat)
    while i < imax:
        lineMat = formatLine(llinesMat[i])
        lvalues = lineMat.split(sep)
        kin = lvalues[0]
        dout[kin] = {}
        j = 0
        jmax = len(lheaders)
        while j < jmax:
            dout[kin][lheaders[j]] = lvalues[j]
            j += 1
        i += 1
    return dout

def formatLine(lineinput, delimitorStr = "\""):

    linein = deepcopy(lineinput)
    linein = linein.replace("\n", "")

    linenew = ""

    imax = len(linein)
    i = 0
    flagchar = 0
    while i < imax:
        if linein[i] == delimitorStr and flagchar == 0:
            flagchar = 1
        elif linein[i] == delimitorStr and flagchar == 1:
            flagchar = 0

        if flagchar == 1 and linein[i] == ",":
            linenew = linenew + " "
        else:
            linenew = linenew + linein[i]
        i += 1

    linenew = linenew.replace(delimitorStr, "")
    return linenew


import multiprocessing
import time

def timeFunction(funct, mol):

    manager = multiprocessing.Manager()
    lout = manager.list()

    p = multiprocessing.Process(target=funct, args=(mol, lout))
    p.start()
    time.sleep(2)

    if p.is_alive():
        p.terminate()
        p.join()
        return "ERROR"
    else:
        p.join()
        #print lout
        return lout[0]

def selectMinimalEnergyLigPrep(psdfin, psdfout):

    # case of only one
    filin = open(psdfin, "r")
    readfile = filin.read()
    filin.close()

    lsdf = readfile.split("$$$$\n")[:-1]


    if len(lsdf) == 1:
        copy(psdfin, psdfout)

    else:
        #find with the lower energy
        lenergy = []
        for sdfin in lsdf:
            energy = sdfin.split("> <r_lp_Energy>\n")[-1].split("\n")[0]
            print (energy)
            lenergy.append(float(energy))

        # take minimal energy
        ibest = lenergy.index(min(lenergy))
        print (ibest)
        filout = open(psdfout, "w")
        filout.write(lsdf[ibest] + "$$$$\n")
        filout.close()

    return psdfout

def parsePadelOut(pfiledesc=""):
    """
    Only case of 2 lines in descriptors file
    :param pfiledesc:
    :return:
    """

    ddesc = {"RDF65s": "NA", "RDF65p": "NA", "RDF65v": "NA", "RDF65u": "NA", "RDF145v": "NA", "RDF65e": "NA",
            "RDF65i": "NA"
        , "RDF65m": "NA", "RDF80v": "NA", "RDF80u": "NA", "RDF80s": "NA", "RDF80p": "NA", "RDF80m": "NA",
            "RDF80i": "NA", "RDF80e"
            : "NA", "RDF145e": "NA", "LOBMAX": "NA", "TDB3e": "NA", "TDB3i": "NA", "TDB3m": "NA", "TDB3r": "NA",
            "TDB3s": "NA",
            "TDB3p": "NA", "TDB3v": "NA", "TDB3u": "NA", "TDB2e": "NA", "RDF105e": "NA", "RDF105i": "NA",
            "RDF105m": "NA",
            "RDF105s": "NA", "RDF105p": "NA", "RDF105v": "NA", "RDF105u": "NA", "RDF15e": "NA", "geomRadius": "NA",
            "RDF15m": "NA", "RDF15i": "NA", "RDF15u": "NA", "RDF15v": "NA", "RDF15p": "NA", "RDF15s": "NA", "P1e": "NA",
            "RDF25v": "NA", "RDF25u": "NA", "RDF25s": "NA", "P1m": "NA", "RDF25p": "NA", "P1p": "NA", "RDF25m": "NA",
            "P1s": "NA",
            "P1u": "NA", "RDF25i": "NA", "RDF25e": "NA", "E2u": "NA", "RDF130e": "NA", "RPCS": "NA", "E2s": "NA",
            "RDF130m": "NA",
            "RDF130i": "NA", "RPCG": "NA", "RDF130v": "NA", "RDF130u": "NA", "RDF130s": "NA", "RDF130p": "NA",
            "E2m": "NA",
            "E2i": "NA", "TDB10p": "NA", "TDB10s": "NA", "TDB10r": "NA", "TDB10u": "NA", "TDB10v": "NA", "L3m": "NA",
            "L3i": "NA",
            "L3v": "NA", "L3u": "NA", "TDB10e": "NA", "L3s": "NA", "L3p": "NA", "TDB10i": "NA", "TDB10m": "NA",
            "DPSA-1": "NA",
            "DPSA-3": "NA", "DPSA-2": "NA", "WPSA-1": "NA", "WPSA-2": "NA", "WPSA-3": "NA", "TDB8v": "NA",
            "TDB8u": "NA",
            "TDB8s": "NA", "TDB8r": "NA", "TDB8p": "NA", "L3e": "NA", "TDB8e": "NA", "TDB8m": "NA", "TDB8i": "NA",
            "RDF100u": "NA",
            "RDF100v": "NA", "RDF100p": "NA", "geomShape": "NA", "RDF100s": "NA", "RDF100m": "NA", "RDF100i": "NA",
            "RDF100e": "NA", "RDF10v": "NA", "RDF10u": "NA", "RDF10s": "NA", "RDF10p": "NA", "RDF10m": "NA",
            "RDF10i": "NA",
            "RDF10e": "NA", "MOMI-XZ": "NA", "RDF40s": "NA", "RDF40p": "NA", "RDF40v": "NA", "RDF40u": "NA",
            "RDF40e": "NA",
            "RDF40i": "NA", "RDF40m": "NA", "RDF150p": "NA", "RDF150s": "NA", "RDF150u": "NA", "RDF150v": "NA",
            "RDF150e": "NA",
            "RDF150i": "NA", "RDF150m": "NA", "TDB4s": "NA", "TDB4r": "NA", "TDB4p": "NA", "TDB4v": "NA", "TDB4u": "NA",
            "TDB4i": "NA", "TDB4m": "NA", "TDB4e": "NA", "TPSA": "NA", "PPSA-2": "NA", "RDF115m": "NA", "P2e": "NA",
            "P2m": "NA",
            "P2i": "NA", "P2u": "NA", "P2v": "NA", "P2p": "NA", "P2s": "NA", "RDF70i": "NA", "RDF70m": "NA",
            "RDF70e": "NA",
            "RDF70p": "NA", "RDF70s": "NA", "RDF70u": "NA", "RDF70v": "NA", "WNSA-2": "NA", "WNSA-3": "NA",
            "WNSA-1": "NA",
            "RDF125m": "NA", "TDB9p": "NA", "TDB9r": "NA", "TDB9s": "NA", "TDB9u": "NA", "TDB9v": "NA", "TDB9e": "NA",
            "TDB9i": "NA", "TDB9m": "NA", "RHSA": "NA", "Ke": "NA", "Ki": "NA", "Km": "NA", "Ks": "NA", "Kp": "NA",
            "Kv": "NA",
            "Ku": "NA", "L1s": "NA", "TDB6v": "NA", "Vm": "NA", "TDB6r": "NA", "Vp": "NA", "Dm": "NA", "Di": "NA",
            "De": "NA",
            "Dv": "NA", "TDB6i": "NA", "Du": "NA", "Ds": "NA", "Dp": "NA", "RDF45m": "NA", "RDF35s": "NA",
            "RDF35p": "NA",
            "RDF45i": "NA", "RDF35u": "NA", "RDF45e": "NA", "E1p": "NA", "E1s": "NA", "E1u": "NA", "E1v": "NA",
            "E1i": "NA",
            "E1m": "NA", "RDF35e": "NA", "RDF45u": "NA", "RDF45v": "NA", "E1e": "NA", "RDF45p": "NA", "RDF45s": "NA",
            "RDF35m": "NA", "MOMI-YZ": "NA", "RDF85e": "NA", "RDF85i": "NA", "RDF85m": "NA", "RDF85p": "NA",
            "RDF85s": "NA"
        , "RDF85u": "NA", "RDF85v": "NA", "RDF110u": "NA", "TDB5e": "NA", "RDF145u": "NA", "RDF110v": "NA",
            "RDF145s": "NA",
            "RDF110p": "NA", "RDF110s": "NA", "RDF145p": "NA", "TDB5m": "NA", "TDB5i": "NA", "RDF110e": "NA",
            "TDB5u": "NA",
            "TDB5v": "NA", "TDB5p": "NA", "TDB5r": "NA", "TDB5s": "NA", "RDF110m": "NA", "RDF145m": "NA",
            "RDF110i": "NA",
            "RDF145i": "NA", "Tp": "NA", "RDF30m": "NA", "RDF30i": "NA", "RDF55i": "NA", "RDF30e": "NA",
            "geomDiameter": "NA",
            "RDF30u": "NA", "RDF30v": "NA", "RDF30p": "NA", "RDF30s": "NA", "RDF55m": "NA", "RDF135i": "NA",
            "PPSA-1": "NA",
            "RDF135m": "NA", "L1u": "NA", "L1v": "NA", "L1p": "NA", "RDF55e": "NA", "RDF135e": "NA", "L1m": "NA",
            "L1i": "NA",
            "RDF135p": "NA", "L1e": "NA", "RDF135s": "NA", "RDF55u": "NA", "RDF55v": "NA", "RDF95m": "NA",
            "RDF95i": "NA",
            "RDF95e": "NA", "FPSA-3": "NA", "FPSA-2": "NA", "FPSA-1": "NA", "RDF95u": "NA", "RDF95v": "NA",
            "RDF95p": "NA",
            "RDF95s": "NA", "TDB1i": "NA", "TDB1m": "NA", "TDB1e": "NA", "RDF55p": "NA", "TDB1p": "NA", "TDB1r": "NA",
            "TDB1s": "NA", "TDB1u": "NA", "TDB1v": "NA", "RDF55s": "NA", "RDF135u": "NA", "Tv": "NA", "RDF135v": "NA",
            "Tu": "NA",
            "RDF60e": "NA", "RDF60i": "NA", "E2v": "NA", "RDF60m": "NA", "Ts": "NA", "RDF60p": "NA", "RDF60s": "NA",
            "RDF60u": "NA", "RDF60v": "NA", "RPSA": "NA", "MOMI-Y": "NA", "Tm": "NA", "P1i": "NA", "PNSA-3": "NA",
            "PNSA-2": "NA",
            "PNSA-1": "NA", "MOMI-XY": "NA", "RDF90e": "NA", "RDF90m": "NA", "RDF90i": "NA", "RDF90v": "NA",
            "RDF90u": "NA",
            "RDF90s": "NA", "RDF90p": "NA", "RDF125u": "NA", "RDF125v": "NA", "RDF125p": "NA", "RDF125s": "NA",
            "RDF115i": "NA",
            "RDF115e": "NA", "RDF125e": "NA", "RDF115v": "NA", "PPSA-3": "NA", "RDF115u": "NA", "RDF125i": "NA",
            "RDF115s": "NA",
            "RDF115p": "NA", "Ve": "NA", "RDF155e": "NA", "TDB6u": "NA", "Vi": "NA", "RDF155i": "NA", "TDB6p": "NA",
            "TDB6s": "NA",
            "RDF155m": "NA", "TDB6m": "NA", "RDF155s": "NA", "RDF155p": "NA", "Vs": "NA", "RDF155v": "NA", "Vu": "NA",
            "Vv": "NA",
            "RDF155u": "NA", "TDB6e": "NA", "GRAV-6": "NA", "GRAV-5": "NA", "GRAV-4": "NA", "GRAV-3": "NA",
            "GRAV-2": "NA",
            "GRAV-1": "NA", "E2p": "NA", "LOBMIN": "NA", "E2e": "NA", "RDF120m": "NA", "RDF120i": "NA", "RDF120e": "NA",
            "RDF120v": "NA", "RDF35v": "NA", "RDF120u": "NA", "RDF120s": "NA", "RDF120p": "NA", "FNSA-1": "NA",
            "FNSA-3": "NA", "FNSA-2": "NA", "TDB2i": "NA", "TDB2m": "NA", "TDB2p": "NA", "TDB2s": "NA", "TDB2r": "NA",
            "TDB2u": "NA", "TDB2v": "NA", "THSA": "NA", "RDF20e": "NA", "RDF75e": "NA", "RDF20m": "NA", "RDF75m": "NA",
            "RDF20i": "NA", "RDF75i": "NA", "RDF75v": "NA", "RDF20u": "NA", "RDF20v": "NA", "RDF75u": "NA",
            "RDF20p": "NA",
            "RDF75s": "NA", "RDF75p": "NA", "RDF20s": "NA", "P1v": "NA", "RDF35i": "NA", "E3v": "NA", "E3u": "NA",
            "E3s": "NA",
            "E3p": "NA", "E3e": "NA", "GRAVH-2": "NA", "GRAVH-3": "NA", "E3m": "NA", "GRAVH-1": "NA", "E3i": "NA",
            "L2e": "NA",
            "L2i": "NA", "L2m": "NA", "L2p": "NA", "L2s": "NA", "L2u": "NA", "L2v": "NA", "RDF50s": "NA",
            "RDF50p": "NA",
            "RDF50v": "NA", "RDF50u": "NA", "MOMI-Z": "NA", "RDF50i": "NA", "MOMI-X": "NA", "RDF50m": "NA", "Ti": "NA",
            "MOMI-R": "NA", "Te": "NA", "RDF50e": "NA", "TDB7v": "NA", "TDB7u": "NA", "TDB7r": "NA", "TDB7s": "NA",
            "TDB7p": "NA",
            "TDB7m": "NA", "TDB7i": "NA", "TDB7e": "NA", "RNCG": "NA", "RDF140p": "NA", "RDF140s": "NA",
            "RDF140u": "NA",
            "RDF140v": "NA", "RDF140i": "NA", "RDF140m": "NA", "RNCS": "NA", "RDF140e": "NA", "Ae": "NA", "Ai": "NA",
            "Am": "NA",
            "Ap": "NA", "As": "NA", "Au": "NA", "Av": "NA"}

    if pfiledesc == "":
        return ddesc
    else:
        filin = open(pfiledesc, "r")
        lline = filin.readlines()
        filin.close()

        dout = {}

        lheaders = lline[0].strip().split(",")

        i = 1
        while i < len(lline):
            lvalues = lline[i].strip().split(",")
            #print lvalues
            j = 1
            dtemp = deepcopy(ddesc)
            while j < len(lheaders):
                dtemp[lheaders[j]] = lvalues[j]
                j += 1
            dout[lvalues[0]] = dtemp
            i += 1
    return dout

def renameHeaderSDF(pfilin):
    """Rename header with name file"""
    namesdf = pfilin.split("/")[-1].split(".")[0]
    filin = open(pfilin, "r")
    llines = filin.readlines()
    filin.close()
    llines[0] = str(namesdf) + "\n"

    filout = open(pfilin, "w")
    filout.write("".join(llines))
    filout.close()


from rdkit import Chem
def convertSMILEStoINCHIKEY(SMILESin):
    molformat = Chem.MolFromSmiles(SMILESin)
    inchi = Chem.inchi.MolToInchi(molformat)
    inchikey = Chem.inchi.InchiToInchiKey(inchi)

    return inchikey

# fast search
from bisect import bisect_left

def binary_search(L, x):
    i = bisect_left(L, x)
    if i == len(L) or L[i] != x:
        return -1
    return i
