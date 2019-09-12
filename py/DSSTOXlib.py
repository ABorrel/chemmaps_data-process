import pathFolder
import toolbox
import runExternalSoft
import DBrequest
import calculate

#=> to import ToxCast librairy
import sys
sys.path.insert(0, "/home/borrela2/development/descriptor/")
#sys.path.insert(0, "C:\\Users\\borrela2\\development\\molecular-descriptors\\")
import Chemical



from os import path, listdir, remove
from copy import deepcopy
from math import sqrt
from rdkit import Chem
from re import search
from random import shuffle


class DSSTOX:

    def __init__(self, plistChem, istart, iend, prDesc, prout):

        self.plistChem = plistChem
        self.istart = istart
        self.iend = iend
        self.plog = prout + "log.txt"
        self.prDesc = prDesc
        self.prout = prout



    def loadlistChem(self):

        prForDB = pathFolder.createFolder(self.prout + "forDB/")
        pfilout = prForDB + "db.csv"
        # try:remove(pfilout)
        # except:pass
        # print(pfilout)
        if path.exists(pfilout):
            dchem = toolbox.loadMatrixToDict(pfilout, sep="\t")
        else:
            dchem = {}
            dchemIn = toolbox.loadMatrixToDict(self.plistChem, sep=",") #rewrtie pfas and tox21 with comma
            for chemIn in dchemIn.keys():
                if "SMILES" in list(dchemIn[chemIn].keys()):
                    SMILES_origin = dchemIn[chemIn]["SMILES"]
                    DTXSID = dchemIn[chemIn]["DTXSID"]
                elif "Original_SMILES" in list(dchemIn[chemIn].keys()):
                    SMILES_origin = dchemIn[chemIn]["Original_SMILES"]
                    DTXSID = dchemIn[chemIn]["dsstox_substance_id"]
                else:
                    print("ERROR")
                    return


                dchem[DTXSID] = {}
                dchem[DTXSID]["DTXSID"] = DTXSID
                dchem[DTXSID]["smiles_origin"] = SMILES_origin

                # prepare ligand
                cchem = Chemical.Chemical(SMILES_origin, self.prDesc)
                cchem.prepChem()
                if cchem.err == 1:
                    qsar_ready = 0
                    cleanSMILES = "NA"
                    inchikey = "NA"
                else:
                    qsar_ready = 1
                    cleanSMILES = cchem.smi
                    inchikey = cchem.generateInchiKey()
                    cchem.writeSMIClean()

                dchem[DTXSID]["smiles_clean"] = cleanSMILES
                dchem[DTXSID]["inchikey"] = inchikey
                dchem[DTXSID]["qsar_ready"] = qsar_ready


            # write table for control -> after open and put in the DB
            filout = open(pfilout, "w", encoding="utf8")
            filout.write("DTXSID\tsmiles_origin\tsmiles_clean\tinchikey\tqsar_ready\n")
            for chem in dchem.keys():
                filout.write("%s\t%s\t%s\t%s\t%s\n" % (
                chem, dchem[chem]["smiles_origin"], dchem[chem]["smiles_clean"], dchem[chem]["inchikey"],
                dchem[chem]["qsar_ready"]))
            filout.close()
        self.dchem = dchem


    def computeDesc (self, insertDB =0, w=0):

        if not "dchem" in self.__dict__:
            self.loadlistChem()

        lchemID = list(self.dchem.keys()) # can be shuffle
        if len(lchemID) < 50000 and self.iend == 0:
            shuffle(lchemID)

        imax = len(lchemID)
        if self.iend == 0 or self.iend > imax:
            iend = imax
        else:
            iend = self.iend

        lchemID = lchemID[self.istart:iend]
        shuffle(lchemID)

        if insertDB == 1:
            cDB = DBrequest.DBrequest()

        pfilout1D2D = self.prout + "1D2D.csv"
        pfilout3D = self.prout + "3D.csv"

        if w == 1:

            ldesc1D2D = Chemical.getLdesc("1D2D")
            ldesc3D = Chemical.getLdesc("3D")

            filout1D2D = open(pfilout1D2D, "w")
            filout1D2D.write("inchikey\t" + "\t".join(ldesc1D2D) + "\n")

            filout3D = open(pfilout3D, "w")
            filout3D.write("inchikey\t" + "\t".join(ldesc3D) + "\n")

        i = 0
        imax = len(lchemID)
        while i < imax:
            if i%1000 == 0:
                print (i)

            SMILESClean = self.dchem[lchemID[i]]["smiles_clean"]
            if SMILESClean == "NA":
                i = i + 1
                continue
            cChem = Chemical.Chemical(SMILESClean, self.prDesc)
            cChem.prepChem()

            # print(SMILESClean)
            if cChem.err == 0:
                # 2D descriptors
                cChem.computeAll2D(update=0)
                if cChem.err == 1:
                    i = i + 1
                    continue
                cChem.writeMatrix("2D")
                if w == 1:
                    filout1D2D.write(
                        "%s\t%s\n" % (cChem.inchikey, "\t".join([str(cChem.all2D[desc]) for desc in ldesc1D2D])))

                # insert in DB
                if insertDB == 1:
                    cDB.verbose = 0
                    out1D2D = cDB.getRow("desc_1d2d", "inchikey='%s'" % (cChem.inchikey))
                    if out1D2D == []:
                        w1D2D = "{" + ",".join(["\"%s\"" % (cChem.all2D[desc1D2D]) for desc1D2D in ldesc1D2D]) + "}"
                        cDB.addElement("desc_1d2d", ["inchikey", "desc_value"], [cChem.inchikey, w1D2D])

                # 3D descriptors
                cChem.set3DChemical()
                # control if 3D generated
                if cChem.err == 0:
                    cChem.computeAll3D(update=0)
                    cChem.writeMatrix("3D")

                    if cChem.err == 1:
                        i = i + 1
                        continue
                    # write master table
                    if w == 1:
                        filout3D.write(
                            "%s\t%s\n" % (cChem.inchikey, "\t".join([str(cChem.all3D[desc]) for desc in ldesc3D])))

                    # put in table descriptors
                    if insertDB == 1:
                        out3D = cDB.getRow("desc_3d", "inchikey='%s'" % (cChem.inchikey))
                        if out3D == []:
                            w3D = "{" + ",".join(["\"%s\"" % (cChem.all3D[desc3D]) for desc3D in ldesc3D]) + "}"
                            cDB.addElement("desc_3d", ["inchikey", "desc_value"], [cChem.inchikey, w3D])
            i = i + 1

        if w == 1:
            filout1D2D.close()
            filout3D.close()

        self.p1D2D = pfilout1D2D
        self.p3D = pfilout3D





    def computepng(self):

        # draw from desc descriptors

        prSMI = self.prDSSTox + "SMI/"
        if not path.exists(prSMI):
            print ("ERROR: CREATE CLEAN SMI FIRST")
            return

        prPNG = self.prDSSTox + "PNG/"
        pathFolder.createFolder(prPNG)

        if "ddesc" in self.__dict__:
            linchikey = list(self.ddesc.keys())
            shuffle(linchikey)
            for inchikey in linchikey:
                ppng = prPNG + inchikey + ".png"
                if path.exists(ppng):
                    continue
                else:
                    DSSTOXid = self.ddesc[inchikey]["DSSTOXid"][0]
                    SMIclean = self.ddesc[inchikey]["SMI"][0]

                    pSMIclean = prPNG + DSSTOXid + ".smi"
                    fSMIcLean = open(pSMIclean, "w")
                    fSMIcLean.write(SMIclean)
                    fSMIcLean.close()

                    runExternalSoft.molconvert(pSMIclean, ppng)
                    remove(pSMIclean)




    def computeCoords(self, corVal, distributionVal, insertDB=1):


        if not "p1D2D" in self.__dict__ and not "p3D" in self.__dict__:
            self.computeDesc(insertDB=0, w=1)

        # create coords
        prmap = pathFolder.createFolder(self.prout + "map_" + str(corVal) + "-" + str(distributionVal) + "/")
        self.prmap = prmap

        pcoordDim1Dim2 = prmap + "coord1D2D.csv"
        pcoordDim3D = prmap + "coord3D.csv"
        if not path.exists(pcoordDim1Dim2) or not path.exists(pcoordDim3D):
            runExternalSoft.RComputeMapFiles(self.p1D2D, self.p3D, prmap, corVal, distributionVal)

        if not path.exists(pcoordDim1Dim2) or not path.exists(pcoordDim3D):
            print("ERROR file map")
            return
        else:
            self.pcoords1D2D = pcoordDim1Dim2
            self.pcoords3D = pcoordDim3D

        if insertDB == 1:
            dcoord1D2D = toolbox.loadMatrixToDict(pcoordDim1Dim2, sep = ",")
            dcoord3D = toolbox.loadMatrixToDict(pcoordDim3D, sep = ",")
            cDB = DBrequest.DBrequest()
            cDB.verbose = 0
            for chem in dcoord1D2D.keys():
                print(chem)
                out1D2D = cDB.getRow("drugbank_coords", "inchikey='%s'" % (chem))
                if out1D2D == []:
                    nbdim1d2d = len(dcoord1D2D[chem].keys()) - 1
                    nbdim3d = len(dcoord3D[chem].keys()) - 1

                    w1D2D = "{" + ",".join(["\"%s\"" % (dcoord1D2D[chem]["DIM" + str(i)]) for i in range(1, nbdim1d2d + 1)]) + "}"
                    w3D = "{" + ",".join(["\"%s\"" % (dcoord3D[chem]["DIM3-" + str(i)]) for i in range(1, nbdim3d + 1)]) + "}"
                    cDB.addElement("drugbank_coords", ["inchikey", "dim1d2d", "dim3d", "origin"], [chem, w1D2D, w3D, "1"])


    def runRprojection(self, corVal, distributionVal):

        if not "p1D2D" in self.__dict__ and not "p3D" in self.__dict__:
            self.computeDesc(insertDB=0)

        # create coords
        prproj = pathFolder.createFolder(self.prout + "proj_" + str(corVal) + "-" + str(distributionVal) + "/")
        runExternalSoft.RComputeCor(self.p1D2D, self.p3D, prproj, corVal, distributionVal)




    def splitMap(self, nbsplit, dim):

        if not "prmap" in self.__dict__:
            print ("Generate the map files first")
            return

        else:
            prout = pathFolder.createFolder(self.prmap + "split_" + str(nbsplit) + "_" + str(dim) + "/")
            self.prmaps = prout

            if len(listdir(self.prmaps)) >= (nbsplit*2):
                print ("Maps already computed")
                return


            coord1D2D = self.prmap + "coord1D2D.csv"
            coord3D = self.prmap + "coord3D.csv"

            d1D2D = toolbox.loadMatrixToDict(coord1D2D, sep = ",")
            d3D = toolbox.loadMatrixToDict(coord3D, sep=",")

            # max and min 1D2D
            maxDim = 0.0
            minDim = 0.0

            nbchem = len(list(d1D2D.keys()))
            nbchembymap = int(nbchem/nbsplit)

            # calibrate max and min
            for chem in d1D2D.keys():

                if dim == 1:
                    dimVal = float(d1D2D[chem]["DIM1"])
                elif dim == 2:
                    dimVal = float(d1D2D[chem]["DIM2"])
                elif dim == 3:
                    dimVal = float(d3D[chem]["DIM3"])

                if dimVal > maxDim:
                    maxDim = dimVal
                if dimVal < minDim:
                    minDim = dimVal


            dmap = {}
            imap = 1
            dmap[imap] = []

            dimVal = minDim
            while dimVal < maxDim:
                dimVal = dimVal + 0.10
                if len(dmap[imap]) > nbchembymap:
                    imap = imap + 1
                    dmap[imap] = []
                ichem = 0
                lchem = list(d1D2D.keys())
                nbchem = len(lchem)
                while ichem < nbchem:

                    if dim == 1:
                        valtemp = float(d1D2D[lchem[ichem]]["DIM1"])
                    elif dim == 2:
                        valtemp = float(d1D2D[lchem[ichem]]["DIM2"])
                    elif dim == 3:
                        valtemp = float(d3D[lchem[ichem]]["DIM3"])


                    if valtemp < dimVal:
                        dmap[imap].append(deepcopy(d1D2D[lchem[ichem]]))
                        del d1D2D[lchem[ichem]]
                        del lchem[ichem]
                        nbchem = nbchem - 1
                        continue
                    else:
                        ichem = ichem + 1


        for d in dmap.keys():
            pfilout1D2D = prout + str(d) + "_map1D2D.csv"
            pfilout3D = prout + str(d) + "_map3D.csv"

            filout1D2D = open(pfilout1D2D, "w")
            filout3D = open(pfilout3D, "w")
            filout1D2D.write("\"\",\"DIM1\",\"DIM2\"\n")
            filout3D.write("\"\",\"DIM3\",\"DIM4\"\n")

            for chem in dmap[d]:
                filout1D2D.write("\"%s\",%s,%s\n"%(chem[""], chem["DIM1"], chem["DIM2"]))
                filout3D.write("\"%s\",%s,%s\n" % (chem[""], d3D[chem[""]]["DIM3"], d3D[chem[""]]["DIM4"]))

            filout1D2D.close()
            filout3D.close()

            print (len(dmap[d]), d)


    def generateMapSplitFile(self, dim):

        if not "prmaps" in self.__dict__:
            print ("Generate Maps first")
            return

        else:
            lfmap = listdir(self.prmaps)
            dCentroid = {}
            lmap = []

            pfileMapChem = self.prmap + str(dim) + "_mapChem.csv"
            pfmapCentroid = self.prmap + str(dim) + "MapCentroid.csv"
            if path.exists(pfileMapChem) and path.exists(pfmapCentroid):
                print ("Already computed")
                return


            fmapChem = open(pfileMapChem, "w")
            fmapChem.write("ID\tMap\n")

            for namefmap in lfmap:
                map = namefmap.split("_")[0]

                if map in lfmap or map == "":
                    continue
                else:
                    pf1D2D = self.prmaps + map + "_map1D2D.csv"
                    pf3D = self.prmaps + map + "_map3D.csv"

                    d1D2D = toolbox.loadMatrixToDict(pf1D2D, sep=",")
                    d3D = toolbox.loadMatrixToDict(pf3D, sep=",")

                    dcoord = {}
                    for chemID in d1D2D.keys():
                        dcoord[chemID] = [d1D2D[chemID]["DIM1"], d1D2D[chemID]["DIM2"], d3D[chemID]["DIM3"]]
                        fmapChem.write("%s\t%s\n"%(chemID, map))

                    coordCentroid = calculate.centroid(dcoord)
                    dCentroid[map] = coordCentroid

                    lmap.append(map)

            fmapChem.close()


            fmapCentroid = open(pfmapCentroid, "w")
            fmapCentroid.write("Map\tX\tY\tZ\n")

            for map in dCentroid.keys():
                lcoord = [str(c) for c in dCentroid[map]]
                fmapCentroid.write("%s\t%s\n"%(map, "\t".join(lcoord)))
            fmapCentroid.close()



    def generateTableProp(self, prDSSTOXPred, pknownSDF, pLD50, pDSSToxmapID=""):

        # list of descriptors to extract
        ldesc = ["inchikey", "SMILES", "GHS_category", "EPA_category", "consensus_LD50", "LD50_mgkg", "MolWeight", "LogOH_pred", "CATMoS_VT_pred", "CATMoS_NT_pred", "CATMoS_EPA_pred",
                 "CATMoS_GHS_pred", "CATMoS_LD50_pred", "CERAPP_Ago_pred", "CERAPP_Anta_pred", "CERAPP_Bind_pred", "Clint_pred", "CoMPARA_Ago_pred", "CoMPARA_Anta_pred",
                 "CoMPARA_Bind_pred", "FUB_pred", "LogHL_pred", "LogKM_pred", "LogKOA_pred", "LogKoc_pred", "LogBCF_pred", "LogD55_pred", "LogP_pred", "MP_pred", "pKa_a_pred",
                 "pKa_b_pred", "ReadyBiodeg_pred", "RT_pred", "LogVP_pred", "LogWS_pred", "BioDeg_LogHalfLife_pred", "BP_pred", "nbLipinskiFailures"]


        if not "prmaps" in self.__dict__:
            print ("Generate map fist")
            return


        # check with DSStox ID in case of the inchikey not find
        dDSSToxmap = toolbox.loadMatrixToDict(pDSSToxmapID, sep = ",")
        print ("Load map CID to XID")


        # to bypass the file creation
        lpmaps = listdir(self.prmaps)
        for pmaps in lpmaps:
            print (pmaps)
            if search("TableProp", pmaps):
                print ("Prop table already computed")
                return

        # intialisation
        dDSSTOX = {}

        # load molecular coord file and update prop table
        dSMILE = {}
        filinDesc2D = open(self.pfdesc["1D2D"], "r")
        filinDesc2D.readline()
        flag = 0
        while flag == 0:
            try:
                l = filinDesc2D.readline().split("\t")
                ID = l[0]
                SMILES = l[1]
                dSMILE[ID] = SMILES
            except:
                flag = 1
        filinDesc2D.close()
        print ("LOAD 2D desc")

        # put in dict out -> initialization to NA
        for IDchem in dSMILE.keys():
            dDSSTOX[IDchem] = {}
            for desc in ldesc:
                dDSSTOX[IDchem][desc] = "NA"
            dDSSTOX[IDchem]["SMILES"] = dSMILE[IDchem]


            # generate inchkey
            chemMol = Chem.MolFromSmiles(dDSSTOX[IDchem]["SMILES"])
            inchi = Chem.inchi.MolToInchi(chemMol)
            inchikey = Chem.inchi.InchiToInchiKey(inchi)
            dDSSTOX[IDchem]["inchikey"] = inchikey

        print ("Dictionnary created")

        # load prediction and update table
        lppred = listdir(prDSSTOXPred)
        for ppred in lppred:
            print (ppred, "Load file")
            if ppred[-3:] == "csv":
                dtemp = toolbox.loadMatrixToDict(prDSSTOXPred + ppred, sep=",")

                for chemIDtemp in dtemp.keys():
                    if search("DTXCID", chemIDtemp):
                        try:
                            chemSID = dDSSToxmap[chemIDtemp]["dsstox_substance_id"]
                            try:
                                for k in dtemp[chemIDtemp].keys():
                                    if k in ldesc:
                                        dDSSTOX[chemSID][k] = dtemp[chemIDtemp][k]
                            except:
                                pass
                        except:
                            pass

                    else:
                        try:
                            for k in dtemp[chemIDtemp].keys():
                                if k in ldesc:
                                    dDSSTOX[chemIDtemp][k] = dtemp[chemIDtemp][k]
                        except:
                            pass

        print ("Prediction loaded")


        #load sdf
        dsdf = loadDB.sdfDB(pknownSDF, "InChI Key_QSARr", self.prmaps)
        dsdf.parseAll()


        #load LD50 file
        dLD50 = toolbox.loadMatrixToDict(pLD50)
        print ("SDF and table LD50 loaded")


        # load MAP
        lmaps = listdir(self.prmaps)
        lmapscompleted = []
        for pmap in lmaps:
            nmap = pmap.split("_")[0]
            if not nmap in lmapscompleted:
                lmapscompleted.append(nmap)
                print (nmap)
                dcoords = toolbox.loadMatrixToDict(self.prmaps + pmap, sep=",")
                pfilout = self.prmaps + str(nmap) + "_TableProp.csv"
                filout = open(pfilout, "w")
                filout.write("ID\t%s\n"%("\t".join(ldesc)))

                for chemID in dcoords.keys():
                    try: tempinchKey = dDSSTOX[chemID]["inchikey"]
                    except:continue
                    # look sdf -> map on the sdf
                    for dchemIDsdf in dsdf.lc:
                        if dchemIDsdf["InChI Key_QSARr"] == tempinchKey:
                            for ksdf in dchemIDsdf.keys():
                                if ksdf in ldesc:
                                    dDSSTOX[chemID][ksdf] = dchemIDsdf[ksdf]

                    # look in LD50 file -> map on the LD50
                    for chemIDLD50 in dLD50.keys():
                        if dLD50[chemIDLD50]["InChI Key_QSARr"] == tempinchKey:
                            for kLD50 in dLD50[chemIDLD50].keys():
                                if kLD50 in ldesc:
                                    dDSSTOX[chemID][kLD50] = dLD50[chemIDLD50][kLD50]

                    filout.write("%s\t%s\n"%(chemID, "\t".join([dDSSTOX[chemID][k] for k in ldesc])))
                filout.close()




    def generateNeighborMatrix(self, nneighbor):

        if not "prmaps" in self.__dict__:
            print ("Generate map fist")
            return

        # to bypass the file creation
        lpmaps = listdir(self.prmaps)
        for pmaps in lpmaps:
            print (pmaps)
            if search("TableNeighbors", pmaps):
                print ("Neighbors table already computed")
                return


        lmaps = listdir(self.prmaps)
        lmapscompleted = []
        for pmap in lmaps:
            if search("map", pmap):
                nmap = pmap.split("_")[0]
                if not nmap in lmapscompleted:
                    lmapscompleted.append(nmap)
                    print (nmap)
                    pcoord1D2D = self.prmaps + str(nmap) + "_map1D2D.csv"
                    pcoord3D = self.prmaps + str(nmap) + "_map3D.csv"

                    d1D2D = toolbox.loadMatrixToDict(pcoord1D2D, sep = ",")
                    d3D = toolbox.loadMatrixToDict(pcoord3D, sep=",")


                    dcor = {}
                    for chemID in d1D2D.keys():
                        x = d1D2D[chemID]["DIM1"]
                        y = d1D2D[chemID]["DIM2"]
                        z = d3D[chemID]["DIM3"]
                        dcor[chemID] = [float(x), float(y), float(z)]

                    ddist = {}
                    for ID in dcor.keys():
                        ddist[ID] = {}
                        for ID2 in dcor.keys():
                            if ID != ID2:
                                ddist[ID][ID2] = sqrt(sum([(xi - yi) ** 2 for xi, yi in zip(dcor[ID], dcor[ID2])]))

                        lID = [i[0] for i in sorted(ddist[ID].items(), key=lambda x: x[1])][:nneighbor]
                        ddist[ID] = lID


                    # write in table
                    ptableneighbor = self.prmaps + str(nmap) + "_TableNeighbors.csv"
                    ftable = open(ptableneighbor, "w")
                    ftable.write("ID\tNeighbors\n")
                    for ID in ddist.keys():
                        ftable.write("%s\t%s\n" % (ID, " ".join(ddist[ID])))
                    ftable.close()












