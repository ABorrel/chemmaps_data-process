import pathFolder
import toolbox
import chemical
import runExternalSoft
import loadDB
import calculate

from os import path, listdir, remove
from copy import deepcopy
from math import sqrt
from rdkit import Chem
from re import search
from random import shuffle


class DSSTOX:

    def __init__(self, plistChem, istart, iend, prDSSTox, prout):

        self.plistChem = plistChem
        self.istart = istart
        self.iend = iend
        self.plog = prout + "log.txt"
        self.prDSSTox = prDSSTox
        self.prSMI = pathFolder.createFolder(self.prDSSTox + "SMI/")
        self.prout = prout



    def loadlistChem(self):

        dchem = toolbox.loadMatrixToDict(self.plistChem, sep=",")
        self.dchem = dchem



    def computeDesc2D3D (self, compute=1):

        pr2D = pathFolder.createFolder(self.prDSSTox + "2D/")
        pr3D = pathFolder.createFolder(self.prDSSTox + "3D/")
        pr3Dsdf = pathFolder.createFolder(self.prDSSTox + "3Dtemp/")

        self.pr2D = pr2D
        self.pr3D = pr3D

        if not "dchem" in self.__dict__:
            self.loadlistChem()

        lchemID = self.dchem.keys()

        i = self.istart
        imax = len(lchemID)
        if self.iend == 0 or self.iend > imax:
            iend = imax
        else:
            iend = self.iend


        # just work with open file
        if compute == 0:
            dout = {}
            while i < iend:
                if i % 1000 == 0:
                    print i
                if "dsstox_substance_id" in self.dchem[lchemID[i]].keys():
                    DSSTOXid = self.dchem[lchemID[i]]["dsstox_substance_id"]
                    pSMI = self.prSMI + DSSTOXid + ".smi"
                    if path.exists(pSMI):

                        fSMI = open(pSMI, "r")
                        llSMI = fSMI.readlines()
                        fSMI.close()
                        if len (llSMI) > 0:
                            lSMI = llSMI[0].strip().split("\t")
                            if len(lSMI) == 3:
                                SMI = lSMI[0]
                                inchikey = lSMI[1]

                                pdesc2D = pr2D + inchikey + ".txt"
                                pdesc3D = pr3D + inchikey + ".txt"

                                if path.exists(pdesc2D) and path.exists(pdesc3D):
                                    if not inchikey in dout.keys():
                                        dout[inchikey] = {}
                                        dout[inchikey]["SMI"] = [SMI]
                                        dout[inchikey]["DSSTOXid"] = [DSSTOXid]
                                    else:
                                        dout[inchikey]["SMI"].append(SMI)
                                        dout[inchikey]["DSSTOXid"].append(DSSTOXid)
                            else:
                                print pSMI
                i = i + 1
            self.ddesc = dout
            return




        dout = {}
        while i < iend:
            if i%1000 == 0:
                print i
            if "dsstox_substance_id" in self.dchem[lchemID[i]].keys():
                DSSTOXid = self.dchem[lchemID[i]]["dsstox_substance_id"]
            else:
                DSSTOXid = self.dchem[lchemID[i]]["DTXSID"]

            if "Original_SMILES" in self.dchem[lchemID[i]].keys():
                smiles = self.dchem[lchemID[i]]["Original_SMILES"].split(" ")[0]
            elif "ORIGINAL  SMILES" in self.dchem[lchemID[i]].keys():
                smiles = self.dchem[lchemID[i]]["ORIGINAL  SMILES"]
            else:
                smiles = self.dchem[lchemID[i]]["SMILES"]

            chem = chemical.chemical(DSSTOXid, smiles)
            chem.prepareChem(self.prSMI)


            # skip error prep
            if chem.err == 1:
                flog = open(self.plog, "a")
                flog.write("\n" + str(DSSTOXid) + "\n")
                flog.write(chem.log)
                flog.close()
                i = i + 1
                continue

            if compute == 1:
                chem.generate3DFromSMILES(pr3Dsdf, "RDKit")
                chem.compute1D2DDesc(pr2D)
                chem.compute3DDesc(pr3D)
                # control if all process is correct

                if chem.err == 0:
                    # write table of desc
                    chem.writeTablesDesc(pr2D, "1D2D")
                    chem.writeTablesDesc(pr3D, "3D")
                else:
                    i = i + 1
                    continue

            # case of everything works
            InchIKey = chem.inchikey
            if path.exists(pr2D + InchIKey + ".txt") and path.exists(pr3D + InchIKey + ".txt"):
                if not InchIKey in dout.keys():
                    dout[InchIKey] = {}
                    dout[InchIKey]["SMI"]=[]
                    dout[InchIKey]["DSSTOXid"] = []
                dout[InchIKey]["SMI"].append(chem.smiclean)
                dout[InchIKey]["DSSTOXid"].append(DSSTOXid)

            i = i + 1

        self.ddesc = dout


    def computepng(self):

        # draw from desc descriptors

        prSMI = self.prDSSTox + "SMI/"
        if not path.exists(prSMI):
            print "ERROR: CREATE CLEAN SMI FIRST"
            return

        prPNG = self.prDSSTox + "PNG/"
        pathFolder.createFolder(prPNG)

        if "ddesc" in self.__dict__:
            linchikey = self.ddesc.keys()
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



    def writeDescMatrix(self, typeDesc):

        # write 2D
        if typeDesc == "1D2D" or typeDesc == "1D" or typeDesc == "2D" :
            prin = self.pr2D
            pfilout = self.prout + "1D2D.csv"
        elif typeDesc == "3D":
            prin = self.pr3D
            pfilout = self.prout + "3D.csv"


        if path.exists(pfilout) and path.getsize(pfilout) > 1000:
            if not "pfdesc" in self.__dict__:
                self.pfdesc = {}
            self.pfdesc[typeDesc] = pfilout
            return pfilout
        else:
            filout = open(pfilout, "w")


        ldesc = []
        for ink in self.ddesc.keys():
            pdesc = prin + ink + ".txt"
            if path.exists(pdesc):
                dink = toolbox.loadMatrixToDict(pdesc, sep="\t")
                if ldesc == []:
                    ldesc = dink[ink].keys()
                    del ldesc[ldesc.index("ID")]
                    filout.write("ID\tSMILES\t" + "\t".join(ldesc) + "\n")

                # put back DSSTOX
                dink[ink]["ID"] = self.ddesc[ink]["DSSTOXid"][0]
                dink[ink]["SMILES"] = self.ddesc[ink]["SMI"][0]
                filout.write("%s\t%s"%(dink[ink]["ID"],  dink[ink]["SMILES"]))
                for desc in ldesc:
                    try: filout.write("\t%s"%(str(dink[ink][desc])))
                    except: filout.write("\tNA")
                filout.write("\n")
        filout.close()

        if not "pfdesc" in self.__dict__:
            self.pfdesc = {}
        self.pfdesc[typeDesc] = pfilout
        return pfilout




    def projection(self, corval, maxquantile):

        prproject = pathFolder.createFolder(self.prout + "projection" + str(corval) + "-" + str(maxquantile) + "/")
        print self.pfdesc["1D2D"], self.pfdesc["3D"], prproject, corval, maxquantile
        runExternalSoft.RComputeCor(self.pfdesc["1D2D"], self.pfdesc["3D"], prproject, corval, maxquantile)


    def generateFileMap(self, corval, maxquantile):

        prmap = pathFolder.createFolder(self.prout + "map_" + str(corval) + "-" + str(maxquantile) + "/")
        if len(listdir(prmap)) > 0:
            self.prmap = prmap
        else:
            runExternalSoft.RComputeMapFiles(self.pfdesc["1D2D"], self.pfdesc["3D"], prmap, corval, maxquantile)

            self.prmap = prmap


    def splitMap(self, nbsplit, dim):

        if not "prmap" in self.__dict__:
            print "Generate the map files first"
            return

        else:
            prout = pathFolder.createFolder(self.prmap + "split_" + str(nbsplit) + "_" + str(dim) + "/")
            self.prmaps = prout

            if len(listdir(self.prmaps)) >= (nbsplit*2):
                print "Maps already computed"
                return


            coord1D2D = self.prmap + "coord1D2D.csv"
            coord3D = self.prmap + "coord3D.csv"

            d1D2D = toolbox.loadMatrixToDict(coord1D2D, sep = ",")
            d3D = toolbox.loadMatrixToDict(coord3D, sep=",")

            # max and min 1D2D
            maxDim = 0.0
            minDim = 0.0

            nbchem = len(d1D2D.keys())
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
                lchem = d1D2D.keys()
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

            print len(dmap[d]), d


    def generateMapSplitFile(self, dim):

        if not "prmaps" in self.__dict__:
            print "Generate Maps first"
            return

        else:
            lfmap = listdir(self.prmaps)
            dCentroid = {}
            lmap = []

            pfileMapChem = self.prmap + str(dim) + "_mapChem.csv"
            pfmapCentroid = self.prmap + str(dim) + "MapCentroid.csv"
            if path.exists(pfileMapChem) and path.exists(pfmapCentroid):
                print "Already computed"
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
            print "Generate map fist"
            return


        # check with DSStox ID in case of the inchikey not find
        dDSSToxmap = toolbox.loadMatrixToDict(pDSSToxmapID, sep = ",")
        print "Load map CID to XID"


        # to bypass the file creation
        lpmaps = listdir(self.prmaps)
        for pmaps in lpmaps:
            print pmaps
            if search("TableProp", pmaps):
                print "Prop table already computed"
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
        print "LOAD 2D desc"

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

        print "Dictionnary created"

        # load prediction and update table
        lppred = listdir(prDSSTOXPred)
        for ppred in lppred:
            print ppred, "Load file"
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

        print "Prediction loaded"


        #load sdf
        dsdf = loadDB.sdfDB(pknownSDF, "InChI Key_QSARr", self.prmaps)
        dsdf.parseAll()


        #load LD50 file
        dLD50 = toolbox.loadMatrixToDict(pLD50)
        print "SDF and table LD50 loaded"


        # load MAP
        lmaps = listdir(self.prmaps)
        lmapscompleted = []
        for pmap in lmaps:
            nmap = pmap.split("_")[0]
            if not nmap in lmapscompleted:
                lmapscompleted.append(nmap)
                print nmap
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
            print "Generate map fist"
            return

        # to bypass the file creation
        lpmaps = listdir(self.prmaps)
        for pmaps in lpmaps:
            print pmaps
            if search("TableNeighbors", pmaps):
                print "Neighbors table already computed"
                return


        lmaps = listdir(self.prmaps)
        lmapscompleted = []
        for pmap in lmaps:
            if search("map", pmap):
                nmap = pmap.split("_")[0]
                if not nmap in lmapscompleted:
                    lmapscompleted.append(nmap)
                    print nmap
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












