import pathFolder
import toolbox
import runExternalSoft
import DBrequest
import calculate
import parseSDF

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


LPROP = ["inchikey", "SMILES", "preferred_name", "GHS_category", "EPA_category", "consensus_LD50", "LD50_mgkg", "MolWeight", "LogOH_pred", "CATMoS_VT_pred", "CATMoS_NT_pred", "CATMoS_EPA_pred",
                 "CATMoS_GHS_pred", "CATMoS_LD50_pred", "CERAPP_Ago_pred", "CERAPP_Anta_pred", "CERAPP_Bind_pred", "Clint_pred", "CoMPARA_Ago_pred", "CoMPARA_Anta_pred",
                 "CoMPARA_Bind_pred", "FUB_pred", "LogHL_pred", "LogKM_pred", "LogKOA_pred", "LogKoc_pred", "LogBCF_pred", "LogD55_pred", "LogP_pred", "MP_pred", "pKa_a_pred",
                 "pKa_b_pred", "ReadyBiodeg_pred", "RT_pred", "LogVP_pred", "LogWS_pred", "BioDeg_LogHalfLife_pred", "BP_pred", "nbLipinskiFailures"]



class DSSTOX:

    def __init__(self, plistChem, nameMap,  istart, iend, prDesc, prout):

        self.plistChem = plistChem
        self.nameMap = nameMap
        self.istart = istart
        self.iend = iend
        self.plog = prout + "log.txt"
        self.prDesc = prDesc
        self.prout = prout



    def loadlistChem(self):

        prForDB = pathFolder.createFolder(self.prout + "forDB/")
        pfilout = prForDB + "db.csv"
        #try:remove(pfilout)
        #except:pass
        #print(pfilout)
        if path.exists(pfilout):
            dchem = toolbox.loadMatrixToDict(pfilout, sep="\t")
        else:
            dchem = {}
            if self.nameMap == "dsstox":
                dchemIn = toolbox.loadMatrixToDict(self.plistChem, sep=",") #rewrite pfas and tox21 with comma
            else:
                dchemIn = toolbox.loadMatrixToDict(self.plistChem, sep="\t")

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
                dchem[DTXSID]["db_id"] = DTXSID
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
            filout.write("db_id\tsmiles_origin\tsmiles_clean\tinchikey\tqsar_ready\t%s\n"%(self.nameMap))
            for chem in dchem.keys():
                filout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
                chem, dchem[chem]["smiles_origin"], dchem[chem]["smiles_clean"], dchem[chem]["inchikey"],
                dchem[chem]["qsar_ready"], 1))
            filout.close()
        self.dchem = dchem




    def pushChemInDB(self):

        if not "dchem" in self.__dict__:
            self.loadlistChem()

        if self.nameMap == "dsstox":
            cDB = DBrequest.DBrequest()
            cDB.verbose = 0
            for chem in self.dchem.keys():
                cDB.addElement("dsstox_chem", ["db_id", "smiles_origin", "smiles_clean", "inchikey", "qsar_ready"],
                                    [chem, self.dchem[chem]["smiles_origin"], self.dchem[chem]["smiles_clean"],
                                    self.dchem[chem]["inchikey"], self.dchem[chem]["qsar_ready"]])
        else:
            cDB = DBrequest.DBrequest()
            ddd
            # have to be wrtie


    def computeDesc (self, insertDB =0, w=0):

        if not "dchem" in self.__dict__:
            self.loadlistChem()

        pfilout1D2D = self.prout + "1D2D.csv"
        pfilout3D = self.prout + "3D.csv"

        if path.exists(pfilout1D2D) and path.exists(pfilout3D) and w == 1 and insertDB == 0:
            self.p1D2D = pfilout1D2D
            self.p3D = pfilout3D
            return


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
        ldesc1D2D = Chemical.getLdesc("1D2D")
        ldesc3D = Chemical.getLdesc("3D")

        if insertDB == 1:
            cDB = DBrequest.DBrequest()

        if w == 1:
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
                        valDesc = [cChem.all2D[desc1D2D] for desc1D2D in ldesc1D2D]
                        valDesc = ['-9999' if desc == "NA" else desc for desc in valDesc]

                        w1D2D = "{" + ",".join(["\"%s\"" % (desc) for desc in valDesc]) + "}"
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
                            valDesc = [cChem.all3D[desc3D] for desc3D in ldesc3D]
                            valDesc = ['-9999' if desc == "NA" else desc for desc in valDesc]

                            w3D = "{" + ",".join(["\"%s\"" % (desc) for desc in valDesc]) + "}"
                            cDB.addElement("desc_3d", ["inchikey", "desc_value"], [cChem.inchikey, w3D])
            i = i + 1

        if w == 1:
            filout1D2D.close()
            filout3D.close()

        self.p1D2D = pfilout1D2D
        self.p3D = pfilout3D



    def pushPropInDB(self):
        
        tableDB = "dsstox_chem"

        if not "dchem" in self.__dict__:
            self.loadlistChem()
        cDB = DBrequest.DBrequest()
        cDB.verbose = 1
        for chem in self.dchem.keys():
            cmdSQL = "SELECT count(*) FROM %s WHERE inchikey = '%s';"%(tableDB, self.dchem[chem]["inchikey"])
            if cDB.execCMD(cmdSQL) > 0:
                cmdSQL = "UPDATE %s SET %s = 1 WHERE inchikey = '%s'"%(tableDB, self.nameMap, self.dchem[chem]["inchikey"])
                cDB.execCMD(cmdSQL)
            dddd
        return 





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
                #out1D2D = cDB.getRow("%s_coords"%(self.nameMap), "inchikey='%s'" % (chem))
                #if out1D2D == []:
                    nbdim1d2d = len(dcoord1D2D[chem].keys()) - 1
                    nbdim3d = len(dcoord3D[chem].keys()) - 1

                    w1D2D = "{" + ",".join(["\"%s\"" % (dcoord1D2D[chem]["DIM" + str(i)]) for i in range(1, nbdim1d2d + 1)]) + "}"
                    w3D = "{" + ",".join(["\"%s\"" % (dcoord3D[chem]["DIM3-" + str(i)]) for i in range(1, nbdim3d + 1)]) + "}"
                    print(w1D2D)
                    ddd
                    cDB.addElement("%s_coords"%(self.nameMap), ["inchikey", "dim1d2d", "dim3d", "in_db"], [chem, w1D2D, w3D, "1"])


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
            prout = pathFolder.createFolder(self.prmap + "split_" + str(nbsplit) + "/")
            self.prmaps = prout

            if not "psplitMap" in self.__dict__:
                self.psplitMap = {}

            # generate only one file with chem and map
            if dim == 1:
                pfilout = prout + "mapx_split.csv"
                
            elif dim == 2: 
                pfilout = prout + "mapy_split.csv"

            else:
                pfilout = prout + "mapz_split.csv"

            self.psplitMap[dim] = pfilout
            if path.exists(pfilout):
                return

            coord1D2D = self.prmap + "coord1D2D.csv"
            coord3D = self.prmap + "coord3D.csv"

            if dim == 1 or dim == 2:
                din = toolbox.loadMatrixCoords(coord1D2D)
            else:
                din = toolbox.loadMatrixCoords(coord3D)

            # max and min 1D2D
            maxDim = 0.0
            minDim = 0.0

            nbchem = len(list(din.keys()))
            nbchembymap = int(nbchem/nbsplit)

            # calibrate max and min
            print("== Initiate calibration ==")
            for chem in din.keys():

                if dim == 1 or dim == 3:
                    dimVal = din[chem][0]
                elif dim == 2:
                    dimVal = din[chem][1]

                if dimVal > maxDim:
                    maxDim = dimVal
                if dimVal < minDim:
                    minDim = dimVal
            print("== End calibration ==")

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
                lchem = list(din.keys())
                nbchem = len(lchem)
                while ichem < nbchem:

                    if dim == 1 or dim == 3:
                        valtemp = din[lchem[ichem]][0]
                    elif dim == 2:
                        valtemp = din[lchem[ichem]][1]


                    if valtemp < dimVal:

                        dmap[imap].append(deepcopy(lchem[ichem]))
                        del din[lchem[ichem]]
                        del lchem[ichem]
                        nbchem = nbchem - 1
                        continue
                    else:
                        ichem = ichem + 1

        print("==== Write output ====")
        filout = open(pfilout, "w")
        filout.write("inchikey\tmap\n")        
        for d in dmap.keys():
            for chem in dmap[d]:
                filout.write("%s\t%s\n"%(chem, d))
        filout.close()


    # have to be optimize
    def generateCentroidFile(self):

        if not "prmaps" in self.__dict__:
            print ("Generate Maps first")
            return

        if not "psplitMap" in self.__dict__:
            print("Generate the split map first")
            return

        else:
            lpfmap = list(self.psplitMap.values())
            #print(lpfmap)

            pfilout = self.prmaps + "centroids.csv"
            #if path.exists(pfilout):
            #    return 
            
            coords1D2D = toolbox.loadMatrixCoords(self.pcoords1D2D)
            coords3D = toolbox.loadMatrixCoords(self.pcoords3D)
            
            dout = {}
            for pmap in lpfmap:
                
                print(pmap)
                nameMap = pmap.split("/")[-1].split("_")[0]
                dmap = toolbox.loadMatrixToDict(pmap)
                print(nameMap)

                i = 1
                while 1:
                    lcoords = []
                    for chem in dmap.keys():
                        if int(dmap[chem]["map"]) == i :
                            lcoords.append([coords1D2D[chem][0], coords1D2D[chem][1], coords3D[chem][0]])
                    if lcoords == []:
                        break
                    else:
                        print(len(lcoords))
                        print(lcoords[0])
                        coordCentroid = calculate.centroid(lcoords)
                        dout[nameMap + "_" + str(i)] = coordCentroid
                    i = i + 1


            #print(dout)
            filout = open(pfilout, "w")
            filout.write("map\tx\ty\tz\n")
            for map in dout.keys():
                filout.write("%s\t%s\t%s\t%s\n"%(map, dout[map][0], dout[map][1], dout[map][2]))
            filout.close()



    def pushDssToxNamePropInDB(self):

        cDB = DBrequest.DBrequest()
        cDB.verbose = 1
        i = 1
        for PROP in LPROP:
            cDB.addElement("dsstox_prop_name", ["id", "name"], [i, PROP])
            i = i + 1

    def generateTablePropAllDSSTOX(self, prDSSTOXPred, pknownSDF, pLD50, pDSSTOXMapOnCID, insertDB=0):

        
        pTableinfo = self.prout + "tablePropForDB.csv"
        if path.exists(pTableinfo) and insertDB == 0:
            self.pTableInAll = pTableinfo
            return

        #print ("LOAD INFO FROM DCHEM")
        #if not "dchem" in self.__dict__:
        #    self.loadlistChem()

        # intialisation ful dictionnary
        dDSSTOX = {}
        dmapCIDtoSID = {}
        print ("LOAD INFO MAP SID to CID")

        filMap = open(pDSSTOXMapOnCID, "r", encoding="utf8", errors="ignore")
        llines = filMap.readlines()
        filMap.close()
        
        lhead = llines[0].replace("\"", "")
        lhead = lhead.strip().split(",")
        #print(lhead)
        iDSSSID = lhead.index("dsstox_substance_id")
        iDSSCID = lhead.index("DSSTox_Structure_Id")
        iname = lhead.index("preferred_name")
        i = 1
        imax = len(llines)
        while i < imax:#####################################
            lineClean = toolbox.formatLine(llines[i])
            lelem = lineClean.strip().split(",")
            try:
                dDSSTOX[lelem[iDSSSID]] = {}
                dDSSTOX[lelem[iDSSSID]]["preferred_name"] = lelem[iname]
                dDSSTOX[lelem[iDSSSID]]["SMILES"] = self.dchem[lelem[iDSSSID]]["smiles_clean"]
                dDSSTOX[lelem[iDSSSID]]["inchikey"] = self.dchem[lelem[iDSSSID]]["inchikey"]
                dmapCIDtoSID[lelem[iDSSCID]] = lelem[iDSSSID]
            except:
                pass
            i = i + 1
        filMap.close()


        print("INIT DICTIONNARY")
        # put in dict out -> initialization to NA
        for chem in dDSSTOX.keys():
            for PROP in LPROP[3:]:
                try: dDSSTOX[chem][PROP] = "NA"
                except: break


        print("LOAD PRED")
        # load prediction and update table
        lppred = listdir(prDSSTOXPred)
        for ppred in lppred:##########################################
            if ppred[-3:] == "csv":
                print (ppred, "Load file")
                dtemp = toolbox.loadMatrixToDict(prDSSTOXPred + ppred, sep=",")
                k1 = list(dtemp.keys())[0]
                #print(dtemp[k1])
                #dddd
                for chemIDtemp in dtemp.keys():
                    DTXCID = dtemp[chemIDtemp]["MoleculeID"]
                    try: DTXSID = dmapCIDtoSID[DTXCID]
                    except: continue
                    for k in dtemp[chemIDtemp].keys():
                        if k in LPROP[3:]:
                            dDSSTOX[DTXSID][k] = dtemp[chemIDtemp][k]

        print ("PRED LOAD")

        print("LOAD SDF AND LD50")
        #load sdf
        dsdf = parseSDF.parseSDF(pknownSDF, "InChI Key_QSARr", self.prout)
        dsdf.parseAll()


        #load LD50 file
        dLD50 = toolbox.loadMatrixToDict(pLD50)
        print ("SDF and table LD50 loaded")

        
        for chem in dDSSTOX.keys():
            tempinchKey = dDSSTOX[chem]["inchikey"]
            # look sdf -> map on the sdf
            for dchemIDsdf in dsdf.lc:
                if dchemIDsdf["InChI Key_QSARr"] == tempinchKey:
                    for ksdf in dchemIDsdf.keys():
                        if ksdf in LPROP[3:]:
                            dDSSTOX[chem][ksdf] = dchemIDsdf[ksdf]

            # look in LD50 file -> map on the LD50
            for chemIDLD50 in dLD50.keys():
                if dLD50[chemIDLD50]["InChI Key_QSARr"] == tempinchKey:
                    for kLD50 in dLD50[chemIDLD50].keys():
                        if kLD50 in LPROP[3:]:
                            dDSSTOX[chem][kLD50] = dLD50[chemIDLD50][kLD50]



        print("WRITE TABLE")
        # load MAP

        filout = open(pTableinfo, "w")
        filout.write("ID\t%s\n"%("\t".join(LPROP)))
        for chem in dDSSTOX.keys():
            filout.write("%s\t%s\n"% (chem, "\t".join([str(dDSSTOX[chem][prop]) for prop in LPROP])))
        filout.close()

        self.pTableInAll = pTableinfo


    def pushTablePropAllInDB(self):

        if not "pTableInAll" in self.__dict__:
            print ("GENERATE TABLE FIRST")

        dtopush = toolbox.loadMatrixToDict(self.pTableInAll)
        cDB = DBrequest.DBrequest()
        cDB.verbose = 0
        i = 0
        i = 500000
        lchem = list(dtopush.keys())
        imax = len(lchem)
        while i < imax:
        #for chem in dtopush.keys():
            chem = lchem[i]
            outChem = cDB.getRow("dsstox_prop", "db_id='%s'" % (chem))
            if outChem == []:
                wprop = "{" + ",".join(["\"%s\"" % (dtopush[chem][PROP].replace("'", "")) for PROP in LPROP]) + "}"
                cDB.addElement("dsstox_prop", ["db_id", "prop_value"], [chem, wprop])
            i = i + 1
        return 





    def generateNeighborMatrix(self, nbNeighbor, lnDim):

        if not "pcoords1D2D" in self.__dict__:
            print("Compute Coord first")
            return 1
        else:
            if self.nameMap == "dsstox":
                # no N dimension because to slow
                prNeighbor = pathFolder.createFolder(self.prout + "Neighbors/")
                #pfilout = prNeighbor + "Table_DIM1D2D-2_1.csv"
                #if path.exists(pfilout):
                #    return
                

                    
                dDim1D2D = toolbox.loadMatrixCoords(self.pcoords1D2D)
                dDim3D = toolbox.loadMatrixCoords(self.pcoords3D)
                lpfmap = self.psplitMap
                lmap = []
                for imap in lpfmap.keys():
                    lmap.append(toolbox.loadMatrixToDict(lpfmap[imap]))
                #print(lmap)

                # from 1D2D coord
                lchem = list(dDim1D2D.keys())
                shuffle(lchem)
                i = 0
                imax = len(lchem)
                while i < imax:
                    inch = lchem[i]
                    
                    pfilout = prNeighbor + inch
                    if path.exists(pfilout):
                        i = i + 1
                        continue
                    filout = open(pfilout, "w")
                    filout.write("ID\tNeighbors\n")

                    # define map where we inspect
                    linmap = []
                    for dmap in lmap:
                        mapin = int(dmap[inch]["map"])
                        for chem in dmap.keys():
                            if int(dmap[chem]["map"]) == mapin or int(dmap[chem]["map"]) == (mapin + 1) or int(dmap[chem]["map"]) == (mapin - 1):
                                if not chem in lmap:
                                    linmap.append(chem)
                        
                    print(len(linmap))
                    ddist = {}
                    ddist[inch] = {}
                    for ID in linmap:
                        if ID != inch:
                            ddist[inch][ID] = sqrt(sum([(xi - yi) ** 2 for xi, yi in zip([dDim1D2D[ID][0], dDim1D2D[ID][1], dDim3D[ID][0]], [dDim1D2D[inch][0], dDim1D2D[inch][1], dDim3D[inch][0]],)]))

                    lID = [i[0] for i in sorted(ddist[inch].items(), key=lambda x: x[1])][:nbNeighbor]
                    filout.write("%s\t%s\n" % (inch, " ".join(lID)))
                    filout.close()
                    i = i + 1
            else:
                # compute all dimension openning withou restriction
                if lnDim == []:
                    dDim1D2D = toolbox.loadMatrixToDict(self.pcoords1D2D, sep = ",")
                    dDim3D = toolbox.loadMatrixToDict(self.pcoords3D, sep=",")
                    chem1 = list(dDim1D2D.keys())[0]
                    n1D2D = len(list(dDim1D2D[chem1].keys())) - 1
                    n3D = len(list(dDim3D[chem1].keys())) - 1
                    lnDim = [n1D2D, n3D]
                
                prNeighbor = pathFolder.createFolder(self.prout + "Neighbors/")
                pfilout = prNeighbor + "Table_DIM1D2D-" + str(lnDim[0]) + "_" + str(lnDim[1]) + ".csv"
                if path.exists(pfilout):
                    return
                else:
                    dDim1D2D = toolbox.loadMatrixToDict(self.pcoords1D2D, sep = ",")
                    dDim3D = toolbox.loadMatrixToDict(self.pcoords3D, sep=",")

                    dcor = {}
                    # from 1D2D coord
                    for inch in dDim1D2D.keys():
                        dcor[inch] = []

                        i = 1
                        while i <= lnDim[0]:
                            dcor[inch].append(float(dDim1D2D[inch]["DIM" + str(i)]))
                            i = i + 1

                        i = 1
                        while i <= lnDim[1]:
                            dcor[inch].append(float(dDim3D[inch]["DIM3-" + str(i)]))
                            i = i + 1

                    ddist = {}
                    for ID in dcor.keys():
                        ddist[ID] = {}
                        for ID2 in dcor.keys():
                            if ID != ID2:
                                ddist[ID][ID2] = sqrt(sum([(xi - yi) ** 2 for xi, yi in zip(dcor[ID], dcor[ID2])]))

                        lID = [i[0] for i in sorted(ddist[ID].items(), key=lambda x: x[1])][:nbNeighbor]
                        ddist[ID] = lID

                    # write in table
                    ftable = open(pfilout, "w")
                    ftable.write("ID\tNeighbors\n")
                    for ID in ddist.keys():
                        ftable.write("%s\t%s\n" % (ID, " ".join(ddist[ID])))
                    ftable.close()


    def pushNeighbors(self):
        prneighbor = pathFolder.createFolder(self.prout + "Neighbors/")
        ptable3Dim = prneighbor + "Table_DIM1D2D-2_1.csv"
        ptableNDim = prneighbor + "Table_DIM1D2D-131_254.csv"
        if path.exists(ptable3Dim) and path.exists(ptableNDim):
            ddist3D = toolbox.loadMatrixToDict(ptable3Dim)
            for chem in ddist3D.keys():
                ddist3D[chem] = ddist3D[chem]["Neighbors"].split(" ")
            ddistND = toolbox.loadMatrixToDict(ptableNDim)
            for chem in ddistND.keys():
                ddistND[chem] = ddistND[chem]["Neighbors"].split(" ")

            cDB = DBrequest.DBrequest()
            cDB.verbose = 0
            for chem in ddist3D.keys():
                # print(chem)
                out1D2D = cDB.getRow("%s_neighbors"%(self.nameMap), "inchikey='%s'" % (chem))
                if out1D2D == []:
                    w3D = "{" + ",".join(["\"%s\"" % (neighbor) for neighbor in ddist3D[chem]]) + "}"
                    wND = "{" + ",".join(["\"%s\"" % (neighbor) for neighbor in ddistND[chem]]) + "}"
                    cDB.addElement("%s_neighbors"%(self.nameMap), ["inchikey", "neighbors_dim3", "neighbors_dimn"], [chem, w3D, wND])










