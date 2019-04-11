import pathFolder
import toolbox
import chemical
import runExternalSoft
import loadDB

from os import path, listdir
from copy import deepcopy
from math import sqrt
from rdkit import Chem
from re import search


class DSSTOX:

    def __init__(self, plistChem, istart, iend, prout):

        self.plistChem = plistChem
        self.prDSSTOX = pathFolder.createFolder(prout)
        self.istart = istart
        self.iend = iend
        self.prSMI = pathFolder.createFolder(self.prDSSTOX + "SMI/")
        self.plog = prout + "log.txt"
        self.prout = prout


    def loadlistChem(self):

        dchem = toolbox.loadMatrixToDict(self.plistChem, sep=",")
        self.dchem = dchem



    def computeDesc2D3D (self, compute=1):
        pr2D = pathFolder.createFolder(self.prDSSTOX + "2D/")
        pr3D = pathFolder.createFolder(self.prDSSTOX + "3D/")
        pr3Dsdf = pathFolder.createFolder(self.prDSSTOX + "3Dtemp/")

        self.pr2D = pr2D
        self.pr3D = pr3D

        if compute == 0:
            return

        if not "dchem" in self.__dict__:
            self.loadlistChem()


        lchemID = self.dchem.keys()

        i = self.istart
        if self.iend == 0:
            iend = len(lchemID)
        else:
            iend = self.iend

        dout = {}
        while i < iend:
            InchIKey = self.dchem[lchemID[i]]["InChI Key_QSARr"]
            DSSTOXid = self.dchem[lchemID[i]]["DSSTox_Structure_Id"]
            smiles = self.dchem[lchemID[i]]["Original_SMILES"].split(" ")[0]

            chem = chemical.chemical(InchIKey, smiles)
            chem.prepareChem(self.prSMI)

            if compute == 1:
                chem.generate3DFromSMILES(pr3Dsdf, "Rdkit")
                chem.compute1D2DDesc(pr2D)
                chem.compute3DDesc(pr3D)

            # control if all process is correct
            if chem.err == 1:
                flog = open(self.plog, "a")
                flog.write("\n" + str(DSSTOXid) + "\n")
                flog.write(chem.log)
                flog.close()
            else:
                # write table of desc
                chem.writeTablesDesc(pr2D, "1D2D")
                chem.writeTablesDesc(pr3D, "3D")

                #print chem.log

                if not InchIKey in dout.keys():
                    dout[InchIKey] = {}
                    dout[InchIKey]["SMI"]=[]
                    dout[InchIKey]["DSSTOXid"] = []
                dout[InchIKey]["SMI"].append(chem.smiclean)
                dout[InchIKey]["DSSTOXid"].append(DSSTOXid)

            i = i + 1

        self.ddesc = dout


    def writeDescMatrix(self, typeDesc):

        # write 2D
        if typeDesc == "1D2D" or typeDesc == "1D" or typeDesc == "2D" :
            prin = self.pr2D
            pfilout = self.prout + "1D2D.csv"
        elif typeDesc == "3D":
            prin = self.pr3D
            pfilout = self.prout + "3D.csv"


        if path.exists(pfilout):
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
        if listdir(prmap) > 0:
            self.prmap = prmap
        else:
            runExternalSoft.RComputeMapFiles(self.pfdesc["1D2D"], self.pfdesc["3D"], prmap, corval, maxquantile)
            self.prmap = prmap


    def splitMap(self, nbsplit):

        if not "prmap" in self.__dict__:
            print "Generate the map files first"
            return

        else:
            prout = pathFolder.createFolder(self.prmap + "split_" + str(nbsplit) + "/")
            self.prmaps = prout

            if len(listdir(self.prmaps)) > 0:
                print "Maps already computed"
                return


            coord1D2D = self.prmap + "coord1D2D.csv"
            coord3D = self.prmap + "coord3D.csv"

            d1D2D = toolbox.loadMatrixToDict(coord1D2D, sep = ",")


            # max and min 1D2D
            maxX = 0.0
            maxY = 0.0
            minX = 0.0
            minY = 0.0

            nbchem = len(d1D2D.keys())
            nbbychem = int(nbchem / nbsplit)

            for chem in d1D2D.keys():
                X = float(d1D2D[chem]["DIM1"])
                Y = float(d1D2D[chem]["DIM2"])
                if X > maxX:
                    maxX = X
                if X < minX:
                    minX = X

                if Y > maxY:
                    maxY = Y
                if Y < minY:
                    minY = Y

            print minX, maxX
            print minY, maxY

            xi = (maxX - minX) / sqrt(nbsplit)
            yi = (maxY - minY) / sqrt(nbsplit)

            print xi,yi
            dmap = {}
            imap = 1
            dmap[imap] = []
            y = minY


            while y < maxY:
                print imap
                y = y + yi

                x = minX
                while x < maxX:
                    x = x + 0.10
                    if len(dmap[imap]) > 10000:
                        imap = imap + 1
                        dmap[imap] = []
                    ichem = 0
                    lchem = d1D2D.keys()
                    nbchem = len(lchem)
                    while ichem < nbchem:
                        if float(d1D2D[lchem[ichem]]["DIM1"]) < x and float(d1D2D[lchem[ichem]]["DIM2"]) < y:
                            dmap[imap].append(deepcopy(d1D2D[lchem[ichem]]))
                            del d1D2D[lchem[ichem]]
                            del lchem[ichem]
                            nbchem = nbchem - 1
                            continue
                        else:
                            ichem = ichem + 1


        d3D = toolbox.loadMatrixToDict(coord3D, sep = ",")

        for d in dmap.keys():
            pfilout1D2D = prout + str(d) + "_map1D2D.csv"
            pfilout3D = prout + str(d) + "_map3D.csv"

            filout1D2D = open(pfilout1D2D, "w")
            filout3D = open(pfilout3D, "w")
            filout1D2D.write("\"\",\"DIM1\",\"DIM2\"\n")
            filout3D.write("\"\",\"DIM1\",\"DIM2\"\n")

            for chem in dmap[d]:
                filout1D2D.write("\"%s\",%s,%s\n"%(chem[""], chem["DIM1"], chem["DIM2"]))
                filout3D.write("\"%s\",%s,%s\n" % (chem[""], d3D[chem[""]]["DIM1"], d3D[chem[""]]["DIM2"]))

            filout1D2D.close()
            filout3D.close()


            print len(dmap[d]), d



    def generateTableProp(self, prDSSTOXPred, pknownSDF, pLD50):

        # list of descriptors to extract
        ldesc = ["inchkey", "SMILES", "GHS_category", "EPA_category", "consensus_LD50", "LD50_mgkg", "MolWeight", "LogOH_pred", "CATMoS_VT_pred", "CATMoS_NT_pred", "CATMoS_EPA_pred",
                 "CATMoS_GHS_pred", "CATMoS_LD50_pred", "CERAPP_Ago_pred", "CERAPP_Anta_pred", "CERAPP_Bind_pred", "Clint_pred", "CoMPARA_Ago_pred", "CoMPARA_Anta_pred",
                 "CoMPARA_Bind_pred", "FUB_pred", "LogHL_pred", "LogKM_pred", "LogKOA_pred", "LogKoc_pred", "LogBCF_pred", "LogD55_pred", "LogP_pred", "MP_pred", "pKa_a_pred",
                 "pKa_b_pred", "ReadyBiodeg_pred", "RT_pred", "LogVP_pred", "LogWS_pred", "BioDeg_LogHalfLife_pred", "BP_pred", "nbLipinskiFailures"]


        if not "prmaps" in self.__dict__:
            print "Generate map fist"
            return


        lpmaps = listdir(self.prmaps)
        for pmaps in lpmaps:
            if search("TableProp", pmaps):
                print "Prop table already computed"
                return

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

        # put in dict out
        for IDchem in dSMILE.keys():
            dDSSTOX[IDchem] = {}
            for desc in ldesc:
                dDSSTOX[IDchem][desc] = "NA"
            dDSSTOX[IDchem]["SMILES"] = dSMILE[IDchem]


            # generate inchkey
            chemMol = Chem.MolFromSmiles(dDSSTOX[IDchem]["SMILES"])
            inchi = Chem.inchi.MolToInchi(chemMol)
            inchikey = Chem.inchi.InchiToInchiKey(inchi)
            dDSSTOX[IDchem]["inchkey"] = inchikey

        print "Dictionnary created"


        # load prediction and update table
        lppred = listdir(prDSSTOXPred)
        for ppred in lppred:
            print ppred
            if ppred[-3:] == "csv":
                dtemp = toolbox.loadMatrixToDict(prDSSTOXPred + ppred, sep=",")

                for chemID in dtemp.keys():
                    for k in dtemp[chemID].keys():
                        if k in ldesc:
                            # if DSSTOX is not already in dict just continue
                            try:
                                dDSSTOX[chemID][k] = dtemp[chemID][k]
                            except:
                                continue


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
                    try: tempinchKey = dDSSTOX[chemID]["inchkey"]
                    except:continue
                    # look sdf
                    for dchemIDsdf in dsdf.lc:
                        if dchemIDsdf["InChI Key_QSARr"] == tempinchKey:
                            for ksdf in dchemIDsdf.keys():
                                if ksdf in ldesc:
                                    dDSSTOX[chemID][ksdf] = dchemIDsdf[ksdf]

                    # look in LD50 file
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
                        z = d3D[chemID]["DIM1"]
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












"""

def Run (psdf, pranalysis, kname, lkinfo=[], corval=0.8, maxquantile=80, control=1, Desc1D2D=1, generation3D = 1, Desc3D=1, projection=1, JS=1, drawPNG=1):

    pathFolder.createFolder(pranalysis)

    ###################################
    # analysis and compute descriptor #
    ###################################

    db = loadDB.sdfDB(psdf, kname, pranalysis)
    db.writeTable("TaleProp.csv")


    # descriptor computation #
    ##########################
    prDesc = pranalysis + "Desc/"
    pathFolder.createFolder(prDesc, clean=0)
    dpfiledesc = computeDB.computeDesc(psdf, prDesc, control=control, Desc1D2D=Desc1D2D, generation3D=generation3D, Desc3D=Desc3D, namek=kname)
    # analyse projection  and compute coordinate #
    ##############################################
    prproject = pathFolder.createFolder(pranalysis + "projection" + str(corval) + "-" + str(maxquantile) + "/")
    if projection == 1:
        print dpfiledesc["1D2D"], dpfiledesc["3D"], prproject, corval, maxquantile
        runExternalSoft.RComputeCor(dpfiledesc["1D2D"], dpfiledesc["3D"], prproject, corval, maxquantile)

    ###################
    # for the website #
    ###################


    # 1. compute png #
    ##################
    if drawPNG == 1:
        prpng = pranalysis + "cpdpng/"
        pathFolder.createFolder(prpng)

        # draw from desc descriptors
        prSMI = prDesc + "SMI/"
        if path.exists(prSMI):
            lsmi = listdir(prSMI)
            # control if nSDF = nPNG
            if len(lsmi) != len(listdir(prpng)):
                for smifile in lsmi:
                    runExternalSoft.molconvert(prSMI + smifile, prpng + smifile.split(".")[0] + ".png")
        else:
            # case where considering only original map
            db.drawMolecules(prpng)


    ##################
    # JS file update #
    ##################
    if JS == 1:
        prforjsupdate = pranalysis + "JS" + str(corval) + "-" + str(maxquantile) + "/"
        pathFolder.createFolder(prforjsupdate)
        pfileDataJS = prforjsupdate + "data.js"

        if path.exists(pfileDataJS):
            remove(pfileDataJS)

        # 2. update JS coords #
        #######################
        pcoordsCombine = prproject + "coordPCAcombined.csv"
        if path.exists(pcoordsCombine):
            createJS.formatCoordinates(pcoordsCombine, pfileDataJS)


        # 3. update JS properties #
        ###########################
        createJS.formatInfo(db, dpfiledesc["1D2D"], lkinfo, pfileDataJS)


        # 4. update neighborhood #
        ##########################
        createJS.extractCloseCompounds(pcoordsCombine, 20, pfileDataJS)


"""