from re import search
from os import listdir

import pathFolder
import runExternalSoft


class parseSDF:

    def __init__(self, psdf, kname, prout):
        self.psdf = psdf
        self.prout = prout
        self.name = kname


    def parseAll (self):

        l_chem = []
        l_prop = []

        filin = open(self.psdf, "r", encoding="utf8")
        handle_read = filin.read()
        filin.close()
        #print handle_read

        l_compound = handle_read.split("$$$$\n")

        for compound in l_compound:
            dcompound = {}
            dcompound["sdf"] = compound + "$$$$"
            llines = compound.split("\n")
            i = 0
            nblines = len(llines)
            while i < nblines:
                if search("> <", llines[i]):
                    kin = llines[i].split("> <")[1]
                    kin = kin.split(">")[0]
                    valuek = llines[i+1].strip()
                    dcompound[kin] = valuek
                    l_prop.append(kin)
                i += 1
            if self.name in list(dcompound.keys()):
                l_chem.append(dcompound)

        #output list of chemicals
        self.lc=l_chem
        l_prop = list(set(l_prop))# remove dup
        self.l_prop = l_prop

    def writeTable(self, filname):

        if not "lc" in dir(self):
            self.parseAll()

        ptable = self.prout + str(filname)
        filout = open(ptable, "w")

        lheader = [self.name]
        for compound in self.lc:
            for k in compound.keys():
                if not k in lheader:
                    lheader.append(k)

        lheader.remove("sdf")
        filout.write("\t".join(lheader) + "\n")

        for compound in self.lc:
            try: filout.write(compound[self.name])
            except:
                continue
            for h in lheader[1:]:
                if h in list(compound.keys()):
                    filout.write("\t" + str(compound[h]))
                else:
                    filout.write("\tNA")
            filout.write("\n")
        filout.close()

    def renameHeader(self):
        if not "lc" in dir(self):
            self.parseAll()
        for compound in self.lc:
            sdfin = compound["sdf"]
            sdfin = sdfin.split("\n")
            sdfin[0] = compound[self.name]
            compound["sdf"] = "\n".join(sdfin)

    def splitSDF(self):

        if not "lc" in dir(self):
            self.parseAll()

        prSDF = self.prout + "cpdsdf/"
        pathFolder.createFolder(prSDF)
        self.prsdf = prSDF

        if len(self.lc) == len(listdir(prSDF)):
            return
        else:
            for compound in self.lc:
                pfilout = self.prsdf + compound[self.name] + ".sdf"
                filout = open(pfilout, "w")
                filout.write(compound["sdf"])
                filout.close()

    def drawMolecules(self, prpng):

        if not "prsdf" in dir(self):
            self.splitSDF()

        pathFolder.createFolder(prpng)
        self.prpng = prpng

        if len(listdir(self.prsdf)) == len(listdir(self.prpng)):
            return

        lfsdf = listdir(self.prsdf)
        for sdf in lfsdf:
            runExternalSoft.molconvert(self.prsdf + sdf, self.prpng + sdf[:-3] + "png")

    def write_table_prop(self, nametable):

        if not "lc" in dir(self):
            self.parseAll()

        pfilout = self.prout + nametable + ".csv"
        filout = open(pfilout, "w", encoding="utf-8")

        # header table
        filout.write("ID" + "\t" + "\t".join(self.l_prop) + "\n")

        for compound in self.lc:
            try:
                filout.write(compound[self.name] + "\t")
            except: continue
            lwrite = []
            for kin in self.l_prop:
                if kin in list(compound.keys()):
                    if kin == "SYNONYMS":
                        lwrite.append(str(compound[kin]).lower())
                    else:
                        lwrite.append(str(compound[kin]))
                else:
                    lwrite.append("NA")
            filout.write("\t".join(lwrite) + "\n")
        filout.close()







