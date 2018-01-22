from pydpi.drug import constitution, topology, connectivity, kappa, bcut, basak, estate, moran, moreaubroto, geary, \
    charge, molproperty, moe, fingerprint
from pydpi import pydrug
from molvs import standardize_smiles, Standardizer


from rdkit.Chem.SaltRemover import SaltRemover
from rdkit import Chem

from copy import deepcopy
from os import path, getcwd, remove, system, listdir
from shutil import copy
from re import search

import toolbox
import pathFolder
import runExternalSoft

LSALT="[Co]"
LSMILESREMOVE=["[C-]#N", "[Al+3]", "[Gd+3]", "[Pt+2]", "[Au+3]", "[Bi+3]", "[Al]", "[Si+4]", "[Fe]", "[Zn]", "[Fe+2]",
               "[Ru+8]", "[Fe+]", "[Sr++]", "[Fe+3]", "[O--]", "[OH-]", "[Mn++]", "[La+3]", "[Lu+3]", "[SH-]", "[Pt+4]",
               "[Fe++]", "[W]", "[Cu+2]", "[Cr+3]", "[Tc+7]", "[Xe]", "[Tl+]", "[Zn+2]", "[F-]", "[C]", "[He]", "N#N",
               "O=O", "Cl[Ra]Cl", "[Mn+2]", "N#[N+][O-]", "II", "[Ga+3]", "[Mo+10]", "[Zn]", "[Fe]", "[Si+4]", "[Al]"]

class Descriptors:

    def __init__(self, dcompound, logfile, prout, writecheck=1):
        self.compound = dcompound
        self.prout = prout
        loader = pydrug.PyDrug()

        # if SMILES, load using SMILES code
        if not "SMILES" in dcompound.keys():
            try:
                smile = runExternalSoft.babelConvertSDFtoSMILE(dcompound["sdf"])
                self.compound["SMILES"] = smile
            except:
                print "ERROR INPUT SDF - l33"
                self.log = "ERROR"
                try:logfile.write(
                    self.compound["DATABASE_ID"] + "\t---\tERROR-SDF ORIGINAL INPUT\n")
                except:pass

                return


        #Standardize smile code
        try: smilestandadized = standardize_smiles(self.compound["SMILES"])
        except:
            logfile.write(self.compound["DATABASE_ID"] + "\t" + str(self.compound["SMILES"]) + "\tERROR-SMILES INPUT"
                                                                                               "\n")
            self.log = "ERROR"
            return

        #Standardize using molvs (http://molvs.readthedocs.io/en/latest/api.html#molvs-fragment)
        s = Standardizer()
        mol = Chem.MolFromSmiles(smilestandadized)
        molstandardized = s.standardize(mol)
        smilestandadized = Chem.MolToSmiles(molstandardized)

        # remove salt
        # 1.default
        remover = SaltRemover()
        mol = Chem.MolFromSmiles(smilestandadized)
        molcleandefault = remover(mol)
        # 2. Personal remover
        homeremover = SaltRemover(defnData=LSALT)
        molclean = homeremover(molcleandefault)
        smilesclean = Chem.MolToSmiles(molclean)
        # 3. SMILES remove other manual salts + fragments -> for fragment take one if exactly same compound
        lelem = smilesclean.split(".")
        if len(lelem) > 1:
            # reduce double, case of several salts are included - 255
            lelem = list(set(lelem))
            for smilesdel in LSMILESREMOVE:
                if smilesdel in lelem:
                    lelem.remove(smilesdel)
            try:lelem.remove("")# case of bad smile
            except:pass
            if len(lelem) == 1:
                smilesclean = str(lelem[0])
            else:
                # 4. Fragments
                #Case of fragment -> stock in log file, check after to control
                logfile.write(self.compound["DATABASE_ID"] + "\t" + str(self.compound["SMILES"]) + "\tFRAGMENT IN INPUT"
                                                                                                   "\n")
                print ".".join(lelem), " - FRAGMENTS - l66"
                self.log = "ERROR"
                return
        else:
            pass


        print self.compound["SMILES"], "SMILES IN - l25 liganddescriptors"
        print smilesclean, "SMILES without salt and standardized"

        # case where only salt are included
        if smilesclean == "":
            logfile.write(self.compound["DATABASE_ID"] + "\t" + str(self.compound["SMILES"]) + "\tEMPTY SMILES AFTER "
                                                                                               "STANDARDIZATION\n")
            print "EMPTY SMILES AFTER STANDARDIZATION - l84"
            self.log = "ERROR"
            return

        self.compound["SMILES"] = smilesclean
        self.log = "OK"

        prCpdSmi = pathFolder.createFolder(self.prout + "SMI/")
        prCpdSDF = pathFolder.createFolder(self.prout + "SDF/")

        if writecheck == 1:
            # SMILES code
            pfileSMILES = prCpdSmi + str(dcompound["DATABASE_ID"]) + ".smi"
            fileSMILES = open(pfileSMILES, "w")
            fileSMILES.write(self.compound["SMILES"])
            fileSMILES.close()

            # SDF input
            pfileSDF = prCpdSDF + str(dcompound["DATABASE_ID"]) + ".sdf"
            fileSDF = open(pfileSDF, "w")
            fileSDF.write(self.compound["sdf"])
            fileSDF.close()

        # read mol
        self.mol = loader.ReadMolFromSmile(self.compound["SMILES"])

    def get_descriptor1D2D(self):


        # 0D and 1D
        try:
            self.consti = constitution.GetConstitutional(self.mol)
        except:
            self.consti = {}
        self.compo = {}
        try:
            self.compo["nheavy"] = self.mol.GetNumHeavyAtoms()
        except:
            self.compo = {}
        try:
            self.molprop = molproperty.GetMolecularProperty(self.mol)
        except:
            self.molprop = {}


        # 2D
        try:
            self.topo = topology.GetTopology(self.mol)
        except:
            self.topo = {}
        try:
            self.connect = connectivity.GetConnectivity(self.mol)
        except:
            self.connect = {}
        try:
            self.kap = kappa.GetKappa(self.mol)
        except:
            self.kap = {}
        try:
            self.burden = bcut.GetBurden(self.mol)
        except:
            self.burden = {}
        try:
            self.basakD = basak.Getbasak(self.mol)
        except:
            self.basakD = {}
        try:
            self.est = estate.GetEstate(self.mol)
        except:
            self.est = {}
        try:
            self.moreauBurto = moreaubroto.GetMoreauBrotoAuto(self.mol)
        except:
            self.moreauBurto = {}
        try:
            self.autcormoran = moran.GetMoranAuto(self.mol)
        except:
            self.autcormoran = {}
        try:
            self.gearycor = geary.GetGearyAuto(self.mol)
        except:
            self.gearycor = {}
        try:
            self.charges = charge.GetCharge(self.mol)
        except:
            self.charges = {}
        try:
            self.MOE = moe.GetMOE(self.mol)
        except:
            self.MOE = {}

        # combine all 1D2D
        self.all1D2D = dict()
        self.all1D2D.update(deepcopy(self.consti))
        self.all1D2D.update(deepcopy(self.compo))
        self.all1D2D.update(deepcopy(self.molprop))
        self.all1D2D.update(deepcopy(self.topo))
        self.all1D2D.update(deepcopy(self.connect))
        self.all1D2D.update(deepcopy(self.kap))
        self.all1D2D.update(deepcopy(self.burden))
        self.all1D2D.update(deepcopy(self.basakD))
        self.all1D2D.update(deepcopy(self.est))
        self.all1D2D.update(deepcopy(self.moreauBurto))
        self.all1D2D.update(deepcopy(self.autcormoran))
        self.all1D2D.update(deepcopy(self.gearycor))
        self.all1D2D.update(deepcopy(self.charges))
        self.all1D2D.update(deepcopy(self.MOE))



        # listdesc
        self.l1D2D = constitution._constitutional.keys()
        self.l1D2D = self.l1D2D + ["nheavy"]
        self.l1D2D = self.l1D2D + molproperty.MolecularProperty.keys()
        self.l1D2D = self.l1D2D + topology._Topology.keys()
        self.l1D2D = self.l1D2D + connectivity._connectivity.keys()
        self.l1D2D = self.l1D2D + kappa._kapa.keys()
        self.l1D2D = self.l1D2D + bcut._bcut
        self.l1D2D = self.l1D2D + basak._basak.keys()
        self.l1D2D = self.l1D2D + estate._estate.keys()
        self.l1D2D = self.l1D2D + moreaubroto._moreaubroto.keys()
        self.l1D2D = self.l1D2D + moran._moran.keys()
        self.l1D2D = self.l1D2D + geary._geary.keys()
        self.l1D2D = self.l1D2D + charge._Charge.keys()
        self.l1D2D = self.l1D2D + moe._moe.keys()



    def get_fingerprints(self):
        # fingerprint
        self.fingerAtomPairs = fingerprint.CalculateAtomPairsFingerprint(self.mol)
        self.fingerDaylight = fingerprint.CalculateDaylightFingerprint(self.mol)
        self.fingerEstate = fingerprint.CalculateEstateFingerprint(self.mol)
        self.fingerFP4 = fingerprint.CalculateFP4Fingerprint(self.mol)
        self.fingerMACCS = fingerprint.CalculateMACCSFingerprint(self.mol)
        self.fingerMorgan = fingerprint.CalculateMorganFingerprint(self.mol)
        self.fingerTorsion = fingerprint.CalculateTopologicalTorsionFingerprint(self.mol)



    def generate3DFromSMILES(self, log):
        """
        Compute descriptors 3D from SMILES code and generate the 3D using ligprep
        :return: dictionary of descriptors in all3D
        """

        # define folder with selected 3D structure
        pr3DSDF = pathFolder.createFolder(self.prout + "SDF3D")

        # clean temp folder - used to compute 3D descriptors
        prtemp = pathFolder.createFolder(self.prout + "temp3D/", clean = 1)
        psdf3Dout = pr3DSDF + self.compound["DATABASE_ID"] + ".sdf"

        # temp SMILES
        pfilesmile = prtemp + "tem.smi"
        filesmile = open(pfilesmile, "w")
        filesmile.write(self.compound["SMILES"])
        filesmile.close()

        pdesc = ""
        # run ligprep
        if not path.exists(psdf3Dout):
            psdf3D = runExternalSoft.runLigprep(psmilin=pfilesmile)

            # case error in ligprep
            if not path.exists(psdf3D) or path.getsize(psdf3D) == 0:
                log.write(self.compound["DATABASE_ID"] + "\t" + self.compound["SMILES"] + "\t" + psdf3D)
            else:
                psdf3Dout = toolbox.selectMinimalEnergyLigPrep(psdfin=psdf3D,
                                                               psdfout=psdf3Dout)
                # take only best energy
                pathFolder.cleanFolder(prtemp)








    def writeTablesDesc(self, prresult):

        # case 1D2D
        if "all1D2D" in self.__dict__:
            if not path.exists(prresult + "1D2D.csv"):
                self.fil1D2D = open(prresult + "1D2D.csv", "w")
                # header
                self.fil1D2D.write("ID\tSMILES\t")
                self.fil1D2D.write("\t".join(self.l1D2D) + "\n")
            else:
                self.fil1D2D = open(prresult + "1D2D.csv", "a")
            self.fil1D2D.write(self.compound['DRUGBANK_ID'] + "\t" + self.compound['SMILES'])

            for desc1D2D in self.l1D2D:
                try:
                    self.fil1D2D.write("\t" + str(self.all1D2D[desc1D2D]))
                except:
                    self.fil1D2D.write("\tNA")
            self.fil1D2D.write("\n")
            self.fil1D2D.close()



def get_descriptor3D(pr3DSDF, pdesc3D):


    pPadel = runExternalSoft.runPadel(pr3DSDF)
    dpadel = toolbox.parsePadelOut(pPadel)

    filout = open(pdesc3D, "w")

    # write descriptor
    ldesc = toolbox.parsePadelOut().keys()
    filout.write("\t".join(ldesc) + "\n")
    for compound in dpadel.keys():
        filout.write(compound)
        for desc in ldesc:
            filout.write("\t" + str(dpadel[desc]))
        filout.write("\n")
    filout.close()

    return pdesc3D