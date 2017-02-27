from pydpi.drug import constitution, topology, connectivity, kappa, bcut, basak, estate, moran, moreaubroto, geary, charge, molproperty, moe, fingerprint
from pydpi import pydrug
from molvs import Standardizer


from copy import deepcopy
from os import path, getcwd, remove, system, listdir

import toolbox
import pathFolder

class Descriptors:

    def __init__(self, dcompound, writecheck=1):
        self.compound = dcompound
        loader = pydrug.PyDrug()

        # if SMILES, load using SMILES code
        if not "SMILES" in dcompound.keys():
            smile = toolbox.babelConvertSDFtoSMILE(dcompound["sdf"])
            self.compound["SMILES"] = smile

        #standardize smile code
        s = Standardizer()
        smilestandadized = s.standardize(self.compound["SMILES"])
        print smilestandadized, "Standardize SMILES - l25 liganddescriptors"
        self.compound["SMILES"] = smilestandadized

        if writecheck == 1:
            #SMILES code
            pfileSMILES = pathFolder.PR_COMPOUNDS + str(dcompound["DATABASE_ID"]) + ".smi"
            fileSMILES = open(pfileSMILES, "w")
            fileSMILES.write(self.compound["SMILES"])
            fileSMILES.close()

            #SDF input
            pfileSDF = pathFolder.PR_COMPOUNDS + str(dcompound["DATABASE_ID"]) + ".sdf"
            fileSDF = open(pfileSDF, "w")
            fileSDF.write(self.compound["sdf"])
            fileSDF.close()


        # read mol
        self.mol = loader.ReadMolFromSmile(self.compound["SMILES"])



    def get_descriptorOD1D(self):
        try: self.consti = constitution.GetConstitutional(self.mol)
        except: self.consti = {}
        self.compo = {}
        try: self.compo["nheavy"] = self.mol.GetNumHeavyAtoms()
        except: self.compo = {}
        try: self.molprop = molproperty.GetMolecularProperty(self.mol)
        except: self.molprop = {}

        # combine all 1D
        self.all1D = {}
        self.all1D.update(self.consti)
        self.all1D.update(self.compo)
        self.all1D.update(self.molprop)

        #listdesc
        self.l1D = constitution._constitutional.keys()
        self.l1D = self.l1D + ["nheavy"]
        self.l1D = self.l1D + molproperty.MolecularProperty.keys()



    def get_descriptor2D(self):
        try: self.topo = topology.GetTopology(self.mol)
        except: self.topo = {}
        try: self.connect = connectivity.GetConnectivity(self.mol)
        except: self.connect = {}
        try: self.kap = kappa.GetKappa(self.mol)
        except: self.kap = {}
        try: self.burden = bcut.GetBurden(self.mol)
        except: self.burden = {}
        try: self.basakD = basak.Getbasak(self.mol)
        except: self.basakD = {}
        try: self.est = estate.GetEstate(self.mol)
        except: self.est = {}
        try: self.moreauBurto = moreaubroto.GetMoreauBrotoAuto(self.mol)
        except: self.moreauBurto = {}
        try: self.autcormoran = moran.GetMoranAuto(self.mol)
        except: self.autcormoran = {}
        try: self.gearycor = geary.GetGearyAuto(self.mol)
        except: self.gearycor = {}
        try: self.charges = charge.GetCharge(self.mol)
        except: self.charges = {}
        try: self.MOE = moe.GetMOE(self.mol)
        except: self.MOE = {}

        # list 2D -> modified in main library !!!!
        self.l2D = topology._Topology.keys()
        self.l2D = self.l2D + connectivity._connectivity.keys()
        self.l2D = self.l2D + kappa._kapa.keys()
        self.l2D = self.l2D + bcut._bcut
        self.l2D = self.l2D + basak._basak.keys()
        self.l2D = self.l2D + estate._estate.keys()
        self.l2D = self.l2D + moreaubroto._moreaubroto.keys()
        self.l2D = self.l2D + moran._moran.keys()
        self.l2D = self.l2D + geary._geary.keys()
        self.l2D = self.l2D + charge._Charge.keys()
        self.l2D = self.l2D + moe._moe.keys()

        # combine all 2D
        self.all2D = dict()
        self.all2D.update(deepcopy(self.topo))
        self.all2D.update(deepcopy(self.connect))
        self.all2D.update(deepcopy(self.kap))
        self.all2D.update(deepcopy(self.burden))
        self.all2D.update(deepcopy(self.basakD))
        self.all2D.update(deepcopy(self.est))
        self.all2D.update(deepcopy(self.moreauBurto))
        self.all2D.update(deepcopy(self.autcormoran))
        self.all2D.update(deepcopy(self.gearycor))
        self.all2D.update(deepcopy(self.charges))
        self.all2D.update(deepcopy(self.MOE))


    def get_fingerprints(self):
        #fingerprint
        self.fingerAtomPairs = fingerprint.CalculateAtomPairsFingerprint(self.mol)
        self.fingerDaylight = fingerprint.CalculateDaylightFingerprint(self.mol)
        self.fingerEstate = fingerprint.CalculateEstateFingerprint(self.mol)
        self.fingerFP4 = fingerprint.CalculateFP4Fingerprint(self.mol)
        self.fingerMACCS = fingerprint.CalculateMACCSFingerprint(self.mol)
        self.fingerMorgan = fingerprint.CalculateMorganFingerprint(self.mol)
        self.fingerTorsion = fingerprint.CalculateTopologicalTorsionFingerprint(self.mol)


    def get_descriptor3D(self):
        """
        Compute descriptors 3D from SMILES code and generate the 3D using ligprep
        :return: dictionary of descriptors in all3D
        """

        # clean temp folder - used to compute 3D descriptors
        prtemp = pathFolder.cleanFolder()

        #temp SMILES
        pfilesmile = prtemp + "tem.smi"
        filesmile = open(pfilesmile, "w")
        filesmile.write(self.compound["SMILES"])
        filesmile.close()


        print self.compound["SMILES"]





    def writeTablesDesc(self, prresult):

        # case 1D
        if "all1D" in self.__dict__:
            if not path.exists(prresult + "1D.csv"):
                self.fil1D = open(prresult + "1D.csv", "w")
                # header
                #self.fil1D.write("ID\tSMILES\t")
                self.fil1D.write("SMILES\t")
                self.fil1D.write("\t".join(self.l1D) + "\n")
            else:
                self.fil1D = open(prresult + "1D.csv", "a")
            self.fil1D.write(self.compound['SMILES'])

            for desc1D in self.l1D:
                try: self.fil1D.write("\t" + str(self.all1D[desc1D]))
                except: self.fil1D.write("\tNA")
            self.fil1D.write("\n")
            self.fil1D.close()

        # case 2D
        if "all2D" in self.__dict__:
            if not path.exists(prresult + "2D.csv"):
                self.fil2D = open(prresult + "2D.csv", "w")
                # header
                #self.fil2D.write("ID\tSMILES\t")
                self.fil2D.write("SMILES\t")
                self.fil2D.write("\t".join(self.l2D) + "\n")
            else:
                self.fil2D = open(prresult + "2D.csv", "a")

            self.fil2D.write(self.compound['SMILES'])
            for desc2D in self.l2D:
                try:
                    self.fil2D.write("\t" + str(self.all2D[desc2D]))
                except:
                    self.fil2D.write("\tNA")
            self.fil2D.write("\n")
            self.fil2D.close()

        # case 3D - not work
        if "all3D" in self.__dict__:
            if not path.exists(prresult + "drugbank3D.csv",):
                self.fil3D = open(prresult + "drugbank3D.csv", "w")
                # header
                self.fil3D.write("ID\tSMILES\t")
                self.fil3D.write("\t".join(self.l3D) + "\n")
            else:
                self.fil3D = open(prresult + "drugbank3D.csv", "a")

            self.fil3D.write(self.compound['DRUGBANK_ID'] + "\t" + self.compound['SMILES'])
            for desc3D in self.l3D:
                try:
                    self.fil3D.write("\t" + str(self.all3D[desc3D]))
                except:
                    self.fil3D.write("\tNA")
            self.fil3D.write("\n")
            self.fil3D.close()



