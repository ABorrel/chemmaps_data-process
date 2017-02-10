from pydpi.drug import constitution, topology, connectivity, kappa, bcut, basak, estate, moran, moreaubroto, geary, charge, molproperty, moe, fingerprint
from pydpi import pydrug

#3D descriptors
from pychem import pychem, geometric, rdf, morse, whim, cpsa

from copy import deepcopy
from os import path, getcwd, remove


class Descriptors:

    def __init__(self, dcompound):
        self.compound = dcompound
        loader = pydrug.PyDrug()

        # if SMILES, load using SMILES code
        if "SMILES" in dcompound.keys():
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

        # based on chempy and generation of 3D using mopac7
        mol3d = pychem.pybel.readstring("smi", self.compound["SMILES"])
        pychem.GetARCFile(mol3d)

        # list of descriptor
        self.l3D = geometric._geometric.keys()
        self.l3D = self.l3D + cpsa._cpsa.keys()
        self.l3D = self.l3D + rdf._rdf.keys()
        self.l3D = self.l3D + morse._morse.keys()
        self.l3D = self.l3D + whim._whim.keys()

        # test if temp file is generated
        if not path.exists(getcwd()+'/temp'):
            self.geo3D = {}
            self.cpsa = {}
            self.rdf = {}
            self.morse = {}
            self.whim = {}

        else:
            self.geo3D = geometric.GetGeometric(mol3d)
            #self.cpsa = cpsa.GetCPSA()
            self.rdf = rdf.GetRDF(mol3d)
            self.morse = morse.GetMoRSE(mol3d)
            self.whim = whim.GetWHIM()


        #combine all
        self.all3D = dict()
        self.all3D.update(self.geo3D)
        #self.all3D.update(self.cpsa)
        self.all3D.update(self.rdf)
        self.all3D.update(self.morse)
        self.all3D.update(self.whim)

        # remove 3D temp
        try:remove(getcwd()+'/temp')
        except: pass


    def writeTablesDesc(self, prresult):

        # case 1D
        if "all1D" in self.__dict__:
            if not path.exists(prresult + "drugbank1D.csv"):
                self.fil1D = open(prresult + "drugbank1D.csv", "w")
                # header
                self.fil1D.write("ID\tSMILES\t")
                self.fil1D.write("\t".join(self.l1D) + "\n")
            else:
                self.fil1D = open(prresult + "drugbank1D.csv", "a")
            self.fil1D.write(self.compound['DRUGBANK_ID'] + "\t" + self.compound['SMILES'])
            for desc1D in self.l1D:
                try: self.fil1D.write("\t" + str(self.all1D[desc1D]))
                except: self.fil1D.write("\tNA")
            self.fil1D.write("\n")
            self.fil1D.close()

        # case 2D
        if "all2D" in self.__dict__:
            if not path.exists(prresult + "drugbank2D.csv"):
                self.fil2D = open(prresult + "drugbank2D.csv", "w")
                # header
                self.fil2D.write("ID\tSMILES\t")
                self.fil2D.write("\t".join(self.l2D) + "\n")
            else:
                self.fil2D = open(prresult + "drugbank2D.csv", "a")

            self.fil2D.write(self.compound['DRUGBANK_ID'] + "\t" + self.compound['SMILES'])
            for desc2D in self.l2D:
                try:
                    self.fil2D.write("\t" + str(self.all2D[desc2D]))
                except:
                    self.fil2D.write("\tNA")
            self.fil2D.write("\n")
            self.fil2D.close()

        # case 3D
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



