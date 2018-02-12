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

    def __init__(self, dcompound, logfile, prout, dfiles, namek="DATABASE_ID"):
        self.compound = dcompound
        self.prout = prout
        self.namek = namek
        self.descFiles = dfiles
        loader = pydrug.PyDrug()

        # folder with SDF and SMI
        prCpdSmi = pathFolder.createFolder(self.prout + "SMI/")
        prCpdSDF = pathFolder.createFolder(self.prout + "SDF/")

        psdf = prCpdSDF + self.compound[namek] + ".sdf"
        psmi = prCpdSmi + self.compound[namek] + ".smi"

        # already clean and SMILES desalt
        if path.exists(psdf) and path.exists(psmi):
            print "Already clean", self.compound[namek]
            fsmile = open(psmi, "r")
            smile = fsmile.readlines()[0].strip()
            fsmile.close()
            self.compound["SMILES"] = smile
            self.log = "OK"
            self.mol = loader.ReadMolFromSmile(self.compound["SMILES"])
        else:

            if not "SMILES" in dcompound.keys():
                try:
                    smile = runExternalSoft.babelConvertSDFtoSMILE(dcompound["sdf"])
                    self.compound["SMILES"] = smile
                    #print smile
                except:
                    print "ERROR INPUT SDF - l33"
                    self.log = "ERROR"
                    try:logfile.write(
                        self.compound[namek] + "\t---\tERROR-SDF ORIGINAL INPUT\n")
                    except:pass

                    return


            #Standardize smile code
            try: smilestandadized = standardize_smiles(self.compound["SMILES"])
            except:
                logfile.write(self.compound[namek] + "\t" + str(self.compound["SMILES"]) + "\tERROR-SMILES INPUT"
                                                                                            "\n")
                self.log = "ERROR"
                return

            #Standardize using molvs (http://molvs.readthedocs.io/en/latest/api.html#molvs-fragment)
            s = Standardizer()
            mol = Chem.MolFromSmiles(smilestandadized)
            molstandardized = s.standardize(mol)
            smilestandadized = Chem.MolToSmiles(molstandardized)

            # remove salt
            # 1.defaultre
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
                try:lelem.remove("") #case of bad smile
                except:pass
                if len(lelem) == 1:
                    smilesclean = str(lelem[0])
                else:
                    # 4. Fragments
                    #Case of fragment -> stock in log file, check after to control
                    logfile.write(self.compound[namek] + "\t" + str(self.compound["SMILES"]) + "\tFRAGMENT IN INPUT"
                                                                                                       "\n")
                    print ".".join(lelem), " - FRAGMENTS - l66"
                    self.log = "ERROR"
                    return
            else:
                pass


            #print self.compound["SMILES"], "SMILES IN - l25 liganddescriptors"
            print smilesclean, "SMILES without salt and standardized"

            # case where only salt are included
            if smilesclean == "":
                logfile.write(self.compound[namek] + "\t" + str(self.compound["SMILES"]) + "\tEMPTY SMILES AFTER "
                                                                                               "STANDARDIZATION\n")
                print "EMPTY SMILES AFTER STANDARDIZATION - l84"
                print self.compound[self.namek]
                self.log = "ERROR"
                return

            self.compound["SMILES"] = smilesclean
            self.log = "OK"


            #write clean SMILES and split sdf
            pfileSMILES = prCpdSmi + str(dcompound[namek]) + ".smi"
            fileSMILES = open(pfileSMILES, "w")
            fileSMILES.write(self.compound["SMILES"])
            fileSMILES.close()

            #SDF input
            pfileSDF = prCpdSDF + str(dcompound[namek]) + ".sdf"
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
        if not "allDesc" in dir(self):
            self.allDesc = dict()
        self.allDesc.update(deepcopy(self.consti))
        self.allDesc.update(deepcopy(self.compo))
        self.allDesc.update(deepcopy(self.molprop))
        self.allDesc.update(deepcopy(self.topo))
        self.allDesc.update(deepcopy(self.connect))
        self.allDesc.update(deepcopy(self.kap))
        self.allDesc.update(deepcopy(self.burden))
        self.allDesc.update(deepcopy(self.basakD))
        self.allDesc.update(deepcopy(self.est))
        self.allDesc.update(deepcopy(self.moreauBurto))
        self.allDesc.update(deepcopy(self.autcormoran))
        self.allDesc.update(deepcopy(self.gearycor))
        self.allDesc.update(deepcopy(self.charges))
        self.allDesc.update(deepcopy(self.MOE))



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
        pr3DSDF = pathFolder.createFolder(self.prout + "SDF3D/")

        # clean temp folder - used to compute 3D descriptors
        prtemp = pathFolder.createFolder(self.prout + "temp3D/", clean = 1)
        psdf3Dout = pr3DSDF + self.compound[self.namek] + ".sdf"

        # control if chemical exists
        #print psdf3Dout
        if path.exists(psdf3Dout):
            return pr3DSDF

        # temp SMILES
        pfilesmile = prtemp + "tem.smi"
        filesmile = open(pfilesmile, "w")
        filesmile.write(self.compound["SMILES"])
        filesmile.close()

        pdesc = ""
        # run ligprep
        if not path.exists(psdf3Dout):
            psdf3D = runExternalSoft.runLigprep(psmilin=pfilesmile)
            #ffff
            # case error in ligprep
            if not path.exists(psdf3D) or path.getsize(psdf3D) == 0:
                log.write(self.compound[self.namek] + "\t" + self.compound["SMILES"] + "\t" + psdf3D)
            else:
                psdf3Dout = toolbox.selectMinimalEnergyLigPrep(psdfin=psdf3D,
                                                               psdfout=psdf3Dout)
                # take only best energy
                pathFolder.cleanFolder(prtemp)
        return pr3DSDF




    def writeTablesDesc(self, kfile):

        # case 1D2D
        if not kfile in self.descFiles.keys():
             return
        else:
            pfilin = self.descFiles[kfile]

        if not path.exists(pfilin):
            filin = open(pfilin, "w")
            # header
            if kfile == "1D2D":
                filin.write("ID\tSMILES\t")
                filin.write("\t".join(self.l1D2D) + "\n")
            # does not work with 3D descriptors, table from Padel formated already
            #elif kfile == "3D":
            #    filin.write("ID\tSMILES\t")
            #    filin.write("\t".join(self.l3D) + "\n")
            #filin.close()

        # write desc
        filin = open(pfilin, "a")
        filin.write(self.compound[self.namek] + "\t" + self.compound['SMILES'])

        if kfile == "1D2D":
            ldescs = self.l1D2D
        #elif kfile == "3D":
        #    ldescs = self.l3D

        for desc in ldescs:
            try:
                filin.write("\t" + str(self.allDesc[desc]))
            except:
                filin.write("\tNA")
        filin.write("\n")
        filin.close()



def get_descriptor3DPadel(pr3DSDF, pdesc3D):


    pPadel = runExternalSoft.runPadel(pr3DSDF)
    dpadel = toolbox.parsePadelOut(pPadel)
    filout = open(pdesc3D, "w")

    # write descriptor
    ldesc = toolbox.parsePadelOut().keys()
    filout.write("ID" + "\t".join(ldesc) + "\n")
    for compound in dpadel.keys():
        #print dpadel[compound]
        filout.write(compound)
        for desc in ldesc:
            if dpadel[compound][desc] == "":
                filout.write("\tNA")
            else:
                filout.write("\t" + str(dpadel[compound][desc]))
        filout.write("\n")
    filout.close()

    return pdesc3D



def get_descriptor3Down():

    # geometric
    #CPSA
    #RDF
    #MoRSE
    #WHIM


    return

