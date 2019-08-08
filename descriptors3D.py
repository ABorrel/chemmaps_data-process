import geo3D
import cpsa3D
import rdf3D
import morse3D
import vector3d
import whim3D

import scipy
from re import search
# compute from CHEMPY
# replace MOPAC file by sdf parsing


####################
# COMPUTE PARSING  #
####################



################################################################################
class Atom:
    """
    #################################################################
    A atom class used for wrapping some properties of atoms.

    Note that Coordinates is the output of the function

    (_ReadCoordinates).
    #################################################################
    """

    def __init__(self, Coordinates):

        self.pos = vector3d.Vector3d()
        self.radius = 0.0
        self.Coordinates = Coordinates
        self.Element = ''

    def SetCoordinates(self):

        temp = self.Coordinates
        self.pos.x = float(temp[0])
        self.pos.y = float(temp[1])
        self.pos.z = float(temp[2])

    def GetCoordinates(self):

        self.SetCoordinates()

        return self.pos

    def SetElement(self):

        temp = self.Coordinates

        self.Element = temp[3]

    def GetElement(self):

        self.SetElement()

        return self.Element

    def SetRadius(self):

        radii = {'H': 1.20, 'N': 1.55, 'Na': 2.27, 'Cu': 1.40, 'Cl': 1.75, 'C': 1.70,
                 'O': 1.52, 'I': 1.98, 'P': 1.80, 'B': 1.85, 'Br': 1.85, 'S': 1.80, 'Se': 1.90,
                 'F': 1.47, 'Fe': 1.80, 'K': 2.75, 'Mn': 1.73, 'Mg': 1.73, 'Zn': 1.39, 'Hg': 1.8,
                 'Li': 1.8, '.': 1.8}

        temp = self.GetElement()

        if temp in radii.keys():
            self.radius = radii[temp]
        else:
            self.radius = radii['.']

    def GetRadius(self):

        self.SetRadius()

        return self.radius




###########################################################################

def GetAtomClassList(Coordinates):
    """
    #################################################################
    Combine all atoms in a molecule into a list form.
    Note that Coordinates is the output of the function (_ReadCoordinates).
    #################################################################
    """
    Atoms = []
    for i in Coordinates:
        atom = Atom(i)
        atom.SetCoordinates()
        atom.SetElement()
        atom.SetRadius()
        Atoms.append(atom)
    return Atoms


def GetAtomCoordinateMatrix(lcoordinates):
    """
    #################################################################
    Get the atom coordinate matrix
    #################################################################
    """
    nAtom = len(lcoordinates)
    CoordinateMatrix = scipy.zeros([nAtom, 3])
    AtomLabel = []

    for i, j in enumerate(lcoordinates):
        CoordinateMatrix[i, :] = [j[0], j[1], j[2]]
        AtomLabel.append(j[3])

    return scipy.matrix(CoordinateMatrix), AtomLabel



def parseSDFfor3D(pfilin):
    """
    Read the coordinates and charge of each atom in molecule from .sdf file.
    """

    dchargeSDF = {7:-3.0, 6:-2.0, 5:-1.0, 0:0.0, 3:1.0, 2:2.0, 1:3.0} # and 4 for radical


    latoms = []

    filin = file(pfilin, 'r')
    llines = filin.readlines()
    filin.close()

    # start at line 5 classical format
    for AtBlock in llines[4:]:
        if len(AtBlock) != 70 and len(AtBlock) != 52:
            break
        else:
            #print "-" + AtBlock[0:10] + "-"
            #print AtBlock
            # remove IND from 3D and protonation issues
            if search("IND", AtBlock):
                #print "INNN"
                continue
            X = float(AtBlock[0:10].replace(" ", ""))
            Y = float(AtBlock[10:20].replace(" ", ""))
            Z = float(AtBlock[20:30].replace(" ", ""))
            elem = AtBlock[31:34].replace(" ", "")
            charge = int(AtBlock[36:39].replace(" ", ""))
            charge = dchargeSDF[charge]

            at = [X, Y, Z, elem, charge]
            #print at
            latoms.append(at)

    return latoms


def get_atomicMass(element):

    atomicMass={'H': 1.0079, 'N': 14.0067, 'Na': 22.9897, 'Cu': 63.546, 'Cl': 35.453, 'C': 12.0107,
                 'O': 15.9994, 'I': 126.9045, 'P': 30.9738, 'B': 10.811, 'Br': 79.904, 'S': 32.065, 'Se': 78.96,
                 'F': 18.9984, 'Fe': 55.845, 'K': 39.0983, 'Mn': 54.938, 'Mg': 24.305, 'Zn': 65.39, 'Hg': 200.59,
                 'Li': 6.941, 'Co': 58.9332, "Si":28.0855, "As": 74.9216, "Te":127.6, "Sr":87.62, "Ag": 107.8682,
                "Cd": 112.411}


    return atomicMass[element]


def get_MW(lcoords, H=1):
    MW = 0
    for coords in lcoords:
        if coords[3] != "H" and H ==0:MW = MW + get_atomicMass(coords[3])
        else:MW = MW + get_atomicMass(coords[3])
    return MW






def get3Ddesc(psdf, geometry=1, cpsa=1, rdf=1, morse=1, whim=1):

    ddesc = {}
    lcoordinates = parseSDFfor3D(psdf)

    try:
        if geometry == 1:
            ddesc['W3DH'] = geo3D.Calculate3DWienerWithH(lcoordinates)
            ddesc['W3D'] = geo3D.Calculate3DWienerWithoutH(lcoordinates)
            ddesc['Petitj3D'] = geo3D.CalculatePetitjean3DIndex(lcoordinates)
            ddesc['Petitj3D'] ="NA"
            ddesc['GeDi'] = geo3D.CalculateGemetricalDiameter(lcoordinates)
            ddesc['GeDi'] = "NA"
            ddesc['grav'] = geo3D.CalculateGravitational3D1(lcoordinates)
            ddesc['rygr'] = geo3D.CalculateRadiusofGyration(lcoordinates)
            ddesc['Harary3D'] = geo3D.CalculateHarary3D(lcoordinates)
            ddesc['AGDD'] = geo3D.CalculateAverageGeometricalDistanceDegree(lcoordinates)
            ddesc['SEig'] = geo3D.CalculateAbsEigenvalueSumOnGeometricMatrix(lcoordinates)


            ddesc['SPAN'] = geo3D.CalculateSPANR(lcoordinates)
            ddesc['ASPAN'] = geo3D.CalculateAverageSPANR(lcoordinates)
            ddesc['MEcc'] = geo3D.CalculateMolecularEccentricity(lcoordinates)
    except:
        return {}

    if cpsa == 1:
        ChargeSA = cpsa3D.GetChargeSA(lcoordinates, RadiusProbe=1.5, n_sphere_point=5000)

        ddesc['ASA'] = cpsa3D.CalculateASA(ChargeSA)
        ddesc['MSA'] = cpsa3D.CalculateMSA(lcoordinates)
        ddesc['PNSA1'] = cpsa3D.CalculatePNSA1(ChargeSA)
        ddesc['PNSA2'] = cpsa3D.CalculatePNSA2(ChargeSA)
        ddesc['PNSA3'] = cpsa3D.CalculatePNSA3(ChargeSA)
        ddesc['PPSA1'] = cpsa3D.CalculatePPSA1(ChargeSA)
        ddesc['PPSA2'] = cpsa3D.CalculatePPSA2(ChargeSA)
        ddesc['PPSA3'] = cpsa3D.CalculatePPSA3(ChargeSA)
        ddesc['DPSA1'] = cpsa3D.CalculateDPSA1(ChargeSA)
        ddesc['DPSA2'] = cpsa3D.CalculateDPSA2(ChargeSA)
        ddesc['DPSA3'] = cpsa3D.CalculateDPSA3(ChargeSA)
        ddesc['FNSA1'] = cpsa3D.CalculateFNSA1(ChargeSA)
        ddesc['FNSA2'] = cpsa3D.CalculateFNSA2(ChargeSA)
        ddesc['FNSA3'] = cpsa3D.CalculateFNSA3(ChargeSA)
        ddesc['FPSA1'] = cpsa3D.CalculateFPSA1(ChargeSA)
        ddesc['FPSA2'] = cpsa3D.CalculateFPSA2(ChargeSA)
        ddesc['FPSA3'] = cpsa3D.CalculateFPSA3(ChargeSA)
        ddesc['WNSA1'] = cpsa3D.CalculateWNSA1(ChargeSA)
        ddesc['WNSA2'] = cpsa3D.CalculateWNSA2(ChargeSA)
        ddesc['WNSA3'] = cpsa3D.CalculateWNSA3(ChargeSA)
        ddesc['WPSA1'] = cpsa3D.CalculateWPSA1(ChargeSA)
        ddesc['WPSA2'] = cpsa3D.CalculateWPSA2(ChargeSA)
        ddesc['WPSA3'] = cpsa3D.CalculateWPSA3(ChargeSA)
        ddesc['TASA'] = cpsa3D.CalculateTASA(ChargeSA)
        ddesc['PSA'] = cpsa3D.CalculateTPSA(ChargeSA)
        ddesc['RASA'] = cpsa3D.CalculateRASA(ChargeSA)
        ddesc['RPSA'] = cpsa3D.CalculateRPSA(ChargeSA)
        ddesc['RNCS'] = cpsa3D.CalculateRNCS(ChargeSA)
        ddesc['RPCS'] = cpsa3D.CalculateRPCS(ChargeSA)
        ddesc['FrTATP'] = cpsa3D.CalculateFractionTATP(ChargeSA)

    if rdf ==1:
        ddesc.update(rdf3D.CalculateUnweightRDF(lcoordinates))
        ddesc.update(rdf3D.CalculateChargeRDF(lcoordinates))
        ddesc.update(rdf3D.CalculateMassRDF(lcoordinates))
        ddesc.update(rdf3D.CalculatePolarizabilityRDF(lcoordinates))
        ddesc.update(rdf3D.CalculateSandersonElectronegativityRDF(lcoordinates))
        ddesc.update(rdf3D.CalculateVDWVolumeRDF(lcoordinates))

    if morse ==1:
        ddesc.update(morse3D.CalculateUnweightMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculateChargeMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculateMassMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculateAtomicNumberMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculatePolarizabilityMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculateSandersonElectronegativityMoRSE(lcoordinates))
        ddesc.update(morse3D.CalculateVDWVolumeMoRSE(lcoordinates))

    if whim ==1:
        CoordinateMatrix, AtomLabel = GetAtomCoordinateMatrix(lcoordinates)
        ddesc['L1u'] = whim3D.GetWHIM1(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['L2u'] = whim3D.GetWHIM2(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['L3u'] = whim3D.GetWHIM3(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['Tu'] = whim3D.GetWHIM4(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['Au'] = whim3D.GetWHIM5(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['Vu'] = whim3D.GetWHIM6(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['P1u'] = whim3D.GetWHIM7(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['P2u'] = whim3D.GetWHIM8(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['Ku'] = whim3D.GetWHIM9(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['E1u'] = whim3D.GetWHIM10(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['E2u'] = whim3D.GetWHIM11(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['E3u'] = whim3D.GetWHIM12(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['Du'] = whim3D.GetWHIM13(CoordinateMatrix, AtomLabel, proname='u')
        ddesc['L1m'] = whim3D.GetWHIM1(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['L2m'] = whim3D.GetWHIM2(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['L3m'] = whim3D.GetWHIM3(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['Tm'] = whim3D.GetWHIM4(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['Am'] = whim3D.GetWHIM5(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['Vm'] = whim3D.GetWHIM6(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['P1m'] = whim3D.GetWHIM7(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['P2m'] = whim3D.GetWHIM8(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['Km'] = whim3D.GetWHIM9(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['E1m'] = whim3D.GetWHIM10(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['E2m'] = whim3D.GetWHIM11(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['E3m'] = whim3D.GetWHIM12(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['Dm'] = whim3D.GetWHIM13(CoordinateMatrix, AtomLabel, proname='m')
        ddesc['L1e'] = whim3D.GetWHIM1(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['L2e'] = whim3D.GetWHIM2(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['L3e'] = whim3D.GetWHIM3(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['Te'] = whim3D.GetWHIM4(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['Ae'] = whim3D.GetWHIM5(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['Ve'] = whim3D.GetWHIM6(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['P1e'] = whim3D.GetWHIM7(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['P2e'] = whim3D.GetWHIM8(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['Ke'] = whim3D.GetWHIM9(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['E1e'] = whim3D.GetWHIM10(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['E2e'] = whim3D.GetWHIM11(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['E3e'] = whim3D.GetWHIM12(CoordinateMatrix, AtomLabel, proname='En')
        ddesc['De'] = whim3D.GetWHIM13(CoordinateMatrix, AtomLabel, proname='En')
        try:ddesc['L1v'] = whim3D.GetWHIM1(CoordinateMatrix, AtomLabel, proname='V')
        except:ddesc['L1v'] = "NA"
        try:ddesc['L2v'] = whim3D.GetWHIM2(CoordinateMatrix, AtomLabel, proname='V')
        except: ddesc['L2v']="NA"
        try:ddesc['L3v'] = whim3D.GetWHIM3(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['L3v'] = "NA"
        try:ddesc['Tv'] = whim3D.GetWHIM4(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['Tv'] = "NA"
        try:ddesc['Av'] = whim3D.GetWHIM5(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['Av'] = "NA"
        try:ddesc['Vv'] = whim3D.GetWHIM6(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['Vv'] = "NA"
        try:ddesc['P1v'] = whim3D.GetWHIM7(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['P1v'] = "NA"
        try:ddesc['P2v'] = whim3D.GetWHIM8(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['P2v'] = "NA"
        try:ddesc['Kv'] = whim3D.GetWHIM9(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['Kv'] = "NA"
        try:ddesc['E1v'] = whim3D.GetWHIM10(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['E1v'] = "NA"
        try:ddesc['E2v'] = whim3D.GetWHIM11(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['E2v'] = "NA"
        try:ddesc['E3v'] = whim3D.GetWHIM12(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['E3v'] = "NA"
        try:ddesc['Dv'] = whim3D.GetWHIM13(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['Dv'] = "NA"
        try:ddesc['L1p'] = whim3D.GetWHIM1(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['L1p'] = "NA"
        try:ddesc['L2p'] = whim3D.GetWHIM2(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['L2p'] = "NA"
        try:ddesc['L3p'] = whim3D.GetWHIM3(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['L3p'] = "NA"
        try:ddesc['Tp'] = whim3D.GetWHIM4(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['Tp'] = "NA"
        try:ddesc['Ap'] = whim3D.GetWHIM5(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['Ap'] = "NA"
        try:ddesc['Vp'] = whim3D.GetWHIM6(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['Vp'] = "NA"
        try:ddesc['P1p'] = whim3D.GetWHIM7(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['P1p'] = "NA"
        try:ddesc['P2p'] = whim3D.GetWHIM8(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['P2p'] = "NA"
        try:ddesc['Kp'] = whim3D.GetWHIM9(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['Kp'] = "NA"
        try:ddesc['E1p'] = whim3D.GetWHIM10(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['E1p'] = "NA"
        try:ddesc['E2p'] = whim3D.GetWHIM11(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['E2p'] = "NA"
        try:ddesc['E3p'] = whim3D.GetWHIM12(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['E3p'] = "NA"
        try:ddesc['Dp'] = whim3D.GetWHIM13(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['Dp'] = "NA"
        try:ddesc['P3p'] = whim3D.GetWHIM14(CoordinateMatrix, AtomLabel, proname='alapha')
        except:
            ddesc['P3p'] = "NA"
        try:ddesc['P3u'] = whim3D.GetWHIM14(CoordinateMatrix, AtomLabel, proname='u')
        except:
            ddesc['P3u'] = "NA"
        try:ddesc['P3m'] = whim3D.GetWHIM14(CoordinateMatrix, AtomLabel, proname='m')
        except:
            ddesc['P3m'] = "NA"
        try:ddesc['P3e'] = whim3D.GetWHIM14(CoordinateMatrix, AtomLabel, proname='En')
        except:
            ddesc['P3e'] = "NA"
        try:ddesc['P3v'] = whim3D.GetWHIM14(CoordinateMatrix, AtomLabel, proname='V')
        except:
            ddesc['P3v'] = "NA"

    return ddesc


