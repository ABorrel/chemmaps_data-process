import geo3D
import cpsa3D
import rdf3D
import morse3D
import vector3d

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


###########################################################################

#def _ReadCoordinates(filename="temp"):
#    """
#    #################################################################
#    Read the coordinates and charge of each atom in molecule from .arc file.
#    #################################################################
#    """
#    res = []

#    f = file(filename, 'r')
#    templine = f.readlines()
#    f.close()

#    for line in range(len(templine)):
#        if templine[line][-7:-1] == "CHARGE":
#            k = line
#            break
#
#    for i in templine[k + 4:len(templine) - 1]:
#        temp = i.split()
#        ElementCoordinate = [string.strip(temp[0]), string.strip(temp[1]),
#                             string.strip(temp[3]), string.strip(temp[5]),
#                             string.strip(temp[10])]
#        res.append(ElementCoordinate)
#
#   return res


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
        if len(AtBlock) != 70:
            break
        else:
            #print "-" + AtBlock[10:20] + "-"
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
                 'Li': 6.941, 'Co': 58.9332}


    return atomicMass[element]


def get_MW(lcoords, H=1):
    MW = 0
    for coords in lcoords:
        if coords[3] != "H" and H ==0:MW = MW + get_atomicMass(coords[3])
        else:MW = MW + get_atomicMass(coords[3])
    return MW






def get3Ddesc(psdf, geometry=0, cpsa=0, rdf=0, morse=1, whim=0):

    ddesc = {}
    lcoordinates = parseSDFfor3D(psdf)

    if geometry == 1:
        ddesc['W3DH'] = geo3D.Calculate3DWienerWithH(lcoordinates)
        ddesc['W3D'] = geo3D.Calculate3DWienerWithoutH(lcoordinates)
        ddesc['Petitj3D'] = geo3D.CalculatePetitjean3DIndex(lcoordinates)
        ddesc['GeDi'] = geo3D.CalculateGemetricalDiameter(lcoordinates)
        ddesc['grav'] = geo3D.CalculateGravitational3D1(lcoordinates)
        ddesc['rygr'] = geo3D.CalculateRadiusofGyration(lcoordinates)
        ddesc['Harary3D'] = geo3D.CalculateHarary3D(lcoordinates)
        ddesc['AGDD'] = geo3D.CalculateAverageGeometricalDistanceDegree(lcoordinates)
        ddesc['SEig'] = geo3D.CalculateAbsEigenvalueSumOnGeometricMatrix(lcoordinates)


        ddesc['SPAN'] = geo3D.CalculateSPANR(lcoordinates)
        ddesc['ASPAN'] = geo3D.CalculateAverageSPANR(lcoordinates)
        ddesc['MEcc'] = geo3D.CalculateMolecularEccentricity(lcoordinates)

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
        print "to do"

    return ddesc


