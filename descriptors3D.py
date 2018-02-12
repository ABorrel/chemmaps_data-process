import geo3D
import cpsa
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




#############################################################################


#def FormatConversion(inputmol):
#    """
#    #################################################################
#    Using Pybel to convert the smi/sdf formats to mop format!
#    #################################################################
#    """
#    # inputmol.removeh()
#    inputmol.addh()
#    inputmol.make3D(forcefield='mmff94', steps=50)  ##Gemetrical optimization
#    ##forcefields = ['uff', 'mmff94', 'ghemical']
#    # make3D(self, forcefield = "mmff94", steps = 50)
#    ##inputmol.localopt(forcefield='mmff94',steps=50)
#    outputmol = pybel.Outputfile('mop', "temp.dat", overwrite=True)
#    outputmol.write(inputmol)
#    outputmol.close()
#    f = file('temp.dat', 'r+')
#    f.write('AM1              ')
#    f.close()


#def RunMOPAC(filename):
#    """
#    #################################################################
#    Run the MOPAC using os.system
#    #################################################################
#    """

#    itest = os.system("run_mopac7" + " " + filename)
#    # time.sleep(1)
#    return itest


############################################################################
#def GetARCFile(inputmol):
#    """
    #################################################################
#    Get ARC file for each molecule
    #################################################################
#    """

#    FormatConversion(inputmol)

#    itest = RunMOPAC('temp')

#    if not itest:
#        print itest, '\t', 'Finshed successfully!'
#    else:
#        print itest, '\t', 'Failed Finished........'

#    os.remove('temp.dat')
#    os.remove('temp.log')
#    os.remove('temp.OUT')
    # os.remove('temp.arc')
#    oldpath = os.getcwd() + '/temp.arc'
#    newpath = os.getcwd() + '/temp'
#    os.rename(oldpath, newpath)


##############################################################################
#if __name__ == "__main__":
#    mol = 'C1C=CCS1'
#    mol = 'SCCC(=O)N1[C@@H](CCC1)C(=O)OCC'
#    inputmol = pybel.readstring('smi', mol)
#    GetARCFile(inputmol)
#    res = _ReadCoordinates('temp')
#    print res



#############################################################################


def get3Ddesc(psdf, psmiles, geometry=1, cpsa=1):

    ddesc = {}
    lcoordinates = parseSDFfor3D(psdf)

    if geometry == 1:
        ddesc['W3DH'] = geo3D.Calculate3DWienerWithH(lcoordinates)
        ddesc['W3D'] = geo3D.Calculate3DWienerWithoutH(lcoordinates)
        ddesc['Petitj3D'] = geo3D.CalculatePetitjean3DIndex(lcoordinates)
        ddesc['GeDi'] = geo3D.CalculateGemetricalDiameter(lcoordinates)
        ddesc['grav'] = geo3D.CalculateGravitational3D1(lcoordinates)
        ddesc['rygr'] = geo3D.CalculateRadiusofGyration(psmiles, lcoordinates)
        ddesc['Harary3D'] = geo3D.CalculateHarary3D(lcoordinates)
        ddesc['AGDD'] = geo3D.CalculateAverageGeometricalDistanceDegree(lcoordinates)
        ddesc['SEig'] = geo3D.CalculateAbsEigenvalueSumOnGeometricMatrix(lcoordinates)


        ddesc['SPAN'] = geo3D.CalculateSPANR(lcoordinates)
        ddesc['ASPAN'] = geo3D.CalculateAverageSPANR(lcoordinates)
        ddesc['MEcc'] = geo3D.CalculateMolecularEccentricity(lcoordinates)

    if cpsa == 1:

        ChargeSA = cpsa.GetChargeSA(lcoordinates, RadiusProbe=1.5, n_sphere_point=5000)

        #ddesc['ASA'] = round(cpsa.CalculateASA(ChargeSA), 3)
        #ddesc['MSA'] = round(cpsa.CalculateMSA(), 3)
        #ddesc['PNSA1'] = round(cpsa.CalculatePNSA1(ChargeSA), 3)
        #ddesc['PNSA2'] = round(cpsa.CalculatePNSA2(ChargeSA), 3)
        #ddesc['PNSA3'] = round(cpsa.CalculatePNSA3(ChargeSA), 3)
        #ddesc['PPSA1'] = round(cpsa.CalculatePPSA1(ChargeSA), 3)
        #ddesc['PPSA2'] = round(cpsa.CalculatePPSA2(ChargeSA), 3)
        #ddesc['PPSA3'] = round(cpsa.CalculatePPSA3(ChargeSA), 3)
        #ddesc['DPSA1'] = round(cpsa.CalculateDPSA1(ChargeSA), 3)
        #ddesc['DPSA2'] = round(cpsa.CalculateDPSA2(ChargeSA), 3)
        #ddesc['DPSA3'] = round(cpsa.CalculateDPSA3(ChargeSA), 3)
        #ddesc['FNSA1'] = round(cpsa.CalculateFNSA1(ChargeSA), 3)
        #ddesc['FNSA2'] = round(cpsa.CalculateFNSA2(ChargeSA), 3)
        #ddesc['FNSA3'] = round(cpsa.CalculateFNSA3(ChargeSA), 3)
        #ddesc['FPSA1'] = round(cpsa.CalculateFPSA1(ChargeSA), 3)
        #ddesc['FPSA2'] = round(cpsa.CalculateFPSA2(ChargeSA), 3)
        #ddesc['FPSA3'] = round(cpsa.CalculateFPSA3(ChargeSA), 3)
        #ddesc['WNSA1'] = round(cpsa.CalculateWNSA1(ChargeSA), 3)
        #ddesc['WNSA2'] = round(cpsa.CalculateWNSA2(ChargeSA), 3)
        #ddesc['WNSA3'] = round(cpsa.CalculateWNSA3(ChargeSA), 3)
        #ddesc['WPSA1'] = round(cpsa.CalculateWPSA1(ChargeSA), 3)
        #ddesc['WPSA2'] = round(cpsa.CalculateWPSA2(ChargeSA), 3)
        #ddesc['WPSA3'] = round(cpsa.CalculateWPSA3(ChargeSA), 3)
        #ddesc['TASA'] = round(cpsa.CalculateTASA(ChargeSA), 3)
        #ddesc['PSA'] = round(cpsa.CalculateTPSA(ChargeSA), 3)
        #ddesc['RASA'] = round(cpsa.CalculateRASA(ChargeSA), 3)
        #ddesc['RPSA'] = round(cpsa.CalculateRPSA(ChargeSA), 3)
        #desc['RNCS'] = round(cpsa.CalculateRNCS(ChargeSA), 3)
        #ddesc['RPCS'] = round(cpsa.CalculateRPCS(ChargeSA), 3)
        #ddesc['FrTATP'] = round(cpsa.CalculateFractionTATP(ChargeSA), 3)


    return ddesc


