from os import system

LIGPREP = "/opt/schrodinger2016-4/ligprep"
PADEL = "/home/aborrel/softwares/padel/PaDEL-Descriptor.jar"


def runLigprep(psmilin, forcefield="OPLS3", stereoisoster=1):

    """Maybe fix more option"""

    if forcefield == "OPLS3":
        bff = "16"
    else:
        bff = "14"

    cmd = LIGPREP + " -ismi " + psmilin + " -osd " + psmilin[0:-4] + ".sdf" + " -bff " + str(bff) + " -epik -s " + str(stereoisoster) + " -WAIT"

    print cmd
    system(cmd)

    return psmilin[0:-4] + ".sdf"

def runPadel(prin):
    """Input include a folder of sdf file"""
    cmd = "java -jar " + PADEL + " -maxruntime 100 -3d -dir " + str(prin) + " -file " + prin[0:-4] + ".desc"
    print cmd
    system(cmd)

    return prin[0:-4] + ".desc"