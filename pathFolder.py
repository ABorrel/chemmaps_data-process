from os import listdir, remove

PR_REF = "/home/aborrel/ChemMap/"
PR_RESULT = "/home/aborrel/ChemMap/results/"
PR_TEMP3D = "/home/aborrel/ChemMap/results/temp3D/"
PR_COMPOUNDS = "/home/aborrel/ChemMap/results/compounds/"



def cleanFolder(prin=PR_TEMP3D):
    lfiles = listdir(prin)
    if len(lfiles) != 0:
        for filin in lfiles:
            # problem with folder
            remove(prin + filin)
    return prin

