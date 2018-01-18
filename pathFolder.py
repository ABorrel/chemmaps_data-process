from os import listdir, remove, makedirs, path
from shutil import rmtree

#PR_REF = "/home/aborrel/ChemMap/"
#PR_RESULT = "/home/aborrel/ChemMap/results/"
#PR_TEMP3D = "/home/aborrel/ChemMap/results/temp3D/"
#PR_COMPOUNDS = "/home/aborrel/ChemMap/results/compounds/"
#PR_ANALYSIS = "/home/aborrel/ChemMap/results/analysis/"





def cleanFolder(prin):

    lfiles = listdir(prin)
    if len(lfiles) != 0:
        for filin in lfiles:
            # problem with folder
            try: remove(prin + filin)
            except: rmtree(prin + filin)

    return prin


def createFolder(prin, clean=0):

    if not path.exists(prin):
        makedirs(prin)

    if clean == 1:
        cleanFolder(prin)

    return prin


#def analyses(psub):

#    if psub == "":
#        return PR_ANALYSIS
#    else:
#        try: makedirs(PR_ANALYSIS + psub + "/")
#        except: pass

#    return PR_ANALYSIS + psub + "/"



