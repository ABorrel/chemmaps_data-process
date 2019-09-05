from os import listdir, path
from random import shuffle

import runExternalSoft



def runPNG(prSMI, prPNG):

    if not path.exists(prSMI):
        print ("ERROR: CREATE CLEAN SMI FIRST")
        return


    lSMI = listdir(prSMI)
    shuffle(lSMI)
    for SMI in lSMI:
        #print prSMI + SMI
        fSMI = open(prSMI + SMI, "r")
        lelem = fSMI.readlines()
        #print(lelem)
        if len(lelem) == 0:
            continue
        lelem = lelem[0].split("\t")
        #if len(lelem) == 1:
        #    print lelem
        #    print prSMI + SMI
        if lelem[0] == "ERROR":
            continue
        else:
            if len(lelem) < 3:
                print (prSMI + SMI)
                #remove(prSMI + SMI)
                continue
            inchikey = lelem[1]
            DSSTOXid = lelem[2]
            SMIclean = lelem[0]

            if not path.exists(prPNG + inchikey + ".png"):
                pSMIclean = prPNG + inchikey + ".smi"
                fSMIcLean = open(pSMIclean, "w")
                fSMIcLean.write(SMIclean)
                fSMIcLean.close()

                runExternalSoft.molconvert(pSMIclean)
                #remove(pSMIclean)
