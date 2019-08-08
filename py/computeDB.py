import chemical
import loadDB
import pathFolder

from os import path, remove

def computeDesc(psdf, prdesc, Desc1D2D=1, generation3D = 1, Desc3D=1, control=0, namek="DATABASE_ID", soft3Dgeneration = "RDKit"):

    # folder
    dout = {}
    dout["1D2D"] = prdesc + "1D2D.csv"
    dout["3D"] = prdesc + "3D.csv"

    # shortcut
    if path.exists(prdesc + "3D.csv") and path.exists(prdesc + "1D2D.csv") and Desc3D == 0 and Desc1D2D == 0:
        return dout



    if control == 2:
        try:remove(dout["1D2D"])
        except:pass

        try:remove(dout["3D"])
        except:pass


    if control == 0:
        if Desc1D2D == 1:
            try: remove(dout["1D2D"])
            except: pass
        if Desc3D == 1:
            try: remove(dout["3D"])
            except: pass


    if Desc1D2D == 1:
        ldesc1D2D = chemical.getLdesc("1D2D")
        dout["1D2D"] = open(dout["1D2D"], "w")
        dout["1D2D"].write("ID\t" + "\t".join(ldesc1D2D) + "\n")

    if Desc3D == 1:
        ldesc3D = chemical.getLdesc("3D")
        dout["3D"] = open(dout["3D"], "w")
        dout["3D"].write("ID\t" + "\t".join(ldesc3D) + "\n")

    # formate database
    DB = loadDB.sdfDB(psdf, namek, prdesc)
    DB.parseAll()
    DB.renameHeader()# change first line of sdf to write the good name

    print (len(DB.lc))

    # folder for the log
    prlog = prdesc + "log/"
    pathFolder.createFolder(prlog)

    i = 0
    for compound in DB.lc:
        nameChem = compound[namek]
        smiles = compound["SMILES"]
        print(i, nameChem)  # for verbose

        prSMIclean = prdesc + "SMIclean/"
        pathFolder.createFolder(prSMIclean)

        # prepare ligand
        chem = chemical.chemical(nameChem, smiles)
        err = chem.prepareChem(prSMIclean)



        #Descriptor 1D and 2D#
        ######################
        if Desc1D2D == 1:
            prDescByChem = prdesc + "1D2DbyChem/"
            pathFolder.createFolder(prDescByChem)
            chem.compute1D2DDesc(prDescByChem)
            err = chem.writeTablesDesc(prDescByChem, "1D2D")

            if err == 0:
                chem.writeDesc(ldesc1D2D, dout["1D2D"], ["1D2D"])


        # descriptor 3D #
        #################
        if Desc3D == 1:
            prSDF3D = prdesc + "SDF3D/"
            pathFolder.createFolder(prSDF3D)
            if generation3D == 1:# generation from the SMILES
                chem.generate3DFromSMILES(prSDF3D, software=soft3Dgeneration)
            else:
                pSDF3D = prSDF3D + nameChem + ".sdf"
                if path.exists(pSDF3D):
                    chem.log = chem.log + "Existing 3D\n"
                    chem.psdf3D = pSDF3D
                elif generation3D == 3:# case where we have already 3D structures
                    fSDF3D = open(pSDF3D, "w")
                    fSDF3D.write(compound["sdf"])
                    fSDF3D.close()
                    chem.psdf3D = pSDF3D
                else:
                    print ("No 3D structure properly computed")
                    err = 1

            pr3DDesc = prdesc + "3DbyChem/"
            pathFolder.createFolder(pr3DDesc)

            err = chem.compute3DDesc(pr3DDesc)
            chem.writeTablesDesc(pr3DDesc, "3D")

            if err == 0:
                chem.writeDesc(ldesc3D, dout["3D"], ["3D"])


        if err ==1:
            chem.writelog(prlog)

        i += 1

    if Desc3D == 1:
        dout["3D"].close()
        dout["3D"] = prdesc + "3D.csv"

    if Desc1D2D == 1:
        dout["1D2D"].close()
        dout["1D2D"] = prdesc + "1D2D.csv"

    return dout





