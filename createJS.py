from os import path
import math



def formatCoordinates(pcoords, pjs):

    if path.exists(pjs):
        js = open(pjs, "a")

    else:
        js = open(pjs, "w")


    # write headers
    js.write("function loadCoordinate(){\n")
    js.write("    var dcoords={")


    fcoords = open(pcoords, "r")
    lcoords = fcoords.readlines()
    fcoords.close()

    lw = []
    for coords in lcoords[1:]:
        lxyz = coords.strip().split(",")
        linew = str(lxyz[0]) + ":[" + str(lxyz[1]) + "," + str(lxyz[2]) + "," + str(lxyz[3]) + "]"
        lw.append(linew)
    js.write(",".join(lw) + "};\n")


    js.write("    return(dcoords);\n};\n\n\n")
    js.close()



def formatInfo(db, pdesc, lkinfo, pjs):


    if path.exists(pjs):
        js = open(pjs, "a")

    else:
        js = open(pjs, "w")

    # write headers
    js.write("function loadInfoDrug(){\n")
    js.write("    var infodrug={")


    # load 1D2D desc
    ddesc = {}
    if path.exists(pdesc):
        fdesc = open(pdesc, "r")
        lcpddesc = fdesc.readlines()
        fdesc.close()

        ldesc = lcpddesc[0].strip().split("\t")
        for cpd in lcpddesc[1:]:
            desc = cpd.strip().split("\t")
            name = desc[0]
            ddesc[name] = {}
            i = 1
            nbdesc = len(ldesc)
            while i < nbdesc:
                ddesc[name][ldesc[i]]= desc[i]
                i += 1

    lw = []
    # write JS
    for cpd in db.lc:
        namecpd = cpd[db.name]
        linfo = []
        for kinfo in lkinfo:
            if kinfo in cpd.keys():
                if cpd[kinfo] != "":
                    linfo.append("\"" + str(cpd[kinfo]) + "\"")
                else:
                    linfo.append("\"NA\"")
            elif namecpd in ddesc.keys() and kinfo in ddesc[namecpd].keys():
                if ddesc[namecpd][kinfo] != "":
                    linfo.append("\"" + str(ddesc[namecpd][kinfo]) + "\"")
                else:
                    linfo.append("\"NA\"")
            else:
                linfo.append("\"NA\"")
        linenew = "\"" + str(namecpd) + "\"" + ":[" + ",".join(linfo) + "]"
        lw.append(linenew)

    js.write(",".join(lw) + "};\n")

    js.write("    return(infodrug);\n};\n\n\n")
    js.close()





def extractCloseCompounds(pcoords, nneighbor, pjs):


    if path.exists(pjs):
        js = open(pjs, "a")

    else:
        js = open(pjs, "w")


    filin = open(pcoords, "r")
    lcords = filin.readlines()
    filin.close()

    dcor = {}

    for cord in lcords[1:]:
        lelem = cord.strip().split(",")
        ID = lelem[0]
        x = lelem[1]
        y = lelem[2]
        z = lelem[3]
        dcor[ID] = [float(x),float(y),float(z)]

    ddist = {}
    for ID in dcor.keys():
        ddist[ID] = {}
        for ID2 in dcor.keys():
            if ID != ID2:
                ddist[ID][ID2] = math.sqrt(sum([(xi-yi)**2 for xi, yi in zip(dcor[ID], dcor[ID2])]))

        lID = [i[0] for i in sorted(ddist[ID].items(), key=lambda x:x[1])][:nneighbor]
        ddist[ID] = lID


    js.write("function loadNeighbors(){\n    var lneighbor = {")

    lwrite = []
    for ID in ddist.keys():
        w = str(ID) + ":[" + ",".join(ddist[ID]) + "]"
        lwrite.append(w)

    js.write(",".join(lwrite) + "};\n    return(lneighbor);\n};\n")
    js.close()

