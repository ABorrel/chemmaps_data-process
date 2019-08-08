from os import path
import math
import toolbox


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



def formatInfo(db, pdesc, lkinfo, pjs, prout):


    if path.exists(pjs):
        js = open(pjs, "a")

    else:
        js = open(pjs, "w")

    # write headers
    js.write("function loadInfoDrug(){\n")
    js.write("    var infodrug={")


    # load 1D2D desc
    if path.exists(pdesc):
        ddesc = toolbox.loadMatrixToDict(pdesc)
    else:
        return

    lw = []
    # write JS
    for cpd in db.lc:
        namecpd = cpd[db.name]
        linfo = []
        for kinfo in lkinfo:
            if kinfo in list(cpd.keys()):
                if cpd[kinfo] != "":
                    linfo.append("\"" + str(cpd[kinfo]) + "\"")
                else:
                    linfo.append("\"NA\"")
            elif namecpd in list(ddesc.keys()) and kinfo in list(ddesc[namecpd].keys()):
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


    pinfo = prout + "tableinfo.csv"
    finfo = open(pinfo, "w")
    finfo.write("ID\t" + "\t".join(lkinfo) + "\n")

    for cpd in db.lc:
        namecpd = cpd[db.name]
        linfo = []
        for kinfo in lkinfo:
            if kinfo in list(cpd.keys()):
                if cpd[kinfo] != "":
                    linfo.append(str(cpd[kinfo]))
                else:
                    linfo.append("NA")
            elif namecpd in list(ddesc.keys()) and kinfo in list(ddesc[namecpd].keys()):
                if ddesc[namecpd][kinfo] != "":
                    linfo.append(str(ddesc[namecpd][kinfo]))
                else:
                    linfo.append("NA")
            else:
                linfo.append("NA")

        finfo.write("%s\t%s\n"%(namecpd, "\t".join(linfo)))

    finfo.close()





def extractCloseCompounds(pcoords, nneighbor, pjs, prout):



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


    #write in table
    ptableneighbor = prout + "tableNeighbor.csv"
    ftable = open(ptableneighbor, "w")
    ftable.write("ID\tNeighbors\n")
    for ID in ddist.keys():
        ftable.write("%s\t%s\n"%(ID, " ".join(ddist[ID])))
    ftable.close()

