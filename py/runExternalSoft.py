from os import system, path, remove, chdir, getcwd, name
from re import search
from time import sleep
import subprocess 

PADEL = "/home/aborrel/softwares/padel/PaDEL-Descriptor.jar"


P_RSCRIPTS = "../R/"
R_BIN = "C://Program Files/R/R-4.2.1/bin/Rscript.exe"


def runRCMD(cmd, out = 0):

    workdir = getcwd()
    chdir(P_RSCRIPTS)
    if name == "nt":
        l_elem = cmd.split(" ")
        cmd_line = [R_BIN] + l_elem
        print(cmd_line)
        p = subprocess.Popen(cmd_line)
        (output, err) = p.communicate() 
        p.wait()
        print(err)
    else:
        print(cmd)
        system(cmd)
    chdir(workdir)


def runPadel(prin=""):
    """Input include a folder of sdf file"""
    if prin == "":
        return "ERROR - Padel Input"
    else:
        cmd = "java -jar " + PADEL + " -maxruntime 10000 -3d -dir " + str(prin) + " -file " + prin + "tem.desc"
        print (cmd)
        system(cmd)

    return prin + "tem.desc"


def babelConvertSDFtoSMILE(sdfread, clean_smi=0, rm_smi=1):

    tempsdf = open("tempsdf.sdf", "w")
    tempsdf.write(sdfread)
    tempsdf.close()

    psmile = "tempsmile.smi"

    cmd_convert = "babel tempsdf.sdf " + psmile + " 2>/dev/null"
    system(cmd_convert)

    try : filin = open (psmile, "r")
    except : return "0"
    l_Fline = filin.readlines ()
    filin.close ()
    try : smile = l_Fline[0].split ("\t")[0]
    except : return "0"

    # rewrite path in filout
    if clean_smi == 1:
        filout = open (psmile, "w")
        filout.write (str (smile))
        filout.close ()

    if rm_smi == 1:
        system("rm " + psmile)


    return smile



def babelConvertMoltoSDF(pmolin, psdfout):

    if not path.exists(psdfout):
        cmd_convert = "babel " + pmolin + " " + psdfout + " 2>/dev/null"
        system(cmd_convert)




def RDrawProjection(pfilin1D2D, pfilin3D, prout, valcor = 0.9, maxquantile=80):
    cmdplotPCA = "./draw_several_projection.R " + str(pfilin1D2D) + " " + str(pfilin3D) + " " + str(prout) + " " + str(valcor) + " " + str(maxquantile)
    runRCMD(cmdplotPCA)



def RComputeMapFiles(p_desc, type_desc, prout, corval, maxquantile):
    cmdMAP = "./compute_coords.R %s %s %s %s %s"%(p_desc, type_desc, prout, corval, maxquantile) 
    runRCMD(cmdMAP)


def molconvert(pfilin, pfilout= ""):
    """Convert with black background"""
    if pfilout == "":
        pfilout = pfilin[:-3] + "png"

    if path.exists(pfilout):
        return pfilout
    cmdconvert = "molconvert \"png:w500,Q100,#00000000\" " + pfilin + " -o " + pfilout
    print(cmdconvert)
    system(cmdconvert)
    return pfilout






