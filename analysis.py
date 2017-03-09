import runExternalSoft
import pathFolder



def PCAplot(pdesc1D, pdesc2D, pdesc3D):

    prout = pathFolder.analyses("PCAs")
    runExternalSoft.PCAplot(pdesc1D, pdesc2D, pdesc3D, prout)
