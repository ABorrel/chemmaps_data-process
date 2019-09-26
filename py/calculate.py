

def centroid(lcoord):

    Cx = 0.0
    Cy = 0.0
    Cz = 0.0

    for coord in lcoord:
        Cx = Cx + float(coord[0])
        Cy = Cy + float(coord[1])
        Cz = Cz + float(coord[2])


    nbchem = len(lcoord)

    lout = [Cx/nbchem, Cy/nbchem, Cz/nbchem]
    return lout

