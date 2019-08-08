

def centroid(dcoord):

    Cx = 0.0
    Cy = 0.0
    Cz = 0.0

    for chemID in dcoord.keys():
        Cx = Cx + float(dcoord[chemID][0])
        Cy = Cy + float(dcoord[chemID][1])
        Cz = Cz + float(dcoord[chemID][2])


    nbchem = len(list(dcoord.keys()))

    lout = [Cx/nbchem, Cy/nbchem, Cz/nbchem]
    return lout

