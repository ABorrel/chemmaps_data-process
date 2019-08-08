from re import search



class drugbank:

    def __init__(self, psdf):
        self.psdf = psdf


    def parseAll (self):

        lout = []

        filin = open(self.psdf, "r")
        handle_read = filin.read()
        filin.close()
        #print handle_read

        l_compound = handle_read.split("$$$$\n")

        for compound in l_compound:
            dcompound = {}
            dcompound["sdf"] = compound + "$$$$"
            llines = compound.split("\n")
            i = 0
            nblines = len(llines)
            while i < nblines:
                if search("> <", llines[i]):
                    kin = llines[i].split("> <")[1]
                    kin = kin.split(">")[0]
                    valuek = llines[i+1].strip()
                    dcompound[kin] = valuek
                i += 1
            lout.append(dcompound)
        self.lc=lout
        print (len(lout))

    def clean(self):
        """TO DO"""


        return




