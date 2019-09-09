import psycopg2
from configparser import ConfigParser


class DBrequest:
    def __init__(self, verbose=1):
        self.dbconfig = "database.ini"
        self.conn = None
        self.verbose = verbose


    def config(self, section='postgresql'):
        parser = ConfigParser()
        parser.read(self.dbconfig)
        dparams = {}
        if parser.has_section(section):
            params = parser.items(section)
            for param in params:
                dparams[param[0]] = param[1]
        else:
            raise Exception('Section {0} not found in the {1} file'.format(section, self.dbconfig))

        self.params = dparams


    def connOpen(self):
        try:
            self.config()
            if self.verbose: print('Connecting to the PostgreSQL database...')
            self.conn = psycopg2.connect(** self.params)

        except (Exception, psycopg2.DatabaseError) as error:
            print(error)

    def connClose(self):
        if self.conn is not None:
            self.conn.close()
            if self.verbose == 1: print('Database connection closed.')


    def addElement(self, nameTable, lcoloumn, lval):
        sqlCMD = "INSERT INTO %s(%s) VALUES(%s);"%(nameTable, ",".join(lcoloumn), ",".join(["\'%s\'"%(val) for val in lval]))
        if self.verbose == 1: print(sqlCMD)
        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(sqlCMD)
                self.conn.commit()
            except (Exception, psycopg2.DatabaseError) as error:
                print(error)
        else:
            print("Open connection first")

    def extractColoumn(self, nameTable, coloumn):
        sqlCMD = "SELECT %s FROM %s LIMIT(10);" % (coloumn, nameTable)
        if self.verbose == 1: print(sqlCMD)
        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(sqlCMD)
                out = cur.fetchall()
                if self.verbose == 1: print(out)
                return out
            except (Exception, psycopg2.DatabaseError) as error:
                print(error)
                return "ERROR: " + error
        else:
            print("Open connection first")
            return "ERROR: open DB"


    def execCMD(self, cmd):
        sqlCMD = "INSERT INTO %s(%s) VALUES(%s);" % (
        nameTable, ",".join(lcoloumn), ",".join(["\'%s\'" % (val) for val in lval]))
        if self.verbose == 1: print(sqlCMD)
        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(cmd)
                #self.conn.commit()
                # print(cur)
                out = cur.fetchall()
                if self.verbose == 1: print(out)
            except (Exception, psycopg2.DatabaseError) as error:
                print(error)
        else:
            print("Open connection first")
        return



#cmd = 'select * from drugbank_chemicals  limit(10)'

#cmd = """INSERT INTO chemmaps_test(dbid, test) VALUES ('test3', 'test6')"""

#dbr = DBrequest()
#dbr.connOpen()
#dbr.addElement("chemmaps_test", ["dbid", "test"], ["tt2", "fff5"])
#out = dbr.extractColoumn("chemmap_1d2d_arr", "data_arr")
#dbr.execCMD(cmd)
#dbr.connClose()


#print(out[0][0][4])


