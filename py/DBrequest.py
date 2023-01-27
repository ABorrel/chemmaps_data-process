import psycopg2
from configparser import ConfigParser

from os import path

class DBrequest:
    def __init__(self, verbose=0):
        self.dbconfig = path.realpath(path.dirname(__file__)) + "\\" + "database.ini"
        self.conn = None
        self.verbose = verbose

    def config(self, section='postgresql'):
        parser = ConfigParser()
        parser.read(self.dbconfig)
        dparams = {}
        if parser.has_section(section):
            params = parser.items(section)
            for param in params:
                if param[0] == "schema":
                    dparams["options"] = "-c search_path=dbo," + param[1]
                else:
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

    def addElement(self, nameTable, lcoloumn, lval, openDB=1):
        if openDB == 1:self.connOpen()
        sqlCMD = "INSERT INTO %s(%s) VALUES(%s);"%(nameTable, ",".join(lcoloumn), ",".join(["\'%s\'"%(val) for val in lval]))
        if self.verbose == 1: print(sqlCMD)
        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(sqlCMD)
                self.conn.commit()
                if openDB == 1:self.connClose()
            except (Exception, psycopg2.DatabaseError) as error:
                print(error)
                if openDB == 1:self.connClose()
        else:
            print("Open connection first")

    def runCMDaddElement(self, sqlCMD):
        self.connOpen()
        if self.verbose == 1: print(sqlCMD)
        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(sqlCMD)
                self.conn.commit()
                self.connClose()
            except (Exception, psycopg2.DatabaseError) as error:
                print(error)
                self.connClose()
        else:
            print("Open connection first")       

    def addMultipleElem(self, nameTable, lcoloumn, lval):
        """Same function than add element without open and close DB
        Have to be open and close before"""
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

    def extractColoumn(self, nameTable, coloumn, condition=""):
        sqlCMD = "SELECT %s FROM %s %s;" % (coloumn, nameTable, condition)
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

    def getRow(self, table, condition):

        self.connOpen()
        sqlCMD = "SELECT * FROM %s WHERE %s;" % (table, condition)
        if self.verbose == 1: print(sqlCMD)
        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(sqlCMD)
                out = cur.fetchall()
                if self.verbose == 1: print(out)
                self.connClose()
                return out
            except (Exception, psycopg2.DatabaseError) as error:
                print(error)
                self.connClose()
                return error
        else:
            self.connClose()
            print("Open connection first")

    def updateTable(self, cmdSQL):
        if self.verbose == 1: print(cmdSQL)
        self.connOpen()
        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(cmdSQL)
                self.conn.commit()
            except (Exception, psycopg2.DatabaseError) as error:
                self.connClose()
                print(error)
                return "Error"
        else:
            print("Open connection first")
            return None
        self.connClose()
        return 1

    def updateTable_run(self, cmdSQL):
        if self.verbose == 1: print(cmdSQL)
        #self.connOpen()
        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(cmdSQL)
                self.conn.commit()
            except (Exception, psycopg2.DatabaseError) as error:
                #self.connClose()
                print(error)
                return "Error"
        else:
            print("Open connection first")
            return None
        #self.connClose()
        return 1

    def updateTableMultiple(self, cmdSQL):
        if self.verbose == 1: print(cmdSQL)
        #self.connOpen()
        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(cmdSQL)
                self.conn.commit()
            except (Exception, psycopg2.DatabaseError) as error:
                #self.connClose()
                print(error)
                return "Error"
        else:
            print("Open connection first")
            return None
        #self.connClose()
        return 1

    def execCMDrun(self, cmdSQL):
        if self.verbose == 1: print(cmdSQL)
        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(cmdSQL)
                #self.conn.commit()
                # print(cur)
                try:out = cur.fetchall()
                except: out = 0
                if self.verbose == 1: print(out)
            except (Exception, psycopg2.DatabaseError) as error:
                print(error)
                return "Error"
        else:
            print("Open connection first")
            return None
        return out

    def execCMD(self, cmdSQL):
        if self.verbose == 1: print(cmdSQL)
        self.connOpen()
        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(cmdSQL)
                #self.conn.commit()
                # print(cur)
                try:out = cur.fetchall()
                except: out = 0
                if self.verbose == 1: print(out)
            except (Exception, psycopg2.DatabaseError) as error:
                self.connClose()
                print(error)
                return "Error"
        else:
            print("Open connection first")
            return None
        self.connClose()
        return out



#cmd = 'select * from drugbank_chemicals  limit(10)'

#cmd = """INSERT INTO chemmaps_test(dbid, test) VALUES ('test3', 'test6')"""

#dbr = DBrequest()
#dbr.connOpen()
#dbr.connClose()
#dbr.connOpen()
#dbr.addElement("chemmaps_test", ["dbid", "test"], ["tt2", "fff5"])
#out = dbr.extractColoumn("chemmap_1d2d_arr", "data_arr")
#dbr.execCMD(cmd)
#dbr.connClose()


#print(out[0][0][4])


