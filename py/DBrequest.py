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

    def execCMD(self, cmd):

        if self.conn != None:
            try:
                cur = self.conn.cursor()
                cur.execute(cmd)
                print(cur)
                out = cur.fetchall()
                if self.verbose == 1: print(out)
            except (Exception, psycopg2.DatabaseError) as error:
                print(error)
        else:
            print("Open connection first")



cmd = 'select * from drugbank_chemicals  limit(10)'

cmd = 'INSERT INTO drugbank_chemicals (drugbank_id, smiles_origin, smiles_clean, inchikey, qsar_ready)\n VALUES ("test1", "test2", "test3", "test4", false); '

dbr = DBrequest()
dbr.connOpen()
dbr.execCMD(cmd)
dbr.connClose()

