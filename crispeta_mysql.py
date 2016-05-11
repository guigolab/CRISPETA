from config import *
from func import *

parser = argparse.ArgumentParser(description = "")
parser.add_argument('-i', '--input',
                    required = True,
                    dest = 'infile',
                    action = 'store',
                    help = 'Path to input file')

parser.add_argument('-u', '--user',
					dest = 'user',
					action = 'store',
					default = 'crispeta',
					help = 'MySQL user name')

parser.add_argument('-p', '--password',
					dest = 'pwd',
					action = 'store',
					default = 'pwd',
					help = 'MySQL user password')

parser.add_argument('-host', '--host',
					dest = 'host',
					action = 'store',
					default = 'localhost',
					help = 'MySQL host')

parser.add_argument('-db', '--database',
					dest = 'db',
					action = 'store',
					default = 'crispeta',
					help = 'MySQL crispeta database name')

parser.add_argument('-t', '--table',
					dest = 'table',
					action = 'store',
					default = 'crispeta',
					help = 'MySQL table name')

op = parser.parse_args()

config['mysql']['user'] = op.user if op.user != config['mysql']['user'] else op.user
config['mysql']['passwd'] = op.pwd if op.pwd != config['mysql']['passwd'] else op.pwd
config['mysql']['host'] = op.host if op.host != config['mysql']['host'] else op.host

config['mysql_db']['user'] = op.user if op.user != config['mysql_db']['user'] else op.user
config['mysql_db']['passwd'] = op.pwd if op.pwd != config['mysql_db']['passwd'] else op.pwd
config['mysql_db']['host'] = op.host if op.host != config['mysql_db']['host'] else op.host
config['mysql_db']['db'] = op.db if op.db != config['mysql_db']['db'] else op.db
table = op.table if op.table != 'crispeta' else op.table


conn = get_conn(config['mysql'])													#connect to mysql
create_database(conn, op.db)														#creates database

conn = get_conn(config['mysql_db'])													#connects to database
create_table(conn, op.table, config['colnames'])									#creates table

print "Loading data to MySQL. This step can take several minutes"
insert_into_table(conn, op.infile, op.table)                            #loads data into table
 
sys.exit("Table '"+op.table+"' has been created succesfully on data base '"+op.db+"'.")   #exit
