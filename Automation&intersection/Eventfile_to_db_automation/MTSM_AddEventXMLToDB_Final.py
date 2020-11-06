import sqlite3
import time
from lxml import etree
import xml.etree.cElementTree as ET
import pandas as pd
import sys, traceback
import os

args = sys.argv[1:]
arglen = len(sys.argv)
print ("Number of arguments: <"+str(arglen)+">")

if arglen >= 3:
     eventXML=args[0]
     dbfile = args[1]

else:
    print("One of the arguments is missing!")
    print("Syntax of command is:")
    print ("> python DFGenTotalSum.py [EVENT_XML_FILE] [DB_FILE] ")
    exit()
start_time = time.time()
print ("\nCurrent working directory")
print(os.getcwd())
os.chdir(os.path.dirname(dbfile))
print ("New working directory3 ")
print(os.getcwd())

print ("\n")
print ("Selected Event file: <"+eventXML+">")
print ("Selected Database file: <"+dbfile+">")


print ("\n**** Please wait while the Database file is generated!! ***")

# input tables
eventTblnm = "event"

# Connect to the database file
conn = sqlite3.connect(dbfile)
c = conn.cursor()

try: 
#initial cleanup
    sqlqry = "DROP TABLE IF EXISTS " + eventTblnm
    print("Running query <"+sqlqry+">. Please wait!")
    c.execute(sqlqry)
  
    print("\nParsing XML file <"+eventXML+">. Please wait!")
    tree = ET.parse(eventXML) 
    root = tree.getroot()
#dfcols = ['person', 'vehicle', 'time', 'link', 'type']
#event_df = pd.DataFrame(columns=dfcols)
    knt = 0
    data=[]

    for i in root.getiterator():
        p_person = i.get('person')
        p_veh = i.get('vehicle')
        p_link =  i.get('link')
  
        if ((p_person is not None) or (p_veh is not None)):
            if (p_link is not None): 
                knt = knt+1
                p_time =  i.get('time')  
                p_type = i.get('type')
                data.append((p_person, p_veh, p_time, p_link, p_type))
      #newrec = pd.DataFrame({'person':[p_person], 'vehicle':[p_veh], 'time':[p_time], 'link':[p_link], 'type':[p_type]})
                if (knt % 5000 == 5):        
        #print (knt," : ",newrec)
                    print(knt," : ",p_person, ": ", p_veh, ": " , p_time, ": " ,p_link, ": " ,p_type) 
      #event_df = event_df.append(newrec)
      
    event_df = pd.DataFrame(data ,columns=['person', 'vehicle', 'time', 'link', 'type'])
    
  #if (knt>10):
  #  break

    event_df = event_df.astype({"person": str})
    event_df = event_df.astype({"vehicle": str})
#convert times from sec to minutes
#event_df['depTime']   = event_df['depTime'] / 60.0
      
    print("\nShape of event_df")
    print(event_df.shape)
    print(event_df.tail(10))

#write table in database 
    event_df.to_sql(eventTblnm, conn, index=False)
    print("\nWriting of SQL table <"+eventTblnm+"> completed!")

except: 
  traceback.print_exc(file=sys.stdout) 
  print ("**** Error encountered. Printing traceback! ****")
  traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
  print("*** Print_exception:")
  # exc_type below is ignored on 3.5 and later
  traceback.print_exception(exc_type, exc_value, exc_traceback,
                              limit=2, file=sys.stdout)
finally:
  #clean up
  del [[event_df]] 
  event_df = pd.DataFrame()
  
  # Close database file
  conn.close()
  
  print ("**** Terminating MTSM_AddEventXMLToDB.pyy ****")
  print("--- %s seconds ---" % (time.time() - start_time))

