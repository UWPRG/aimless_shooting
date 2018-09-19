import os
import re

logfile = "commit.log"
search_term = "COMMITED"
status = None
for line in open(logfile, 'r'):
    if re.search(search_term, line):
        status = line.rstrip()
        #status = ''.join(line.splitlines())
if not status: 
    print ("None")
    print (type(status))
else: 
    print (type(status))
    print (status[-1])
