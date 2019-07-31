"""

print out info on a .r data file 

now includes .spine files!


usage: trackinfo.py filename [filename2]


basically just reads it in and runs "printinfo" from the trackdata
object


djo, 2/10


"""


# ------------------------- imports -------------------------
from __future__ import with_statement

import sys


from lib import spinedata, trackdata
from lib import util


# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    if len(sys.argv) == 1:
        print "usage: trackinfo.py filename [filename2]"
        sys.exit(2)
    
    filename = sys.argv[1]
    
    if len(sys.argv) == 2:
        # one filename = show all data:
        if "spine" in filename:
            with util.Timer("time to read") as t:
                sd = spinedata.readspinefile(filename)
            sd.printinfo()
        else:
            with util.Timer("time to read") as t:
                td = trackdata.readfile(filename)
            td.printinfo()
    elif len(sys.argv) > 2:
        filename2 = sys.argv[2]
        # going to compare data; just read, don't print
        if "spine" in filename:
            with util.Timer("time to read") as t:
                sd1 = spinedata.readspinefile(filename)
            with util.Timer("time to read") as t:
                sd2 = spinedata.readspinefile(filename2)
            print "read two files into sd1, sd2"
        else:
            with util.Timer("time to read") as t:
                td1 = trackdata.readfile(filename)
            with util.Timer("time to read") as t:
                td2 = trackdata.readfile(filename2)
            print "read two files into td1, td2"
        

