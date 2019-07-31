"""

this script reads data from uncombined data; basically replaces
trackinfo.py 

usage: trackinfo.py chore-results-folder scalar


djo, 6/11

"""



# ------------------------- imports -------------------------
from __future__ import with_statement

import sys


from lib import spinedata, trackdata
from lib import util


# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print "usage: trackinfo.py chore-dir [scalar]"
        sys.exit(2)
    
    choredir = sys.argv[1]
    if len(sys.argv) > 2:
        scalar = sys.argv[2]
    else:
        scalar = None
    
    
    if "spine" in choredir:
        with util.Timer("time to read") as t:
            sd = spinedata.readuncombinedspinerfile(choredir)
        sd.printinfo()
    else:
        if scalar is None:
            print "usage: trackinfo.py chore-dir [scalar]"
            sys.exit(2)
        with util.Timer("time to read") as t:
            td = trackdata.readuncombinedchoredir(choredir, scalar)
        td.printinfo()



