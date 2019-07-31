"""

this script changes filenames of data so that I can run it in my 
"test" trackers

now adjusted for pipeline branching of early 2011

djo, 2/10

"""



# ------------------------- imports -------------------------
import os
import shutil
import sys


from lib import constants as const
from lib import util


# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    if len(sys.argv) < 3:
        sys.exit("usage: retrack.py oldtracker newtracker (always branch rd)")
    
    oldtracker = sys.argv[1]
    newtracker = sys.argv[2]
    
    # new location:
    triples = util.walktracker(["rd", newtracker], location="source")
    for tr, gene, stim in triples:
        basedir = util.getdirectory("source", "rd", tr, gene, stim)
        for filename in os.listdir(basedir):
            if oldtracker in filename and filename.endswith('.r'):
                newfilename = filename.replace(oldtracker, newtracker)
                os.rename(os.path.join(basedir, filename), os.path.join(basedir, newfilename))
    
    

