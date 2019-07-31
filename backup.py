"""

duplicate all output data for a given branch, tracker, gene, and stim

usage: python backup.py branch tracker geneeff stimno


meant to create backup data for comparison testing


djo, 1/10


"""

# ------------------------- imports -------------------------
import os
import shutil
import sys


# from twocof lib:
import lib.constants as const
from lib import util


# ------------------------- constants -------------------------
suffix = "-safe"

locationlist = [
    "analysis",
    "derived",
    "histograms",
    "plots",
    "statistics",
    "tables",
    ]

productiontrackers = [
    "t1",
    "t2",
    "t3",
    "t4",
    "t5",
    "t6",
    ]



# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    if len(sys.argv) < 5:
        print "usage: backup.py branch tracker geneeff stimno"
        sys.exit(2)
    
    branch = sys.argv[1]
    tracker = sys.argv[2]
    geneeff = sys.argv[3]
    stimno = sys.argv[4]
    
    if tracker in productiontrackers and "--sure" not in sys.argv:
        print "you must use the --sure flag for production trackers!"
        sys.exit(1)
    
    for location in locationlist:
        oldpath = util.getdirectory(location, branch, tracker, geneeff, stimno)
        
        newpath = oldpath + suffix
        
        if os.path.exists(oldpath):
            print "copying %s" % location
            shutil.copytree(oldpath, newpath)
        else:
            print "%s has no files; skipping" % location
        

