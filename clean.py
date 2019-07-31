"""

delete all data that I produce in a given tracker

usage: python clean.py branch tracker [--sure]

-- if tracker is a production tracker, you must provide the "--sure" flag!
-- if tracker == dev: deletes data in all test trackers (see below for list)



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
productiontrackers = [
    "t1",
    "t2",
    "t3",
    "t4",
    "t5",
    "t6",
    ]

testtrackers = [
    "t1-don-test", 
    "t2-don-test",
    "t3-don-test",
    ]

locationlist = [
    "analysis",
    "derived",
    "histograms",
    "plots",
    "statistics",
    "tables",
    ]
folderlist = [const.locations[loc] for loc in locationlist]


# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    if len(sys.argv) < 3:
        print "usage: clean.py branch tracker"
        sys.exit(2)
    
    branch = sys.argv[1]
    tracker = sys.argv[2]
    if tracker == "dev":
        trackerlist = testtrackers
    elif tracker in productiontrackers:
        if "--sure" in sys.argv:
            trackerlist = [tracker]
        else:
            print "you must use the --sure flag for production trackers!"
            sys.exit(1)
    else:
        trackerlist = [tracker]
    
    # delete all files in the above directories
    for tracker in trackerlist:
        for location in locationlist:
            folderpath = util.getdirectory(location, branch, tracker)
            if os.path.exists(folderpath):
                filelist = os.listdir(folderpath)
                print "deleting %d items from %s" % (len(filelist), folderpath)
                for filename in filelist:
                    filepath = os.path.join(folderpath, filename)
                    if os.path.isdir(filepath):
                        shutil.rmtree(filepath)
                    else:
                        os.remove(filepath)
    
