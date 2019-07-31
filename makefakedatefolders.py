"""

this script takes a current test data set and creates empty directories
in chore-results that will allow the new "per run stats" to be run


djo, 4/11


"""


# ------------------------- imports -------------------------
import optparse
import os
import shutil
import sys


from lib import constants as const
from lib import trackdata
from lib import util






# ------------------------- constants -------------------------

branch = "rd"



# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    # parse input
    usage = "usage: %prog tracker geneeff stimno1 [options] (always branch rd)"
    parser = optparse.OptionParser(usage=usage, version="%%prog %s" % const.__version__)
    parser.add_option("--verbose", action="store_true", dest="verbose", default=False,
        help="increase level of output")
    (options, args) = parser.parse_args()
    
    
    if len(args) < 3:
        parser.print_usage()
        sys.exit(2)
    
    tracker = args[0]
    geneeff = args[1]
    stimno = args[2]
    
    os.umask(const.umask)
    
    # just need one .r file, doesn't matter which:
    sourcefilename = "%s@%s@%s.x.r" % (geneeff, tracker, stimno)
    sourcedir = util.getdirectory("source", branch, tracker, geneeff, stimno)
    sourcepath = os.path.join(sourcedir, sourcefilename)
    
    td = trackdata.readfile(sourcepath)
    
    
    choredir = util.getdirectory("choreography", branch, tracker, geneeff, stimno)
    
    for datestamp in td.getdatestamps():
        datefolder = os.path.join(choredir, datestamp)
        
        if options.verbose:
            print "creating %s" % datefolder
        os.makedirs(datefolder)
    
    
