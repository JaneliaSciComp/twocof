"""

this is the driver script for the processing; it will run a subscript
many times for each set of files we plan to process

this version will run the scripts sequentially on a single computer; a
different version will submit jobs to the cluster


usage:

process0seq.py tracker [gene@effector [stim@animalno]] [--overwrite] [--branch branchname]

tracker = t1, etc.

tracker, gene@effector, stim@animalno: 
-- if all three present, process that dataset
-- if stim@animalno missing, do all stim in that gene@effector
-- if gene@effector missing, do all in tracker




djo, 1/10

"""


# ------------------------- imports -------------------------
# std lib
import commands
import optparse
import os
import sys
import time



# our stuff
import lib.constants as const
from lib import util


# ------------------------- constants -------------------------
# use a relative path to get the process1.py that is from the
#   same revision as this process0seq.py:
scriptpath = os.path.join(os.path.dirname(__file__), "process1.py")
aggscriptpath = os.path.join(os.path.dirname(__file__), "aggregate-stats.py")


# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    
    # set the file creation mask:
    os.umask(const.umask)
    
    # parse input
    usage = "usage: %prog tracker [trialspec] [options]"
    parser = optparse.OptionParser(usage=usage, version="%%prog %s" % const.__version__)
    parser.add_option("--branch", action="store", dest="branch", default=const.defaultbranch,
        help="branch to process (default %s)" % const.defaultbranch)
    parser.add_option("--overwrite", action="store_true", dest="overwrite", default=True,
        help="overwrite output files if they exist")
    parser.add_option("--nooverwrite", action="store_false", dest="overwrite", 
        help="do not overwrite output files if they exist")
    parser.add_option("--all", action="store_true", dest="processall", default=False, 
        help="process all experiments even if no new data")
    (options, args) = parser.parse_args()
    
    # at a minimum, we need a tracker and branch:
    if len(args) < 1:
        parser.print_usage()
        sys.exit(2)
    tracker = args[0]
    branch = options.branch
    
    # does it exist:
    trackerpath = util.getdirectory("choreography", branch, tracker)
    if not os.path.exists(trackerpath):
        sys.exit("file %s not found" % trackerpath)
    if not os.path.isdir(trackerpath):
        sys.exit("%s does not seem to be a directory" % trackerpath)
    
    
    # as before, the method is going to be to first generate a list of 
    #   gene/stim to be done and then run them either sequentially or
    #   on the cluster
    
    tasklist = []
    
    
    
    if len(args) == 1:
        # do everything in tracker:
        for item in os.listdir(trackerpath):
            geneeffpath = os.path.join(trackerpath, item)
            if os.path.isdir(geneeffpath):
                # assume it's one of our folders (gene@effector)
                for item2 in os.listdir(geneeffpath):
                    stimnopath = os.path.join(geneeffpath, item2)
                    if os.path.isdir(stimnopath):
                        # should be stimulus@animalno:
                        tasklist.append((item, item2))
    elif len(args) > 1:
        geneeff = args[1]
        geneeffpath = os.path.join(trackerpath, geneeff)
        # check that it exists:
        if not os.path.exists(geneeffpath):
            sys.exit("file %s not found" % geneeffpath)
        if not os.path.isdir(geneeffpath):
            sys.exit("%s does not seem to be a directory" % geneeffpath)
        
        if len(args) == 2:
            # no stim@animalno, do everything in that gene@effector:
            for item in os.listdir(geneeffpath):
                stimnopath = os.path.join(geneeffpath, item)
                if os.path.isdir(stimnopath):
                    # should be stimulus@animalno:
                    tasklist.append((geneeff, item))
        else:
            # it's a specific combination:
            stimno = args[2]
            # check that it exists:
            stimnopath = os.path.join(geneeffpath, stimno)
            if not os.path.exists(stimnopath):
                sys.exit("file %s not found" % stimnopath)
            if not os.path.isdir(stimnopath):
                sys.exit("%s does not seem to be a directory" % stimnopath)
            
            tasklist = [(geneeff, stimno)]
            
    
    starttime = time.time()
    print
    print "process0seq.py beginning at %s" % time.asctime()
    print "args: %s" % " ".join(sys.argv[1:])
    
    for geneeff, stimno in tasklist:
        
        # check whether to process or not based on presence 
        #   of new data:
        if (not options.processall and 
            not util.isnewdata(branch, tracker, geneeff, stimno)):
            continue
        
        
        # need to pass on options here; this is clunky, but I
        #   can't see how to get actual options from optparse:
        
        
        if options.overwrite:
            overwriteoption = "--overwrite"
        else:
            overwriteoption = ""
        
        cmdstring = "python %s %s %s %s %s %s" % (scriptpath, branch, tracker, 
            geneeff, stimno, overwriteoption)
        print
        print "executing: %s" % cmdstring
        ans = commands.getoutput(cmdstring)
        print "output: %s" % ans
        print
        # (do logging based on ans)
        
    # if we did the whole tracker, also run the stats aggregation
    #   script
    if len(args) == 1:
        cmdstring = "python %s %s %s" % (aggscriptpath, branch, tracker)
        print "executing: %s" % cmdstring
        ans = commands.getoutput(cmdstring)
        print "output: %s" % ans
        print
    
    
    print "process0seq.py ending at %s" % time.asctime()
    endtime = time.time()
    print "elapsed time: %ss" % (endtime - starttime)
    
    
