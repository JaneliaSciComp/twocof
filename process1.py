#!/bin/env python
"""

this script does the processing on a single "stem" in a single folder;
it is anticipated that it will be called over and over again from another
script; it is designed so that it may be parallelized on the cluster on
a per-stem basis

usage:

process1.py branch tracker gene@effector stim@animalno [--overwrite]

branch, tracker, gene@effector, and stim@animalno are also folders 
in the hierarchy containing the .r files to be worked on


basic plan: given specification of folder:
-- generate list of tasks and output files
-- check if any output files are missing; build list of 
    tasks for missing files (or all files, if overwriting)
-- if no tasks: exit
-- examine tasks; read in data needed for tasks at outset
-- loop through tasks, in order, and perform them



djo, 1/10

"""


# ------------------------- imports -------------------------
import optparse
import os
import shutil
import sys
import time
import traceback


# json not in Python pre-2.5; try the 3rd party lib (mostly for
#   me on my dev box):
try:
    import json
except ImportError:
    try:
        import simplejson as json
    except ImportError:
        raise ImportError("no JSON library found!")

# suppress a warning regarding the matplotlib backend; I set
#   it in two modules, and I don't need to be told that the
#   second time has no effect
# I couldn't get the filters on the specific message to work,
#   so I'm suppressing all UserWarnings, ugh
import warnings
warnings.simplefilter('ignore', UserWarning)
# warnings.filterwarnings('ignore',
#     message=r'.*This call to matplotlib.use() has no effect.*', append=True)


# suppress warnings due to NaN's in numpy:
import numpy
numpy.seterr(invalid='ignore')


# from twocof
import lib.constants as const
from lib import tasks, trackdata, util


# ------------------------- deletederived() -------------------------
def deletederived():
    """
    delete derived .r files
    """
    
    shutil.rmtree(util.getdirectory("derived", branch, tracker, geneeff, stimno))
    
    # end deletederived()

# ------------------------- deleteprocessingstatus() -------------------------
def deleteprocessingstatus():
    """
    delete the processing status file
    """
    
    statusdir = util.getdirectory(const.processingstatuslocation, 
        branch, tracker, geneeff, stimno) 
    statuspath = os.path.join(statusdir, const.processingstatusfilename)    
    
    if os.path.exists(statuspath):
        os.remove(statuspath)
    
    # end deleteprocessingstatus()

# ------------------------- reportscalars() -------------------------
def reportscalars(scalarlist):
    """
    report on list of scalars
    """
    
    print
    print "measurements:"
    for s in sorted(scalarlist):
        print s
    print "%d total measurements" % len(scalarlist)
    
    
    # end reportscalars()

# ------------------------- reportoutputfiles() -------------------------
def reportoutputfiles(tasklistlist):
    """
    report on what files will be output by the tasks
    
    input: list of list of tasks
    output: none (prints to screen)
    """
    
    print
    print "listing output files"
    pathlist = []
    for tasklist in tasklistlist:
        for task in tasklist:
            pathlist.extend(task.outputfilepaths())
    pathlist.sort()
    for path in pathlist:
        print path
    print "%d total output files" % len(pathlist)
    
    # end reportoutputfiles()

# ------------------------- writeprocessingstatus() -------------------------
def writeprocessingstatus(branch, tracker, geneeff, stimno):
    """
    write a file indicating what's been successfully processed
    """
    
    statusdir = util.getdirectory(const.processingstatuslocation, 
        branch, tracker, geneeff, stimno) 
    statuspath = os.path.join(statusdir, const.processingstatusfilename)
    
    
    data = {}
    data["date"] = util.getdatestamp()
    data["status"] = "successful"
    data["dates processed"] = util.getchoreographydatestamps(branch, tracker, geneeff, stimno)
    
    f = open(statuspath, 'wt')
    json.dump(data, f)
    f.close()
    
    # end writeprocessingstatus()

# ------------------------- writesummary() -------------------------
def writesummary():
    """
    write a short summary
    """
    
    print ("%d tasks processed; %d warnings, %d no data, %d errors, %d exceptions" % 
        (sum(len(tasklist) for tasklist in tasklistlist), nwarnings, nnodata, nerrors, nexceptions))
    endtime = time.time()
    print "elapsed time: %ss" % (endtime - starttime)
    # check memory usage (crude):
    print trackdatacache.memcheck()
    
    # this should always be the last line (it will be relied on in
    #   future log parsing)
    print "process1.py ending at %s" % time.asctime()
    
    # end writesummary()

# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    # set the file creation mask:
    os.umask(const.umask)
    
    # parse input
    usage = "usage: %prog branch tracker gene@effector stim@animalno [options]"
    parser = optparse.OptionParser(usage=usage, version="%%prog %s" % const.__version__)
    parser.add_option("--keep-derived", action="store_true", dest="keepderived", default=False,
        help="do not delete derived data at the end of processing")
    parser.add_option("--listoutput", action="store_true", dest="listoutput", default=False,
        help="list all potential output files but do not process")
    parser.add_option("--overwrite", action="store_true", dest="overwrite", default=True,
        help="overwrite output files if they exist")
    parser.add_option("--nooverwrite", action="store_false", dest="overwrite", 
        help="do not overwrite output files if they exist")
    parser.add_option("--raise", action="store_true", dest="raiseexceptions", default=False,
        help="raise exceptions and halt on error (for debugging)")
    parser.add_option("--verbose", action="store_true", dest="verbose", default=False,
        help="increase level of output")
    (options, args) = parser.parse_args()
    
    if len(args) != 4:
        parser.print_usage()
        sys.exit(2)
    
    # input OK on the surface, let's get going
    branch, tracker, geneeff, stimno = args
    
    
    starttime = time.time()
    print
    # this should always be the first line; it'll be relied on
    #   in future log parsing
    print "process1.py beginning at %s" % time.asctime()
    print "args: %s" % " ".join(sys.argv[1:])
    print "twocof version: %s" % const.__version__
    
    
    # if we're on the cluster, log the job ID, etc:
    #if "SGE_O_HOST" in os.environ:
    #print "running on host %s" % os.environ["HOSTNAME"]
    #print "job ID: %s" % os.environ["JOB_ID"]
    #print "job name: %s" % os.environ["JOB_NAME"]
    #print "reserved slots: %s" % os.environ["NSLOTS"]
        
    # if verbose, add lots of version stuff:
    if options.verbose:
        print "Python version %s" % sys.version.replace('\n', ' ')
        print "numpy version %s" % numpy.__version__
        
    
    # look at stem; do all expected stem files exist?  if no, error?
    #   or just warning in logs
    
    
    
    
    # load phase: set up the caching track data thing; we need to give
    #   this to the tasks for later, even though they will not load
    #   data until when (if) the tasks are run
    trackdatacache = trackdata.TrackDataCache()
    
    # generate lists of tasks:
    
    
    # derived file generation tasks; these tasks produce more
    #   .r files that will need to be processed in the main processing
    #   loop
    generationtasklist = [taskclass(branch, tracker, geneeff, stimno, 
            trackdatacache, options.overwrite) 
            for taskclass in tasks.rfilegenerationtasks[tracker]]
    
    # ask all the generation task classes what their scalars are
    derivedscalardata = []
    for taskclass in tasks.rfilegenerationtasks[tracker]:
        derivedscalardata.extend(taskclass.getscalars(tracker))
    
    # r-file processing tasks (primary and derived scalars)
    processingtasklist = [
        taskclass(branch, tracker, geneeff, stimno, scalar, trackdatacache, options.overwrite)
            for scalar in (const.primaryscalardata[tracker] + derivedscalardata)
            for taskclass in tasks.rfileprocessingtasks[tracker]
        ]
    
    # miscellaneous analysis tasks; these tasks produce no data that 
    #   will be further processed, but depend on the availability of 
    #   at least all generated data; by putting them after the main
    #   processing tasks, they can also read plots/stats/etc if desired
    analysistasklist = [taskclass(branch, tracker, geneeff, stimno, 
            trackdatacache, options.overwrite) 
            for taskclass in tasks.analysistasks[tracker]]
    
    # need to go through all tasks, but in defined order, so 
    #   loop over lists of tasks, then over tasks:
    tasklistlist = [generationtasklist, processingtasklist, analysistasklist]
    nwarnings = 0
    nerrors = 0
    nnodata = 0
    nexceptions = 0
    
    
    # at this point we have all the tasks that could be run; this is
    #   the time to report the potential output, if that's what the 
    #   user has specified
    if options.listoutput:
        reportscalars(const.primaryscalardata[tracker] + derivedscalardata)
        reportoutputfiles(tasklistlist)
        writesummary()
        sys.exit(0)
    
    
    # if we reach here, we're actually doing the processing
    
    # start by deleting the processing status file:
    deleteprocessingstatus()
    
    # if we're overwriting files, actually clean out the directories (this is
    #   so old, "obsolete" files aren't left behind for the db loader
    #   to find); obviously, if the dir doesn't exist, it doesn't need to 
    #   be emptied!
    if options.overwrite:
        outputdirset = set()
        for tasklist in tasklistlist:
            for task in tasklist:
                outputdirset.update(os.path.dirname(path) for path in task.outputfilepaths())
        for outputdir in outputdirset:
            # I think it's slightly more efficient to delete the entire dir and
            #   recreate it...the OS may have internal optimizations for removing
            #   whole dir
            if os.path.exists(outputdir):
                shutil.rmtree(outputdir)
                os.mkdir(outputdir)
    
    # now do the work:
    for tasklist in tasklistlist:
        if not options.overwrite:
            #   check for what we have to actually do; if we aren't overwriting
            #   all files, check for existing output and filter out tasks:
            tasklist = [task for task in tasklist if not task.outputexists()]
        for task in tasklist:
            if options.verbose:
                print "working on %s with target %s" % (task.gettaskname(), task.gettasktarget())
            try:
                # in the future, may have some reason to parse
                #   the string and identify errors; for now
                #   just print it
                status = task.performtask()
                if status.startswith("error"):
                    nerrors += 1
                    # could enable this later, but for now, don't
                    #   output more info on errors (which I generate)
                    #   only on exceptions (which I don't)
                    # status += ("\n  in task %s\n  with target %s" %
                    #     (task.gettaskname(), task.gettasktarget()))
                elif status.startswith("warning"):
                    nwarnings += 1
            except util.NoDataError, e:
                nnodata += 1
                status = "error: no data error"
                status += ("\n  in task %s\n\twith target %s" %
                    (task.gettaskname(), task.gettasktarget()))
                status += "\n  traceback: %s" % traceback.format_exc()
                
                if options.raiseexceptions:
                    # for debugging:
                    raise
            
            except Exception, e:
                nexceptions += 1
                status = "error: caught exception"
                status += ("\n  in task %s\n\twith target %s" %
                    (task.gettaskname(), task.gettasktarget()))
                status += "\n  traceback: %s" % traceback.format_exc()
                
                if options.raiseexceptions:
                    # for debugging:
                    raise
            
            if not options.verbose:
                # process the output a little; in particular, we write so many
                #   files that one line per written file is overwhelming; so
                #   prune them out (for now) (later may allow a command-line
                #   switch):
                statuslines = status.split('\n')
                # suppress lines that report writing a file:
                statuslines = [line for line in statuslines if not line.startswith("wrote")]
                # suppress lines that report *not* writing a file:
                statuslines = [line for line in statuslines if "skipping" not in line]
                status = '\n'.join(statuslines)
            
            if status:
                if options.verbose:
                    print time.asctime()
                    print "\t%s" % status
                    
                    # try to monitor memory, too:
                    meminfo = util.getmeminfo()
                    print "\tfree memory: %s M" % (meminfo["free"] // 1000) 
                    
                    # not sure if this will help with the buffering issue or not:
                    sys.stdout.flush()
                else:
                    print status
    
    # clean up section
    
    # delete all the derived .r files (save space, don't really need
    #   them once stats etc. are produced:
    if not options.keepderived:
        deletederived()
    
    
    writesummary()
    
    # only write this file if we have no exceptions!
    if nexceptions == 0:
        writeprocessingstatus(branch, tracker, geneeff, stimno)
    
    
