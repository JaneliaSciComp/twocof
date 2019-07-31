"""

this is the driver script for the processing; it will run a subscript
many times for each set of files we plan to process

this version must be run on the cluster head node; it will
submit jobs to the cluster


usage:

process0cluster.py tracker [gene@effector [stim@animalno]] [other options] 

tracker = t1, etc.

tracker, gene@effector, stim@animalno: 
-- if all three present, process that dataset
-- if stim@animalno missing, do all stim in that gene@effector
-- if gene@effector missing, do all in tracker

other options:

--branch branchname: specify branch of pipeline (default = ??)

--email address: address to email at end of summary job (not clear if 
    this is working on our cluster?)

--overwrite: overwrite existing files

--short: submit to "short" queue on cluster (for short jobs, testing)

note on short queue: in May 2016, short queue was changed; now
you specify a time < 1 hour and your job will be killed if it
exceeds that time; the jobs I run in short have always been
< 10 minutes, so I set it to 30 min to be safe





(derived from process0seq.py)

djo, 3/10

"""


# ------------------------- imports -------------------------
# std lib
import commands
import optparse
import os
import sys
import time



# mine
import lib.constants as const
from lib import util


# ------------------------- constants -------------------------
# executables; using our own location, create absolute paths for the
#   Python scripts we submit:
mydir = os.path.dirname(__file__)
wrapperpath = os.path.abspath(os.path.join(mydir, "runsce.sh"))
process1path = os.path.abspath(os.path.join(mydir, "process1.py"))
aggregatepath = os.path.abspath(os.path.join(mydir, "aggregate-stats.py"))
summarize0path = os.path.abspath(os.path.join(mydir, "summarize0cluster.py"))

# output filenames:
summaryoutputfilename = "summarize0.out"
aggstatsoutputfilename = "aggregate-stats.out"


# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    
    # check that we're on the cluster; see if SGE_ROOT is defined:
    #if "SGE_ROOT" not in os.environ:
    #    print "this script must be run on the cluster head node; I don't see the SGE_ROOT variable, so you're probably not there"
    #    sys.exit(0)
    
    
    # set the file creation mask:
    os.umask(const.umask)
    
    # note: building the task list ought to be refactored out of the 
    #   two process0 script; or, alternately, make the cluster/single machine
    #   choice a command-line option in a unified process0.py?
    
    
    
    # parse input
    usage = "usage: %prog tracker [trialspec] [options]"
    parser = optparse.OptionParser(usage=usage, version="%%prog %s" % const.__version__)
    parser.add_option("--branch", action="store", dest="branch", default=const.defaultbranch,
        help="branch to process (default %s)" % const.defaultbranch)
    parser.add_option("--charge", action="store_true", dest="charge", default=False,
        help="charge for processing (default off)")
    parser.add_option("--nocharge", action="store_false", dest="charge",
        help="do not charge for processing")
    parser.add_option("--email", action="store", dest="email",
        help="(this feature is disabled)")
    parser.add_option("--overwrite", action="store_true", dest="overwrite", default=True,
        help="overwrite output files if they exist")
    parser.add_option("--nooverwrite", action="store_false", dest="overwrite", 
        help="do not overwrite output files if they exist")
    parser.add_option("--all", action="store_true", dest="processall", default=False, 
        help="process all experiments even if no new data")
    parser.add_option("--nosubmit", action="store_true", dest="nosubmit", default=False,
        help="don't submit jobs; print submit commands instead (for testing)")
    parser.add_option("--short", action="store_true", dest="short", default=False,
        help="submit to the 'short' cluster queue (for testing)")
    parser.add_option("--verbose", action="store_true", dest="verbose", default=False,
        help="increase level of output")
    (options, args) = parser.parse_args()
    
    # at a minimum, we need a tracker and branch:
    if len(args) < 1:
        parser.print_usage()
        sys.exit(2)
    branch = options.branch
    tracker = args[0]
    
    # does it exist:
    trackerpath = util.getdirectory("choreography", branch, tracker)
    
    if not os.path.exists(trackerpath):
        print "file %s not found" % trackerpath
        sys.exit(1)
    if not os.path.isdir(trackerpath):
        print "%s does not seem to be a directory" % trackerpath
        sys.exit(1)
    
    
    
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
            print "file %s not found" % geneeffpath
            sys.exit(1)
        if not os.path.isdir(geneeffpath):
            print "%s does not seem to be a directory" % geneeffpath
            sys.exit(1)
        
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
                print "file %s not found" % stimnopath
                sys.exit(1)
            if not os.path.isdir(stimnopath):
                print "%s does not seem to be a directory" % stimnopath
                sys.exit(1)
            
            tasklist = [(geneeff, stimno)]
            
    
    
    starttime = time.time()
    # the date/time stamp on all output is the start of the run:
    print
    print "process0cluster.py beginning at %s" % time.asctime()
    print "args: %s" % " ".join(sys.argv[1:])
    
    
    # logging:
    # logfiles will go in a central folder, nested dates, based on start date:
    # eg, logs/2011/01/23/....
    logfiledir = util.getdirectory("logs", branch, tracker)
    
    datetimestamp = time.strftime("%Y%m%d%H%M")
    yearstring = datetimestamp[:4]
    monthstring = datetimestamp[4:6]
    datestring = datetimestamp[6:8]
    trackerdatetime = "%s-%s-%s" % (branch, tracker, datetimestamp)
    datedlogfiledir = os.path.join(logfiledir, "twocof", yearstring, monthstring, datestring)
    templogfiledir = os.path.join(datedlogfiledir, "%s-%s-%s" % (branch, tracker, datetimestamp))
    if not os.path.exists(templogfiledir) and not options.nosubmit:
        # recursively create all needed directories:
        os.makedirs(templogfiledir)
    
    # continuing effort to fix TWOCOF-33...specify a cluster-specific
    #   matplotlib config dir so it doesn't conflict with my own
    #   personal one (which is in ~/.matplotlib); there are font cache
    #   issues between the two environments, and we can't let them
    #   interfere with each other; this is only needed when we run on
    #   the cluster, so it's in this file rather than, eg, process1.py
    # the "-V" in the command template means "pass along any env vars",
    #   so we need merely set the variable in this process:
    os.environ["MPLCONFIGDIR"] = const.matplotlibconfigdir 
    
    
    # note: -pe batch # --> reserves number of slots per node
    #   (there are 8 available on current cluster; using # = 8 reserves 
    #   whole node)
    #cmdtemplate = 'qsub -N %(process1jobname)s -j y -pe batch %(slots)s -o %(outputpath)s -b y -V %(sgeoptions)s %(wrapper)s %(process1)s %(process1args)s'
    cmdtemplate = 'bsub -J %(process1jobname)s -n %(slots)s -o %(outputpath)s %(lsfoptions)s %(wrapper)s %(process1)s %(process1args)s'
    #print cmdtemplate
    
    # dict for holding arguments to qsub and the processing job
    jobargs = {}
    jobargs["process1jobname"] = "%s-%s-process1.py" % (branch, tracker)
    jobargs["wrapper"] = wrapperpath
    jobargs["process1"] = process1path
    
    njobs = 0
    for geneeff, stimno in tasklist:
        
        # check whether to process or not based on presence 
        #   of new data:
        if (not options.processall and 
            not util.isnewdata(branch, tracker, geneeff, stimno)):
            #print "Skipping.."
            continue
        
        outputfilename = "%s-%s-%s.txt" % (tracker, geneeff, stimno)
        print "output: " + outputfilename
        jobargs["outputpath"] = os.path.join(templogfiledir, outputfilename)
        
        # prepare the arguments; need to escape the '#' and then
        #   double-quote it for filenames to arrive correctly
        #   at process1.py:
        arguments= [branch, tracker, geneeff, stimno]
        jobargs["process1args"] = ['"%s"' % item.replace('#', '#') for item in arguments]
        
        # optional arguments (to process1):
        if options.overwrite:
            jobargs["process1args"].append('--overwrite')
        if options.verbose:
            jobargs["process1args"].append("--verbose")
        
        # concatenate the process1 args into a string:
        jobargs["process1args"] = " ".join(jobargs["process1args"])
        
        # optional options (to SGE):
        jobargs["sgeoptions"] = ""
        jobargs["lsfoptions"] = ""
        if options.short:
            # if you're using the short queue, don't reserver more
            #   than one slot:
            # sge arg
            #jobargs["sgeoptions"] += ' -l h_rt=1800' 
            # lsf arg
            jobargs["lsfoptions"] += ' -W 30'
            jobargs["slots"] = 1
        else:
            # choose number of slots based on data size:
            if util.islargedata(branch, tracker, geneeff, stimno):
                jobargs["slots"] = const.slotsperlargedataset
            else:
                jobargs["slots"] = const.slotspersmalldataset
        if options.charge:
            # charge the appropritae account for the computer time:
            #jobargs["sgeoptions"] += ' -A %s' % const.accountstring
            # lsf billing arg
            jobargs["lsfoptions"] += ' -P %s' % const.accountstring
        
        cmdstring = cmdtemplate % jobargs
        print "Here" + cmdstring
        if options.nosubmit:
            print cmdstring
        else:
            ans = commands.getoutput(cmdstring)
        njobs += 1
    
    
    # the next set of jobs are only submitted if any new data has
    #   been processed; we use njobs to tell us if that's the case:
    
    if njobs == 0:
        # remove the log file dir that will forever otherwise be empty:
        if os.path.isdir(templogfiledir) and len(os.listdir(templogfiledir)) == 0:
            os.rmdir(templogfiledir)
        print "all processing up-to-date!"
        
    else:
        # if we did the whole tracker, also submit the stats aggregator;
        #   like summary, hold until rest is done; also, use the short queue
        if len(args) == 1:
            jobargs["outputpath"] = os.path.join(datedlogfiledir, "%s-%s" % (trackerdatetime, aggstatsoutputfilename))
            jobargs["aggjobname"] = "%s-%s-aggregate-stats.py" % (branch, tracker)
            jobargs["aggregatepath"] = aggregatepath
            #summarytemplate = 'qsub -N %(aggjobname)s -j y -o %(outputpath)s -l h_rt=1800 -b y -hold_jid %(process1jobname)s -V %(options)s %(wrapper)s %(aggregatepath)s %(branch)s %(tracker)s'
            summarytemplate = 'bsub -J %(aggjobname)s -o %(outputpath)s -W 30 -w \'ended(\"%(process1jobname)s\")\' %(options)s %(wrapper)s %(aggregatepath)s %(branch)s %(tracker)s'
            print summarytemplate
            # add more options as desired
            jobargs["options"] = ""
            if options.charge:
                # charge the appropritae account for the computer time:
                jobargs["options"] += ' -P %s' % const.accountstring
            
            jobargs["branch"] = branch
            jobargs["tracker"] = tracker
            
            cmdstring = summarytemplate % jobargs
            print cmdstring
            if options.nosubmit:
                print cmdstring
            else:
                ans = commands.getoutput(cmdstring)
            njobs += 1
        
        
        # submit a summary generating job here that will
        #   only run after the others have run (and they have a known job name)
        # use the short queue
        jobargs["outputpath"] = os.path.join(datedlogfiledir, "%s-%s" % (trackerdatetime, summaryoutputfilename))
        jobargs["summaryjobname"] = "%s-%s-summarize0cluster.py" % (branch, tracker)
        jobargs["summarize0path"] = summarize0path
        jobargs["templogfiledir"] = templogfiledir
        #summarytemplate = 'qsub -N %(summaryjobname)s -j y -o %(outputpath)s -l h_rt=1800 -b y -hold_jid %(process1jobname)s -V %(options)s %(wrapper)s %(summarize0path)s %(templogfiledir)s'
        summarytemplate = 'bsub -J %(summaryjobname)s -o %(outputpath)s -W 30 -w \'ended(\"%(process1jobname)s\")\' %(options)s %(wrapper)s %(summarize0path)s %(templogfiledir)s'
        #print summarytemplate
        
        # add more options as desired
        jobargs["options"] = ""
        if options.email:
            # email on abort or end
            # NOTE: this option is not available on our cluster at this time!
            pass
            # optionaloptions += " -m ae -M %s" % options.email
        if options.charge:
            # charge the appropritae account for the computer time:
            #jobargs["options"] += ' -A %s' % const.accountstring
            jobargs["options"] += ' -P %s' % const.accountstring
        
        cmdstring = summarytemplate % jobargs
        print cmdstring
        if options.nosubmit:
            print cmdstring
        else:
            ans = commands.getoutput(cmdstring)
        njobs += 1
    
    # # of jobs: one per tasks + aggregate stats + summary:
    if options.nosubmit:
        print "%d jobs not submitted" % (njobs)
    else:
        print "%d jobs submitted" % (njobs)
    print "process0cluster.py ending at %s" % time.asctime()
    endtime = time.time()
    print "elapsed time: %ss" % (endtime - starttime)
    
