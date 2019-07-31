#!/usr/bin/python
"""

look through output files and present a summary of what happened;
run by "process0cluster.py" after all actual processing jobs are 
done

runs under python 2.4 so I can submit to cluster without a
shell wrapper to set up a different python environment


djo, 4/10

"""



# ------------------------- imports -------------------------
import getpass
import os
import shutil
import sys


# ------------------------- constants -------------------------
# name of aggregate file ("parallel" to input dated folder)
# test name:
summaryfilenametemplate = "%s-summary.txt"



# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print "usage: summarize0cluster.py datedoutputdir"
        sys.exit(2)
    
    datedoutputpath = os.path.abspath(sys.argv[1])
    
    outputdir = os.path.dirname(datedoutputpath)
    datedoutputdirname = os.path.basename(datedoutputpath)
    
    # open up output file
    summaryfilepath = os.path.join(outputdir, summaryfilenametemplate % datedoutputdirname)
    outputfilelist = os.listdir(datedoutputpath)
    # filter out non-.txt, .DS_Store (etc)
    outputfilelist = sorted([fn for fn in outputfilelist 
        if fn[0] != '.' and fn.endswith('.txt')])
    f = open(summaryfilepath, 'wt')
    f.write("\nsummary file for %s\n" % datedoutputdirname)
    
    f.write("\nrun by username %s\n" % getpass.getuser())
    
    stored = {}
    for fn in outputfilelist:
        filepath = os.path.join(datedoutputpath, fn)
        data = open(filepath, 'rt').readlines()
        
        # write whole thing to aggregate file
        f.write("\n")
        f.writelines(data)
        
        # now look for specific info we want
        # from the end, look for "elapsed", "tasks processed", and "args"
        
        # also: add alternate output if usual info can't be found; this
        #   will catch, eg, if the SGE system kills the job
        
        timefound = False
        errorsfound = False
        argsfound = False
        for line in reversed(data):
            if not argsfound and line.startswith("args:"):
                argsline = line
                argsfound = True
            if not timefound and "elapsed time" in line:
                # line looks like: elapsed time: 251.506650925s
                elapsedtime = float(line.split()[-1][:-1])
                timefound = True
            if not errorsfound and "tasks processed" in line:
                # line looks like: 153 tasks processed; 0 warnings, 0 no data, 0 errors, 0 exceptions
                temp = line.split()
                tasks = int(temp[0])
                warnings = int(temp[3])
                nodata = int(temp[5])
                errors = int(temp[8])
                exceptions = int(temp[10])
                errorsfound = True
            
            if timefound and errorsfound and argsfound:
                break
            
        if timefound and errorsfound and argsfound:
            stored[fn] = [argsline, tasks, warnings, nodata, errors, exceptions, elapsedtime]
        else:
            # couldn't find all data; consider it an exception and output the
            #   filename, which will specify the experiment:
            stored[fn] = ["%s\n" % fn, 1, 0, 0, 0, 1, 0.0]
        
    
    # write out the summary data
    f.write("\nsummary information\n")
    f.write("experiments with errors or exceptions:\n\n")
    fnlist = sorted(stored.keys())
    totaltime = 0.0
    totalwarnings = 0
    totalnodata = 0
    totalerrors = 0
    totalexceptions = 0
    totaltasks = 0
    maxtime = 0
    for fn in fnlist:
        argsline, tasks, warnings, nodata, errors, exceptions, elapsedtime = stored[fn]
        if errors or exceptions or nodata:
            f.write("%s" % argsline)
            f.write("\t%d warnings, %d no data, %d errors, %d exceptions\n" % (warnings, nodata, errors, exceptions))
        totaltasks += tasks
        totalwarnings += warnings
        totalnodata += nodata
        totalerrors += errors
        totalexceptions += exceptions
        totaltime += elapsedtime
        maxtime = max(maxtime, elapsedtime)
    f.write("\n%d total datasets\n" % len(outputfilelist))
    f.write("%s total tasks; %d total warnings, %d total no data, %d total errors; %d total exceptions\n" %
        (totaltasks, totalwarnings, totalnodata, totalerrors, totalexceptions))
    f.write("longest run: %f (%.1f hours)\n" % (maxtime, maxtime / 3600.))
    f.write("total elapsed time: %f (%.1f hours)\n\n" % (totaltime, totaltime / 3600.))
    
    
    f.close()
    
    
    # clean up: delete the files!
    shutil.rmtree(datedoutputpath)
    # print "pretending to remove %s" % datedoutputpath
