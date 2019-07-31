#!/bin/env python
"""

this script reads and parses a summary file written by process0cluster.py

it produces a variety of information regarding how long each job ran, etc

works for v1.2.1 and later (relies on that logging format)



djo, 3/11


"""


# ------------------------- imports -------------------------
# std lib
import os
import sys
import time


# twocof
from lib import constants as const
from lib import util


# ------------------------- constants -------------------------
# for strptime; it's derived from: Mon Apr 25 22:03:14 2011
timestampformat = "%a %b %d %H:%M:%S %Y"



# ------------------------- findglobalspan() -------------------------
def findglobalspan(loglines):
    """
    given lines from a summary file, find the overall span of time
    for which the jobs ran; ie, what's the length of time between
    the start of the first job and the end of the last?
    
    input: list of lines
    output: (first, last) times
    """
    
    times = []
    for line in loglines:
        if line.startswith("process1.py"):
            line = line.strip()
            i = line.find("ing at") + 7
            temptime = time.strptime(line[i:], timestampformat)
            times.append(temptime)
    
    return min(times), max(times)
    
    # end findglobalspan()

# ------------------------- parsesummary() -------------------------
def parsesummary(loglines):
    """
    parse lines from a summary file
    
    lines look like:
    
    args: rd t1-don-test kolodkin1m@n b_100Hz500mV_20s1x30s0s#n#n#n@1 --overwrite
    twocof version: 1.2.2.x
    719 tasks processed; 0 warnings, 0 errors, 0 exceptions
    elapsed time: 433.622341871s
    memory: VmSize:	  253048 kB	VmRSS:	   80932 kB	
    
    
    input: list of lines
    output: dict: {"screen tracker geneeff stimno": (elapsed time (s), 
        VmSize (kB), VmRSS (kB))}
    """
    
    data = {}
    
    for stanza in stanzaiter(loglines):
        for line in stanza:
            if line.startswith("args"):
                # grab all the args as our identifier
                name = line[6:]
            if line.startswith("elapsed"):
                items = line.split()
                time = float(items[2][:-1])
            if line.startswith("memory"):
                items = line.split()
                memory1 = int(items[2])
                memory2 = int(items[5])
        data[name] = time, memory1, memory2
    
    return data
    
    # end parsesummary()

# ------------------------- stanzaiter() -------------------------
def stanzaiter(lines):
    """
    iterator over "stanzas" in log file
    
    NOTE: will fail on logs that have errors from within process1.py!
        (in those cases, process1.py won't log its own errors like it
        logs task errors, so the log won't be parsed correctly)
    
    input: list of lines
    output: iterator returning list of lines between "process1.py beginning/ending"
    """
    
    lineiter = iter(lines)
    
    line = lineiter.next()
    stanzalines = []
    while not line.startswith("summary information"):
        if line.startswith("process1.py beginning"):
            # new stanza:
            stanzalines = []
        elif line.startswith("process1.py ending"):
            # done with this stanza, output:
            yield stanzalines
        else:
            # continue existing stanza:
            stanzalines.append(line)
        
        line = lineiter.next()
    
    # end stanzaiter()

# ------------------------- datantiles() -------------------------
def datantiles(datadict, n, column):
    """
    quick n-tile finder
    
    input: data and how many divisions; which column
    output: n + 1 division borders
    """
    
    data = getsorteddata(datadict, column)
    ntiles = [data[i * len(data) // n] for i in range(n)]
    ntiles.append(data[-1])
    
    return ntiles
    
    # end datantiles()

# ------------------------- getsorteddata() -------------------------
def getsorteddata(datadict, column):
    """
    get data from dict in column and sort it
    """
    
    return sorted([item[column] for item in datadict.values()])
    
    # end getsorteddata()

# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    
    filename = sys.argv[1]
    
    if not os.path.exists(filename):
        sys.exit("can't find input file %s" % filename)
    loglines = open(filename).readlines()
    
    data = parsesummary(loglines)
    
    n = 4
    print
    print "time %d-tiles" % n
    print datantiles(data, n, 0)
    
    print
    print "VmSize %d-tiles" % n
    print datantiles(data, n, 1)
    
    print
    print "VmRSS %d-tiles" % n
    print datantiles(data, n, 2)
    
    print
    
    # how many runs longer than twice the median?
    timelist = getsorteddata(data, 0)
    nexpt = len(timelist)
    median = timelist[nexpt // 2]
    print "median time = %s" % median
    
    print "5 longest =", timelist[-5:]
    
    for mult in [2, 3, 5, 10]:
        morethan = len([t for t in timelist if t > mult * median])
        print "%d times of %d longer than %d x median" % (morethan, nexpt, mult)
    
    
    mint, maxt = findglobalspan(loglines)
    print
    print "overall start time: %s" % time.strftime(timestampformat, mint)
    print "overall end time: %s" % time.strftime(timestampformat, maxt)
