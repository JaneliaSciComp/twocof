"""

this is a script for aggregating some statistics; it's
a per-tracker task, and that doesn't fit into our pipeline
right now; however, we need the data, so I'll run it by hand
in the short term and integrate it later


given a tracker:
- read in all stats files
- given a set of variables
- for each time window, protocol, and variable, aggregate all stats, one row
    per gene eff
- protocol = stimulus/animal number piece


output looks like:

statistics/
    tracker/
        aggregated-stats/
            variable1/
                stimno/
                    timewindow1
                    timewindow2
            variable2/
                stimno/
                    timewindow1
                    timewindow2


djo, 11/10

"""


# ------------------------- imports -------------------------
import os
import shutil
import sys


# from twocof:
from lib import constants as const
from lib import util


# ------------------------- constants -------------------------
# output file: this is a bit awkward since t1, t2 both have
#   decimal points in them:
# tracker.variable.t1-t2.stats.txt
filenamepattern = "%s.%s.%s-%s.stats.txt"


# ------------------------- createdirectories() -------------------------
def createdirectories(branch, tracker, clean=True):
    """
    create the upper level of output directories; 
    
    input: tracker; clean = True means delete existing first
    """
    
    trackerstatspath = util.getdirectory("statistics", branch, tracker)
    aggregatedstatspath = os.path.join(trackerstatspath, const.aggregatedstatsdir)
    
    if os.path.exists(aggregatedstatspath):
        if clean:
            # delete all contents:
            shutil.rmtree(aggregatedstatspath)
            os.mkdir(aggregatedstatspath)
    else:
        os.mkdir(aggregatedstatspath)
    
    
    for scalar in const.aggregatedstatsscalarlist[tracker]:
        scalarpath = os.path.join(aggregatedstatspath, scalar)
        if not os.path.exists(scalarpath):
            os.mkdir(scalarpath)
    
    # end createdirectories()

# ------------------------- readalldata() -------------------------
def readalldata(branch, targettracker):
    """
    read all data from a tracker
    
    input: branch, tracker
    output: dict with rawdata[scalar][expt][interval] = stats
    """
    
    rawdata = {}
    for scalar in const.aggregatedstatsscalarlist[targettracker]:
        rawdata[scalar] = {}
        for tracker, geneeff, stimno in util.walktracker([branch, targettracker], "statistics"):
            filedir = util.getdirectory("statistics", branch, tracker, geneeff, stimno)
            filename = "%s@%s@%s.%s-stats-perobj.txt" % (geneeff, tracker, stimno, scalar)
            filepath = os.path.join(filedir, filename)
            
            if os.path.exists(filepath):
                # eventually put in try-except with reporting here?
                rawdata[scalar][(geneeff, stimno)] = readstatsfile(filepath)
    return rawdata
    
    # end readalldata()

# ------------------------- readstatsfile() -------------------------
def readstatsfile(filepath):
    """
    read a stats file and return data
    
    input: path to stats file
    output: dict {(t1, t2): [ave, stddev, N, sqsum, sumsqares]
    """
    
    data = open(filepath).readlines()
    
    results = {}
    # skip header line:
    for line in data[1:]:
        line = line.strip()
        items = [float(x) for x in line.split()]
        # to normalize the time intervals, use string formatting to
        #   truncate (not round) to tenths of a second:
        t1 = "%.1f" % items[0]
        t2 = "%.1f" % items[1]
        results[t1, t2] = items[2:]
    
    return results
    
    # end readstatsfile()

# ------------------------- transformdata() -------------------------
def transformdata(rawdata, tracker):
    """
    here, interval = (t1, t2)
    
    input: dict with rawdata[scalar][(geneeff, stimno)][interval] = stats
           tracker
    output: dict with sliceddata[scalar][stimno][interval][geneeff] = stats
    """
    
    sliceddata = {}
    for scalar in const.aggregatedstatsscalarlist[tracker]:
        sliceddata[scalar] = {}
        for geneeff, stimno in rawdata[scalar]:
            if stimno not in sliceddata[scalar]:
                sliceddata[scalar][stimno] = {}
            for interval in rawdata[scalar][geneeff, stimno]:
                if interval not in sliceddata[scalar][stimno]:
                    sliceddata[scalar][stimno][interval] = {}
                sliceddata[scalar][stimno][interval][geneeff] = rawdata[scalar][geneeff, stimno][interval]
    return sliceddata
    
    # end transformdata()

# ------------------------- writeaggregatedstatsfile() -------------------------
def writeaggregatedstatsfile(filepath, data):
    """
    write out the file
    
    input:  filepath, dict: {geneeff: stats list}
            stats list = ave, stddev, N, sq sum, sum sqrs
    output: none (writes file)
    """
    
    f = open(filepath, 'wt')
    
    # header lines
    f.write("# %s\n" % filepath)
    f.write("gene\tave\tstddev\tsquaredsum\tsumsquares\tt-value\tp-value\n")
    
    for geneeff in sorted(data.keys()):
        genestring = "%s" % geneeff
        linedata = tuple([genestring] + data[geneeff])
        f.write("%s\t%f\t%f\t%d\t%f\t%f\t\t\n" % linedata)
    
    f.close()
    
    # end writeaggregatedstatsfile()

# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    # parse input
    if len(sys.argv) < 3:
        sys.exit("usage: python aggregate-stats.py branch tracker")
    branch = sys.argv[1]
    tracker = sys.argv[2]
    if tracker not in const.stim2col:
        sys.exit("unknown tracker %s" % tracker)
    if branch not in const.branches:
        sys.exit("unknown branch %s" % branch)
    
    # check that anything needs to be done:
    if not const.aggregatedstatsscalarlist[tracker]:
        sys.exit("no scalars to aggregate")
    
    # set the file permissions
    os.umask(const.umask)
    
    # read all the data in:
    rawdata = readalldata(branch, tracker)
    
    # re-organize:
    sliceddata = transformdata(rawdata, tracker)
    
    
    # set up directories
    createdirectories(branch, tracker)
    
    # iterate over scalars and time intervals:
    for scalar in sliceddata:
        for stimno in sliceddata[scalar]:
            for t1, t2 in sliceddata[scalar][stimno]:
                filename = filenamepattern % (tracker, scalar, t1, t2)
                dirpath = util.getdirectory("statistics", branch,
                    tracker, const.aggregatedstatsdir)
                dirpath = os.path.join(dirpath, scalar, stimno)
                if not os.path.exists(dirpath):
                    os.mkdir(dirpath)
                filepath = os.path.join(dirpath, filename)
                writeaggregatedstatsfile(filepath, sliceddata[scalar][stimno][(t1, t2)])
    
    


