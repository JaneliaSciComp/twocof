"""

this script runs periodically via Hudson and checks whether any twocof
jobs have finished recently; if so, it signals Hudson via its exit status,
and console output will contain some minimal summary information; Hudson
will be set up to include that output in notification emails that it will
send upon the appropriate exist status


initial details:
- run daily, early am?
- use exit = success to signal "recent job completion"
- recent = previous 24 hrs (or: input number of hours/days to look back?)
- detected by looking for new summary files in log area
- only looks back some fixed amount (current and previous month)
- summary should include minimally: 
-- which branch, tracker
-- number of errors/etc. (for devs)
- summary could include:
-- list of expts processed
-- total/longest CPU time
-- list of expts with errors (for devs)



other notes:
- annoyance: logs are stored by date started, so you'll have to
    go through a bunch of folders to see what's finished recently;
    probably need to set a limit to how far back to go? 
    current and previous month is enough?
- creation of log dirs is done inline (no function)
- filename spec for summary file is hard-coded in summarize script, ugh;
    move it into constants



djo, 12/11

"""

# ------------------------- imports -------------------------
# stdlib
import datetime
import os
import time
import sys

# twocof
import lib.constants as const
from lib import util


# ------------------------- constants -------------------------

# number of days to look back, if not specified by the user:
defaultndays = 7


# ------------------------- filtersummaries() -------------------------
def filtersummaries(summarylist, ndays):
    """
    filter a list of summary files; return subset that were created/modified
    in the last ndays (these files are created all at once and not modified)
    
    input: list of paths; number of days
    output: list of paths
    """
    
    now = datetime.datetime.now()
    delta = datetime.timedelta(days=ndays)
    results = []
    for summarypath in summarylist:
        mtime = os.path.getmtime(summarypath)
        mtime = datetime.datetime.fromtimestamp(mtime)
        if mtime > now - delta:
            results.append(summarypath)
    
    return results
    
    # end filtersummaries()

# ------------------------- formatsummaryinfo() -------------------------
def formatsummaryinfo(info):
    """
    input: summary info struct
    output: string to pring
    """
    
    lines = []
    
    lines.append("Tracker %s/%s was run\n" % (info.branch, info.tracker))
    
    lines.extend(info.lastfivelines)
    
    result = "".join(lines)
    
    return result
    
    # end formatsummaryinfo()

# ------------------------- getsummarypaths() -------------------------
def getsummarypaths(monthdir):
    """
    for a given month dir, find all the summaries for all dates
    
    input: path to month dir
    output: list of paths to summary files in that month dir
    """
    
    results = []
    
    # try to ignore files like .DS_Store throughout
    
    if not os.path.exists(monthdir):
        # the month dir for this month hasn't been created yet
        return results

    datelist = [fn for fn in os.listdir(monthdir) if not fn.startswith('.')]
    for date in datelist:
        datepath = os.path.join(monthdir, date)
        filenames = os.listdir(datepath)
        for fn in filenames:
            if not fn.startswith('.') and fn.endswith("-summary.txt"):
                results.append(os.path.join(datepath, fn))
    
    return results
    
    # end getsummarypaths()


# ------------------------- geteligiblemonths() -------------------------
def geteligiblemonths():
    """
    input: none
    output: list of (month, year) to look at
    """
    
    now = time.localtime()
    
    curryear = now.tm_year
    currmonth = now.tm_mon
    
    if currmonth == 1:
        # previous month is previous year:
        prevyear = curryear - 1
        prevmonth = 12
    else:
        prevyear = curryear
        prevmonth = currmonth - 1
        
    return [(prevmonth, prevyear), (currmonth, curryear)]
    
    # end geteligiblemonths()


# ------------------------- logpathfrommonth() -------------------------
def logpathfrommonth(branch, tracker, month, year):
    """
    input: branch, tracker, month, year
    output: path to log file dir for given year and month
    """
    
    return os.path.join(util.getdirectory("logs", branch, tracker), 
        "twocof", "%4d" % year, "%02d" % month)
    
    # end logpathfrommonth()

# ------------------------- parsesummaryfile() -------------------------
def parsesummaryfile(summaryfilepath):
    """
    input: path to a summary file
    output: struct containing info
    """
    
    # desired info:
    #   minimal: branch, tracker, date finished 
    #   highly desired: ntrials, nerrors
    #   also: list of expts processed, list of expts with error/exception,
    #       total CPU, longest CPU
    
    
    info = util.struct(path=summaryfilepath)
    
    # extract info from the filename and the file contents:
    
    info = parsesummaryfilename(info)
    
    info = parsesummaryfilecontents(info)
    
    
    return info
    
    # end parsesummaryfile()

# ------------------------- parsesummaryfilecontents() -------------------------
def parsesummaryfilecontents(info):
    """
    parse info from contents of summary file
    
    input: info struct
    output: info struct (with even more info)
    """
    
    # for a start, just take the last five lines of the 
    #   summary file (includes a final empty line)
    
    data = open(info.path, 'rt').readlines()
    
    info.lastfivelines = data[-5:]
    
    return info
    
    # end parsesummaryfilecontents()

# ------------------------- parsesummaryfilename() -------------------------
def parsesummaryfilename(info):
    """
    parse info out of the summary file filename
    
    input: info struct
    output: info struct (with more info)
    """
    
    fn = os.path.basename(info.path)
    info.filename = fn
    
    # I want to grab branch and tracker out of the filename; form is:
    #   rd-t1-don-test-201103091350-summary.txt
    # unfortunately, I allow '-' in the tracker names...
    
    if fn.count("-") < 3:
        # improperly formed name, should never happen:
        raise ValueError("bad filename: %s" % fn)
    
    # grab the branch from the left
    branch, sep, rest = fn.partition("-")
    
    # peel off the suffix and date from the right:
    temp = rest.replace("-summary.txt", "")
    tracker, sep, datestamp = temp.rpartition("-")
    
    info.branch = branch
    info.tracker = tracker
    
    return info
    
    # end parsesummaryfilename()


# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    if len(sys.argv) > 1:
        ndays = int(sys.argv[1])
    else:
        ndays = defaultndays
    
    
    
    
    # gather list of (potential) branches/trackers, and months/years;
    #   assemble list of (branch, tracker, month, year) to examine
    brtrlist = [(branch, tracker) for branch in const.branches for tracker in const.trackers]
    moylist = geteligiblemonths()
    
    btmylist = [(branch, tracker, month, year) for branch, tracker in brtrlist for month, year in moylist]
    
    
    # loop over list: for each combination, get unique list of 
    #   month/year log file dirs
    monthdirs = set(logpathfrommonth(branch, tracker, month, year) 
        for branch, tracker, month, year in btmylist)
    
    
    # get summary files:
    allsummarylist = []
    for monthdir in monthdirs:
        allsummarylist.extend(getsummarypaths(monthdir))
    
    
    # filter summary files to those created/modified in desired time window:
    newsummarylist = filtersummaries(allsummarylist, ndays)
    
    
    
    # print header
    print
    print "Report generated %s for processing finished in the previous %d days" % (time.asctime(), ndays)
    print
    
    # at this point, if there are no new runs, newsummarylist will be empty;
    #   if so, we don't want to notify or send email; use script success/failure
    #   as a signal to Hudson to send or not send email (which will be the 
    #   console output of this script)
    
    if not newsummarylist:
        # no jobs, no notify = failure
        print "(none)"
        sys.exit(1)
        
    else:
        # loop over files and print info
        for summaryfilepath in newsummarylist:
            summaryinfo = parsesummaryfile(summaryfilepath)
            print formatsummaryinfo(summaryinfo)
        
        # success means send email:
        sys.exit(0)
        
