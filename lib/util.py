"""

handy routines useful to all the lib

originally part of the monolithic "twocof.py"


djo, 3/10


"""


# ------------------------- imports -------------------------
# std lib:
import os
import re
import time


import numpy


# json not in Python pre-2.5; try the 3rd party lib (mostly for
#   me on my dev box):
try:
    import json
except ImportError:
    try:
        import simplejson as json
    except ImportError:
        raise ImportError("no JSON library found!")


# from twocof:
import constants as const


# ------------------------- constants -------------------------
perobjectabsoluteintervals = [
    # list edited per Marta's email of Feb. 19, 2011 
    #   (original sources noted below)
    
    # dictated to me by Marta, July 15, 2010
    (14.0, 29.0),
    (30.5, 33.5),
    (31.0, 32.0),
    (33.5, 38.5), 
    (40.5, 43.5),
    (41.0, 42.0),
    (43.5, 48.5), 
    (50.5, 53.5),
    (51.0, 52.0),
    (53.5, 58.5), 
    (60.5, 63.5),
    (61.0, 62.0),
    (63.5, 68.5), 
    (70.5, 73.5),
    (71.0, 72.0),
    (73.5, 78.5), 
    (80.5, 83.5),
    (81.0, 82.0),
    (83.5, 88.5), 
    (90.5, 93.5),
    (91.0, 92.0),
    (93.5, 98.5), 
    (100.5, 103.5),
    (101.0, 102.0),
    (103.5, 108.5), 
    (140.5, 143.5), 
    (141.0, 142.0),
    (143.5, 148.5),
    # added from email from Marta, Oct. 22, 2010:
    (38, 40),
    (45, 59),
    (46, 47),
    (49, 50),
    (53, 55),
    # dupe
    # (61, 62),
    (62, 63),
    (64, 65),
    (68, 70),
    (73, 75),
    (76, 77),
    (79, 80),
    (83, 85),
    (90, 95),
    (96, 101),
    (106, 107),
    (110, 115),
    (121, 122),
    (136, 137),
    (151, 152),
    (150, 154),
    (155, 160),
    (160, 165),
    # added from email from Marta, Nov. 1, 2010:
    (31, 35),
    (31, 41),
    (35, 37),
    (37, 39),
    (47, 49),
    (48, 50),
    (50, 52),
    (51, 53),
    (52, 54),
    # dupe
    # (53, 55),
    (55, 56),
    (55.5, 58.5),
    (56, 57),
    (58, 59.5),
    (57.5, 59.5),
    (61, 63),
    (62, 65),
    (65, 67), 
    (66, 68),
    (67, 69),
    (70, 71),
    (73, 74.5),
    # dupe
    # (76, 77),
    (78, 80),
    (80, 82),
    (81, 83),
    (82, 84),
    # dupe
    # (83, 85),
    (85, 86),
    (85.5, 88.5),
    (86, 87),
    (88, 89.5),
    (93, 95),
    (95, 97),
    (96, 98),
    (97, 99),
    (98, 100),
    (100, 101),
    (103, 104.5),
    # dupe
    # (106, 107),
    (108, 110),
    (110, 112),
    (111, 113),
    (112, 114),
    (113, 115),
    (115, 116),
    (115.5, 118.5),
    (116, 117),
    (118, 119.5),
    # dupe
    # (121, 122),
    (123, 125),
    (125, 127),
    (126, 128),
    (127, 129),
    (128, 130),
    (130, 131),
    (130.5, 133.5),
    (131, 132),
    (133, 134.5),
    (138, 140),
    (140, 142),
    (141, 143),
    (142, 144),
    (143, 145),
    (145, 146),
    (145.5, 148.5),
    (146, 147),
    (148, 149.5),
    (164, 166),
    # added from email from Marta Feb 19, 2011
    (14, 220),
    (75, 85),
    (100, 130),
    (105, 115),
    (155, 165),
    (190, 220),
    (210, 220),
    (265, 275),
    (445, 455),
    ]

# ------------------------- class NoDataError -------------------------
class NoDataError(Exception):
    """
    exception to be raised when a track or spine data object/file
    doesn't have any data or any usable data
    """
    
    pass
    
    # end class NoDataError

# ------------------------- adjuststimuli() -------------------------
def adjuststimuli(td, sourcename, t):
    """
    given a trackdata and a source, and a new time series, reassign
    the stimuli onsets in the trackdata to the nearest time point in
    the new series
    
    NOTE: new time series must be evenly spaced! 
    
    input: trackdata; source name; new time series
    output: list of four new stimulus columns
    """
    
    # look at the tmin/max/dt for this track's output; knowing
    #   the grid that implies, try to reassign the original
    #   stimuli points to the nearest point on the grid (this
    #   glosses over the whole "0.5s window" part)
    tmin = t[0]
    tmax = t[-1]
    dt = t[1] - t[0]
    
    # this is actually too big, but convenient for
    #   maintaining nice indices
    stride = td.getdatastride()
    output = [0] * stride
    
    data = td.getdatablockbysource(sourcename)
    wherefinite = numpy.isfinite(data[:, td.valuecolumn])
    for j in td.stimcolumnrange:
        stimin = data[:, j][wherefinite]
        stimout = numpy.zeros((len(t)))
        
        for stimt in data[:, 0][stimin == 1]:
            if tmin <= stimt <= tmax:
                index = int(numpy.round((stimt - tmin) / dt))
                stimout[index] = 1
        output[j] = stimout
    
    # strip off the unused entries and return: 
    return output[2:]
    
    # end adjuststimuli()

# ------------------------- getchoreographydatestamps() -------------------------
def getchoreographydatestamps(branch, tracker, geneeff, stimno):
    """
    get list of date stamps based on what's in choreography-results
    """
    
    choreoutputdir = getdirectory("choreography", branch,
        tracker, geneeff, stimno)
    
    # should this be an error, or return []?
    if not os.path.exists(choreoutputdir):
        raise ValueError("could not find Choreography output folder %s" % choreoutputdir)
    
    possiblelist = [fn for fn in os.listdir(choreoutputdir) if re.search(const.datestampregex, fn)]
    
    # check each directory; if there is only one file instead of many, it's
    #   probably a run where the so-called "goodnumer" is zero, implying 
    #   that no objects were tracked; if so, do not include that date stamp
    # side note: in a totally successful run, there will be 22 files (as of May 2011):
    checkedlist = []
    for datestamp in possiblelist:
        filelist = [fn for fn in os.listdir(os.path.join(choreoutputdir, datestamp))
            if not fn.startswith('.')]
        # could be one or two goodnumber files (.dat and/or .r)
        # NOTE: Marta is going to discontinue output of goodnumber files, so
        #   failed run ought to have nothing; on the other hand, successful run
        #   should have more than two files anyway, so leave this check in for now
        if len(filelist) > 2:
            checkedlist.append(datestamp)
    
    return checkedlist
    
    # end getchoreographydatestamps()

# ------------------------- getdatestamp() -------------------------
def getdatestamp():
    """
    return today's date stamp in our usual format
    """
    
    return time.strftime("%Y%m%d")
    
    # end getdatestamp()

# ------------------------- getdirectory() -------------------------
def getdirectory(kind, branch, tracker, geneeff=None, stimno=None):
    """
    construct a directory path from pieces
    
    input:  kind = (eg) plots, stats, source, etc (see constants.py for list)
            branch, tracker required
            geneeff optional; stimno optional (ignore if geneeff is None)
    output: path to directory; if geneeff, stimno omitted, those parts 
                not concatenated
    """
    
    if branch not in const.branches:
        raise ValueError("unknown branch %s" % branch)
    
    if kind not in const.locations:
        raise ValueError("unknown location kind %s" % kind)
    
    if kind == "logs":
        # ignore tracker for logs:
        temp = os.path.join(const.trackerroots[tracker], const.pipelinedir, 
            branch, const.locations[kind])
    else:
        temp = os.path.join(const.trackerroots[tracker], const.pipelinedir, 
            branch, const.locations[kind], tracker)
    
    if geneeff is not None:
        temp = os.path.join(temp, geneeff)
        if stimno is not None:
            temp = os.path.join(temp, stimno)
    #print "directory is: " + temp
    return temp
    
    # end getdirectory()

# ------------------------- checkintervallist() -------------------------
def checkintervallist(td, intervallist):
    """
    because Marta wants us to keep running even in the face of
    errors, this routine has become more of a "fix it" than
    a check
    
    fix:
    -- all times before ignoretime become ignoretime
    -- times later than data end time become data end time
    -- remove degenerate intervals (t1 == t2)
    
    exceptions:
    -- start time > end time
    
    input: trackdata object, interval list
    output: laundered, checked interval list
    """
    
    # first: go through all times, and any that are earlier than
    #   the ignore time, move them up to the ignore time:
    intervallist = [(max(t1, const.ignoretime), max(t2, const.ignoretime)) 
        for t1, t2 in intervallist]
    
    # any times after end of data become end of data:
    tmin, tmax = td.findtimerange()
    intervallist = [(min(t1, tmax), min(t2, tmax)) for t1, t2 in intervallist]
    
    # next, remove any degenerate intervals (shouldn't be an error):
    intervallist = [(t1, t2) for t1, t2 in intervallist if abs(t1 - t2) > const.intervalepsilon]
    
    # lastly, check that the remaining invervals are well-ordered:
    for t1, t2 in intervallist:
        if t1 >= t2:
            raise ValueError("interval (%f, %f) is not properly ordered" %
                (t1, t2))
    
    return intervallist
    
    # end checkintervallist()

# ------------------------- getintervallist() -------------------------
def getintervallist(td, kind, check=True):
    """
    return a list of time intervals for various purposes
    
    NOTE: if you change this method, also update (if needed)
    GenerateTableTasks.outputfilepaths() in tasks.py; 
    those filenames need to be in one-to-one 
    correspondance with the list of intervals returned 
    for "basic stats"
    
    input:  trackdata; which kind; flag whether to check the intervals
            (used for debugging)
    output: list of (t1, t2) intervals (sorted)
    """
    
    # first cut: per-object stats, absolute doesn't care about
    #   number of stimuli:
    # (yes, I know this is insane)
    # (yes, I know there are duplicates in Marta's lists, but 
    #   I'm not interested in tracking them down)
    if kind == "per-object stats, absolute":
        intervallist = perobjectabsoluteintervals 
    elif kind == "unified rules 1":
        # note that if there are stimuli, there is also an interval
        #   that depends on something above the individual stim level;
        #   it's "last offset" to + 29, and last offset must be 
        #   determined here (for nstim > 0)
        nstim = td.getnstimuli()
        if nstim == 0:
            lastoffset, intervallist = getintervalunified1(None, td)
            lastoffsetlist = []
        elif nstim == 1:
            stimulus = td.filenamedata.stimuli[0]
            lastoffset, intervallist = getintervalunified1(stimulus, td)
            lastoffsetlist = [lastoffset]
        else:
            intervalset = set() 
            lastoffsetlist = []
            for stimulus in td.filenamedata.stimuli:
                lastoffset, templist = getintervalunified1(stimulus, td) 
                intervalset.update(templist)
                lastoffsetlist.append(lastoffset)
            intervallist = list(intervalset)
        # usually we ignore 's', but if it's there, add some intervals:
        if 's' in td.filenamedata.ignoredstimuli:
            intervallist.append((14, 220))
            intervallist.append((100, 130))
            intervallist.append((190, 220))
        if lastoffsetlist:
            lastoffset = max(lastoffsetlist)
            intervallist.append((lastoffset, lastoffset + 29))
    else:
        # pick up with everything else:
        nstim = td.getnstimuli()
        tmin, tmax = td.findtimerange()
        
        if nstim == 0:
            if kind == "basic stats":
                # no stimulus; one big window
                intervallist = [
                    (0.0, tmax),
                    ]
            elif kind == "xy-motion":
                # no stimulus; if there is a "scratch" ignored stimulus,
                #   still need windows (hardcoded):
                if 's' in td.filenamedata.ignoredstimuli:
                    intervallist = [
                        (0.0, 15.0),
                        (15.0, 45.0),
                        (45.0, 70.0),
                        (0.0, tmax),
                        ]
                else:
                    intervallist = [
                        (0.0, tmax),
                        ]
            elif kind == "line fits":
                # (this is currently the same as xy-motion, but I keep
                #   it separate to make it easier to change in the future)
                # no stimulus; if there is a "scratch" ignored stimulus,
                #   still need windows (hardcoded):
                if 's' in td.filenamedata.ignoredstimuli:
                    intervallist = [
                        (0.0, 15.0),
                        (15.0, 45.0),
                        (45.0, 70.0),
                        (0.0, tmax),
                        ]
                else:
                    intervallist = [
                        (0.0, tmax),
                        ]
            else:
                raise ValueError("unknown kind %s" % kind)
        elif nstim == 1:
            # single stim; pulsed or not?
            stimulus = td.filenamedata.stimuli[0]
            if not stimulus.pulsed:
                intervallist = getintervalsingle(stimulus, td, kind)
            else:
                # single stim, pulsed stimulus
                intervallist = getintervalpulsed(stimulus, td, kind)
        else:
            # multiple stimuli; need to grab all the interval lists
            #   and union them:
            intervalset = set() 
            for stimulus in td.filenamedata.stimuli:
                if not stimulus.pulsed:
                    intervallist = getintervalsingle(stimulus, td, kind)
                else:
                    intervallist = getintervalpulsed(stimulus, td, kind)
                intervalset.update(intervallist)
            intervallist = list(intervalset)
    
    if check:
        intervallist = checkintervallist(td, intervallist)
    intervallist.sort()
    return intervallist
    
    # end getintervallist()

# ------------------------- getintervalsingle() -------------------------
def getintervalsingle(stimulus, td, kind):
    """
    get the interval list for a single stimulus
    
    input: stimulus struct; trackdata; kind of stimulus
    output: list of (t1, t2) intervals
    """
    
    tmin, tmax = td.findtimerange()
    
    stimonset = stimulus.onset
    stimduration = stimulus.duration
    stimoffset = stimonset + stimduration
    offsetduration = tmax - (stimonset + stimduration)
    
    if kind == "basic stats":
        intervallist = [
            (0.0, stimonset),
            (stimonset + 1, stimonset + 4),
            (stimonset, stimonset + 0.5 * stimduration),
            (stimonset + 0.5 * stimduration, stimoffset),
            (stimoffset, stimoffset + 0.5 * offsetduration),
            (stimoffset + 0.5 * offsetduration, tmax),
            ]
    elif kind == "line fits":
        intervallist = [
            (stimonset, stimonset + 2.0),
            (stimonset + 2.0, stimonset + 4.0),
            (stimonset + 0.5 * stimduration, stimoffset),
            (stimoffset, tmax),
            (stimonset + 0.5 * stimduration, tmax),
            ]
        # this interval is only valid for long enough durations:
        if stimonset + 4.0 < stimonset + 0.5 * stimduration:
            intervallist.append((stimonset + 4.0, stimonset + 0.5 * stimduration))
    elif kind == "xy-motion":
        intervallist = [
            (0.0, stimonset),
            (stimonset, stimoffset),
            (stimonset, stimonset + 30),
            (stimonset + 30, tmax),
            ]
    else:
        raise ValueError("unknown interval kind %s" % kind)
    
    return intervallist
    
    # end getintervalsingle()

# ------------------------- getintervalpulsed() -------------------------
def getintervalpulsed(stimulus, td, kind):
    """
    get the intervals for a pulsed stimulus
    
    input: stimulus struct; trackdata; kind of stimulus
    output: list of (t1, t2) intervals
    """
    
    tmin, tmax = td.findtimerange()
    
    stimcount = stimulus.count
    stimduration = stimulus.duration
    stimonset = stimulus.onset
    stiminterval = stimulus.interval
    stimdelta = stimulus.onsetdelta
    
    if kind == "basic stats":
        intervallist = [
            (0.0, stimonset),
        ]
        for i in range(stimcount):
            # Marta specified these intervals in email,
            #   Feb 26, 2010
            currentonset = stimonset + i * stimdelta
            intervallist.append((currentonset, currentonset + 3))
            intervallist.append((currentonset, currentonset + stimdelta))
            intervallist.append((currentonset, currentonset + 0.5 * stimdelta))
            intervallist.append((currentonset + 0.5 * stimdelta, currentonset + stimdelta))
        # one more to cover the end; remember that last stimulus ends
        #   partway through last delta:
        intervallist.append((stimonset + (stimcount - 1) * stimdelta + stimduration, tmax))
    elif kind == "line fits":
        # mimic the stats stuff for now:
        intervallist = []
        for i in range(stimcount):
            # Marta specified these intervals in email,
            #   Feb 26, 2010
            currentonset = stimonset + i * stimdelta
            intervallist.append((currentonset, currentonset + 3))
            intervallist.append((currentonset, currentonset + stimdelta))
            intervallist.append((currentonset, currentonset + 0.5 * stimdelta))
            intervallist.append((currentonset + 0.5 * stimdelta, currentonset + stimdelta))
        # one more to cover the end; remember that last stimulus ends
        #   partway through last delta:
        intervallist.append((stimonset + (stimcount - 1) * stimdelta + stimduration, tmax))
    elif kind == "xy-motion":
        # not clear that this is useful; during intervals are too short to
        #   be interesting?
        # also: did NOT add the "stim to stim + 30" and "stim+30 to end"
        #   intervals for all the pulses
        intervallist = [
            (0.0, stimonset),
        ]
        for i in range(stimcount):
            currentonset = stimonset + i * stimdelta
            intervallist.append((currentonset, currentonset + stimduration))
    else:
        raise ValueError("unknown interval kind %s" % kind)    
    
    return intervallist
    
    # end getintervalpulsed()

# ------------------------- getintervalunified1() -------------------------
def getintervalunified1(stimulus, td):
    """
    get the intervals for a stimulus under the unified interval rules; 
    given suffix "1" in case there are major modifications in the future
    
    specs from Marta in email, Feb. 19 2011
    
    input: stimulus struct (may be None); trackdata
    output: last offset time, list of (t1, t2) intervals
    """
    
    # start with hardcoded values
    intervallist = [
        # pre-stim:
        (14.0, 29.0),
        # apparently removed, re: Gennady's email of Feb. 23, 2011
        # end of experiment, allegedly needs to be hardcoded:
        # (75, 85),
        # (105, 115),
        # (155, 165),
        # (210, 220),
        # (265, 275),
        # (445, 455),
        ]
    
    # if we don't have a stimulus, stop here:
    if stimulus is None:
        return None, intervallist
    else:
        stimonset = stimulus.onset
    
    # hardcoded values for "scratch" stim in t5:
    if stimulus.stimcode == 's':
        intervallist.extend([
            (14, 220),
            (100, 130),
            (190, 220),
            ])
    
    # unlike previously, we don't differentiate between "single" 
    #   and "pulsed"; we distinguish by shorter or longer than 15s
    #   (there are short singles which I need to treat like I had
    #   previously treated pulsed)
    if stimulus.duration > 15:
        stimoffset = stimonset + stimulus.duration
        intervallist.extend([
            # added Nov. 2011, see email Nov. 1-2, and Tomoko's
            #   request email Oct. 31, ~4:30 pm:
            (stimonset + 0.25, stimonset + 0.75),
            # from Feb. 2011:
            (stimonset + 1, stimonset + 2),
            (stimonset + 1, stimonset + 5),
            (stimonset + 1, stimonset + 10),
            (stimonset + 5, stimonset + 15),
            (stimonset + 15, stimonset + 29),
            (stimoffset - 6, stimoffset - 1),
            (stimoffset + 2, stimoffset + 5),
            (stimoffset + 2, stimoffset + 10),
            (stimoffset + 10, stimoffset + 20),
            (stimoffset + 19, stimoffset + 29),
            ])
        lastoffset = stimoffset
    else:
        # for pulsed, need to find first and last pulse and work from there
        firstonset = stimonset
        lastonset = stimonset + (stimulus.count - 1) * stimulus.onsetdelta
        
        # first:
        intervallist.extend([
            (firstonset + 0.5, firstonset + 1.5),
            (firstonset + 1, firstonset + 2),
            (firstonset + 1, firstonset + 5),
            (firstonset + 2, firstonset + 4),
            (firstonset + 4, firstonset + 5),
            (firstonset + 9, firstonset + 11),
            (firstonset + 9, firstonset + 15),
            (firstonset + 13, firstonset + 15),
            ])
        
        # last: 
        intervallist.extend([
            # apparently removed, re: Gennady's email of Feb. 23, 2011
            # (lastonset + 0.5, lastonset + 1.5),
            # (lastonset + 1, lastonset + 2),
            # (lastonset + 1, lastonset + 5),
            # (lastonset + 2, lastonset + 4),
            # (lastonset + 4, lastonset + 5),
            # (lastonset + 9, lastonset + 11),
            # (lastonset + 9, lastonset + 15),
            # (lastonset + 13, lastonset + 15),
            (lastonset + 15, lastonset + 30),
            ])
        lastoffset = lastonset + stimulus.duration
    
    return lastoffset, intervallist
    
    # end getintervalunified1()

# ------------------------- getmeminfo() -------------------------
def getmeminfo():
    """
    try to return some memory info 
    
    input: none
    output: dict with memory info
    """
    
    data = open("/proc/meminfo", 'rt').readlines()
    
    memdict = {}
    
    # let's start with the free memory:
    for line in data:
        if line.startswith("MemFree"):
            items = line.split()
            memdict["free"] = int(items[1])
            break
    
    return memdict
    
    # end getmeminfo()

# ------------------------- getresourcelocation() -------------------------
def getresourcelocation(resource):
    """
    returns the location of a given type of resource
    
    input: resource type; allowed: "executable",
    output: file path of directory containing resource
    """
    
    # TODO: figure out best way to do this
    
    libdir = os.path.dirname(__file__)
    packagedir = os.path.dirname(libdir)
    if resource == "executable":
        return os.path.join(packagedir, "bin")
    else:
        raise ValueError("unknown resource type %s" % resource)
    
    # end getresourcelocation()

# ------------------------- islargedata() -------------------------
def islargedata(branch, tracker, geneeff, stimno):
    """
    is the given dataset "large"?  in this context, "large" means
    "needing special handling"; most likely this will be used to
    assign the right number of slots on the cluster to the job
    that processes the dataset
    
    in the future, you could also stat all the files and count
    bytes, but not sure we have enough experience to interpret
    the results
    
    input: branch, tracker, geneeff, stimno 
    output: boolean
    """
    
    # simple test: compare number of choreography output folders
    #   to a threshold:
    
    datestamplist = getchoreographydatestamps(branch, tracker, geneeff, stimno)
    
    return len(datestamplist) > const.largedatathreshold
    
    # end islargedata()

# ------------------------- isnewdata() -------------------------
def isnewdata(branch, tracker, geneeff, stimno):
    """
    has the data in the tracker changed since the last
    recorded succussful run?  

    NOTE: as the name suggests, used to only check for new data;
        now also checks for removed data
    
    output: True available data has changed since last run 
        (data added or removed), or if tracking file is missing
    """
    
    statusdir = getdirectory(const.processingstatuslocation, 
        branch, tracker, geneeff, stimno) 
    statuspath = os.path.join(statusdir, const.processingstatusfilename)
    print statuspath
    if not os.path.exists(statuspath):
        return True
    
    data = json.load(open(statuspath, 'rt'))
    
    if data["status"] != "successful":
        return True
    
    processedset = set(data["dates processed"])
    currentset = set(getchoreographydatestamps(branch, tracker, geneeff, stimno))
    
    # has the set of dates processed changed (anything added or removed)?
    return currentset != processedset
    
    # end isnewdata()

# ------------------------- rectangularize() -------------------------
def rectangularize(inputarray):
    """
    take a list of list of arrays and create one big rectangular
    array; used for creating .r files
    
    input: data list: [(t, val, stim1...stim4), ...]
    output: numpy array of above, rectangular
    """
    
    stride = len(inputarray[0])
    ntracks = len(inputarray)
    
    # assemble into rectangular array
    longest = max(len(data[0]) for data in inputarray)
    outputarray = numpy.zeros((longest, stride * ntracks), dtype=numpy.float)
    outputarray.fill(numpy.NaN)
    
    for i, data in enumerate(inputarray):
        ny = len(data[0])
        for j in range(len(data)):
            # this basically loops over t, value, stim columns:
            outputarray[:ny, i * stride + j] = data[j]
    
    return outputarray
    
    # end rectangularize()

# ------------------------- class rFileWriter -------------------------
class rFileWriter(object):
    """
    a class for writing "version 2" .r files (repeated stanzas of
    source location then columnar tab-separated data)
    """
    # ......................... __init__ .........................
    def __init__(self, filename, NaNstring=None):
        """
        
        input: filename; string to substitue for NaN in data
        """
        
        self.filename = filename
        self.NaNstring = NaNstring
        
        self.f = open(filename, "wt")
        
        # end __init__()
    
    # ......................... close() .........................
    def close(self):
        """
        close the file
        """
        
        self.f.close()
        
        # end close()
    
    # ......................... writearray() .........................
    def writearray(self, array):
        """
        write an array
        """
        
        if len(array.shape) != 2:
            raise ValueError("array must be two-dimensional")
        
        ny, nx = array.shape
        format = '\t'.join(["%g"] * nx) + '\n'
        
        if self.NaNstring is None:
            for row in array:
                self.f.write(format % tuple(row))
        else:
            originalNaN = "%g" % numpy.NaN
            for row in array:
                line = format % tuple(row)
                self.f.write(line.replace(originalNaN, self.NaNstring))
        
        # end writearray()
    
    # ......................... writesourceline() .........................
    def writesourceline(self, source):
        """
        write a source line
        """
        
        self.f.write("# %s\n" % source)
        
        # end writesourceline()
    
    # end class rFileWriter

# ------------------------- class rFileWriterDummy -------------------------
class rFileWriterDummy(object):
    """
    dummy class; can be used in place of an rFileWriter, but doesn't
    acutally do anything
    """
    def __init__(self, filename, NaNstring=None): pass
    
    def close(self): pass
    
    def writearray(self, array): pass
    
    def writesourceline(self, source): pass
    
    # end class rFileWriterDummy

# ------------------------- savetabseparated() -------------------------
def savetabseparated(array, filename, headerlines=None, NaNstring=None):
    """
    save the given array as a tab-separated file; will overwrite
    existing files
    
    input: numpy array and filename; optional list of strings to prepend 
            (strings should include \n's); optional string to 
            substitute for NaNs in data
    output: none (writes file)
    """
    
    if len(array.shape) != 2:
        raise ValueError("array must be two-dimensional")
    
    ny, nx = array.shape
    format = '\t'.join(["%f"] * nx) + '\n'
    
    f = open(filename, 'wt')
    if headerlines is not None:
        for line in headerlines:
            f.write(line)
    if NaNstring is None:
        for row in array:
            f.write(format % tuple(row))
    else:
        originalNaN = "%f" % numpy.NaN
        for row in array:
            line = format % tuple(row)
            f.write(line.replace(originalNaN, NaNstring))
    f.close()
    
    # end savetabseparated()

# ------------------------- class struct -------------------------
class struct(object):
    """
    a simple structure class; just add fields as you need 
    them...this is one of the simpler possible implementations
    
    example:
    
    >>> myproject = struct(name="mine")
    >>> myproject.path = "/path/to/files"
    >>> print myproject.name, myproject.path
    mine /path/to/files
    """
    # ......................... __init__ .........................
    def __init__(self, **entries):
        """
        create a new structure with given entries 
        input:keyword entries & values
        """
        
        self.__dict__.update(entries)
        
        # end __init__()
    
    # ......................... __repr__() .........................
    def __repr__(self):
        """
        same as str
        """
        
        return self.__str__()
        
        # end __repr__()
    
    # ......................... __str__() .........................
    def __str__(self):
        """
        called when an instance needs to be printed
        """
        
        # this formatting is almost certainly too clever...
        
        nentries = len(self.__dict__)
        if nentries == 0:
            return "<struct with no keys>"
        elif nentries <= 3:
            return "<struct with keys %s>" % ", ".join(self.__dict__.keys())
        else:
            return "<struct with keys %s, ...>" % ", ".join(self.__dict__.keys()[:3])
        
        # end __str__()
    
    # ......................... add() .........................
    def add(self, **entries):
        """
        adds a bunch of keyword=value pairs to the struct
        """
        
        self.__dict__.update(entries)
        
        # end add()
    
    # ......................... __contains__() .........................
    def __contains__(self, item):
        """
        is the key contained in the struct?
        """
        
        return item in self.__dict__
        
        # end __contains__()
    
    # ......................... clear() .........................
    def clear(self):
        """
        clears the structure
        """
        
        self.__dict__.clear()
        
        # end clear()
    
    # ......................... copy() .........................
    def copy(self):
        """
        returns a shallow copy of the struct
        """
        
        newcopy = struct()
        newcopy.update(self.entries())
        return newcopy
        
        # end copy()
    
    # ......................... entries() .........................
    def entries(self):
        """
        returns a dictionary with the struct's entries (shallow copy)
        """
        
        return self.__dict__.copy()
        
        # end entries()
    
    # ......................... remove() .........................
    def remove(self, item):
        """
        remove item from struct
        """
        
        del self.__dict__[item]
        
        # end remove()
    
    # ......................... update() .........................
    def update(self, entrydict):
        """
        adds the contents of the input dictionary to the struct
        """
        
        self.__dict__.update(entrydict)
        
        # end update()
    
    # end class struct

# ------------------------- class Timer -------------------------
class Timer(object):
    """
    adapted from examples by David Beazley and Phil Winston
    
    usage:
    
    with Timer() as t:
        do stuff
    print "duration = %s" % t.seconds
    
    """
    def __init__(self, message=None): self.message = message
    
    def __enter__(self): 
        self._start = time.time()
        return self
    
    def __exit__(self, type, value, traceback): 
        self._end = time.time()
        self.seconds = self._end - self._start
        if self.message is not None:
            print "%s: %d" % (self.message, self.seconds)
    
    # end class Timer

# ------------------------- walktracker() -------------------------
def walktracker(args, location=const.locations["source"]):
    """
    given a list of command-line args, identify a tracker and
    optionally a gene/effector and stimuls/animalno component;
    returns a list of (tracker, geneeff, stimno) triples for
    all components not specified
    
    or, basically it'll walk the directory and find all the
    stuff there for all levels below the ones you specify
    
    NOTE: output should probably be either (geneeff, stimno)
        alone, or all four: (branch, tracker, geneeff, stimno)
    
    input:  list of command-line args; optional which location
            of our structure to walk (default = combiner)
    output: list of (tracker, geneeff, stimno)
    """
    
    # at a minimum, we need a tracker and branch
    if len(args) < 2:
        raise ValueError("need a branch and tracker!")
    branch = args[0]
    tracker = args[1]
    
    # does it exist:
    trackerpath = getdirectory(location, branch, tracker)
    if not os.path.exists(trackerpath):
        raise ValueError("%s not found" % trackerpath)
    if not os.path.isdir(trackerpath):
        raise ValueError("%s does not seem to be a directory" % trackerpath)
    
    triplets = []
    
    if len(args) == 2:
        # do everything in tracker:
        for item in os.listdir(trackerpath):
            geneeffpath = os.path.join(trackerpath, item)
            if os.path.isdir(geneeffpath):
                # assume it's one of our folders (gene@effector)
                for item2 in os.listdir(geneeffpath):
                    stimnopath = os.path.join(geneeffpath, item2)
                    if os.path.isdir(stimnopath):
                        # should be stimulus@animalno:
                        triplets.append((tracker, item, item2))
    elif len(args) > 2:
        geneeff = args[2]
        geneeffpath = os.path.join(trackerpath, geneeff)
        # check that it exists:
        if not os.path.exists(geneeffpath):
            raise ValueError("file %s not found" % geneeffpath)
        if not os.path.isdir(geneeffpath):
            raise ValueError("%s does not seem to be a directory" % geneeffpath)
        
        if len(args) == 3:
            # no stim@animalno, do everything in that gene@effector:
            for item in os.listdir(geneeffpath):
                stimnopath = os.path.join(geneeffpath, item)
                if os.path.isdir(stimnopath):
                    # should be stimulus@animalno:
                    triplets.append((tracker, geneeff, item))
        else:
            # it's a specific combination:
            stimno = args[3]
            # check that it exists:
            stimnopath = os.path.join(geneeffpath, stimno)
            if not os.path.exists(stimnopath):
                raise ValueError("file %s not found" % stimnopath)
            if not os.path.isdir(stimnopath):
                raise ValueError("%s does not seem to be a directory" % stimnopath)
            
            triplets = [(tracker, geneeff, stimno)]    
    
    return triplets
    
    # end walktracker()

