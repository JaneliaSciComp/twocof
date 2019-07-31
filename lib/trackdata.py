"""

classes and functions for dealing with the actual data

originally part of the monolithic "twocof.py"


djo, 3/10


"""


# ------------------------- imports -------------------------
# std lib
import gzip
import os
import pprint
import random
import re
import sys


# ***** NOTE *****
# these imports are left here for a reason!  even though this
#   module no longer contains plots, it tends to be imported
#   early, and we need to get this cleanly imported with our
#   preferred backend as soon as possible; leaving it here
#   helps ensure that this is the case


# numerical and graphical
import numpy
import matplotlib

# backend choice (if specified) must precede import of pyplot;
#   check for SGE environment variable to detect if we're on a node
#if "SGE_O_HOST" in os.environ:
#    # we're on a node
matplotlib.use('Agg')

#else:
#    # use the default, probably TkAgg
#    pass

import matplotlib.pyplot as plt


# inside the lib
import constants as const
import util


# ------------------------- constants -------------------------

# set the seed for reproducibility (my office phone #)
numpy.random.seed(2094656)


# tolerance for comparing onset times (seconds):
onsettolerance = 0.1

# minimum separation for onset clusters (seconds):
onsetclusterseparation = 1.0


# ------------------------- parsefilename() -------------------------
def parsefilename(filename):
    """
    grab info out of filename; return in a dictionary
    
    basically: 
    (genotype)@(effector)@(tracker)@(stimulus1)#(stimulus2)#(stimulus3)#(stimulus4)@(animalno).(scalar).r
    
    with:
    (type)_(specs)_(onset)s(number)x(duration)s(interval)s
    
    onset, duration, interval may be floating point (ie, may include '.')
        (this is slightly inconvenient for me, but better than alternatives)
    
    input: filename
    output: struct with filename info
    """
    
    results = util.struct()
    results.filename = filename
    
    results.branch = "unknown"
    results.geneeff = "not set"
    results.stimno = "not set"
    
    # first, grab the branch from the full path (it's not in the
    #   filename, but we need it)
    filepathitems = filename.split('/')
    for item in filepathitems:
        if item in const.branches:
            results.branch = item
            break
    
    # strip off the directory part, then find the parts; parse as much
    #   as you can; unlike earlier versions, try to be more strict;
    #   raise exceptions rather than returning incomplete results
    
    filename = os.path.basename(filename)
    
    # first split is on @:
    atparts = filename.split('@')
    if len(atparts) != 5:
        return results
    
    results.genotype = atparts[0]
    results.effector = atparts[1]
    results.geneeff = "%s@%s" % (results.genotype, results.effector)
    
    results.tracker = atparts[2]
    
    # skip to the last part next:
    animalno, scalar, ext = atparts[4].split('.')
    results.animalno = animalno
    results.stimno = "%s@%s" % (atparts[3], results.animalno)
    results.scalar = scalar
    
    
    # there are four stimulus blocks; need to parse each and 
    #   store data; only parses timed stimuli, not stimuli
    #   that span the whole experiment
    stimlist = atparts[3].split('#')
    pattern = r"([\d\.]+)s(\d+)x([\d\.]+)s([\d\.]+)s"
    if len(stimlist) != 4:
        raise ValueError("expected four stimulus fields, found %s" % stimlist)
    
    
    # first make a determination if all four stimulus blocks are 
    #   populated by stimuli we pay attention to and 
    #   store the results
    results.nstimuli = len([stimfield for stimfield in stimlist 
        if stimfield[0] not in const.ignoredstimuli])
    
    
    # process each block and store results:
    results.stimuli = []
    results.ignoredstimuli = []
    for i, stimfield in enumerate(stimlist):
        if stimfield[0] not in const.ignoredstimuli:
            # store data in a stimulus structure:
            stimulus = util.struct()
            stimcode, ignored2, times = stimfield.split('_')
            
            # we're now in regex territory:
            ans = re.search(pattern, times)
            onset, count, duration, interval = ans.groups()
            
            # stimulus column is fixed per tracker:
            stimulus.column = const.stim2col[results.tracker][stimcode]
            stimulus.stimcode = stimcode
            stimulus.stimfield = stimfield
            stimulus.onset = float(onset)
            stimulus.count = int(count)
            stimulus.pulsed = (stimulus.count > 1)
            stimulus.duration = float(duration)
            
            # the interpretation of "interval" has been an issue since the 
            #   beginning; my notes online and my memory, plus notes on
            #   paper from May 5, 2010, say the interval is the time between
            #   stim offset and next onset (for pulsed stim); however, 
            #   users often thought of it as time between successive onsets,
            #   and by early 2011, the new time intervals specified by
            #   Marta and implemented by Gennady reflect that change in 
            #   definition; my solution: record the "onset delta" here,
            #   so I can quickly adapt to whichever definition is the
            #   current fashion
            stimulus.interval = float(interval)
            stimulus.onsetdelta = stimulus.interval
            # the original, now deprecated alternative:
            # stimulus.onsetdelta = stimulus.interval + stimulus.duration
            
            results.stimuli.append(stimulus)
        else:
            # record the "ignored" stimuli (except 'n'), just
            #   in case you need to know later on:
            if stimfield[0] != 'n':
                results.ignoredstimuli.append(stimfield[0])
    
    return results
    
    # end parsefilename()

# ------------------------- parseseqconcat() -------------------------
def parseseqconcat(data, filenamedata):
    """
    parse data from a sequentially concatenated file
    
    note that source line has the form
    
    input: data, filenamedata structure
    output: trackdata object
    """
    
    filenamedata.readformat = "sequential concatenated"
    
    stanzas = {}
    for filepath, datalines in stanzaiterator(data):
        # filepath is not quite right for a source identifier; need
        #   to strip out the scalar/variable name and leave a generic;
        #   the generic is a handy placeholder will we later replace
        #   when writing out our own derived .r files (in effect,
        #   creating a dummy source fitting our usual naming scheme)
        sourcename = filepath.strip()
        sourcename = sourcename.replace(".%s." % filenamedata.scalar, ".scalar.")
        
        numbers = [line.rstrip('\n').split() for line in datalines]
        stanzas[sourcename] = numpy.array(numbers, dtype=numpy.float)
    
    return TrackData(stanzas, filenamedata=filenamedata)
    
    # end parseseqconcat()

# ------------------------- readfile() -------------------------
def readfile(filename, combinerflag=True):
    """
    read in a data file; it's expected to be in a specific 
    format (tab sep values with gaps, specific column structure)
    
    we support multiple historical formats at this time, although
    the older formats are not guaranteed to work for all tasks
    anymore
    
    if combiner flag is not none, try to read uncombined data
    (details below); this is to help transition to uncombined
    data over .r files
    
    input:  filename; combiner flag
    output: TrackData object 
    """
    
    filenamedata = parsefilename(filename)
    
    # derived scalars are still being written as .r files, remember;
    #   try to read uncombined first (if requested), then fall 
    #   back to .r file method:
    
    if combinerflag:
        try:
            return readuncombinedrfile(filename)
        except:
            # kind of dangerous, but we don't care the details of the error
            pass
    
    if not os.path.exists(filename):
        raise util.NoDataError("filename %s doesn't exist" % filename)
    
    data = open(filename).readlines()
    
    # some .r files can be empty; raise an exception
    if len(data) == 0:
        raise NoDataError("file %s contains no data!" % filename)
    
    # examine the first line to determine which format it is:
    if data[0].startswith('#'):
        # these .dat files have been concatentated sequentially,
        #   with separating lines: '# filepath' 
        return parseseqconcat(data, filenamedata)
    else:
        # these are assumed to be the old way, in which .dat files
        #   are concatenated into one big grid of numbes:
        
        # no longer supported!
        raise ValueError("old columnar .r files not supported!")
    
    # end readfile()

# ------------------------- readheader() -------------------------
def readheader(filename):
    """
    read in the header info for a data file; returns an empty
    TrackData object with no data except for the filename data
    (use at your own risk!  some methods will work, not all)
    
    input: filename of .r file
    output: "empty" TrackData object
    """
    
    filenamedata = parsefilename(filename)
    return TrackData(None, filenamedata=filenamedata)
    
    # end readheader()

# ------------------------- readuncombinedchoredir() -------------------------
def readuncombinedchoredir(choredirpath, scalar):
    """
    read uncombined data from a choreography directory
    
    input: path to choreagraphy output; scalar
    output: TrackData object
    """
    
    # this function isn't used now (May 2011), but it is preferable that
    #   we move toward using it in the future
    
    # untested!
    
    if not os.path.isabs(choredirpath):
        raise ValueError("input path %s must be an absolute path!" % choredirpath)
    
    # parse dirname into its parts an construct a fake .r 
    #   filepath; this is fragile and not a good idea in the long run
    
    # recall that abs path starts with '/', and items[0] is empty string
    items = choredirpath.split('/')
    fakerfilename = "%s@%s@%s.%s.r" % (items[8], items[7], items[9], scalar)
    
    # again, fragile: simple replacement to get right directory path:
    fakerfiledir = choredirpath.replace(const.locations["choreography"],
        const.locations["source"])
    fakerfilepath = os.path.join(fakerfiledir, fakerfilename)
    
    return _readuncombineddata(fakerfilepath, choredirpath)
    
    # end readuncombinedchoredir()

# ------------------------- _readuncombineddata() -------------------------
def _readuncombineddata(fakerfilepath, choredirpath):
    """
    read data from uncombined data
    
    input: (possibly fake) path to .r file; path to choreography results
    output: TrackData object
    """
    
    filenamedata = parsefilename(fakerfilepath)
    filenamedata.readformat = "uncombined"
    
    # get datestamp folders
    datestamplist = [fn for fn in os.listdir(choredirpath) if re.search(const.datestampregex, fn)]
    
    
    data = {}
    for datestamp in datestamplist:
        datestampfolder = os.path.join(choredirpath, datestamp)
        
        # data may be text or gzipped; check for each:
        filename = "%s@%s@%s@%s.%s.dat" % (datestamp, filenamedata.geneeff, filenamedata.tracker,
            filenamedata.stimno, filenamedata.scalar)
        filepath = os.path.join(datestampfolder, filename)
        if os.path.exists(filepath):
            filedata = open(filepath).readlines()
        else:
            # if gzipped:
            gzipfilepath = filepath + ".gz"
            if os.path.exists(gzipfilepath):
                filedata = gzip.open(gzipfilepath).readlines()
            else:
                # for goodnumber = 0 runs, the datestamp folder
                #   exists, but has no usable .dat files, so this isn't an exception
                continue
        
        numbers = [line.rstrip('\n').split()[1:] for line in filedata]
        datestampdata = numpy.array(numbers, dtype=numpy.float)
        
        idlist = numpy.unique(datestampdata[:, 0])
        
        for idnumber in idlist:
            # NOTE: this source name looks like a file path, but it isn't; at
            #   one time it was, but I kept the form because it's (a) backward
            #   compatible, and (b) a nice unique string anyway
            # thus this string does not need to be adjusted when we change
            #   file schemes; contrariwise, it's of less use for tracking data
            sourcename = "/groups/zlatic/zlaticlab/Data/results/%s/%s/%s/%s/%s@%s@%s@%s@.scalar.%d.dat" % (filenamedata.tracker,
                filenamedata.geneeff, filenamedata.stimno, datestamp, datestamp,
                filenamedata.geneeff, filenamedata.tracker, filenamedata.stimno,
                idnumber)
            data[sourcename] = datestampdata[datestampdata[:, 0] == idnumber][:, 1:]
    
    # at this point, check whether we have any data at all:
    if data:
        return TrackData(data, filenamedata=filenamedata)
    else:
        raise util.NoDataError("%s had no usable data" % fakerfilepath)
    
    # end _readuncombineddata()

# ------------------------- readuncombinedrfile() -------------------------
def readuncombinedrfile(rfilepath):
    """
    read the data corresponding to a .r file path, but read from the
    uncombined data instead of the .r file (which in fact need not
    exist)
    
    input: file path to real or fictional .r file
    output: TrackData from uncombined data
    """
    
    # construct the directory and call the other routine:
    
    if not os.path.isabs(rfilepath):
        raise ValueError("input path %s must be an absolute path!" % rfilepath)
    
    # could just replace "combiner-results" with "choreography-results" in
    #   this file's dir path, but let's be a little more careful
    filenamedata = parsefilename(rfilepath)
    choredirname = util.getdirectory("choreography", filenamedata.branch,
        filenamedata.tracker, filenamedata.geneeff, filenamedata.stimno)
    
    return _readuncombineddata(rfilepath, choredirname)
    
    # end readuncombinedrfile()

# ------------------------- stanzaiterator() -------------------------
def stanzaiterator(data):
    """
    auxilliary routine for parsing files that are sequentially 
    concatenated; returns data one stanza at a time
    
    input: whole data file (list of lines)
    output: iterator that returns (filepath, list of data lines)
    """
    
    data = iter(data)
    firstline = data.next()
    # path is prefaced by '# '
    currentpath = firstline[2:]
    currentlist = []
    
    # having initialized things by consuming the first line, now 
    #   iterate through the rest:
    for line in data:
        if line.startswith('#'):
            # next file; empty out this one and prepare for the next:
            yield currentpath, currentlist
            currentpath = line[2:]
            currentlist = []
        else:
            currentlist.append(line)
        
    yield currentpath, currentlist
    
    # end stanzaiterator()

# ------------------------- class TrackData -------------------------
class TrackData(object):
    """
    class for holding data from multiple tracks
    """
    # ......................... constants .........................
    # minimum increment for interpolation
    dtincrement = 0.05
    
    # column assignments
    # time column assumed first (0)
    valuecolumn = 1
    stimcolumnrange = [2, 3, 4, 5]
    
    # stimwidth cutoff; if stimulus onset range is
    #   greater than this, assume it's "pulsed" stimuli
    #   (ie, multiple onsets)
    stimonsetwidthcutoff = 5.0
    
    # ......................... __init__ .........................
    def __init__(self, data, filenamedata=None):
        """
        note that without the filenamedata struct, lots of things
        won't work
        
        note that without the sources dictionary, calculations 
        requiring correlation of different .r files won't be done
        (eg, vel_bias calcs)
        
        input:  dict of data arrays, keyed by source
                optional filenamedata struct
        """
        
        # raw input data
        self._rawdata = data
        
        # this is convenient; allow for the "header only" case:
        if self._rawdata is not None:
            self.sources = self._rawdata.keys()
        else:
            self.sources = []
        
        # this holds a struct with info from the parsed filename
        self.filenamedata = filenamedata
        
        # will hold resampled data
        self._grid = None
        
        
        # memoize some values:
        self._tmin = None
        self._tmax = None
        
        
        '''
        # consistancy check: not enabled yet!  Marta wants things to run
        #   even if inaccurate at this point
        if not self.checkstimonsets():
            raise ValueError("%s failed stim onset consistancy check!" % self.filenamedata.filename)
        '''
        
        # parse the sources in the file:
        self.parsesources()
        
        # end __init__()
    
    # ......................... checkstimonsets() .........................
    def checkstimonsets(self, verbose=False):
        """
        check for potential errors; try to determine whether the 
        stimulus onsets from the filename are consistent with the
        onsets reported in the data (ie, reported by the tracker)
        
        input:  verbose flag; in verbose mode, print output to screen
                (meant to assist debugging, diagnosing data problems)
        output: boolean: consistant or not?  (True = consistant)
        """
        
        # first, arrange the filename stimuli by column:
        filenamestims = dict((i, None) for i in self.stimcolumnrange)
        for stim in self.filenamedata.stimuli:
            filenamestims[stim.column] = stim
        
        
        # check consistancy for each column 
        overallmatch = True
        for stimcol in self.stimcolumnrange:
            if verbose:
                print
                print "examining column %d" % stimcol
            
            onsetlist = self.findstimonsetlist(stimcol)
            stim = filenamestims[stimcol]
            
            datahasstim = (len(onsetlist) > 0)
            filenamehasstim = (stim is not None)
            
            if not datahasstim and not filenamehasstim:
                # no problems yet, do nothing an move on to the 
                #   next column
                if verbose:
                    print "no stimuli in either"
            elif not datahasstim and filenamehastim:
                if verbose:
                    print "filename has stimulus, data doesn't"
                overallmatch = False
            elif datahasstim and not filenamehasstim:
                if verbose:
                    print "data has stimulus, filename doesn't"
                overallmatch = False
            else:
                # both have stim, time to compare:
                
                # data onsets should be in clusters; need to identify
                #   them, find a center, and compare to expected values
                #   from filename
                
                # onset list from filename:
                filenameonsets = [stim.onset + stim.onsetdelta * i 
                    for i in range(stim.count)]
                
                if verbose:
                    print "filename onsets =", filenameonsets
                
                # debug
                # import pdb; pdb.set_trace()
                
                # clusters from data; note that these are intervals:
                onsetlist = self.findstimonsetlist(stimcol)
                dataonsets = list(self._onsetclusterfinder(onsetlist, onsetclusterseparation)) 
                
                if verbose: 
                    print "data onset ranges:"
                    for t1, t2 in dataonsets:
                        print "%s - %s" % (t1, t2)
                
                # brute force: iterate over both and matches; if anything 
                #   doesn't match, we fail; yes, we're testing all matched
                #   times twice, that's why it's brute force
                
                for t in filenameonsets:
                    localmatch = False
                    for t1, t2 in dataonsets:
                        if t1 - onsettolerance <= t <= t2 + onsettolerance:
                            localmatch = True
                            break
                    if not localmatch:
                        if verbose:
                            print "filename onset %s has no match in data" % t
                        overallmatch = False
                        break
                
                for t1, t2 in dataonsets:
                    localmatch = False
                    for t in filenameonsets:
                        if t1 - onsettolerance <= t <= t2 + onsettolerance:
                            localmatch = True
                            break
                    if not localmatch:
                        if verbose:
                            print "data onset interval %s, %s has no match in filename" % (t1, t2)
                        overallmatch = False
                        break
        
        if verbose:
            print
            if overallmatch:
                print "everything matches"
            else:
                print "mismatches detected"
        return overallmatch
        
        # end checkstimonsets()
    
    # ......................... findaverageperobject() .........................
    def findaverageperobject(self, tmin, tmax):
        """
        calculate a scalar average over the given time window on
        a per-object basis
        
        - uses the raw data
        - switch to 1/(N-1) weighting
        
        input: time endpoints (need not be grid aligned)
        output: average, stddev, ntracks, sum of squares
        """
        
        # loop over tracks, find average, then combine averages
        #   as if they were new data points
        
        total = 0
        sumsquares = 0
        numtracks = 0
        
        # loop over the tracks and see if the data is in the window;
        #   if so, do stats!
        for sourcename in self.sources:
            t, value = self.gettrackbysource(sourcename, resampled=False)
            data = value[numpy.logical_and(t > tmin, t < tmax)]
            
            if data.any():
                trackave = numpy.average(data)
                total += trackave
                sumsquares += (trackave ** 2)
                numtracks += 1
        
        # use 1/(N-1) weighting; calculate with 1/N then
        #   multiply by N/(N-1): 
        if numtracks == 0:
            ave = numpy.NaN
            stddev = numpy.NaN
        elif numtracks == 1:
            ave = total / numtracks
            stddev = 0
        else:
            ave = total / numtracks
            stddev = numpy.sqrt(sumsquares / numtracks - ave * ave)
            stddev *= numpy.sqrt(numtracks / (numtracks - 1))
        
        return ave, stddev, numtracks, sumsquares
        
        # end findaverageperobject()
    
    # ......................... findaverageminmax() .........................
    def findaverageminmax(self, tmin, tmax, kind):
        """
        calculate the average in the given time window of the
        min or max of each track
        
        - uses the raw data
        - switch to 1/(N-1) weighting
        
        input: time endpoints (need not be grid aligned); kind = 'min' or 'max'
        output: average, stddev, ntracks, sum of squares
        """
        
        # which is it?
        if kind == "min":
            extremafunction = numpy.min
        elif kind == "max":
            extremafunction = numpy.max
        else:
            raise ValueError("unknown kind %s" % kind)
        
        # loop over tracks, find min or max, then averages
        total = 0
        sumsquares = 0
        numtracks = 0
        
        # loop over the tracks and see if the data is in the window;
        #   if so, stats!
        for sourcename in self.sources:
            t, value = self.gettrackbysource(sourcename, resampled=False)
            data = value[numpy.logical_and(t > tmin, t < tmax)]
            
            if data.any():
                extrema = extremafunction(data)
                total += extrema
                sumsquares += (extrema ** 2)
                numtracks += 1
        
        # use 1/(N-1) weighting; calculate with 1/N then
        #   multiply by N/(N-1): 
        if numtracks == 0:
            ave = numpy.NaN
            stddev = numpy.NaN
        elif numtracks == 1:
            ave = total / numtracks
            stddev = 0
        else:
            ave = total / numtracks
            stddev = numpy.sqrt(sumsquares / numtracks - ave * ave)
            stddev *= numpy.sqrt(numtracks / (numtracks - 1))
        
        return ave, stddev, numtracks, sumsquares
        
        # end findaverageminmax()
    
    # ......................... findaveragetrack() .........................
    def findaveragetrack(self):
        """
        find the averages over the resampled data
        
        input: none
        output: average, stddev, n (count) arrays, sum of squares, squared sum
        """
        
        # this task is much easier with the large single array; without it,
        #   we need to keep a bunch of accumulator arrays and iterate over
        #   the sources, finding where each track's time series should be
        #   added to the accumulators:
        
        tfullrange = self.gettimegrid()
        count = numpy.zeros(tfullrange.shape, dtype=numpy.int)
        total = numpy.zeros(tfullrange.shape, dtype=numpy.float)
        sumsquares = numpy.zeros(tfullrange.shape, dtype=numpy.float)
        
        for sourcename in self.sources:
            t, value = self.gettrackbysource(sourcename, resampled=True)
            
            # in rare instances, this can fail to return data, so check:
            if len(t) > 0:
                # since we're on the resampled grid, finding the first element
                #   is enough to identify where we are in the full time series:
                i = numpy.where(tfullrange == t[0])[0][0]
                n = len(t)
                count[i:i + n] += 1
                total[i:i + n] += value
                sumsquares[i:i + n] += value * value
        
        # along with the usual stats, throw in sum of squares and square of sum
        #   for convenience
        
        squaredsum = total * total
        ave = total / count
        stddev = numpy.sqrt(sumsquares / count - ave * ave)
        
        return ave, stddev, count, sumsquares, squaredsum
        
        # end findaveragetrack()
    
    # ......................... findaveragewindowed() .........................
    def findaveragewindowed(self, tmin, tmax):
        """
        calculate a scalar average over the given time window
        
        uses the resampled data
        
        input: time endpoints (need not be grid aligned)
        output: average, stddev, count, sum sqrs
        """
        
        total = 0.0
        sumsquares = 0.0
        count = 0
        
        for sourcename in self.sources:
            t, value = self.gettrackbysource(sourcename, resampled=True)
            
            # only accumulate if we have data:
            if len(value) > 0:
                valuewindowed = value[numpy.logical_and(t >= tmin, t <= tmax)]
                
                total += valuewindowed.sum()
                sumsquares += (valuewindowed ** 2).sum()
                count += len(valuewindowed)
        
        if count:
            ave = total / count
            stddev = numpy.sqrt(sumsquares / count - ave * ave)
        else:
            ave = numpy.NaN
            stddev = numpy.NaN
        
        return ave, stddev, count, sumsquares
        
        # end findaveragewindowed()
    
    # ......................... findaveragewindowedraw() .........................
    def findaveragewindowedraw(self, tmin, tmax):
        """
        calculate a scalar average over the given time window
        
        uses the raw data
        
        input: time endpoints (need not be grid aligned)
        output: average, stddev, npoints, ntracks, sum of squares
        """
        
        total = 0
        sumsquares = 0
        numpoints = 0
        numtracks = 0
        
        # no other way to do this than grab each time column
        #   and extract corresponding elements that fall in time
        #   window (checking that they are finite or not!)
        for sourcename in self.sources:
            t, value = self.gettrackbysource(sourcename, resampled=False)
            wherefinite = numpy.isfinite(value) 
            t = t[wherefinite]
            value = value[wherefinite]
            data = value[numpy.logical_and(t > tmin, t < tmax)]
            
            total += data.sum()
            sumsquares += (data ** 2).sum()
            
            # count total number of data points, plus
            #   number of tracks which has any data points
            #   included at all
            numpoints += len(data)
            if len(data) > 0:
                numtracks += 1
        
        ave = total / numpoints
        stddev = numpy.sqrt(sumsquares / numpoints - ave * ave)
        
        return ave, stddev, numpoints, numtracks, sumsquares
        
        # end findaveragewindowedraw()
    
    # ......................... findbaselineaverage() .........................
    def findbaselineaverage(self, stimulus):
        """
        find the average before the stimulus; time returned is
        middle of pre-stim range
        
        input: stimulus
        output: t, ave, stddev, number
        """
        
        # use filename, check for onset < 1s
        if stimulus.onset <= 1.0:
            # no data = no averages
            return 0.0, numpy.NaN, numpy.NaN, 0
        else:
            ave, stddev, count, sumsquares = self.findaveragewindowed(0, stimulus.onset - 1)
            t = 0.5 * (stimulus.onset - 1)
            return t, ave, stddev, count
        
        # end findbaselineaverage()
    
    # ......................... findbinnedaverage() .........................
    def findbinnedaverage(self, tmin, tmax, tstep, kind='summed'):
        """
        find averages in windows
        
        returned time vector is currently left end of interval!
        
        input: time start/stop/step; kind = 'summed' or 'per object'
        output: t, ave, stdev, count vectors
        """
        
        if kind == "summed":
            averagefunction = self.findaveragewindowed
        elif kind == "per object":
            averagefunction = self.findaverageperobject
        else:
            raise ValueError("unknown kind %s" % kind)
        
        if tstep <= 0:
            raise ValueError("tstep must be positive")
        if tmax <= tmin:
            raise ValueError("tmax must be greater than tmin")
        
        timelist = []
        avelist = []
        stddevlist = []
        countlist = []
        
        t = tmin
        halfstep = 0.5 * tstep
        while t < tmax:
            ave, stddev, count, sumsquares = averagefunction(t, t + tstep)
            timelist.append(t)
            # could (should?) use center of interval:
            # timelist.append(t + halfstep)
            avelist.append(ave)
            stddevlist.append(stddev)
            countlist.append(count)
            t += tstep
        
        return timelist, avelist, stddevlist, countlist
        
        # end findbinnedaverage()
    
    # ......................... findlineraw() .........................
    def findlineraw(self, tmin, tmax):
        """
        find the least-squares best line through the raw data in
        the given time interval
        
        input: time interval
        output: (slope, intercept, corr coeff squared, count) of best line
        """
        
        # let time = x and value = y, then do standard fit;
        #   probably a clever numpy way to do this that I
        #   don't see (eg, numpy.polyfit(x, y, 1), if I could get
        #   all my data into one vector easily)
        sumx = 0
        sumy = 0
        sumx2 = 0
        sumy2 = 0
        sumxy = 0
        count = 0
        
        # no other way to do this than grab each time column
        #   and extract corresponding elements that fall in time
        #   window (checking that they are finite or not!)
        for sourcename in self.sources:
            t, value = self.gettrackbysource(sourcename, resampled=False)
            wherefinite = numpy.isfinite(value) 
            t = t[wherefinite]
            value = value[wherefinite]
            inrange = numpy.logical_and(t > tmin, t < tmax)
            t = t[inrange]
            data = value[inrange]
            
            sumx += t.sum()
            sumx2 += (t * t).sum()
            sumy += data.sum()
            sumy2 += (data * data).sum()
            sumxy += (t * data).sum()
            count += len(data)
        
        if count:
            xave = sumx / count
            yave = sumy / count
        else:
            xave = numpy.NaN
            yave = numpy.NaN
        
        # definitions taken from mathworld.wolfram.com
        slope = (sumxy - count * xave * yave) / (sumx2 - count * xave * xave)
        intercept = yave - slope * xave 
        
        temp = (sumx2 - count * xave * xave) * (sumy2 - count * yave * yave)
        corrcoeff2 = (sumxy - count * xave * yave) ** 2 / temp
        
        return slope, intercept, corrcoeff2, count
        
        # end findlineraw()
    
    # ......................... findstimcolumns() .........................
    def findstimcolumns(self):
        """
        find which stimulus columns have actual stimuli in them,
        as reported by the tracker
        
        input: 
        output: list of which columns have stimulus in any of the data
        """
        
        stimfound = {}
        for col in self.stimcolumnrange:
            stimfound[col] = False
        for sourcename in self.sources:
            for col in self.stimcolumnrange:
                if not stimfound[col]:
                    stimfound[col] = self._rawdata[sourcename][:, col].any()
        return [col for col in stimfound if stimfound[col]] 
        
        # end findstimcolumns()
    
    # ......................... findstimonset() .........................
    def findstimonset(self, column):
        """
        find the time range of stimulus onset in a given stimulus column
        
        input: column (allowed range currently 2 <= column <= 5)
        output: tmin, tmax
        """
        
        stimonset = self.findstimonsetlist(column)
        if not stimonset:
            raise ValueError("no stimulus onsets found in data for column %s" % column)
        return min(stimonset), max(stimonset)
        
        # end findstimonset()
    
    # ......................... findstimonsetlist() .........................
    def findstimonsetlist(self, column):
        """
        return a list of all the onsets found in a given stimulus column
        
        input: column (allowed range currently 2 <= column <= 5)
        output: list of times (sorted)
        """
        
        # using tricky numpy indexing here!  find values in time column
        #   for which a relation on the stimuls column is true:
        
        onsetlist = []
        for sourcename in self.sources:
            onsetlist.extend(self._rawdata[sourcename][:, 0][self._rawdata[sourcename][:, column] == 1])
        return sorted(onsetlist)
        
        # end findstimonsetlist()
    
    # ......................... findtimerange() .........................
    def findtimerange(self):
        """
        locate the min/max time stamps
        
        input: 
        output: tmin, tmax
        """
        
        if self._tmin is None:
            # rely on correct ordering of the time column:
            self._tmin = min(data[0][0] for data in self._rawdata.values())
            self._tmax = max(data[-1][0] for data in self._rawdata.values())
        
        return self._tmin, self._tmax
        
        # end findtimerange()
    
    # ......................... findvaluerange() .........................
    def findvaluerange(self):
        """
        find the min/max of the values
        
        input:
        output: valuemin, valuemax
        """
        
        minimum = min(min(data[:, self.valuecolumn]) for data in self._rawdata.values())
        maximum = max(max(data[:, self.valuecolumn]) for data in self._rawdata.values())
        
        return minimum, maximum
        
        # end findvaluerange()
    
    # ......................... getchoreographydatestamps() .........................
    def getchoreographydatestamps(self):
        """
        get list of date stamps based on the list in the choreography output
        """
        
        return util.getchoreographydatestamps(self.filenamedata.branch,
            self.filenamedata.tracker, self.filenamedata.geneeff,
            self.filenamedata.stimno)
        
        # end getchoreographydatestamps()
    
    # ......................... getdatablockbysource() .........................
    def getdatablockbysource(self, source):
        """
        get data block by source name in its usual tabular form (always raw data)
        
        input: source name
        output: array
        """
        
        return self._rawdata[source]
        
        # end getdatablockbysource()
    
    # ......................... getdatastride() .........................
    def getdatastride(self):
        """
        current this is a constant
        
        NOTE: given how I use this, it would be getter to attach
        it to self.stride and have that call this automagically;
        or, make self.stride a constant that I *could* attach to
        a property later if I need to
        """
        
        return 6
        
        # end getdatastride()
    
    # ......................... getdatawindowedraw() .........................
    def getdatawindowedraw(self, tmin, tmax):
        """
        returns an array of the data, unresampled; in successivepairs of
        columns, you will have (t, value); NaNs will be all over the place
        
        input: time range
        output: array
        """
        
        # NOTE: this routine is currently only used in output of tables,
        #   which is (a) disabled and (b) not a good idea to start with;
        #   don't bother updating this routine until it's needed!
        
        raise NotImplementedError("not fixed for internal non-rectangular data yet")
        
        # no other way to do this than grab each time column
        #   and extract corresponding elements that fall in time window
        
        # data won't be rectangular, ugh, so need to go through it
        #   once to find the max length and how many columns have data
        #   in the given range
        # length per columns depends on both (a) is time in range,
        #   and (b) is data finite (ie, not NaN)?
        ntracks = self.getntracks()
        stride = self.getdatastride()
        
        maxlength = 0
        totaltracks = 0
        for i in range(ntracks):
            t = self._rawdata[:, stride * i]
            value = self._rawdata[:, stride * i + self.valuecolumn]
            inrange = numpy.logical_and(t > tmin, t < tmax)
            isfinite = numpy.isfinite(value)
            both = numpy.logical_and(inrange, isfinite)
            npts = len(t[both])
            maxlength = max(maxlength, npts)
            if npts > 0:
                totaltracks += 1
        
        # create new array; default value is NaN
        data = numpy.zeros((maxlength, 2 * totaltracks))
        data[:, :] = numpy.NaN
        currenttrack = 0
        for i in range(ntracks):
            t = self._rawdata[:, stride * i]
            value = self._rawdata[:, stride * i + self.valuecolumn]
            inrange = numpy.logical_and(t > tmin, t < tmax)
            isfinite = numpy.isfinite(value)
            both = numpy.logical_and(inrange, isfinite)
            t = t[both]
            if len(t) > 0:
                data[:len(t), 2 * currenttrack] = t
                data[:len(t), 2 * currenttrack + 1] = value[both]
                currenttrack += 1
        
        return data
        
        # end getdatawindowedraw()
    
    # ......................... getdatestamps() .........................
    def getdatestamps(self):
        """
        returns a list of date stamps (strings of form 20101123_121827)
        """
        
        return self.sourcesbydatestamp.keys()
        
        # end getdatestamps()
    
    # ......................... getnaivedt() .........................
    def getnaivedt(self):
        """
        try to guess a good dt without looking at stimulus at all
        
        note: may be overridden by constant dt flag
        """
        
        if const.constantresampleddt:
            return const.resampleddt
        else:
            # try to make a good guess; look at some columns
            #   and see what kind of time interval we see:
            mindt = 999
            for sourcename in self.sources:
                t, value = self.gettrackbysource(sourcename, resampled=False)
                delta = (t - numpy.roll(t, 1))[1:]
                mindt = min(mindt, delta.min())
            return self.dtincrement * numpy.ceil(mindt / self.dtincrement)
        
        # end getnaivedt()
    
    # ......................... getnstimuli() .........................
    def getnstimuli(self):
        """
        return how many not-ignored stimuli we have (should be
        0 to 4)
        """
        
        return self.filenamedata.nstimuli
        
        # end getnstimuli()
    
    # ......................... getntracks() .........................
    def getntracks(self):
        """
        returns the number of tracks (which is the number tracked
        objects)
        """
        
        return len(self._rawdata)
        
        # end getntracks()
    
    # ......................... gettrackdatabydatestamp() .........................
    def gettrackdatabydatestamp(self, datestamp):
        """
        generate a new TrackData object that only contains data from the
        requested date stamp
        
        NOTE: only data shared between full TrackData and reduced 
        one is the filename data
        
        input: string date stamp (format as in .dat filenames & dirs)
        output: TrackData object with that source's data only
        """
        
        # create a "header-only" trackdata to start, then fill it in ourselves
        td = TrackData(None, filenamedata=self.filenamedata)
        
        td._rawdata = {}
        
        for source in self.sourcesbydatestamp[datestamp]:
            td._rawdata[source] = self._rawdata[source].copy()
        td.sources = td._rawdata.keys()
        td.parsesources()
        
        return td
        
        # end gettrackdatabydatestamp()
    
    # ......................... getrandompoints() .........................
    def getrandompoints(self, npts):
        """
        return a group of random points from the track
        
        distribution is not entirely uniform; a source is
        chosen randomly, then a point within the source; this
        will weight the short sources disproportionately
        
        input: npts
        output: t, val arrays
        """
        
        npersource = {}
        for i in range(npts):
            sample = random.choice(self.sources)
            if sample in npersource:
                npersource[sample] += 1
            else:
                npersource[sample] = 1
        
        tsampled = numpy.zeros((npts, ), dtype=numpy.int)
        valsampled = numpy.zeros((npts, ), dtype=numpy.float)
        
        currindex = 0
        for source in npersource:
            for i in range(npersource[source]):
                ny, nx = self._rawdata[source].shape
                randomindex = random.randrange(ny)
                tsampled[currindex] = self._rawdata[source][randomindex, 0]
                valsampled[currindex] = self._rawdata[source][randomindex, self.valuecolumn]
                currindex += 1
        
        return tsampled, valsampled
        
        # end getrandompoints()
    
    # ......................... getsuggesteddt() .........................
    def getsuggesteddt(self):
        """
        calculate a good dt for reinterpolation
        
        if there is only one stim and it's not pulsed, look
        at the range of stimulus onset to choose; otherwise,
        look at the time spacing of the data
        
        (in practice, it will amount to much the same due
        to the discrete steps I impose)
        
        note: may be overridden by constant dt flag
        
        input: none
        output: dt
        """
        
        if const.constantresampleddt:
            return const.resampleddt
        else:
            if self.getnstimuli() == 1:
                stimulus = self.filenamedata.stimuli[0] 
                if not stimulus.pulsed:
                    # try to estimate from onset:
                    tmin, tmax = self.findstimonset(stimulus.column)
                    stimrange = tmax - tmin
                    if stimrange:
                        return self.dtincrement * numpy.ceil(stimrange / self.dtincrement)
                    else:
                        return self.dtincrement
            # basically, any other case is the naive case:
            return self.getnaivedt()
        
        # end getsuggesteddt()
    
    # ......................... gettimegrid() .........................
    def gettimegrid(self, tmin=None, tmax=None):
        """
        return an array of time points; if min/max not given: 
        
        minimum is 0
        maximum is chosen from data (one time step beyond max)
        
        input: optional time endpoints
        output: time array
        """
        
        if tmin is None:
            tmin = 0
        
        if tmax is None:
            ignored, tmax = self.findtimerange()
        
        dt = self.getsuggesteddt()
        
        newtminindex = numpy.floor(tmin / dt)
        newtmaxindex = numpy.ceil(tmax / dt)
        return numpy.arange(newtminindex, newtmaxindex + 1) * dt
        
        # end gettimegrid()
    
    # ......................... gettrackbysource() .........................
    def gettrackbysource(self, source, resampled=True):
        """
        retrieve a track by source name
        
        input: source; optional resampled flag
        output: tlist, valuelist
        """
        
        if self.sources is None:
            raise ValueError("this dataset doesn't have source information")
        
        origt = self._rawdata[source][:, 0]
        origval = self._rawdata[source][:, self.valuecolumn]
        
        if len(origt) == 0:
            return [ ], [ ]
        
        if resampled:
            tmin, tmax = min(origt), max(origt)
            t = self.gettimegrid(tmin, tmax)
            
            # only use finite values
            wherefinite = numpy.isfinite(origval)
            origt = origt[wherefinite]
            origval = origval[wherefinite]
            
            # TWOCOF-31: numpy 1.4 and 1.6 react differently here
            #   if origt and origval have no finite values, so check;
            #   note above that we have a null value defined that we
            #   can return
            if not wherefinite.any():
                return [ ], [ ]
            val = numpy.interp(t, origt, origval, right=numpy.NaN, left=numpy.NaN) 
            
        else:
            t = origt
            val = origval
        
        # don't forget to strip NaN
        wherefinite = numpy.isfinite(val)
        return t[wherefinite], val[wherefinite]
        
        # end gettrackbysource()
    
    # ......................... hassources() .........................
    def hassources(self):
        """
        does the data contain the source for each track?
        """
        
        return self.sources is not None
        
        # end hassources()
    
    # ......................... hasstimulus() .........................
    def hasstimulus(self):
        """
        is there stimulus in any column (as reported by the tracker)?
        
        note that "continuous" stimuli are ignored; check the 
        filename for those
        """
        
        return len(self.findstimcolumns()) > 0
        
        # end hasstimulus()
    
    # ......................... iscontrol() .........................
    def iscontrol(self):
        """
        is this data a control or not?
        
        NOTE: never finished, unused!
        """
        
        raise NotImplementedError
        
        # examine the genotype and apply rules:
        genotype = self.filenamedata.genotype
        
        # compare case-insensitive:
        if genotype.lower().startswith(const.controls):
            return True
        else:
            return False
        
        # end iscontrol()
    
    # ......................... _onsetclusterfinder() .........................
    def _onsetclusterfinder(self, onsetlist, dtmin):
        """
        return an iterator over a list of intervals (centers? means?) 
        from a list of onsets
        
        input: list of times; minimum separation of clusters
        output: iterator; iterator returns [tmin, tmax] of cluster
        """
        
        onsetlist = iter(onsetlist)
        firsttime = onsetlist.next()
        
        currentcluster = [firsttime]
        
        for t in onsetlist:
            if t - currentcluster[-1] > dtmin:
                # output current cluster, start new
                yield min(currentcluster), max(currentcluster)
                currentcluster = [t]
            else:
                currentcluster.append(t)
        
        yield min(currentcluster), max(currentcluster)
        
        # end _onsetclusterfinder()
    
    # ......................... parsesources() .........................
    def parsesources(self):
        """
        iterate over the sources and generate secondary data structures;
        specifically, look for date-time stamp and separate sources by
        them (each "run")
        """
        
        self.sourcesbydatestamp = {}
        
        for source in self.sources:
            # the source is an absolute path, prefaced by '# "; 
            #   exact form varies based on a number of factors, so
            #   look for the datestamp by regex; it looks like
            #   this: 20101123_121827
            sourceitems = source.split('/')
            datetimestamp = ""
            for item in sourceitems:
                if re.search(const.datestampregex, item):
                    datetimestamp = item
                    break
            
            if not datetimestamp:
                raise ValueError("coudn't find date stamp in source %s" % source)
            
            if datetimestamp in self.sourcesbydatestamp:
                self.sourcesbydatestamp[datetimestamp].append(source)
            else:
                self.sourcesbydatestamp[datetimestamp] = [source]
        
        # end parsesources()
    
    # ......................... printinfo() .........................
    def printinfo(self):
        """
        print basic info to screen about the data
        """
        
        print "basic data info:"
        print "stride = %d" % self.getdatastride()
        print "number of tracks = %d" % self.getntracks() 
        
        tmin, tmax = self.findtimerange()
        print "time range = [%f, %f]" % (tmin, tmax)
        
        vmin, vmax = self.findvaluerange()
        print "value range = [%f, %f]" % (vmin, vmax)
        
        stimcolumns = self.findstimcolumns()
        print "data has stimuli in columns: %s" % stimcolumns 
        
        if stimcolumns:
            for column in stimcolumns:
                stimtmin, stimtmax = self.findstimonset(column)
                print ("stimulus onset range for column %d is [%f, %f]" %
                    (column, stimtmin, stimtmax))
        
        print "naive dt is %f" % self.getnaivedt()
        print "suggested dt is %f" % self.getsuggesteddt()
        
        print "filename data:"
        pprint.pprint(self.filenamedata.entries())
        
        if self.filenamedata.stimuli:
            print "stimulus data:"
            for stim in self.filenamedata.stimuli:
                pprint.pprint(stim.entries())
                print
        
        if self.filenamedata.ignoredstimuli:
            print "ignored stimuli:"
            print self.filenamedata.ignoredstimuli
        
        # disable the check when outputting intervals:
        print "intervals"
        '''
        intervaltypelist = ["basic stats", "line fits", "per-object stats, absolute",
            "xy-motion", "unified rules 1"]
        '''
        intervaltypelist = ["unified rules 1"]
        for kind in intervaltypelist:
            intervallist = util.getintervallist(self, kind, check=False)
            print "%s (%d):" % (kind, len(intervallist))
            print intervallist
            print
        
        # end printinfo()
    
    # end class TrackData

# ------------------------- class TrackDataCache -------------------------
class TrackDataCache(object):
    """
    this class acts as a caching container for multiple trackdata objects;
    originally kept all data cached, but in Nov. 2010, we discovered that
    we were having memory issues, so additional caching schemes were added
    """
    # ......................... __init__ .........................
    def __init__(self, scheme="none"):
        """
        
        input: none
        """
        
        # bookkeeping
        if scheme in ["all", "none"]:
            self.scheme = scheme
        else:
            raise ValueError("unknown scheme %s" % scheme)
        
        
        # internal storage
        self._data = {}
        
        # end __init__()
    
    # ......................... _getdatacacheall() .........................
    def _getdatacacheall(self, filename):
        """
        return data, cache everything (original scheme)
        """
        
        if not os.path.isabs(filename):
            filename = os.path.abspath(filename)
        
        if filename not in self._data:
            self._data[filename] = readfile(filename)
        return self._data[filename]
        
        # end _getdatacacheall()
    
    # ......................... _getdatacachenone() .........................
    def _getdatacachenone(self, filename):
        """
        return data, no caching
        """
        
        if not os.path.isabs(filename):
            filename = os.path.abspath(filename)
        
        return readfile(filename)
        
        # end _getdatacachenone()
    
    # ......................... gettrackdata() .........................
    def gettrackdata(self, filename):
        """
        retrieve track data from file; use cached copy if it exists
        
        input: filename (should be absolute path; will attempt to be
                converted if not)
        output: trackdata
        """
        
        if self.scheme == "all":
            return self._getdatacacheall(filename)
        elif self.scheme == "none":
            return self._getdatacachenone(filename)
        
        # end gettrackdata()
    
    # ......................... memcheck() .........................
    def memcheck(self):
        """
        update memory display; Linux only, since I intend
        to get the info out of the /proc directories
        """
        
        data = open("/proc/%d/status" % os.getpid()).readlines()
        
        # format looks like:
        # VmSize:	   77412 kB
        
        # I'm just going to output exactly what was read (brute force):
        output = ""
        for tag in ["VmSize", "VmRSS"]:
            for line in data:
                line = line.strip()
                if line.startswith(tag):
                    output += "%s\t" % line
        
        return "memory: %s" % output
        
        # end memcheck()
    
    # end class TrackDataCache

# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    sys.exit("nothing to see here")
    
    
