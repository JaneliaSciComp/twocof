"""

class and functions for working with spine data


djo, 5/10

"""


# ------------------------- imports -------------------------
# stdlib
import gzip
import os
import pprint
import re


# numerics
import numpy


# from twocof
import constants as const
import trackdata
import util


# ------------------------- constants -------------------------
# number of nodes in spine and columns in the spine file (time plus
#   x, y per node)
nspinenodes = 11
nspinecolumns = 1 + 2 * nspinenodes


# ------------------------- readspinefile() -------------------------
def readspinefile(filename, combinerflag=True):
    """
    read in a data file; it's expected to be in a specific 
    format (tab sep values with gaps, specific column structure)
    
    input:  filename; combiner flag
    output: SpineData object 
    """
    
    filenamedata = trackdata.parsefilename(filename)
    if filenamedata.scalar != "spine":
        raise ValueError("%s does not appear to be a spine file!" % filenamedata.scalar)
    
    if combinerflag:
        return readuncombinedspinerfile(filename)
    
    if not os.path.exists(filename):
        raise util.NoDataError("filename %s doesn't exist" % filename)
    
    data = open(filename).readlines()
    
    # some .r files can be empty; raise an exception
    if len(data) == 0:
        raise util.NoDataError("file %s contains no data!" % filename)
    
    # examine the first line to determine which format it is:
    if data[0].startswith('#'):
        stanzas = {}
        for filepath, datalines in trackdata.stanzaiterator(data):
            # store by "source"; we adjust the source line for later transformation; 
            #   note that this is not quite the same format as for
            #   the other scalars! we're going to fake it to look
            #   like it came from .dat files
            # this involves (1) changing .spine to .dat, and
            #   (2) inserting a ".scalar" where it should be
            # so it starts as: "(stuff).00001.spine", becomes
            #   "(stuff).scalar.00001.dat", to match what comes
            #   out of the "normal" .r files (so we can correlate
            #   with them by source)
            sourcename = filepath.strip()
            sourcename = sourcename.replace(".spine", ".dat")
            # looks like the "@." combination will only occur at
            #   the point I want to insert:
            sourcename = sourcename.replace("@.", "@.scalar.")
            
            numbers = [line.rstrip('\n').split() for line in datalines]
            stanzas[sourcename] = numpy.array(numbers, dtype=numpy.float)
        
        return SpineData(stanzas, filenamedata=filenamedata)
        
    else:
        # don't support the old format for spine files
        raise ValueError("%s appears to be in the old format; not supported for spine files" 
            % filenamedata.scalar)
    
    # end readspinefile()

# ------------------------- _readuncombinedspinedata() -------------------------
def _readuncombinedspinedata(fakerfilepath, choredirpath):
    """
    read data from uncombined spine data
    
    input: (possibly fake) path to .r file; path to choreography results
    output: TrackData object
    """
    
    filenamedata = trackdata.parsefilename(fakerfilepath)
    
    # get datestamp folders
    datestamplist = [fn for fn in os.listdir(choredirpath) if re.search(const.datestampregex, fn)]
    
    
    data = {}
    for datestamp in datestamplist:
        datestampfolder = os.path.join(choredirpath, datestamp)
        
        # data may be text or gzip, check:
        filename = "%s@%s@%s@%s.%s" % (datestamp, filenamedata.geneeff, filenamedata.tracker,
            filenamedata.stimno, filenamedata.scalar)
        filepath = os.path.join(datestampfolder, filename)
        if os.path.exists(filepath):
            filedata = open(filepath).readlines()
        else:
            # gzipped:
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
            sourcename = "/groups/zlatic/zlaticlab/Data/results/%s/%s/%s/%s/%s@%s@%s@%s@.scalar.%d.dat" % (filenamedata.tracker,
                filenamedata.geneeff, filenamedata.stimno, datestamp, datestamp,
                filenamedata.geneeff, filenamedata.tracker, filenamedata.stimno,
                idnumber)
            data[sourcename] = datestampdata[datestampdata[:, 0] == idnumber][:, 1:]
    
    # at this point, check whether we have any data at all:
    if data:
        return SpineData(data, filenamedata=filenamedata)
    else:
        raise util.NoDataError("%s had no usuable data" % fakerfilepath)
    
    # end _readuncombinedspinedata()

# ------------------------- readuncombinedspinerfile() -------------------------
def readuncombinedspinerfile(rfilepath):
    """
    read spine data from uncombined data
    
    input: file path to real or ficticious spine.r file
    output: SpineData from uncombined data
    """
    
    # construct the directory and call the other routine:
    
    if not os.path.isabs(rfilepath):
        raise ValueError("input path %s must be an absolute path!" % rfilepath)
    
    # could just replace "combiner-results" with "choreography-results" in
    #   this file's dir path, but let's be a little more careful
    filenamedata = trackdata.parsefilename(rfilepath)
    choredirname = util.getdirectory("choreography", filenamedata.branch,
        filenamedata.tracker, filenamedata.geneeff, filenamedata.stimno)
    
    return _readuncombinedspinedata(rfilepath, choredirname)
    
    # end readuncombinedspinerfile()

# ------------------------- class SpineData -------------------------
class SpineData(object):
    """
    holds spine data
    """
    # ......................... __init__ .........................
    def __init__(self, data, filenamedata=None):
        """
        
        input:  dictionary of data arrays, keyed by source name
                optional filenamedata struct
        """
        
        # raw input data
        self._rawdata = data
        self.sources = self._rawdata.keys()
        
        # this holds a struct with info from the parsed filename
        self.filenamedata = filenamedata
        
        # end __init__()
    
    # ......................... anglebetweenbysource() .........................
    def anglebetweenbysource(self, source, seg1, seg2):
        """
        find the angle between two segments (clockwise = positive)
        for a given source
        
        segment = (node1, node2) = index of nodes (0 offset)
        
        input: two segments (0-offset); source string
        output: list of [t, angle] (2 x 1-D arrays, not 1 x 2D array) 
        """
        
        # nodes = 0 to 10 index of which node in spine
        node1, node2 = seg1
        node3, node4 = seg2
        
        
        sourcedata = self._rawdata[source]
        t = sourcedata[:, 0]
        
        # previously, we filtered to remove NaN's which we 
        #   added; don't think they exist anymore, but I'm 
        #   always surprised by things:
        tfinite = numpy.isfinite(t)
        t = t[tfinite]
        
        # the dx1, dy1, etc. are the coordinates of the
        #   difference vector for each segment; careful with
        #   the indexing here!
        dx1 = (sourcedata[:, 1 + 2 * node2] - 
            sourcedata[:, 1 + 2 * node1])[tfinite]
        dy1 = (sourcedata[:, 2 + 2 * node2] - 
            sourcedata[:, 2 + 2 * node1])[tfinite]
        dx2 = (sourcedata[:, 1 + 2 * node4] - 
            sourcedata[:, 1 + 2 * node3])[tfinite]
        dy2 = (sourcedata[:, 2 + 2 * node4] - 
            sourcedata[:, 2 + 2 * node3])[tfinite]
        
        # calculate the absolute angle using the dot product
        theta = numpy.arccos((dx1 * dx2 + dy1 * dy2) / 
            numpy.sqrt((dx1 * dx1 + dy1 * dy1) * (dx2 * dx2 + dy2 * dy2)))
        
        # get the sign of the angle from the cross product; according to
        #   Nicholas Swierczek, one of the worm tracker authors,
        #   their coordinate system is with the origin in the top left,
        #   so their natural positive angle is clockwise (x-->y axis);
        #   that's what we want to report, so the sign off the cross 
        #   product is OK:
        sign = numpy.where(dx1 * dy2 - dy1 * dx2 > 0, 1, -1)
        
        return [t, sign * theta]
        
        # end anglebetweenbysource()
    
    # ......................... findtimerange() .........................
    def findtimerange(self):
        """
        locate the min/max time stamps
        
        input: 
        output: tmin, tmax
        """
        
        # rely on correct ordering of the time column:
        minimum = min(data[0][0] for data in self._rawdata.values())
        maximum = max(data[-1][0] for data in self._rawdata.values())
        
        return minimum, maximum
        
        # end findtimerange()
    
    # ......................... getntracks() .........................
    def getntracks(self):
        """
        returns the number of tracks, which is the number of sources
        """
        
        return len(self._rawdata)
        
        # end getntracks()
    
    # ......................... printinfo() .........................
    def printinfo(self):
        """
        print some data
        """
        
        print "basic data info:"
        print "number of tracks = %d" % self.getntracks() 
        
        tmin, tmax = self.findtimerange()
        print "time range = [%f, %f]" % (tmin, tmax)
        
        print "filename data:"
        pprint.pprint(self.filenamedata.entries())
        
        # end printinfo()
    
    # end class SpineData

