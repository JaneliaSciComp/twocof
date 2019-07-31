"""

constants for the twocof library; contains:

-- control parameters
-- time and variable/measurement constants
-- file path stuff
-- experimental data: stimulus mappings, controls 
-- plot related: colors, format, etc.


NOTE: to add a tracker:
- in this file:
-- add lists for primary, normalized, and aggregated scalars
--- add new variables if necessary (see below)
-- add stim2col for the new tracker
-- add tracker root location
- in tasks/__init__.py:
-- add lists for which tasks the tracker uses
- in Hudson:
-- if processed as t1-6 (every two weeks):
--- create jobs
--- add tracker to "manual start" list
- on disk:
-- create directories for tracker in derived, plots, statistics
    for rd and screen; chmod g+w if needed


NOTE: to add a variable:
- in this file:
-- add to scalar lists for appropriate trackers
-- add scale info
-- add unit info


NOTE: to add a task:
- code it within the tasks package, following the existing
    class hierarchy and API
- add task to trackers in tasks/__init__.py


djo, 3/09

"""

# ------------------------- imports -------------------------

from version import __version__


# ------------------------- control -------------------------
# these constants control how the scripts are run (esp. on 
#   the cluster)

# when we run jobs on the cluster, we need to reserve the
#   appropriate number of CPUs; the details are cluster-
#   dependent; currently we segregate into large and small
#   datasets (usually control and non-control)
# note that our cluster has 8 CPU per node and 24G per node;
#   given those values, we have chosen:
slotsperlargedataset = 8
slotspersmalldataset = 1

# related: what constitutes a large dataset?  if you have
#   more than some number of datestamp folders in Choreography
#   output, that's large:
# (see script countdates.py; survey of current data shows
#   5 to be a good, though potentially conservative, threshold 
largedatathreshold = 5


# who pays?  when the appropriate command-line parameter
#   is added, this string is passed as an account code
#   using -A to the qsub command:
accountstring = "zlatic"


# ------------------------- times -------------------------

# flag whether to use constant dt for resampled data or not:
constantresampleddt = True

# if we use a constant dt, what is it? (in seconds)
resampleddt = 0.05 


# seconds to (usually) ignore at the beginning of a trial
ignoretime = 5.0

# minimum difference between interval endpoints to consider 
#   the interval non-degenerate
intervalepsilon = 0.01


# ------------------------- variables -------------------------
# Marta calls the data types "variables"; I call them "scalars"
#   they're also called "measurements"
# NOTE: I delete and/or add scalars at Marta's request; we are
#   dependent on what Choreography outputs, so this needs to be
#   done in coordination with the rest of the pipeline
# NOTE 2: added dependence on tracker relatively late in the
#   project

# the variables are divided into "primary" (.r files
#   provided by Rob from Choreography by the combiner), and
#   "derived", which are .r files that I generate

# Marta's original group:
basevariablelist = [
    "area", 
    "bias", 
    "crabspeed", 
    "curve", 
    "dir", 
    "kink", 
    "midline", 
    "morpwidth", 
    "speed085", 
    "x", 
    "y",
    ]
primaryscalardata = {}
primaryscalardata["t2"] = basevariablelist
primaryscalardata["t3"] = basevariablelist
primaryscalardata["t4"] = basevariablelist
primaryscalardata["t5"] = basevariablelist
primaryscalardata["t6"] = basevariablelist
primaryscalardata["t14"] = basevariablelist
primaryscalardata["t15"] = basevariablelist
primaryscalardata["t16"] = basevariablelist
primaryscalardata["t150"] = basevariablelist
primaryscalardata["t159"] = basevariablelist

# t1, t111 are different:
primaryscalardata["t1"] = [
    "area", 
    "bias", 
    "cast", 
    "crabspeed",
    "curve",
    "dir",
    "headx", 
    "heady", 
    "kink",
    "midline", 
    "morpwidth", 
    "speed", 
    "speed085",
    "tailvecx", 
    "tailvecy", 
    "x", 
    "y", 
    ]
primaryscalardata["t111"] = primaryscalardata["t1"]

# larval guests; t7-t12:
baselarvallist = [
    "area", 
    "cast", 
    "flux", 
    "headx", 
    "heady", 
    "midline", 
    "morpwidth", 
    "nosex", 
    "nosey", 
    "pathlength", 
    "speed", 
    "tailline", 
    "tailvecx", 
    "tailvecy", 
    "tailx", 
    "taily", 
    "x", 
    "y", 
    ]
primaryscalardata["t7"] = baselarvallist
# primaryscalardata["t8"] = baselarvallist
# Nov 2016: t8 reused
primaryscalardata["t8"] = basevariablelist
primaryscalardata["t9"] = baselarvallist
primaryscalardata["t10"] = baselarvallist
primaryscalardata["t11"] = baselarvallist

# t12 is different
primaryscalardata["t12"] = ["area", "cast", "midline", "speed", "x", "y"]


# my test trackers:
# note my test t1 is old, like t2
primaryscalardata["t1-don-test"] = primaryscalardata["t2"]
primaryscalardata["t2-don-test"] = primaryscalardata["t2"]
primaryscalardata["t3-don-test"] = primaryscalardata["t3"]


# we specify explicitly which ones we *want* to normalize for each tracker:
normalizedscalars = {}
basenormlist = ["area", "midline", "morpwidth", "speed085"]
normalizedscalars["t2"] = basenormlist
normalizedscalars["t3"] = basenormlist
normalizedscalars["t4"] = basenormlist
normalizedscalars["t5"] = basenormlist
normalizedscalars["t6"] = basenormlist
normalizedscalars["t14"] = basenormlist
normalizedscalars["t15"] = basenormlist
normalizedscalars["t16"] = basenormlist
normalizedscalars["t150"] = basenormlist
normalizedscalars["t159"] = basenormlist

# t1, t111 are different:
normalizedscalars["t1"] = ["area", "midline", "morpwidth", "speed"]
normalizedscalars["t111"] = ["area", "midline", "morpwidth", "speed"]

# larval guests, t7-12:
baselarvalnormlist = ["area", "midline", "morpwidth", "speed"]
normalizedscalars["t7"] = baselarvalnormlist
# Nov 2016: t8 reused
# normalizedscalars["t8"] = baselarvalnormlist
normalizedscalars["t8"] = basenormlist
normalizedscalars["t9"] = baselarvalnormlist
normalizedscalars["t10"] = baselarvalnormlist
normalizedscalars["t11"] = baselarvalnormlist

normalizedscalars["t12"] = ["area", "midline", "speed"]


# my test trackers:
# note my test t1 is old, like t2
normalizedscalars["t1-don-test"] = normalizedscalars["t2"]
normalizedscalars["t2-don-test"] = normalizedscalars["t2"]
normalizedscalars["t3-don-test"] = normalizedscalars["t3"]


# derived scalars are no longer listed here; they are generated
#   by the tasks (that is, each task knows what new scalars it
#   generates, and it tells process1.py what they are so it
#   can process them)

# I also no longer keep track of all the other files that
#   I don't process


# scalars that we aggregate statistics for:
aggregatedstatsscalarlist = {}
baseaggregatedlist = [
    "area", "areanorm",
    "crabspeed",
    "curve", 
    "dir",
    "midline", "midlinenorm",
    "speed085", "speed085norm",
    ]
aggregatedstatsscalarlist["t2"] = baseaggregatedlist
aggregatedstatsscalarlist["t3"] = baseaggregatedlist
aggregatedstatsscalarlist["t4"] = baseaggregatedlist
aggregatedstatsscalarlist["t5"] = baseaggregatedlist
aggregatedstatsscalarlist["t6"] = baseaggregatedlist
aggregatedstatsscalarlist["t14"] = baseaggregatedlist
aggregatedstatsscalarlist["t15"] = baseaggregatedlist
aggregatedstatsscalarlist["t16"] = baseaggregatedlist
aggregatedstatsscalarlist["t150"] = baseaggregatedlist
aggregatedstatsscalarlist["t159"] = baseaggregatedlist

# t1, t111:
aggregatedstatsscalarlist["t1"] = [
    "area", "areanorm",
    "midline", "midlinenorm",
    "speed", "speednorm",
    ]
aggregatedstatsscalarlist["t111"] = aggregatedstatsscalarlist["t1"]

# larval guests, t7-t12:
aggregatedstatsscalarlist["t7"] = [
    "area", "areanorm",
    "midline", "midlinenorm",
    "speed", "speednorm",
    ]
# Nov 2016: t8 reused
# aggregatedstatsscalarlist["t8"] = aggregatedstatsscalarlist["t7"]
aggregatedstatsscalarlist["t8"] = baseaggregatedlist
aggregatedstatsscalarlist["t9"] = aggregatedstatsscalarlist["t7"]
aggregatedstatsscalarlist["t10"] = aggregatedstatsscalarlist["t7"]
aggregatedstatsscalarlist["t11"] = aggregatedstatsscalarlist["t7"]

aggregatedstatsscalarlist["t12"] = [
    "area", "areanorm",
    "midline", "midlinenorm",
    "speed", "speednorm",
    ]


# my test trackers:
# note my test t1 is old, like t2
aggregatedstatsscalarlist["t1-don-test"] = aggregatedstatsscalarlist["t2"]
aggregatedstatsscalarlist["t2-don-test"] = aggregatedstatsscalarlist["t2"]
aggregatedstatsscalarlist["t3-don-test"] = aggregatedstatsscalarlist["t3"]


# plots are scaled to a more-or-less constant resolution;
#   this is the table of the standard value range per scalar:

# numbers were chosen by looking at a couple datasets and
#   making reasonable guesses; typically, units are:
#   s for time
#   mm, sq mm, mm/s for length, area, speed
#   degrees for angles
#   unitless for biases, cosines, etc.

scalarscaledata = {
    # primary
    "t": [0, 90],
    "area": [2, 5],
    "areanorm": [0.5, 1.5],
    "bias": [-1.5, 1.5],
    "cast": [-120, 120],
    "crabspeed": [-1, 5],
    "crabspeednorm": [0, 1.5],
    "curve": [0, 40],
    "curvenorm": [0.5, 1.5],
    "dir": [-1, 2],
    "dirS": [0, 1],
    "flux": [-1, 1],
    "headx": [-40, 40],
    "heady": [-40, 40],
    "kink": [0, 120],
    "kinknorm": [0, 1.0],
    "midline": [1, 5],
    "midlinenorm": [0.5, 1.5],
    "morpwidth": [0, 1],
    "nosex": [0, 2000],
    "nosey": [0, 2000],
    "orientation": [-180, 180],
    "pathlength": [-100, 100],
    "speed": [0, 3],
    "speednorm": [0.5, 1.5],
    "speed085": [0, 3],
    "speed085norm": [0, 2],
    "tailline": [-1, 1],
    "tailvecx": [-50, 50],
    "tailvecy": [-50, 50],
    "tailx": [0, 2000],
    "taily": [0, 2000],
    "vel_x": [-5, 5],
    "vel_y": [-5, 5],
    "x": [0, 60],
    "y": [0, 60],
    # derived
    "abs_accel_x": [0, 10],
    "abs_accel_y": [0, 10],
    "abscast": [0, 120],
    "abs_vel_x": [0, 5],
    "abs_vel_y": [0, 5],
    "angle": [-45, 45],
    "costheta12": [-1, 1],
    "norm": [-2, 2],
    "theta12": [0, 180],
    "vel_x_bias": [-1, 1],
    "vel_y_bias": [-1, 1],
    "vel_xy_bias": [-1, 1],
    }
def getscalarscale(scalar):
    """
    retrieve scale data for a scalar
    
    -- "norm" is special; if xxxxnorm has no scale,
        use norm scale instead (not xxxx)
    -- "normperrun" is special; treat it as "norm"
    
    """
    if scalar in scalarscaledata:
        return scalarscaledata[scalar]
    elif scalar.endswith("normperrun"):
        return getscalarscale(scalar.replace("normperrun", "norm"))
    else:
        items = scalar.split("-")
        if scalar.endswith("norm"):
            return scalarscaledata["norm"]
        elif scalar in scalarscaledata:
            return scalarscaledata[scalar]
        else:
            raise ValueError("couldn't find scale for scalar %s" % scalar) 

# units of each scalar, so we can have nice axes labels:
scalarunits = {
    "t": "(s)",
    "abs_accel_x": "(mm/s^2)",
    "abs_accel_y": "(mm/s^2)",
    "abs_vel_x": "(mm/s)",
    "abs_vel_y": "(mm/s)",
    "area": "(sq mm)",
    "bias": "",
    "cast": "(degrees)",
    "crabspeed": "(mm/s)",
    "curve": "(degrees)",
    "dir": "",
    "dirS": "",
    "flux": "",
    "headx": "",
    "heady": "",
    "kink": "(degrees)",
    "midline": "(mm)",
    "morpwidth": "(mm)",
    "norm": "",
    "nosex": "",
    "nosey": "",
    "orientation": "(degrees)",
    "pathlength": "(mm)",
    "speed": "(mm/s)",
    "speed085": "(mm/s)",
    "tailline": "",
    "tailvecx": "",
    "tailvecy": "",
    "tailx": "",
    "taily": "",
    "vel_x": "(mm/s)",
    "vel_y": "(mm/s)",
    "x": "(mm)",
    "y": "(mm)",
    "abscast": "(degrees)",
    "angle": "(degrees)",
    "costheta12": "",
    "theta12": "(degrees)",
    "vel_x_bias": "",
    "vel_y_bias": "",
    "vel_xy_bias": "",
    }
def getscalarunit(scalar):
    """
    retrieve scale data for a scalar
    
    -- "norm" (and "normperrun") is special; if present, it takes precedence
    
    """
    if "norm" in scalar:
        return scalarunits["norm"]
    elif scalar in scalarunits:
        return scalarunits[scalar]
    else:
        raise ValueError("couldn't find scale for scalar %s" % scalar) 
    

# ------------------------- file  paths -------------------------
# root directories; all directory hierarchies are assumed to
#   be identical, but the root can change based on tracker;
#   tracker names must be unique

# Marta's lab share--original location of all processing
zlaticlabdir = "/nrs/zlatic/zlaticlab"
		
# larval olympiad area
larvalolympiaddir = "/groups/larvalolympiad/larvalolympiad"

# my test area in my home directory:
olbrisdir = "/home/olbrisd/projects/zlatic"

# for now I assume that both branches of a tracker are in the
#   same location
trackerroots = {
    # Marta's trackers
    "t1": zlaticlabdir,
    "t2": zlaticlabdir,
    "t3": zlaticlabdir,
    "t4": zlaticlabdir,
    "t5": zlaticlabdir,
    "t6": zlaticlabdir,
    "t14": zlaticlabdir,
    "t15": zlaticlabdir,
    "t16": zlaticlabdir,
    "t111": zlaticlabdir,
    "t150": zlaticlabdir,
    "t159": zlaticlabdir,
    
    # larval guests:
    "t7": larvalolympiaddir,
    # Nov 2016: t8 reused
    # "t8": larvalolympiaddir,
    "t8": zlaticlabdir,
    "t9": larvalolympiaddir,
    "t10": larvalolympiaddir,
    "t11": larvalolympiaddir,
    "t12": larvalolympiaddir,
    
    # my test trackers
    "t1-don-test": olbrisdir,
    "t2-don-test": olbrisdir,
    "t3-don-test": olbrisdir,
    }


# this is the top of our file hierarchy
pipelinedir = "pipeline" 

branches = [
    "screen",
    "rd",
    "memory",
    ]
defaultbranch = branches[0]

# these are individual places whose directories can be asked for

locations = {
    "choreography": "choreography-results",
    "source": "combiner-results",
    "derived": "derived",
    "plots": "plots",
    "statistics": "statistics",
    "tables": "tables",
    "histograms": "histograms",
    "analysis": "targeted-analysis",
    "logs": "logs",
    }

# this is a specialized directory within the stats directory
aggregatedstatsdir = "aggregated-statistics"


# umask to use when writing files & dirs; need to be sure other
#   people can run the pipeline and overwrite existing data;
#   this opens things up to the group:
umask = 0002


# regular expression for our date/time stamps:
datestampregex = "^\d{8}_\d{6}$"


# processing status file; location is a key
#   into the above locations dictionary
processingstatusfilename = "_processing-status.json"
processingstatuslocation = "plots"


# location of matplotlib config dir, when run on the cluster:
matplotlibconfigdir = "/groups/zlatic/zlaticlab/code/twocof-production/matplotlibconfig"


# ------------------------- experimental data -------------------------

# stimulus mappings
# which columns in the data correspond to which stimuli for each
#   tracker
# these column numbers used to be 0-indexed in each block of 6 in the
#   .r files, with columns 0 and 1 being time and value of variable;
#   however, when we started working with the uncombined .dat files, the
#   columns shifted over two (to accommodate the time/date stamp and animal
#   ID); since the data structures below work on a .r file that is 
#   reconstructed in memory, the column numbers are technically correct,
#   but they are misleading if you look at the .dat files on disk

stim2col = {}
stim2col["t1"] = {
    'r': 2,     # channelrhodopsin
    'p': 3,     # puff
    'b': 4,     # buzz
    'f': 5,     # (tentative) second puff
    'l': 5,     # light
    'g': 5,     # green
    'a': 5,     # blue ("azure", since b = buzz)
    }
stim2col["t111"] = stim2col["t1"]

stim2col["t2"] = {
    'f': 2,     # (tentative) second puff
    'p': 3,     # puff
    'b': 4,     # buzz
    'i': 5,     # IR laser
    }
stim2col["t3"] = {
    'r': 2,     # channelrhodopsin
    'l': 2,     # light
    'a': 2,     # blue ("azure")
    'p': 3,     # puff
    'h': 4,     # hot
    'c': 4,     # cold
    'b': 4,     # buzz
    }

# note: this may change!  b shouldn't map to two columns; however,
#   since the stim onset error checking is turned off, may not
#   matter
stim2col["t4"] = {
    'p': 3,     # puff
    'b': 4,     # buzz
    'b': 5,     # buzz
    }

stim2col["t5"] = {
    't': 2,     # tap ('t' will likely be changed; we use 't' for temperature)
    'p': 3,     # puff
    'b': 4,     # buzz
    }

stim2col["t6"] = {
    'r': 2,     # channelrhodopsin
    'p': 3,     # puff
    'b': 4,     # buzz
    'g': 5,     # green light (probable)
    }

stim2col["t14"] = {
    'r': 2,     # channelrhodopsin
    'p': 3,     # puff
    'b': 4,     # buzz
    'g': 5,     # green light (probable)
    }

stim2col["t15"] = {
    'r': 2,     # channelrhodopsin
    'p': 3,     # puff
    'b': 4,     # buzz
    'g': 5,     # green light (probable)
    }

stim2col["t150"] = {
    'r': 2,     # channelrhodopsin
    'p': 3,     # puff
    'b': 4,     # buzz
    'g': 5,     # green light (probable)
    }

stim2col["t159"] = stim2col["t15"]
stim2col["t16"] = stim2col["t15"]

# larval olympiad guests, t7-t12: appear to have no stimuli other
#   than from the ignored list:
stim2col["t7"] = { }
# Nov 2016: t8 reused
# stim2col["t8"] = { }
stim2col["t8"] = stim2col["t15"]
stim2col["t9"] = { }
stim2col["t10"] = { }
stim2col["t11"] = { }
stim2col["t12"] = { }


# get the tracker list before adding the test trackers
trackers = stim2col.keys()


# these stimuli are continuous throughout the whole
#   trial, and I do not take any action for them; may 
#   also include stimuli I ignore for other reasons
#   (note: includes 'n' = no stimuli at all)

ignoredstimuli = [
    'n',        # none
    'y',        # yeast
    'o',        # hole
    's',        # scratch
    't',        # temperature
    ]


# ...and for testing:
stim2col["t1-don-test"] = stim2col["t1"]
stim2col["t2-don-test"] = stim2col["t2"]
stim2col["t2-don-test"]['r'] = 2
stim2col["t3-don-test"] = stim2col["t1-don-test"]

testtrackers = ["t1-don-test", "t2-don-test", "t3-don-test"]

# which genotypes are controls?  any which start with the following
#   strings (all converted to lower case when compared):
# NOTE: never finished, unused!
controls = [
    "w1118",
    "pbd",
    ]



# ------------------------- plot related -------------------------

# file format for plots, expressed as file extension; matplotlib
#   handles most common types easily, including svg
# plotextension = ".png"
plotextension = ".svg"


# we try to keep a constant scale for the plots, but sometimes
#   data varies so much that some plots end up distorted and
#   entirely too tall and narrow to be useful; put a limit on
#   those (this is y/x):
aspectratiolimit = 5.0


# constants for plot stimulus markers; when stimulus is short,
#   you need more saturation to see it, but when stimulus is long,
#   large saturation is too obtrusive (my subjective opinion)
basealpha = 0.3
thinalpha = 0.4


# order of colors for multiple stimuli:
stimuluscolorlist = ['r', 'g', 'b', 'y']

# note that 30 is darker than 50; use 30 so it
#   shows up better now that I've raised the alpha
#   on the stimulus shading
gray30 = (0.3, 0.3, 0.3, 1.0)
gray50 = (0.5, 0.5, 0.5, 1.0)
gray70 = (0.7, 0.7, 0.7, 1.0)

graylinecolor = gray30
majorgridcolor = gray30
minorgridcolor = gray70


