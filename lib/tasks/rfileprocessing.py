"""

this file contains routines for .r file processing--tasks that
are run on each and every .r file


djo, 12/10


"""



# ------------------------- imports -------------------------
# std lib
import os

# numerical and graphical
import numpy
import matplotlib

# backend choice (if specified) must precede import of pyplot;
#   check for SGE environment variable to detect if we're on a node
#if "SGE_O_HOST" in os.environ:
# we're on a node

matplotlib.use('Agg')

#else:
# use the default, probably TkAgg
#    pass

import matplotlib.pyplot as plt


# inside the lib
import common
from lib import constants as const
from lib import plots, reports
from lib import trackdata
from lib import util



# ------------------------- class ProcessRFileTask -------------------------
class ProcessRFileTask(common.ProcessingTask):
    """
    
    """
    # ......................... __init__ .........................
    def __init__(self, branch, tracker, geneeff, stimno, scalar, trackdatacache, 
        overwrite=False):
        """
        
        input:  branch: eg, screen
                tracker: eg, t1 or t1-don-text
                stem: eg, w1118@none@t1@buzz@1v100hz15s10x5s25
                scalar: eg, aspect
                trackdata cache instance
                flag whether to overwrite existing or not
        """
        
        # housekeeping
        self.branch = branch
        self.tracker = tracker
        self.geneeff = geneeff
        self.stimno = stimno
        self.scalar = scalar
        self.trackdatacache = trackdatacache
        self.overwrite = overwrite
        
        # generate the stem from the above:
        self.stem = "%s@%s@%s" % (self.geneeff, self.tracker, self.stimno)
        
        # end __init__()
    
    
    # ......................... gettasktarget() .........................
    def gettasktarget(self):
        """
        return target
        """
        
        return "%s -- %s -- %s -- %s" % (self.tracker, self.geneeff, 
            self.stimno, self.scalar)
        
        # end gettasktarget()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        for most r-file processing, this is well defined (can go in
        superclass)
        
        input: none
        output: list of file paths the task will write
        """
        
        filename = "%s.%s%s" % (self.stem, self.scalar, self.suffix)
        directory = util.getdirectory(self.outputlocation, self.branch, self.tracker, 
            self.geneeff, self.stimno) 
        return [os.path.join(directory, filename)]
        
        # end outputfilepaths()
    
    # ......................... readdata() .........................
    def readdata(self):
        """
        read in and return the data
        """
        
        # construct the filename and request the data from the cache 
        
        sourcefilename = "%s.%s.r" % (self.stem, self.scalar)
        if self.scalar in const.primaryscalardata[self.tracker]:
            sourcelocation = "source"
        else:
            sourcelocation = "derived"
        sourcedir = util.getdirectory(sourcelocation, self.branch,
                self.tracker, self.geneeff, self.stimno)
        sourcepath = os.path.join(sourcedir, sourcefilename)
        
        return self.trackdatacache.gettrackdata(sourcepath)        
        
        # end readdata()
    
    # end class ProcessRFileTask

# ------------------------- class NormalizedPerRunReader -------------------------
class NormalizedPerRunReader(object):
    """
    this is a mixin; use for any per-run calculations that should, when
    performed on normalized variables, use per-datestamp normalization
    instead of global normalization
    
    this must precede other classes in the inheritance hierarchy so
    its readdata() method is called
    """
    
    # ......................... readdata() .........................
    def readdata(self):
        """
        read in data; we want to read normalized per-run data,
        so if the scalar ends in "norm", we instead read the
        corresponding "normperrun" data
        """
        
        if not self.scalar.endswith("norm"):
            # fall back to old method:
            return ProcessRFileTask.readdata(self)
        else:
            scalar = self.scalar.replace("norm", "normperrun")
            sourcefilename = "%s.%s.r" % (self.stem, scalar)
            sourcelocation = "derived"
            sourcedir = util.getdirectory(sourcelocation, self.branch,
                    self.tracker, self.geneeff, self.stimno)
            sourcepath = os.path.join(sourcedir, sourcefilename)
            
            return self.trackdatacache.gettrackdata(sourcepath)        
        
        # end readdata()
    
    # end class NormalizedPerRunReader

# ------------------------- class PlotDataSaver -------------------------
class PlotDataSaver(object):
    """
    this is a mixin; use it for any task that generates plots for
    which you want to save the data as well; must precede other
    task classes in the inheritance hierarchy so its outputfilepaths()
    method is called before ProcessRFileTask's
    
    """
    
    # ......................... _getdatafilepath() .........................
    def _getdatafilepath(self):
        """
        returns the data file path
        """
        
        filename = "%s.%s%s" % (self.stem, self.scalar, self.datasuffix)
        return os.path.join(self._getplotsdir(), filename)
        
        # end _getdatafilepath()
    
    # ......................... _getplotsdir() .........................
    def _getplotsdir(self):
        """
        returns the plots dir we're going to output to
        """
        
        return util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        
        # end _getplotsdir()
    
    # ......................... _getplotfilepath() .........................
    def _getplotfilepath(self):
        """
        returns the plot file path
        """
        
        filename = "%s.%s%s" % (self.stem, self.scalar, self.plotsuffix)
        return os.path.join(self._getplotsdir(), filename)
        
        # end _getplotfilepath()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        returns list of output files
        """
        
        return [self._getplotfilepath(), self._getdatafilepath()]
        
        # end outputfilepaths()
    
    # end class PlotDataSaver

# ------------------------- class PerRunPlotDataSaver -------------------------
class PerRunPlotDataSaver(object):
    """
    see PlotDataSaver; this mixin is the same, but for per-run processing
    """
    
    # ......................... _getdatafilepath() .........................
    def _getdatafilepath(self, datestamp):
        """
        returns the data file path
        """
        
        filename = "%s.%s%s" % (self.stem, self.scalar, self.datasuffix)
        return os.path.join(self._getplotsdir(datestamp), filename)
        
        # end _getdatafilepath()
    
    # ......................... _getplotsdir() .........................
    def _getplotsdir(self, datestamp):
        """
        returns the plots dir we're going to output to
        """
        
        baseoutputdir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        return os.path.join(baseoutputdir, datestamp)
        
        # end _getplotsdir()
    
    # ......................... _getplotfilepath() .........................
    def _getplotfilepath(self, datestamp):
        """
        returns the plot file path
        """
        
        filename = "%s.%s%s" % (self.stem, self.scalar, self.plotsuffix)
        return os.path.join(self._getplotsdir(datestamp), filename)
        
        # end _getplotfilepath()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        return list of output this task will write
        """
        
        # look in choreography output to see which datestamps we need
        #   to output for; this way we don't need to read the actual
        #   data file
        
        # primary scalar data must exist
        scalar = const.primaryscalardata[self.tracker][0]
        sourcefilename = "%s.%s.r" % (self.stem, scalar)
        sourcedir = util.getdirectory("source", self.branch,
                self.tracker, self.geneeff, self.stimno)
        sourcepath = os.path.join(sourcedir, sourcefilename)
        
        td = trackdata.readheader(sourcepath)
        plotlist = [self._getplotfilepath(datestamp) for datestamp in td.getchoreographydatestamps()]
        datalist = [self._getdatafilepath(datestamp) for datestamp in td.getchoreographydatestamps()]
        
        return plotlist + datalist
        
        # end outputfilepaths()
    
    # end class PerRunPlotDataSaver

# ------------------------- class GenerateAveErrPlotTask -------------------------
class GenerateAveErrPlotTask(PlotDataSaver, ProcessRFileTask):
    """
    generate plot of ave track with error bars
    """
    # ......................... constants .........................
    plotsuffix = "-ave-err%s" % const.plotextension
    datasuffix = "-ave-err-data.txt"
    outputlocation = "plots"
    
    # ......................... performtask() .........................
    def performtask(self):
        
        plotfilepath = self._getplotfilepath()
        datafilepath = self._getdatafilepath()
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateAveErrPlot!"
            
            td = self.readdata()
            print "Plot File Path: " + plotfilepath            
            fig, statsdata = plots.plotaverageerror(td, showN=True, show=False, 
                filename=plotfilepath)
            # fig.clf() not enough to release memory!  need to use close()
            plt.close(fig)
            
            # now write the data
            header = ["data for %s\n" % os.path.basename(plotfilepath),
                "time\taverage\tstddev\tcount\tsumsqr\tsqrsum\n",
                ]
            dataarray = numpy.column_stack(statsdata)
            util.savetabseparated(dataarray , datafilepath, headerlines=header)
            
            return "wrote %s\nwrote %s" % (os.path.basename(plotfilepath), 
                os.path.basename(datafilepath))
        else:
            return "skipping ave-err plots...output exists"
        
        # end performtask()
    
    # end class GenerateAveErrPlotTask

# ------------------------- class GenerateAveErrPerRunPlotTask -------------------------
class GenerateAveErrPerRunPlotTask(NormalizedPerRunReader, PerRunPlotDataSaver, ProcessRFileTask):
    """
    generate plot of ave in bins using per-object statistics; be
    aware of stimulus onset (one bin should start at onset, previous
    bin should not overlap onset)
    """
    # ......................... constants .........................
    plotsuffix = "-ave-err%s" % const.plotextension
    datasuffix = "-ave-err-data.txt"
    outputlocation = "plots"
    
    # ......................... performtask() .........................
    def performtask(self):
        
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateAveErrPerRunPlotTask!"
            
            td = self.readdata()
            
            # loop over datestamps:
            messagelist = []
            for datestamp in td.getdatestamps():
                plotfilepath = self._getplotfilepath(datestamp)
                datafilepath = self._getdatafilepath(datestamp)
                
                tdreduced = td.gettrackdatabydatestamp(datestamp)
                
                fig, statsdata = plots.plotaverageerror(tdreduced, showN=True, show=False, 
                    filename=plotfilepath)
                # fig.clf() not enough to release memory!  need to use close()
                plt.close(fig)
                
                # now write the data
                header = ["data for %s\n" % os.path.basename(plotfilepath),
                    "time\taverage\tstddev\tcount\tsumsqr\tsqrsum\n",
                    ]
                dataarray = numpy.column_stack(statsdata)
                util.savetabseparated(dataarray , datafilepath, headerlines=header)
                
                messagelist.append("wrote %s\nwrote %s" % (os.path.basename(plotfilepath), 
                        os.path.basename(datafilepath)))
                
            return "\n".join(messagelist)
        else:
            return "skipping ave-err per run plots...output exists"
        
        # end performtask()
    
    # end class GenerateAveErrPerRunPlotTask

# ------------------------- class GenerateAveFitsPlotTask -------------------------
class GenerateAveFitsPlotTask(ProcessRFileTask):
    """
    generate plot of ave track with fit lines
    """
    # ......................... constants .........................
    suffix = "-ave-linefits%s" % const.plotextension
    outputlocation = "plots"
    
    # ......................... performtask() .........................
    def performtask(self):
        
        # there's only one output file path for r file processing:
        filepath = self.outputfilepaths()[0]
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateAveFitsPlot!"
            
            td = self.readdata()
            
            if td.hasstimulus():
                fig = plots.plotaveragefits(td, filename=filepath, show=False)
                plt.close(fig)
                
                return "wrote %s" % os.path.basename(filepath)
            else:
                return "skipping %s...no stimulus" % os.path.basename(filepath)
        else:
            return "skipping %s...output exists" % os.path.basename(filepath)
        
        # end performtask()
    
    # end class GenerateAveFitsPlotTask

# ------------------------- class GenerateBinnedAvePlotTask -------------------------
class GenerateBinnedAvePlotTask(ProcessRFileTask):
    """
    generate plot of ave in bins, with pre-stim interval
    """
    # ......................... constants .........................
    suffix = "-ave-binned%s" % const.plotextension
    outputlocation = "plots"
    
    # ......................... performtask() .........................
    def performtask(self):
        
        # there's only one output file path for r file processing:
        filepath = self.outputfilepaths()[0]
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateBinnedAvePlot!"
            
            td = self.readdata()
            
            if td.hasstimulus():
                # plot the average in bins of 5s, including the pre-stimulus
                #   average, using resampled data (this is a bit messy)
                # if multiple stimuli, use the earliest onset
                binsize = 5
                if td.getnstimuli() > 0:
                    tmin, tmax = td.findtimerange()
                    stimulus = min(td.filenamedata.stimuli,
                        key=lambda x: getattr(x, "onset"))
                    stimonset = stimulus.onset
                    stimcolumn = stimulus.column
                    # prestimvalues: needs to be a list for ave/std/etc,
                    #   not a single value; also, check that there is
                    #   actually data there (count > 0, no NaN's); throw
                    #   in some 0 values if so
                    prestimvalues = td.findbaselineaverage(stimulus)
                    if prestimvalues[-1] == 0 or numpy.isnan(prestimvalues[1]):
                        prestimvalues = [[0], [0], [0], [0]]
                    else:
                        prestimvalues = [[item] for item in prestimvalues]
                    fig, statsdata = plots.plotbinnedaverage(td, stimonset, tmax, 
                        binsize, showN=True,
                        prestimvalues=prestimvalues,
                        filename=filepath, show=False)
                    plt.close(fig)
                return "wrote %s" % os.path.basename(filepath)
            else:
                return "skipping %s...no stimulus" % os.path.basename(filepath)
        else:
            return "skipping %s...output exists" % os.path.basename(filepath)
        
        # end performtask()
    
    # end class GenerateBinnedAvePlotTask

# ------------------------- class GenerateBinnedAvePerObjPlotTask -------------------------
class GenerateBinnedAvePerObjPlotTask(PlotDataSaver, ProcessRFileTask):
    """
    generate plot of ave in bins using per-object statistics; be
    aware of stimulus onset (one bin should start at onset, previous
    bin should not overlap onset)
    """
    # ......................... constants .........................
    plotsuffix = "-ave-binned-perobj%s" % const.plotextension
    datasuffix = "-ave-binned-perobj-data.txt"
    outputlocation = "plots"
    
    # ......................... performtask() .........................
    def performtask(self):
        
        plotfilepath = self._getplotfilepath()
        datafilepath = self._getdatafilepath()
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateBinnedAvePerObjPlot!"
            
            td = self.readdata()
            
            # used to use td.hassstimulus(), but it turns out some 
            #   'ignored' stimuli still have data in the stim columns
            #   in the .r files!  so check the filename:
            if td.filenamedata.stimuli:
                # plot the average in bins of 5s, including the pre-stimulus
                #   average, using resampled data (this is a bit messy)
                # if multiple stimuli, use the earliest onset
                binsize = 5
                tmin, tmax = td.findtimerange()
                stimulus = min(td.filenamedata.stimuli,
                    key=lambda x: getattr(x, "onset"))
                stimonset = stimulus.onset
                stimcolumn = stimulus.column
                if stimonset > binsize:
                    # prestim: starts at t=0, last bin should not overlap 
                    #   stim onset, so that would be max t = stimonset - tstep
                    prestimvalues = list(td.findbinnedaverage(
                        0, stimonset - binsize, binsize, kind="per object"))
                else:
                    # stimulus is too early; no prestim values:
                    prestimvalues = None
                fig, statsdata = plots.plotbinnedaverage(td, 
                    stimonset, tmax, binsize, 
                    kind="per object",
                    prestimvalues=prestimvalues,
                    showN=True, 
                    filename=plotfilepath, 
                    show=False,
                    )
                plt.close(fig)
                
                # now write the data
                # NOTE: writing two dummy columns in preparation for future
                #   sum of squares and squared sum output; it's convenient
                #   to have them now for SAGE loading, and fill with data later
                zerodata = numpy.zeros_like(statsdata[0])
                statsdata = list(statsdata)
                statsdata.extend([zerodata, zerodata])
                
                header = ["data for %s\n" % os.path.basename(plotfilepath),
                    "time\taverage\tstddev\tcount\tsumsqr\tsqrsum\n",
                    ]
                dataarray = numpy.column_stack(statsdata)
                util.savetabseparated(dataarray, datafilepath, headerlines=header)
                
                return "wrote %s\nwrote %s" % (os.path.basename(plotfilepath), 
                    os.path.basename(datafilepath))
            else:
                return "skipping ave-binned-perobj plots...no stimulus"
        else:
            return "skipping ave-binned-perobj plots...output exists"
        
        # end performtask()
    
    # end class GenerateBinnedAvePerObjPlotTask

# ------------------------- class GenerateBinnedAvePerObjPerRunPlotTask -------------------------
class GenerateBinnedAvePerObjPerRunPlotTask(NormalizedPerRunReader, PerRunPlotDataSaver, ProcessRFileTask):
    """
    generate plot of ave in bins using per-object statistics; be
    aware of stimulus onset (one bin should start at onset, previous
    bin should not overlap onset)
    """
    # ......................... constants .........................
    plotsuffix = "-ave-binned-perobj%s" % const.plotextension
    datasuffix = "-ave-binned-perobj-data.txt"
    outputlocation = "plots"
    
    # ......................... performtask() .........................
    def performtask(self):
        
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateBinnedAvePerObjPerRunPlotTask!"
            
            td = self.readdata()
            
            # loop over datestamps:
            messagelist = []
            for datestamp in td.getdatestamps():
                plotfilepath = self._getplotfilepath(datestamp)
                datafilepath = self._getdatafilepath(datestamp)
                
                tdreduced = td.gettrackdatabydatestamp(datestamp)
                
                # used to use td.hassstimulus(), but it turns out some 
                #   'ignored' stimuli still have data in the stim columns
                #   in the .r files!  so check the filename:
                if tdreduced.filenamedata.stimuli:
                    # plot the average in bins of 5s, including the pre-stimulus
                    #   average, using resampled data (this is a bit messy)
                    # if multiple stimuli, use the earliest onset
                    binsize = 5
                    tmin, tmax = tdreduced.findtimerange()
                    stimulus = min(tdreduced.filenamedata.stimuli,
                        key=lambda x: getattr(x, "onset"))
                    stimonset = stimulus.onset
                    stimcolumn = stimulus.column
                    if stimonset > binsize:
                        # prestim: starts at t=0, last bin should not overlap 
                        #   stim onset, so that would be max t = stimonset - tstep
                        prestimvalues = list(tdreduced.findbinnedaverage(
                            0, stimonset - binsize, binsize, kind="per object"))
                    else:
                        # stimulus is too early; no prestim values:
                        prestimvalues = None
                    fig, statsdata = plots.plotbinnedaverage(tdreduced, 
                        stimonset, tmax, binsize, 
                        kind="per object",
                        prestimvalues=prestimvalues,
                        showN=True, 
                        filename=plotfilepath, 
                        show=False,
                        )
                    plt.close(fig)
                    
                    # now write the data
                    # NOTE: writing two dummy columns in preparation for future
                    #   sum of squares and squared sum output; it's convenient
                    #   to have them now for SAGE loading, and fill with data later
                    zerodata = numpy.zeros_like(statsdata[0])
                    statsdata = list(statsdata)
                    statsdata.extend([zerodata, zerodata])
                    
                    header = ["data for %s\n" % os.path.basename(plotfilepath),
                        "time\taverage\tstddev\tcount\tsumsqr\tsqrsum\n",
                        ]
                    dataarray = numpy.column_stack(statsdata)
                    util.savetabseparated(dataarray, datafilepath, headerlines=header)
                    
                    messagelist.append("wrote %s\nwrote %s" % (os.path.basename(plotfilepath), 
                        os.path.basename(datafilepath)))
                else:
                    return "skipping ave-binned-perobj per run plots...no stimulus"
            return "\n".join(messagelist)
        else:
            return "skipping ave-binned-perobj per run plots...output exists"
        
        # end performtask()
    
    # end class GenerateBinnedAvePerObjPerRunPlotTask

# ------------------------- class GenerateHistogramsTask -------------------------
class GenerateHistogramsTask(ProcessRFileTask):
    """
    generate plot of histograms of mean values (per track) in a time window
    """
    # ......................... constants .........................
    plotsuffixpattern = "-ave-hist-%.2f-%.2f" + const.plotextension
    datasuffixpattern = "-ave-hist-%.2f-%.2f-data.txt"
    outputlocation = "histograms"
    
    # ......................... _getdatafilepath() .........................
    def _getdatafilepath(self, tmin, tmax):
        """
        returns the data file path
        """
        
        suffix = self.datasuffixpattern % (tmin, tmax)
        filename = "%s.%s%s" % (self.stem, self.scalar, suffix)
        outputdir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        return os.path.join(outputdir, filename)
        
        # end _getdatafilepath()
    
    # ......................... _getplotfilepath() .........................
    def _getplotfilepath(self, tmin, tmax):
        """
        returns the plot file path
        """
        
        suffix = self.plotsuffixpattern % (tmin, tmax)
        filename = "%s.%s%s" % (self.stem, self.scalar, suffix)
        outputdir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        return os.path.join(outputdir, filename)
        
        # end _getplotfilepath()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        return list of output this task will write
        """
        
        # histograms are unusual in that their output file names
        #   depend on the input data; in this case, time intervals
        #   appear in the output filename, and they are a function
        #   of the stimuli; fortunately, we can get that by parsing
        #   the .r filename (but we don't need to read the data)
        
        # just pick one scalar rather than use the correct scalar, because
        #   for "derived", the correct scalar may point you to a file that
        #   does not exist yet!  (the intervals will all be the same anyway):
        scalar = const.primaryscalardata[self.tracker][0]
        sourcefilename = "%s.%s.r" % (self.stem, scalar)
        sourcedir = util.getdirectory("source", self.branch,
                self.tracker, self.geneeff, self.stimno)
        sourcepath = os.path.join(sourcedir, sourcefilename)
        
        # get an empty trackdata object with real filename info;
        #   that's all the interval routine needs, for the assumed
        #   interval scheme; do NOT check the list, because that
        #   does require a data read (to find tmin, tmax)
        # NOTE: the side effect of this is that since too-long intervals
        #   aren't trimmed, we'll probably always trigger recalculation;
        #   better that be done and the data read once in the task than
        #   read it here to be sure, though
        td = trackdata.readheader(sourcepath)
        intervallist = util.getintervallist(td, "unified rules 1", check=False)
        
        return ([self._getdatafilepath(tmin, tmax) for tmin, tmax in intervallist]
            + [self._getplotfilepath(tmin, tmax) for tmin, tmax in intervallist])
        
        # end outputfilepaths()
    
    # ......................... performtask() .........................
    def performtask(self):
        
        resultstrings = []
        td = self.readdata()
        skipped = 0
        for tmin, tmax in util.getintervallist(td, "unified rules 1"):
            plotfilepath = self._getplotfilepath(tmin, tmax)
            datafilepath = self._getdatafilepath(tmin, tmax)
            if (self.overwrite or not os.path.exists(plotfilepath) or
                not os.path.exists(datafilepath)):
                
                if not self.checkcreateoutputfolders():
                    return "error: could not create output folders for GenerateHistograms!"
                
                fig, numbers, bins = plots.plotavehistogram(td, tmin, tmax, 
                    show=False, filename=plotfilepath)
                # fig.clf() not enough to release memory!  need to use close()
                plt.close(fig)
                
                # construct the data to write out; array of bins = edges
                #   of the bins, so I'll output them explicitly (which
                #   is redundant)
                
                header = ["data for %s\n" % os.path.basename(plotfilepath),
                    "bin-L\tbin-R\tnumber (normalized)\n",
                    ]
                if len(numbers) > 0:
                    dataarray = numpy.column_stack([bins[:-1], bins[1:], numbers])
                else:
                    dataarray = numpy.array([[0, 0, 0]])
                util.savetabseparated(dataarray , datafilepath, headerlines=header)
                
                
                resultstrings.append("wrote %s" % os.path.basename(plotfilepath))
            else:
                # see notes in outputfilepath() above; we're going to be running
                #   this a lot without needing to write any files, so ease up
                #   on the "skipped" messages:
                skipped += 1
                # resultstrings.append("skipping %s...output exists" % os.path.basename(plotfilepath))
        
        if skipped:
            resultstrings.append("skipped %d files...output exists" % skipped)
        return "\n".join(resultstrings)
        
        # end performtask()
    
    # end class GenerateHistogramsTask

# ------------------------- class GenerateLinearFitsTask -------------------------
class GenerateLinearFitsTask(ProcessRFileTask):
    """
    generate linear fits
    """
    # ......................... constants .........................
    suffix = "-linfit.txt"
    outputlocation = "statistics"
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write statistics
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        # there's only one output file path for r file processing:
        filepath = self.outputfilepaths()[0]
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateLinearFits!"
            
            td = self.readdata()
            
            if td.hasstimulus():
                reports.report1lf(td, filepath)
                return "wrote %s" % os.path.basename(filepath)
            else:
                return "skipping %s...no stimulus" % os.path.basename(filepath)
        else:
            return "skipping %s...output exists" % os.path.basename(filepath)
        
        # end performtask()
    
    # end class GenerateLinearFitsTask


# ------------------------- class GenerateAllDataPlotTask -------------------------
class GenerateAllDataPlotTask(ProcessRFileTask):
    """
    generate scatter plot with all the data
    """
    
    # ......................... constants .........................
    suffix = "-alldata%s" % const.plotextension
    outputlocation = "plots"
    
    # ......................... performtask() .........................
    def performtask(self):
        
        # there's only one output file path for r file processing:
        filepath = self.outputfilepaths()[0]
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateAllDataPlot!"
            
            td = self.readdata()
            
            fig = plots.plotalldata(td, show=False)
            fig.set_size_inches(40, 30)
            try:
                memoryproblem = False
                fig.savefig(filepath, dpi=300)
            except RuntimeError:
                # memory issue; try to go on
                memoryproblem = True
            finally:
                # fig.clf() not enough to release memory!  need the close()
                plt.close(fig)
            
            if memoryproblem:
                return "error: runtime error trying to save figure %s" % filepath
            else:
                return "wrote %s" % os.path.basename(filepath)
        else:
            return "skipping %s...output exists" % os.path.basename(filepath)
        
        # end performtask()    
    
    # end class GenerateAllDataPlotTask

# ------------------------- class GenerateScatterFitsPlotTask -------------------------
class GenerateScatterFitsPlotTask(ProcessRFileTask):
    """
    generate plot of randomly chosen scatter plot with fit lines
    """
    # ......................... constants .........................
    suffix = "-scatter-linefits%s" % const.plotextension
    outputlocation = "plots"
    npoints = 1000
    
    # ......................... performtask() .........................
    def performtask(self):
        
        # there's only one output file path for r file processing:
        filepath = self.outputfilepaths()[0]
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateScatterFitsPlot!"
            
            td = self.readdata()
            
            if td.hasstimulus():
                fig = plots.plotrandomraw(td, self.npoints, linefit=True,
                    filename=filepath, show=False)
                plt.close(fig)
                
                return "wrote %s" % os.path.basename(filepath)
            else:
                return "skipping %s...no stimulus" % os.path.basename(filepath)
        else:
            return "skipping %s...output exists" % os.path.basename(filepath)
        
        # end performtask()
    
    # end class GenerateScatterFitsPlotTask

# ------------------------- class GenerateStatsTask -------------------------
class GenerateStatsTask(ProcessRFileTask):
    """
    generate statistics
    """
    
    # ......................... constants .........................
    suffix = "-stats.txt"
    outputlocation = "statistics"
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write statistics
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        # there's only one output file path for r file processing:
        filepath = self.outputfilepaths()[0]
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateStats!"
            
            td = self.readdata()
            
            reports.report1a(td, filepath)
            return "wrote %s" % os.path.basename(filepath)
        else:
            return "skipping %s...output exists" % os.path.basename(filepath)
        
        # end performtask()
    
    # end class GenerateStatsTask

# ------------------------- class GenerateMaxStatsTask -------------------------
class GenerateMaxStatsTask(ProcessRFileTask):
    """
    generate statistics
    """
    
    # ......................... constants .........................
    suffix = "-stats-max.txt"
    outputlocation = "statistics"
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write statistics
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        # there's only one output file path for r file processing:
        filepath = self.outputfilepaths()[0]
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateMaxStats!"
            
            td = self.readdata()
            
            reports.report1max(td, filepath)
            return "wrote %s" % os.path.basename(filepath)
        else:
            return "skipping %s...output exists" % os.path.basename(filepath)
        
        # end performtask()
    
    # end class GenerateMaxStatsTask

# ------------------------- class GenerateMinStatsTask -------------------------
class GenerateMinStatsTask(ProcessRFileTask):
    """
    generate statistics
    """
    
    # ......................... constants .........................
    suffix = "-stats-min.txt"
    outputlocation = "statistics"
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write statistics
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        # there's only one output file path for r file processing:
        filepath = self.outputfilepaths()[0]
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateMinStats!"
            
            td = self.readdata()
            
            reports.report1min(td, filepath)
            return "wrote %s" % os.path.basename(filepath)
        else:
            return "skipping %s...output exists" % os.path.basename(filepath)
        
        # end performtask()
    
    # end class GenerateMinStatsTask

# ------------------------- class GenerateStatsPerObjectTask -------------------------
class GenerateStatsPerObjectTask(ProcessRFileTask):
    """
    generate statistics per object
    """
    
    # ......................... constants .........................
    suffix = "-stats-perobj.txt"
    outputlocation = "statistics"
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write statistics
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        # there's only one output file path for r file processing:
        filepath = self.outputfilepaths()[0]
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateStatsPerObject!"
            
            td = self.readdata()
            
            reports.report1po(td, filepath)
            return "wrote %s" % os.path.basename(filepath)
        else:
            return "skipping %s...output exists" % os.path.basename(filepath)
        
        # end performtask()
    
    # end class GenerateStatsPerObjectTask

# ------------------------- class GenerateStatsPerObjectPerRunTask -------------------------
class GenerateStatsPerObjectPerRunTask(NormalizedPerRunReader, ProcessRFileTask):
    """
    generate statistics per object, per run (= per datestamp), but only
    for controls
    """
    
    # ......................... constants .........................
    suffix = "-stats-perobj.txt"
    outputlocation = "statistics"
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        return list of output this task will write
        """
        
        # this task is unusual in that it only produces output for
        #   some .r files--those belonging to controls; so
        #   pick a scalar that we know exists, read the header, and
        #   determine whether to output files or not (and what their
        #   filenames will be)
        
        # primary scalar data must exist
        scalar = const.primaryscalardata[self.tracker][0]
        sourcefilename = "%s.%s.r" % (self.stem, scalar)
        sourcedir = util.getdirectory("source", self.branch,
                self.tracker, self.geneeff, self.stimno)
        sourcepath = os.path.join(sourcedir, sourcefilename)
        td = trackdata.readheader(sourcepath)
        
        # look in choreography output to see which datestamps we need
        #   to output for; this way we don't need to read the actual
        #   data file
        
        filename = "%s.%s%s" % (self.stem, self.scalar, self.suffix)
        baseoutputdir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        
        return [os.path.join(baseoutputdir, datestamp, filename)
            for datestamp in td.getchoreographydatestamps()]
        
        # end outputfilepaths()
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write statistics
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateStatsPerObjectPerRun!"
            
            filename = "%s.%s%s" % (self.stem, self.scalar, self.suffix)
            baseoutputdir = util.getdirectory(self.outputlocation, self.branch, 
                self.tracker, self.geneeff, self.stimno)
            
            # we're going to be looping over datestamps:
            td = self.readdata()
            messagelist = []
            for datestamp in td.getdatestamps():
                filepath = os.path.join(baseoutputdir, datestamp, filename)
                tdreduced = td.gettrackdatabydatestamp(datestamp)
                reports.report1po(tdreduced, filepath)
                messagelist.append("wrote %s" % filepath)
            return '\n'.join(messagelist)
        else:
            return "skipping per-object, per-run stats...all output exists"
        
        # end performtask()
    
    # end class GenerateStatsPerObjectPerRunTask

# ------------------------- class GenerateTablesTask -------------------------
class GenerateTablesTask(ProcessRFileTask):
    """
    generates tabular output (ie, tab-delimited data files); we output
    raw data separated into intervals as for statistics
    """
    # ......................... constants .........................
    suffixpattern = "-raw-%s-%s.txt"
    outputlocation = "tables"
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        unlike most r-file processing, this task has multiple output
        files; return the list of filenames in a list
        
        NOTE: this method must be kept in careful sync with
        the "util.getintervallist()" function!  the list of filenames
        below must be in one-to-one correspondance with the
        list of intervals returned for "basic stats"
        
        input: none
        output: list of file paths the task will write
        """
        
        # the construction of these filenames is a bit messy; by the
        #   end of it, the end of the filename will look like (for example):
        #   (long filename).area-raw-on-st.txt,
        # where "-on" and "-st" are mnemonics for the start and end
        #   times for the data window (which match the "basic stats" window)
        filenamepattern = "%s.%s%s" % (self.stem, self.scalar, self.suffixpattern)
        outputdir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        filepathpattern = os.path.join(outputdir, filenamepattern)
        
        # one of these for each stimulus variation; since we want to
        #   generate these filenames to be able to check for their 
        #   presence before data is read, must get stim details from
        #   the filename:
        sourcefilename = "%s.%s.r" % (self.stem, self.scalar)
        sourcefilenamedata = trackdata.parsefilename(sourcefilename)
        
        if sourcefilenamedata.nstimuli == 1:
            # single stimulus: pulsed or not?
            stimulus = sourcefilenamedata.stimuli[0]
            stimcount = stimulus.count
            if stimcount == 1:
                # not pulsed
                return [
                    filepathpattern % ("st", "on"),
                    filepathpattern % ("on1", "on4"),
                    filepathpattern % ("on", "onhalf"),
                    filepathpattern % ("onhalf", "off"),
                    filepathpattern % ("off", "offhalf"),
                    filepathpattern % ("offhalf", "end"),
                    ]
            else:
                # pulsed (yes, this is a mess, and duplicated other code...)
                stimduration = stimulus.duration
                stimonset = stimulus.onset
                
                outputlist = [filepathpattern % ("st", "on")]
                for i in range(stimcount):
                    currentonset = stimonset + i * stimulus.onsetdelta
                    outputlist.append(filepathpattern % ("on%d" % i, "on%d+3" % i))
                    outputlist.append(filepathpattern % ("on%d" % i, "off%d" % i))
                    outputlist.append(filepathpattern % ("on%d" % i, "on%dhalf" % i))
                    outputlist.append(filepathpattern % ("on%dhalf" % i, "off%d" % i))
                outputlist.append(filepathpattern % ("off", "end"))
                
                return outputlist
        else:
            # multiple, can't handle this, but the error shouldn't 
            #   be coming from here; therefore return a dummy filename
            #   which shouldn't exist so that later routines will
            #   attempt processing but log errors (I'm not sure 
            #   this is a great idea, but I'm going with it for now)
            filepath = util.getdirectory(self.outputlocation, self.branch, 
                self.tracker, self.geneeff, self.stimno)
            filepath = os.path.join(filepath, "dummyfilenamethatshouldnotexist.txt")
            return [filepath]
        
        # end outputfilepaths()
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write out a table
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateTables!"
            
            td = self.readdata()
            
            nstimuli = td.getnstimuli()
            if nstimuli == 0:
                return "no stimulus...skipping tables"
            if nstimuli > 1:
                return "warning: found %d stimulus columns...skipping tables" % nstimuli
            
            # grab interval list and filename list; verify they 
            #   are the same length:
            intervallist = util.getintervallist(td, "basic stats")
            filepathlist = self.outputfilepaths()
            if len(intervallist) != len(filepathlist):
                return "error: interval list is different length than filename list"
            
            
            # loop over intervals
            outputlist = []
            for interval, filepath in zip(intervallist, filepathlist):
                t1, t2 = interval
                
                tempdata = td.getdatawindowedraw(t1, t2)
                util.savetabseparated(tempdata, filepath, NaNstring="NaN")
                
                outputlist.append("wrote %s" % os.path.basename(filepath))
            return "\n".join(outputlist)
        else:
            return "skipping tables...output exists"
        
        # end performtask()
    
    # end class GenerateTablesTask

