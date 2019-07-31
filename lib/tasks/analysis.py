"""

this file contains routines for tasks that run on a unit other than
"once for each and every .r file"

(or, in other words, AnalysisTask and its subclasses)


djo, 12/10

"""



# ------------------------- imports -------------------------
# std lib
import os
import tempfile

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
from lib import util


# ------------------------- class AnalysisTask -------------------------
class AnalysisTask(common.ProcessingTask):
    """
    an analysis task is similar to a GenerateSecondaryDataTask in that
    it operates on an entire experiment and not a single r-file, but
    it differs in that it produces nothing that needs further processing;
    it's an end-point step
    """
    # ......................... constants .........................
    outputlocation = "analysis"
    
    # ......................... __init__ .........................
    def __init__(self, branch, tracker, geneeff, stimno, trackdatacache, overwrite=False):
        """
        
        input:  branch: eg, screen
                tracker: eg, t1 or t1-don-text
                gene/effector: eg, w1118@n
                stimulus/animalno: eg, b_100Hz3V_20s1x30s0s#n#n#n@1
                trackdata cache instance
                flag whether to overwrite existing or not
        """
        
        # housekeeping
        self.branch = branch
        self.tracker = tracker
        self.geneeff = geneeff
        self.stimno = stimno
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
        
        return "%s -- %s -- %s" % (self.tracker, self.geneeff, self.stimno)
        
        # end gettasktarget()
    
    # end class AnalysisTask


# ------------------------- class AnalyzeXYMotionTask -------------------------
class AnalyzeXYMotionTask(AnalysisTask):
    """
    analyze tracks to probe whether there is preferential motion along
    a scratch in the plate (which is usually oriented on the y-axis,
    but we do x as well for comparison)
    """
    # ......................... constants .........................
    # tolerance for determining "along axis", in degrees:
    angletolerance = 10.0
    
    # y-axis:
    # if we perform atan2(dy, dx), we want angles that fall in
    #   the range [90 - tol, 90 + tol] and [270 - tol, 270 + tol] degrees;
    #   fortunately, those correspond to the same range in tangent, and
    #   since the tangent is singular along the y-axis, and anti-symmetric
    #   around it, we need only check that the absolute magnitude of the
    #   tangent (= dy / dx) be greater than:
    
    ytangenttolerance = numpy.tan((numpy.pi / 180.) * (90 - angletolerance))
    
    # x-axis:
    # we want angles in [-tol, +tol] and [180 - tol, 180 + tol]; so
    #   in a similar way, we want the abs mag of the tangent to
    #   be less than:
    
    xtangenttolerance = numpy.tan((numpy.pi / 180.) * angletolerance)
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        output file paths
        """
        
        basedir = util.getdirectory("analysis", self.branch, 
            self.tracker, self.geneeff, self.stimno)
        
        self.xsumspath = os.path.join(basedir, "%s.xmotion-sums.txt" % self.stem)
        self.ysumspath = os.path.join(basedir, "%s.ymotion-sums.txt" % self.stem)
        
        return [self.xsumspath, self.ysumspath]
        
        # end outputfilepaths()
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write the files
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for AnalyzeXYMotion!"
            
            # we need two of the .r files:
            sourcedir = util.getdirectory("source", self.branch, 
                self.tracker, self.geneeff, self.stimno)
            
            sourcefilename = "%s.%s.r" % (self.stem, "x")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            tdx = self.trackdatacache.gettrackdata(sourcepath)        
            
            sourcefilename = "%s.%s.r" % (self.stem, "y")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            tdy = self.trackdatacache.gettrackdata(sourcepath)        
            
            # identify and loop over intervals, writing one line per interval
            # do for x and y:
            
            messages = []
            f = open(self.ysumspath, 'wt')
            f.write("t1\tt2\ttime on y\ttotal time\trel time\tpath on y\ttotal path\trel path\tntracks\tave t on y\tstd dev\tsum sqrs\tave path on y\tstd dev\tsum sqrs\n")
            for t1, t2 in util.getintervallist(tdx, "xy-motion"):
                (ytime, totaltime, ypath, totalpath, ntracks, 
                    timeave, timestd, timesumsqrs, 
                    ypathave, ypathstd, ypathsumsqrs) = self.windowedcalc(tdx, tdy, t1, t2, 'y')
                # note that numpy.divide will emit NaN or Inf, appropriate for 
                #   tabular output
                f.write("%.1f\t%.1f\t%.1f\t%.1f\t%.3f\t%.1f\t%.1f\t%.3f\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % 
                    (t1, t2, ytime, totaltime, numpy.divide(ytime, totaltime), 
                    ypath, totalpath, numpy.divide(ypath, totalpath), ntracks, 
                    timeave, timestd, timesumsqrs, ypathave, ypathstd, ypathsumsqrs))
            f.close()
            messages.append('wrote %s' % self.ysumspath)
            
            f = open(self.xsumspath, 'wt')
            f.write("t1\tt2\ttime on x\ttotal time\trel time\tpath on x\ttotal path\trel path\tntracks\tave t on x\tstd dev\tsum sqrs\tave path on x\tstd dev\tsum sqrs\n")
            for t1, t2 in util.getintervallist(tdx, "xy-motion"):
                (ytime, totaltime, ypath, totalpath, ntracks, 
                    timeave, timestd, timesumsqrs, 
                    ypathave, ypathstd, ypathsumsqrs) = self.windowedcalc(tdx, tdy, t1, t2, 'x')
                # note that numpy.divide will emit NaN or Inf, appropriate for 
                #   tabular output
                f.write("%.1f\t%.1f\t%.1f\t%.1f\t%.3f\t%.1f\t%.1f\t%.3f\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % 
                    (t1, t2, ytime, totaltime, numpy.divide(ytime, totaltime), 
                    ypath, totalpath, numpy.divide(ypath, totalpath), ntracks, 
                    timeave, timestd, timesumsqrs, ypathave, ypathstd, ypathsumsqrs))
            f.close()
            messages.append('wrote %s' % self.xsumspath)
            
            return '\n'.join(messages)
            
        else:
            return "skipping xy-motion analysis...all output exists"
        
        # end performtask()
    
    # ......................... windowedcalc() .........................
    def windowedcalc(self, tdx, tdy, t1, t2, axis):
        """
        calculate values for a window
        
        input: beginning and end time; axis = 'y' or 'x'
        output: time on y-axis, total time, y-projection of path, total path,
                ntracks, time on y-axix ave per track, std dev (1/N-1), 
                ypath ave per track, std dev (1/N-1)
        """
        
        # loop over the data and add stuff up!
        
        # NOTE: initially written only for y-axis, so variables all
        #   are named that way; only difference is in the tangent 
        #   test for "close to axis"
        
        # these are overall:
        totaltime = 0.0
        timeonyaxis = 0.0
        totalpath = 0.0
        ypath = 0.0
        
        # these are per source:
        timeratiosum = 0.0
        timeratiosumsquares = 0.0
        ypathratiosum = 0.0
        ypathratiosumsquares = 0.0
        rationumber = 0
        
        
        valuecolumn = tdx.valuecolumn
        for source in tdx.sources:
            datax = tdx.getdatablockbysource(source)
            datay = tdy.getdatablockbysource(source)
            
            t = datax[:, 0]
            x = datax[:, valuecolumn]
            y = datay[:, valuecolumn]
            
            # filter down to window:
            x = x[numpy.logical_and(t > t1, t < t2)]
            y = y[numpy.logical_and(t > t1, t < t2)]
            t = t[numpy.logical_and(t > t1, t < t2)]
            
            if len(t) > 1:
                # need to do these per source
                totaltimesource = 0.0
                timeonyaxissource = 0.0
                totalpathsource = 0.0
                ypathsource = 0.0
                
                # loop over time; we're looking at deltas, so 
                #   roll/subtract/strip off spurious endpoint:
                dt = numpy.roll(t, -1) - t
                dx = numpy.roll(x, -1) - x
                dy = numpy.roll(y, -1) - y
                dt = dt[:-1]
                dx = dx[:-1]
                dy = dy[:-1]
                
                for deltat, deltax, deltay in zip(dt, dx, dy):
                    # check if the motion is along an axis in the 
                    #   time step is using tangent trick (see
                    #   notes in constants for this class)
                    # probably a better way to do it than this
                    #   repeatedly executed if clause...
                    if axis == 'y':
                        if abs(deltay / deltax) > self.ytangenttolerance:
                            timeonyaxis += deltat
                            timeonyaxissource += deltat
                        ypath += abs(deltay)
                        ypathsource += abs(deltay)
                    else:
                        # axis == 'x':
                        if abs(deltay / deltax) < self.xtangenttolerance:
                            timeonyaxis += deltat
                            timeonyaxissource += deltat
                        ypath += abs(deltax)
                        ypathsource += abs(deltax)
                    totaltime += deltat
                    totaltimesource += deltat
                    
                    # sum path and projection:
                    totalpath += numpy.hypot(deltax, deltay)
                    totalpathsource += numpy.hypot(deltax, deltay)
                timeratio = timeonyaxissource / totaltimesource
                timeratiosum += timeratio
                timeratiosumsquares += timeratio * timeratio
                
                ypathratio = ypathsource / totalpathsource
                ypathratiosum += ypathratio
                ypathratiosumsquares += ypathratio * ypathratio
                
                rationumber += 1
        
        if rationumber > 1:
            # need 2 items to do stats:
            timeratioave = timeratiosum / rationumber
            timeratioerr = timeratiosumsquares / rationumber - timeratioave * timeratioave
            timeratioerr = numpy.sqrt(timeratioerr * rationumber / (rationumber - 1))
            
            ypathratioave = ypathratiosum / rationumber
            ypathratioerr = ypathratiosumsquares / rationumber - ypathratioave * ypathratioave
            ypathratioerr = numpy.sqrt(ypathratioerr * rationumber / (rationumber - 1))
        else:
            timeratioave = numpy.NaN
            timeratioerr = numpy.NaN
            
            ypathratioave = numpy.NaN
            ypathratioerr = numpy.NaN
        
        return (timeonyaxis, totaltime, ypath, totalpath, rationumber, 
            timeratioave, timeratioerr, timeratiosumsquares,
            ypathratioave, ypathratioerr, ypathratiosumsquares)
        
        # end windowedcalc()
    
    # end class AnalyzeXYMotionTask

# ------------------------- class GenerateAveErrMulti1PlotTask -------------------------
class GenerateAveErrMulti1PlotTask(AnalysisTask):
    """
    generate a stacked 5-plot
    """
    # ......................... constants .........................
    outputlocation = "plots"
    
    # scalars to plot (with plot format):
    scalardict = {
        "dir": 'ro-',
        "curve": 'bo-',
        "speed085": 'go-',
        "midline": 'co-',
        "area": 'ko-',
        }
    
    # we cut back on the number of points shown:
    stride = 10
    
    # ......................... __init__ .........................
    def __init__(self, branch, tracker, geneeff, stimno, trackdatacache, overwrite=False):
        """
        
        input:  branch: eg, screen
                tracker: eg, t1 or t1-don-text
                gene/effector: eg, w1118@n
                stimulus/animalno: eg, b_100Hz3V_20s1x30s0s#n#n#n@1
                trackdata cache instance
                flag whether to overwrite existing or not
        """
        
        # housekeeping
        self.branch = branch
        self.tracker = tracker
        self.geneeff = geneeff
        self.stimno = stimno
        self.trackdatacache = trackdatacache
        self.overwrite = overwrite
        
        # generate the stem from the above:
        self.stem = "%s@%s@%s" % (self.geneeff, self.tracker, self.stimno)
        
        # filename pattern for input files:
        self.filenamepattern = "%s@%s@%s.%s.r"
        
        # end __init__()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        output file paths
        """
        
        basedir = util.getdirectory(self.outputlocation, self.branch,
            self.tracker, self.geneeff, self.stimno)
        return[os.path.join(basedir, "%s.combo1-ave-err%s" % (self.stem, const.plotextension))]
        
        # end outputfilepaths()
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        do it!
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        # only one file:
        outputfilepath = self.outputfilepaths()[0]
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateAveErrMulti1PlotTask!"
            
            # read data
            td = {}
            for scalar in self.scalardict:
                filename = self.filenamepattern % (self.geneeff, self.tracker, 
                    self.stimno, scalar)
                sourcedir = util.getdirectory("source", self.branch,
                    self.tracker, self.geneeff, self.stimno)
                sourcefilepath = os.path.join(sourcedir, filename)
                td[scalar] = self.trackdatacache.gettrackdata(sourcefilepath)
                
            # not very modularized yet...
            fig = plt.figure()
            
            axis = fig.add_subplot(5, 1, 1)
            axis.set_title(self.stem)
            self.setupgridticks(axis, "area", 0.5, majorlabels=False)
            self.plotscalar(td, "area", axis, elinewidth=0.5, markersize=2)
            
            axis = fig.add_subplot(5, 1, 2)
            self.setupgridticks(axis, "curve", 10, majorlabels=False)
            self.plotscalar(td, "curve", axis, elinewidth=0.5, markersize=2)
            
            axis = fig.add_subplot(5, 1, 3)
            self.setupgridticks(axis, "dir", 0.5, majorlabels=False)
            self.plotscalar(td, "dir", axis, elinewidth=0.5, markersize=2)
            
            axis = fig.add_subplot(5, 1, 4)
            self.setupgridticks(axis, "midline", 0.5, majorlabels=False)
            self.plotscalar(td, "midline", axis, elinewidth=0.5, markersize=2)
            
            axis = fig.add_subplot(5, 1, 5)
            axis.set_xlabel('time (s)')
            self.setupgridticks(axis, "speed085", 0.5, majorlabels=True)
            self.plotscalar(td, "speed085", axis, elinewidth=0.5, markersize=2)
            
            # default (saved) is 800 x 600 pixels at 100 dpi = 8 x 6 inches;
            #   double each length:
            fig.set_size_inches(16., 12.)
            fig.savefig(outputfilepath)
            
            # fig.clf() not enough to release memory!  need to use close()
            plt.close(fig)
            
            return "wrote %s" % os.path.basename(outputfilepath)
        else:
            return "skipping %s...output exists" % os.path.basename(outputfilepath)
        
        # end performtask()
    
    # ......................... plotscalar() .........................
    def plotscalar(self, tddict, scalar, axis, **kw):
        """
        plot a scalar
        
        input: dict of data, scalar, axis; keyword args to pass to plot
        output: none (plots)
        """
        
        t = tddict[scalar].gettimegrid()
        ave, stddev, count, sumsquares, squaredsum = tddict[scalar].findaveragetrack()
        axis.errorbar(t[::self.stride], ave[::self.stride], yerr=stddev[::self.stride], 
            fmt=self.scalardict[scalar], **kw)
        
        # end plotscalar()
    
    # ......................... setupplot() .........................
    def setupgridticks(self, axis1, scalar, ytickinterval, majorlabels=True):
        """
        set up the grid the major/minor ticks the axis
        """
        
        axis1.set_ylabel("%s %s" % (scalar, const.getscalarunit(scalar)))
        
        axis1.minorticks_on()
        # turn on major and minor, but make the minor grid fainter
        axis1.grid(True, which='major', linewidth=1, linestyle='-', color=const.majorgridcolor)
        # axis1.grid(True, which='minor', linewidth=0.5, linestyle='-', color=const.minorgridcolor)
        axis1.grid(True, which='minor', linewidth=2, linestyle='-', color='k')
        
        # Marta wants some specific times with darker lines; make those minor,
        #   since we want the labels at the usual places (major)
        xmajor = matplotlib.ticker.MultipleLocator(5)
        ticklist = [i * 15 for i in range(2, 10)]
        xminor = matplotlib.ticker.FixedLocator(ticklist)
        axis1.xaxis.set_major_locator(xmajor)
        axis1.xaxis.set_minor_locator(xminor)
        
        # get rid of time labels on all plots but bottom?
        if not majorlabels:
            axis1.xaxis.set_major_formatter(matplotlib.ticker.NullFormatter())
        
        # y-axis: 
        ymajor = matplotlib.ticker.MultipleLocator(ytickinterval)
        yminor = matplotlib.ticker.NullLocator()
        axis1.yaxis.set_major_locator(ymajor)
        axis1.yaxis.set_minor_locator(yminor)
        
        # grid should be on bottom:
        axis1.set_axisbelow(True)
        
        # end setupgridticks()
    
    # end class GenerateAveErrMulti1PlotTask

# ------------------------- class PeriodicAnalysis -------------------------
class PeriodicAnalysis(AnalysisTask, common.ExecutableRunner):
    """
    run an external prog to do some periodic analysis
    """
    # ......................... constants .........................
    executablename = "periodogram.py"
    
    
    # ......................... __init__ .........................
    def __init__(self, branch, tracker, geneeff, stimno, trackdatacache, overwrite=False):
        """
        
        input:  branch: eg, screen
                tracker: eg, t1 or t1-don-text
                gene/effector: eg, w1118@n
                stimulus/animalno: eg, b_100Hz3V_20s1x30s0s#n#n#n@1
                trackdata cache instance
                flag whether to overwrite existing or not
        """
        
        # housekeeping
        self.branch = branch
        self.tracker = tracker
        self.geneeff = geneeff
        self.stimno = stimno
        self.trackdatacache = trackdatacache
        self.overwrite = overwrite
        
        # generate the stem from the above:
        self.stem = "%s@%s@%s" % (self.geneeff, self.tracker, self.stimno)
        
        # filename pattern for input files:
        self.filenamepattern = "%s@%s@%s.%s.r"
        
        # end __init__()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        output file paths
        """
        
        basedir = util.getdirectory("analysis", self.branch,
            self.tracker, self.geneeff, self.stimno)
        return[os.path.join(basedir, "%s.periodic.txt" % self.stem)]
        
        # end outputfilepaths()
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        do it!
        
        input: none
        output: results string; should begin with "error" if an
                error occurred        
        """
        
        # only one file:
        outputfilepath = self.outputfilepaths()[0]
        
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for PeriodicAnalysis!"
            
            # generate input filename:
            sourcedir = util.getdirectory("source", self.branch, 
                self.tracker, self.geneeff, self.stimno)
            filename = "%s.speed.r" % self.stem
            sourcefilepath = os.path.join(sourcedir, filename)
            
            td = self.trackdatacache.gettrackdata(sourcefilepath)
            
            
            # need to pass the intervals to the program, so write out
            #   to a temporary file:
            intervalpath = self.writeintervalfile(td)
            
            
            # run the command and return its output; note that the
            #   "self.runcommand" method ought to trap all exceptions,
            #   so no need to try/except here
            execpath = os.path.join(util.getresourcelocation("executable"), 
                self.executablename)
            commandstring = "python %s %s %s %s" % (execpath, inputfilepath, 
                intervalpath, outputfilepath) 
            results = self.runcommand(commandstring)
            
            # clean up
            os.remove(intervalpath)
            
            return results
        
        # end performtask()
    
    # ......................... writeintervalfile() .........................
    def writeintervalfile(self, td):
        """
        write out a file with the list of intervals to process
        
        input: trackdata object (for getting intervals)
        output: the path of the (temporary) interval file
        """
        
        fd, intervalpath = tempfile.mkstemp()
        
        # at this point the file is open, and we can wrap the descriptor
        #   in a file object:
        f = os.fdopen(fd, 'w')
        for t1, t2 in util.getintervallist(td, "per-object stats, absolute"):
            f.write("%s\t%s\n" % (t1, t2))
        f.close()
        
        return intervalpath
        
        # end writeintervalfile()
    
    # end class PeriodicAnalysis

