"""

this file contains the plotting routines, both the specific routines
for specific plots and all the supporting routines that
(eg) place labels, format axes, etc.


(this stuff used to be in trackdata.py, a holdover from the days
when the whole project was one monolithic file; moving plot
routines to this file is part of the modularization effort)


djo, 12/10

"""


# ------------------------- imports -------------------------
# std lib
import os
import sys


# numerical and graphical
import numpy
import matplotlib

# backend choice (if specified) must precede import of pyplot;
#   check for SGE environment variable to detect if we're on a node
#if "SGE_O_HOST" in os.environ:
# we're on a node

matplotlib.use('Agg')

#else:
#    # use the default, probably TkAgg
#    pass

import matplotlib.pyplot as plt


# inside the lib
import constants as const
import trackdata
import util


# ------------------------- plotalldata() -------------------------
def plotalldata(td, filename=None, show=True):
    """
    plots all the raw data as a scatter plot 
    
    input: track data object, flags 
    output: figure
    """
    
    fig = plt.figure()
    axes = fig.add_subplot(111)
    setfirstaxislabel(td, axes)
    setupgridticks(axes)
    setstimulusmarkers(td, axes)
    
    
    # this may not be maximally efficient...
    tylist = [td.gettrackbysource(source, resampled=False) for source in td.sources]
    tlist, ylist = zip(*tylist)
    t = numpy.hstack(tlist)
    y = numpy.hstack(ylist)
    
    # print "n data points =", len(t)
    
    # note that '.' = small disk, ',' = pixel
    axes.plot(t, y, ',')
    
    
    if show:
        plt.show()
    
    if filename is not None:
        fig.savefig(filename)
    
    return fig
    
    # end plotalldata()

# ------------------------- plotavehistogram() -------------------------
def plotavehistogram(td, tmin, tmax, filename=None, show=True):
    """
    plot histogram of mean values per track in a time window;
    use relative frequencies; also returns histogram data
    
    input:  track data object; time interval; filename to save; show flag
    output: figure, numbers, bins
    """
    
    fig = plt.figure()
    axes1 = fig.add_subplot(111)
    sethistaxislabel(td, axes1, relative=True)
    
    # turn on grid, but don't adjust ticks
    axes1.grid(True, which='major', linewidth=0.5, linestyle='-', color=const.minorgridcolor)
    axes1.set_axisbelow(True)
    
    avelist = []
    
    # loop over the tracks and see if the data is in the window;
    #   if so, grab the average
    for sourcename in td.sources:
        t, value = td.gettrackbysource(sourcename, resampled=False)
        data = value[numpy.logical_and(t > tmin, t < tmax)]
        
        if data.any():
            avelist.append(numpy.average(data))
    
    if avelist:
        # try to tame the range:
        scalarrange = const.getscalarscale(td.filenamedata.scalar)
        n, bins, patches = axes1.hist(avelist, normed=True, range=scalarrange)
    else:
        # no data!
        n, bins, patches = [], [], []
    
    if show:
        plt.show()
    
    if filename is not None:
        fig.savefig(filename)
    
    return fig, n, bins
    
    # end plotavehistogram()

# ------------------------- plotaverage() -------------------------
def plotaverage(td, filename=None, show=True):
    """
    plots the average, no error, no count, no "stride"
    
    input: track data object, flags 
    output: figure
    """
    
    fig = plt.figure()
    axes1 = fig.add_subplot(111)
    setfirstaxislabel(td, axes1)
    setupgridticks(axes1)
    setstimulusmarkers(td, axes1)
    
    t = td.gettimegrid()
    ave, stddev, count, sumsquares, squaredsum = td.findaveragetrack()
    
    axes1.plot(t, ave)
    
    if show:
        plt.show()
    
    if filename is not None:
        fig.savefig(filename)
    
    return fig
    
    # end plotaverage()

# ------------------------- plotaveragefits() -------------------------
def plotaveragefits(td, filename=None, show=True):
    """
    plots the average, no error, no count, no "stride"
    
    adds the linear regression line
    
    input: track data object, flags 
    output: figure
    """
    
    fig = plt.figure()
    axes1 = fig.add_subplot(111)
    setfirstaxislabel(td, axes1)
    setupgridticks(axes1)
    setstimulusmarkers(td, axes1)
    
    t = td.gettimegrid()
    ave, stddev, count, sumsquares, squaredsum = td.findaveragetrack()
    
    # in order to see the superimposed plots, make these 
    #   lines a lighter gray
    axes1.plot(t, ave, color=const.graylinecolor)
    
    
    # add fit lines; kick up the line width for visibility
    for t1, t2 in util.getintervallist(td, "line fits"):
        times = numpy.array([t1, t2])
        m, b, r2, count = td.findlineraw(t1, t2)
        line = m * times + b
        axes1.plot(times, line, '-', linewidth=3)
    
    setplotsize(ave, td, fig, axes1)
    
    if show:
        plt.show()
    
    if filename is not None:
        fig.savefig(filename)
    
    return fig
    
    # end plotaveragefits()

# ------------------------- plotaverageerror() -------------------------
def plotaverageerror(td, showN=False, filename=None, show=True):
    """
    plot the average track, plus errors
    
    input:  track data object
            flag to also plot N
            output: figure, (arrays of t, ave, std, count, sumsquares, squaredsum)
    """
    
    fig = plt.figure()
    axes1 = fig.add_subplot(111)
    setfirstaxislabel(td, axes1)
    
    t = td.gettimegrid()
    ave, stddev, count, sumsquares, squaredsum = td.findaveragetrack()
    
    # cut back on number of data points
    stride = 10
    axes1.errorbar(t[::stride], ave[::stride], yerr=stddev[::stride], 
        fmt='bo-')
    
    if showN:
        axes2 = axes1.twinx()
        axes2.set_ylabel('count', color='g')
        axes2.plot(t[::stride], count[::stride], 'go-')
        setupgridticks(axes1, axes2)
    else:
        setupgridticks(axes1)
    setstimulusmarkers(td, axes1)
    setplotsize([ave[::stride] - stddev[::stride], ave[::stride] + stddev[::stride]], 
        td, fig, axes1)
    
    if show:
        plt.show()
    
    if filename is not None:
        fig.savefig(filename)
    
    # return copies, not slices
    return fig, (numpy.array(t[::stride]), numpy.array(ave[::stride]),
        numpy.array(stddev[::stride]), numpy.array(count[::stride]),
        numpy.array(sumsquares[::stride]), numpy.array(squaredsum[::stride]),)
    
    # end plotaverageerror()

# ------------------------- plotbinnedaverage() -------------------------
def plotbinnedaverage(td, tmin, tmax, tstep, kind="summed", 
        showN=False, prestimvalues=None, filename=None, show=True):
    """
    plot a binned average; optionally include the prestim value 
    as well
    
    (uses resampled data)
    
    input:  the data; time start/stop/step
            kind = "summed" or "per object"
            flag to show the number of tracks averaged
            optional prestimulus value = (time, ave, std, count)
            filename for save, flag for show
    output: figure, t data, ave data
    """
    
    fig = plt.figure()
    axes1 = fig.add_subplot(111)
    setfirstaxislabel(td, axes1)
    
    # t, ave, stddev, count = trackdata.findbinnedaverage(tmin, tmax, tstep)
    avelist = list(td.findbinnedaverage(tmin, tmax, tstep, kind))
    
    if prestimvalues is not None:
        for i in range(len(avelist)):
            avelist[i] = numpy.concatenate([prestimvalues[i], avelist[i]])
    
    t, ave, stddev, count= avelist
    axes1.errorbar(t, ave, yerr=stddev, fmt='bo-')
    
    
    if showN:
        axes2 = axes1.twinx()
        axes2.set_ylabel('count', color='g')
        axes2.plot(t, count, 'go-')
        setupgridticks(axes1, axes2)
    else:
        setupgridticks(axes1)
    setstimulusmarkers(td, axes1)
    
    # ave and stddev are lists at this point; find the min/max 
    #   ave -/+ stddev for the plot limits:
    plotmin = numpy.nanmin([a - s for a, s in zip(ave, stddev)])
    plotmax = numpy.nanmax([a + s for a, s in zip(ave, stddev)])
    setplotsize([plotmin, plotmax], td, fig, axes1)
    
    if show:
        plt.show()
    
    if filename is not None:
        fig.savefig(filename)
    
    return fig, avelist
    
    # end plotbinnedaverage()

# ------------------------- plotrandomraw() -------------------------
def plotrandomraw(td, npts, linefit=False, 
    filename=None, show=True):
    """
    plot random raw data points
    
    input: trackdata; how many points to plot; flags
    output: figure
    """
    
    fig = plt.figure()
    axes = fig.add_subplot(111)
    setfirstaxislabel(td, axes)
    setupgridticks(axes)
    setstimulusmarkers(td, axes)
    
    t, y = td.getrandompoints(npts)
    
    if linefit:
        # lighter gray if superpimposing line fites:
        axes.plot(t, y, '.', color=const.graylinecolor)
    else:
        axes.plot(t, y, '.')
    
    
    # possible bug (that we possibly rely on): no check for
    #   "linefit" flag here!
    
    # maybe plot some line fits:
    for t1, t2 in util.getintervallist(td, "line fits"):
        times = numpy.array([t1, t2])
        m, b, r2, count = td.findlineraw(t1, t2)
        line = m * times + b
        axes.plot(times, line, '-', linewidth=3)
    
    setplotsize(y, td, fig, axes)
    
    if show:
        plt.show()
    
    if filename is not None:
        fig.savefig(filename)
    
    return fig
    
    
    # end plotrandomraw()

# ------------------------- plottracks() -------------------------
def plottracks(td, stride=1, resampled=True,
    filename=None, style=',', show=True):
    """
    test routine--plots data points with a given style
    
    use stride > 1 to cut down on density of wiggles
    
    input:  trackdata object; flags
            note style = matplotlib style; default ',' is a dot
    output: figure
    """
    
    fig = plt.figure()
    axes = fig.add_subplot(111)
    setfirstaxislabel(td, axes)
    setupgridticks(axes)
    setstimulusmarkers(td, axes)
    
    for sourcename in td.sources:
        t, y = td.gettrackbysource(sourcename, resampled=resampled)
        axes.plot(t[::stride], y[::stride], style)
    
    if show:
        plt.show()
    
    if filename is not None:
        fig.savefig(filename)
    
    return fig
    
    # end plottracks()

# ------------------------- setfirstaxislabel() -------------------------
def setfirstaxislabel(td, axis, title=True):
    """
    set up the label on an axis; despite the name, you can use
    this for the second axis, too, just set title=False if you
    set the title with the first axis
    
    input:  trackdata object; axis object; flag to set title
    output: none (formats axis)
    """
    
    axis.set_xlabel('time (s)')
    if td.filenamedata is not None:
        if title:
            axis.set_title(os.path.basename(td.filenamedata.filename))
        axis.set_ylabel("%s %s" % (td.filenamedata.scalar, 
            const.getscalarunit(td.filenamedata.scalar)))
    else:
        axis.set_ylabel('quantity')
    
    # end setfirstaxislabel()

# ------------------------- sethistaxislabel() -------------------------
def sethistaxislabel(td, axis, relative=False):
    """
    set up the labels for histograms
    
    input: trackdata object; axis object; flag for relative freq.
    output: none (formats axis)
    """
    
    if relative:
        axis.set_ylabel('relative frequency')
    else:
        axis.set_ylabel('count')
    if td.filenamedata is not None:
        axis.set_title(os.path.basename(td.filenamedata.filename))
        axis.set_xlabel("%s %s" % (td.filenamedata.scalar, 
            const.getscalarunit(td.filenamedata.scalar)))
    else:
        axis.set_xlabel('quantity')
    
    # end sethistaxislabel()

# ------------------------- setplotsize() -------------------------
def setplotsize(data, td, fig, axis):
    """
    set the plot to a standard size
    
    input: data to be plotted (any shape), trackdata, figure, axis
    output: none (formats axis)
    """
    
    # first, sanity check: if all data isn't finite, leave at default 
    #   plot size
    # could be any shape:
    testdata = numpy.array(data)
    if numpy.all(numpy.logical_or(numpy.isnan(testdata), numpy.isinf(testdata))):
        return
    
    tmin, tmax = td.findtimerange()
    datamin, datamax = numpy.nanmin(data), numpy.nanmax(data)
    
    scalar = td.filenamedata.scalar
    
    # initial size (inches):
    plotorigx, plotorigy = fig.get_size_inches()
    
    timemin, timemax = const.getscalarscale("t")
    valuemin, valuemax = const.getscalarscale(scalar)
    
    # set the limits; for time, the min is zero regardless
    # set_xbound and set_xlim seem to do the same thing?    
    axis.set_autoscale_on(False)
    if tmax <= timemax:
        # smaller; set the max, keep image size:
        axis.set_xbound(0, timemax)
        plotfinalx = plotorigx
    else:
        # bigger; expand image size
        axis.set_xbound(0, tmax)
        plotfinalx = plotorigx * tmax / float(timemax)
    
    # for the data, need to fiddle with upper and lower limits:
    finalymin = min(datamin, valuemin)
    finalymax = max(datamax, valuemax)
    axis.set_ybound(finalymin, finalymax)
    
    if finalymin != valuemin or finalymax != valuemax:
        plotfinaly = plotorigy * (finalymax - finalymin) / float(valuemax - valuemin)
    else:
        plotfinaly = plotorigy
    
    # final check; if aspect ratio is too large, just
    #   go back to the default size
    if plotfinaly / plotfinalx > const.aspectratiolimit:
        plotfinalx, plotfinaly = plotorigx, plotorigy
    fig.set_size_inches(plotfinalx, plotfinaly)
    
    # end setplotsize()

# ------------------------- setstimulusmarkers() -------------------------
def setstimulusmarkers(td, axis, kind="duration"):
    """
    mark the stimuli on the plot
    
    input: the data and the axis object; kind = "onset" or "duration"
    output: none (formats axis)
    """
    
    # temporarily disabled; not correct for all runs,
    #   and Marta needs plots for a paper
    return
    
    if kind not in ["onset", "duration"]:
        raise ValueError("unknown kind %s" % kind)
    
    nstim = td.getnstimuli()
    if nstim == 1:
        stimulus = td.filenamedata.stimuli[0]
        color = const.stimuluscolorlist[0]
        if not stimulus.pulsed:
            setstimulusmarkerssingle(axis, stimulus, kind=kind, color=color)
        else:
            setstimulusmarkerspulsed(axis, stimulus, kind=kind, color=color)
    elif nstim > 1:
        # multiple stimuli
        '''
        #   try dividing the vertical space, with same color for all
        for i, stimulus in enumerate(td.filenamedata.stimuli):
            # I need floats, so turn nstim into a float:
            nstim = 1.0 * nstim
            ymin = i / nstim
            ymax = (i + 1) / nstim
            if not stimulus.pulsed:
                setstimulusmarkerssingle(axis, stimulus, 
                    ymin=ymin, ymax=ymax, kind=kind)
            else:
                setstimulusmarkerspulsed(axis, stimulus, 
                    ymin=ymin, ymax=ymax, kind=kind)
        '''
        # don't divide vertical space, different colors:
        for i, stimulus in enumerate(td.filenamedata.stimuli):
            color = const.stimuluscolorlist[i]
            if not stimulus.pulsed:
                setstimulusmarkerssingle(axis, stimulus, kind=kind, color=color)
            else:
                setstimulusmarkerspulsed(axis, stimulus, kind=kind, color=color)
    else:
        # no stimulus = no markers
        pass
    
    # end setstimulusmarkers()

# ------------------------- setstimulusmarkerssingle() -------------------------
def setstimulusmarkerssingle(axis, stimulus, ymin=0, ymax=1, kind='duration', color='r'):
    """
    mark stimulus on a plot
    
    input: axis; stimulus; vertical span
    output: none
    """
    
    if kind == "onset":
        axis.axvline(x=stimulus.onset, color=color)
    elif kind == "duration":
        endtime = stimulus.onset + stimulus.duration
        axis.axvspan(stimulus.onset, endtime, ymin=ymin, ymax=ymax, 
            facecolor=color, edgecolor='none', alpha=const.basealpha)    
    
    # end setstimulusmarkerssingle()

# ------------------------- setstimulusmarkerspulsed() -------------------------
def setstimulusmarkerspulsed(axis, stimulus, ymin=0, ymax=1, kind='duration', color='r'):
    """
    set stimulus markers, pulsed stim
    
    input: axis; stimulus; vertical span
    output: none
    """
    
    stimcount = stimulus.count
    stimduration = stimulus.duration
    stimonset = stimulus.onset
    delta = stimulus.onsetdelta
    stimonsetlist = [stimonset + i * delta for i in range(stimcount)]
    for t in stimonsetlist:
        if kind == "onset":
            axis.axvline(x=t, color=color)
        elif kind == "duration":
            # really thin bars need heavier alpha:
            if stimduration < 1.0:
                alpha = const.thinalpha
            else:
                alpha = const.basealpha
            axis.axvspan(t, t + stimduration, ymin=ymin, ymax=ymax, 
                facecolor=color, edgecolor='none', alpha=alpha)
    
    # end setstimulusmarkerspulsed()

# ------------------------- setupgridticks() -------------------------
def setupgridticks(axis1, axis2=None):
    """
    set up the grid the major/minor ticks for one or (optionally)
    two axes
    """
    
    axis1.minorticks_on()
    # put grid on axis1 instead of 2 (makes more sense; can do both
    #   but they superimpose = confusing, ugly); turn on major and 
    #   minor, but make the minor grid fainter
    axis1.grid(True, which='major', linewidth=1, linestyle='-', color=const.majorgridcolor)
    axis1.grid(True, which='minor', linewidth=0.5, linestyle='-', color=const.minorgridcolor)
    
    
    if axis2 is not None:
        # note that x-axis takes its ticks from axis2, not axis1,
        #   if both are present
        axis2.minorticks_on()
    
    # add time ticks on 1s intervals; add to "front" axis, 
    #   which is axis 2 if present;  on the time axis,
    #   shift the minor ticks --> major spacing so that grid is less dense
    xmajor = matplotlib.ticker.MultipleLocator(5)
    xminor = matplotlib.ticker.MultipleLocator(5)
    if axis2 is not None:
        axis2.xaxis.set_major_locator(xmajor)
        axis2.xaxis.set_minor_locator(xminor)
    else:
        axis1.xaxis.set_major_locator(xmajor)
        axis1.xaxis.set_minor_locator(xminor)
    
    # grid should be on bottom:
    axis1.set_axisbelow(True)
    if axis2 is not None:
        axis2.set_axisbelow(True)
    
    # end setupgridticks()

