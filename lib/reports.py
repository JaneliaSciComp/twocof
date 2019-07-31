"""

this file contains functions that write reports to files based on
trackdata

migrated from methods in rfileprocessing tasks; originally they were
functions in the trackdata file


djo, 3/11

"""



# ------------------------- imports -------------------------
# std lib


# inside the lib
from lib import util


# ------------------------- report1a() -------------------------
def report1a(td, filename):
    """
    tab-delimited file with ave, std dev, count, etc. for 
    a set of data windows
    
    smart enough to detect control vs stimulus trials
    
    input: data and output filename
    output: none (writes file)
    """
    
    intervallist = util.getintervallist(td, "basic stats")
    
    f = open(filename, 'wt')
    f.write("time1\ttime2\tave\tstddev\tNpoints\tNtracks\tsumsquares\n")
    for t1, t2 in intervallist:
        a, s, npts, ntracks, sqrs = td.findaveragewindowedraw(t1, t2)
        f.write("%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.2f\n" % (t1, t2, a, s, npts, ntracks, sqrs))
    f.close()
    
    # end report1a()

# ------------------------- report1lf() -------------------------
def report1lf(td, filename):
    """
    do some linear fits; output range slope intercept corr coeff
    
    input: data and filename
    output: none (writes file)
    """
    
    if not td.hasstimulus():
        raise ValueError("input data must have a stimulus!")
    
    intervallist = util.getintervallist(td, "line fits")
    
    f = open(filename, 'wt')
    f.write("time1\ttime2\tslope\tintercept\tr2\tcount\n")
    for t1, t2 in intervallist:
        m, b, r2, count = td.findlineraw(t1, t2)
        f.write("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n" % (t1, t2, m, b, r2, count))
    f.close()
    
    # end report1lf()

# ------------------------- report1max() -------------------------
def report1max(td, filename):
    """
    statistics on set of minimum values per track in a time window
    
    input: data and output filename
    output: none (writes file)
    """
    
    intervallist = util.getintervallist(td, "unified rules 1")
    
    f = open(filename, 'wt')
    f.write("time1\ttime2\tave\tstddev\tNtracks\tsquaredsum\tsumsquares\n")
    for t1, t2 in intervallist:
        a, s, ntracks, sumsqrs = td.findaverageminmax(t1, t2, 'max')
        sqrsum = (a * ntracks) ** 2
        f.write("%.2f\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t%.2f\n" % 
            (t1, t2, a, s, ntracks, sqrsum, sumsqrs))
    f.close()
    
    # end report1max()

# ------------------------- report1min() -------------------------
def report1min(td, filename):
    """
    statistics on set of minimum values per track in a time window
    
    input: data and output filename
    output: none (writes file)
    """
    
    intervallist = util.getintervallist(td, "unified rules 1")
    
    f = open(filename, 'wt')
    f.write("time1\ttime2\tave\tstddev\tNtracks\tsquaredsum\tsumsquares\n")
    for t1, t2 in intervallist:
        a, s, ntracks, sumsqrs = td.findaverageminmax(t1, t2, 'min')
        sqrsum = (a * ntracks) ** 2
        f.write("%.2f\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t%.2f\n" % 
            (t1, t2, a, s, ntracks, sqrsum, sumsqrs))
    f.close()
    
    # end report1min()

# ------------------------- report1po() -------------------------
def report1po(td, filename):
    """
    tab-delimited file with ave, std dev, count, etc. for 
    a set of data windows; this is for per-object stats
    
    input: data and output filename
    output: none (writes file)
    """
    
    intervallist = util.getintervallist(td, "unified rules 1")
    
    f = open(filename, 'wt')
    f.write("time1\ttime2\tave\tstddev\tNtracks\tsquaredsum\tsumsquares\n")
    for t1, t2 in intervallist:
        a, s, ntracks, sumsqrs = td.findaverageperobject(t1, t2)
        # calculate the squared sum of x, which Marta would like 
        #   output for convenience (likewise for min/max):
        sqrsum = (a * ntracks) ** 2
        f.write("%.2f\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t%.2f\n" % 
            (t1, t2, a, s, ntracks, sqrsum, sumsqrs))
    f.close()
    
    # end report1po()


