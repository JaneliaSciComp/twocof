"""

this is the __init__ file for the tasks subpackage

note that this init imports the contents of all the modules, so
you can just "import tasks" and get everything


djo, 12/10

"""

# get all the tasks into this namespace:

from common import *

from analysis import *

from rfileprocessing import *

from secondarydata import *


# aggregate the tasks into the groups that process1 needs; 
#   with the arrival of Larval Olympiad guest projects, 
#   it's necessary to make the list of tasks dependent on tracker

rfilegenerationtasks = {}
rfileprocessingtasks = {}
analysistasks = {}

# ---- Marta's list
# this is also the master list, in that all tasks appear in this
#   list, commented out if they are off

baserfilegenerationtasks = [
    GenerateAbsCastTask,
    # GenerateAnglesTask,
    # GenerateAccelerationTask,
    # GenerateDeltaxyTask,
    GenerateNormalizedTask,
    GenerateNormalizedPerRunTask,
    # GenerateKinksTask,
    # GenerateVelBiasTask,
    ]

baserfileprocessingtasks = [
    # stats, line fits:
    # GenerateStatsTask,
    GenerateStatsPerObjectTask,
    GenerateStatsPerObjectPerRunTask,
    # GenerateMaxStatsTask,
    # GenerateMinStatsTask,
    # GenerateLinearFitsTask,
    
    # plots:
    GenerateAveErrPlotTask,
    GenerateAveErrPerRunPlotTask,
    # GenerateAllDataPlotTask,
    # GenerateBinnedAvePlotTask,
    # GenerateBinnedAvePerObjPlotTask,
    # GenerateBinnedAvePerObjPerRunPlotTask,
    # GenerateAveFitsPlotTask,
    # GenerateScatterFitsPlotTask,
    
    # histograms:
    # GenerateHistogramsTask,
    
    # tabular data:
    # GenerateTablesTask,
    ]

baseanalysistasks = [
    # for Tihana's scratch tests:
    # AnalyzeXYMotionTask,
    
    # multiple plots on one chart:
    GenerateAveErrMulti1PlotTask,
    
    # periodic analysis:
    # PeriodicAnalysis,
    ]

# most of Marta's trackers are the same:
for tracker in ["t2", "t3", "t4", "t5", "t6", "t14", "t15", "t150","t159"]:
    rfilegenerationtasks[tracker] = baserfilegenerationtasks
    rfileprocessingtasks[tracker] = baserfileprocessingtasks
    analysistasks[tracker] = baseanalysistasks

# t1, t111 are a bit different; no multiplot:
for tracker in ["t1", "t111"]:
    rfilegenerationtasks[tracker] = baserfilegenerationtasks
    rfileprocessingtasks[tracker] = baserfileprocessingtasks
    analysistasks[tracker] = []


# larval guests, t7-t12; seem to be the usual with no analysis:
for tracker in ["t7", "t8", "t9", "t10", "t11", "t12"]:
    rfilegenerationtasks[tracker] = baserfilegenerationtasks
    rfileprocessingtasks[tracker] = baserfileprocessingtasks
    analysistasks[tracker] = []



# my tests follow the same-named trackers:
for tracker in ["t1", "t2", "t3"]:
    rfilegenerationtasks["%s-don-test" % tracker] = rfilegenerationtasks[tracker]
    rfileprocessingtasks["%s-don-test" % tracker] = rfileprocessingtasks[tracker]
    analysistasks["%s-don-test" % tracker] = analysistasks[tracker]
    


