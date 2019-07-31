"""

this file contains routines that generate secondary .r files from
the primary .r files


djo, 12/10

"""




# ------------------------- imports -------------------------
# std lib
import os

# numerical and graphical
import numpy

# inside the lib
import common
from lib import constants as const
from lib import plots
from lib import spinedata
from lib import trackdata
from lib import util


# ------------------------- class CreateSecondaryDataTask -------------------------
class CreateSecondaryDataTask(common.ProcessingTask):
    """
    create derived data (typically r-files) that will need further processing
    """
    # ......................... constants .........................
    outputlocation = "derived"
    
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
    
    # ......................... getscalars() .........................
    @staticmethod
    def getscalars(tracker):
        """
        return a list of what scalars you generate; eg, "vel_x_bias";
        this will be used by process1.py to schedule the .r files
        for those scalars to be processed (without having to maintain
        a hard-coded list, as was done previously)
        
        NOTE: therefore, if for some reason you do *not* want a scalar 
            to be processed, don't return it from this method!
        
        NOTE: this is a static method, to be called on the class
        
        NOTE: many if not most subclasses will probably ignore "tracker"
        
        input: tracker
        output: list of scalars (strings)
        """
        
        raise NotImplementedError("subclass must provide!")
        
        # end getscalars()
    
    # ......................... gettasktarget() .........................
    def gettasktarget(self):
        """
        return target
        """
        
        return "%s -- %s -- %s" % (self.tracker, self.geneeff, self.stimno)
        
        # end gettasktarget()
    
    # end class CreateSecondaryDataTask
    
# ------------------------- class GenerateAbsCastTask -------------------------
class GenerateAbsCastTask(CreateSecondaryDataTask):
    """
    generate abs(cast) from cast primary data
    """
    
    # file suffix (also the scalar name):
    
    suffix = "abscast"
    scalar = suffix
    
    # ......................... getscalars() .........................
    @staticmethod
    def getscalars(tracker):
        """
        returns list of scalars we generate
        """
        
        return [GenerateAbsCastTask.suffix]
        
        # end getscalars()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        we have one
        """
        
        # semi-hard coded:
        basedir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        return [
            os.path.join(basedir, "%s.%s.r" % (self.stem, self.scalar)), 
            ]
        
        # end outputfilepaths()
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write the abscast files
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateAbsCastTask!"
            
            sourcedir = util.getdirectory("source", self.branch, 
                self.tracker, self.geneeff, self.stimno)
            sourcefilename = "%s.%s.r" % (self.stem, "cast")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            td = self.trackdatacache.gettrackdata(sourcepath)
            outputdir = util.getdirectory(self.outputlocation, self.branch, 
                self.tracker, self.geneeff, self.stimno)                    
            outputfilepath = os.path.join(outputdir, 
                "%s.%s.r" % (self.stem, self.scalar))
            writer = util.rFileWriter(outputfilepath, NaNstring="NaN")
            
            # this is simple: we want the absolute value of the cast variable:
            for source in td.sources:
                # normalize each block as we go, then write out:
                data = td.getdatablockbysource(source)
                data[:, td.valuecolumn] = abs(data[:, td.valuecolumn])
                
                # source needs to be transformed for each variable:
                writer.writesourceline(source.replace(".scalar.", ".%s." % self.scalar))
                writer.writearray(data)
                
            
            writer.close()
            
            return "wrote abscast file"
        else:
            return "skipping abscast generation...all output exists"
        
        # end performtask()
    
    # end class GenerateAbsCastTask

# ------------------------- class GenerateAccelerationTask -------------------------
class GenerateAccelerationTask(CreateSecondaryDataTask):
    """
    generate the new .r files for second derivative of x, y (acceleration);
    we report the absolute value of the acceleration
    """
    
    # file suffixes (also the scalar names):
    
    accelxsuffix = "abs_accel_x"
    accelysuffix = "abs_accel_y"
    
    # ......................... getscalars() .........................
    @staticmethod
    def getscalars(tracker):
        """
        returns list of scalars we generate
        """
        
        return [GenerateAccelerationTask.accelxsuffix, GenerateAccelerationTask.accelysuffix]
        
        # end getscalars()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        we have two...
        """
        
        # semi-hard coded:
        basedir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        return [
            os.path.join(basedir, "%s.%s.r" % (self.stem, self.accelxsuffix)), 
            os.path.join(basedir, "%s.%s.r" % (self.stem, self.accelysuffix)), 
            ]
        
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
                return "error: could not create output folders for GenerateAcceleration!"
            
            # we need two of the .r files:
            sourcedir = util.getdirectory("source", self.branch, 
                self.tracker, self.geneeff, self.stimno)
            
            sourcefilename = "%s.%s.r" % (self.stem, "x")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            tdx = self.trackdatacache.gettrackdata(sourcepath)        
            
            sourcefilename = "%s.%s.r" % (self.stem, "y")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            tdy = self.trackdatacache.gettrackdata(sourcepath)        
            
            outputdir = util.getdirectory(self.outputlocation, self.branch,
                self.tracker, self.geneeff, self.stimno)
            outputfilepathx = os.path.join(outputdir, "%s.%s.r" % (self.stem, self.accelxsuffix))
            outputfilepathy = os.path.join(outputdir, "%s.%s.r" % (self.stem, self.accelysuffix))
            
            # gotta make sure both are written:
            if self.overwrite or not os.path.exists(outputfilepathx):
                writerx = util.rFileWriter(outputfilepathx, NaNstring="NaN")
            else:
                writerx = util.rFileWriterDummy(outputfilepathx, NaNstring="NaN")
            if self.overwrite or not os.path.exists(outputfilepathy):
                writery = util.rFileWriter(outputfilepathy, NaNstring="NaN")
            else:
                writery = util.rFileWriterDummy(outputfilepathy, NaNstring="NaN")
            
            for source in tdx.sources:
                # since we need to adjust stimuli, use resampled data:
                tx, x = tdx.gettrackbysource(source)
                ty, y = tdy.gettrackbysource(source)
                
                thalfx = 0.5 * (numpy.roll(tx, -1) + tx)
                dtx = numpy.roll(tx, -1) - tx
                thalfy = 0.5 * (numpy.roll(ty, -1) + ty)
                dty = numpy.roll(ty, -1) - ty
                dx = numpy.roll(x, -1) - x
                dy = numpy.roll(y, -1) - y
                
                # rolling the array for use as a time offset makes the
                #   last point wrong, so drop it:
                thalfx = thalfx[:-1]
                thalfy = thalfy[:-1]
                dtx = dtx[:-1]
                dty = dty[:-1]
                dx = dx[:-1]
                dy = dy[:-1]
                
                # up to here, the calculation is the same as for
                #   the Deltaxy tasks, which is basically the velocity
                #   (since we have vel_x/y, I don't do that anymore)
                # just repeat the procedure to do another derivative:
                vx = dx / dtx
                vy = dy / dty
                
                thalfx = 0.5 * (numpy.roll(thalfx, -1) + thalfx)
                dtvx = numpy.roll(thalfx, -1) - thalfx
                thalfy = 0.5 * (numpy.roll(thalfy, -1) + thalfy)
                dtvy = numpy.roll(thalfy, -1) - thalfy
                dvx = numpy.roll(vx, -1) - vx
                dvy = numpy.roll(vy, -1) - vy
                
                # we want the absolute value:
                tempx = numpy.zeros((len(thalfx), 6), dtype=numpy.float)
                tempx[:, 0] = thalfx
                tempx[:, 1] = abs(dvx / dtvx)
                stimx = util.adjuststimuli(tdx, source, thalfx)
                for i in range(4):
                    tempx[:, 2 + i] = stimx[i]
                
                tempy = numpy.zeros((len(thalfy), 6), dtype=numpy.float)
                tempy[:, 0] = thalfy
                tempy[:, 1] = abs(dvy / dtvy)
                stimy = util.adjuststimuli(tdy, source, thalfy)
                for i in range(4):
                    tempy[:, 2 + i] = stimy[i]
                
                # drop the last time point again:
                tempx = tempx[:-1]
                tempy = tempy[:-1]
                
                # source needs to be transformed for each variable:
                writerx.writesourceline(source.replace(".scalar.", ".%s." % self.accelxsuffix))
                writerx.writearray(tempx)
                
                writery.writesourceline(source.replace(".scalar.", ".%s." % self.accelysuffix))
                writery.writearray(tempy)
            
            # close files, write messages
            writerx.close()
            writery.close()
            
            # at least one of these is written or we'd have skipped 
            #   this whole task:
            messagelist = []
            if isinstance(writerx, util.rFileWriter): 
                messagelist.append("wrote %s" % outputfilepathx)
            if isinstance(writery, util.rFileWriter): 
                messagelist.append("wrote %s" % outputfilepathy)
            return '\n'.join(messagelist)
        else:
            return "skipping acceleration generation...all output exists"
        
        # end performtask()
    
    # end class GenerateAccelerationTask

# ------------------------- class GenerateAnglesTask -------------------------
class GenerateAnglesTask(CreateSecondaryDataTask):
    """
    create .r files for angles and cosines
    """
    
    # file names (actually the scalar names):
    anglesuffix = "angle"
    cossuffix = "costheta12"
    theta12suffix = "theta12"
    
    # ......................... findcosinetheta() .........................
    def findcosinetheta(self, t, x, y, offset):
        """
        find the cosine (and other things)
        
        input:  three vectors, time, x, y (resampled, uniform spacing); 
                offset = integer = how many units to displace for the 
                    time between points
        output: t, dx, dy, dr, and cos theta
        """
        
        thalf = 0.5 * (t + numpy.roll(t, -offset))
        dx = numpy.roll(x, -offset) - x
        dy = numpy.roll(y, -offset) - y
        dr = numpy.hypot(dx, dy)
        cos = dx / dr
        
        # chop off the parts at the end that are not properly rolled:
        return (thalf[:-offset], dx[:-offset], dy[:-offset], 
            dr[:-offset], cos[:-offset])
        
        # end findcosinetheta()
    
    # ......................... findcostheta12() .........................
    def findcostheta12(self, t, dx, dy, offset):
        """
        find the cosine of the angle between two successive dr vectors    
        
        input: time, dx, dy; offset
        output: time, costheta12
        """
        
        thalf = 0.5 * (t + numpy.roll(t, -offset))
        dr = numpy.hypot(dx, dy)
        
        dx2 = numpy.roll(dx, -offset)
        dy2 = numpy.roll(dy, -offset)
        dr2 = numpy.hypot(dx2, dy2)
        
        # basically the dot product, normalized:
        costheta12 = (dx * dx2 + dy * dy2) / (dr * dr2)
        
        # chop off the parts at the end that are not properly rolled:
        return thalf[:-offset], costheta12[:-offset]
        
        # end findcostheta12()
    
    # ......................... finddtheta() .........................
    def finddtheta(self, t, theta, offset):
        """
        find dtheta from one time point to the next
        
        input: t, theta; offset
        output: t, dtheta
        """
        
        thalf = 0.5 * (t + numpy.roll(t, -offset))
        dtheta = numpy.roll(theta, -offset) - theta
        
        # chop off the parts at the end that are not properly rolled:
        return thalf[:-offset], dtheta[:-offset]
        
        # end finddtheta()
    
    # ......................... getscalars() .........................
    @staticmethod
    def getscalars(tracker):
        """
        returns list of scalars we generate
        """
        
        return [GenerateAnglesTask.anglesuffix, GenerateAnglesTask.cossuffix, 
            GenerateAnglesTask.theta12suffix]
        
        # end getscalars()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        output files
        """
        
        # semi-hard coded:
        basedir = util.getdirectory(self.outputlocation, self.branch,
            self.tracker, self.geneeff, self.stimno)
        return [
            os.path.join(basedir, "%s.%s.r" % (self.stem, self.anglesuffix)), 
            os.path.join(basedir, "%s.%s.r" % (self.stem, self.cossuffix)), 
            os.path.join(basedir, "%s.%s.r" % (self.stem, self.theta12suffix)), 
            ]
        
        # end outputfilepaths()
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write the angle files
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateAngles!"
            
            # I don't think I have any stimuli restrictions; I'm just
            #   copying the data from the original file, however many
            #   stimuli are marked there
            
            
            # we need two of the .r files:
            sourcedir = util.getdirectory("source", self.branch, 
                self.tracker, self.geneeff, self.stimno)
            
            sourcefilename = "%s.%s.r" % (self.stem, "x")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            tdx = self.trackdatacache.gettrackdata(sourcepath)        
            
            sourcefilename = "%s.%s.r" % (self.stem, "y")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            tdy = self.trackdatacache.gettrackdata(sourcepath)        
            
            # we need to correlate tracks, so if they don't have sources,
            #   we can't continue:
            if tdx.sources is None:
                return "error: x file has no sources!"
            if tdy.sources is None:
                return "error: y file has no sources!"
            
            outputdir = util.getdirectory(self.outputlocation, self.branch,
                self.tracker, self.geneeff, self.stimno)
            outputfilepathangle = os.path.join(outputdir, "%s.%s.r" % (self.stem, self.anglesuffix))
            outputfilepathcos = os.path.join(outputdir, "%s.%s.r" % (self.stem, self.cossuffix))
            outputfilepaththeta12 = os.path.join(outputdir, "%s.%s.r" % (self.stem, self.theta12suffix))
            
            if self.overwrite or not os.path.exists(outputfilepathangle):
                writerangle = util.rFileWriter(outputfilepathangle, NaNstring="NaN")
            else:
                writerangle = util.rFileWriterDummy(outputfilepathangle, NaNstring="NaN")
            if self.overwrite or not os.path.exists(outputfilepathcos):
                writercos = util.rFileWriter(outputfilepathcos, NaNstring="NaN")
            else:
                writercos = util.rFileWriterDummy(outputfilepathcos, NaNstring="NaN")
            if self.overwrite or not os.path.exists(outputfilepaththeta12):
                writertheta12 = util.rFileWriter(outputfilepaththeta12, NaNstring="NaN")
            else:
                writertheta12 = util.rFileWriterDummy(outputfilepaththeta12, NaNstring="NaN")
            
            # find time points that are at least 0.5s apart; since the 
            #   most common dt = 0.05s, I expect the points will 
            #   usually be 10 slots apart
            
            t, x = tdx.gettrackbysource(tdx.sources[0])
            actualdt = t[1] - t[0]
            offset = int(numpy.ceil(0.5 / actualdt))
            
            # we don't store this quantity anymore, so calculate it:
            stride = tdx.getdatastride()
            
            for source in tdx.sources:
                # use resampled data
                t, x = tdx.gettrackbysource(source)
                t, y = tdy.gettrackbysource(source)
                
                thalf_costheta, dx, dy, dr, costheta = self.findcosinetheta(t, 
                    x, y, offset)
                
                # the absolute angle between successive time points is what 
                #   we really want
                theta = numpy.arccos(costheta)
                thalf_dtheta, dtheta = self.finddtheta(thalf_costheta, theta, 1) 
                
                angletracks = [0] * stride
                angletracks[0] = thalf_dtheta
                angletracks[1] = abs(dtheta * 180 / numpy.pi)
                
                # insert stimuli at appropriate places and write:
                angletracks[2:] = util.adjuststimuli(tdx, source, thalf_dtheta)
                writerangle.writesourceline(source.replace(".scalar.", ".angle."))
                writerangle.writearray(util.rectangularize([angletracks]))
                
                # cos theta 12: do as for angles
                thalf_costheta12, costheta12 = self.findcostheta12(thalf_costheta, 
                    dx, dy, 1)
                costheta12tracks = [0] * stride
                costheta12tracks[0] = thalf_costheta12
                costheta12tracks[1] = costheta12
                
                costheta12tracks[2:] = util.adjuststimuli(tdx, source, thalf_costheta12)
                writercos.writesourceline(source.replace(".scalar.", ".costheta12."))
                writercos.writearray(util.rectangularize([costheta12tracks]))
                
                # theta12 = arccos of previous quantity:
                theta12 = numpy.arccos(costheta12)
                theta12tracks = [0] * stride 
                theta12tracks[0] = thalf_costheta12
                theta12tracks[1] = theta12 * 180 / numpy.pi
                
                theta12tracks[2:] = util.adjuststimuli(tdx, source, thalf_costheta12)
                writertheta12.writesourceline(source.replace(".scalar.", ".theta12."))
                writertheta12.writearray(util.rectangularize([theta12tracks]))
                
            
            # close files, write messages
            writerangle.close()
            writercos.close()
            writertheta12.close()
            
            # at least one of these is written or we'd have skipped 
            #   this whole task:
            messagelist = []
            if isinstance(writerangle, util.rFileWriter): 
                messagelist.append("wrote %s" % outputfilepathangle)
            if isinstance(writercos, util.rFileWriter): 
                messagelist.append("wrote %s" % outputfilepathcos)
            if isinstance(writertheta12, util.rFileWriter): 
                messagelist.append("wrote %s" % outputfilepaththeta12)
            return '\n'.join(messagelist)
        else:
            return "skipping angle files generation...all output exists"
        
        # end performtask()
    
    # end class GenerateAnglesTask

# ------------------------- class GenerateDeltaxyTask -------------------------
class GenerateDeltaxyTask(CreateSecondaryDataTask):
    """
    generate the new .r files for the delta x, y (= derivatives);
    we take the absolute value
    """
    
    # file suffixes (also the scalar names):
    
    deltaxsuffix = "abs_vel_x"
    deltaysuffix = "abs_vel_y"
    
    # ......................... getscalars() .........................
    @staticmethod
    def getscalars(tracker):
        """
        returns list of scalars we generate
        """
        
        return [GenerateDeltaxyTask.deltaxsuffix, GenerateDeltaxyTask.deltaysuffix]
        
        # end getscalars()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        we have two...
        """
        
        # semi-hard coded:
        basedir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        return [
            os.path.join(basedir, "%s.%s.r" % (self.stem, self.deltaxsuffix)), 
            os.path.join(basedir, "%s.%s.r" % (self.stem, self.deltaysuffix)), 
            ]
        
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
                return "error: could not create output folders for GenerateDeltaxy!"
            
            # we need two of the .r files:
            sourcedir = util.getdirectory("source", self.branch, 
                self.tracker, self.geneeff, self.stimno)
            sourcefilename = "%s.%s.r" % (self.stem, "x")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            tdx = self.trackdatacache.gettrackdata(sourcepath)        
            
            sourcefilename = "%s.%s.r" % (self.stem, "y")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            tdy = self.trackdatacache.gettrackdata(sourcepath)        
            
            outputdir = util.getdirectory(self.outputlocation, self.branch,
                self.tracker, self.geneeff, self.stimno)
            outputfilepathx = os.path.join(outputdir, "%s.%s.r" % (self.stem, self.deltaxsuffix))
            outputfilepathy = os.path.join(outputdir, "%s.%s.r" % (self.stem, self.deltaysuffix))
            
            # not sure if this is all necessary for delta x/y, but I just
            #   modified the vel_bias task
            if self.overwrite or not os.path.exists(outputfilepathx):
                writerx = util.rFileWriter(outputfilepathx, NaNstring="NaN")
            else:
                writerx = util.rFileWriterDummy(outputfilepathx, NaNstring="NaN")
            if self.overwrite or not os.path.exists(outputfilepathy):
                writery = util.rFileWriter(outputfilepathy, NaNstring="NaN")
            else:
                writery = util.rFileWriterDummy(outputfilepathy, NaNstring="NaN")
            
            for source in tdx.sources:
                # since we need to adjust stimuli, use resampled data:
                tx, x = tdx.gettrackbysource(source)
                ty, y = tdy.gettrackbysource(source)
                
                thalfx = 0.5 * (numpy.roll(tx, -1) + tx)
                dtx = numpy.roll(tx, -1) - tx
                thalfy = 0.5 * (numpy.roll(ty, -1) + ty)
                dty = numpy.roll(ty, -1) - ty
                dx = numpy.roll(x, -1) - x
                dy = numpy.roll(y, -1) - y
                
                # rolling the array for use as a time offset makes the
                #   last point wrong, so drop it:
                thalfx = thalfx[:-1]
                thalfy = thalfy[:-1]
                dtx = dtx[:-1]
                dty = dty[:-1]
                dx = dx[:-1]
                dy = dy[:-1]
                
                # recall we want the absolute value:
                tempx = numpy.zeros((len(thalfx), 6), dtype=numpy.float)
                tempx[:, 0] = thalfx
                tempx[:, 1] = abs(dx / dtx)
                stimx = util.adjuststimuli(tdx, source, thalfx)
                for i in range(4):
                    tempx[:, 2 + i] = stimx[i]
                
                tempy = numpy.zeros((len(thalfy), 6), dtype=numpy.float)
                tempy[:, 0] = thalfy
                tempy[:, 1] = abs(dy / dty)
                stimy = util.adjuststimuli(tdy, source, thalfy)
                for i in range(4):
                    tempy[:, 2 + i] = stimy[i]
                
                
                # source needs to be transformed for each variable:
                writerx.writesourceline(source.replace(".scalar.", ".%s." % self.deltaxsuffix))
                writerx.writearray(tempx)
                
                writery.writesourceline(source.replace(".scalar.", ".%s." % self.deltaysuffix))
                writery.writearray(tempy)
            
            # close files, write messages
            writerx.close()
            writery.close()
            
            # at least one of these is written or we'd have skipped 
            #   this whole task:
            messagelist = []
            if isinstance(writerx, util.rFileWriter): 
                messagelist.append("wrote %s" % outputfilepathx)
            if isinstance(writery, util.rFileWriter): 
                messagelist.append("wrote %s" % outputfilepathy)
            return '\n'.join(messagelist)
        else:
            return "skipping delta x, y generation...all output exists"
        
        # end performtask()
    
    # end class GenerateDeltaxyTask

# ------------------------- class GenerateKinksTask -------------------------
class GenerateKinksTask(CreateSecondaryDataTask):
    """
    generate the new .r files for the angles in the spine file
    """
    
    # pattern for scalars:
    # (double check that this will parse correctly!)
    scalarpattern = "kink-%d_%dto%d_%d"
    
    # list of angles to calculate; numbers refer to the index into
    #   the list of nodes on the spine
    
    '''
    # working list (1-offset) (partial, testing)
    anglelist = [
        ((1, 2), (2, 3)),
        ((1, 3), (3, 4)),
        ]
    
    '''
    # final list:
    anglelist = [
        ((1, 2), (2, 3)), 
        ((1, 3), (3, 4)),
        ((1, 4), (4, 5)),
        ((1, 5), (5, 6)),
        ((1, 6), (6, 7)),
        ((1, 7), (7, 8)),
        ((1, 8), (8, 9)),
        ((1, 9), (9, 10)),
        ((1, 10), (10, 11)),
        ((1, 2), (5, 7)),
        ((1, 3), (5, 7)),
        ((1, 4), (5, 7)),
        ((1, 5), (5, 7)),
        ((1, 2), (9, 11)), 
        ((1, 3), (9, 11)),
        ((1, 4), (9, 11)),
        ((1, 5), (9, 11)),
        ((1, 6), (9, 11)),
        ((1, 7), (9, 11)),
        ((1, 8), (9, 11)),
        ((1, 9), (9, 11)),
        ((1, 2), (8, 11)),
        ((1, 3), (8, 11)),
        ((1, 4), (8, 11)),
        ((1, 5), (8, 11)),
        ((1, 6), (8, 11)),
        ((1, 7), (8, 11)),
        ((1, 8), (8, 11)),
        # added 5-oct-10:
        ((2, 3), (3, 4)),
        ((3, 4), (4, 5)),
        ((4, 5), (5, 6)),
        ((5, 6), (6, 7)),
        ((2, 4), (4, 5)),
        ((3, 5), (5, 6)),
        ((4, 6), (6, 7)),
        ((5, 7), (7, 8)),
        ((6, 8), (8, 9)),
        ((7, 9), (9, 10)),
        ((8, 10), (10, 11)),
        ((2, 5), (5, 6)),
        ((3, 6), (6, 7)),
        ((4, 7), (7, 8)),
        ((5, 8), (8, 9)),
        ((6, 9), (9, 10)),
        ((7, 10), (10, 11)),
        ]
    
    # ......................... getscalars() .........................
    @staticmethod
    def getscalars(tracker):
        """
        returns list of scalars we generate
        """
        
        # unlike other tasks, we'll generate this programmatically
        
        return [GenerateKinksTask.scalarpattern % (a, b, c, d) 
            for (a, b), (c, d) in GenerateKinksTask.anglelist]
        
        # end getscalars()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        we have far too many
        """
        
        # generate from scalar list:
        basedir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        return [os.path.join(basedir, "%s.%s.r" % (self.stem, scalar))
            for scalar in GenerateKinksTask.getscalars(tracker)]
        
        # end outputfilepaths()
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write the kink files
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateKinks!"
            
            # we need the spine file:
            sourcedir = util.getdirectory("source", self.branch,
                self.tracker, self.geneeff, self.stimno) 
            sourcefilename = "%s.%s.r" % (self.stem, "spine")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            sd = spinedata.readspinefile(sourcepath)
            
            # we also need a trackdata file to grab the stimuli from;
            #   let's use the "kink" one (seems appropriate)
            sourcefilename2 = "%s.%s.r" % (self.stem, "kink")
            sourcepath2 = os.path.join(sourcedir, sourcefilename2)
            td = trackdata.readfile(sourcepath2)
            
            
            # the spine data object will get us the angles, we just
            #   have to correlate with the stimuli; and we have
            #   to do it over and over and over for each kink...
            
            messagelist = []
            basedir = util.getdirectory(self.outputlocation, self.branch,
                self.tracker, self.geneeff, self.stimno)
            for (a1, a2), (b1, b2) in self.anglelist:
                scalar = self.scalarpattern % (a1, a2, b1, b2) 
                kinkfilename = os.path.join("%s.%s.r" % (self.stem, scalar))
                kinkfilepath = os.path.join(basedir, kinkfilename) 
                if self.overwrite or not os.path.exists(kinkfilepath):
                    writer = util.rFileWriter(kinkfilepath, NaNstring="NaN")
                else:
                    writer = util.rFileWriterDummy(kinkfilepath, NaNstring="NaN")
                
                # angledata = repeating (t, angle) columns, one per track: 
                for source in sd.sources:
                    
                    # some .spine files have more sources than the corresponding
                    #   .kink file; if we don't have stimuli for a source in 
                    #   the .kink file, ignore that source
                    if source not in td.sources:
                        continue
                    
                    # we actually want the greater of the angles for the
                    #   specified nodes, and for their "complement"--the
                    #   same pattern of nodes indexed from the other
                    #   end of the animal; segment 1-->2 becomes 11-->10,
                    #   and so on (the tracker can't distinguish head from
                    #   tail, so we calculate both and report biggest)
                    #   (yes, this is what Marta wants; only the head angle
                    #   gets big, so if angle is small, doesn't matter if
                    #   it's from the head or tail)
                    # note the switch to zero-offset
                    angledata1 = sd.anglebetweenbysource(source, (a1 - 1, a2 - 1), (b1 - 1, b2 - 1))
                    c1 = spinedata.nspinenodes + 1 - a1
                    c2 = spinedata.nspinenodes + 1 - a2
                    d1 = spinedata.nspinenodes + 1 - b1
                    d2 = spinedata.nspinenodes + 1 - b2
                    angledata2 = sd.anglebetweenbysource(source, (c1 - 1, c2 - 1), (d1 - 1, d2 - 1))
                    # we want the absolute value, and the larger of the 
                    #   two, per time point:
                    angledata1[1] = abs(angledata1[1])
                    angledata2[1] = abs(angledata2[1])
                    angledata = [angledata1[0], numpy.where(angledata1[1] > angledata2[1], 
                        angledata1[1], angledata2[1]) * 180 / numpy.pi]
                    
                    # stimulus adjustment stuff here...tack the stimuli
                    #   columns onto the list, then convert to array
                    
                    # adjust stimuli; note that we can only adjust stims on
                    #    a uniform time grid (to save our sanity); thus,
                    #   resample the kink angles; use the same time grid
                    #   as for the trackdata
                    t = angledata[0]
                    tmin = numpy.nanmin(t)
                    tmax = numpy.nanmax(t)
                    dt = td.getsuggesteddt()
                    nsteps = (tmax - tmin) / dt + 1
                    tgrid = numpy.arange(nsteps) * dt + tmin
                    
                    anglesgrid = numpy.interp(tgrid, t, angledata[1],
                        right=numpy.NaN, left=numpy.NaN)
                    
                    angledata = [tgrid, anglesgrid]
                    angledata.extend(util.adjuststimuli(td, source, tgrid))
                    anglestowrite = numpy.array(angledata).transpose()
                    
                    writer.writesourceline(source.replace(".scalar.", ".%s." % scalar))
                    writer.writearray(anglestowrite)
                    
                writer.close()
                messagelist.append("wrote %s" % kinkfilepath)
            return '\n'.join(messagelist)
        else:
            return "skipping kinks generation...all output exists"
        
        # end performtask()
    
    # end class GenerateKinksTask

# ------------------------- class GenerateNormalizedTask -------------------------
class GenerateNormalizedTask(CreateSecondaryDataTask):
    """
    generate the new .r files for many scalars; normalize the 
    input by some initial data window
    """
    
    scalarpattern = "%snorm"
    
    # ......................... getscalars() .........................
    @staticmethod
    def getscalars(tracker):
        """
        returns list of scalars we generate
        """
        
        # generate them on the fly:
        return [GenerateNormalizedTask.scalarpattern % item 
            for item in const.normalizedscalars[tracker]]
        
        # end getscalars()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        
        basedir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        return [os.path.join(basedir, "%s.%s.r" % (self.stem, scalar)) 
            for scalar in GenerateNormalizedTask.getscalars(self.tracker)]
        
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
                return "error: could not create output folders for GenerateNormalized!"
            
            
            # loop over scalars; read data, identify interval,
            #   find mean in interval, normalize data, write file
            nfiles = 0
            for scalar in const.normalizedscalars[self.tracker]:
                sourcedir = util.getdirectory("source", self.branch, 
                    self.tracker, self.geneeff, self.stimno)
                sourcefilename = "%s.%s.r" % (self.stem, scalar)
                sourcepath = os.path.join(sourcedir, sourcefilename)
                td = self.trackdatacache.gettrackdata(sourcepath)
                outputdir = util.getdirectory(self.outputlocation, self.branch, 
                    self.tracker, self.geneeff, self.stimno)                    
                outputfilepath = os.path.join(outputdir, 
                    "%s.%s.r" % (self.stem, self.scalarpattern % scalar))
                writer = util.rFileWriter(outputfilepath, NaNstring="NaN")
                
                # identify baseline interval; this is arbitrary/temporary,
                #   and I'll later fix it to match what Marta wants
                # find the ave (ignore the other stats), and watch for 
                #   means near zero (can't really normalize those?)
                '''
                if td.getnstimuli() > 0:
                    stimulus = min(td.filenamedata.stimuli,
                        key=lambda x: getattr(x, "onset"))
                    stimonset = stimulus.onset
                    t1 = const.ignoretime
                    t2 = stimonset - 5
                else:
                    t1 = const.ignoretime
                    t2 =  25
                '''
                # hard code for now:
                t1 = 15.0
                t2 = 29.0
                
                ave, std, npts, ntracks, sumsq = td.findaveragewindowedraw(t1, t2)
                if ave < 0.1:
                    # warning here?
                    pass
                
                
                for source in td.sources:
                    # normalize each block as we go, then write out:
                    data = td.getdatablockbysource(source)
                    data[:, td.valuecolumn] /= ave
                    
                    # need a clip here?  didn't help--still have
                    #   one plot exceeding allowed size!
                    # numpy.clip(data[:, td.valuecolumn], -10, 10, data[:, td.valuecolumn])
                    
                    # source needs to be transformed for each variable:
                    writer.writesourceline(source.replace(".scalar.", ".%s." % (self.scalarpattern % scalar)))
                    writer.writearray(data)
                    
                
                writer.close()
                nfiles += 1
            
            return "wrote %d GenerateNormalized files" % nfiles
        else:
            return "skipping normalized generation...all output exists"
        
        # end performtask()
    
    # end class GenerateNormalizedTask

# ------------------------- class GenerateNormalizedPerRunTask -------------------------
class GenerateNormalizedPerRunTask(CreateSecondaryDataTask):
    """
    generate the new .r files for many scalars; normalize the 
    input by some initial data window, but do so on a per-run
    (per datestamp) basis
    
    (compare GenerateNormalizedTask, which normalizes everything)
    """
    
    scalarpattern = "%snormperrun"
    
    # ......................... getscalars() .........................
    @staticmethod
    def getscalars(tracker):
        """
        returns list of scalars we generate; in this case, we will
        choose in other tasks whether we want norm or normperrun,
        so here we lie and say we don't produce any scalars
        """
        
        return [ ]
        
        # end getscalars()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        
        basedir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        return [os.path.join(basedir, "%s.%s.r" % (self.stem, self.scalarpattern % scalar)) 
            for scalar in const.normalizedscalars[self.tracker]]
        
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
                return "error: could not create output folders for GenerateNormalized!"
            
            
            # loop over scalars; read data, identify interval,
            #   find mean in interval, normalize data, write file
            nfiles = 0
            for scalar in const.normalizedscalars[self.tracker]:
                sourcedir = util.getdirectory("source", self.branch, 
                    self.tracker, self.geneeff, self.stimno)
                sourcefilename = "%s.%s.r" % (self.stem, scalar)
                sourcepath = os.path.join(sourcedir, sourcefilename)
                td = self.trackdatacache.gettrackdata(sourcepath)
                outputdir = util.getdirectory(self.outputlocation, self.branch, 
                    self.tracker, self.geneeff, self.stimno)                    
                outputfilepath = os.path.join(outputdir, 
                    "%s.%s.r" % (self.stem, self.scalarpattern % scalar))
                writer = util.rFileWriter(outputfilepath, NaNstring="NaN")
                
                # identify baseline interval; hard coded for now:
                t1 = 15.0
                t2 = 29.0
                
                for datestamp in td.getdatestamps():
                    tdreduced = td.gettrackdatabydatestamp(datestamp)
                    
                    ave, std, npts, ntracks, sumsq = tdreduced.findaveragewindowedraw(t1, t2)
                    if ave < 0.1:
                        # warning here?
                        pass
                    
                    for source in tdreduced.sources:
                        # normalize each block as we go, then write out:
                        data = tdreduced.getdatablockbysource(source)
                        data[:, tdreduced.valuecolumn] /= ave
                        
                        # source needs to be transformed for each variable:
                        writer.writesourceline(source.replace(".scalar.", ".%s." % (self.scalarpattern % scalar)))
                        writer.writearray(data)
                
                writer.close()
                nfiles += 1
            
            return "wrote %d GenerateNormalizedPerRun files" % nfiles
        else:
            return "skipping normalized per run generation...all output exists"
        
        # end performtask()
    
    # end class GenerateNormalizedPerRunTask

# ------------------------- class GenerateVelBiasTask -------------------------
class GenerateVelBiasTask(CreateSecondaryDataTask):
    """
    generate the new .r files for the velocity biases
    """
    
    # file suffixes (also the scalar names):
    
    velxsuffix = "vel_x_bias"
    velysuffix = "vel_y_bias"
    velxysuffix = "vel_xy_bias"
    
    # ......................... getscalars() .........................
    @staticmethod
    def getscalars(tracker):
        """
        returns list of scalars we generate
        """
        
        return [GenerateVelBiasTask.velxsuffix, GenerateVelBiasTask.velysuffix, 
            GenerateVelBiasTask.velxysuffix]
        
        # end getscalars()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        we have three...
        """
        
        # semi-hard coded:
        basedir = util.getdirectory(self.outputlocation, self.branch, 
            self.tracker, self.geneeff, self.stimno)
        return [
            os.path.join(basedir, "%s.%s.r" % (self.stem, self.velxsuffix)), 
            os.path.join(basedir, "%s.%s.r" % (self.stem, self.velysuffix)), 
            os.path.join(basedir, "%s.%s.r" % (self.stem, self.velxysuffix)), 
            ]
        
        # end outputfilepaths()
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        write the vel_bias files
        
        input: none
        output: results string; should begin with "error" if an
                error occurred
        """
        
        if self.overwrite or not self.outputexists():
            # this branch should always be taken in the process I
            #   have planned, but I don't want to rely on that
            
            if not self.checkcreateoutputfolders():
                return "error: could not create output folders for GenerateVelBias!"
            
            # we need two of the .r files:
            sourcedir = util.getdirectory("source", self.branch, 
                self.tracker, self.geneeff, self.stimno)
            
            sourcefilename = "%s.%s.r" % (self.stem, "vel_x")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            tdvx = self.trackdatacache.gettrackdata(sourcepath)        
            
            sourcefilename = "%s.%s.r" % (self.stem, "vel_y")
            sourcepath = os.path.join(sourcedir, sourcefilename)
            tdvy = self.trackdatacache.gettrackdata(sourcepath)        
            
            # we need to match up corresponding tracks from vx and vy 
            #   and calculate track-by-track, outputting as we go;
            #   (a) both datasets must have sources, and (b) use
            # dummy objects to simplify the output process:
            if tdvx.sources is None:
                return "error: vel_x file has no sources!"
            if tdvy.sources is None:
                return "error: vel_y file has no sources!"
            outputdir = util.getdirectory(self.outputlocation, self.branch,
                self.tracker, self.geneeff, self.stimno) 
                
            outputfilepathx = os.path.join(outputdir, "%s.%s.r" % (self.stem, self.velxsuffix))
            outputfilepathy = os.path.join(outputdir, "%s.%s.r" % (self.stem, self.velysuffix))
            outputfilepathxy = os.path.join(outputdir, "%s.%s.r" % (self.stem, self.velxysuffix))
            
            if self.overwrite or not os.path.exists(outputfilepathx):
                writerx = util.rFileWriter(outputfilepathx, NaNstring="NaN")
            else:
                writerx = util.rFileWriterDummy(outputfilepathx, NaNstring="NaN")
            if self.overwrite or not os.path.exists(outputfilepathy):
                writery = util.rFileWriter(outputfilepathy, NaNstring="NaN")
            else:
                writery = util.rFileWriterDummy(outputfilepathy, NaNstring="NaN")
            if self.overwrite or not os.path.exists(outputfilepathxy):
                writerxy = util.rFileWriter(outputfilepathxy, NaNstring="NaN")
            else:
                writerxy = util.rFileWriterDummy(outputfilepathxy, NaNstring="NaN")
            
            valuecolumn = tdvx.valuecolumn
            for source in tdvx.sources:
                datax = tdvx.getdatablockbysource(source)
                datay = tdvy.getdatablockbysource(source)
                
                vx = datax[:, valuecolumn]
                vy = datay[:, valuecolumn]
                v = numpy.sqrt(vx ** 2 + vy ** 2)
                
                xbias = datax.copy()
                xbias[:, valuecolumn] /= v
                
                ybias = datay.copy()
                ybias[:, valuecolumn] /= v
                
                xybias = datax.copy()
                xybias[:, valuecolumn] = (abs(vx) - abs(vy)) / v
                
                # source needs to be transformed for each variable:
                writerx.writesourceline(source.replace(".scalar.", ".vel_x_bias."))
                writerx.writearray(xbias)
                
                writery.writesourceline(source.replace(".scalar.", ".vel_y_bias."))
                writery.writearray(xbias)
                
                writerxy.writesourceline(source.replace(".scalar.", ".vel_xy_bias."))
                writerxy.writearray(xbias)
            
            # close files, write messages
            writerx.close()
            writery.close()
            writerxy.close()
            
            # at least one of these is written or we'd have skipped 
            #   this whole task:
            messagelist = []
            if isinstance(writerx, util.rFileWriter): 
                messagelist.append("wrote %s" % outputfilepathx)
            if isinstance(writery, util.rFileWriter): 
                messagelist.append("wrote %s" % outputfilepathy)
            if isinstance(writerxy, util.rFileWriter): 
                messagelist.append("wrote %s" % outputfilepathxy)
            return '\n'.join(messagelist)
        else:
            return "skipping vel_bias generation...all output exists"
        
        # end performtask()
    
    # end class GenerateVelBiasTask

