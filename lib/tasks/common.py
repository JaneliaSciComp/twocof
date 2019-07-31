"""

this file contains the top task superclass and other stuff
useful to all the tasks (eg, mixins)


djo, 12/10


"""


# ------------------------- imports -------------------------
# std lib
import commands
import os
import traceback


# inside the lib
# import constants as const



# ------------------------- class ProcessingTask -------------------------
class ProcessingTask(object):
    """
    base class for our processing tasks
    """
    # ......................... __init__ .........................
    def __init__(self):
        """
        
        input:
        """
        
        raise NotImplementedError("subclass must implement method")
        
        # end __init__()
    
    # ......................... checkcreateoutputfolders() .........................
    def checkcreateoutputfolders(self):
        """
        create all output folders if they don't exist
        
        input: none
        output: boolean: successful or not?
        """
        
        folderset = set(os.path.dirname(filepath) for filepath 
            in self.outputfilepaths())
        for folder in folderset:
            if not os.path.exists(folder):
                # this call will recursively create all 
                #   intermediate directories:
                try:
                    os.makedirs(folder)
                except:
                    return False
            # not sure if this check is necessary, if I'm catching
            #   the errors from above:
            if not os.path.exists(folder):
                return False 
        return True
        
        # end checkcreateoutputfolders()
    
    # ......................... gettaskname() .........................
    def gettaskname(self):
        """
        return the name of the task
        """
        
        # at a minimum, it's the class name:
        return self.__class__.__name__
        
        # end gettaskname()
    
    # ......................... gettasktarget() .........................
    def gettasktarget(self):
        """
        return the target of the task (eg, the tracker/geneeff/stimno)
        """
        
        return "target not defined"
        
        # end gettasktarget()
    
    # ......................... outputexists() .........................
    def outputexists(self):
        """
        does the output of this task exists or not?
        
        input: none
        output: boolean
        """
        
        return all(os.path.exists(filename) for filename 
            in self.outputfilepaths()) 
        
        # end outputexists()
    
    # ......................... outputfilepaths() .........................
    def outputfilepaths(self):
        """
        return the output file paths this tasks will generate
        
        input: none
        output: list of file paths (absolute)
        """
        
        raise NotImplementedError("subclass must implement method")
        
        # end outputfilepaths()
    
    # ......................... performtask() .........................
    def performtask(self):
        """
        execute the task; raise exceptions if there are problems;
        they will be caught and logged
        
        input: none
        output: results string, allowed to contains multiple lines separated
                by \n; this string is parsed for status purposes:
                - string should begin with "error" if an error occurred
                - string should begin with "warning" if a warning should be noted
                - should be one \n separated line per file written, each
                    line beginning with "wrote"
        """
        
        raise NotImplementedError("subclass must implement method")
        
        
        
        # typical sequence:
        
        # read the data
        
        # examine stimulus types etc.; take appropriate action for
        #   those we don't handle (error or skip)
        
        # perform the task
        
        # return a string describing what you did (for logging
        #   and reporting purposes); eg, "wrote xxxx-stats.txt"
        # should begin with "error" if an error occurred
        
        
        # end performtask()
    
    # end class ProcessingTask

# ------------------------- class ExecutableRunner -------------------------
class ExecutableRunner(object):
    """
    this class is actually a mixin; a task should inherit from this class
    as well as its primary superclass if it intends to run an external
    command-line program; this class provides some functionality to
    support that
    """
    
    # ......................... runcommand() .........................
    def runcommand(self, command):
        """
        run the given command
        
        input: string containing a command to be run
        output: string containing the stdout (stderr?) results from the command
        """
        
        # note: this is very bare-bones right now, but I wanted to
        #   abstract it out in case we wanted to exert more 
        #   control in the future (eg, capturing separate stdout and
        #   stderr streams, killing the job after a timeout, etc.)
        
        try:
            ans = commands.getoutput(command)
        except Exception:
            return "error: external command caused an exception; details: %s" % traceback.format_exc()
        
        return ans
        
        # end runcommand()
    
    # end class ExecutableRunner


