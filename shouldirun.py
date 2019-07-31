"""

this script determines whether processing should be run or not

prints to standard output: "1" if you should run, "0" if not


expected to be run by shell snippet in Hudson, roughly:

if ((`python shouldirun.py` == "1"))
then
    echo "run processing"
else
    echo "do not run processing"
fi


djo, 2/12

"""


# ------------------------- imports -------------------------

import time


# print 1 if we want to run, 0 if not:
def run(shouldrun=True):
    if shouldrun:
        print 1
    else:
        print 0

def dontrun():
    run(False)


# ------------------------- script starts here -------------------------
if __name__ == "__main__":
    
    
    # ----- biweekly:
    
    # get week of year, with Sunday as first day; convert to int
    w = int(time.strftime("%U"))
    
    # run on even weeks
    if w % 2 == 0:
        run()
    else: 
        dontrun()
    
