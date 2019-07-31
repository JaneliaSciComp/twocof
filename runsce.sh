#!/bin/bash

# this is a wrapper script for running python programs in the SCE environment

# all the input arguments are passed to Python, so the first one ought to
#   be a python script!

# generalized from a set of less general scripts; by forcing the caller
#   to provide the script path, you (a) allow the caller to specify the
#   correct absolute path (which varies for dev/prod, and local/hudson),
#   (b) put more intelligence in the calling script than in what should
#   be a thin wrapper, and (c) replace three scripts with one

# djo, 12/11


# set up the environment
source /usr/local/SCE/SCE/build/Modules-3.2.6/Modules/3.2.6/init/bash
module use /usr/local/SCE/SCE/build/COTS/
module load cse-build
module load cse-tools
module load cse/matplotlib

python $*



