#!/usr/bin/env python

#            QMcBeaver
#
#         Constructed by 
#
#     Michael Todd Feldmann 
#              and 
#   David Randall "Chip" Kent IV
#
# Copyright 2000.  All rights reserved.
#
# drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

import sys
import string
import os

# FLAGS
# r - run

if len(sys.argv) < 4: 
	print "llnl_run.py <filename> <number of nodes> <time limit in h> <flags>"
	sys.exit(0)

filename   = sys.argv[1]
nodes      = sys.argv[2]
processors = str(string.atoi(nodes)*4)
time       = sys.argv[3]
flags      = sys.argv[4]

file = open(filename[:len(filename)-4]+'run','w')

file.write("#!/bin/csh         # Sets your shell\n")
file.write("#PSUB -s /bin/csh  # Sets your shell in batch\n")
file.write("#PSUB -c 100Mb     # Select a 100Mb memory limit\n")
file.write("#PSUB -ln " + nodes + "        # Number of nodes you want to use\n")
file.write("#PSUB -g " + processors + "         # Number of tasks to use, 4 tasks per node\n")
file.write("#PSUB  -eo         # Send std error & std out to the same file\n")
file.write("#PSUB -tM " + time + "h       # Select your time limit. The default time limit\n")
file.write("                   # is only 30 minutes! Time can be HH:MM:SS or HH:MM\n")
file.write("#PSUB -b caltech\n")
file.write("#PSUB -me          # mail at end of calc\n")
file.write("#PSUB -mb          # mail at begining of calc\n")
file.write("\n")
file.write("cd " + os.getcwd() + "\n")
file.write("poe /g/g20/drk3548/code/QMcBeaver/bin/pQMcBeaver.llnl " + filename + "\n")
file.close()

if string.find(flags,'r') != -1:
	print "Submitting " + filename
	os.system("psub " + filename[:len(filename)-4]+'run')

