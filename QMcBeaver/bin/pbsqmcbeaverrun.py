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

if len(sys.argv) < 2: 
	print "pbsqmcbeaverrun.py <filename> <number of nodes> <processors per node> <queue name>"
	sys.exit(0)

filename   = sys.argv[1]
nodes      = sys.argv[2]
ppn        = sys.argv[3]
queue      = sys.argv[4]

processors = string.atoi(nodes)*string.atoi(ppn)

file = open(filename[:len(filename)-4]+'run','w')

#file.write("#!/bin/bash    \n")
file.write("#PBS -l nodes=" + nodes + ":ppn=" + ppn + " \n")
file.write("#PBS -me        \n")
file.write("#PBS -mb        \n")
file.write("#PBS -q " + queue + " \n")
file.write("\n")
file.write("#!/bin/bash     \n")
file.write("cd " + os.getcwd() + "\n")
file.write("echo " + filename + "\n")

if processors == 1:
	file.write("QMcBeaver " + filename + " >& " + filename[:-4] + \
		   "output \n")
else:
	file.write("mpirun -np " + processors + \
		   " -machinefile $PBS_NODEFILE pQMcBeaver " + filename + \
		   " >& " + filename[:-4] + "output \n")	
file.close()

print "Submitting " + filename
command = "qsub " + filename[:len(filename)-4]+'run'
print command
os.system("qsub " + filename[:len(filename)-4]+'run')

