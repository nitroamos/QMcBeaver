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
	print "pbsqmcbeaverrun.py <exe> <filename> <number of nodes> <processors per node> <queue name>"
	sys.exit(0)

exe        = sys.argv[1]
filename   = sys.argv[2]   

if len(sys.argv) > 3:
	nodes      = sys.argv[3]
else:
	nodes      = "1"

if len(sys.argv) > 4:
	ppn        = sys.argv[4]
else:
	ppn        = "2"

if len(sys.argv) > 5:
	queue      = sys.argv[5]
else:
	queue      = "default"
				
processors = string.atoi(nodes)*string.atoi(ppn)

file = open(filename[:len(filename)-4]+'run','w')

#file.write("#!/bin/bash    \n")
file.write("#PBS -l nodes=" + nodes + ":ppn=" + ppn + " \n")
file.write("#PBS -N " + filename + " \n")
file.write("#PBS -me        \n")
file.write("#PBS -mb        \n")
if queue != "default":
	file.write("#PBS -q " + queue + " \n")
file.write("\n")
file.write("#!/bin/bash     \n")
file.write("cd " + os.getcwd() + "\n")
file.write("echo " + filename + "\n")

if processors == 1:
	file.write(exe + " " + filename + " >& " + filename[:-4] + \
		   "output \n")
else:
	file.write("mpirun -np " + str(processors) + \
		   " -machinefile $PBS_NODEFILE " + exe + " " + filename + \
		   " >& " + filename[:-4] + "output \n")	
file.close()

print "Submitting " + filename
command = "qsub " + filename[:len(filename)-4]+'run'
print command
os.system("qsub " + filename[:len(filename)-4]+'run')

