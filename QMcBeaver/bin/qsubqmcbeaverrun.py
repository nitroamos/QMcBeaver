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

default_nodes = 1
default_ppn   = 1
default_hours = 96

emailaddress = "your@address.edu"

if len(sys.argv) < 2: 
	print "qsubqmcbeaverrun.py <exe> <filename> " + \
		"[number of nodes="+str(default_nodes)+"] " \
		"[processors per node="+str(default_ppn)+"] " \
		"[hours="+str(default_hours)+"] "
	sys.exit(0)

exe        = sys.argv[1]
filename   = sys.argv[2]   

if len(sys.argv) > 3:
	nodes      = string.atoi(sys.argv[3])
else:
	nodes      = default_nodes

if len(sys.argv) > 4:
	ppn        = string.atoi(sys.argv[4])
else:
	ppn        = default_ppn

if len(sys.argv) > 5:
      hours        = string.atoi(sys.argv[5])
else:
      hours        = default_hours

processors = nodes*ppn

base = filename[:-4]
file = open(base+'run','w')
ckmf_file = open(filename,'r')
ckmf_data = ckmf_file.readlines()
ckmf_file.close()
jobname = base[:-1]

for i in range(len(ckmf_data)):
	if string.find(ckmf_data[i],"temp_dir") == 0:
		temp_dir = ckmf_data[i+1][:-1]
		break

print "Using temp space: " + temp_dir

debug = 0
#debug = 1

file.write("#!/bin/csh    \n")
file.write("#PBS -N " + jobname +        " \n")
file.write("#PBS -l nodes=" + str(nodes) + ":ppn=" + str(ppn) + ",walltime=" + str(hours) + ":20:00\n")
file.write("#PBS -m ea        \n") #b for begin, e for end, a for abort
file.write("#PBS -j oe -k n         \n")# -j oe merge std out and std err, -k n don't keep the file
file.write("\n")
file.write("cd " + os.getcwd() + "\n")

if processors == 1:
	file.write(exe + " " + filename + " >& " + base + "out \n")
else:
	file.write("cat $PBS_NODEFILE > " + base + "nodes\n")
	file.write("mpirun -machinefile $PBS_NODEFILE -np " + str(processors) + \
		   " " + exe + " " + filename + \
		   " >& " + base + "out\n")

if emailaddress != "your@address.edu":
	file.write("csh -c \"cat " + base + "run; cat " + base + "out;\" | /bin/mailx -s \"[qsc] " + base + "out\" " + emailaddress + "\n\n\n")
file.close()

handle = os.popen("echo `hostname`")
hostname = handle.readline()
hostname = hostname.rstrip()

print "Submitting " + jobname + " to " + str(processors) + " processors on " + hostname + " for " + str(hours) + " hours"
command = "qsub " + base +"run"
print command
os.system(command)





