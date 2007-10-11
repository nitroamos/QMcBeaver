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
default_ppn   = 4
default_hours = 12
default_queue = "largeq"

emailaddress = "your@address.edu"

if len(sys.argv) < 2: 
	print "bsubqmcbeaverrun.py <exe> <filename> " + \
		"[number of nodes="+str(default_nodes)+"] " \
		"[processors per node="+str(default_ppn)+"] " \
		"[queue name=" + default_queue + "]"
	sys.exit(0)

exe        = sys.argv[1]
filename   = sys.argv[2]   

if not os.path.isfile(filename):
    print "This is not a valid input file: " + filename
    sys.exit(0)
if not os.path.isfile(exe):
    print "This is not a valid executable: " + exe
    sys.exit(0)
    
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

if len(sys.argv) > 6:
	queue      = sys.argv[5]
else:
	queue      = default_queue

processors = nodes*ppn

base = filename[:-4]
file = open(base+'run','w')
ckmf_file = open(filename,'r')
ckmf_data = ckmf_file.readlines()
ckmf_file.close()

for i in range(len(ckmf_data)):
	if string.find(ckmf_data[i],"temp_dir") == 0:
		temp_dir = ckmf_data[i+1][:-1]
		break

print "Using temp space: " + temp_dir

debug = 0
#debug = 1

file.write("#!/bin/csh    \n")
file.write("#BSUB -J " + filename +        " \n")
file.write("#BSUB -n " + str(processors) + " \n")
if debug == 1:
	file.write("#BSUB -I                         \n")
	file.write("#BSUB -W 0:30                    \n")
else:
	file.write("#BSUB -B                         \n")
	#file.write("#BSUB -W 8:00                    \n")
	file.write("#BSUB -W "+ str(hours) + ":00\n")
file.write("#BSUB -q " + queue + "           \n")
file.write("#BSUB -wa 'URG' -wt '10'        \n")
file.write("\n")
file.write("cd " + os.getcwd() + "\n")
if debug == 1:
	file.write("module load totalview_default\n")

if processors == 1:
	file.write(exe + " " + filename + " >& " + base + "out \n")
else:
	file.write("mpirun  -np " + str(processors) + \
		   " " + exe + " " + filename + \
		   " >& " + base + "out\n")

if emailaddress != "your@address.edu":
	file.write("csh -c \"cat " + base + "run; cat " + base + "out;\" | /bin/mailx -s \"[qsc] " + base + "out\" " + emailaddress + "\n\n\n")
file.close()

print "Submitting " + filename + " to " + str(processors) + " processors on " + queue + " for " + str(hours) + " hours"
command = "bsub < " + base +"run"
print command
os.system(command)





