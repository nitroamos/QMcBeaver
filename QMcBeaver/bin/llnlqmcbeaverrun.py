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
default_ppn   = 8
default_bank  = "caltech"

emailaddress = "your@address.edu"

if len(sys.argv) < 2: 
	print "llnlqmcbeaverrun.py <exe> <filename> " + \
		"[number of nodes="+str(default_nodes)+"] " \
		"[processors per node="+str(default_ppn)+"] " \
		"[bank name=" + default_bank + "]"
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
	queue      = sys.argv[5]
else:
	queue      = default_bank

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

file.write("#!/bin/csh         # Sets your shell\n")
#file.write("#PSUB -s /bin/csh  # Sets your shell in batch\n")
#file.write("#PSUB -c 100Mb     # Select a 100Mb memory limit\n")
file.write("#PSUB -ln " + str(nodes) + "        # Number of nodes you want to use\n")
file.write("#PSUB -g " + str(processors) + "         # Number of tasks to use\n")
file.write("#PSUB  -eo         # Send std error & std out to the same file\n")
file.write("#PSUB  -ro\n")
file.write("#PSUB -tM 12:00       # Select your time limit.\n")
file.write("#PSUB -b " + queue + "\n")
file.write("#PSUB -me          # mail at end of calc\n")
file.write("#PSUB -mb          # mail at begining of calc\n")
file.write("\n")
file.write("cd " + os.getcwd() + "\n")
file.write("rm " + base + "out\n")
file.write("poe " + exe + " " + filename + " >& " + base + "out -coredir none -labelio no \n")
#file.write("./" + exe + " " + filename + " >& " + base + "out\n")

#if emailaddress != "your@address.edu":
#	file.write("csh -c \"cat " + base + "run; cat " + base + "out;\" | /bin/mailx -s \"[up] " + base + "out\" " + emailaddress + "\n\n\n")
file.close()

print "Submitting " + filename + " to " + str(processors) + " processors."
command = "psub " + base +"run"
print command
os.system("echo; cat " + base + "run")
os.system(command)
