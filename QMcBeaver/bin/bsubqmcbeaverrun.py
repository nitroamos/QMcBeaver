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

#this script will email you <input>out when it's done if you
#change this address
emailaddress = "your@address.edu"

if len(sys.argv) < 2: 
	print "bsubqmcbeaverrun.py <exe> <filename> <number of nodes> <processors per node> <queue name>"
	print "(you can edit this script for it to email you the output file!)"
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
	ppn        = "4"

if len(sys.argv) > 5:
	queue      = sys.argv[5]
else:
#	queue      = "largeq"
	queue      = "smallq"

processors = string.atoi(nodes)*string.atoi(ppn)

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
	file.write("#BSUB -W 8:00                    \n")
file.write("#BSUB -q " + queue + "           \n")
file.write("#BSUB -wa 'URG' -wt '10'        \n")
file.write("\n")
file.write("cd " + os.getcwd() + "\n")
if debug == 1:
	file.write("module load totalview_default\n")

if processors == 1:
	file.write(exe + " " + filename + " >& " + base + "out \n")
else:
	if debug == 1:
		file.write("totalview prun -a -t " + " " + exe + " " + filename +  " >& " + base + "out\n")
	else:
		file.write("prun " + " " + exe + " " + filename +  " >& " + base + "out\n")
if emailaddress != "your@address.edu":
	file.write("csh -c \"cat " + base + "run; cat " + base + "out;\" | /bin/mailx -s \"[qsc] " + base + "out\" " + emailaddress + "\n\n\n")
file.close()

print "Submitting " + filename + " to " + str(processors) + " processors on " + queue
command = "bsub < " + base +"run"
print command
os.system(command)





