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

def getParam(filename, param, default):
    file = open(filename,'r')
    value = ""
    file.seek(0)
    for line in file:
	if string.find(line,param) == 0:
	    value = file.next().strip()    	    
	    file.close()
	    return value
    file.close()
    return default

def setParam(filename, param, new, after="run_type"):
    file = open(filename,'r').readlines()
    
    #look for where the parameter is used
    for i in range(len(file)):
	if string.find(file[i],param) == 0:
	    file[i+1] = " " + str(new) + "\n"
	    open(filename,'w').writelines(file)
	    return

    #ok, it's not used, so write the parameter right after 'after'
    for i in range(len(file)):
	if string.find(file[i],after) == 0:
	    file.insert(i+2,param + "\n " + str(new) + "\n")
	    open(filename,'w').writelines(file)
	    return
	
    print "Warning: can not set parameter "+param
    return

default_nodes = 1
default_ppn   = 4
default_hours = 15
max_hours     = 15
default_queue = "largeq"

emailaddress = "your@address.edu"

if len(sys.argv) < 2: 
	print "bsubqmcbeaverrun.py <exe> <filename> " + \
		"[number of nodes="+str(default_nodes)+"] " \
		"[processors per node="+str(default_ppn)+"] " \
		"[hours=" + str(default_hours) + "] " \
		"[queue name=" + default_queue + "] "
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

if hours > max_hours:
    print "You requested ",hours," which is higher than the max of ",max_hours
    hours = max_hours

if len(sys.argv) > 6:
	queue      = sys.argv[5]
else:
	queue      = default_queue

processors = nodes*ppn

base = filename[:-4]
file = open(base+'run','w')
    
use_check = int(getParam(filename,"checkpoint",0))
use_avail = int(getParam(filename,"use_available_checkpoints",0))
use_check = (use_check == 1) or (use_avail == 1)

if use_check or use_avail:
    check_in = getParam(filename,"checkpoint_input_name",base[:-1])
    check_out = getParam(filename,"checkpoint_output_name",base[:-1])

    if use_check:
	print "Using checkpoint output directory: " + check_out    
    if os.path.isdir(check_in) and use_avail:
	print "Using checkpoint  input directory: " + check_in
    
    temp_dir = "/scratch2/"+str(os.getlogin())+"/"+check_out
    setParam(filename,"temp_dir",temp_dir,"checkpoint")
    print "Setting temp_dir: " + temp_dir
else:
    print "We are neither reading nor writing checkpoints."

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
    if hours == 0:
	#5 minutes minimum
	file.write("#BSUB -W 0:05\n")
    else:
	file.write("#BSUB -W "+ str(hours) + ":00\n")
file.write("#BSUB -q " + queue + "           \n")
file.write("#BSUB -u nitroamos@gmail.com           \n")
file.write("#BSUB -wa 'URG' -wt '10'        \n")
file.write("\n")

file.write("set CWD="+os.getcwd()+"\n")
if use_check == 1 or int(getParam(filename,"optimize_Psi",0)) == 1:
    file.write("mkdir $CWD/"+check_out+"\n")
file.write("cd $CWD\n")

if use_check == 1:
    #set the temp directory for checkpoint files
    file.write("set TD="+temp_dir+"\n")
    file.write("/bin/rm -r $TD\n")
    file.write("mkdir $TD\n")
    if os.path.isdir(check_in):
	file.write("cp "+check_in+"/* $TD/\n")
    file.write("\n")

if debug == 1:
    file.write("module load totalview_default\n")

if processors == 1:
    file.write(exe + " " + filename + " >& " + base + "out \n")
else:
    file.write("mpirun -np " + str(processors) + \
	       " " + exe + " " + filename + \
	       " >& " + base + "out\n")

if use_check == 1:
    #copy all the files back to our local storage
    file.write("\n")    
    file.write("cp $TD/* $CWD/"+check_out+"\n")
    file.write("\n")
    
if emailaddress != "your@address.edu":
	file.write("csh -c \"cat " + base + "run; cat " + base + "out;\" | /bin/mailx -s \"[cy] " + base + "out\" " + emailaddress + "\n\n\n")
file.close()

print "Submitting " + filename + " to " + str(processors) + " processors on " + queue + " for " + str(hours) + " hours"
command = "bsub < " + base +"run"
print command
os.system(command)




