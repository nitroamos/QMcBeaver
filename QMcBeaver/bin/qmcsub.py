#!/usr/bin/env python
# QMcBeaver queue submission script

import sys
import string
import os

######## Set up your default settings here
default_procs = 1
default_hours = 100
default_queue = ""
ppn           = 2
mpi           = "LAMPI"
qserver       = "qsub"
max_hours     = 1000
scratch       = "/temp1"

# Set to 0 for no email, 1 for queue emails, 2 for all emails
emailme       = 1

# Search this script for how this is used
debug = 0
#debug = 1    

######## Set up defaults for machines you use
handle   = os.popen("echo `hostname`")
host     = handle.readline()
host     = host.strip()
hostname = (host.split('.'))[0]

if string.find(hostname, "hive") != -1:
    mpi           = "OpenMPI"
elif string.find(hostname, "cy") != -1:
    ppn           = 4
    mpi           = "LAMPI"
    qserver       = "bsub"
    default_hours = 15
    max_hours     = 15
    default_queue = "largeq"
    scratch       = "/scratch2"
    emailme       = 1
elif string.find(hostname, "uP") != -1:
    ppn           = 8
    mpi           = "poe"
    qserver       = "psub"
    default_hours = 12
    max_hours     = 12
    default_queue = "caltech"
    emailme       = 1
    
######## Start reading input parameters
if len(sys.argv) <= 2:
    print "Script usage:"
    print sys.argv[0] +" <exe> <ckmf file> " + \
	  "[processors="+str(default_procs)+"] " \
	  "[hours=" + str(default_hours) + "] " \
	  "[queue name=" + default_queue + "] "
    sys.exit(0)
else:
    exe        = sys.argv[1]
    filename   = sys.argv[2]   
    processors = default_procs
    hours      = default_hours
    queue      = default_queue

if len(sys.argv) > 3:
    processors = string.atoi(sys.argv[3])
else:
    if processors == 1:
	mpi = "none"

if len(sys.argv) > 4:
    hours = string.atoi(sys.argv[4])
    
if len(sys.argv) > 5:
    queue = sys.argv[5]

######## Validate job
if not os.path.isfile(filename):
    print "This is not a valid input file: " + filename
    sys.exit(0)
if filename[-5:] != ".ckmf":
    print "Input file is not recognized as a ckmf file: " + filename
    sys.exit(0)
if not os.path.isfile(exe):
    print "This is not a valid executable: " + exe
    sys.exit(0)

base = filename[:-4]
jobname = base[:-1]

username = os.getlogin()
emailaddress = username + "@" + host

#change this line to override username@host
#emailaddress = "your@address.edu"

if emailme > 0:
    print "Sending emails to " + emailaddress
    
######## Set up machine parameters
if processors <= ppn:
    nodes = 1
    ppn   = processors
else:
    if processors%ppn != 0:
	print "Warning: rounding "+str(processors)+" processors to fill nodes."
    nodes = processors / ppn
    processors = nodes * ppn

if hours > max_hours:
    print "You requested ",hours," which is higher than the max of ",max_hours
    hours = max_hours

if hours == 0 or debug:
    time = "00:10"
else:
    time = str(hours) + ":00"

######## Start writing queue submission script
file = open(base+'run','w')
file.write("#!/bin/csh    \n")
if qserver == "qsub":
    file.write("#PBS -N " + jobname + " \n")
    file.write("#PBS -j oe -k n         \n") # -j oe merge std out and std err, -k n don't keep the file
    file.write("#PBS -l nodes=" + str(nodes) + ":ppn=" + str(ppn) + ",walltime=" + time+":00\n")

    if emailme > 0:
	file.write("#PBS -M "+emailaddress+" \n")
	file.write("#PBS -m ea               \n") #b for begin, e for end, a for abort

    if default_queue != "":
	file.write("#PBS -q " + queue + " \n")

elif qserver == "bsub":
    file.write("#BSUB -J " + filename +        " \n")
    file.write("#BSUB -n " + str(processors) + " \n")
    file.write("#BSUB -W "+ time + "\n")
    file.write("#BSUB -wa 'URG' -wt '10'         \n")

    if default_queue != "":
	file.write("#BSUB -q " + queue + "       \n")
	
    if debug == 1:
	file.write("#BSUB -I                    \n")
    else:
	file.write("#BSUB -B                    \n")

    if emailme > 0:
	file.write("#BSUB -u "+emailaddress+"   \n")

elif qserver == "psub":
    file.write("#PSUB -ln " + str(nodes) + "    # Number of nodes you want to use\n")
    file.write("#PSUB -g " + str(processors) + "# Number of tasks to use\n")
    file.write("#PSUB  -eo                      # Send std error & std out to the same file\n")
    file.write("#PSUB  -ro\n")
    file.write("#PSUB -tM "+time+"              # Select your time limit.\n")

    if emailme > 0:
	file.write("#PSUB -me                   # mail at end of calc\n")
	file.write("#PSUB -mb                   # mail at begining of calc\n")

    if default_queue != "":
	#this actually selects the name of a bank, not a queue
	file.write("#PSUB -b " + queue +"\n")
else:
    print "I don't know anything about qserver == ",qserver
    sys.exit(0)
file.write("\n")

######## Start looking in the ckmf file for some details regarding
######## directory usage
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

is_opt    = int(getParam(filename,"optimize_Psi",0))
use_check = int(getParam(filename,"checkpoint",0))
use_avail = int(getParam(filename,"use_available_checkpoints",0))
use_check = (use_check == 1) or (use_avail == 1)

if use_check or use_avail or is_opt:
    check_in = getParam(filename,"checkpoint_input_name",base[:-1])
    check_out = getParam(filename,"checkpoint_output_name",base[:-1])

    if use_check or is_opt:
	print "Using output directory: " + check_out    
    if os.path.isdir(check_in) and use_avail:
	print "Using  input directory: " + check_in
    
    temp_dir = scratch + "/" + username +"/"+ check_out
    setParam(filename,"temp_dir",temp_dir,"checkpoint")
    print "Setting temp_dir: " + temp_dir
else:
    print "We are not using any auxilliary directories"

######## Set up the local directory
file.write("set CWD="+os.getcwd()+"\n")
if use_check or is_opt:
    file.write("mkdir $CWD/"+check_out+"\n")
file.write("cd $CWD\n")

######## Set up the scratch directory, if we're using it
if use_check == 1:
    file.write("set TD="+temp_dir+"\n")
    file.write("/bin/rm -r $TD\n")
    file.write("mkdir $TD\n")
    if os.path.isdir(check_in):
	file.write("cp "+check_in+"/* $TD/\n")
    file.write("\n")

######## This line will remove old files. This is helpful if we use the same input multiple times.
os.system("/bin/rm -f " + base + "o*\n")
os.system("/bin/rm -f " + base + "nodes\n")

######## The executable line
if mpi == "none":
    file.write(exe + " " + filename + " >& " + base + "out \n")
elif mpi == "OpenMPI":
    file.write("cat $PBS_NODEFILE > " + base + "nodes\n")
    file.write("mpirun -machinefile $PBS_NODEFILE -np " + str(processors) + \
	       " " + exe + " " + filename + \
	       " >& " + base + "out\n")
elif mpi == "LAMPI":
    file.write("lamboot\n")
    file.write("mpirun -np " + str(processors) + \
	       " " + exe + " " + filename + \
	       " >& " + base + "out\n")
    file.write("lamhalt\n")
elif mpi == "poe":
    file.write("poe " + exe + " " + filename + " >& " + base + "out -coredir none -labelio no \n")
else:
    print "Error: I don't know anything about mpi == ",mpi
    sys.exit(0)

######## If we wrote checkpoint files, then we want to copy them back from the scratch directory
if use_check == 1:
    file.write("\n")    
    file.write("cp $TD/* $CWD/"+check_out+"\n")
    file.write("\n")

######## Do you want to be emailed the output file?
if emailme == 2:
    file.write("csh -c \"cat " + base + "run; cat " + base + "out;\" | /bin/mailx -s \"["+hostname+"] " + base + "out\" " + emailaddress + "\n\n\n")

######## Wrap up and submission
file.close()
print "Submitting " + jobname + " to " + str(processors) + " processors on " + hostname + " for " + time + " hours"
command = qserver + " <  " + base +"run"
print command
os.system("cat "+base+"run")
os.system(command)





