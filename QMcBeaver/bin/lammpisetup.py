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

import os
import sys
import string

Machines = sys.argv[1:]
if len(Machines) < 1:
    print "Useage: lammpisetup.py machine1 <machine2> <machine3> ..."
    sys.exit(0);

print "MAKING LAMHOSTS*********************************************"
hosts=open("lamhosts","w")
hosts.write("# this is the lamhost file which needs to list all machines which will be accessed\n")

for machine in Machines:
    hosts.write(machine+".wag.caltech.edu\n")
hosts.close()


print "RECON*********************************************"
os.system("recon -v lamhosts")
os.system("recon -v lamhosts > recon.mpiout")
print "LAMBOOT"
os.system("lamboot -v lamhosts")
os.system("lamboot -v lamhosts > lamboot.mpiout")
lamboot=open("lamboot.mpiout","r")
data_lamboot=lamboot.readlines()
if((-1)==string.find(data_lamboot[len(data_lamboot)-1],"topology") or
   (-1)==string.find(data_lamboot[len(data_lamboot)-1],"done")):
    print "ERROR in **lamboot -v lamhosts** command"
    print data_lamboot[len(data_lamboot)-1]
    exit(1)


print "TPING*********************************************"
os.system("tping -c1 N")
os.system("tping -c1 N >tping.mpiout")
tping=open("tping.mpiout","r")
data_tping=tping.readlines()
if((-1)==string.find(data_tping[len(data_tping)-1],"roundtrip")):
    print "ERROR in **tping -c1 N** command"
    print data_tping[len(data_tping)-1]
    exit(1)



#print "RUNNING MPIRUN*********************************************"
#np=sys.argv[1]
#jobname=sys.argv[2]
#command="mpirun -np "+np+" "+jobname
#print "EXECUTING COMMAND:"
#print command
#os.system(command)
#print "COMMAND SUBMITTED"


print "DONE"

