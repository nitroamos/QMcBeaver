#!/usr/bin/env python

import string
import sys
import glob

if len(sys.argv) < 3:
    print "Useage: % configs2xyz.py <QMcBeaver Input> <QMcBeaver Checkpoint>"
    sys.exit(0)


input_name=argv[1]
config_name=argv[2]

#read in the *.ckmf geometry
CKMF = open(input_name)
ckmflines = CKMF.readlines()

#put atoms in geom list
geom=[]
for i in range(len(ckmflines)):
    if(1<=string.count(ckmflines[i],"&geometry")):
        j=1
        while(0==string.count(ckmflines[i+j],"&")):
            sline=string.split(ckmflines[i+j])
            atom=[]
            #atom type
            atom.append(sline[0])
            #atom charge
            atom.append(sline[1])
            #xyz
            atom.append(sline[2])
            atom.append(sline[3])
            atom.append(sline[4])
            j=j+1
            geom.append(atom)

#read in the config lines
CONFIG = open(config_name)
configlines = CONFIG.readlines()

#make a list of configs xyz
configs=[]

#fill each element in the config list
for i in range(len(configlines)):
    sline=string.split(configlines[i])
    if("R"==sline[0]):
        n_electrons=string.atoi(sline[1])
        config=[]
        for j in range(1,n_electrons+1):
            electron=[]
            sline=string.split(configlines[i+j])
            #xyz
            electron.append(sline[1])
            electron.append(sline[2])
            electron.append(sline[3])
            config.append(electron)
        configs.append(config)

#print out in a format Mark would like

#test which just prints out the geom and first four configs
for atom in geom:
    print atom

for i in range(4):
    print ""
    for electron in configs[i]:
        print electron






