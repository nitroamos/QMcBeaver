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

# This program will plot .qmc files using gnuplot
# syntax:  plotQMC.py file.qmc

import sys
import string
import os

def GetEnergy(QMCfile):
    try:
        ckmfFile = open(QMCfile[:len(QMCfile)-3]+'ckmf')
        data = ckmfFile.readlines()
        ckmfFile.close()
        EnergyPos = -1
        for i in range(len(data)):
            if data[i][:-1] == 'energy':
                print data[i]
                EnergyPos = i+1
        if(len(data) < 1 or EnergyPos == -1):  return 'noEnergy'
        return data[EnergyPos]
    except IOError, e:
        return 'noEnergy'


QMCfile = sys.argv[1]
GNUfile = QMCfile[:len(QMCfile)-3]+'gnp'

Loop = 1

while Loop:
    Printfile = open(GNUfile,'w')
    Printfile.write("#plot %s\n" %(QMCfile))
    Printfile.write('plot \"'+ QMCfile + '\" u 1:2 w l')
    Printfile.write(' ,"'+ QMCfile + '\" u 1:4 w l')
    Printfile.write(' , \"'+ QMCfile + '\" u 1:5 w l')
    Energy = GetEnergy(QMCfile)
    print Energy
    if(Energy == 'noEnergy'):
        Printfile.write('\n')
    else:
        Printfile.write(' , ' + GetEnergy(QMCfile) + '\n' )
    Printfile.write('pause -1 "Type Return to Continue"')
    Printfile.close()

    Command = "gnuplot " + GNUfile
    print Command
    os.system(Command)

    Command = "rm " + GNUfile
    os.system(Command)
    
    Quit = raw_input("Type \'x\' to exit or return to continue:")
    if (Quit == 'x'): Loop = 0

