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
from double_factorial import *

Infile=sys.argv[1]
IN=open(Infile,'r')
gamesdata=IN.readlines()
IN.close()

Outfile=sys.argv[1][0:len(sys.argv[1])-4]+".ckmf"
OUT=open(Outfile,'w')

ANGtoBOHRconversion = 1.0/0.529177249

##################  EXTRACT GEOMETRY: START #########################

# determine if single point or geometry opt
geopt = 0
for line in gamesdata:
    if string.find(line,'RUNTYP=OPTIMIZE') != -1: geopt = 1

if geopt == 0:
    # extract the initial geometry in BOHR
    for i in range(len(gamesdata)):
        if string.find(gamesdata[i], 'RUN TITLE') != -1:
            start = i
        if string.find(gamesdata[i], 'INTERNUCLEAR DISTANCES') != -1:
            end = i
            break
    data = gamesdata[start:end]
    geometry = []
    start = 0
    for line in data:
        if start : geometry = geometry + [line]
        if string.find(line,'CHARGE') != -1: start = 1
    geometry = geometry[:len(geometry)-1]

    # split up the data
    for i in range(len(geometry)):
        geometry[i] = string.split(geometry[i])
        for j in range(2,5):
            geometry[i][j] = string.atof(geometry[i][j])

if geopt == 1:
    # extract the optimized geometry in ANG
    start = 10000000; end = 0
    for i in range(len(gamesdata)):
        if string.find(gamesdata[i], 'EQUILIBRIUM GEOMETRY LOCATED') != -1:
            start = i
        if string.find(gamesdata[i], 'INTERNUCLEAR DISTANCES') != -1 \
           and i > start:
            end = i
            break
    data = gamesdata[start:end]
    geometry = []
    start = 0
    for line in data:
        print line
        if start == 2 : geometry = geometry + [line]
        if string.find(line,'CHARGE') != -1: start = start + 1
    geometry = geometry[1:len(geometry)-1]
    print geometry
    
    # split up the data
    for i in range(len(geometry)):
        print i
        geometry[i] = string.split(geometry[i])
        for j in range(2,5):
            print j
            geometry[i][j] = string.atof(geometry[i][j])

    # convert from ANGS to BOHR if necessary
    changeunits = 0
    for line in data:
        print line
        if string.find(line,'(ANGS)') != -1: changeunits = 1
    if changeunits == 1:
        for i in range(len(geometry)):
            print i
            for j in range(2,5):
                print j
                geometry[i][j] = geometry[i][j] * ANGtoBOHRconversion

##################  EXTRACT GEOMETRY: END  #########################

##################  EXTRACT BASIS SET: BEGIN #######################

start = 0
end = 0
for i in range(len(gamesdata)):
    if string.find(gamesdata[i], 'ATOMIC BASIS SET') != -1:
        start = i
    if string.find(gamesdata[i], '$CONTRL OPTIONS') != -1:
        end = i
        break
data = gamesdata[start:end]

end = 0
for i in range(len(data)) :
    if string.find(data[i],'TOTAL NUMBER OF SHELLS') != -1 :
        end = i
        break
    if string.find(data[i],'TOTAL NUMBER OF BASIS SET SHELLS') != -1 :
        end = i
        break
basisdata = data[7:end]

basis = []
atom = []
bf = []
atomnumber = -1
oldbfnumber = 0
for line in basisdata:
    if line[0] != '\n' and line[1] != ' ':
        if bf != [] :
            atom = atom + [bf]
            bf = []
        if atom != [] : basis = basis + [atom]
        atom = []
        line = string.split(line)
        atom = atom + [line]

    elif line != '\n':
        line = string.replace(line,')',' ')
        line = string.replace(line,'(',' ')
        line = string.split(line)
        newbfnumber = string.atoi(line[0])
        temp = [line[1]] + line[3:]
        if len(line) > 6 : temp = temp + [line[6]]
        line = temp
        if newbfnumber > oldbfnumber :
            if bf != [] : atom = atom + [bf]
            bf = [line]
            oldbfnumber = newbfnumber
        else :
            bf = bf + [line]
atom = atom + [bf]
basis = basis + [atom]

# extract some flags data from this section
for line in data:
    if string.find(line,'TOTAL NUMBER OF BASIS FUNCTIONS') !=-1 :
        nbasisfunc = string.split(line)[6]
    if string.find(line,'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS') !=-1 :
        nbasisfunc = string.split(line)[7]
    if string.find(line,'CHARGE OF MOLECULE') !=-1 :
        charge = string.split(line)[4]
    if string.find(line,'TOTAL NUMBER OF ATOMS') !=-1 :
        atoms = string.split(line)[5]
    if string.find(line,'NUMBER OF OCCUPIED ORBITALS (ALPHA)') != -1 :
        nalpha = string.atoi(string.split(line)[6])
    if string.find(line,'NUMBER OF OCCUPIED ORBITALS (BETA )') != -1 :
        nbeta = string.atoi(string.split(line)[7])

##################  EXTRACT BASIS SET: END #######################

##################  EXTRACT WAVEFUNCTION: BEGIN ##################

start = -1
end = -1
for i in range(len(gamesdata)):
    if string.find(gamesdata[i],'MOLECULAR ORBITALS') !=-1 :
        start = i
    if string.find(gamesdata[i],'EIGENVECTORS') !=-1 :
        start = i
    if string.find(gamesdata[i],'ENERGY COMPONENTS') != -1 :
        end = i - 2
data = gamesdata[start+2:end]
for i in range(len(data)) :
    if string.find(data[i],'......') !=-1 :
        end = i
        data = data[:end]
        break
data = data + ['\n']

temp = []
oldbreak = 0
for i in range(len(data)) :
    if data[i] == '\n' :
        newbreak = i
        temp = temp + [data[oldbreak+4:newbreak]]
        oldbreak = newbreak

for i in range(len(temp)) :
    for j in range(len(temp[i])) :
        temp[i][j] = string.split(temp[i][j])
        if temp[i][j][2] == 'XXX' or \
           temp[i][j][2] == 'YYY' or \
           temp[i][j][2] == 'ZZZ' or \
           temp[i][j][2] == 'XXY' or \
           temp[i][j][2] == 'XXZ' or \
           temp[i][j][2] == 'YYX' or \
           temp[i][j][2] == 'YYZ' or \
           temp[i][j][2] == 'ZZX' or \
           temp[i][j][2] == 'ZZY' or \
           temp[i][j][2] == 'XYZ' :
            end = 3
        else : end = 4
        temp[i][j] = temp[i][j][end:]
temp = temp[1:]

wavefunction = []
line = []
for i in range(string.atoi(nbasisfunc)) :
    for j in range(len(temp)) :
        line = line + temp[j][i]
    wavefunction = wavefunction + [line]
    line = []

norbitals = len(wavefunction[0])
        
##################  EXTRACT WAVEFUNCTION: END   ##################

##################  EXTRACT ENERGY: BEGIN  #######################

start = -1
for i in range(len(gamesdata)):
    if string.find(gamesdata[i],'ENERGY COMPONENTS') !=-1 :
        start = i
for line in gamesdata[start:start+20]:
    if string.find(line,'TOTAL ENERGY') !=-1 :
        energy = string.split(line)[3]
        
##################  EXTRACT ENERGY: END   ########################

##################  SET RUN TYPE: BEGIN   ########################

run_type = 'variational'

##################  SET RUN TYPE: END     ########################

##################  PRINT FLAGS: BEGIN    ########################

OUT.write('&flags\n')
OUT.write('atoms\n %s\n'%atoms)
OUT.write('charge\n %s\n'%charge)
OUT.write('energy\n %s\n'%energy)
OUT.write('norbitals\n %s\n'%norbitals)
OUT.write('nbasisfunc\n %s\n'%nbasisfunc)
OUT.write('run_type\n %s\n'%run_type)
OUT.write('chip_and_mike_are_cool\n Yea_Baby!\n')
OUT.write('&\n')

##################  PRINT FLAGS: END      #######################        

##################  PRINT GEOMETRY: BEGIN #######################

OUT.write('&geometry\n')
for line in geometry:
    OUT.write('%s\t%i\t%f\t%f\t%f\n'%(line[0],string.atof(line[1]),line[2],line[3],line[4]))
OUT.write('&\n')

##################  PRINT GEOMETRY: END  ########################

##################  PRINT BASIS: BEGIN   ########################

# calculate the number of basis functions for atom and maximum gaussians
# in any basis function
for i in range(len(basis)):
    label = basis[i][0][0]
    basis[i] = basis[i][1:]
    nbf = 0
    maxgaussian = 0
    for bf in basis[i]:
        if bf[0][0] == 'S' : nbf = nbf + 1
        if bf[0][0] == 'P' : nbf = nbf + 3
        if bf[0][0] == 'D' : nbf = nbf + 6
        if bf[0][0] == 'F' : nbf = nbf + 10
        if bf[0][0] == 'L' : nbf = nbf + 4
        if len(bf) > maxgaussian : maxgaussian = len(bf)
    basis[i] = [[label,nbf,maxgaussian]] + basis[i]

OUT.write('&basis\n')
for atom in geometry :
    for ATOM in basis:
        if atom[0] == ATOM[0][0] :
            atomicbasis = ATOM
    head = atomicbasis[0]
    atomicbasis = atomicbasis[1:]
    OUT.write('%s\t%i\t%i\n'%(head[0],head[1],head[2]))
    for pbf in atomicbasis :
        if pbf[0][0] == 'S' :
            OUT.write('\t%i\t%s\n'%(len(pbf),'s'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,0,0,gaussian[2],gaussian[1])))
        if pbf[0][0] == 'P' :
            OUT.write('\t%i\t%s\n'%(len(pbf),'px'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(1,0,0,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'py'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,1,0,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'pz'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,0,1,gaussian[2],gaussian[1])))
        if pbf[0][0] == 'D' :
            OUT.write('\t%i\t%s\n'%(len(pbf),'dxx'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(2,0,0,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'dyy'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,2,0,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'dzz'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,0,2,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'dxy'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(1,1,0,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'dxz'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(1,0,1,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'dyz'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,1,1,gaussian[2],gaussian[1])))
        if pbf[0][0] == 'F' :
            OUT.write('\t%i\t%s\n'%(len(pbf),'fxxx'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(3,0,0,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'fyyy'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,3,0,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'fzzz'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,0,3,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'fxxy'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(2,1,0,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'fxxz'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(2,0,1,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'fyyx'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(1,2,0,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'fyyz'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,2,1,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'fzzx'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(1,0,2,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'fzzy'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,1,2,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'fxyz'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(1,1,1,gaussian[2],gaussian[1])))
        if pbf[0][0] == 'L' :
            OUT.write('\t%i\t%s\n'%(len(pbf),'s'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,0,0,gaussian[2],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'px'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(1,0,0,gaussian[3],gaussian[1])))
            OUT.write('\t%i\t%s\n'%(len(pbf),'py'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,1,0,gaussian[3],gaussian[1])))   
            OUT.write('\t%i\t%s\n'%(len(pbf),'pz'))
            for gaussian in pbf :
                OUT.write('\t\t%s\t%s\n'%(gaussian[1], \
                             normalize(0,0,1,gaussian[3],gaussian[1])))
                
OUT.write('&\n')

##################  PRINT BASIS: END     ########################

##################  PRINT WAVEFUNCTION: BEGIN ###################

OUT.write('&wavefunction\n')

for i in range(norbitals) :
    if i < nalpha : OUT.write('1 ')
    else : OUT.write('0 ')
    if i < nbeta : OUT.write('1\n')
    else : OUT.write('0\n')
    for j in range(string.atoi(nbasisfunc)) :
        OUT.write('\t%s'%wavefunction[j][i])
        if (j+1)%3 == 0  : OUT.write('\n')
        else : OUT.write('\t')
    OUT.write('\n\n')

OUT.write('&\n')

##################  PRINT WAVEFUNCTION: END   ###################


def periodic_table(i):
    mass=0.0
    i=int(i)
    if(i==1):
        mass=1.0
    elif(i==2):
        mass=4.0
    elif(i==3):
        mass=7.0
    elif(i==4):
        mass=9.0
    elif(i==5):
        mass=11.0
    elif(i==6):
        mass=12.0
    elif(i==7):
        mass=14.0
    elif(i==8):
        mass=16.0
    elif(i==9):
        mass=19.0
    elif(i==10):
        mass=20.0
    elif(i==11):
        mass=23.0
    elif(i==12):
        mass=24.0
    elif(i==13):
        mass=27.0
    elif(i==14):
        mass=28.0
    elif(i==15):
        mass=31.0
    elif(i==16):
        mass=32.0
    elif(i==17):
        mass=35.0
    elif(i==18):
        mass=40.0
    elif(i==19):
        mass=39.0
    elif(i==20):
        mass=40.0
    elif(i==21):
        mass=45.0
    elif(i==22):
        mass=48.0
    elif(i==23):
        mass=51.0
    elif(i==24):
        mass=52.0
    elif(i==25):
        mass=55.0
    elif(i==26):
        mass=56.0
    elif(i==27):
        mass=59.0
    elif(i==28):
        mass=59.0
    elif(i==29):
        mass=64.0
    elif(i==30):
        mass=65.0
    elif(i==31):
        mass=70.0
    elif(i==32):
        mass=73.0
    elif(i==33):
        mass=75.0
    elif(i==34):
        mass=79.0
    elif(i==35):
        mass=80.0
    elif(i==36):
        mass=84.0
    elif(i==37):
        mass=85.0
    elif(i==38):
        mass=88.0
    elif(i==39):
        mass=89.0
    elif(i==40):
        mass=91.0
    elif(i==41):
        mass=93.0
    elif(i==42):
        mass=96.0
    elif(i==43):
        mass=98.0
    elif(i==44):
        mass=101.0
    elif(i==45):
        mass=103.0
    elif(i==46):
        mass=106.0
    elif(i==47):
        mass=108.0
    elif(i==48):
        mass=112.0
    elif(i==49):
        mass=115.0
    elif(i==50):
        mass=119.0
    elif(i==51):
        mass=122.0
    elif(i==52):
        mass=128.0
    elif(i==53):
        mass=127.0
    elif(i==54):
        mass=131.0
    elif(i==55):
        mass=133.0
    elif(i==56):#this has become tiresome, finish later if needed
        mass=138.0
    else:
        print "Need to pick a real atom type in periodic_table in GamesToCkMf.py"
        exit(1)
    conversion=1833.15038419
    return mass*conversion



################## PRINT JASTROW: BEGIN ##############################

# Make a list of all the different atom types
atom_types = []
atom_type_charges = []
for atom in geometry:
    is_in_list = 0;
    for atom_type in atom_types:
        if atom[0] == atom_type:
            is_in_list = 1;
    if not is_in_list:
        atom_types = atom_types + [atom[0]]
        atom_type_charges = atom_type_charges + [atom[1]]

# write out the jastrow
OUT.write('\n&Jastrow\n\n')

# up down jastrow
if nalpha > 0 and nbeta > 0:
    OUT.write('ParticleTypes: Electron_Up Electron_Down\n')
    OUT.write('CorrelationFunctionType: FixedCuspPade\n')
    OUT.write('NumberOfParameterTypes: 2\n')
    OUT.write('NumberOfParametersOfEachType: 0 1\n')
    OUT.write('Parameters: 3.0\n')
    OUT.write('NumberOfConstantTypes: 1\n')
    OUT.write('NumberOfConstantsOfEachType: 1\n')
    OUT.write('Constants: 0.5\n')      
    OUT.write('\n')

# up up jastrow
if nalpha > 1:
    OUT.write('ParticleTypes: Electron_Up Electron_Up\n')
    OUT.write('CorrelationFunctionType: FixedCuspPade\n')
    OUT.write('NumberOfParameterTypes: 2\n')
    OUT.write('NumberOfParametersOfEachType: 0 1\n')
    OUT.write('Parameters: 100.0\n')
    OUT.write('NumberOfConstantTypes: 1\n')
    OUT.write('NumberOfConstantsOfEachType: 1\n')
    OUT.write('Constants: 0.25\n')      
    OUT.write('\n')
    
# down down jastrow
if nbeta > 1:
    OUT.write('ParticleTypes: Electron_Down Electron_Down\n')
    OUT.write('CorrelationFunctionType: FixedCuspPade\n')
    OUT.write('NumberOfParameterTypes: 2\n')
    OUT.write('NumberOfParametersOfEachType: 0 1\n')
    OUT.write('Parameters: 100.0\n')
    OUT.write('NumberOfConstantTypes: 1\n')
    OUT.write('NumberOfConstantsOfEachType: 1\n')
    OUT.write('Constants: 0.25\n')      
    OUT.write('\n')
    
# up nuclear jastrow
if nalpha > 0:
    for i in range(len(atom_types)):
        OUT.write('ParticleTypes: Electron_Up ' + atom_types[i] + '\n')
        OUT.write('CorrelationFunctionType: FixedCuspPade\n')
        OUT.write('NumberOfParameterTypes: 2\n')
        OUT.write('NumberOfParametersOfEachType: 0 1\n')
        OUT.write('Parameters: 100.0\n')
        OUT.write('NumberOfConstantTypes: 1\n')
        OUT.write('NumberOfConstantsOfEachType: 1\n')
        OUT.write('Constants: -' + atom_type_charges[i] + '\n')      
        OUT.write('\n')

# down nuclear jastrow
if nbeta > 0:
    for i in range(len(atom_types)):
        OUT.write('ParticleTypes: Electron_Down ' + atom_types[i] + '\n')
        OUT.write('CorrelationFunctionType: FixedCuspPade\n')
        OUT.write('NumberOfParameterTypes: 2\n')
        OUT.write('NumberOfParametersOfEachType: 0 1\n')
        OUT.write('Parameters: 100.0\n')
        OUT.write('NumberOfConstantTypes: 1\n')
        OUT.write('NumberOfConstantsOfEachType: 1\n')
        OUT.write('Constants: -' + atom_type_charges[i] + '\n')      
        OUT.write('\n')

OUT.write('&Jastrow\n')

################## PRINT JASTROW: END ################################

