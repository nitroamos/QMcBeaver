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
from atomic_symbol_to_Z import *
#from double_factorial import *

PI = 3.14159265359
a0 = 0.529177257507

def normalize( xexp, yexp, zexp, precoeff, expcoeff):
    precoeff = string.atof(precoeff)
    expcoeff = string.atof(expcoeff)

    value = 0
    if xexp == 0 and yexp == 0 and zexp == 0:
        value = 2.0 * expcoeff**(3.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)
    elif xexp == 1 and yexp == 0 and zexp == 0:
        value = 4.0 * expcoeff**(5.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)
    elif xexp == 0 and yexp == 1 and zexp == 0:
        value = 4.0 * expcoeff**(5.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)
    elif xexp == 0 and yexp == 0 and zexp == 1:
        value = 4.0 * expcoeff**(5.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)
    elif xexp == 2 and yexp == 0 and zexp == 0:
        value = 8.0*expcoeff**(7.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                3.0**(-0.5)
    elif xexp == 0 and yexp == 2 and zexp == 0:
        value = 8.0*expcoeff**(7.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                3.0**(-0.5)
    elif xexp == 0 and yexp == 0 and zexp == 2:
        value = 8.0*expcoeff**(7.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                3.0**(-0.5)
    elif xexp == 1 and yexp == 1 and zexp == 0:
        value = 8.0*expcoeff**(7.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)
    elif xexp == 1 and yexp == 0 and zexp == 1:
        value = 8.0*expcoeff**(7.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)
    elif xexp == 0 and yexp == 1 and zexp == 1:
        value = 8.0*expcoeff**(7.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)
    elif xexp == 3 and yexp == 0 and zexp == 0:
        value = 16.0*expcoeff**(9.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                15.0**(-0.5)
    elif xexp == 0 and yexp == 3 and zexp == 0:
        value = 16.0*expcoeff**(9.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                15.0**(-0.5)
    elif xexp == 0 and yexp == 0 and zexp == 3:
        value = 16.0*expcoeff**(9.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                15.0**(-0.5)
    elif xexp == 2 and yexp == 1 and zexp == 0:
        value = 16.0*expcoeff**(9.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                15.0**(-0.5)
    elif xexp == 2 and yexp == 0 and zexp == 1:
        value = 16.0*expcoeff**(9.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                15.0**(-0.5)
    elif xexp == 1 and yexp == 2 and zexp == 0:
        value = 16.0*expcoeff**(9.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                15.0**(-0.5)
    elif xexp == 0 and yexp == 2 and zexp == 1:
        value = 16.0*expcoeff**(9.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                15.0**(-0.5)
    elif xexp == 1 and yexp == 0 and zexp == 2:
        value = 16.0*expcoeff**(9.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                15.0**(-0.5)
    elif xexp == 0 and yexp == 1 and zexp == 2:
        value = 16.0*expcoeff**(9.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)* \
                15.0**(-0.5)
    elif xexp == 1 and yexp == 1 and zexp == 1:
        value = 16.0*expcoeff**(9.0/4.0)/PI**(3.0/4.0)/2.0**(1.0/4.0)
    else:
        print "ERROR in normalize!"
        value = 0
    return value * precoeff

Infile=sys.argv[1]
Outfile=sys.argv[1][0:len(sys.argv[1])-6]+".ckmf"
JagOutfile=sys.argv[1][0:len(sys.argv[1])-6]+".out"
JagBasisfile=sys.argv[1][0:len(sys.argv[1])-6]+".bas"

IN=open(Infile,'r')
jaguardata=IN.readlines()
IN.close()
Jagoutput=open(JagOutfile,'r')
jagoutdata=Jagoutput.readlines()
Jagoutput.close()
OUT=open(Outfile,'w')

Charge=0
Energy=0
for i in range(len(jaguardata)):
    if string.find(jaguardata[i],"molchg=") != -1:
        Charge=string.atoi(jaguardata[i][7:len(jaguardata[i])])
        
for i in range(len(jaguardata)):
    if string.find(jaguardata[i],"&zmat") == 0:
        Start = i+1
        break
for i in range(Start,len(jaguardata)):
    if jaguardata[i][0] == '&':
        End = i
        break

#
tempspot = -1
crashmark=""
Error = 0
for j in range(len(jagoutdata)):
    if string.find(jagoutdata[j],"SCFE") != -1:
        tempspot = j
    if string.find(jagoutdata[j],"ERROR") != -1:
        Error = Error + 1
if tempspot < 0:
    print "Could Not Find Energy!"
energyline=string.split(jagoutdata[tempspot])
Energy=energyline[4]    
#
    
# Get Geometry
Natoms = End-Start
Qmcdata=[]
for i in range(Start,End):
    Qmcdata=Qmcdata+[string.split(jaguardata[i])]

for i in range(len(Qmcdata)):
    Qmcdata[i][0]=string.replace(Qmcdata[i][0],"_","")
    for j in range(0,10):
        Qmcdata[i][0]=string.replace(Qmcdata[i][0],repr(j),"")

for i in range(len(Qmcdata)):
    Qmcdata[i] = [ Qmcdata[i][0], str(symbol_to_Z(Qmcdata[i][0])) ] \
                 + [string.atof(Qmcdata[i][1])/a0,string.atof(Qmcdata[i][2]) \
                    /a0,string.atof(Qmcdata[i][3])/a0] 

geometry = Qmcdata
# End Get Geometry

# Get Wavefunction
for i in range(len(jaguardata)):
    if string.find(jaguardata[i],"&guess") == 0:
        Start = i+1
        break
for i in range(Start,len(jaguardata)):
    if jaguardata[i][0] == '&':
        End = i
        break
Wavefunction=[]
Orbital=[]
LongLine=[]
Occupation=-1
for i in range(Start,End):
    Line=string.split(jaguardata[i])
    Index=string.find(jaguardata[i],"Occupation")
    if Index != -1:
        if i != 0:
             Wavefunction=Wavefunction + [[[Occupation],LongLine]]
        Occupation=Line[5]
        LongLine=[]
        continue
    LongLine=LongLine+Line
Wavefunction=Wavefunction + [[[Occupation],LongLine]]
Wavefunction=Wavefunction[1:End]
# End Get Wavefunction

# Get Basis Set

IN = open(JagBasisfile,'r')
data = IN.readlines()
IN.close()

# break into atoms
atoms = []
atom  = []
for line in data:
    if string.find(line,"!") != -1 :
        atoms = atoms + [atom]
        # print atom
        atom  = []
    else:
        atom = atom + [line[:-1]]
atoms = atoms + [atom]
atoms = atoms[1:]

for i in range(len(atoms)):
    atoms[i] = atoms[i][:1]+atoms[i][2:]

for i in range(len(atoms)):
    numbers = atoms[i][1]
    numbers = string.split(numbers)
    for j in range(len(numbers)):
        numbers[j] = string.atoi(numbers[j])

    COEFFS = atoms[i][2:]
    SPLITCOEFFS = []
    for k in range(len(numbers)):
        SPLITCOEFFS = SPLITCOEFFS + [COEFFS[:numbers[k]]]
        COEFFS = COEFFS[numbers[k]:]

    for k in range(len(SPLITCOEFFS)):
        field = SPLITCOEFFS[k]

        for j in range(len(field)):
            stuff = [field[j][:11][0]] + [string.replace(field[j][11:25], \
                         'D','E')] + [string.replace(field[j][25:],'D','E')]
            field[j] = stuff

        type = field[0][0]

        for j in range(len(field)):
            field[j] = [field[j][1]] + [field[j][2]]

        MAXPBF = 0
        for number in numbers:
            if MAXPBF < number: MAXPBF = number

        SPLITCOEFFS[k] = [type,len(field),field]

    atoms[i] = [atoms[i][0],MAXPBF,SPLITCOEFFS]

# End Get Basis Set

# Write out data
OUT.write("&flags\n")
OUT.write("atoms\n %i\n"%(Natoms))
OUT.write("charge\n %i\n"%(Charge))
OUT.write("energy\n %s\n"%(Energy))
OUT.write("norbitals\n %i\n"%(len(Wavefunction)))
OUT.write("nbasisfunc\n %i\n"%(len(Wavefunction[0][1])))
OUT.write("ndeterminants\n 1\n")
OUT.write("run_type\n variational\n")
OUT.write("chip_and_mike_are_cool\n Yea_Baby!\n")
OUT.write("&\n")

# Write Geometry
OUT.write("&geometry\n")
for i in range(len(Qmcdata)):
    for j in range(len(Qmcdata[i])):
        OUT.write(str(Qmcdata[i][j]))
        OUT.write("\t")
    OUT.write("\n")
OUT.write("&\n")
# End Write Geometry

# Write Basis Set

OUT.write("&basis\n")

for atom in atoms:
    NumberBF = 0
    for BF in atom[2]:
        if   BF[0] == 'S' : NumberBF = NumberBF + 1
        elif BF[0] == 'P' : NumberBF = NumberBF + 3
        elif BF[0] == 'D' : NumberBF = NumberBF + 6
        elif BF[0] == 'F' : NumberBF = NumberBF + 10
        else : print "ERROR: Improper basis function type!"

    OUT.write("%s\t%i\t%s\n" % ( atom[0], NumberBF, atom[1] ) )

    for BF in atom[2]:
        if BF[0] == 'S' :
            OUT.write("\t%i\t%s\n" % (len(BF[2]),'s'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize(0,0,0,C[1],C[0])))
                
        if BF[0] == 'P' :
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'px'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize(1,0,0,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'py'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize(0,1,0,C[1],C[0])))
            OUT.write("\t%i\t%s\n" % (len(BF[2]),'pz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize(0,0,1,C[1],C[0])))
                
        if BF[0] == 'D' :
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'dxx'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize(2,0,0,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'dyy'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize(0,2,0,C[1],C[0])))
            OUT.write("\t%i\t%s\n" % (len(BF[2]),'dzz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize(0,0,2,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'dxy'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize(1,1,0,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'dxz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize(1,0,1,C[1],C[0])))
            OUT.write("\t%i\t%s\n" % (len(BF[2]),'dyz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize(0,1,1,C[1],C[0]))) 

        if BF[0] == 'F' :
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fxxx'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize(3,0,0,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fyyy'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize(0,3,0,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fzzz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize(0,0,3,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fyyx'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize(1,2,0,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fxxy'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize(2,1,0,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fxxz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize(2,0,1,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fzzx'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize(1,0,2,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fzzy'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize(0,1,2,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fyyz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize(0,2,1,C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fxyz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize(1,1,1,C[1],C[0])))
                
OUT.write("&\n")

# End Write Basis Set

# Write Wavefunction
nalpha = 0
nbeta  = 0
OUT.write("&wavefunction\n")
for i in range(len(Wavefunction)):
    # Get Occupation
    if (string.atof(Wavefunction[i][0][0]) == 1.0):
        nalpha = nalpha + 1
        nbeta  = nbeta + 1
    elif (string.atof(Wavefunction[i][0][0]) == 0.5):
        nalpha = nalpha + 1
    elif (string.atof(Wavefunction[i][0][0]) == 0.0):
        continue
    else:
        OUT.write("ERROR\n\t")

    # Write Orbital
for i in range(len(Wavefunction)):
    for j in range(len(Wavefunction[i][1])):
        OUT.write("\t%s" % (Wavefunction[i][1][j]))
        if( j%3 == 2 ): OUT.write("\n")
    OUT.write("\n\n")

OUT.write("Alpha Occupation\n")
for i in range(nalpha):
    OUT.write("1\t")
for i in range(len(Wavefunction)-nalpha):
    OUT.write("0\t")
OUT.write("\n\n")

OUT.write("Beta Occupation\n")
for i in range(nbeta):
    OUT.write("1\t")
for i in range(len(Wavefunction)-nbeta):
    OUT.write("0\t")
OUT.write("\n\n")

OUT.write("CI Coeffs\n")
OUT.write("1.0\n")
OUT.write("\n\n")

OUT.write("&\n")

# End Write Wavefunction

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




