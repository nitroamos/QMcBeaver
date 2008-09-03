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
import time
import os
from atomic_symbol_to_Z import *
from utilities import *

if len(sys.argv) < 2: 
    print "jaguar2qmcbeaver.py <jaguar restart>.01.in"
    sys.exit(0)

Infile=sys.argv[1]
if string.find(Infile,'.in') == -1:
    #check for .01.in ??
    print "You need to use a Jaguar restart file, but you used ", Infile, ".\n"
    sys.exit(0)
    
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
trialFunctionType="restricted"
for j in range(len(jagoutdata)):
    if string.find(jagoutdata[j],"SCFE") != -1 and string.find(jagoutdata[j],"GVB") != -1:
        trialFunctionType = "gvb"
    if string.find(jagoutdata[j],"SCFE") != -1:
        tempspot = j
    if string.find(jagoutdata[j],"ERROR") != -1:
        Error = Error + 1
if tempspot < 0:
    print "Could Not Find Energy!"
energyline=string.split(jagoutdata[tempspot])
Energy=energyline[4]

for i in range(len(jaguardata)):
    if string.find(jaguardata[i],"iuhf=1") == 0:
        trialFunctionType = "unrestricted"
        break
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

for i in range(len(jaguardata)):
    if string.find(jaguardata[i],"iunit") != -1:
	if string.find(jaguardata[i],"0") != -1 or string.find(jaguardata[i],"2") != -1:
	    #units are already bohr, so we don't need to convert
	    a0 = 1.0
	    break

for i in range(len(Qmcdata)):
    Qmcdata[i] = [ Qmcdata[i][0], str(symbol_to_Z(Qmcdata[i][0])) ] \
                 + [string.atof(Qmcdata[i][1])/a0,string.atof(Qmcdata[i][2]) \
                    /a0,string.atof(Qmcdata[i][3])/a0] 

geometry = Qmcdata
# End Get Geometry

# Get Wavefunction
Wavefunction = []
AlphaOrbitals = []
BetaOrbitals = []
nalpha = 0
nbeta  = 0

if trialFunctionType == "restricted" or trialFunctionType == "gvb":
    for i in range(len(jaguardata)):
        if string.find(jaguardata[i],"&guess") == 0:
            Start = i+1
            break
    for i in range(Start,len(jaguardata)):
        if jaguardata[i][0] == '&':
            End = i
            break
    LongLine=[]
    Occupation=-1
    
    for i in range(Start,End):
        Line=string.split(jaguardata[i])
        Index=string.find(jaguardata[i],"Occupation")
        if Index != -1:
            if i != Start and string.atof(Occupation) > 0.0:
		#we found a new orbital, so that means the old one was completed
		Wavefunction=Wavefunction + [[[Occupation],LongLine]]
            Occupation=Line[5]
            LongLine=[]
            continue
        LongLine=LongLine+Line
    #add the last orbital
    if string.atof(Occupation) > 0.0:
	Wavefunction=Wavefunction + [[[Occupation],LongLine]]

    norbitals=len(Wavefunction)

    for i in range(norbitals):
	# Get Occupation
	if (string.atof(Wavefunction[i][0][0]) == 1.0):
	    nalpha = nalpha + 1
	    nbeta  = nbeta + 1
	elif (string.atof(Wavefunction[i][0][0]) == 0.5):
	    nalpha = nalpha + 1
	elif (string.atof(Wavefunction[i][0][0]) == 0.0):
	    continue

    ndeterminants=1
    CI_coeffs = [1]

    AlphaOccupation = [0]
    BetaOccupation = [0]
    AlphaOccupation[0] = range(norbitals)
    BetaOccupation[0] = range(norbitals)
    
    for i in range(nalpha):
	AlphaOccupation[0][i] = 1
    for j in range(nalpha,norbitals):
	AlphaOccupation[0][j] = 0

    for i in range(nbeta):
	BetaOccupation[0][i] = 1
    for j in range(nbeta,norbitals):
	BetaOccupation[0][j] = 0

elif trialFunctionType == "unrestricted":
    for i in range(len(jaguardata)):
        if string.find(jaguardata[i],"Alpha Molecular Orbitals") == 0:
            start_alpha = i+1
            break
    for i in range(start_alpha,len(jaguardata)):
        if string.find(jaguardata[i],"Beta Molecular Orbitals") == 0:
            end_alpha = i
            start_beta = i+1
            break
    for i in range(start_beta,len(jaguardata)):
        if jaguardata[i][0] == '&':
            end_beta = i
            break

    LongLine = []
    Occupation = -1
    for i in range(start_alpha,end_alpha):
        Line = string.split(jaguardata[i])
        Index = string.find(jaguardata[i],"Occupation")
        if Index != -1:
            if i!= 0:
                AlphaOrbitals = AlphaOrbitals + [[[Occupation],LongLine]]
            Occupation = Line[6]
            LongLine = []
            continue
        LongLine = LongLine + Line
    AlphaOrbitals = AlphaOrbitals + [[[Occupation],LongLine]]
    AlphaOrbitals = AlphaOrbitals[1:]

    LongLine = []
    Occupation = -1
    for i in range(start_beta,end_beta):
        Line = string.split(jaguardata[i])
        Index = string.find(jaguardata[i],"Occupation")
        if Index != -1:
            if i!= 0:
                BetaOrbitals = BetaOrbitals + [[[Occupation],LongLine]]
            Occupation = Line[6]
            LongLine = []
            continue
        LongLine = LongLine + Line
    BetaOrbitals = BetaOrbitals + [[[Occupation],LongLine]]
    BetaOrbitals = BetaOrbitals[1:]

    for i in range(len(AlphaOrbitals)):
        if (string.atof(AlphaOrbitals[i][0][0]) == 1.0):
            nalpha = nalpha + 1
        elif (string.atof(AlphaOrbitals[i][0][0]) == 0.0):
            continue
        else:
            OUT.write("ERROR in getting Alpha Occupation\n\t")
    for i in range(len(BetaOrbitals)):
        if (string.atof(BetaOrbitals[i][0][0]) == 1.0):
            nbeta = nbeta + 1
        elif (string.atof(BetaOrbitals[i][0][0]) == 0.0):
            continue
        else:
            OUT.write("ERROR in getting Beta Occupation\n\t")


    ndeterminants=1
    CI_coeffs = [1]

    AlphaOccupation = [0]
    BetaOccupation = [0]
    AlphaOccupation[0] = range(len(AlphaOrbitals))
    BetaOccupation[0] = range(len(BetaOrbitals))

    for i in range(nalpha):
        AlphaOccupation[0][i] = 1
    for j in range(len(AlphaOrbitals)-nalpha): 
        AlphaOccupation[0][j] = 0

    for i in range(nbeta):
        BetaOccupation[0][i] = 1
    for j in range(len(BetaOrbitals)-nbeta):
        BetaOccupation[0][j] = 0

    norbitals = len(AlphaOrbitals)
else:
    print "Unknown trial function type: ", trialFunctionType
    sys.exit(0)

if trialFunctionType == "gvb":
    ncore = 0
    for i in range(norbitals):
        if abs(string.atof(Wavefunction[i][0][0])-1.0) < 1e-15:
            ncore += 1
    
    for j in range(len(jagoutdata)):
        if string.find(jagoutdata[j],'first natural orbital') != -1:
            start_ci = j+4
            break
    
    for i in range(start_ci,len(jagoutdata)):
	ci_line = string.split(jagoutdata[i])
	if len(ci_line) != 11:
	    end_ci = i
	    break

    nalpha = ncore + (end_ci-start_ci)
    nbeta  = ncore + (end_ci-start_ci)
    AlphaOccupation = range(1) 
    CI_coeffs = range(1)
    CI_coeffs[0] = 1.0
   
    for i in range(len(AlphaOccupation)):
        AlphaOccupation[i] = range(norbitals)

    for i in range(len(AlphaOccupation)):
        for j in range(norbitals):
	    if j < ncore:
		AlphaOccupation[i][j] = 1
	    else:
		AlphaOccupation[i][j] = 0
    
    # we need to expand out the geminal pairs into separate determinants
    # the start in the form (c1 + c2) * (c3 + c4) * (c5 + c6)...
    # expanding to (c1*c3 + c1*c4 + c2*c3 + c2*c4)*(c5 * c6)*(...)
    # I've implemented a recursion in a loop.
    print "core=",ncore,"start=",start_ci, " end=", end_ci
    index = 3
    idet = 0
    for p in range(end_ci-start_ci):
	old = len(CI_coeffs)
	ci_line = string.split(jagoutdata[p+start_ci])
	coef1 = string.atof(ci_line[4])
	index += 1
	coef2 = string.atof(ci_line[8])
	index += 3
	
	orb1 = string.atoi(ci_line[1])-1
	orb2 = string.atoi(ci_line[5])-1

	pairCI    = range(old*2)
	pairAlpha = range(old*2)
	for i in range(len(pairAlpha)):
	    pairAlpha[i] = range(norbitals)
	    
	for orb in range(len(pairCI)):
	    for i in range(len(AlphaOccupation[orb%old])):
		pairAlpha[orb][i] = AlphaOccupation[orb%old][i]
	    
	    if orb < old:
		pairCI[orb] = CI_coeffs[orb%old] * coef1
		pairAlpha[orb][orb1] = 1
	    else:
		pairCI[orb] = CI_coeffs[orb%old] * coef2
		pairAlpha[orb][orb2] = 1

	print "Adding geminal pair %i (%i %20.10f / %i %20.10f)"%(p+1,orb1,coef1,orb2,coef2)

	#print "Determinant CI coefficients = ",pairCI	

	CI_coeffs = pairCI
	AlphaOccupation = pairAlpha
    BetaOccupation = AlphaOccupation
    ndeterminants = pow(2,(end_ci - start_ci))
    assert(ndeterminants == len(CI_coeffs))

# End Get Wavefunction

# Get Basis Set
try:
    IN = open(JagBasisfile,'r')
except:
    print "Did you run the calculation with ip164=2 to get the MQM .bas file?\n"
    sys.exit(0)

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
##################  PRINT FLAGS: BEGIN    ########################
#
# We'll use ckmf template files (extension ckmft) to help make
# our ckmf file. The reason is to save us from having to go through
# and manually choose all our parameters. These template files are meant
# to be pretty close to what we'll end up wanting.
my_path, my_name = os.path.split(__file__)

#A couple default placed to look for "ckmft" files
templatedir = [".","..","../..","../examples","/ul/amosa/ckmf_origs",my_path]

templates = []
for dir in templatedir:
    if os.path.exists(dir):
	for file in os.listdir(dir):
	    if file[-5:] == "ckmft":
		templates.append(dir+"/"+file)

print "\nAvailable ckmf templates:"
for i in range(len(templates)):
    print " %3i : " % i, templates[i]

try:
    choice = string.atoi(raw_input("Your choice [0]:"))
except:
    choice = 0

print "You chose ckmf template: ", templates[choice]
myStandardFlags=open(templates[choice],'r')

OUT.write('# Created on %s\n'%(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
OUT.write('# Using Jaguar restart file:  %s\n'% os.path.abspath(Infile))
OUT.write('# Using ckmft template file: %s\n'% templates[choice])

#let's save some of the important info from a GAMESS calculation
for i in range(len(jagoutdata)):
    line = -1;
    if string.find(jagoutdata[i],'SCFE: SCF energy:') != -1:
	    line = i
    if string.find(jagoutdata[i],'Molecular Point Group:') != -1:
	    line = i
    if string.find(jagoutdata[i],'Stoichiometry:') != -1:
	    line = i
    if string.find(jagoutdata[i],'basis set:') != -1:
	    line = i
    if line > 0:
	    OUT.write('#%s' % jagoutdata[line])
	    print '#%s' % jagoutdata[line],

	    
OUT.write("\n")
OUT.write(myStandardFlags.read());
myStandardFlags.close()

# Write out data
OUT.write("atoms\n %i\n"%(Natoms))
OUT.write("charge\n %i\n"%(Charge))
OUT.write("energy\n %s\n"%(Energy))

if trialFunctionType == "restricted":
    OUT.write("norbitals\n %i\n"%(norbitals))
    OUT.write("nbasisfunc\n %i\n"%(len(Wavefunction[0][1])))
    OUT.write("trial_function_type\n %s\n"%(trialFunctionType))
elif trialFunctionType == "unrestricted":
    OUT.write("norbitals\n %i\n"%(len(AlphaOrbitals)))
    OUT.write("nbasisfunc\n %i\n"%(len(AlphaOrbitals[0][1])))
    OUT.write("trial_function_type\n %s\n"%(trialFunctionType))
elif trialFunctionType == "gvb":
    OUT.write("norbitals\n %i\n"%(norbitals))
    OUT.write("nbasisfunc\n %i\n"%(len(Wavefunction[0][1])))
    OUT.write("trial_function_type\n %s\n"%("restricted"))
    
OUT.write("ndeterminants\n %i\n"%(ndeterminants))
OUT.write("chip_and_mike_are_cool\n Yea_Baby!\n")
OUT.write("&\n")

# Write Geometry
OUT.write("&geometry\n")
for i in range(len(Qmcdata)):
    OUT.write("%s %5s"%(str(Qmcdata[i][0]),str(Qmcdata[i][1])))
    for j in range(2,len(Qmcdata[i])):
        OUT.write("%20s"%str(Qmcdata[i][j]))
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
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize("s",C[1],C[0])))
                
	elif BF[0] == 'P' :
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'px'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize("px",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'py'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize("py",C[1],C[0])))
            OUT.write("\t%i\t%s\n" % (len(BF[2]),'pz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize("pz",C[1],C[0])))
                
	elif BF[0] == 'D' :
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'dxx'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize("dxx",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'dyy'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize("dyy",C[1],C[0])))
            OUT.write("\t%i\t%s\n" % (len(BF[2]),'dzz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize("dzz",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'dxy'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize("dxy",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'dxz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize("dxz",C[1],C[0])))
            OUT.write("\t%i\t%s\n" % (len(BF[2]),'dyz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n" %(C[0], normalize("dyz",C[1],C[0]))) 

	elif BF[0] == 'F' :
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fxxx'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize("fxxx",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fyyy'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize("fyyy",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fzzz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize("fzzz",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fyyx'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize("fyyz",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fxxy'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize("fxxy",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fxxz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize("fxxz",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fzzx'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize("fzzx",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fzzy'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize("fyzz",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fyyz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize("fyyz",C[1],C[0])))
            OUT.write("\t%i\t%s\n" %(len(BF[2]),'fxyz'))
            for C in BF[2]:
                OUT.write("\t\t%s\t%s\n"%(C[0], normalize("fxyz",C[1],C[0])))

	else:
	    print "Unknown basis function type: ",BF[0]
OUT.write("&\n")

# End Write Basis Set

# Write Wavefunction
OUT.write("&wavefunction\n")

if trialFunctionType == "restricted" or trialFunctionType == "gvb":
    # Write Orbital
    for i in range(norbitals):
        for j in range(len(Wavefunction[i][1])):
            OUT.write("\t%s" % (Wavefunction[i][1][j]))
            if( j%3 == 2 ): OUT.write("\n")
        OUT.write("\n\n")

elif trialFunctionType == "unrestricted":            
    # Write Alpha Orbitals
    OUT.write("Alpha Molecular Orbitals\n\n")
    for i in range(len(AlphaOrbitals)):
        for j in range(len(AlphaOrbitals[i][1])):
            OUT.write("\t%s" % (AlphaOrbitals[i][1][j]))
            if( j%3 == 2 ): OUT.write("\n")
        OUT.write("\n\n")

    # Write Beta Orbitals
    OUT.write("Beta Molecular Orbitals\n\n")
    for i in range(len(BetaOrbitals)):
        for j in range(len(BetaOrbitals[i][1])):
            OUT.write("\t%s" % (BetaOrbitals[i][1][j]))
            if( j%3 == 2 ): OUT.write("\n")
        OUT.write("\n\n")    
    
# Write CI Coefficients
OUT.write("Alpha Occupation\n")
for i in range(ndeterminants):
    nume = 0
    for j in range(norbitals):
        OUT.write('%i '%AlphaOccupation[i][j])
	print '%i '%AlphaOccupation[i][j],
	nume += string.atof(AlphaOccupation[i][j])
    OUT.write('\n')
    print ""
    if abs(nalpha-nume) > 1e-10:
	print "nalpha ",nalpha, " does not match electron count ",nume
OUT.write('\n')

OUT.write("Beta Occupation\n")
for i in range(ndeterminants):
    nume = 0
    for j in range(norbitals):
        OUT.write('%i '%BetaOccupation[i][j])
	nume += string.atof(BetaOccupation[i][j]) 
    if abs(nbeta-nume) > 1e-10:
	print "nbeta ",nbeta, " does not match electron count ",nume
    OUT.write('\n')
OUT.write('\n')

OUT.write("CI Coeffs\n")
for i in range(ndeterminants):
    OUT.write('%25s\n'%CI_coeffs[i])
OUT.write('\n')

print str(ndeterminants) + " CI determinant(s) used, with coefficients:"
cum = 0
for i in range(ndeterminants):
    ci = string.atof(CI_coeffs[i])
    cum += ci*ci
    print "%3i) %15.7f has percentage %15.10e, cumulative remaining %15.10e" % (i+1,ci,ci*ci,1.0-cum)

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


OUT.write('&Jastrow\n')

################## PRINT JASTROW: END ################################




