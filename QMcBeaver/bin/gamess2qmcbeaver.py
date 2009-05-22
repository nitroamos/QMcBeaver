#!/usr/bin/env python

# This script will convert the output from a GAMESS calculation into an input
# file for QMcBeaver.
# * It will copy all the basis function data and orbitals
# * It will look for the energies calculated in GAMESS
#   and add them as comments.
# * To find a good set of QMcBeaver flags, it will look for a "ckmft" file
#   in a few directories (see "templatedir" variable below)
#   to copy a good set of defaults.
#
# Usage:
# gamess2qmcbeaver.py <GAMESS output file> [determinant cutoff = 0.0]
# * We recognize .log and .inp.out as GAMESS output extensions.
# * If the absolute value of the CI coefficient is below the determinant
#   cutoff, then it will not be included in the ckmf file.
#
# Permissible RUNTYP = ENERGY and OPTIMIZE
# Permissible SCFTYP = anything other than MCSCF
#
# To use SCFTYP=MCSCF:
# 1) Run the MCSCF calculation in GAMESS. This is the hardest part... Look in
#    the GAMESS manual and the "Further Information" document for hints.
# 2) Make a 2nd GAMESS input file with a $VEC section from the natural orbitals
#    and minimized geometry of the MCSCF run. Specify SCFTYP=NONE and CITYP=ALDET.
#    This will produce a CI expansion in these orbitals. You might be able to use
#    other CITYP, but we haven't programmed them.
# 3) This script can read the ALDET output file, and will find as many determinants
#    as were printed out, and put them in the ckmf file. You might need to modify
#    PRTTOL in the $DET section to get more determinants.
#
# NOTE: check your ALDET runs... I've found that the occupations don't always match
# the orbitals printed! For one of my runs, it sorted the natural orbitals according to occupation,
# which was different from the input order.

import re
import sys
import copy
import math
import string
import time
import os
from utilities import *

if len(sys.argv) < 2: 
    print "gamess2qmcbeaver.py <filename>[.log, .inp.out] [detcutoff=0.0]"
    sys.exit(0)
    
Infile = sys.argv[1]
IN = open(Infile,'r')
gamess_output = IN.readlines()
IN.close()

filebase = ""
if string.find(Infile,'.inp.out') != -1:
	filebase = sys.argv[1][0:len(sys.argv[1])-7]
elif string.find(Infile,'.log') != -1:
	filebase = sys.argv[1][0:len(sys.argv[1])-3]
else:
	print "The file ", Infile, " is not recognized as a GAMESS log file!"
	sys.exit(0)
	
Datafile = filebase + "dat"	
IN2 = open(Datafile,'r')
gamess_data = IN2.readlines()
IN2.close()

Outfile = filebase + "ckmf"
OUT = open(Outfile,'w')

run_type  = "ENERGY"
scf_type  = "RHF"
ci_type   = "NONE"
pp_type   = "NONE"
spin_mult = 1
istate    = 1

detcutoff = 1e-10
if len(sys.argv) == 3:
    detcutoff = string.atof(sys.argv[2])
    print "Removing all determinants with coefficients less than ",detcutoff
	
# Get run type and scf type

for i in range(len(gamess_output)):
    if string.find(gamess_output[i],'$CONTRL OPTIONS') != -1:
	k = i
	while string.find(gamess_output[k],'$SYSTEM OPTIONS') == -1:
	    line = re.split('[\s=]+',gamess_output[k])
	    for j in range(len(line)):
		if string.find(line[j],'SCFTYP') != -1:
		    scf_type = line[j+1]
		if string.find(line[j],'VBTYP') != -1:
		    if string.find(line[j+1],'NONE') == -1:
			scf_type = line[j+1]
		if string.find(line[j],'RUNTYP') != -1:
		    run_type = line[j+1]
		if string.find(line[j],'CITYP') != -1:
		    ci_type = line[j+1]
		if string.find(line[j],'MULT') != -1:
		    spin_mult = string.atoi(line[j+1])
		if string.find(line[j],'PP') != -1:
		    pp_type = line[j+1]
	    k += 1

if ci_type == "GENCI":
    #These are effectively the same kind of calculation
    #Just different lists of determinants
    ci_type = "ALDET"

#################### EXTRACT GEOMETRY: START ####################

# Find where the geometry is stored.

if run_type == "ENERGY" or run_type == "HESSIAN":
    for i in range(len(gamess_output)):
        if string.find(gamess_output[i], 'RUN TITLE') != -1:
            start_geometry = i
        if string.find(gamess_output[i], 'INTERNUCLEAR DISTANCES') != -1:
            end_geometry = i
            break

elif run_type == "OPTIMIZE":
    for i in range(len(gamess_output)):
        if string.find(gamess_output[i],'EQUILIBRIUM GEOMETRY LOCATED') != -1:
            start_geometry = i
            for j in range(i,len(gamess_output)):
                if string.find(gamess_output[j],'INTERNUCLEAR DISTANCES') !=-1:
                    end_geometry = j-1
                    break
                elif string.find(gamess_output[j],'INTERNAL COORDINATES') !=-1:
                    end_geometry = j-3
                    break
                elif string.find(gamess_output[j],'SUBSTITUTED Z-MATRIX') !=-1:
                    end_geometry = j-1
                    break
            break    

try:
    geom_data = gamess_output[start_geometry:end_geometry]
except:
    print "Failed to find geometry for run_type = ", run_type
    raise
geometry = []

start = 0

if run_type == "ENERGY":
    for line in geom_data:
        if start: geometry = geometry + [line]
        if string.find(line,'CHARGE') != -1: start = 1
    geometry = geometry[:len(geometry)-1]

elif run_type == "OPTIMIZE":
    for line in geom_data:
        if start == 2: geometry = geometry + [line]
        if string.find(line,'CHARGE') != -1: start = start + 1
    geometry = geometry[1:]

#split up the data
for i in range(len(geometry)):
    geometry[i] = string.split(geometry[i])
    for j in range(2,5):
        geometry[i][j] = string.atof(geometry[i][j])

#convert from ANGs to BOHR if necessary

ANGtoBOHRconversion = 1.0/0.529177249

for line in geom_data:
    if string.find(line,'(ANGS)') != -1:
        for i in range(len(geometry)):
            for j in range(2,5):
                geometry[i][j] = geometry[i][j] * ANGtoBOHRconversion
        break

#################### EXTRACT GEOMETRY: END ######################

#################### EXTRACT BASIS SET: BEGIN ###################

start_basis = 0
end_basis = 0
for i in range(len(gamess_output)):
    if string.find(gamess_output[i], 'ATOMIC BASIS SET') != -1:
        start_basis = i
    if string.find(gamess_output[i], '$CONTRL OPTIONS') != -1:
        end_basis = i
        break
basisdata = gamess_output[start_basis:end_basis]

end = 0
for i in range(len(basisdata)):
    if string.find(basisdata[i],'TOTAL NUMBER OF SHELLS') != -1 :
        end = i
        break
    if string.find(basisdata[i],'TOTAL NUMBER OF BASIS SET SHELLS') != -1 :
        end = i
        break
basisdata = basisdata[7:end]

basis = []
atom = []
bf = []
atomnumber = -1
bfnumber = 0
for line in basisdata:
    if line != '\n':
        line = string.replace(line,')',' ')
        line = string.replace(line,'(',' ')
        line = string.split(line)
        if len(line) == 1:
            # We are starting a new atom
            if bf != []:
                # We add the old contracted basis function to the atom and clear the temp
                atom = atom + [bf]
                bf = []
            if atom != []:
                # We add the old atom to the basis and clear the temp space
                basis = basis + [atom]
                atom = []
            # We start the new atom with the label
            atom = atom + [line]
            
        elif len(line) > 1 and line[0] != str(bfnumber):
            bfnumber = string.atoi(line[0])
            # We are starting a new contracted basis function for this atom
            if bf != []:
                # We add the old contracted basis function to the atom and clear the temp
                atom = atom + [bf]
                bf = []
            # We start the new contracted basis function
            temp = [line[1]] + line[3:]
            if len(line) > 6:
                temp = temp + [line[6]]
            line = temp
            bf = [line]

        elif len(line) > 1 and line[0] == str(bfnumber):
            # We are continuing to add primitive basis functions to the contracted one
            temp = [line[1]] + line[3:]
            if len(line) > 6:
                temp = temp + [line[6]]
            line = temp
            bf = bf + [line]
atom = atom + [bf]
basis = basis + [atom]

# extract some flags data from this section
for line in gamess_output:
    if string.find(line,'TOTAL NUMBER OF BASIS FUNCTIONS') !=-1 :
        nbasisfunc = string.atoi(string.split(line)[6])
    if string.find(line,'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS') !=-1 :
        nbasisfunc = string.atoi(string.split(line)[7])
    if string.find(line,'CHARGE OF MOLECULE') !=-1 :
        charge = string.atoi(string.split(line)[4])
#we'll rely on orbital occupations to indicate charge
#	charge = 0
    if string.find(line,'TOTAL NUMBER OF ATOMS') !=-1 :
        atoms = string.atoi(string.split(line)[5])
    if string.find(line,'NUMBER OF OCCUPIED ORBITALS (ALPHA)') != -1 :
	nalpha = string.atoi(line.split('=')[1])
    if string.find(line,'NUMBER OF OCCUPIED ORBITALS (BETA )') != -1 :
	nbeta = string.atoi(line.split('=')[1])

#################### EXTRACT BASIS SET: END #######################

#################### EXTRACT WAVEFUNCTION: BEGIN ##################


# First, we want to load in all the $VEC .. $END sections we can find.
# We look in both the .dat file, and in the .inp file, and we save
# the name that GAMESS gave it.
collecting = 0
Inputfile = filebase + "inp"	
INPFILE = open(Inputfile,'r')
input_data = INPFILE.readlines()
INPFILE.close()

name = "MOREAD orbitals from " + Inputfile + "\n"
raw_orbitals = []
orbital_vecs = []
orbital_name = []
iorder_vecs  = []
norder_vecs  = 0
for i in range(len(input_data)):
    m = re.search('norder\s*=(\d+)',input_data[i],re.I)
    if m:
	norder_vecs = int(m.group(1))
    m = re.search('iorder\((\d+)\)=([\d,]+)',input_data[i],re.I)
    if m:
	iorder_vecs.append(int(m.group(1)))
	iorder_vecs += [int(k) for k in m.group(2).split(',')]
	
    if string.find(input_data[i],'$END') != -1 and collecting == 1:
	collecting = 0
	orbital_vecs = orbital_vecs + [raw_orbitals]
	orbital_name = orbital_name + [name]
	raw_orbitals = []
    if string.find(input_data[i],'$VEC') != -1:
	collecting = 1
    elif collecting == 1:
	raw_orbitals = raw_orbitals + [input_data[i]]

if norder_vecs == 1:
    print "Notice: found IORDER section for MOREAD orbitals: ",iorder_vecs
       
for i in range(len(gamess_data)):
    if string.find(gamess_data[i],'NO-S OF CI STATE') != -1 or \
	   string.find(gamess_data[i], 'GVB ORBITALS') != -1 or \
	   string.find(gamess_data[i], 'LOCALIZED') != -1 or \
	   string.find(gamess_data[i], 'OPEN SHELL ORBITALS') != -1 or \
	   string.find(gamess_data[i], 'CLOSED SHELL ORBITALS') != -1 or \
	   string.find(gamess_data[i], 'MP2 NATURAL ORBITALS') != -1 or \
	   string.find(gamess_data[i], 'OPTIMIZED MCSCF') != -1 or \
	   string.find(gamess_data[i], 'NATURAL ORBITALS OF MCSCF') != -1:
	name = gamess_data[i]
    if string.find(gamess_data[i],'$END') != -1 and collecting == 1:
	collecting = 0
	orbital_vecs = orbital_vecs + [raw_orbitals]
	orbital_name = orbital_name + [name]
	raw_orbitals = []
    if string.find(gamess_data[i],'$VEC') != -1:
	collecting = 1
    elif collecting == 1:
	raw_orbitals = raw_orbitals + [gamess_data[i]]
    m=re.search('^E\([\w\-]+\)=\s*([\d\-\.]+)',gamess_data[i])
    if m:
	energy = m.group(1)
    m=re.search('CI STATE\s+\d+\sE=\s*([\d\-\.]+)',gamess_data[i])
    if m:
	energy = m.group(1)
	
if scf_type == "RHF":
    default_orb_string = 'CLOSED SHELL ORBITALS'    
elif scf_type == "ROHF":
    default_orb_string = 'OPEN SHELL ORBITALS'
elif scf_type == "UHF":
    default_orb_string = ''
elif scf_type == "GVB":
    default_orb_string = 'GVB ORBITALS'
elif scf_type == "NONE" and ci_type == "ALDET":
    default_orb_string = 'MOREAD'
elif scf_type == "VB2000":
    cicoef=[]
    for i in range(len(gamess_output)):
        if string.find(gamess_output[i],'Normalized structure coefficients') != -1:
	    line = gamess_output[i+1]
            cicoef += string.split(line)	    
        if string.find(gamess_output[i],'ENERGY AND DIFF OF MACROITER') != -1:
	    line = gamess_output[i]
            energy = (string.split(line))[7]

    pcum = 0
    for i in range(len(cicoef)):
	cicoef[i] = string.atof(cicoef[i])
	pcum += cicoef[i]*cicoef[i]
    print "VB Coeff = ",cicoef, "\nnorm = ",pcum
    print "VB Energy = ",energy

    print "Fix VB2000 to get the right orbitals!";
    sys.exit(0)
    Datafile = filebase + "vec"	
    IN2 = open(Datafile,'r')
    gamess_data = IN2.readlines()
    IN2.close()
else:
    print "SCFTYP", scf_type, "is not supported."
    sys.exit(0) 

orb_choice = len(orbital_name)-1
print "We found %i orbital sets:"%len(orbital_name)
for i in range(len(orbital_name)):
    lines_per_orb = int(nbasisfunc/5)+1
    print "%i) %g orbitals"%(i,((len(orbital_vecs[i]))/float(lines_per_orb))),
    print "for: ", orbital_name[i],
    if string.find(orbital_name[i],default_orb_string) != -1:
	orb_choice = i
try:
    orb_choice = string.atoi(raw_input("Your choice [%i]:"%orb_choice))
except:
    print "",

# Get the wavefunction parameters from the .dat file.
orbital_number = []
orbital_coeffs = [] 
wavefunction   = []
current_index  = 1
for n in range(len(orbital_vecs[orb_choice])):
    orbital_index = string.atoi(orbital_vecs[orb_choice][n][0:2])

    #if the index changed, then we completed the old orbital, so we add it to the wf
    if orbital_index != current_index:
	wavefunction.append(orbital_coeffs)
	current_index = orbital_index
	orbital_coeffs = []

    #turn the text into numbers: index coeff1 coeff2 coeff3 coeff4 coeff5
    len_line = len(orbital_vecs[orb_choice][n])
    number_of_entries = len_line/15
    line_data = range(number_of_entries)
    for i in range(number_of_entries):
        line_data[number_of_entries-i-1] = \
                   orbital_vecs[orb_choice][n][len_line-15*(i+1)-1:len_line-15*i-1]
    for j in range(len(line_data)):
        line_data[j] = string.atof(line_data[j])

    #append the current line to the current orbital
    orbital_coeffs = orbital_coeffs + line_data

#Add that last orbital
wavefunction.append(orbital_coeffs)    
norbitals = len(wavefunction)

if re.search("MOREAD",orbital_name[orb_choice],re.I) and len(iorder_vecs) > 0 and norder_vecs == 1:
    print "Reordering according to IORDER: ",
    new_orbitals = []
    num = len(iorder_vecs)-1
    print num, ", first = ",iorder_vecs[0]
    for i in range(iorder_vecs[0]-1):
	print i+1,
	new_orbitals.append(wavefunction[i])
    for i in range(num):
	print iorder_vecs[i+1],
	new_orbitals.append(wavefunction[iorder_vecs[i+1]-1])
    for i in range(iorder_vecs[0]-1+num,norbitals):
	print i+1,
	new_orbitals.append(wavefunction[i])
    print "\n"
    wavefunction = new_orbitals

################ Set up the occupation and CI coefficent arrays.
if scf_type == "RHF" or scf_type == "ROHF" or scf_type == "UHF":
    AlphaOcc    = [0]
    BetaOcc     = [0]
    AlphaOcc[0] = range(norbitals)
    BetaOcc[0]  = range(norbitals)

    for j in range(0,norbitals):
        AlphaOcc[0][j] = 0
        BetaOcc[0][j]  = 0
    
    for i in range(nalpha):
        AlphaOcc[0][i] = 1

    #The VEC section is double in size... The first half are the alpha
    #orbitals, and the second half are the beta orbitals... If that statement
    #isn't always true, then this will have problems.
    beta_start = 0
    if scf_type == "UHF":
	if norbitals % 2 == 0:
	    beta_start = norbitals/2
	else:
	    print "Error: unexpected problem reading UHF wavefunction...\n"
	    sys.exit(0)

    for i in range(beta_start,nbeta+beta_start):
        BetaOcc[0][i] = 1

    ncore = 0
    ndeterminants = 1
    CI = [1]

elif scf_type == "GVB":
    #the CI Coefficients in the dat file are more precise, so we want them
    #this might have a problem if too many GVB pairs are used
    cicoef=[]
    for i in range(len(gamess_data)):
        if string.find(gamess_data[i],'CICOEF') != -1:
	    line = gamess_data[i]
	    p = re.compile("\(\s+")
	    line = p.sub("(",line)
	    line = string.replace(line,',',' ')
            cicoef += string.split(line)	    

    core_line_number = -1
    start_ci = 0
    end_ci = 0
    for i in range(len(gamess_output)):
        if string.find(gamess_output[i],'ROHF-GVB INPUT PARAMETERS') != -1:
            core_line_number = i+3
            core_line = string.split(gamess_output[core_line_number])
	    norb = string.atoi(core_line[2])
            ncore = string.atoi(core_line[5])
	    pair_line = string.split(gamess_output[core_line_number+1]) 
	    npair = string.atoi(pair_line[2])
	    nseto = string.atoi(pair_line[5])
	    odegen = 0
	    if re.search("NO",gamess_output[core_line_number+2],re.I):
		no_line = string.split(gamess_output[core_line_number+2])
		odegen  = string.atoi(no_line[2])
		
	    print "GVB settings: mult=",spin_mult,"ncore=",ncore,"norb=",norb,"npair=",npair,"nseto=",nseto

	    #if odegen > 0:
	#	print "Error: open shell too complicated, no = ",odegen
	#	sys.exit(0)
        #    break

    AlphaOcc = range(1)
    CI = range(1)
    CI[0] = 1.0
    
    for i in range(len(AlphaOcc)):
        AlphaOcc[i] = range(norbitals)

    for i in range(len(AlphaOcc)):
        for j in range(norbitals):
	    if j < ncore:
		AlphaOcc[i][j] = 1
	    else:
		AlphaOcc[i][j] = 0

    if npair == 0:
	ndeterminants = 1
	ncore = 0
	
    if npair > 0:
	#The two perfect paired electrons are spin coupled into a singlet
	# See Eq 56 from "SCF Equations for GVB" by Bobrowicz and Goddard
	# WF = anti[(c1 11 - c2 22)ab] = c1 anti[11ab] - c2 anti[22ab]

	for j in range(core_line_number,len(gamess_output)):
	    if string.find(gamess_output[j],'CI COEFFICIENTS') != -1:
		start_ci = j+2
		break

	for i in range(start_ci,len(gamess_output)):
	    ci_line = string.split(gamess_output[i])
	    if len(ci_line) != 9:
		end_ci = i
		break

	# we need to expand out the geminal pairs into separate determinants
	# the start in the form (c1 + c2) * (c3 + c4) * (c5 + c6)...
	# expanding to (c1*c3 + c1*c4 + c2*c3 + c2*c4)*(c5 * c6)*(...)
	# I've implemented a recursion in a loop.
	index = 2
	idet = 0
	for p in range(npair):
	#for p in range(npair-1,-1,-1):
	    old = len(CI)
	    coef1 = string.atof(cicoef[index])
	    index += 1
	    coef2 = string.atof(cicoef[index])
	    index += 2
	    
	    ci_line = string.split(gamess_output[p+start_ci])
	    orb1 = string.atoi(ci_line[1]) - 1
	    orb2 = string.atoi(ci_line[2]) - 1
	    
	    pairCI    = range(old*2)
	    pairAlpha = range(old*2)
	    for i in range(len(pairAlpha)):
		pairAlpha[i] = range(norbitals)
	    
	    for ci in range(len(pairCI)):
		branch1 = ci%old
		branch2 = old-1-ci%old
		#branch2 = branch1
		for i in range(len(AlphaOcc[ci%old])):
		    if ci < old:
			pairAlpha[ci][i] = AlphaOcc[branch1][i]
		    else:
			pairAlpha[ci][i] = AlphaOcc[branch2][i]
			
		if ci < old:
		    pairCI[ci] = CI[branch1] * coef1
		    pairAlpha[ci][orb1] = 1
		else:
		    pairCI[ci] = CI[branch2] * coef2
		    pairAlpha[ci][orb2] = 1

	    #for ci in range(len(pairAlpha)):
	    #	for o in range(2*(p+1)):
	    #	    if o % 2 == 0:
	    #		print pairAlpha[ci][o+ncore],
	    #	print
	    print "Geminal ",p," with coeffs ",coef1, ", ",coef2, " uses orbitals ", orb1, " and ", orb2
	    #print "Determinant CI coefficients = ",pairCI	
	    
	    CI = pairCI
	    AlphaOcc = pairAlpha
	ndeterminants = pow(2,(end_ci - start_ci))	
	#end npair > 0

    BetaOcc = copy.deepcopy(AlphaOcc)

    if nseto > 0 and nseto != 2 and npair > 0:
	print "\n\nWarning: the script maybe doesn't know how to handle npair=",npair, " with nseto=",nseto
    if nseto > 2:
	print "\n\nWarning: the script does not handle nseto=",nseto," correctly!!!"

    #if nseto == 2 and spin_mult == 3:
    if nseto > 0 and spin_mult > 1:
	# This case is easy, since the NSETO orbitals are alpha
	# WF = anti[12aa]
	for det in range(len(AlphaOcc)):
	    for orb in range(len(AlphaOcc[det])):
		if orb >= ncore and orb < ncore+nseto:
		    AlphaOcc[det][orb] = 1
	# If you want the other triplet (even though there doesn't appear to be a difference):
	# WF = anti[12(ab + ba)] = anti[12ab] - anti[21ab]
	# then you need to use spinCouple

    elif nseto == 2:
	#The two electrons are spin coupled into a singlet
	# See Eq 33a from "SCF Equations for GVB" by Bobrowicz and Goddard
	# WF = anti[12(ab - ba)] = anti[(12 + 21)ab] = anti[12ab] + anti[21ab]
	for j in range(core_line_number,len(gamess_output)):
	    if string.find(gamess_output[j],'OPEN SHELL ORBITALS') != -1:
		start_ci = j+1
		break

	ci_line = string.split(gamess_output[start_ci])
	orb1 = string.atoi(ci_line[4])-1
	ci_line = string.split(gamess_output[start_ci+1])
	orb2 = string.atoi(ci_line[4])-1
	(AlphaOcc,BetaOcc,CI) = spinCouple(orb1,orb2,AlphaOcc,BetaOcc,CI,spin_mult)
	
	ndeterminants = len(CI)

    assert(ndeterminants == len(CI))

elif scf_type == "VB2000":
    rumer = []
    for i in range(len(gamess_output)):
        if string.find(gamess_output[i],'GENERAL CONTROLS ($GENCTL)') != -1:
            core_line_number = i+8
            core_line = string.split(gamess_output[core_line_number])
            ncore = string.atoi(core_line[3])
	    print "VB2000 settings: ncore=",ncore
	if string.find(gamess_output[i],'RUMER PATTERN') != -1:
	    for r in range(len(cicoef)):
		line = gamess_output[i+r+1]
		#the first number is just the index
		rumer = rumer + [(string.split(line))[1:]]

    coreO = range(1)
    for i in range(len(coreO)):
        coreO[i] = range(norbitals)

    for i in range(len(coreO)):
        for j in range(norbitals):
	    if j < ncore:
		coreO[i][j] = 1
	    else:
		coreO[i][j] = 0

    AlphaOcc = []
    BetaOcc = []
    CI = []
    for r in range(len(rumer)):
	print "rumer =", rumer[r]
	tempA = copy.deepcopy(coreO)
	tempB = copy.deepcopy(coreO)
	tempC = range(1)
	tempC[0] = cicoef[r]
	pcum += tempC[0]*tempC[0]
	for p in range(len(rumer[r])/2):
	    orb1 = ncore-1 + string.atoi(rumer[r][2*p])
	    orb2 = ncore-1 + string.atoi(rumer[r][2*p+1])	    
	    (tempA,tempB,tempC) = spinCouple(orb1,orb2,tempA,tempB,tempC,1) 
	#for ci in range(len(tempA)):
	#    print tempC[ci], ": ",
	#    for o in range(ncore,norbitals):
	#	print tempA[ci][o],
	#    print
	AlphaOcc = AlphaOcc + tempA
	BetaOcc  = BetaOcc + tempB
	CI = CI + tempC

    ndeterminants = len(CI)
    #sys.exit(0)
elif scf_type == "NONE" and ci_type == "ALDET":
    core_line_number = -1
    nstates = 1
    for i in range(len(gamess_output)):
        if string.find(gamess_output[i],'NUMBER OF CORE ORBITALS') != -1:
            core_line_number = i
            core_line = string.split(gamess_output[core_line_number])
            ncore = string.atoi(core_line[5])
	if string.find(gamess_output[i],'NUMBER OF CI STATES REQUESTED') != -1:
            line = string.split(gamess_output[i])
            nstates = string.atoi(line[6])
	if string.find(gamess_output[i],'PARTIAL TWO ELECTRON INTEGRAL TRANSFORMATION') != -1:
            break

    print ""
    for i in range(core_line_number,len(gamess_output)):
        if string.find(gamess_output[i],'ENERGY=') != -1 and \
	       string.find(gamess_output[i],'CONVERGED') == -1:
	    print gamess_output[i],
	    
    if nstates > 1:
	istate = string.atoi(raw_input("Choose which CI state you want [1 to %i]: "%nstates))

    for i in range(core_line_number,len(gamess_output)):
        if string.find(gamess_output[i],'ENERGY=') != -1 and \
	       string.find(gamess_output[i],'CONVERGED') == -1:
	    line = string.split(gamess_output[i])
	    if string.atoi(line[1]) == istate:
		start_mc_data = i
		break

    end_ci = -1
    start_ci = -1
    for j in range(start_mc_data,len(gamess_output)):
	m=re.search('STATE\s+\d+\s+ENERGY=\s*([\d\-\.]+)',gamess_output[j])
	if m:
	    energy = m.group(1)

	if string.find(gamess_output[j],'ALPH') != -1:
            start_ci = j+2
	if start_ci != -1 and len(gamess_output[j]) == 1:
	    end_ci = j-1
	    break

    for i in range(start_ci,end_ci):
	try:
	    ci_line = string.split(gamess_output[i])
	    coeffs = string.atof(ci_line[4])
	    if abs(coeffs) < detcutoff:
		end_ci = i-1
		break
	except:
	    print "Error extracting ALDET state ", istate, ":"
	    print "First det line = ", start_ci
	    print "Last  det line = ", end_ci
	    print "ENERGY= line   = ", start_mc_data
	    print "Cur det index  = ", i
	    print "Cur det data   = ", gamess_output[i]
	    raise

	if string.find(gamess_output[i+1],'DONE WITH DETERMINANT CI') != -1 or string.find(gamess_output[i+1],'DONE WITH GENERAL CI') != -1:
	    end_ci = i
	    break

    ndeterminants = end_ci - start_ci + 1

    AlphaOcc = range(ndeterminants)
    BetaOcc = range(ndeterminants)
    CI = range(ndeterminants)

    for i in range(ndeterminants):
        AlphaOcc[i] = range(norbitals)
        BetaOcc[i] = range(norbitals)

    for i in range(ndeterminants):
        for j in range(ncore):
            AlphaOcc[i][j] = 1
            BetaOcc[i][j] = 1

    for i in range(ndeterminants):
        ci_line = string.split(gamess_output[i+start_ci])
        CI[i] = string.atof(ci_line[4])
        alpha_occ = ci_line[0]
        beta_occ = ci_line[2]
        for j in range(len(alpha_occ)):
            AlphaOcc[i][j+ncore] = string.atoi(alpha_occ[j])
            BetaOcc[i][j+ncore] = string.atoi(beta_occ[j])
        for k in range(len(alpha_occ)+ncore,norbitals):
            AlphaOcc[i][k] = 0
            BetaOcc[i][k] = 0

######### Use a cutoff criteria to decide which CSFs to include
for i in range(ndeterminants-1,-1,-1):
    try:
	if abs(CI[i]) < detcutoff:
	    #print "Removing",i,"with",CI[i]
	    del CI[i]
	    del AlphaOcc[i]
	    del BetaOcc[i]
    except:
	continue
ndeterminants = len(CI)

######## Remove any orbitals that aren't used

for i in range(norbitals-1,-1,-1):
    keep_this_orbital = 0
    for j in range(ndeterminants):
	if AlphaOcc[j][i] == 1 or BetaOcc[j][i] == 1:
	    keep_this_orbital = 1
	    break
    if(keep_this_orbital == 1):
	continue
    del wavefunction[i]
    for j in range(ndeterminants):
	del AlphaOcc[j][i]
	del BetaOcc[j][i]
norbitals = len(wavefunction)

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
    choice = i

try:
    choice = string.atoi(raw_input("Your choice [%i]:"%choice))
except:
    print "",

myStandardFlags=open(templates[choice],'r')

OUT.write('# Created on %s\n'%(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
OUT.write('# Using gamess output file:  %s\n'% os.path.abspath(Infile))
OUT.write('# Using ckmft template file: %s\n'% templates[choice])
OUT.write('# Orbitals are: %s'%orbital_name[orb_choice])
#let's save some of the important info from a GAMESS calculation
for i in range(len(gamess_output)):
    line = -1;
    if string.find(gamess_output[i],'RUN TITLE') != -1:
	line = i+2
    if string.find(gamess_output[i]," STATE   %i"%istate) != -1 and \
	   string.find(gamess_output[i],'ENERGY=') != -1 and \
	   string.find(gamess_output[i],'SYM=') != -1:
	line = i
    if string.find(gamess_output[i],'CCSD(T) ENERGY:') != -1:
	line = i
    if string.find(gamess_output[i],'CCSD[T] ENERGY:') != -1:
	line = i
    if string.find(gamess_output[i],'CCSD') != -1 and \
	   string.find(gamess_output[i],'ENERGY:') != -1 and \
	   string.find(gamess_output[i],'CORR.') != -1:
	line = i
    if string.find(gamess_output[i],'MBPT(2) ENERGY:') != -1:
	line = i
    if string.find(gamess_output[i],'CORR.') != -1 and \
	   string.find(gamess_output[i],'CR-CC') != -1:
	line = i
    if string.find(gamess_output[i],'FINAL') != -1:
	line = i
    if string.find(gamess_output[i],'$BASIS') != -1:
	line = i
    if string.find(gamess_output[i],'ITER:') != -1:
	line = -1
	    
    if line > 0:
	    OUT.write('#%s' % gamess_output[line])
	    print '#%s' % gamess_output[line],

	    
OUT.write("\n")
OUT.write(myStandardFlags.read());
myStandardFlags.close()

OUT.write('atoms\n %i\n'%atoms)
OUT.write('charge\n %i\n'%charge)
OUT.write('energy\n %s\n'%energy)

if string.atof(energy) >= 0.0:
    print "\nEnergy",energy," didn't converge!!! Quitting.\n"
    sys.exit(0)
    
OUT.write('norbitals\n %i\n'%norbitals)
OUT.write('nbasisfunc\n %i\n'%nbasisfunc)
OUT.write('ndeterminants\n %i\n'%ndeterminants)
OUT.write('&\n')

##################  PRINT FLAGS: END      #######################

##################  PRINT GEOMETRY: BEGIN #######################

OUT.write('&geometry\n')
p = re.compile("[0-9]+")
for line in geometry:
    atom = line[0]
    # if the atom title has a number in it, then we need to remove it
    # so that QMC believes that all the atoms are the same,
    # since Jastrows are specific to the label.
    atom = p.sub("",atom)
    OUT.write('%s\t%i\t%f\t%f\t%f\n'\
	      %(atom,string.atof(line[1]),line[2],line[3],line[4]))
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
        elif bf[0][0] == 'P' : nbf = nbf + 3
        elif bf[0][0] == 'D' : nbf = nbf + 6
        elif bf[0][0] == 'F' : nbf = nbf + 10
	elif bf[0][0] == 'G' : nbf = nbf + 15
	elif bf[0][0] == 'H' : nbf = nbf + 21
	elif bf[0][0] == 'I' : nbf = nbf + 28
        elif bf[0][0] == 'L' : nbf = nbf + 4
	else:
	    print "Error: we don't know about basis function type: ",bf[0][0]
	    sys.exit(0)
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
	# There are a few special basis function types, and you have to
	# program them individually
        if pbf[0][0] == 'L' :
	    mterms = getM('S')
	    for m in mterms:
		OUT.write('\t%i\t%s\n'%(len(pbf),m))
		for gs in pbf:
		    OUT.write('\t\t%s\t%s\n'%(gs[1],normalize(m,gs[2],gs[1])))
	    mterms = getM('P')
	    for m in mterms:
		OUT.write('\t%i\t%s\n'%(len(pbf),m))
		for gs in pbf:
		    OUT.write('\t\t%s\t%s\n'%(gs[1],normalize(m,gs[3],gs[1])))
	else:
	    mterms = getM(pbf[0][0])
	    for m in mterms:
		OUT.write('\t%i\t%s\n'%(len(pbf),m))
		for gs in pbf:
		    OUT.write('\t\t%s\t%s\n'%(gs[1],normalize(m,gs[2],gs[1])))
		    

OUT.write('&\n')

##################  PRINT BASIS: END     ########################

##################  PRINT WAVEFUNCTION: BEGIN ###################

OUT.write('&wavefunction\n\n')
print "charge     = %i"%charge
print "norbitals  = %d\nnbasisfunc = %d\nenergy     = %s\n" % (norbitals,nbasisfunc,energy)


for i in range(norbitals):
    if len(wavefunction[i]) != nbasisfunc:
	print "Error: Orbital",i,"has",len(wavefunction[i]),"basisfunctions, instead of the expected",nbasisfunc
	sys.exit(0)
    for j in range(nbasisfunc):
	try:
	    OUT.write('%20s'%wavefunction[i][j])
	except:
	    print "Error:\nnorbitals = %d\nnbasisfunc = %d\ni = %d\nj = %d\n" % (norbitals,nbasisfunc,i,j)
	    print "wavefunction is %d by %d\n" % (len(wavefunction),len(wavefunction[i]))
	    sys.exit(1)
	if (j+1)%5 == 0:
	    OUT.write('\n')
    OUT.write('\n\n')

print "Alpha Occupation:"
OUT.write("Alpha Occupation\n")
for i in range(ndeterminants):
    nume = 0
    for j in range(norbitals):
	nume += AlphaOcc[i][j]
        OUT.write('%i '%AlphaOcc[i][j])
	if j >= ncore:
	    sys.stdout.write('%i'%AlphaOcc[i][j])
    OUT.write('\n')
    print " (=",nume,")"
OUT.write('\n')
print ""
print "Beta Occupation:"
OUT.write("Beta Occupation\n")
for i in range(ndeterminants):
    nume = 0
    for j in range(norbitals):
	nume += BetaOcc[i][j]
        OUT.write('%i '%BetaOcc[i][j])
	if j >= ncore:
	    sys.stdout.write('%i'%BetaOcc[i][j])
    OUT.write('\n')
    print " (=",nume,")"
OUT.write('\n')
print ""
constraints = []
OUT.write("CI Coeffs\n")
for i in range(ndeterminants):
    match = 0
    constraints.append(-1)
    if ci_type == "ALDET" or 1:
	for j in range(i):
	    ratio = string.atof(CI[i])/string.atof(CI[j])
	    if abs(abs(ratio)-1.0)< 1e-5:
		match = 1
		# There are some couplings that are required to get the correct spin function.
		# We want to include the constraints so that QMC knows which are free to optimize.
		# You'll get ratio = 1 for singlet (ab-ba), and ratio = -1 for triplet (ab+ba)
		#print "Using CI constraint: Det[%i] = %5.3f * Det[%i]"%(i,ratio,j)
		constraints[i] = j
		OUT.write('c %i %5.3f\n'%(j, ratio))
		break
    if match == 0:
	OUT.write('%s\n'%CI[i])
OUT.write('\n')

print str(ndeterminants) + " CI determinant(s) used, with coefficients:"
cum = 0
for i in range(ndeterminants):
    try:
	ci = string.atof(CI[i])
	cum += ci*ci
	print "%3i) %25.7e has percentage %15.8f, cumulative remaining %15.8e" % (i+1,ci,ci*ci*100,1.0-cum),
	if constraints[i] == -1:
	    print ""
	else:
	    rel_diff = ci/string.atof(CI[constraints[i]])-1.0
	    print ", constrained to %2i %25.7e"%(constraints[i]+1,rel_diff)
    except:
	print "%3i) %25s has percentage %15.10e, cumulative remaining %15.10e" % (i+1,CI[i],0,1.0-cum)



OUT.write('&\n')

##################  PRINT WAVEFUNCTION: END   ###################

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

if 0:
    # up down jastrow
    if nalpha > 0 and nbeta > 0:
	OUT.write('ParticleTypes: Electron_Up Electron_Down\n')
	OUT.write('CorrelationFunctionType: Cambridge2\n')
	OUT.write('NumberOfParameterTypes: 2\n')
	OUT.write('NumberOfParametersOfEachType: 1 8\n')
	OUT.write('Parameters: 0.30 0.3\n')
	OUT.write('NumberOfConstantTypes: 2\n')
	OUT.write('NumberOfConstantsOfEachType: 1  1\n')
	OUT.write('Constants: 0.5  3\n')
	OUT.write('\n')
	
	# up up jastrow
	if nalpha > 1:
	    OUT.write('ParticleTypes: Electron_Up Electron_Up\n')
	    OUT.write('CorrelationFunctionType: Cambridge2\n')
	    OUT.write('NumberOfParameterTypes: 2\n')
	    OUT.write('NumberOfParametersOfEachType: 1 8\n')
	    OUT.write('Parameters: 0.30 0.1\n')
	    OUT.write('NumberOfConstantTypes: 2\n')
	    OUT.write('NumberOfConstantsOfEachType: 1  1\n')
	    OUT.write('Constants: 0.25  3\n')
	    OUT.write('\n')
	    
	# down down jastrow
	if nbeta > 1:
	    OUT.write('ParticleTypes: Electron_Down Electron_Down\n')
	    OUT.write('CorrelationFunctionType: Cambridge2\n')
	    OUT.write('NumberOfParameterTypes: 2\n')
	    OUT.write('NumberOfParametersOfEachType: 1 8\n')
	    OUT.write('Parameters: 0.30 0.1\n')
	    OUT.write('NumberOfConstantTypes: 2\n')
	    OUT.write('NumberOfConstantsOfEachType: 1  1\n')
	    OUT.write('Constants: 0.25  3\n')
	    OUT.write('\n')

	# up nuclear jastrow
	if nalpha > 0:
	    for i in range(len(atom_types)):
		OUT.write('ParticleTypes: Electron_Up ' + atom_types[i] + '\n')
		OUT.write('CorrelationFunctionType: Cambridge2\n')
		OUT.write('NumberOfParameterTypes: 2\n')
		OUT.write('NumberOfParametersOfEachType: 1 8\n')
		OUT.write('Parameters: 0.30 -0.3\n')
		OUT.write('NumberOfConstantTypes: 2\n')
		OUT.write('NumberOfConstantsOfEachType: 1  1\n')
		OUT.write('Constants: 0  3\n')
		#        OUT.write('Constants: -' + atom_type_charges[i] + '\n')
		OUT.write('\n')
		
	# down nuclear jastrow
       	if nbeta > 0:
	    for i in range(len(atom_types)):
		OUT.write('ParticleTypes: Electron_Down ' + atom_types[i] + '\n')
		OUT.write('CorrelationFunctionType: Cambridge2\n')
		OUT.write('NumberOfParameterTypes: 2\n')
		OUT.write('NumberOfParametersOfEachType: 1 8\n')
		OUT.write('Parameters: 0.30 -0.3\n')
		OUT.write('NumberOfConstantTypes: 2\n')
		OUT.write('NumberOfConstantsOfEachType: 1  1\n')
		OUT.write('Constants: 0  3\n')
		OUT.write('\n')
else:
    print "\nDont forget to add jastrows!"
    
OUT.write('&Jastrow\n')

################## PRINT JASTROW: END ################################

################## PRINT PSEUDOPOTENTIAL: BEGIN  #####################
if pp_type != "NONE":
    Inpfile = filebase + "inp"	
    INP = open(Inpfile,'r')
    gamess_input = INP.readlines()
    INP.close()

    #Copy the PP right from the input file. QMcBeaver is programmed to use
    #exactly the same format, except it needs to be a GEN PP
    OUT.write('&pseudopotential\n')
    for i in range(len(gamess_input)):
	if string.find(gamess_input[i].upper(),'$ECP') != -1:
	    k = i+1
	    while string.find(gamess_input[k].upper(),'$END') == -1:
		line = string.split(gamess_input[k])
		if len(line) == 4 and line[1].upper() != "GEN" and line[1].upper() != "NONE":
		    print "Pseudoptential for",line[0], "is",line[1],
		    print ": is unknown. It needs to be GEN or NONE.\n";
		OUT.write(gamess_input[k])
		k += 1
    OUT.write('&\n')
################## PRINT PSEUDOPOTENTIAL: END  #######################
    
print "\nFinished writing file ", Outfile

