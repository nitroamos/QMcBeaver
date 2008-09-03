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

import re
import sys
import copy
import math
import string
import time
import os

if len(sys.argv) < 2: 
    print "gamess2qmcbeaver.py <filename>[.log, .inp.out] [detcutoff=0.0]"
    sys.exit(0)

PI = 3.14159265359

#double factorial
def fact2( i ) :
    if i == 1 or i == 0 or i == -1: return 1.0
    else : return i * fact2( i-2 )

#calculate the normalization coefficient for a GTO (see pg 280 in the HLR book)
def normalize( m, pf, ef ) :
    prefactor = string.atof(pf)
    expfactor = string.atof(ef)
    a = m.count('x')
    b = m.count('y')
    c = m.count('z')
    temp = (fact2(2*a-1)*fact2(2*b-1)*fact2(2*c-1))**(-0.5)
    temp *= (2.0 * expfactor / PI)**(0.75)
    temp *= (4.0 * expfactor)**( 0.5*(a+b+c) )
    return temp*prefactor

# This function will make all the basis functions for m
# It is critical to get the order of these correct.
def getM(type):
    type = string.lower(type)
    
    if type == 's':
	return ["s"]
    elif type == 'p':
	return ["px","py","pz"]
    elif type == 'd':
	return ["dxx","dyy","dzz","dxy","dxz","dyz"]
    elif type == 'f':
	return ["fxxx","fyyy","fzzz","fxxy","fxxz","fxyy","fyyz","fxzz","fyzz","fxyz"]
    elif type == 'g':
	return ["gxxxx","gyyyy","gzzzz","gxxxy","gxxxz","gxyyy","gyyyz","gxzzz","gyzzz",
		"gxxyy","gxxzz","gyyzz","gxxyz","gxyyz","gxyzz"]
    else:
	# If you're looking for some of the hybrid types,
	# then don't program them in this function.
	print "Unknown basis function type: ", type
	sys.exit(0)

    return []

def spinCouple(orb1,orb2,A,B, CI, mult):
    #orb1 = string.atoi(orb1)
    #orb2 = string.atoi(orb2)
    if orb1 == orb2:
	for i in range(len(A)):
	    A[i][orb1] = 1
	    B[i][orb1] = 1	
	return (A,B,CI)
    
    old = len(A)
    newCI = range(old*2)
    newA = range(old*2)
    newB = range(old*2)

    for i in range(len(A)):
	ratio = math.sqrt(1.0/2.0)
	ratio = 1.0
	
	newCI[2*i]      = ratio * CI[i]
	newA[2*i]       = copy.deepcopy(A[i])
	newB[2*i]       = copy.deepcopy(B[i])
	newA[2*i][orb1] = 1
	newB[2*i][orb2] = 1
	
	if mult == 3:
	    newCI[2*i+1]      = "c %i -1.0"%(2*i)
	elif mult == 1:
	    newCI[2*i+1]      = "c %i 1.0"%(2*i)
	else:
	    print "Don't know how to couple mult=",mult
	    sys.exit(0)
	newA[2*i+1]       = copy.deepcopy(A[i])
	newB[2*i+1]       = copy.deepcopy(B[i])
	newA[2*i+1][orb2] = 1
	newB[2*i+1][orb1] = 1
    
    #print "Singlet paired orbitals", orb1, "and", orb2
    #print "--> Remember: these CI can't be optimized or you will mess up the spin function!! <--"
    return (newA,newB,newCI)
    
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

detcutoff = 0
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

if run_type == "ENERGY":
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

geom_data = gamess_output[start_geometry:end_geometry]
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
for line in gamess_output:
    if string.find(line,'TOTAL NUMBER OF BASIS FUNCTIONS') !=-1 :
        nbasisfunc = string.atoi(string.split(line)[6])
    if string.find(line,'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS') !=-1 :
        nbasisfunc = string.atoi(string.split(line)[7])
    if string.find(line,'CHARGE OF MOLECULE') !=-1 :
        charge = string.atoi(string.split(line)[4])
#we'll rely on orbital occupations to indicate charge
	charge = 0
    if string.find(line,'TOTAL NUMBER OF ATOMS') !=-1 :
        atoms = string.atoi(string.split(line)[5])
    if string.find(line,'NUMBER OF OCCUPIED ORBITALS (ALPHA)') != -1 :
	nalpha = string.atoi(line.split('=')[1])
    if string.find(line,'NUMBER OF OCCUPIED ORBITALS (BETA )') != -1 :
	nbeta = string.atoi(line.split('=')[1])

#################### EXTRACT BASIS SET: END #######################

#################### EXTRACT WAVEFUNCTION: BEGIN ##################

# Find where the wavefunction is in the .dat file.

energy_line = -1
start_wavefunction = -1

if scf_type == "RHF":
    for i in range(len(gamess_data)):
        if string.find(gamess_data[i],'E(RHF)=') != -1:
            energy_line = i
        elif string.find(gamess_data[i],'E(R-B3LYP)=') != -1:
            energy_line = i
    energy_data = string.split(gamess_data[energy_line])
    for j in range(len(energy_data)):
        if (energy_data[j] == "E(RHF)="):
            energy = energy_data[j+1]
            energy = energy[0:len(energy)-1]
            break
        elif (energy_data[j] == "E(R-B3LYP)="):
            energy = energy_data[j+1]
            energy = energy[0:len(energy)-1]
            break
    for i in range(energy_line,len(gamess_data)):
        if string.find(gamess_data[i],'LOCALIZED') != -1:
            start_wavefunction = i+2
            break
    if start_wavefunction == -1:
        for i in range(energy_line,len(gamess_data)):
	    if string.find(gamess_data[i],'$VEC') != -1:
	        start_wavefunction = i+1
		break
    for j in range(start_wavefunction,len(gamess_data)):
        if string.find(gamess_data[j],'$END') != -1:
            end_wavefunction = j
            break

elif scf_type == "ROHF":
    for i in range(len(gamess_data)):
        if string.find(gamess_data[i],'E(ROHF)=') != -1:
            energy_line = i
    energy_data = string.split(gamess_data[energy_line])
    for j in range(len(energy_data)):
        if (energy_data[j] == "E(ROHF)="):
            energy = energy_data[j+1]
            energy = energy[0:len(energy)-1]
            break
    for i in range(energy_line,len(gamess_data)):
        if string.find(gamess_data[i],'$VEC') != -1:
            start_wavefunction = i+1
            break
    for j in range(start_wavefunction,len(gamess_data)):
        if string.find(gamess_data[j],'$END') != -1:
            end_wavefunction = j
            break

elif scf_type == "UHF":
    for i in range(len(gamess_data)):
        if string.find(gamess_data[i],'E(UHF)=') != -1:
            energy_line = i
    energy_data = string.split(gamess_data[energy_line])
    for j in range(len(energy_data)):
        if (energy_data[j] == "E(UHF)="):
            energy = energy_data[j+1]
            energy = energy[0:len(energy)-1]
            break
    for i in range(energy_line,len(gamess_data)):
        if string.find(gamess_data[i],'$VEC') != -1:
            start_wavefunction = i+1
            break
    for j in range(start_wavefunction,len(gamess_data)):
        if string.find(gamess_data[j],'$END') != -1:
            end_wavefunction = j
            break

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

    for i in range(len(gamess_data)):
        if string.find(gamess_data[i],'E(GVB)=') != -1:
            energy_line = i
    energy_data = string.split(gamess_data[energy_line])
    for j in range(len(energy_data)):
        if (energy_data[j] == "E(GVB)="):
            energy = energy_data[j+1]
            energy = energy[0:len(energy)-1]
            break
    for i in range(energy_line,len(gamess_data)):
        if string.find(gamess_data[i],'$VEC') != -1:
            start_wavefunction = i+1
            break
    for j in range(start_wavefunction,len(gamess_data)):
        if string.find(gamess_data[j],'$END') != -1:
            end_wavefunction = j
            break

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
    Datafile = filebase + "vec"	
    IN2 = open(Datafile,'r')
    gamess_data = IN2.readlines()
    IN2.close()

    for i in range(energy_line,len(gamess_data)):
        if string.find(gamess_data[i],'$VEC') != -1:
            start_wavefunction = i+1
            break
    for j in range(start_wavefunction,len(gamess_data)):
        if string.find(gamess_data[j],'$END') != -1:
            end_wavefunction = j
            break

elif scf_type == "NONE" and ci_type == "ALDET":

    for i in range(len(gamess_data)):
        if string.find(gamess_data[i],'NO-S OF CI') != -1:
            energy_line = i
    energy_data = string.split(gamess_data[energy_line])
    for j in range(len(energy_data)):
        if (energy_data[j] == "E="):
            energy = energy_data[j+1]
            break
    for i in range(energy_line,len(gamess_data)):
        if string.find(gamess_data[i],'$VEC') != -1:
            start_wavefunction = i+1
            break
    for j in range(start_wavefunction,len(gamess_data)):
        if string.find(gamess_data[j],'$END') != -1:
            end_wavefunction = j
            break

else:
    print "SCFTYP", scf_type, "is not supported."
    sys.exit(0) 

# Get the wavefunction parameters from the .dat file.

orbital_coeffs = [] 
for n in range(start_wavefunction, end_wavefunction):
    len_line = len(gamess_data[n])
    number_of_entries = len_line/15
    line_data = range(number_of_entries)
    for i in range(number_of_entries):
        line_data[number_of_entries-i-1] = \
                   gamess_data[n][len_line-15*(i+1)-1:len_line-15*i-1]
    for j in range(len(line_data)):
        line_data[j] = string.atof(line_data[j])
    orbital_coeffs.append(line_data)

wavefunction = []
alpha_orbitals = []
beta_orbitals = []

trialFunctionType = "restricted"

if scf_type != "UHF":

    current_index = 1
    m = start_wavefunction+1
    temp_coeffs = orbital_coeffs[0]
    while 1:
        wavefunction_line = string.split(gamess_data[m])
        if string.atoi(wavefunction_line[0]) == current_index % 100:
            temp_coeffs = temp_coeffs + orbital_coeffs[m-start_wavefunction]
        else:
            wavefunction.append(temp_coeffs)
            current_index = current_index+1
            temp_coeffs = orbital_coeffs[m-start_wavefunction]
        if m == end_wavefunction-1:
            wavefunction.append(temp_coeffs)
            break
        else:
            m = m+1
    norbitals = len(wavefunction)

elif scf_type == "UHF":

    trialFunctionType = "unrestricted"

    current_index = 1

    alpha_only = 0

    m = start_wavefunction+1
    temp_coeffs = orbital_coeffs[0]

    while 1:
        wavefunction_line = string.split(gamess_data[m])
        if wavefunction_line[0] == "$END":
            alpha_only = 1
            break
        elif string.atoi(wavefunction_line[0]) == current_index:
            temp_coeffs = temp_coeffs + orbital_coeffs[m-start_wavefunction]
        elif string.atoi(wavefunction_line[0]) > current_index:
            alpha_orbitals.append(temp_coeffs)
            current_index = current_index+1
            temp_coeffs = orbital_coeffs[m-start_wavefunction]
        if current_index != 1 and string.atoi(wavefunction_line[0]) == 1:
            alpha_orbitals.append(temp_coeffs)
            break
        else:
            m = m+1

    if alpha_only == 0:
        current_index = 1
        m = m+1
        temp_coeffs = orbital_coeffs[m-start_wavefunction-1]
    while 1:
        if alpha_only == 1:
            break
        wavefunction_line = string.split(gamess_data[m])
        if string.atoi(wavefunction_line[0]) == current_index:
            temp_coeffs = temp_coeffs + orbital_coeffs[m-start_wavefunction]
        else:
            beta_orbitals.append(temp_coeffs)
            current_index = current_index+1
            temp_coeffs = orbital_coeffs[m-start_wavefunction]
        if m == end_wavefunction-1:
            beta_orbitals.append(temp_coeffs)
            break
        else:
            m = m+1
    norbitals = len(alpha_orbitals)

# Get the occupation and CI coefficents.

# RHF, ROHF, and UHF have one determinant.
if scf_type == "RHF" or scf_type == "ROHF" or scf_type == "UHF":
    AlphaOcc = [0]
    BetaOcc = [0]
    AlphaOcc[0] = range(norbitals)
    BetaOcc[0] = range(norbitals)
    
    for i in range(nalpha):
        AlphaOcc[0][i] = 1
    for j in range(nalpha,norbitals):
        AlphaOcc[0][j] = 0

    for i in range(nbeta):
        BetaOcc[0][i] = 1
    for j in range(nbeta,norbitals):
        BetaOcc[0][j] = 0

    ncore = norbitals
    ndeterminants = 1
    CI = [1]

elif scf_type == "GVB":
    core_line_number = -1
    start_ci = 0
    end_ci = 0
    for i in range(len(gamess_output)):
        if string.find(gamess_output[i],'ROHF-GVB INPUT PARAMETERS') != -1:
            core_line_number = i+3
            core_line = string.split(gamess_output[core_line_number])
            ncore = string.atoi(core_line[5])
	    norb = string.atoi(core_line[2])
	    pair_line = string.split(gamess_output[core_line_number+1]) 
	    npair = string.atoi(pair_line[2])
	    nseto = string.atoi(pair_line[5]) 
	    print "GVB settings: mult=",spin_mult,"ncore=",ncore,"norb=",norb,"npair=",npair,"nseto=",nseto
            break

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
	BetaOcc = copy.deepcopy(AlphaOcc)
	ndeterminants = pow(2,(end_ci - start_ci))
	
	#end npair > 0

    if nseto > 0 and nseto != 2 and npair > 0:
	print "\n\nWarning: the script maybe doesn't know how to handle npair=",npair, " with nseto=",nseto
    if nseto > 2:
	print "\n\nWarning: the script does not handle nseto=",nseto," correctly!!!"

    if nseto == 2 and spin_mult == 3:
	# This case is easy, since the NSETO orbitals are both alpha
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

    if nstates > 1:
	istate = string.atoi(raw_input("Choose which CI state you want [1 to %i]: "%nstates))

    for i in range(core_line_number,len(gamess_output)):
        if string.find(gamess_output[i],'ENERGY=') != -1:
	    line = string.split(gamess_output[i])
	    if string.atoi(line[1]) == istate:
		start_mc_data = i
		break

    end_ci = -1
    start_ci = -1
    for j in range(start_mc_data,len(gamess_output)):
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

# GAMESS usually creates many more orbitals than are occupied.  We get rid of
# the unoccupied orbitals because we don't need them.

total_orbitals = norbitals
for i in range(total_orbitals):
    index = total_orbitals-i-1
    keep_this_orbital = 0
    for j in range(len(AlphaOcc)):
	if AlphaOcc[j][index] == 1 or BetaOcc[j][index] == 1:
	    keep_this_orbital = 1
	    break
    if (keep_this_orbital == 0):
	norbitals = norbitals-1
    elif (keep_this_orbital == 1):
	break

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
    #choice = 1
    choice = string.atoi(raw_input("Your choice [0]:"))
except:
    choice = 0

print "You chose ckmf template: ", templates[choice]
myStandardFlags=open(templates[choice],'r')

OUT.write('# Created on %s\n'%(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())))
OUT.write('# Using gamess output file:  %s\n'% os.path.abspath(Infile))
OUT.write('# Using ckmft template file: %s\n'% templates[choice])

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

if abs(string.atof(energy)) < 1e-10:
    print "\nEnergy",energy," didn't converge!!!\n"
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
print "norbitals  = %d\nnbasisfunc = %d\n" % (norbitals,nbasisfunc)

if scf_type != "UHF":

    for i in range(norbitals):
        for j in range(nbasisfunc):
	    try:
	        OUT.write('%20s'%wavefunction[i][j])
	    except:
	        print "Error:\nnorbitals = %d\nnbasisfunc = %d\ni = %d\nj = %d\n" % (norbitals,nbasisfunc,i,j)
		print "wavefunction is %d by %d\n" % (len(wavefunction),len(wavefunction[i]))
		sys.exit(1)
            if (j+1)%3 == 0:
                OUT.write('\n')
        OUT.write('\n\n')

elif scf_type == "UHF":

    OUT.write("Alpha Molecular Orbitals\n\n")
    for i in range(norbitals):
        for j in range(nbasisfunc):
            OUT.write('%20s'%alpha_orbitals[i][j])
            if (j+1)%3 == 0:
                OUT.write('\n')
        OUT.write('\n\n')

    if alpha_only == 0:
        OUT.write("Beta Molecular Orbitals\n\n")
        for i in range(norbitals):
            for j in range(nbasisfunc):
                OUT.write('%20s'%beta_orbitals[i][j])
                if (j+1)%3 == 0:
                    OUT.write('\n')
            OUT.write('\n\n')

    elif alpha_only == 1:
        OUT.write("Beta Molecular Orbitals\n\n")
        for i in range(norbitals):
            for j in range(nbasisfunc):
                OUT.write('%20s'%alpha_orbitals[i][j])
                if (j+1)%3 == 0:
                    OUT.write('\n')
            OUT.write('\n\n')        

print "Alpha Occupation:"
OUT.write("Alpha Occupation\n")
for i in range(ndeterminants):
    nume = 0
    for j in range(norbitals):
	nume += AlphaOcc[i][j]
        OUT.write('%i\t'%AlphaOcc[i][j])
	if j >= ncore:
	    print '%i'%AlphaOcc[i][j],
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
        OUT.write('%i\t'%BetaOcc[i][j])
	if j >= ncore:
	    print '%i'%BetaOcc[i][j],	
    OUT.write('\n')
    print " (=",nume,")"
OUT.write('\n')
print ""
OUT.write("CI Coeffs\n")
for i in range(ndeterminants):
    match = 0
    if ci_type == "ALDET":
	for j in range(i):
	    ratio = string.atof(CI[i])/string.atof(CI[j])
	    if abs(ratio) == 1.0:
		match = 1
		# There are some couplings that are required to get the correct spin function.
		# We want to include the constraints so that QMC knows which are free to optimize.
		# You'll get ratio = 1 for singlet (ab-ba), and ratio = -1 for triplet (ab+ba)
		print "Using CI constraint: Det[%i] = %5.3f * Det[%i]"%(i,ratio,j)
		OUT.write('c %i %5.3f\n'%(j, ratio))
    if match == 0:
	OUT.write('%s\n'%CI[i])
OUT.write('\n')

print str(ndeterminants) + " CI determinant(s) used, with coefficients:"
cum = 0
for i in range(ndeterminants):
    try:
	ci = string.atof(CI[i])
	cum += ci*ci
	print "%3i) %25.7f has percentage %15.10e, cumulative remaining %15.10e" % (i+1,ci,ci*ci,1.0-cum)
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

