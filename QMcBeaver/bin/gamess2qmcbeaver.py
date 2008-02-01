#!/usr/bin/env python

# Get the SCF type, run type, and CI type of the calculation.  Usable run types
# are ENERGY and OPTIMIZE.    Usable SCF types are RHF, UHF, ROHF, and NONE.

#This script requires 1 argument, which is the output from a gamess run
#due to this script's method of searching for the rest of the files, the
#extension of the script needs to be .inp.out (7 characters)

#MCSCF notes:
#This script chooses the number of determinants that it will use based on
#the number of determinants chosen with the tolerance PRTTOL in the $DET section
#A QMC wavefunction does not have to be normalized, so the coefficients used
#do not need to be normalized.

# To make a multideterminant trial function, do a MCSCF calculation, then do a
# calculation with SCFTYP=NONE and CITYP=ALDET with the natural orbitals of the
# MCSCF as the read in $VEC.  This will provide the right CI expansion
# coefficients for the MCSCF wavefunction.

#to come up with an MCSCF input file that is based on a molecular minimum:
#1) run gamess with some method like mp2 to generate an initial guess wavefunction
# this is necessary since any CI method will require a $VEC section.
#
#2) using the new $VEC and minimized geometry, make a new gamess input file
# to use the CI parameters of your choice. Run an MCSCF calculation. This calculation
# will not provide the expansion coefficients.
#
#3) make a 3rd gamess input file with a $VEC section from the natural orbitals
# and minimized geometry of the previous run. Specify SCFTYPE=NONE and CITYP=ALDET
# meaning that this will be a single point calculation, using the provided wavefunction.

import sys
import string
import time
import os

if len(sys.argv) < 2: 
    print "gamess2qmcbeaver.py <filename>[.log, .inp.out] [detcutoff]"
    sys.exit(0)

PI = 3.14159265359

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

run_type = "ENERGY"
scf_type = "RHF"
ci_type  = "NONE"

detcutoff = 0
if len(sys.argv) == 3:
	detcutoff = string.atof(sys.argv[2])
	
# Get run type and scf type

for line in gamess_output:
    if string.find(line,'$CONTRL') != -1:
        control_line = string.split(line)
        for element in control_line:
            if string.find(element,'SCFTYP=') != -1:
                scf_type = element[7:len(element)]
            elif string.find(element,'RUNTYP=') != -1:
                run_type = element[7:len(element)]
            elif string.find(element,'CITYP=') != -1:
                ci_type = element[6:len(element)]

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
bfnumber = 0
for line in basisdata:
    if line != '\n':
        line = string.replace(line,')',' ')
        line = string.replace(line,'(',' ')
        line = string.split(line)
        if len(line) == 1 and line[0] != str(bfnumber) and line[0] != str(bfnumber+1):
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
            
        elif len(line) > 1 and line[0] == str(bfnumber+1):
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
            bfnumber = bfnumber+1

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
	charge = 0
    if string.find(line,'TOTAL NUMBER OF ATOMS') !=-1 :
        atoms = string.atoi(string.split(line)[5])
    if string.find(line,'NUMBER OF OCCUPIED ORBITALS (ALPHA)') != -1 :
        nalpha = string.atoi(string.split(line)[6])
    if string.find(line,'NUMBER OF OCCUPIED ORBITALS (BETA )') != -1 :
        nbeta = string.atoi(string.split(line)[7])

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
        orbital_coeffs.append(line_data[j])

wavefunction = []
alpha_orbitals = []
beta_orbitals = []

trialFunctionType = "restricted"

if scf_type != "UHF":

    norbitals = len(orbital_coeffs)/nbasisfunc

    wavefunction = range(norbitals)
    temp_coeffs = []
    
    for i in range(norbitals):
        temp_coeffs = range(nbasisfunc)
        for j in range(nbasisfunc):
            temp_coeffs[j] = orbital_coeffs[i*nbasisfunc+j]
        wavefunction[i] = temp_coeffs

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
if ci_type != "ALDET" and scf_type != "GVB":
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

    ndeterminants = 1
    CI_coeffs = [1]

elif scf_type == "GVB":
    core_line_number = -1
    for i in range(len(gamess_output)):
        if string.find(gamess_output[i],'ROHF-GVB INPUT PARAMETERS') != -1:
            core_line_number = i+3
            core_line = string.split(gamess_output[core_line_number])
            ncore = string.atoi(core_line[5])
            break
		
    for j in range(core_line_number,len(gamess_output)):
        if string.find(gamess_output[j],'CI COEFFICIENTS') != -1:
            start_ci = j+2
            break
    
    for i in range(start_ci,len(gamess_output)):
	ci_line = string.split(gamess_output[i])
	if len(ci_line) != 9:
	    end_ci = i
	    break
    
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
    index = 3
    idet = 0
    for p in range(end_ci-start_ci):
	old = len(CI_coeffs)
	coef1 = string.atof(cicoef[index])
	index += 1
	coef2 = string.atof(cicoef[index])
	index += 3

	ci_line = string.split(gamess_output[p+start_ci])
	orb1 = string.atoi(ci_line[1]) - 1
	orb2 = string.atoi(ci_line[2]) - 1

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

	print "Geminal ",p," with coeffs ",coef1, ", ",coef2, " uses orbitals ", orb1, " and ", orb2
	print "Determinant CI coefficients = ",pairCI	

	CI_coeffs = pairCI
	AlphaOccupation = pairAlpha
    BetaOccupation = AlphaOccupation
    ndeterminants = 2*(end_ci - start_ci)
    assert(ndeterminants == len(CI_coeffs))

elif scf_type == "NONE" and ci_type == "ALDET":
    core_line_number = -1
    for i in range(len(gamess_output)):
        if string.find(gamess_output[i],'NUMBER OF CORE ORBITALS') != -1:
            core_line_number = i
            core_line = string.split(gamess_output[core_line_number])
            ncore = string.atoi(core_line[5])
            break

    for i in range(core_line_number,len(gamess_output)):
        if string.find(gamess_output[i],'ENERGY=') != -1:
            start_mc_data = i
            break

    for j in range(start_mc_data,len(gamess_output)):
        if string.find(gamess_output[j],'ALPH') != -1:
            start_ci = j+2
            break
    
    for i in range(start_ci,len(gamess_output)):
	    try:
		    ci_line = string.split(gamess_output[i])
		    coeffs = string.atof(ci_line[4])
		    if abs(coeffs) < detcutoff:
			    end_ci = i
			    break
	    except:
		    print "i        = ", i
		    print "start_ci = ", start_ci
		    print "line     = ", gamess_output[i]
		    raise
	    if string.find(gamess_output[i+1],'DONE WITH DETERMINANT CI') != -1:
		    end_ci = i
		    break

    ndeterminants = end_ci - start_ci

    AlphaOccupation = range(ndeterminants)
    BetaOccupation = range(ndeterminants)
    CI_coeffs = range(ndeterminants)

    for i in range(ndeterminants):
        AlphaOccupation[i] = range(norbitals)
        BetaOccupation[i] = range(norbitals)

    for i in range(ndeterminants):
        for j in range(ncore):
            AlphaOccupation[i][j] = 1
            BetaOccupation[i][j] = 1

    for i in range(ndeterminants):
        ci_line = string.split(gamess_output[i+start_ci])
        CI_coeffs[i] = ci_line[4]
        alpha_occ = ci_line[0]
        beta_occ = ci_line[2]
        for j in range(len(alpha_occ)):
            AlphaOccupation[i][j+ncore] = string.atoi(alpha_occ[j])
            BetaOccupation[i][j+ncore] = string.atoi(beta_occ[j])
        for k in range(len(alpha_occ)+ncore,norbitals):
            AlphaOccupation[i][k] = 0
            BetaOccupation[i][k] = 0

# GAMESS usually creates many more orbitals than are occupied.  We get rid of
# the unoccupied orbitals because we don't need them.

total_orbitals = norbitals

if (ci_type == "ALDET" or scf_type == "GVB"):
    total_orbitals = norbitals
    for i in range(total_orbitals):
        index = total_orbitals-i-1
        keep_this_orbital = 0
        for j in range(len(AlphaOccupation)):
            if (AlphaOccupation[j][index] == 1):
                keep_this_orbital = 1
                break
            if (BetaOccupation[j][index] == 1):
                keep_this_orbital = 1
                break
        if (keep_this_orbital == 0):
            norbitals = norbitals-1
        elif (keep_this_orbital == 1):
            break
elif (ci_type != "ALDET"):
    if (nalpha >= nbeta):
        norbitals = nalpha
    elif (nbeta > nalpha):
        norbitals = nbeta

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
OUT.write('# Using gamess output file:  %s\n'% os.path.abspath(Infile))
OUT.write('# Using ckmft template file: %s\n'% templates[choice])

#let's save some of the important info from a GAMESS calculation
for i in range(len(gamess_output)):
    line = -1;
    if string.find(gamess_output[i],'RUN TITLE') != -1:
	    line = i+2
    if string.find(gamess_output[i],'REFERENCE ENERGY:') != -1:
	    line = i
    if string.find(gamess_output[i],' STATE') != -1 and string.find(gamess_output[i],'ENERGY') != -1:
	    line = i
    if string.find(gamess_output[i],'CCSD(T) ENERGY:') != -1:
	    line = i
    if string.find(gamess_output[i],'FINAL') != -1:
	    line = i
    if string.find(gamess_output[i],'$BASIS') != -1:
	    line = i
    if line > 0:
	    OUT.write('#%s' % gamess_output[line])
	    print '#%s' % gamess_output[line],

	    
OUT.write("\n")
OUT.write(myStandardFlags.read());
myStandardFlags.close()

OUT.write('atoms\n %i\n'%atoms)
OUT.write('charge\n %i\n'%charge)
OUT.write('energy\n %s\n'%energy)
OUT.write('norbitals\n %i\n'%norbitals)
OUT.write('nbasisfunc\n %i\n'%nbasisfunc)
OUT.write('ndeterminants\n %i\n'%ndeterminants)
OUT.write('&\n')

##################  PRINT FLAGS: END      #######################

##################  PRINT GEOMETRY: BEGIN #######################

OUT.write('&geometry\n')
for line in geometry:
    OUT.write('%s\t%i\t%f\t%f\t%f\n'\
                       %(line[0],string.atof(line[1]),line[2],line[3],line[4]))
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

OUT.write("Alpha Occupation\n")
for i in range(ndeterminants):
    for j in range(norbitals):
        OUT.write('%i\t'%AlphaOccupation[i][j])
	print '%i\t'%AlphaOccupation[i][j],
    OUT.write('\n')
    print '\n',
OUT.write('\n')

OUT.write("Beta Occupation\n")
for i in range(ndeterminants):
    for j in range(norbitals):
        OUT.write('%i\t'%BetaOccupation[i][j])
    OUT.write('\n')
OUT.write('\n')

OUT.write("CI Coeffs\n")
for i in range(ndeterminants):
    OUT.write('%s\n'%CI_coeffs[i])
OUT.write('\n')

print str(ndeterminants) + " CI determinant(s) used, with coefficients:"
cum = 0
for i in range(ndeterminants):
    ci = string.atof(CI_coeffs[i])
    cum += ci*ci
    print "%3i) %15.7f has percentage %15.10e, cumulative remaining %15.10e" % (i+1,ci,ci*ci,1.0-cum)


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

print "\nFinished writing file ", Outfile

