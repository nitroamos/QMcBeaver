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

import string

PI = 3.14159265359
a0 = 0.529177257507

#double factorial
def fact2( i ) :
    if i == 1 or i == 0 or i == -1: return 1.0
    else : return i * fact2( i-2 )

#calculate the normalization coefficient for a GTO (see pg 280 in the HLR book)
def normalize( m, pf, ef ) :
    a = m.count('x')
    b = m.count('y')
    c = m.count('z')
    prefactor = string.atof(pf)
    expfactor = string.atof(ef)
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
