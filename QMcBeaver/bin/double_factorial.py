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

from string import *

def double_factorial( i ) :
    if i <= 1 : return 1.0
    else : return i * double_factorial( i-2 )
    

def normalize( xexp, yexp, zexp, pf, ef ) :
    prefactor = atof(pf)
    expfactor = atof(ef)

    temp = (double_factorial(2*xexp-1)*double_factorial(2*yexp-1)* \
            double_factorial(2*zexp-1)/ \
            double_factorial(2*(xexp+yexp+zexp)-1))**(-0.5)

    return temp*prefactor
