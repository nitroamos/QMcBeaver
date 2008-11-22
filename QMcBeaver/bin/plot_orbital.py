#! usr/bin/env python

from Gnuplot import Gnuplot,Data
from math import exp
import sys,string
from Numeric import arange

# Usage: ./plot_orbital.py outfile,
# Where outfile is a QMcBeaver output file created with the
# plot_replacement_orbitals flag set to 1.  This will plot the original and
# replacement orbitals for visualization.

g = Gnuplot()
g.reset()

number = 0
atom = -1
orbital = -1

file = sys.argv[1]
IN = open(file,'r')
data = IN.readlines()
IN.close()

print "Read file", file

rc = 0.0
n_alpha = 5
alpha = range(5)
for element in alpha:
    element = 0.0
C = 0.0
sgn = 1

end = -1
n_coeffs = -1

while 1:
    if string.find(data[number],'***ElectronNucleusCuspParameters***') != -1:
        line = string.split(data[number-1])
        atom = line[5]
        orbital = line[7]
        line = string.split(data[number+1])
        rc = string.atof(line[2])
        line = string.split(data[number+3])
        for i in range(5):
            alpha[i] = string.atof(line[i])
        line = string.split(data[number+5])
        n_ideal = len(line)
        ideal = range(n_ideal)
        for i in range(n_ideal):
            ideal[i] = string.atof(line[i])
        line = string.split(data[number+7])
        C = string.atof(line[2])
        line = string.split(data[number+8])
        sgn = string.atoi(line[2])
        end = number+14
        n_coeffs = 0
        while string.find(data[end],'local') == -1:
            n_coeffs = n_coeffs + 1
            end = end+1
        exponents = range(n_coeffs)
        coeffs = range(n_coeffs)
        for i in range(n_coeffs):
            line = string.split(data[number+14+i])
            exponents[i] = string.atof(line[0])
            coeffs[i] = string.atof(line[1])
        print "Rep orbital for atom", atom, "orbital", orbital, "rc", rc

        max_dist = rc*1.1
        dr = max_dist/100
        r = range(101)
        rep = range(101)
        orig = range(101)
        for i in range(101):
            r[i] = dr*i
            rep[i] = 0.0
            orig[i] = 0.0
        for i in range(101):
            R = alpha[0]+r[i]*(alpha[1]+r[i]*(alpha[2]+r[i]*(alpha[3]+r[i]*alpha[4])))
            rep[i] = C+sgn*exp(R)
            orig[i] = 0.0
            for j in range(n_coeffs):
                orig[i] = orig[i]+coeffs[j]*exp(-exponents[j]*r[i]*r[i])
        d1 = Data(r[:],rep[:],with='lines')
        d2 = Data(r[:],orig[:],with='lines')
        g.plot(d1,d2)
        raw_input("Press enter to continue...")

    number = number+1
            
                            
                       
