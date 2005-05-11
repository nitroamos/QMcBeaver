#!/usr/bin/env python

# Plot a QMcBeaver output file on the fly

try:
    # If the package has been installed correctly, this should work:
    import Gnuplot, Gnuplot.funcutils
except ImportError:
    # It may be that the user is just testing out the package by
    # running 'python demo.py' it the package's directory.  If that is
    # the case, the following should work:
    import __init__
    Gnuplot = __init__
    import funcutils
    Gnuplot.funcutils = funcutils

def getRawData(qmcfile):
    # Get the step, energy, variance from the qmc file

    import string

    qmcdata = open(qmcfile,'r').readlines()

    Step = []
    Energy = []
    Variance = []

    ReadingActualData = 0
    for line in qmcdata:
        # Determine when the header has been read in
        if not ReadingActualData and string.find(line,'----') != -1:
            ReadingActualData = 1

        if ReadingActualData and string.find(string.lower(line),'nan') == -1:
            splitline = string.split(line)

            if len(splitline) > 1:
                Step.append(string.atoi(splitline[0]))
                Energy.append(string.atof(splitline[1]))
                Variance.append(string.atof(splitline[2]))

    return Step, Energy, Variance

def getGnuplotData(qmcfile):
    # Formate the raw qmc data so that gnuplot can use it
    Step, Energy, Variance = getRawData(qmcfile)

    EnergyMeanGnuplot = []
    EnergyUpperBoundGnuplot = []
    EnergyLowerBoundGnuplot = []

    for i in range(len(Step)):
        StDev = Variance[i]**0.5
        
        EnergyMeanGnuplot.append([Step[i],Energy[i]])
        EnergyUpperBoundGnuplot.append([Step[i],Energy[i]+StDev])
        EnergyLowerBoundGnuplot.append([Step[i],Energy[i]-StDev])


    return EnergyMeanGnuplot, EnergyUpperBoundGnuplot, EnergyLowerBoundGnuplot


def plotQMCData(qmcfile,refreshtime):
    # plot the qmc data

    import time

    # create and initialize the gnuplot
    g = Gnuplot.Gnuplot(debug=1)
    g.title(qmcfile)
    g.xlabel('Simulation Time')
    g.ylabel('Energy (a.u.)')
    g('set data style linespoints')

    while 1:
        E_mean, E_upper, E_lower = getGnuplotData(qmcfile)
        if len(E_mean) > 1:
            g.plot(E_mean,E_upper,E_lower)
        time.sleep(refreshtime)
        
    raw_input('Please press return to continue...\n')
    

# when executed get the command line parameters and run plot
if __name__ == '__main__':
    import sys
    import string

    if len(sys.argv) != 3:
        print '% newQMcPlot.py <qmcfile> <refresh time in seconds>'
        sys.exit(0)

    qmcfile = sys.argv[1]
    refreshtime = string.atof(sys.argv[2])
    plotQMCData(qmcfile,refreshtime)



