#!/usr/bin/env python

# This program will shorten a config file to the input number of configs
# or the maximum number of configs in the file, whichever is smaller

import sys
import string

if len(sys.argv) != 3:
    print "USAGE: % shortenconfigs.py <ConfigFile> <MaxNumberOfConfigs> "
    sys.exit(0)

ConfigFileName = sys.argv[1]
OutputFileName = ConfigFileName + ".Short"
MaxNumber = string.atoi(sys.argv[2])

ConfigFile = open(ConfigFileName)
OutputFile = open(OutputFileName,'w')

Counter = -1

while Counter < MaxNumber:
    line = ConfigFile.readline()
#    OutputFile.write(line)

    if string.find(line,'&') != -1:
        Counter = Counter + 1

    if Counter < MaxNumber:
        OutputFile.write(line)
