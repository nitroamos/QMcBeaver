#!/usr/bin/env python

import sys
import string

if len(sys.argv) < 3:
    print "Useage: % changeqmcbeaverflags.py <file> <flag=value> ... <flag=value>"
    sys.exit(0)

FileName = sys.argv[1]

FlagsAndValues = []
for item in sys.argv[2:]:
    FlagsAndValues.append(string.split(item,'='))

data = open(FileName).readlines()

for i in range(len(data)):
    for FlagAndValue in FlagsAndValues:
        if string.find(data[i],FlagAndValue[0]) != -1:
            data[i+1] = " " + FlagAndValue[1] + "\n"

File = open(FileName,'w')

for line in data:
    File.write(line)

File.close()

