#!/usr/bin/env python

'''

Check to see if the checkpoint files are properly structured.

'''

import glob
import string

def checkFile(fileName):
    file = open(fileName)

    fileOK = 0


    line = file.readline()
    while line:
        if string.find(line,"</walkers>") != -1:
            fileOK = 1
        line = file.readline()


    file.close()

    return fileOK

def checkAllFiles():
    print "Bad Files:"

    fileNames = glob.glob("*.checkpoint.*")

    for fileName in fileNames:
        fileOK = checkFile( fileName )

        if not fileOK:
            print fileName


if __name__ == "__main__":
    checkAllFiles()
