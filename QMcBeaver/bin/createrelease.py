#!/usr/bin/env python

import os
import os.path
import string

''' 
Delete the CVS directories.
'''
def delete_cvs_directories(directory):

    if os.path.isdir(directory+"CVS"):
        print "delete "+directory+"CVS"
        os.system("rm -rf " + directory+"CVS")

    everything = os.listdir(directory)

    for thing in everything:
        if os.path.isdir(thing):
            delete_cvs_directories(directory+thing+"/")

def set_version_number(version):
    configure = open("./QMcBeaver/configure.py").readlines()

    out = open("./QMcBeaver/configure.py","w")

    for line in configure:
        if string.find(line,"self.VER = ") != -1:
            line = string.replace(line, "str( getVersionNumber() )", \
                                  "str("+str(version)+")")

        out.write(line)

    out.close()

def generate_zip_file(version):
    command = "tar -zcvf QMcBeaver-%d.tar.gz"%(version)
    os.system(command)

if __name__ == '__main__':

    import configure
    
    delete_cvs_directories("./QMcBeaver/")

    versionNumber = configure.getVersionNumber()

    set_version_number(versionNumber)

    generate_zip_file(version)
