#!/usr/bin/env python

import os
import os.path
import string

def execute(command):
    print command
    os.system(command)

def copy_qmcbeaver(version):
    command = "cp -r QMcBeaver QMcBeaver-%d"%version
    execute(command)

''' 
Delete the CVS directories.
'''
def delete_cvs_directories(directory):

    if os.path.isdir(directory+"CVS"):
        execute("rm -rf " + directory+"CVS")

    everything = os.listdir(directory)

    for thing in everything:
        if os.path.isdir(thing):
            delete_cvs_directories(directory+thing+"/")

def delete_cvs(version):
    delete_cvs_directories("./QMcBeaver-%d/"%version)


def set_version_number(version):
    configure = open("./QMcBeaver-%d/configure.py"%version).readlines()

    out = open("./QMcBeaver-%d/configure.py"%version,"w")

    for line in configure:
        if string.find(line,"self.VER = ") != -1:
            line = string.replace(line, "str( getVersionNumber() )", \
                                  "str("+str(version)+")")

        out.write(line)

    out.close()

def generate_zip_file(version):
    command = "tar -zcvf QMcBeaver-%d.tar.gz QMcBeaver-%d"%(version,version)
    execute(command)

def generate_docs(version):
    # make doxygen documents
    command = "cd QMcBeaver-%d/docs && doxygen | tee doxy.log"%version
    execute(command)
    # build doxygen API manual
    command = "cd QMcBeaver-%d/docs/latex && latex refman.tex"%version
    execute(command)
    command = "cd QMcBeaver-%d/docs/latex && latex refman.tex"%version
    execute(command)
    command = "cd QMcBeaver-%d/docs/latex && dvipdf refman.dvi"%version
    execute(command)
    command = "mv QMcBeaver-%d/docs/latex/refman.pdf QMcBeaver-API_manual-%d.pdf"%(version,version)
    execute(command)
    # build user guide
    command = "cd QMcBeaver-%d/docs && pdflatex beaver"%version
    execute(command)
    command = "mv QMcBeaver-%d/docs/beaver.pdf QMcBeaver-user_guide-%d.pdf"%(version,version)
    execute(command)
    

if __name__ == '__main__':
    import sys

    sys.path.append("./QMcBeaver/")

    import configure

    versionNumber = configure.getVersionNumber()

    copy_qmcbeaver(versionNumber)
    
    delete_cvs(versionNumber)

    set_version_number(versionNumber)

    generate_zip_file(versionNumber)

    generate_docs(versionNumber)
