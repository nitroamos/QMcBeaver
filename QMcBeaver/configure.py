#!/usr/bin/env python

#            QMcBeaver
#
#         Constructed by 
#
#     Michael Todd Feldmann 
#              and 
#   David Randall "Chip" Kent IV
#
# Copyright 2000-2005.  All rights reserved.
#
# drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net


#/**************************************************************************
# This SOFTWARE has been authored or contributed to by an employee or 
# employees of the University of California, operator of the Los Alamos 
# National Laboratory under Contract No. W-7405-ENG-36 with the U.S. 
# Department of Energy.  The U.S. Government has rights to use, reproduce, 
# and distribute this SOFTWARE.  Neither the Government nor the University 
# makes any warranty, express or implied, or assumes any liability or 
# responsibility for the use of this SOFTWARE.  If SOFTWARE is modified 
# to produce derivative works, such modified SOFTWARE should be clearly 
# marked, so as not to confuse it with the version available from LANL.   
#
# Additionally, this program is free software; you can distribute it and/or 
# modify it under the terms of the GNU General Public License. Accordingly, 
# this program is  distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A  PARTICULAR PURPOSE.  See the GNU General Public License 
# for more details. 
#**************************************************************************/

# A class containing all of the information necessary to build the code

import os
import os.path
import string

####################### Define Compilers ####################

class Compiler:
    def __init__(self):
        self.CXX = ''
        self.FLAGS = ''
        self.INCLUDE = ''
        self.DEPENDENCY = ''
        self.OPTIMIZATION = ''
        self.DEBUG = ''
        self.PROFILE = ''
        self.LINK = ''

    def printString(self):
        text = 'CXX = ' + self.CXX + '\n'
        text += 'FLAGS = ' + self.FLAGS + '\n'
        text += 'INCLUDE = ' + self.INCLUDE + '\n'
        text += 'DEPENDENCY = ' + self.DEPENDENCY + '\n'
        text += 'OPTIMIZATION = ' + self.OPTIMIZATION + '\n'
        text += 'DEBUG = ' + self.DEBUG + '\n'
        text += 'PROFILE = ' + self.PROFILE + '\n'
        text += 'LINK = ' + self.LINK + '\n'
        return text

class CompilerAIX(Compiler):
    def __init__(self,optimize,debug,profile):
        self.CXX = 'mpCC'
        self.FLAGS = ''
        self.INCLUDE = ''
        self.DEPENDENCY = '-M -qsyntaxonly'

        if optimize:
            self.OPTIMIZATION = '-O3 -qnostrict'
        else:
            self.OPTIMIZATION = ''

        if debug:
            self.DEBUG = '-g'
        else:
            self.DEBUG = ''

        if profile:
            self.PROFILE = '-p'
        else:
            self.PROFILE = ''
            
        self.LINK = '-lm'


class CompilerCCMALLOC(Compiler):
    def __init__(self,optimize,debug,profile):
        self.CXX = 'ccmalloc g++'
        self.FLAGS = '-Wall'
        self.INCLUDE = ''
        self.DEPENDENCY = '-MM'

        if optimize:
            self.OPTIMIZATION = '-O3 -ffast-math'
        else:
            self.OPTIMIZATION = ''

        if debug:
            self.DEBUG = '-g'
        else:
            self.DEBUG = ''

        if profile:
            self.PROFILE = '-pg'
        else:
            self.PROFILE = ''
            
        self.LINK = '-lm -lccmalloc'


class CompilerCrayX1(Compiler):
    def __init__(self,optimize,debug,profile):
        self.CXX = 'CC -h new_for_init -h fp0 -DCRAY_X1 '
        self.FLAGS = '-h new_for_init -h fp0 -DCRAY_X1'
        self.INCLUDE = ''
        self.DEPENDENCY = '-M'

        if optimize:
            self.OPTIMIZATION = '-02'
        else:
            self.OPTIMIZATION = ''

        if debug:
            self.DEBUG = '-g'
        else:
            self.DEBUG = ''

        if profile:
            self.PROFILE = '-p'
        else:
            self.PROFILE = ''
            
        self.LINK = '-lm'


class CompilerGCC(Compiler):
    def __init__(self,optimize,debug,profile):
        self.CXX = 'g++'
        self.FLAGS = '-Wall -Wno-unused -Wno-unused-function -Wno-format'
        self.INCLUDE = ''
        self.DEPENDENCY = '-MM'
        self.LINK = '-lm'
		
        if optimize:
            self.OPTIMIZATION = '-O3 -fexpensive-optimizations -ffast-math #-mtune=???'
        else:
            self.OPTIMIZATION = ''

        if debug:
            self.DEBUG = '-gdwarf-2'
        else:
            self.DEBUG = ''

        if profile:
            self.PROFILE = '-pg'
        else:
            self.PROFILE = ''

class CompilerInsure(Compiler):
    def __init__(self,optimize,debug,profile):
        self.CXX = 'insure'
        self.FLAGS = '-Wall'
        self.INCLUDE = ''
        self.DEPENDENCY = '-MM'

        if optimize:
            self.OPTIMIZATION = '-O3'
        else:
            self.OPTIMIZATION = ''

        if debug:
            self.DEBUG = '-g'
        else:
            self.DEBUG = ''

        if profile:
            self.PROFILE = '-p'
        else:
            self.PROFILE = ''
            
        self.LINK = '-lm'

class CompilerIntel(Compiler):
    def __init__(self,optimize,debug,profile):
        self.CXX = 'icpc'
        self.FLAGS = '-Wall -w1 -Wcheck -ansi'
        self.INCLUDE = ''
        self.DEPENDENCY = '-MM'

        if optimize:
            #The -fast option automatically selects -xP, which means SSE3
            #Using these compile paramters against GCC default parameters,
            #icc outperforms GCC by 7 percent on a harmonic oscillator test.
            self.OPTIMIZATION = '-O3 -ipo -static -parallel -openmp -xN'
        else:
            self.OPTIMIZATION = ''

        if debug:
            self.DEBUG = '-g -debug extended'
        else:
            self.DEBUG = ''

        if profile:
            self.PROFILE = '-p'
        else:
            self.PROFILE = ''
            
        self.LINK = '-lm -lstdc++ -static -lcxaguard'


class CompilerMPICC(Compiler):
    def __init__(self,optimize,debug,profile):
        self.CXX = 'mpiCC'
        self.FLAGS = '-Wall '
        self.INCLUDE = ''
        self.DEPENDENCY = '-MM'

        if optimize:
            self.OPTIMIZATION = '-O3 -ffast-math'
        else:
            self.OPTIMIZATION = ''

        if debug:
            self.DEBUG = '-g'
        else:
            self.DEBUG = ''

        if profile:
            self.PROFILE = '-pg'
        else:
            self.PROFILE = ''
            
        self.LINK = '-lm -lstdc++'
        
class CompilerPGI(Compiler):
    def __init__(self,optimize,debug,profile):
        self.CXX = 'pgCC'
        self.FLAGS = '-Minform=warn'
        self.INCLUDE = ''
        self.DEPENDENCY = '-M'

        if optimize:
            self.OPTIMIZATION = '-fastsse -Minline=levels:10 -Mipa=fast'
        else:
            self.OPTIMIZATION = ''

        if debug:
            self.DEBUG = '-g -Mdwarf2'
        else:
            self.DEBUG = ''

        if profile:
            self.PROFILE = '-pg'
        else:
            self.PROFILE = ''
            
        self.LINK = '-lm'


class CompilerSgiIrix(Compiler):
    def __init__(self,optimize,debug,profile):
        self.CXX = 'CC'
        self.FLAGS = '-64 -LANG:std -ansi'
        self.INCLUDE = ''
        self.DEPENDENCY = '-M'

        if optimize:
            self.OPTIMIZATION = '-O3'
        else:
            self.OPTIMIZATION = ''

        if debug:
            self.DEBUG = '-g'
        else:
            self.DEBUG = ''

        if profile:
            self.PROFILE = '-p'
        else:
            self.PROFILE = ''
            
        self.LINK = '-lm'


class CompilerTru64(Compiler):
    def __init__(self,optimize,debug,profile):
        self.CXX = 'cxx'
        self.FLAGS = '-ieee -underflow_to_zero -D__USE_STD_IOSTREAM -D__FUNCTION__ -DUSING_QSC'
        self.INCLUDE = ''
        self.DEPENDENCY = '-M'

        if optimize:
            self.OPTIMIZATION = '-fast'
        else:
            self.OPTIMIZATION = ''

        if debug:
            if optimize:
                self.DEBUG = '-g3'
            else:
                self.DEBUG = '-g'
        else:
            self.DEBUG = ''

        if profile:
            self.PROFILE = '-p'
        else:
            self.PROFILE = ''
            
        self.LINK = '-lm'


compilers = {
    'AIX'     : CompilerAIX,
    'CrayX1'  : CompilerCrayX1,
    'CCMALLOC' : CompilerCCMALLOC,
    'GCC'     : CompilerGCC,
    'Insure'  : CompilerInsure,
    'Intel'   : CompilerIntel,
    'MPICC'   : CompilerMPICC,
    'PGI'     : CompilerPGI,
    'SgiIrix' : CompilerSgiIrix,
    'Tru64'   : CompilerTru64,
    }


####################### Define MPIs ####################

class Mpi:
    def __init__(self):
        self.PARALLEL = '-DPARALLEL'
        self.FLAGS = ''
        self.INCLUDE = ''
        self.LINK = ''

class MpiNone(Mpi):
    def __init__(self):
        self.PARALLEL = ''
        self.FLAGS = ''
        self.INCLUDE = ''
        self.LINK = ''

class MpiLAMPI(Mpi):
    def __init__(self):
        self.PARALLEL = '-DPARALLEL'
        self.FLAGS = ''
        self.INCLUDE = ' -I$(MPI_ROOT)/include'
        self.LINK = ' -L$(MPI_ROOT)/lib -lmpi -lstdc++'

class MpiOpenMPI(Mpi):
    def __init__(self):
        self.PARALLEL = '-DPARALLEL'
        self.FLAGS = ''
        self.INCLUDE = ' -I$(MPI_ROOT)/include -I$(MPI_ROOT)/include/openmpi/ompi'
        self.LINK = ' -L$(MPI_ROOT)/lib -lmpi'

class MpiSgiIrix(Mpi):
    def __init__(self):
        self.PARALLEL = '-DPARALLEL'
        self.FLAGS = ' -DMPI_NO_CPPBIND'
        self.INCLUDE = ''
        self.LINK = ' -lmpi'

class MpiTru64(Mpi):
    def __init__(self):
        self.PARALLEL = '-DPARALLEL'
        self.FLAGS = ''
        self.INCLUDE = ' $(MPI_COMPILE_FLAGS)'
        self.LINK = ' $(MPI_LD_FLAGS) -lmpi'

mpis = {
    'None'    : MpiNone,
    'LAMPI'   : MpiLAMPI,
    'MPICH'   : Mpi,
    'OpenMPI' : MpiOpenMPI,
    'SgiIrix' : MpiSgiIrix,
    'Tru64'   : MpiTru64,
    }


####################################################################
# Don't Change Below This Line Unless You Know What You Are Doing! #
####################################################################


####################### Define Optional Args ####################

optionalarguments = {
    '--optimize': 'self.optimize = 1',
    '--nooptimize': 'self.optimize = 0',
    '--debug' : 'self.debug = 1',
    '--nodebug' : 'self.debug = 0',
    '--profile' : 'self.profile = 1',
    '--noprofile' : 'self.profile = 0',
    '--hdf5': 'self.hdf5 = 1',
    '--atlas': 'self.atlas = 1',
    '--goto': 'self.goto = 1',
    '--lapack': 'self.lapack = 1',
    '--sprng': 'self.sprng = 1',
    '--float': 'self.float = 1',
    '--gpu': 'self.gpu = 1',
    '--vtune': 'self.vtune = 1',
    }    

#########################


class CommandLineArgs:
    def __init__(self):
        self.optimize = 1
        self.debug = 0
        self.profile = 0
        self.parallel = 0
        self.hdf5 = 0
        self.atlas = 0
	self.goto = 0
        self.lapack = 0
        self.sprng = 0
        self.float = 0
        self.gpu = 0
        self.vtune = 0
        self.TAG = ''
        self.LINK = ''
        self.INCLUDE = ''
        self.EXE = 'QMcBeaver.$(LABEL).$(VERSION).x'
        self.MAKE = self._getMake()

        if len(sys.argv) > 1 and sys.argv[1][0:2] != '--':
            self.COMPILER = sys.argv[1]
        else:
            self.printHelp()
            sys.exit(0)

        if len(sys.argv) > 2 and sys.argv[2][0:2] != "--":
            self.MPI = sys.argv[2]
        else:
            self.MPI = 'None'
            sys.argv.insert(2,'None')

        if self.COMPILER == 'MPICC':
            self.MPI = 'MPICH'

        if self.MPI != 'None':
            self.parallel = 1
        
        for arg in sys.argv[3:]:
            temp = string.split(arg,"=")
            if arg[0:2] == "--":
                self._testArgValidity(temp[0])
                exec optionalarguments[temp[0]]
            else:
                if self.TAG == '':
                    self.TAG = '_'
                self.TAG += temp[0];

            if len(temp) == 2:
                libpath = os.path.expanduser(temp[1])
                includepath = libpath + "/include/"
                libpath += "/lib/"
                if os.path.exists(libpath):
                    self.LINK += " -L" + libpath
                if os.path.exists(includepath):
                    self.INCLUDE += " -I" + includepath
                
        self._testInputValidity()

        self.LABEL = self._getLabel()
    
    def printHelp(self):
        print 'Syntax: ./configure.py compiler mpi [option1 [option2 etc]]'
        print 'Compiler Options'
        compilerKeys = compilers.keys()
        compilerKeys.sort()
        for compiler in compilerKeys:
            print '\t', compiler
        print

        print 'MPI Options'
        mpiKeys = mpis.keys()
        mpiKeys.sort()
        for mpi in mpiKeys:
            print '\t', mpi
        print
        print 'Optional Arguments'
        optArgKeys = optionalarguments.keys()
        optArgKeys.sort()
        for arg in optArgKeys:
            print '\t', arg
        print 'Any Optional Arguments without the -- will be treated as a tag.'
        print 'The vtune option is to create a position-independent exe'
        print 'necessary for vtune\'s Call Graph.\n'
        print 'Some options require special libraries.'
        print 'See lib/README for setup details.\n'
      
    def _getMake(self):
        paths = os.environ['PATH'].split(":")

        for path in paths:
            if os.path.isfile(path+'/gmake'):
                return 'gmake'

        return 'make'

    def _testArgValidity(self, arg):
        if not optionalarguments.has_key( arg ):
            print 'Incorrect argument specified! (%s)'%arg
            self.printHelp()
            sys.exit(1)

    def _testInputValidity(self):
        if not compilers.has_key( self.COMPILER ):
            print 'Incorrect compiler specified! (%s)'%self.COMPILER
            self.printHelp()
            sys.exit(1)

        if not mpis.has_key( self.MPI ):
            print 'Incorrect MPI specified! (%s)'%self.MPI
            self.printHelp()
            sys.exit(1)

    def _getLabel(self):

        extra = [self.COMPILER]

        #the order here should facilitate selecting the exe via tabs
#        if self.optimize:
#            extra.append("optimize")
        if self.hdf5:
            extra.append("hdf5")
        if self.atlas or self.goto:
            extra.append("blas")
        if self.parallel:
            extra.append("parallel")
        if self.profile:
            extra.append("profile")
        if self.float:
            extra.append("float")
        if self.debug:
            extra.append("debug")
        if self.gpu:
           extra.append("gpu")
        if self.vtune:
           extra.append("vtune")
        if self.lapack:
           extra.append("lapack")
        if self.sprng:
           extra.append("sprng")

        return '_'.join(extra)


class MakeConfigBuilder:
    def __init__(self):
        self._inputs = CommandLineArgs()

        self._compiler = compilers[ self._inputs.COMPILER ](
            self._inputs.optimize,
            self._inputs.debug,
            self._inputs.profile)

        self._mpi = mpis[ self._inputs.MPI ]()
        
        self.VERSION = str( getVersionNumber() )
        self.HOME = os.getcwd()        
        self.INCLUDE = '-I$(HOME)/include'

    def printString(self):
        text =  '# Makefile.config is generated by configure.py.  Either use \n'
        text += '# one of the configurations configure.py generates, \n'
        text += '# edit the script to generate your configuration, \n'
        text += '# or edit Makefile.config for your configuration \n'
	text += '# Built with: '
	for arg in sys.argv:
	    text += arg + ' '
	text += '\n\n'
        text += 'TAG = ' + self._inputs.TAG + '\n'
        text += 'COMPILER = ' + self._compiler.CXX + '\n'
        text += 'LABEL = ' + self._inputs.LABEL + self._inputs.TAG + '\n'
        text += 'VERSION = ' + self.VERSION + '\n'
        text += 'HOME = ' + self.HOME + '\n'
        if self._inputs.hdf5:
            text += 'H5 = /ul/amosa/lib/hdf5/lib/\n'
        text += 'INCLUDE = ' + self.INCLUDE + self._inputs.INCLUDE + " " + self._mpi.INCLUDE
        if self._inputs.hdf5:
            text += ' -I/ul/amosa/lib/hdf5/include -I/ul/amosa/lib/szip/szip2.0-linux-enc/include'
        text += '\n'
        text += 'FLAGS = ' + self._compiler.FLAGS + self._mpi.FLAGS
        if self._inputs.hdf5:
            text += ' -DUSEHDF5 '
        if self._inputs.atlas or self._inputs.goto:
            text += ' -DUSEBLAS '
        if self._inputs.atlas and self._inputs.goto:
	    print "You can't specify both ATLAS and Goto!"
	    sys.exit(1)
        if self._inputs.lapack:
            text += ' -DUSELAPACK '
        if self._inputs.sprng:
            text += ' -DUSESPRNG '
        if self._inputs.float:
            text += ' -DSINGLEPRECISION'
        if self._inputs.gpu:
            text += ' -DQMC_GPU'
	if self._inputs.debug:
	    text += ' -DQMC_DEBUG'
        text += '\n'
        text += 'DIROBJ = obj_$(LABEL)\n'
        text += 'MAKE = ' + self._inputs.MAKE + '\n'
        text += 'EXE = ' + self._inputs.EXE + '\n'
        text += 'DEPENDENCY = ' + self._compiler.DEPENDENCY + '\n'
        text += 'OPTIMIZATION = ' + self._compiler.OPTIMIZATION + '\n'
        if self._inputs.vtune:
            self._compiler.DEBUG = '-gstabs -g3'
        text += 'DEBUG = ' + self._compiler.DEBUG + '\n'
        text += 'PROFILE = ' + self._compiler.PROFILE + '\n'
        text += 'PARALLEL = ' + self._mpi.PARALLEL + '\n'
        text += 'LINK = ' + self._compiler.LINK + self._inputs.LINK + self._mpi.LINK + ' -L$(HOME)/lib'
        if self._inputs.hdf5:
            text += ' -L/ul/amosa/lib/szip/szip2.0-linux-enc/lib -static -lsz $(H5)/libhdf5_cpp.a $(H5)/libhdf5.a -lz'
        if self._inputs.sprng:
            text += ' -lsprng'
        if self._inputs.lapack:
            text += ' -llapack -lg2c'
        if self._inputs.atlas:
            text += ' -lf77blas -latlas'
	if self._inputs.goto:
            text += ' -lgoto'

        if self._inputs.gpu:
            text+= ' -lkernel32 -lgdi32 -luser32 -lopengl32 -lglut32 -lversion -lglew32 -lglu32 -lcg -lcgGL'
        if self._inputs.vtune:
            text += ' -Wl,-pie -pie'
        text += '\n'
        text += 'CXX = $(COMPILER) -DVERSION=$(VERSION) $(FLAGS) $(DEBUG) $(PROFILE) $(PARALLEL)\n'
        return text



'''
Gets the most recent date a file was modified in this directory.  The date
is returned in the QMcBeaver versioning format (yyyymmdd).
'''

def __getMostRecentCVSDateLocalDirectory(directory):
    if not os.path.isdir(directory+"CVS"):        
        return 0

    mostRecentVersionNumber = 0

    # get the most recent file date from the CVS files
    data = open( directory + "CVS/Entries" ).readlines()

    for datum in data:
	temp = string.split(datum,"/")
	if len(temp) > 3 and string.find(datum,'dummy') == -1:

	    date = temp[-3]
	    date = string.split(date)
                
	    if len(date) > 1:
		try:
		    day   = string.atoi(date[2])
		    month = date[1]
		    year  = string.atoi(date[4])
		except :
		    print "Error: directory = ",directory, " date =",date, " datum = ",datum
		    continue

		if month == 'Jan': month = 1
		elif month == 'Feb': month = 2
		elif month == 'Mar': month = 3
		elif month == 'Apr': month = 4
		elif month == 'May': month = 5
		elif month == 'Jun': month = 6
		elif month == 'Jul': month = 7
		elif month == 'Aug': month = 8
		elif month == 'Sep': month = 9
		elif month == 'Oct': month = 10
		elif month == 'Nov': month = 11
		elif month == 'Dec': month = 12
		else: print "ERROR: Unknown Month (", month,")!"
	    
	    versionNumber = year * 10000 + month * 100 + day
	else:
	    versionNumber = 0

	if versionNumber > mostRecentVersionNumber:
	    mostRecentVersionNumber = versionNumber


    return mostRecentVersionNumber



'''
Gets the most recent date a file was modified in all subdirectories.
The date is returned in the QMcBeaver versioning format (yyyymmdd).
'''

def __getMostRecentCVSDate(directory):

    # get the version information for the local directory
    mostRecentVersionNumber = \
                            __getMostRecentCVSDateLocalDirectory(directory)

    directoryContents = ['src','include']
    for item in directoryContents:
	versionNumber = __getMostRecentCVSDateLocalDirectory( directory + item + "/" )

        if versionNumber > mostRecentVersionNumber:
            mostRecentVersionNumber = versionNumber

    return mostRecentVersionNumber


'''
Gets the correct version number of the software.
'''

def getVersionNumber():
    return __getMostRecentCVSDate("./")



if __name__ == '__main__':
    import sys

    dat = MakeConfigBuilder()
    print dat.printString()

    file = open('Makefile.config','w')
    file.write(dat.printString())
    file.close()




