#!/usr/bin/env python

# A class containing all of the information necessary to build the code

import os
import os.path
import string

class ControlMake:
    # self.MAKE -- make program
    # self.EXE  -- executable name
    # self.CXX  -- compiler
    # self.DEP  -- flag used to get dependancy
    # self.OPT  -- optimization flags
    # self.DBG  -- debug flags
    # self.PFL  -- profiling flags
    # self.PLL  -- parallel flags
    # self.LNK  -- linking flags
    # self.VER  -- version number of the software
    
    def printHelp(self):
        print 'System Options'
        print '\tdefault'
        print '\tlinux'
        print '\tinsure'
        print '\tsgi'
        print '\tllnl'
        print '\tlanl'
	print '\tcplant'
        print '\tccmalloc'
        print
        print 'Compilation Options'
        print '\t--{no}optimize'
        print '\t--{no}debug'
        print '\t--{no}profile'
        print '\t--{no}parallel'
        print '\t--exe=EXECUTABLE_NAME'
        print '\t--cxx=COMPILER'
        print '\t--make=MAKE_PROGRAM'
        print '\t--optimize-flags=FLAGS'
        print '\t--debug-flags=FLAGS'
        print '\t--profile-flags=FLAGS'
        print '\t--link-flags=FLAGS'
    
    def __init__(self,inputflags):
        if len(inputflags) > 1 and inputflags[1][0:2] != '--':
            self.SYS = inputflags[1]
        else:
            self.printHelp()
            sys.exit(0)

        self.testSysValidity()
        
        import os
        self.HOME = os.getcwd()
        
        self.setINC()
        self.setBase()
        self.parseInput(inputflags)
        self.setEXE()
        self.VER = str( self.getVersionNumber() )

    def testSysValidity(self):
        if self.SYS != 'default' and self.SYS != 'linux' and \
           self.SYS != 'insure' and self.SYS != 'sgi' and \
           self.SYS != 'llnl' and self.SYS != 'lanl' and \
	   self.SYS != 'cplant' and \
           self.SYS != 'ccmalloc':
            print 'Incorrect system specified!'
            self.printHelp()
            sys.exit(0)

    def setINC(self):
        self.INC = '-I'+self.HOME+'/include'

    def setEXE(self):
        if self.EXE == '':
            self.EXE = 'QMcBeaver.' + self.SYS
            if self.PLL != '':
                self.EXE = 'p' + self.EXE

    def setBase(self):
        if   self.SYS == 'linux':
            self.MAKE = 'gmake'
            self.EXE  = ''
            self.CXX  = 'g++ -static -Wall'
            self.DEP  = '-MM'
            self.setOptimize()
            self.DBG  = ''
            self.PFL  = ''
            self.PLL  = ''
            self.LNK  = '-lm'
            
        elif   self.SYS == 'insure':
            self.MAKE = 'gmake'
            self.EXE  = ''
            self.CXX  = 'insure'
            self.DEP  = '-MM'
            self.OPT  = ''
            self.DBG  = '-g'
            self.PFL  = ''
            self.PLL  = ''
            self.LNK  = '-lm'
            
        elif self.SYS == 'sgi':
            self.MAKE = 'gmake'
            self.EXE  = ''
            self.CXX  = 'CC -LANG:std -ansi'
            self.DEP  = '-M'
            self.setOptimize()
            self.DBG  = ''
            self.PFL  = ''
            self.PLL  = ''
            self.LNK  = '-lm'
            
        elif self.SYS == 'llnl':
            self.MAKE = 'gmake'
            self.EXE  = ''
            self.CXX  = 'newmpCC'
            self.DEP  = '-M'
            self.setOptimize()
            self.DBG  = ''
            self.PFL  = ''
            self.setParallel()
            self.LNK  = '-lm'

        elif self.SYS == 'lanl':
            self.MAKE = 'gmake'
            self.EXE  = ''
            self.CXX  = 'CC -64 -LANG:std -ansi'
            self.DEP  = '-M'
            self.setOptimize()
            self.DBG  = ''
            self.PFL  = ''
            self.PLL  = ''
            self.LNK  = '-lm'
            self.setParallel()

        elif self.SYS == 'cplant':
            self.MAKE = 'gmake'
            self.EXE  = ''
            self.CXX  = '/usr/local/cplant/west/current/bin/c++'
            self.DEP  = '-MM'
            self.setOptimize()
            self.DBG  = ''
            self.PFL  = ''
            self.PLL  = ''
            self.LNK  = '-lm'
            self.setParallel()

        elif self.SYS == 'ccmalloc':
            self.MAKE = 'gmake'
            self.EXE  = ''
            self.CXX  = '/ul/chip/TEMP/MemDebugging/ccmalloc/ccmalloc-0-3-4/' \
                        + 'bin/ccmalloc g++'
            self.DEP  = '-MM'
            self.OPT  = ''
            self.setDebug()
            self.PFL  = ''
            self.PLL  = ''
            self.LNK  = '-lm -L/ul/chip/TEMP/MemDebugging/ccmalloc/' \
                        + 'ccmalloc-0-3-4/lib -lccmalloc'
            
        elif self.SYS == 'default':
            self.MAKE = 'gmake'
            self.EXE  = ''
            self.CXX  = 'g++'
            self.DEP  = '-MM'
            self.OPT  = ''
            self.DBG  = ''
            self.PFL  = ''
            self.PLL  = ''
            self.LNK  = '-lm'
            
        else:
            print 'ERROR: Unknown system type!'
            print self.SYS
            sys.exit(0)
            

    def setOptimize(self):
        if self.SYS == 'linux':
            self.OPT = '-O3 -felide-constructors -fnonnull-objects -ffast-math'

        elif self.SYS == 'insure':
            self.OPT = '-O3'
            
        elif self.SYS == 'sgi':
            self.OPT = '-O3'
            
        elif self.SYS == 'llnl':
            self.OPT = '-O3'

        elif self.SYS == 'lanl':
            self.OPT = '-O3'

	elif self.SYS == 'cplant':
	    self.OPT = '-O3'

        elif self.SYS == 'ccmalloc':
            self.OPT = '-O3 -felide-constructors -fnonnull-objects -ffast-math'
            
        elif self.SYS == 'default':
            self.OPT = '-O3'
            
        else:
            print 'ERROR: Unknown system type!'
            print self.SYS
            sys.exit(0)

    def setDebug(self):
        if self.SYS == 'linux':
            self.DBG = '-g'

        elif self.SYS == 'insure':
            self.DBG = '-g'
        
        elif self.SYS == 'sgi':
            self.DBG = '-g'
            
        elif self.SYS == 'llnl':
            self.DBG = '-g'

        elif self.SYS == 'lanl':
            self.DBG = '-g'

	elif self.SYS == 'cplant':
	    self.DBG = '-g'

        elif self.SYS == 'ccmalloc':
            self.DBG = '-g'
            
        elif self.SYS == 'default':
            self.DBG = '-g'
            
        else:
            print 'ERROR: Unknown system type!'
            print self.SYS
            sys.exit(0)

    def setParallel(self):
        self.PLL = '-DPARALLEL'

        if self.SYS == 'linux':
            self.CXX = 'mpiCC'

        elif self.SYS == 'insure':
            self.LNK = self.LNK + ' -lm -lmpi++ -lmpio -lmpi -ltstdio ' \
                       + '-ltrillium -largs -lt'
            self.INC = self.INC + ' -I/usr/local/lam-6.3.2/include/mpi2c++ ' \
                       + '-I/usr/local/lam-6.3.2/include -L/usr/local/' \
                       + 'lam-6.3.2/lib'

        elif self.SYS == 'sgi':
            self.LNK = self.LNK + ' -lmpi'

        elif self.SYS == 'ccmalloc':
            self.INC = self.INC + ' -I/usr/local/lam-6.3.2/include/mpi2c++ ' \
                       + '-I/usr/local/lam-6.3.2/include'
            self.LNK = self.LNK + ' -L/usr/local/lam-6.3.2/lib -lmpi++ ' \
                       + '-lmpio -lmpi -ltstdio -ltrillium -largs -lt'

        elif self.SYS == 'lanl':
            self.PLL = self.PLL + ' -DMPI_NO_CPPBIND'
            self.LNK = self.LNK + ' -lmpi'

        elif self.SYS == 'cplant':
            self.INC = self.INC + ' -I/usr/local/cplant/west/current/include'
	    self.LNK = self.LNK + ' -lmpi'

    def setProfile(self):
        if self.SYS == 'linux':
            self.PFL = '-p'

        elif self.SYS == 'insure':
            self.PFL = '-p'
            
        elif self.SYS == 'sgi':
            self.PFL = '-p'
            
        elif self.SYS == 'llnl':
            self.PFL = '-p'

        elif self.SYS == 'lanl':
            self.PFL = '-p'

	elif self.SYS == 'cplant':
	    self.PFL = '-p'

        elif self.SYS == 'ccmalloc':
            self.PFL = '-p'
            
        elif self.SYS == 'default':
            self.PFL = '-p'
            
        else:
            print 'ERROR: Unknown system type!'
            print self.SYS
            sys.exit(0)

    def printString(self):
        text = 'SYS = ' + self.SYS + '\n'
        text = 'VER = ' + self.VER + '\n'
        text = text + 'HOME = ' + self.HOME + '\n'
        text = text + 'INC = ' + self.INC + '\n'
        text = text + 'MAKE = ' + self.MAKE + '\n'
        text = text + 'EXE = ' + self.EXE + '\n'
        text = text + 'DEP = ' + self.DEP + '\n'
        text = text + 'OPT = ' + self.OPT + '\n'
        text = text + 'DBG = ' + self.DBG + '\n'
        text = text + 'PFL = ' + self.PFL + '\n'
        text = text + 'PLL = ' + self.PLL + '\n'
        text = text + 'LINK = ' + self.LNK + '\n'
        text = text + 'CXX = ' + self.CXX + \
               ' -DVERSION=$(VER) $(OPT) $(DBG) $(PFL) $(PLL)\n'
        return text

    def parseInput(self,inputflags):
        for flag in inputflags[2:]:
            if flag == '--optimize':
                self.setOptimize()

            elif flag == '--nooptimize':
                self.OPT = ''

            elif flag == '--debug':
                self.setDebug()

            elif flag == '--nodebug':
                self.DBG = ''

            elif flag == '--profile':
                self.setProfile()

            elif flag == '--noprofile':
                self.PFL = ''

            elif flag == '--parallel':
                self.setParallel()

            elif flag == '--noparallel':
                self.PLL = ''

            elif string.find(flag,'--exe=') != -1:
                self.EXE = flag[6:]

            elif string.find(flag,'--cxx=') != -1:
                self.CXX = flag[6:]

            elif string.find(flag,'--make=') != -1:
                self.MAKE = flag[7:]

            elif string.find(flag,'--optimize-flags=') != -1:
                self.OPT = flag[17:]

            elif string.find(flag,'--debug-flags=') != -1:
                self.DBG = flag[14:]

            elif string.find(flag,'--profile-flags=') != -1:
                self.PFL = flag[16:]

            elif string.find(flag,'--link-flags=') != -1:
                self.LNK = flag[13:]

            else:
                print 'ERROR: Unknown flag ' + flag
                print
                self.printHelp()
                sys.exit(0)


    '''
    Gets the most recent date a file was modified in this directory.  The date
    is returned in the QMcBeaver versioning format (yyyymmdd).
    '''

    def __getMostRecentCVSDateLocalDirectory(self,directory):

        if os.path.isdir(directory):
            directoryContents = os.listdir(directory)
        else:
            directoryContents = []
            
        mostRecentVersionNumber = 0

        # see if there is a CVS directory
    
        if "CVS" in directoryContents:

            # get the most recent file date from the CVS files

            data = open( directory + "CVS/Entries" ).readlines()

            for datum in data:
                temp = string.split(datum,"/")

                if len(temp) > 3:

                    date = temp[-3]
                    date = string.split(date)
                
                    if len(date) > 1:
                        day   = string.atoi(date[2])
                        month = date[1]
                        year  = string.atoi(date[4])

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

        else:
            # no CVS directory so return a date of zero
            pass

        return mostRecentVersionNumber



    '''
    Gets the most recent date a file was modified in all subdirectories.
    The date is returned in the QMcBeaver versioning format (yyyymmdd).
    '''

    def __getMostRecentCVSDate(self,directory):

        # get the version information for the local directory
        mostRecentVersionNumber = \
                        self.__getMostRecentCVSDateLocalDirectory(directory)

        if os.path.isdir(directory):
            directoryContents = os.listdir(directory)
        else:
            directoryContents = []

        for item in directoryContents:
            versionNumber = \
                          self.__getMostRecentCVSDate( directory + item + "/" )

            if versionNumber > mostRecentVersionNumber:
                mostRecentVersionNumber = versionNumber

        return mostRecentVersionNumber


    '''
    Gets the correct version number of the software.
    '''

    def getVersionNumber(self):
        return self.__getMostRecentCVSDate("./")



if __name__ == '__main__':
    import sys

    dat = ControlMake(sys.argv)
    print dat.printString()

    file = open('Makefile.config','w')
    file.write(dat.printString())
    file.close()

    file = open('src/Makefile.dep','w')
    file.close()


















