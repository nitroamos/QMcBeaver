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




include Makefile.config

# generic make procedure
all: dep 
	@echo "Compiling..."
	-cd src; $(MAKE) 

# this works on some stubborn machines
bruteforce:
	@echo "Brute force make procedure..."
	-cd src; $(CXX) $(OPT) $(INC) -o ../bin/$(EXE) *.cpp $(LINK)


.PHONY:dep
dep: src/Makefile.dep
	@echo "Resolving dependencies..."
	@mkdir -p src/d.files
	-cd src; $(MAKE) Makefile.dep

test: all
	-cd test; $(MAKE)

clean:
	-cd src; $(MAKE) clean
	rm -f *~

execlean:
	rm -f bin/QMcBeaver.* bin/pQMcBeaver.*

testclean: 
	-cd test; $(MAKE) clean

cleanall: clean execlean testclean 

