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
	@echo "Updating object files..."
	@-cd src; $(MAKE) 

# this works on some stubborn machines
bruteforce:
	@echo "Brute force make procedure..."
	-cd src; $(CXX) $(OPTIMIZATION) $(INCLUDE) -o ../bin/$(EXE) *.cpp $(LINK)

.PHONY:dep
dep: src/$(DIROBJ)/Makefile.dep
	@echo "Updating dependency files..."
	@mkdir -p src/$(DIROBJ)
	@-cd src; $(MAKE) $(DIROBJ)/Makefile.dep

src/$(DIROBJ)/Makefile.dep:

#when executables are created the first time, the link command
#was saved. this one command allows one to update all the different
#executables created on the computer from all the different options
#(e.g. debug, parallel) assuming that the dependency files haven't
#changed.
.PHONY:updateall
updateall:
	@echo "Updating all $(SYS) executables..."
	@-cd src; $(MAKE) updateall

test: all
	-cd test; $(MAKE)

#this target will remove the object files and dependency files
#associated with the current Makefile.config
.PHONY:clean
clean:
	-cd src; $(MAKE) clean
	rm -f *~
	rm -f bin/$(EXE)

.PHONY:execlean
execlean:
	-cd bin
	-rm -f *.$(LABEL)

.PHONY:testclean
testclean: 
	@cd test; $(MAKE) clean

cleanall: clean execlean testclean 
	rm -rf src/obj_*


