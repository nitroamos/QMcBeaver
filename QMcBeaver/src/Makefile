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
include ../Makefile.config

OBJECTS  := $(patsubst %.cpp,$(DIROBJ)/%.o,$(wildcard *.cpp))
DEPFILES := $(patsubst %.cpp,$(DIROBJ)/%.d,$(wildcard *.cpp))

all: $(DIROBJ).e

#if we selected the updateall target, then we can't allow gmake to recreate the
#dependency files because the config dependent options are all different
ifneq ($(MAKECMDGOALS),updateall)

$(DIROBJ)/%.d: %.cpp ../include/%.h
	@rm -f $@
	@echo "Writing $*.d"
	@$(CXX) $(INCLUDE) $(DEPENDENCY) $< | sed -e "s/$*\.o/$(DIROBJ)\/$*\.o/g" >> $@
	@echo "	@echo Compiling $*.o" >> $@
	@echo "	@$(CXX) $(OPTIMIZATION) $(INCLUDE) -c -o $(patsubst %.cpp,$(DIROBJ)/%.o,$<) $*.cpp" >> $@
	@echo "" >> $@
	@echo "" >> $@

$(DIROBJ)/QMcBeaver.d: QMcBeaver.cpp
	@rm -f $@
	@echo "Writing QMcBeaver.d"
	@$(CXX) $(INCLUDE) $(DEPENDENCY) $< | sed -e "s/QMcBeaver\.o/$(DIROBJ)\/QMcBeaver\.o/g" >> $@
	@echo "	@echo Compiling QMcBeaver.o" >> $@
	@echo "	@$(CXX) $(OPTIMIZATION) $(INCLUDE) -c -o $(DIROBJ)/QMcBeaver.o QMcBeaver.cpp" >> $@
	@echo "" >> $@
	@echo "" >> $@

$(DIROBJ)/Makefile.dep: $(DEPFILES) ../Makefile.config
	@cat $(DEPFILES) > $(DIROBJ)/Makefile.dep
	@echo "$(DIROBJ).e: $(OBJECTS)" >> $@
	@echo "	@echo Creating $(EXE)..." >> $@
	@echo "	@$(CXX) $(OPTIMIZATION) $(INCLUDE) -o ../bin/$(EXE) $(OBJECTS) $(LINK)" >> $@

#This if block makes sure not to update the .d files right before we delete them
ifneq ($(MAKECMDGOALS),clean)
-include $(DIROBJ)/Makefile.dep
endif

else

DIRSET := $(wildcard *$(TAG))
include $(addsuffix /Makefile.dep, $(DIRSET))

endif

.PHONY:updateall
updateall: $(addsuffix .e,$(DIRSET))

.PHONY:clean
clean:
	rm -f *.o *~ *#
	rm -rf ii_files d.files cxx_repository $(DIROBJ)



