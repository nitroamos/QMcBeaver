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

include ../Makefile.config

OBJECTS  := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
DEPFILES := $(patsubst %.cpp,d.files/%.d,$(wildcard *.cpp))

all: QMcBeaver.e

QMcBeaver.e: $(OBJECTS) $(LIBS)
	$(CXX) $(INC) -o ../bin/$(EXE) $^ $(LINK)

Makefile.dep: $(DEPFILES)
	@cat $(DEPFILES) > Makefile.dep

include Makefile.dep

.cpp.o:
	$(CXX) $(INC) -c $<

d.files/%.d: %.cpp ../include/%.h
	@$(CXX) $(INC) $(DEP) $< >> $@
	@echo "" >> $@
	@echo "" >> $@

d.files/QMcBeaver.d: QMcBeaver.cpp
	@$(CXX) $(INC) $(DEP) $< >> $@
	@echo "" >> $@
	@echo "" >> $@

clean:
	rm -f *.o *~ 
	rm -rf ii_files d.files






