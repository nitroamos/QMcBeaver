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

include Makefile.config

# generic make procedure
all: dep 
	@echo "Compiling..."
	-cd src; $(MAKE) 

# this works on some stubborn machines
bruteforce:
	@echo "Brute force make procedure..."
	-cd src; $(CXX) $(INC) -o ../bin/$(EXE) *.cpp $(LINK)


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
	rm -f bin/QMcBeaver.*

libclean:
	-cd optimization; $(MAKE) cleanall

testclean: 
	-cd test; $(MAKE) clean

cleanall: clean execlean testclean #libclean

