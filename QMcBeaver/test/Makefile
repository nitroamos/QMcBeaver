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

CKMF_FILES  :=  $(patsubst %.ckmf,%.qmc,$(wildcard *.ckmf))

all: clean $(CKMF_FILES)

%.qmc: %.ckmf
	@echo $*
	../bin/$(EXE) $<

clean:
	rm -f *.01.* *.02.* *.qmc *.rslts *.gnp *~ *.run *.o*

