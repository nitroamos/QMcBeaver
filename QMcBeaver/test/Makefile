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

CKMF_FILES  :=  $(wildcard *.ckmf)

all: clean runall_results.txt

runall_results.txt :
	../bin/runqmcbeavers.pl $(CKMF_FILES)

clean:
	rm -f *.01.* *.02.* *.energy.* *.qmc *.rslts *.force *.gnp *~ *.run *.o* runall_results.txt

