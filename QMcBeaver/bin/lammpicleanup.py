#!/usr/bin/env python

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

import os
import sys
import string


print "LAMCLEAN*********************************************"
os.system("lamclean -v")



print "WIPE*********************************************"
os.system("wipe -v lamhosts")


print "DONE"

