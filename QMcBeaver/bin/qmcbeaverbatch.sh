#!/bin/sh

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

# Be sure to set the QMCBEAVER_HOME and QMCBEAVER_EXE environment variables

for file do
    echo "file: $file"
    logfile=$file'.log'
    echo "log:$logfile"

    echo "Job started: `date `" > $logfile
    echo "Machine: $HOST" >> $logfile
    echo "" >> $logfile
    echo "Files To Run" >> $logfile
    echo "------------" >> $logfile
    echo "" >> $logfile

    echo "file: $file"

    for item in `cat $file` 
    do
	echo "$item"
	echo "$item" >> $logfile
    done

    echo "" >> $logfile
    echo "" >> $logfile
    echo "Running QMcBeaver" >> $logfile
    echo "-----------------" >> $logfile
    echo "" >> $logfile

    for item in `cat $file` 
    do
	echo "$item: Running" >> $logfile
	$QMCBEAVER_HOME/$QMCBEAVER_EXE $item
	echo "$item: Done" >> $logfile
    done
    echo "" >> $logfile
    echo "Job Finished: `date`" >> $logfile
done



