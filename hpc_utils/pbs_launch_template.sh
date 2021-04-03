#!/bin/bash

#PBS -l place=free
#PBS -V
#PBS -l walltime=24:00:00

## This is an example script to be modified for a specific task

workUnitPrefix="VDBM";
mem_req="32gb";
outFile_extension='vsc.tsv';
workingFile_extension='wk';
target_DataType='reads'

if [ ! $ncores ]; then ncores="2"; fi;

NODE_TEMPFOLDER=/tmp/
SERVER_BASEFOLDER=${oDir}/${dn}/

# Temp Folder Selection

FREESPACE=$(df --output=avail $NODE_TEMPFOLDER | tail -n1);
WORKSIZE=$(du -L -ck ${prefix}/${dn}/${target_DataType}/${fn}/*${extension} | tail -n1 | cut -f1)
REQSIZE=$((WORKSIZE*6))

if [[ $FREESPACE -lt $REQSIZE ]]; then 
	TEMPFOLDER= #SELECT 1
else
	TEMPFOLDER= #SELECT 2
fi

if [ -d ${TEMPFOLDER} ]; then 
	rm -r ${TEMPFOLDER}/;
fi;

mkdir -p $TEMPFOLDER/ 
mkdir -p $SERVER_BASEFOLDER/

# Command to launch

$uncompress_cmd ${prefix}/${dn}/${target_DataType}/${fn}/*${extension} > ${TEMPFOLDER}/<OUT>


# Remove working file: the task is completed 

if [ -f ${TEMPFOLDER}/${fn}.${outFile_extension} ]; then 
	rm ${SERVER_BASEFOLDER}/${dn}__${fn}.${workingFile_extension}
fi;

# Move the data to final location

mv ${TEMPFOLDER}/*.tsv ${SERVER_BASEFOLDER}/;

# Clean

rm -r ${TEMPFOLDER};












