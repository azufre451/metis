#!/bin/bash

## This is a custom script for UniTN HPC PBS
##
## Multi Scheduler. This is basically a while loop on a list of samples 
## organized as in the CIBIOCM Pipelines schema:
## <PREFIX>/<DATASET_NAME>/<target_DataType>/<SAMPLE_N>/...
## a scipt "launch.sh" (has to be written for each task!) is scheduled on PBS
## for each folder (i.e. sample) in the dataset
##
##
## Intended use on things like:
## <PREFIX>/HMP_2012/reads/SRS0000001/ (e.g launch metaphlan on each)


curDir=$(realpath $(dirname $0));
my_user=$(whoami);

## can be used to define packets of data to work on
#  source ${curDir}/workunit_data.sh;

workUnitPrefix="VDBM";
mem_req="32gb";
outFile_extension='vsc.tsv';
workingFile_extension='wk';
target_DataType='reads'

##############################################
 
mode=$1;

if [ $mode == "U" ]; then
	extension='.fastq'
	extraction_cmd='cat'
fi;

if [ $mode == "bzip" ]; then
	extension='.fastq.bz2'
	extraction_cmd='bzcat'
fi;
	
if [ $mode == "gzip" ]; then
	extension='.fastq.gz'
	extraction_cmd='zcat'
fi;

if [ $mode == "fqgz" ]; then
	extension='.fq.gz'
	extraction_cmd='zcat'
fi;

# build the list of folders to schedule

for et in $etp; do
	for utp in $(find ${prefix}/${et}/${targetDataType} -maxdepth 1 -mindepth 1 -type d); do
		e+=($utp);
	done;

	mkdir -p ${logFolder}/${et}/;
	mkdir -p ${odir}/${et}/

done;

while true; do
	
	bt=0;
	b_short=$(qstat -u ${my_user} | grep short | grep ${workUnitPrefix} | wc -l);
	b_common=$(qstat -u ${my_user} | grep common | grep ${workUnitPrefix} | wc -l);
	b_cibio=$(qstat -u ${my_user} | grep CIBIO_ | grep ${workUnitPrefix} | wc -l);
	b_cibiocm=$(qstat -u ${my_user} | grep CIBIOCM | grep ${workUnitPrefix} | wc -l);

	
	for folder in ${e[@]}; do
		
		sample=$(basename $folder);
		dataset=$(basename $(dirname $(dirname $folder)));
		odirE=${odir}/${dataset}/


		if [ ! -f ${odir}/${dataset}/${sample}.${outFile_extension} ] && [ ! -f ${odir}/${dataset}/${dataset}__${sample}.${workingFile_extension} ]; then
 
			if [ $b_short -lt 30 ]; then 

				echo "QUEUE (short)" ${dataset} ${sample};
				qsub -q short_cpuQ -v prefix=\"${base}\",dn=\"${dataset}\",fn=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"8\",oDir=\"${odir}\" -N ${workUnitPrefix}_${dataset}_${sample} -o ${logFolder}/${dataset}/${workUnitPrefix}_${dataset}_${sample}.o.log -e ${logFolder}/${dataset}/${workUnitPrefix}_${dataset}_${sample}.e.log -l select=1:ncpus=8:mem=${mem_req} ${curDir}/launch.sh;
				touch ${odir}/${dataset}/${dataset}__${sample}.wk
				b_short=$((b_short+1));
				bt=$((bt+1));

			elif [ $b_common -lt 60 ]; then 

				echo "QUEUE (common)" ${dataset} ${sample};
				qsub -q common_cpuQ -v prefix=\"${base}\",dn=\"${dataset}\",fn=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"2\",oDir=\"${odir}\" -N ${workUnitPrefix}_${dataset}_${sample} -o ${logFolder}/${dataset}/${workUnitPrefix}_${dataset}_${sample}.o.log -e ${logFolder}/${dataset}/${workUnitPrefix}_${dataset}_${sample}.e.log -l select=1:ncpus=2:mem=${mem_req} ${curDir}/launch.sh;
				touch ${odir}/${dataset}/${dataset}__${sample}.wk
				b_common=$((b_common+1));
				bt=$((bt+1));

			elif [ $b_cibio -lt 25 ]; then 

				echo "QUEUE (common)" ${dataset} ${sample};
				qsub -q CIBIO_cpuQ -v prefix=\"${base}\",dn=\"${dataset}\",fn=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"2\",oDir=\"${odir}\" -N ${workUnitPrefix}_${dataset}_${sample} -o ${logFolder}/${dataset}/${workUnitPrefix}_${dataset}_${sample}.o.log -e ${logFolder}/${dataset}/${workUnitPrefix}_${dataset}_${sample}.e.log -l select=1:ncpus=2:mem=${mem_req} ${curDir}/launch.sh;
				touch ${odir}/${dataset}/${dataset}__${sample}.wk
				b_cibio=$((b_cibio+1));
				bt=$((bt+1));


			elif [ $b_cibiocm -lt 25 ]; then 

				echo "QUEUE (common)" ${dataset} ${sample};
				qsub -q CIBIOCM_cpuQ -v prefix=\"${base}\",dn=\"${dataset}\",fn=\"${sample}\",uncompress_cmd=\"${extraction_cmd}\",extension=\"${extension}\",ncores=\"2\",oDir=\"${odir}\" -N ${workUnitPrefix}_${dataset}_${sample} -o ${logFolder}/${dataset}/${workUnitPrefix}_${dataset}_${sample}.o.log -e ${logFolder}/${dataset}/${workUnitPrefix}_${dataset}_${sample}.e.log -l select=1:ncpus=2:mem=${mem_req} ${curDir}/launch.sh;
				touch ${odir}/${dataset}/${dataset}__${sample}.wk
				b_cibiocm=$((b_cibiocm+1));
				bt=$((bt+1));

			fi;
#		else
#			echo "SKIP" ${dataset} ${sample};
		fi;

	done;

	a=$(qstat -u ${my_user} | grep ${workUnitPrefix} | grep " R " | wc -l);
	q=$(qstat -u ${my_user} | grep ${workUnitPrefix} | grep " Q " | wc -l);
	ete=${#e[@]}
	ft=$(find ${odir}/ -name "*."${outFile_extension} | wc -l);
	
	##tg=#$(($ete-$ft));
	echo $(date) "| " $bt "jobs launched, " $a "running ", $q "queued " $ft" done. See you in one minute." ;
	sleep 60
	b=0;
done;
