#!/bin/bash

# Check data directory
DATA_DIR=data
if [[ ! -d ${DATA_DIR}/ ]]
then
	echo "Data directory '${DATA_DIR}/' not found!"
	exit
fi

# Check script directory
SCRIPT_DIR=scripts
if [[ ! -d ${SCRIPT_DIR}/ ]]
then
	echo "Scripts directory '${SCRIPT_DIR}/' not found!"
	exit
fi

# Get Working directory
WD=`pwd`
echo ${WD}

cd ${DATA_DIR}/
for d in 100 300 999
do
	DS=DREAM5_SysGenA${d}
	if [[ ! -d ${DS} ]]
	then
		echo "DREAM5 directory '${DS}' not found!"
		exit
	fi
	pushd ${DS}
	for i in {1..5}
	do
		NET=${DS}_Network${i}
		# Prepare input and run tools
		Rscript ${WD}/${SCRIPT_DIR}/analysis_script.R \
			${NET}_Genotype.tsv \
			${NET}_Expression.tsv \
			${NET} \
			${WD}/${DATA_DIR}/${DS}/tmp.output${i}/ \
			${WD}/${SCRIPT_DIR}/ \
		# Join Results in a matrix for Machine Learning
		python ${WD}/${SCRIPT_DIR}/joinRes.py \
			-a ${NET}_MatrixEQTL_ANOVA.tsv \
			-l ${NET}_MatrixEQTL_LINEAR.tsv \
			-m ${NET}_mRMR.tsv \
			-e ${NET}_rqtl_em.csv \
			-k ${NET}_rqtl_hk.csv \
			-t Truth/${DS}_Edges_Network${i}.tsv \
			-o ${NET}_matrix.csv
		gzip ${NET}_matrix.csv
	done
	popd
done
