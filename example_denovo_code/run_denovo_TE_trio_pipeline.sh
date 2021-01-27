#!/bin/bash

#SBATCH -c 1                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-00:10                         # Runtime in D-HH:MM format
#SBATCH -p priority                           # Partition to run in
#SBATCH --mem=22GB                          # Memory total in MB (for all cores)

sh ${PATH_LOCAL_FILES}/example_denovo_code/denovo_TE_trio_pipeline.sh -i ${INPUT_FILE} -n ${OUTDIR} -f ${RSCRIPTS_DIR} -t ${CALLSDIR} -s ${SNAPSHOT_DIR} -e ${SESSIONS_DIR} -m ${PARENTAL_RAM} -p ${PARENTAL_CLIP_TE} -a ${PARENTAL_CLIP_ALL} -c ${MIN_CLIP} -r ${MIN_RAM} -k ${REF_KNR_DIR} -b ${REF_ALL_DIR} -d ${GENE_ANNOT_PATH} -g ${SFARI_ANNOT_PATH} -l ${PLI_ANNOT_PATH} -j ${SFARI_SCORES_PATH}

