#!/usr/bin/bash
##############################################
##############################################
# Code to design primers for MEI validation
# Rebeca Borges Monroy

##############################################
##############################################

################### Parameters ################
###### To change
# Can be hg38 or hg37
Genome_version=hg38
#Genome_version=hs37d5
# I already ran R code Get_bed12.R and sorted.
PRIMER_DESIGN_PATH=""
OUTDIR=""

INPUT_FILES_DIR=${PRIMER_DESIGN_PATH}/input_files
INPUT_CODE_DIR=${PRIMER_DESIGN_PATH}/code
Input_sorted_bed=${INPUT_FILES_DIR}/HG002_ManuallyInspected_Polymorphic_10_4_20_masked_sorted.bed12

Path_fasta=${PRIMER_DESIGN_PATH}/input_files/fasta/hg38.fa.masked
# primer3 settings file
PRIMER3_INPUT_PATH=${INPUT_FILES_DIR}/primer3_input.txt
#isPCR path
isPCR_path=${PRIMER_DESIGN_PATH}/code/InSilico

################################################
############ O2 modules to load ################
################################################

module load primer3
module load gcc/6.2.0
module load blat/35

#primer3_core [-format_output] [-default_version=1|-default_version=2] [-io_version=4] [-p3_settings_file=<file_path>] [-echo_settings_file] [-strict_tags] [-output=<file_path>] [-error=<file_path>] [input_file]
# current version:
# primer3_core --about
# libprimer3 release 2.3.7

###############################################
### Get Fasta Sequence

# Obtain up and downstream regions from breakpoints (-50 to -800 5' and +50 to +800 3')
# Uses the following input:
# FAMILY	SAMPLE_ID	TYPE	CHR	START	END	PHASE	FILTER_DENOVO	CONFIDENCE	Gene	Pli	SAFARI
#HG002	NA24385	L1	chr14	30681602	30681617	NA	two_side_tprt_both	2	NA	NA	NA
#HG002	NA24385	L1	chr5	90154952	90154962	NA	two_side_tprt_both	2	NA	NA	NA
#HG002	NA24385	L1	chr3	38584574	38584591	NA	two_side_tprt_both	2	NA	NA	NA

# Rscript ${INPUT_CODE_DIR}/Get_bed12.R

# Output example:
#chr14	30680802	30682417	HG002_NA24385_L1_chr14_30681602	0	+	30681602	30681602	NA	2	700,700	0,915
#chr5	90154152	90155762	HG002_NA24385_L1_chr5_90154952	0	+	90154952	90154952	NA	2	700,700	0,910
#chr3	38583774	38585391	HG002_NA24385_L1_chr3_38584574	0	+	38584574	38584574	NA	2	700,700	0,917

# Change input in code to get your coordinates.

# Get fasta coordinates, read each line of input file.

mkdir -p ${OUTDIR}/fasta
mkdir -p ${OUTDIR}/Primer3
mkdir -p ${OUTDIR}/isPCR
mkdir -p ${OUTDIR}/Blat

while read line
do
	export ID=`echo $line | awk '{print $4}'`
	export CHR=`echo $line | awk '{print $1}'`

	echo -e "# Working on ${ID}"

	############################ get fasta sequence for coordinates  ##########################
	echo $line | awk '{print $0}' | sed 's/ /\t/g' > ${OUTDIR}/fasta/${ID}.bed12.input

	bedtools getfasta -fi ${Path_fasta} -bed ${OUTDIR}/fasta/${ID}.bed12.input -split  > ${OUTDIR}/fasta/${ID}.fasta

	# Check if fasta file found
	if [[ -f ${OUTDIR}/fasta/${ID}.fasta ]] ; then
		:
	else
		echo -e "#ERROR ${OUTDIR}/fasta/${ID}.fasta empty "
	continue
	fi ;

	#rm ${OUTDIR}/fasta/${ID}.bed12.input

	############################ Run Primer3 to get primers  ##########################
	 # Primer3 input format:
	 #SEQUENCE_ID=chr1test
	 #SEQUENCE_TEMPLATE=SequenceFromGetFasta
	 # Output in this format for isPCR:

	 # Name FWDPrimer REVERSEPrimer
	 #TEST_PRIMER	GGTATCTTTATCACCGCGACT	aacctccctttcccaggttc

	export SEQ_ID_INPUT_PRIMER3=${ID}
	# Fasta code with region sequence for up and downstream breakpoint
	export SEQ_INPUT_PRIMER3=`tail -n 1  ${OUTDIR}/fasta/${ID}.fasta`

	# save ID and sequence in file
	sh ${INPUT_CODE_DIR}/primer3_header_input.sh  > ${OUTDIR}/Primer3/${ID}_primer3_input_header.out
	cat ${OUTDIR}/Primer3/${ID}_primer3_input_header.out ${PRIMER3_INPUT_PATH}  > ${OUTDIR}/Primer3/${ID}_primer3.input

	primer3_core ${OUTDIR}/Primer3/${ID}_primer3.input --output=${OUTDIR}/Primer3/${ID}_primer3.out

	rm ${OUTDIR}/Primer3/${ID}_primer3_input_header.out
	rm ${OUTDIR}/Primer3/${ID}_primer3.input

	# Check if primer output found
	if [[ -f ${OUTDIR}/Primer3/${ID}_primer3.out ]] ; then
		:
	else
		echo -e "#ERROR ${OUTDIR}/Primer3/${ID}_primer3.out empty "
	continue
	fi ;

	############################ Run isPCR  ##########################
	# Make isPCR input file from Primer3 output

	for i in {0..4}
	do
		#echo $i
		FORWARD=`grep "PRIMER_LEFT_${i}_SEQUENCE=" ${OUTDIR}/Primer3/${ID}_primer3.out | sed  "s/PRIMER_LEFT_${i}_SEQUENCE=//g"`
		REVERSE=`grep "PRIMER_RIGHT_${i}_SEQUENCE=" ${OUTDIR}/Primer3/${ID}_primer3.out | sed  "s/PRIMER_RIGHT_${i}_SEQUENCE=//g"`

		[[ -z "$FORWARD" || -z "$REVERSE" ]] && echo "# ${ID} Primers ${i} NA" || echo -e ${ID}_${i}' \t '${FORWARD}' \t '${REVERSE}

	done > ${OUTDIR}/isPCR/${ID}_isPCR_input_com.txt

	# remove comments for no primer
	grep -v "#" ${OUTDIR}/isPCR/${ID}_isPCR_input_com.txt > ${OUTDIR}/isPCR/${ID}_isPCR.input

	# Check if primers output found
	if [[ -f ${OUTDIR}/isPCR/${ID}_isPCR.input ]] ; then
		:
	else
		echo -e "#ERROR ${OUTDIR}/isPCR/${ID}_isPCR.input empty "
	continue
	fi ;

	# run isPCR
	${isPCR_path}/isPcr ${isPCR_path}/twoBit/${Genome_version}_${CHR}.2bit ${OUTDIR}/isPCR/${ID}_isPCR.input ${OUTDIR}/isPCR/${ID}_isPCR.output

	# Check if isPCR output exists
	if [[ -f ${OUTDIR}/isPCR/${ID}_isPCR.output ]] ; then
		:
	else
		echo -e "#ERROR ${OUTDIR}/isPCR/${ID}_isPCR.output empty "
	continue
	fi ;

	# check whether there is only one match
	while read line
	do

		export PRIMER_NAME=`echo $line | awk '{print $1}'`
		export IS_MATCHES=`grep "${PRIMER_NAME}" ${OUTDIR}/isPCR/${ID}_isPCR.output | wc -l`

		if [ "$IS_MATCHES" == 1 ]; then
			echo ${line}
		fi;

	done < ${OUTDIR}/isPCR/${ID}_isPCR.input > ${OUTDIR}/isPCR/${ID}_isPCR_uniqueMatches.output

	# Check if output exists
	if [[ -f ${OUTDIR}/isPCR/${ID}_isPCR_uniqueMatches.output ]] ; then
		:
	else
		echo -e "#ENDED ${OUTDIR}/isPCR/${ID}_isPCR_uniqueMatches.output does not have unique matches "
	continue
	fi ;

	############################ Run Blat  ##########################

	# make input fasta file for blat.
	awk -F' ' '{ print ">"$1"_left""\n"$2"\n"">"$1"_right""\n"$3}' ${OUTDIR}/isPCR/${ID}_isPCR_uniqueMatches.output > ${OUTDIR}/isPCR/${ID}_isPCR_uniqueMatches.output.fasta

	# Run blat
	blat -stepSize=5 -minScore=20 -minIdentity=80 ${isPCR_path}/twoBit/${Genome_version}.2bit ${OUTDIR}/isPCR/${ID}_isPCR_uniqueMatches.output.fasta ${OUTDIR}/Blat/${ID}_blat.output.psl

	# Check if output exists
	if [[ -f ${OUTDIR}/Blat/${ID}_blat.output.psl ]] ; then
		:
	else
		echo -e "#ERROR ${OUTDIR}/Blat/${ID}_blat.output.psl is empty "
	continue
	fi ;

	# Check if there are primers with only one blat match.

	while read line
	do

		export PRIMER_NAME=`echo $line | awk '{print $1}'`
		export IS_MATCHES_LEFT=`grep "${PRIMER_NAME}_left" ${OUTDIR}/Blat/${ID}_blat.output.psl | wc -l`
		export IS_MATCHES_RIGHT=`grep "${PRIMER_NAME}_right" ${OUTDIR}/Blat/${ID}_blat.output.psl | wc -l`

		if [[ "$IS_MATCHES_LEFT" == 1 && "$IS_MATCHES_RIGHT" == 1 ]]; then
			echo ${line}
		fi;

	done < ${OUTDIR}/isPCR/${ID}_isPCR.input > ${OUTDIR}/Blat/${ID}_Blat_uniqueMatches.output

done < ${Input_sorted_bed}




