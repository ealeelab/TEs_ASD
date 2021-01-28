# Supporting Code for "Whole-genome analysis of de novo and polymorphic retrotransposon insertions in Autism Spectrum Disorder "

## Primer Design

The Primer_design folder contains the folder code, which has the following 3 files:

1. Get_bed12.R : This code is used to convert insertion candidates to bed12 format. See and edit code for paths and input file format. Sort output with bedtools sort before continuing. 
2. Primer_design.sh : Contains the main code used to design primers. Edit the parameters section with the appropriate paths. 
3. primer3_header_input.sh : Is a supporting file used by Primer_design.sh

The input folder contains the following files:

1. Primer3_misprimingLibrary_human.txt : Is used by primer3 as a repeat or mispriming library. 
2. primer3_input.txt : Are the primer3 settings, change these if you want to change any parameters used. 

## Merged TE intervals and Allele Frequencies

The SSC_TEs folder contains the merged intervals for parental insertions and all insertions (including unaffected and affected children) in the SSC cohort. These were used to count how many unique insertions were in the cohort. 
The format for the TEFAMILY_mergedRegions.txt.gz files are:

chromosome, interval start, interval end, samples with insertions in these regions, families with insertions in these regions, estimated insertions size, and genotype (1 for heterozygous and 2 for homozygous)
chr1	776261	776341	SSC09486	13709	101	1

Parental merged regions also have the population allele frequency based on the number of chromosomes carrying the insertion in the population. 

## De novo Insertions Code

The example_denovo_code folder contains the code used:

1. Parameters_run_code.txt : Contains the paths needed to run the code. 
2. run_denovo_TE_trio_pipeline.sh : Used in SLURM to run the main code
3. denovo_TE_trio_pipeline.sh : Main pipeline

The input_files folder are input files used for de novo calling:

1. Example_input.txt : One line input example, other samples would be on an additional line in the same format. 
2. ReferenceYoung_KNR_1000G_SSCParents_*.hg38.bed : Contains the intervals used to filter our young reference insertions, KNR insertions from 1000 genome, previous studies, and a subset of parental samples from this cohort
3. SFARI_gene_score.txt : Contains the scores for confidence levels for SFARI genes used, and whether they are syndromic genes.
4. SFARI_genes_hg38.bed : Contains the list of SFARI genes and hg38 positions
5. denovo_breakpointmargins.R : Is the code used to obtain breakpoint margins and bed format for insertions
6. hg38.fa.out.*.chr.bed : Reference TE annotations coordinates
7. hg38_refseq_geneannot_sort.txt.gz : Refseq gene annotation used
8. pLI.txt : pLI scores used
9. xTea_parent_merged_clip.R : Supporting code for main pipeline to merge parental clipped reads. 

The DENOVO_EXAMPLE_L1_01_26_21 folder contains the output of the input example provided. 

The "p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF_Original.txt" file contains the 2 de novo candidates which do not have clipped reads or discordant reads in parents at the threshold provided, and do not overlap with reference insertions. In this case one seemed like a parental false negative insertion and the other one like a true candidate upon manual inspection on IGV.  

The "p1_HighConfInsertions_noDISCparent_noCLIPparent_Original.txt" file is the same as described above, but includes insertions which overlap older reference elements that were not present in the reference and KNR filer file.
