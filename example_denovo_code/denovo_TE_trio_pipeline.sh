#!/usr/bin/bash
## Code by Rebeca Borges-Monroy

#July 26, 2018
#Code for de novo calling based on xTea output
#September 17, 2018 Added clipped and discordant reads filter and filtering out reference and knr calls
#September 20, 2018 Added check for sample names
#Semptember 21, 2018 edit code to take name as input and not bam path
# February 20, 2019 edit code to save entire coordinates with new coordinates, to make it easier to trace back calls. Also now using midpoint of clipped position instead of just a margin
# March 8, 2019 CHanged snapshot output
# March 11, 2019 Code for High confidence, changed margins for plotting, and added grep for original call
# March 18, 2019 Added gene annotation and now reference is not filtered, but annotated at the end.
# January 25, 2021 Edited code for sharing, changed to smaller margins used in paper, deleted some commented sections.

usage="$(basename "$0") [-h] [-i file -n path -f path -t path -s path -e path -k path -c n -r n -m n -p n -b path -d path -g path -l path -j path] -- program to call de novo te insertions in trios, v2 for new tea output August 2017

where:
    -h  show this help text
    -i input file with samples to run
    -n output path for de novo calling
    -f path with R scripts
    -t path where tea ran
    -s path for saving snapshots
    -e path for saving sessions in xml for IGV
    -k path of file with knr/reference calls
    -c number of clipped reads minimum in call
    -r number of ram/discordant reads minimum in call
    -m parental ram: number of allowed discordant reads total for parents in interval of progeny call
    -p parental clipped TE reads: number of allowed clipped reads total that map to TEs for parents in interval of progeny call
    -a parental clipped TE reads: number of allowed clipped reads total for parents in interval of progeny call
    -b reference annotation path with annotated reference sites
    -d gene annotation path with annotation of genes and genic regions in bed format
    -g coordinates for SFARI genes
    -l pLI scores
    -j SFARI scores file"

while getopts ":i:n:f:t:s:e:k:c:r:m:p:a:b:d:g:l:j:h" opt; do
  case $opt in
    i) input="$OPTARG"
    ;;
    n) dn_path="$OPTARG"
    ;;
    f) rfiles_path="$OPTARG"
    ;;
    t) calls_path="$OPTARG"
    ;;
    s) snapshot_path="$OPTARG"
    ;;
    e) sessions_path="$OPTARG"
    ;;
    k) knrref_path="$OPTARG"
    ;;
    m) parental_disc="$OPTARG"
    ;;
    p) parental_TEclip="$OPTARG"
    ;;
    a) parental_Allclip="$OPTARG"
    ;;
    c) clip_req="$OPTARG"
    ;;
    r) ram_req="$OPTARG"
    ;;
    b) ref_all_path="$OPTARG"
    ;;
    d) gene_annot_path="$OPTARG"
    ;;
    g) SFARI_annot_path="$OPTARG"
    ;;
    l) pLI_annot_path="$OPTARG"
    ;;
    j) SFARI_scores_path="$OPTARG"
    ;;

    h) echo "$usage"
    exit
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
    :)
    echo "Option -$OPTARG requires an argument." >&2
    exit 1
    ;;
  esac
done



while read line
do
export familydname=`echo $line | awk '{print $2}'`

mkdir -p "$dn_path/$familydname"
#echo -e "$familydname"
#echo -e "$dn_path/$familydname"
#export bamASDpath=`echo $line | awk '{print $3}'`
export bamASDpath=`echo $line | awk '{print $1}'`
#echo -e "$bamASDpath"
export bamASDsampleName=`echo $line | awk '{print $1}'`
#echo -e "$bamASDsampleName"

#obtain sample name

#link ASD bai file to raw files path
# add father bam name
#export bamFATHERpath=`echo $line | awk '{print $8}'`
export bamFATHERpath=`echo $line | awk '{print $3}'`
#echo -e "$bamFATHERpath"
export bamFATHERsampleName=`echo $line | awk '{print $3}'`
#echo -e "$bamFATHERsampleName"
#filename without extension
#export bamFATHERsampleName=`basename $bamFATHERpath .final.bam`

# add mother bam and bai file
#export bamMOTHERpath=`echo $line | awk '{print $9}'`
export bamMOTHERpath=`echo $line | awk '{print $4}'`
#echo -e "$bamMOTHERpath"
export bamMOTHERsampleName=`echo $line | awk '{print $4}'`
#echo -e "$bamMOTHERsampleName"
#filename without extension
#export bamMOTHERsampleName=`basename $bamMOTHERpath .final.bam`


#Check if names are NA
if [[ $bamASDsampleName=="NA"  ]] ; then
:
else
echo -e "#Error bamASDsampleName is NA "
continue
fi ;

if [[ $bamFATHERsampleName=="NA"  ]] ; then
:
else
echo -e "#Error bamFATHERsampleName is NA "
continue
fi ;

if [[ $bamMOTHERsampleName=="NA"  ]] ; then
:
else
echo -e "#Error bamMOTHERsampleName is NA "
continue
fi ;


# check if xtea files were run


if [[ -f $calls_path/$bamASDsampleName.final/candidate_disc_filtered_cns.txt.high_confident.post_filtering2.txt ]] ; then
:
else
echo -e "#ERROR $calls_path/$bamASDsampleName.final/candidate_disc_filtered_cns.txt.high_confident.post_filtering2.txt not found "
continue
fi ;

if [[ -f $calls_path/$bamFATHERsampleName.final/discordant_reads_tmp0  ]] ; then
:
else
echo -e "#ERROR $calls_path/$bamFATHERsampleName.final/discordant_reads_tmp0 not found "
continue
fi ;

if [[ -f $calls_path/$bamMOTHERsampleName.final/discordant_reads_tmp0 ]] ; then
:
else
echo -e "#ERROR $calls_path/$bamMOTHERsampleName.final/discordant_reads_tmp0 not found "
continue
fi ;

if [[ -f $calls_path/$bamFATHERsampleName.final/clip_reads_tmp0   ]] ; then
:
else
echo -e "#ERROR $$calls_path/$bamFATHERsampleName.final/clip_reads_tmp0 not found "
continue
fi ;

if [[ -f $calls_path/$bamMOTHERsampleName.final/clip_reads_tmp0 ]] ; then
:
else
echo -e "#ERROR $calls_path/$bamMOTHERsampleName.final/clip_reads_tmp0 not found "
continue
fi ;

#Create bed, first obtain pbp and nbp with margins considering if onesided or not then keep only pbp and nvp


# High Confidence Calls
## Change this if you want to run on all xTea calls, or another filtered file. Col. 6 and 7 are number of left and right clipped reads, 8 and 9 number of left and right discordant reads
cat $calls_path/$bamASDsampleName.final/candidate_disc_filtered_cns.txt.high_confident.post_filtering2.txt | awk -v ram_req="$ram_req" -v clip_req="$clip_req" -F"\t" '($6 + $7 >= clip_req) && ($8 + $9 >= ram_req) { print $0}' > $dn_path/$familydname/p1_HighConfInsertions_calls.txt


#r script is in main folder
# this creates a start and end interval based on lclip and rclip position and adds a 40bp plus/minus margin
#if onesided add 40 plus minus the only site given, if s is then negative make 1

cd  $dn_path/$familydname/

Rscript $rfiles_path/denovo_breakpointmargins.R

sortBed -i $dn_path/$familydname/p1_HighConfInsertions_calls.bed > $dn_path/$familydname/p1_HighConfInsertions_calls_sort.bed

rm $dn_path/$familydname/p1_HighConfInsertions_calls.bed

#check if empty output then go to next family
if [[ -s $dn_path/$familydname/p1_HighConfInsertions_calls_sort.bed ]] ; then
:
else
continue
fi ;


#filter based on reference and knr knrref_path

bedtools intersect -wa -a $dn_path/$familydname/p1_HighConfInsertions_calls_sort.bed -b $knrref_path -v > $dn_path/$familydname/p1_HighConfInsertions_calls_sort_reffilt.bed

#check if empty output then go to next family
if [[ -s $dn_path/$familydname/p1_HighConfInsertions_calls_sort_reffilt.bed ]] ; then
:
else
echo -e "# Error No calls in $dn_path/$familydname/p1_HighConfInsertions_calls_sort_reffilt.bed "
continue
fi ;

#########################################
######## Discordant Reads Parents #######
#########################################

#Test mother and father  disc files to test if there are discordant reads in the insertion intervals

cat $calls_path/$bamFATHERsampleName.final/discordant_reads_tmp0  |  awk '{for(i=1;i<=($3+$4);i++)if ($2 > 50) print $1,$2-50,$2+50; else print $1,0,50;}'  OFS="\t" > $dn_path/$familydname/fa.disc.bed
cat $calls_path/$bamMOTHERsampleName.final/discordant_reads_tmp0 | awk '{for(i=1;i<=($3+$4);i++)if ($2 > 50) print $1,$2-50,$2+50; else print $1,0,50;}'  OFS="\t" > $dn_path/$familydname/mo.disc.bed


sortBed -i  $dn_path/$familydname/fa.disc.bed >  $dn_path/$familydname/fa.disc_sort.bed
sortBed -i $dn_path/$familydname/mo.disc.bed > $dn_path/$familydname/mo.disc_sort.bed


rm $dn_path/$familydname/fa.disc.bed
rm $dn_path/$familydname/mo.disc.bed


#Error if parent files are empty

if [[ -s $dn_path/$familydname/fa.disc_sort.bed  ]] ; then
:
else
echo -e "# Error $dn_path/$familydname/fa.disc_sort.bed  is empty"
continue
fi ;

if [[ -s $dn_path/$familydname/mo.disc_sort.bed ]] ; then
:
else
echo -e "# Error $dn_path/$familydname/mo.disc_sort.bed is empty"
continue
fi ;


#test if there are discordant reads in fa and mo intervals, only print intervals without disc reads here
#add regions where parental files have 1 disc reads or clipped reads
#use those with counts 0 or 1.

bedtools coverage -counts -a $dn_path/$familydname/p1_HighConfInsertions_calls_sort_reffilt.bed  -b $dn_path/$familydname/fa.disc_sort.bed $dn_path/$familydname/mo.disc_sort.bed > $dn_path/$familydname/p1_HighConfInsertions_sort_DISCparentIntersectCount.bed

awk -v parental_disc="$parental_disc"  -F"\t" '($4 <= parental_disc) { print $1"\t"$2"\t"$3 }' $dn_path/$familydname/p1_HighConfInsertions_sort_DISCparentIntersectCount.bed > $dn_path/$familydname/p1_HighConfInsertions_sort_noDISCparent.bed


#check if empty output then go to next family
if [[ -s $dn_path/$familydname/p1_HighConfInsertions_sort_noDISCparent.bed  ]] ; then
:
else
echo -e "# Ended $dn_path/$familydname/p1_HighConfInsertions_sort_noDISCparent.bed is empty"
continue
fi ;

#####################################
########### Clipped Reads Parents ###
#####################################

#now test if there are clipped reads in interval

#for reads aligned to TE
cat $calls_path/$bamFATHERsampleName.final/clip_reads_tmp0  |  cut -f 1,2,6,7 | awk '{for(i=1;i<=($3+$4);i++)if ($2 > 50) print $1,$2-50,$2+50; else print $1,0,50;}' OFS="\t" > $dn_path/$familydname/fa.softclip.bed
cat $calls_path/$bamMOTHERsampleName.final/clip_reads_tmp0 |  cut -f 1,2,6,7 | awk '{for(i=1;i<=($3+$4);i++)if ($2 > 50) print $1,$2-50,$2+50; else print $1,0,50;}' OFS="\t" > $dn_path/$familydname/mo.softclip.bed

#for all clipped reads

#added check for whether the coordinate is below
cat $calls_path/$bamFATHERsampleName.final/clip_reads_tmp0  |  cut -f 1,2,3,4 | awk '{for(i=1;i<=($3+$4);i++)if ($2 > 50) print $1,$2-50,$2+50; else print $1,0,50;}' OFS="\t" > $dn_path/$familydname/fa.softclipall.bed
cat $calls_path/$bamMOTHERsampleName.final/clip_reads_tmp0 |  cut -f 1,2,3,4| awk '{for(i=1;i<=($3+$4);i++)if ($2 > 50) print $1,$2-50,$2+50; else print $1,0,50;}' OFS="\t" > $dn_path/$familydname/mo.softclipall.bed

sortBed -i  $dn_path/$familydname/fa.softclip.bed >  $dn_path/$familydname/fa.softclip_sort.bed
sortBed -i $dn_path/$familydname/mo.softclip.bed > $dn_path/$familydname/mo.softclip_sort.bed

sortBed -i  $dn_path/$familydname/fa.softclipall.bed >  $dn_path/$familydname/fa.softclipall_sort.bed
sortBed -i $dn_path/$familydname/mo.softclipall.bed > $dn_path/$familydname/mo.softclipall_sort.bed

rm $dn_path/$familydname/fa.softclip.bed
rm $dn_path/$familydname/mo.softclip.bed

rm $dn_path/$familydname/fa.softclipall.bed
rm $dn_path/$familydname/mo.softclipall.bed

#Error if parent files are empty

if [[ -s $dn_path/$familydname/fa.softclip_sort.bed  ]] ; then
:
else
echo -e "# Error $dn_path/$familydname/fa.softclip_sort.bed  is empty"
continue
fi ;

if [[ -s $dn_path/$familydname/mo.softclip_sort.bed ]] ; then
:
else
echo -e "# Error $dn_path/$familydname/mo.softclip_sort.bed is empty"
continue
fi ;
if [[ -s $dn_path/$familydname/fa.softclipall_sort.bed  ]] ; then
:
else
echo -e "# Error $dn_path/$familydname/fa.softclipall_sort.bed  is empty"
continue
fi ;

if [[ -s $dn_path/$familydname/mo.softclipall_sort.bed ]] ; then
:
else
echo -e "# Error $dn_path/$familydname/mo.softclipall_sort.bed is empty"
continue
fi ;

#test if there are clipped reads that map to TEs in fa and mo intervals, only print intervals without clipped reads here

bedtools coverage -counts -a $dn_path/$familydname/p1_HighConfInsertions_sort_noDISCparent.bed  -b $dn_path/$familydname/fa.softclip_sort.bed $dn_path/$familydname/mo.softclip_sort.bed > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_TECLIPparentcounts.bed
bedtools coverage -counts -a $dn_path/$familydname/p1_HighConfInsertions_sort_noDISCparent.bed  -b $dn_path/$familydname/fa.softclipall_sort.bed $dn_path/$familydname/mo.softclipall_sort.bed > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_ALLCLIPparentcounts.bed


#keep those with the limit of te mapping clipped reads in parents
awk -v parental_TEclip="$parental_TEclip" -F"\t" '($4 <= parental_TEclip) { print $1"\t"$2"\t"$3}' $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_TECLIPparentcounts.bed  > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noTECLIPparent.bed

#keep those with the limit of all clipped reads in parents
awk -v parental_Allclip="$parental_Allclip" -F"\t" '($4 <= parental_Allclip) { print $1"\t"$2"\t"$3}' $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_ALLCLIPparentcounts.bed  > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noALLCLIPparent.bed


#end if no call passes
if [[ -s $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noTECLIPparent.bed ]] ; then
:
else
echo -e "# Ended $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noTECLIPparent.bed is empty"
continue
fi ;

if [[ -s $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noALLCLIPparent.bed ]] ; then
:
else
echo -e "# Ended $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noALLCLIPparent.bed is empty"
continue
fi ;

#merge these to keep only those that are found in both files meaning there are no clipped reads in parent or te mapping clipped reads in parents

Rscript $rfiles_path/xTea_parent_merged_clip.R


# Annotate reference and filter based on all reference and based on young reference

if [[ -s $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent.bed  ]] ; then
    #create folder to save files snapshots later in batch job
    #[[ -d $snapshot_path/$familydname ]] || mkdir $snapshot_path/$familydname
    #print in output file
:
else
echo -e "# Ended $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent.bed is empty"
continue
fi ;


bedtools intersect -wa -wb -loj -a $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent.bed  -b $ref_all_path | sort -k 1,1 -k2,2n | uniq  | sed 's/\./NA/g' > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_RefAnnot.bed


grep "NA" $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_RefAnnot.bed | cut -f 1,2,3 > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF.bed


############################# Print for igv manual inspection ###############################

if [[ -s $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF.bed  ]] ; then
    #create folder to save files snapshots later in batch job
    #[[ -d $snapshot_path/$familydname ]] || mkdir $snapshot_path/$familydname
    #print in output file
    echo -e "new\nsnapshotDirectory $snapshot_path/$familydname\ngenome hg38\nmaxPanelHeight 3000\nload $sessions_path/$familydname.xml\nsquish"
else
continue
fi ;


#print coordinates for batch job IGV
#add a wider range

while read chrom st end
do
    newst=$(expr $st - 400)
    newend=$(expr $end + 400)
    printf "%s\n" "goto $chrom:$newst-$newend"
    printf "%s\n" "squish"
    printf "%s\n" "snapshot"

done < $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF.bed
printf "%s\n"

############################################################################
# Obtain original calls that are de novo
while read chrom st end
do
    grep -w $chrom $dn_path/$familydname/p1_HighConfInsertions_calls_Original.txt  | grep -w $st | grep -w $end
done < $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent.bed > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_Original.txt

while read chrom st end
do
    grep -w $chrom $dn_path/$familydname/p1_HighConfInsertions_calls_Original.txt  | grep -w $st | grep -w $end
done < $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF.bed > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF_Original.txt


########################## Intersect Genes  ###############################

# Obtain gene annotation
bedtools intersect -wa -wb -a $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF.bed  -b $gene_annot_path | sort -k 1,1 -k2,2n | uniq  > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF_geneAnnot.bed

# all those that do intersect with reference
bedtools intersect -wa -wb -a $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent.bed  -b $gene_annot_path | sort -k 1,1 -k2,2n | uniq > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_geneAnnot.bed

# obtain just gene name for pli
bedtools intersect -wa -wb -a $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF.bed  -b $gene_annot_path | sed -e 's/_/\t/g' |  sort -k 1,1 -k2,2n | uniq | cut -f 1,2,3,7   > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF_geneAnnot_gene.bed
bedtools intersect -wa -wb -a $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent.bed  -b $gene_annot_path | sed -e 's/_/\t/g' |  sort -k 1,1 -k2,2n | uniq | cut -f 1,2,3,7   > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_geneAnnot_gene.bed

# intersect with pLI
while read chrom st end GENE
do
    grep -w $GENE $pLI_annot_path
    echo "$pLI_annot_path"
done < $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF_geneAnnot_gene.bed > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF_pLI.bed

while read chrom st end GENE
do
    grep -w $GENE $pLI_annot_path
    echo "$pLI_annot_path"
done < $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_geneAnnot_gene.bed > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_pLI.bed


# Obtain autism gene overlap
bedtools intersect -wa -wb  -a $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF.bed  -b $SFARI_annot_path | sort -k 1,1 -k2,2n | uniq > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF_SFARI.bed

# all those that do intersect with reference
bedtools intersect -wa -wb -a $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent.bed  -b $SFARI_annot_path | sort -k 1,1 -k2,2n | uniq > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_SFARI.bed

# intersect with scores:

while read chrom st end  chrom2 st2 end2 NAME
do
    grep -w $NAME $SFARI_scores_path
done < $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF_SFARI.bed > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_noREF_SFARI_scores.bed

while read chrom st end  chrom2 st2 end2 NAME
do
    grep -w $NAME $SFARI_scores_path
done < $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_SFARI.bed > $dn_path/$familydname/p1_HighConfInsertions_noDISCparent_noCLIPparent_SFARI_scores.bed

rm fa.disc_sort.bed
rm fa.softclipall_sort.bed
rm fa.softclip_sort.bed
rm mo.disc_sort.bed
rm mo.softclipall_sort.bed
rm mo.softclip_sort.bed

done < "$input"



