#!/bin/bash

#######################################################################################################
# A pipeline mapping sequences using bwa, samtools and gatk ReAlignerTargetCreator and indelReAligner.  
# Based on a script originally provided by Carla Coll Costa.
# Run this script by providing the population and individual/sample name.
#######################################################################################################

#By pair of fastq (i.e. read groups) 
#Step1: MAPPING EACH READ PAIR -> bwa mem
#Step2a: FILL IN MATE COORDINATES AND INSERT SIZE FIELDS -> samtools fixmate -m
#Step2b: SORT READS BY COORDINATES -> samtools sort (by coordinates)
#Step2c: INDEX SAM/BAM FILES -> samtools index 

#By sample
#Step3: MERGE MULTIPLE FILES IN ONE OUTPUT -> samtools merge
#Step4a: LIST OF REGIONS TO REALIGN -> gatk ReAlignerTargetCreator 
#Step4b: REALIGN AROUND INDELS -> gatk indelReAligner 
#Step5: samtools markdup

debug=0

time=`date +%F-%T`
reference=/path_to_ref_genome/Gasterosteus_aculeatus.BROADS1.dna.toplevel.fa # path to indexed reference genome
path=/path_to_folder_with_fasta_files/

pop=$1 # population given as an argument on command line
sample=$2 # individual given as an argument on command line
BAMMAP=/BAM/$pop/bamFileList # tab-separated file to store paths of the output bamfiles, providing path to folder where output BAM files will be saved

outdir=/BAM/$pop/$sample # output directory
tmpdir=/tmp/$pop/$sample # temporal directory
sortdir=/tmp/sort/ # directory for sorting in temporal directory

ls $path/$pop/$sample/*1.f*.gz 2>/dev/null


FILES1=( `ls $path/$pop/$sample/*1.f*.gz 2>/dev/null`" " )  # lists fastq files (mate) sequenced at BGI and MACROGEN
FILES2=( `ls $path/$pop/$sample/*2.f*.gz 2>/dev/null`" " ) # lists fastq files (pair) sequenced at BGI and MACROGEN
echo $FILES1
echo $FILES2


mkdir -p $outdir
mkdir -p $tmpdir
	
cd $tmpdir
	
	
for ((i=0;i<${#FILES1[@]};++i)); do ## for loop over the read groups (each mate and pair is one readgroup)
    ID=`zcat ${FILES1[i]} |head -1|cut -f1 -d " "`
    
	echo $ID
    if true; then			
	
	logfile=$tmpdir/$sample""_analysis.log.$time
	
	out1=$tmpdir/$sample""$ID""_aligned.bam # bwa mem out
	out2a=$tmpdir/$sample""$ID""_fixmate.bam # fixmate out 
	out2b=$tmpdir/$sample""$ID""_sorted.bam # sorted bam out
	out2c=$tmpdir/$sample""$ID""_sorted.bai # indexed bam out
	out3=$tmpdir/$sample""_merged.bam # bam after merging
	out3b=$tmpdir/$sample""_merged.bai # index for merged bam
	out4a=$tmpdir/$sample""_target_intervals.list # targets for realignment
	out4b=$tmpdir/$sample""_realigned.bam # realigned bam
	out4c=$tmpdir/$sample""_realigned.bai # index for realigned bam
	out5=$tmpdir/$sample""_markdup.bam # mark duplicates - don't remove them
	out5a=$tmpdir/$sample""_markdup.bai # index duplicate marked bam

	echo $out1
	echo $out2a
	echo $out2b
	echo $out2c 
	# Step1: bwa
	
	date >> $logfile
	echo -e "\n3s: bwa\nRead group: "$ID >> $logfile
	echo -e "${FILES1[i]}\n${FILES2[i]}" >> $logfile
	
	if [ ! -e $out1 ] && [ ! -e $out2a ] && [ ! -e $out2b ] && [ ! -e $out2c ] && [ ! -e $out3 ]; then
	
	  cat ${FILES1[i]} > file1.fq.gz
	  cat ${FILES2[i]} > file2.fq.gz
	  
	  ( bwa mem -t 10 -M -R "@RG\tID:$ID\tSM:$sample\tPL:ILLUMINA\tLB:LIB-1" \
		$reference \
		file1.fq.gz \
		file2.fq.gz  \
		| samtools view -h -b -o $out1 - ) 2>> $logfile
		
	fi
	
	# Step2a: samtools fixmate
	echo >> $logfile
	date >> $logfile
	echo -e "\n3s: samtools fixmate\nRead group: "$ID >> $logfile
	
	if [ -e $out1 ] && [ ! -e $out2a ] && [ ! -e $out2b ] && [ ! -e $out2c ] && [ ! -e $out3 ]; then
	
	  samtools fixmate -m $out1 $out2a &>> $logfile  
	  if [ -e $out2a ] && [ $debug -eq 0 ]; then 
		  rm -f $out1 file1.fq.gz file2.fq.gz
	  fi
	fi
	
	
	# Step2b: samtools sort
	echo >> $logfile
	date >> $logfile
	echo -e "\n3s: samtools sort\nRead group: "$ID >> $logfile
		
	if [ -e $out2a ] && [ ! -e $out2b ] && [ ! -e $out2c ] && [ ! -e $out3 ]; then
	  
	  samtools sort -@ 10 -T $sortdir/$sample""_$ID -O bam -o $out2b $out2a &>> $logfile
	  if [ -e $out2b ] && [ $debug -eq 0 ]; then 
		  rm -f $out2a 
	  fi
	fi
	
	# Step2c: samtools index
	echo >> $logfile
	date >> $logfile
	echo -e "\n3s: samtools index\nRead group: "$ID >> $logfile
	
	if [ -e $out2b ] && [ ! -e $out2c ] && [ ! -e $out3 ]; then
	  samtools index  $out2b $out2c &>> $logfile
	fi 
    
    fi ## End of the debug statement
done ## end of the for loop over read groups

if true; then

debug=0

#By sample
#Step3: samtools merge
#Step4a: gatk ReAlignerTargetCreator
#Step4b:gatk indelReAligner
#Step5:samtools markdup
#Step5a: samtools index


# Step3 samtools merge
echo >> $logfile
date >> $logfile
echo -e "\n3s: samtools merge\n" >> $logfile
if [ ! -e $out3 ] && [ ! -e $out3b ] && [ ! -e $out4a ] && [ ! -e $out4b ]  ; then
  samtools merge $out3 $tmpdir/*_sorted.bam
  if [ -e $out3 ] && [ $debug -eq 0 ]; then
    rm -f $tmpdir/$sample""*_sorted.ba*
  fi
fi

# Step 3b samtools index
echo >> $logfile
date >> $logfile
echo -e "\n3s: samtools index\n" >> $logfile
if [ -e $out3 ] && [ ! -e $out3b ] && [ ! -e $out4a ] && [ ! -e $out4b ] ; then
    samtools index $out3 $out3b &>> $logfile
fi

# Step4a  gatk ReAlignerTargetCreator 
echo >> $logfile
date >> $logfile
echo -e "\n3s: gatk RealignerTargetCreator\n" >> $logfile

if [ -e $out3 ] && [ -e $out3b ] && [ ! -e $out4a ] && [ ! -e $out4b ]  ; then

  java -Xmx16G -XX:MaxHeapSize=2560m -jar /path_to_GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	-R $reference \
	-I $out3 \
	-o $out4a \
	-nt 10 \
	&>> $logfile
	  
fi

# Step4: gatk indelReAligner
echo >> $logfile
date >> $logfile
echo -e "\n3s: gatk indelReAligner\n" >> $logfile

if [ -e $out3 ] && [ -e $out3b ] && [ -e $out4a ] && [ ! -e $out4b ] && [ ! -e $out4c ]; then

  java -Xmx16G -XX:MaxHeapSize=2560m -jar /path_to_GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R $reference \
	-I $out3 \
	-targetIntervals $out4a \
	-o $out4b \
	&>> $logfile

fi

#Step5:samtools markdup
 #DUPLICATES ARE NOT REMOVED, THEY ARE ONLY MARKED!
echo >> $logfile
date >> $logfile
echo -e "\n3s: samtools markdup\n" >> $logfile

if [ -e $out4b ] && [ -e $out4c ] && [ ! -e $out5 ] && [ ! -e $out5a ]; then

  samtools markdup $out4b $out5

fi
#Step5a: samtools index
echo >> $logfile
date >> $logfile
echo -e "\n3s: samtools index\n" >> $logfile

if [ -e $out5 ] && [ ! -e $out5a ]; then

    samtools index $out5 $out5a &>> $logfile
fi

##	below lines to copy the output and remove files.
if [ -e $out4b ] && [ -e $out4c ] && [ -e $out5 ] && [ -e $out5a ]; then
  rm -f $out4b $out4c
  echo -e "$sample\t$outdir/$sample""_markdup.bam" >> $BAMMAP
  mv $out5 $outdir
  mv $out5a $outdir
fi
echo >> $logfile
date >> $logfile
echo -e "\n3s: Pipeline finished\n" >> $logfile

#read logfile
#mv "$logfile" "$outdir"

fi 	

#Output: 1 bam and 1 bai file per sample!
