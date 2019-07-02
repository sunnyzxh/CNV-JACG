# CNV-JACG
CNV-JACG is a random forest based framework for assessing the accuracy of CNVs detected based on paired-end whole genome sequencing data

## Installation
CNV-JACG is written in perl and R

just make sure you have already install perl(v5.22.0) and R(3.5.1)

### Required perl model
Statistics::Basic (https://metacpan.org/pod/distribution/Statistics-Basic/lib/Statistics/Basic.pod)
### Required R packages
randomForest (https://cran.r-project.org/web/packages/randomForest/randomForest.pdf)
### Other required program
Bedtools(v2.25.0)

Samtools(1.3.1)

## Usage
Please prepare

1,tab-separated file containg the coordinates of putative CNVs, "Chr\tStart\tEnd\tDEL/DUP\tother" ($inputCNV)

2,the vcf file (optional) ($vcf)

3,the bam file ($bam)

### Step1: Get the probability of the hetrozygosity of putative CNVs given the common SNPs from 1000 Genomes Project
perl geno.pl

Function
    
    Calculate het SNPs probability in given genome regions

Usage
    
    perl geno.pl -f ref.fa -m maf_f.gz -r region_f -i in.bam [Options] -o out

Options
   
    -h|-help                help
    
    -i|-inbam         [s]   inbam file
    
    -f|-fasta         [s]   reference.fa
    
    -m|-mafFile       [s]   bgziped maf file (chr, pos, ref, alt, maf)
    
    -r|-regionFile    [s]   region file (chr, start, end)
    
    -b|-baseQuaCutoff [i]   base quality cutoff [20]
    
    -d|-debug               debug mode
    
    -s|-seed          [i]   seed for srand [time]
    
    -o|-outfile       [s]   outfile

### Step2: Get the overlapping with repeat region
bedtools intersect -wao -a $inputCNV -b Repeat/SD_N_region_RepeatMasker_SimpleRepeat_merge.bed Repeat/N_region > $inputCNV.repeat

### Step3: Get the depth and GC of putative CNV regions
perl get_depth_gc_samtools.pl $bam $inputCNV

the output of this step is $inputCNV.depth.gc

### Step4: Get other features and assessment result
perl CNV-JACG.clear.pl

-sample <String>    The sample name;

-precnv <String>    The coordinates of to-be-assessed CNVs "Chr\tStart\tEnd\tDEL/DUP\tCalled_by";

-vcf    <String>    The vcf file;

-outdir <String>    The output dir;

-bam    <String>    The bam file;

the output of this step includes

1, $inputCNV.judge.txt

2, $inputCNV.RFM.tmp

3, $inputCNV.RFM.txt

4, $inputCNV.final.judge.txt

### Note
1, The bam file should contain the “SA” tag, which could be realized by BWA “mem”

2, The prefix of the the file $inputCNV should keep the same in all four steps.

# Support
Please feel free to contanct sunnyzxh@connect.hku.hk should you have any question

