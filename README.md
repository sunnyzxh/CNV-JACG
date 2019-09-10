# CNV-JACG
CNV-JACG is a random forest based framework for assessing the accuracy of CNVs detected base on paired-end whole genome sequencing data

## Getting Started
CNV-JACG is written in perl and R

### Prerequisites

#### perl(v5.22.0) and following perl models
Statistics::Basic (https://metacpan.org/pod/distribution/Statistics-Basic/lib/Statistics/Basic.pod)

Bio::DB::HTS (https://github.com/Ensembl/Bio-DB-HTS)

Bio::DB::HTS::Tabix (https://metacpan.org/pod/Bio::DB::HTS::Tabix)

If you encounter error like "perl: symbol lookup error: */perl5/lib/perl5/x86_64-linux-thread-multi//auto/Clone/Clone.so: undefined symbol: Perl_xs_handshake"
it means the last two packages were not successfully installed in your current used perl.

#### R(3.5.1) and the following R package
randomForest (https://cran.r-project.org/web/packages/randomForest/randomForest.pdf)

If you encounter error "Error in library(randomForest) : there is no package called 'randomForest' Execution halted"
it indicates the failure of installing the package.

#### Other required software
Bedtools(v2.25.0)

Samtools(1.3.1)

### Installing
Download CNV-JACG.tar.gz from 

### Usage
Please prepare

1, The tab-separated file containg the coordinates of putative CNVs, "Chr\tStart\tEnd\tDEL/DUP\tother" ($inputCNV)

2, The vcf file (optional) ($vcf)

3, The bam file ($bam)

### Step1: Get the probability of the heterozygosity of putative CNV regions given the common SNPs from the 1000 Genomes Project
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

### Step2: Get the overlapping of putative CNV regions with repeat region
bedtools intersect -wao -a $inputCNV -b Repeat/SD_N_region_RepeatMasker_SimpleRepeat_merge.bed Repeat/N_region > $inputCNV.repeat

### Step3: Get the depth and GC content of putative CNV regions
perl get_depth_gc_samtools.pl $bam $inputCNV

The output of this step is $inputCNV.depth.gc

### Step4: Get other features and assessment results
perl CNV-JACG.clear.pl

Options

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
Should you have any question, please feel free to contanct sunnyzxh@connect.hku.hk 
