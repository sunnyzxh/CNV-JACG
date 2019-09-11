# CNV-JACG (v1.0)
CNV-JACG is a random forest based framework for assessing the accuracy of CNVs detected base on paired-end whole genome sequencing data

## Getting Started
CNV-JACG is written in perl and R

### Prerequisites

#### perl(v5.22.0) and following perl models
* Statistics::Basic (https://metacpan.org/pod/distribution/Statistics-Basic/lib/Statistics/Basic.pod)

- Bio::DB::HTS (https://github.com/Ensembl/Bio-DB-HTS)

* Bio::DB::HTS::Tabix (https://metacpan.org/pod/Bio::DB::HTS::Tabix)

If you encounter error like "perl: symbol lookup error: */perl5/lib/perl5/x86_64-linux-thread-multi//auto/Clone/Clone.so: undefined symbol: Perl_xs_handshake"
it means the last two packages were not successfully installed in your current used perl.

#### R(3.5.1) and the following R package
* randomForest (https://cran.r-project.org/web/packages/randomForest/randomForest.pdf)

If you encounter error "Error in library(randomForest) : there is no package called 'randomForest' Execution halted"
it indicates the failure of installing the package.

#### Other required softwares
* Bedtools(v2.25.0)

- Samtools(1.3.1)

***
### Installing
    
    wget https://github.com/sunnyzxh/CNV-JACG/archive/master.zip
    unzip master.zip
    cd CNV-JACG-master

#### Running the test
    cd bin
    chmod 777 samtools
    chmod 777 bedtools
    cd ../example
    sh test.sh
    cd test.result
    
If "ALL JOBS RUN SUCCESSFULLY!", there should be 5 files under test.result

Note: Since it is a demo, and due to the file size limitation of Github, we only provide a small partial of the following files. For you real data, you need to replace/supply the following files to the genome-wide. You could also contact sunnyzxh@connect.hku.hk to send you the whole files if you can't find the files youself.

* human reference genome (hg19)
- maf > 5% SNPs from 1000 Genomes Project
* Repeat masker, Segmental duplication coordinates

***

### Usage
#### Please prepare the following inputs

* The tab-separated file containg the coordinates of putative CNVs, "Chr\tStart\tEnd\tDEL/DUP\tother" ($preCNV)

- The bam file

Note: 
* Please **DO NOT** mix Deletion and Duplication within a file, it will lead to using deletion classifier to predict duplication, and vice versa. 
- The bam file should contain the **“SA”** tag, which could be generated by BWA “mem”

#### Running

perl CNV-JACG.pl

Function
    
    Assess the accuracy of putative CNVs from illumina pair-end WGS data

Usage
    
    perl CNV-JACG.pl -s example -p example.precnv -b example.bam -r ref.fa -o out

Options
   
    -h|-help            help
    
    -s|-sample   [s]    sample name
    
    -p|-precnv   [s]    putative CNVs file (chr,start,end,type(DEL|DUP),other. tab-separated)
    
    -b|-bam      [s]    bam file
    
    -r|-ref      [s]    reference genome
    
    -o|-outdir   [s]    output dir

#### Outputs

* $preCNV.het.prob

- $preCNV.repeat

* $preCNV.feature

- $preCNV.RFM.txt

* $preCNV.RFM.prediction (This is the final result, containing the prediction result in the first column, and other infomation)

***

### Support
Should you have any question, please feel free to contanct sunnyzxh@connect.hku.hk 
