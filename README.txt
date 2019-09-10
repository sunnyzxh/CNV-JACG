### Before Downloading CNV-JACG! ###
PLEASE make sure you have successfully installed the following Perl and R Packages!

Perl:
1) Statistics::Basic (https://metacpan.org/pod/distribution/Statistics-Basic/lib/Statistics/Basic.pod)
2) Bio::DB::HTS (https://github.com/Ensembl/Bio-DB-HTS)
3) Bio::DB::HTS::Tabix (https://metacpan.org/pod/Bio::DB::HTS::Tabix)
If you encounter error like "perl: symbol lookup error: **/perl5/lib/perl5/x86_64-linux-thread-multi//auto/Clone/Clone.so: undefined symbol: Perl_xs_handshake", it means the last two packages were not successfully installed in your current used perl.

R:
randomForest (https://cran.r-project.org/web/packages/randomForest/randomForest.pdf)
If you encounter error "Error in library(randomForest) : there is no package called 'randomForest' Execution halted", it indicates the failure of installing the package.

### Note! ###
DO NOT mix Deletion and Duplication in a file, it will lead to using deletion classifier to predict duplication, and vice versa.

### Installation ###
Download CNV-JACG.tar.gz from .....
$tar -xvzf CNV-JACG.tar.gz
$cd CNV-JACG
$cd example
$cat test.sh
change the -ref parameter by using the reference genome of your direction
run the changed script


PLEASE feel free to contact sunnyzxh@connect.hku.hk should you have any questions!

To do list:
1) Add Merge program
2) Population level
