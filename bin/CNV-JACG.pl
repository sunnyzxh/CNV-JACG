#!/usr/bin/perl -w
use strict;
use File::Spec;
use File::Basename;
use Statistics::Basic qw(:all);
use Getopt::Long;

my ( $Help, $sample, $precnv, $bam, $ref, $outdir );

GetOptions(
    'help|?'     => \$Help,
    'sample=s'   => \$sample,
    'precnv=s'   => \$precnv,
    'bam=s'      => \$bam,
    'ref=s'      => \$ref,
    'outdir=s'   => \$outdir
);

die `pod2text $0` if ( $Help or !$sample or !$ref or !$precnv or !$outdir or !$bam );

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Configuration
my $path_curf = File::Spec->rel2abs(__FILE__);
my ( $vol, $dirs, $file ) = File::Spec->splitpath($path_curf);
my $home      = dirname($dirs);
my $repeat    = "$home/lib/Repeat/repeat_chr";
my $snps      = "$home/lib/maf.sort.txt.gz";
my $fixregion = "$home/lib/fixed.region.noDGV.txt";
my $samtools  = "$home/bin/samtools";
my $bedtools  = "$home/bin/bedtools";
my $precnvbase = basename($precnv);

# Calculate het probability
&showLog("Calculating het probability...");
`/home/yerui/miniconda2/bin/perl $home/bin/Het-prob.pl -f $ref -m $snps -r $precnv -i $bam -o $outdir/$precnvbase.het.prob`;
my $hetprob = "$outdir/$precnvbase.het.prob";

# calculate depth in fixed region (for normalization)
&showLog("Calculating mean depth in fixed region...");
open IND, "$samtools bedcov $fixregion $bam|" or die $!;
my @all = ();
while (<IND>) {
    chomp;
    my @infod = split;
    next if ( $infod[-1] == 0 );
    my $depth = $infod[-1] / ( $infod[2] - $infod[1] + 1 );
    push @all, $depth;
}
my $normaldep = mean(@all);
&showLog("Done! Sample:$sample\tDepth:$normaldep");
close IND;

# get repeat percent for each CNV
&showLog("Calculating repeat percentage for each CNV...");
`$bedtools intersect -wao -a $precnv -b $repeat/SD_N_region_RepeatMasker_SimpleRepeat_merge.bed $repeat/N_region > $outdir/$precnvbase.repeat`;
open INR, "$outdir/$precnvbase.repeat";
my %repeats;
my %repeats_all;
while (<INR>) {
    chomp;
    my @info = split /\t/;
    next if ( $info[-3] == -1 );
    if ( $info[-5] == 1 ) {
        $repeats_all{"$info[0]\t$info[1]\t$info[2]"} += $info[-1];
    }
    elsif ( $info[-5] == 2 ) {
        $repeats{"$info[0]\t$info[1]\t$info[2]"} += $info[-1];
    }
    my $r = "$info[0]\t$info[1]\t$info[2]";
}
close INR;
&showLog("Done!");

# get depth and gc for each CNV
&showLog("Calculateing depth and GC for each CNV...");
open ING, "$samtools bedcov $precnv $bam|" or die $!;
my %depth;
my %gc;
while (<ING>) {
    chomp;
    my @info  = split;
    my $dna   = `$samtools faidx $ref $info[0]:$info[1]-$info[2]`;
    my $gcnum = ( $dna =~ s/[gcGC]//g );
    my $gcper = ( $gcnum / ( $info[2] - $info[1] + 1 ) ) * 100;
    $gcper = sprintf( "%.3f", $gcper );
    my $d = $info[-1] / ( $info[2] - $info[1] + 1 );
    $d = sprintf( "%.3f", $d );

    $depth{"$info[0]\t$info[1]\t$info[2]"} = $d;
    $gc{"$info[0]\t$info[1]\t$info[2]"}    = $gcper;
}
close ING;
&showLog("Done!");

&showLog("Getting other features...");

# Configure input and output files
open IN, "$precnv";
open OUT, ">$outdir/$precnvbase.feature";
print OUT
"CNV_Type\tSample\tChr\tStart\tRepeat_region?\tEnd\tRepeat_region?\tLength\tAverage_Depth($normaldep)\tFlag_Left\tFlag_Right\tFlag_SA_Reads\t#SA_Reads\tMicrohomo\tStart\tLeft_SA_Cluster\tEnd\tRight_SA_Cluster\tGC\t#Common_P_SNP\t#Het_P_SNP\tHet_Prob\tRepeat.Pct\tCalled_by\n";

open OUT1, ">$outdir/$precnvbase.RFM.txt";
print OUT1
"Chr\tStart\tEnd\tLeft.Repeat\tRight.Repeat\tRepeat.Pct\tMean.depth\tLeft.LL.RR\tLeft.RL\tRight.LL.RR\tRight.RL\tDI\tIN\tLeft.SCR\tRight.SCR\tLeft.Right.SCR\tLeft.SCR.cluster\tRight.SCR.cluster\tMicrohomo\tGC\tLength\tCommon_SNP\tHet_SNP\tHet_Prob\tCalled_by\n";
my $type = "Ini";

# Other features
while (<IN>) {
    chomp;
    my @info1 = split /\t/;
    next if ( $info1[2] - $info1[1] > 100000000 );
    if ( $type ne "Ini" and $type ne $info1[3] ) { die "ERROR!! Please Do NOT MIX DUP and DEL in a file!!\n"; }
    $type = $info1[3];

    # Check if the breakpoint located on repeat region: repeatmasker, simplerepeat ####
    my ( $l_repeat,      $r_repeat )      = ( 0, 0 );
    my ( $l_repeat_flag, $r_repeat_flag ) = ( 1, 1 );
    $l_repeat =
`less $repeat/$info1[0].repeat.txt.gz | awk '$info1[1] > \$2 && $info1[1] < \$3'`;
    $r_repeat =
`less $repeat/$info1[0].repeat.txt.gz | awk '$info1[2] > \$2 && $info1[2] < \$3'`;
    chomp $l_repeat;
    chomp $r_repeat;
    $l_repeat = join ",", ( split /\t/, $l_repeat );
    $r_repeat = join ",", ( split /\t/, $r_repeat );
    if ( !$l_repeat ) { $l_repeat = "no"; $l_repeat_flag = 0; }
    if ( !$r_repeat ) { $r_repeat = "no"; $r_repeat_flag = 0; }

    # Get repeat percentage
    my $region = "$info1[0]\t$info1[1]\t$info1[2]";
    $repeats_all{$region} = 0 if ( !$repeats_all{$region} );
    my $rp = $repeats_all{$region} / ( $info1[2] - $info1[1] + 1 );
    $rp = sprintf( "%.3f", $rp );

    # Get CG content
    my $gcper = $gc{$region};

    # Get the het probability
    my $snp =
`less $hetprob | awk '\$1 == "$info1[0]" && \$2 == $info1[1] && \$3 == $info1[2]'`;
    chomp $snp;
    my ( $snp_num, $het_num, $het_prob ) = ( split /\t/, $snp )[ -3, -2, -1 ];
    $het_prob = 0 if ( $het_prob eq "-" );

    # Orientation of reads on the 150bp upstream/downstream of left/right side
    my $s_l        = $info1[1] - 150;
    my $e_l        = $info1[1] + 150;
    my $s_r        = $info1[2] - 150;
    my $e_r        = $info1[2] + 150;
    my $flag_left  = orientation_flag( $info1[0], $s_l, $e_l );
    my $flag_right = orientation_flag( $info1[0], $s_r, $e_r );

    # Cluster SA (chimeric) reads around breakpoint (Exclude Secondary and PCR-induced Dup reads)
    open IN1,
"$samtools view -F 1280 $bam $info1[0]:$s_l-$e_l $info1[0]:$s_r-$e_r | grep SA: |"
      or die $!;

    my $sr        = 0;
    my %sa_flag   = ( "DI" => 0, "IN" => 0 );
    my @micro     = ();
    my @left_b    = ();
    my @right_b   = ();
    my @left_pos  = ();
    my @right_pos = ();
    my $comp      = 0;

    while (<IN1>) {
        chomp;
        my @info = split;
        next
          if ( $info[5] !~ /[S]/
            or $info[5] !~ /[M]/
            or ( $info[5] =~ tr/[SM]/[SM]/ ) != 2 );
        my ( $cigar2, $chr, $pos, $strand ) =
          ( split /[:,]/, ( ( split /;/, $info[11] )[0] ) )[ 5, 2, 3, 4 ];
        next
          if ( $cigar2 !~ /[S]/
            or $cigar2 !~ /[M]/
            or ( $cigar2 =~ tr/[SM]/[SM]/ ) != 2 );
        my $cigar1_s = ( $info[5] =~ /(\d+)S/ )[0];
        my $cigar1_m = ( $info[5] =~ /(\d+)M/ )[0];
        my $cigar2_s = ( $cigar2  =~ /(\d+)S/ )[0];
        my $cigar2_m = ( $cigar2  =~ /(\d+)M/ )[0];

        if ( abs( $cigar1_s - $cigar2_m ) == abs( $cigar1_m - $cigar2_s ) ) {
            $sr++;
            $comp = 1;
        }
        if (    ( $info[2] eq $info1[0] )
            and ( abs( $info[3] - $info1[1] ) <= 100 )
            and ( $comp == 1 ) )
        {
            push @left_b,   abs( $cigar1_s - $cigar2_m );
            push @left_pos, "$chr,$pos";
            push @micro,    abs( $cigar1_s - $cigar2_m );
        }
        elsif ( ( $info[2] eq $info1[0] )
            and ( abs( $info[3] - $info1[2] ) <= 100 )
            and ( $comp == 1 ) )
        {
            push @right_b,   abs( $cigar1_s - $cigar2_m );
            push @right_pos, "$chr,$pos";
            push @micro,     abs( $cigar1_s - $cigar2_m );
        }

        # Orientation of chimeric reads around breakpoints (all of them span breakpoint ?)
        $sa_flag{"DI"}++
          if ( ( ( $info[1] & 0x10 ) == 0 and $strand eq "+" )
            || ( ( $info[1] & 0x10 ) == 16 and $strand eq "-" ) );
        $sa_flag{"IN"}++
          if ( ( ( $info[1] & 0x10 ) == 0 and $strand eq "-" )
            || ( ( $info[1] & 0x10 ) == 16 and $strand eq "+" ) );
    }
    close IN1;

    my ( $left_t, $right_t, $micros ) = ( 0, 0, 0 );
    $left_t  = pos_cluster(@left_pos)  if ( $#left_pos >= 0 );
    $right_t = pos_cluster(@right_pos) if ( $#right_pos >= 0 );
    $micros  = micro_cluster(@micro)   if ( $#micro >= 0 );
    my $len = $info1[2] - $info1[1] + 1;

    # Get mean coverage of this CNV
    my $av_dep = $depth{$region};

    # Get detail info.
    my @n_micros   = split /;/,     $micros;
    my @n_left_t   = split /;/,     $left_t;
    my @n_left_ts  = split /[,;:]/, $left_t;
    my @n_right_t  = split /;/,     $right_t;
    my @n_right_ts = split /[,;:]/, $right_t;
    my @flag_lefts  = ( split /[:;]/, $flag_left );
    my @flag_rights = ( split /[:;]/, $flag_right );

    # Normalization
    $av_dep        = $av_dep / $normaldep;
    $av_dep        = sprintf( "%.3f", $av_dep );
    $sa_flag{"DI"} = $sa_flag{"DI"} / $normaldep;
    $sa_flag{"DI"} = sprintf( "%.3f", $sa_flag{"DI"} );
    $sa_flag{"IN"} = $sa_flag{"IN"} / $normaldep;
    $sa_flag{"IN"} = sprintf( "%.3f", $sa_flag{"IN"} );
    $sr            = $sr / $normaldep;
    $sr = sprintf( "%.3f", $sr );

    # Output
    my $diflag = $sa_flag{"DI"};
    my $inflag = $sa_flag{"IN"};
    print OUT
"$type\t$sample\t$info1[0]\t$info1[1]\t$l_repeat\t$info1[2]\t$r_repeat\t$len\t$av_dep\t$flag_left\t$flag_right\tDI:$diflag;IN:$inflag\t$sr\t$micros\t$info1[1]\t$left_t\t$info1[2]\t$right_t\t$gcper\t$snp_num\t$het_num\t$het_prob\t$rp\t$info1[-1]\n";

    my ( $left_scr, $right_scr, $left_scr_cluster, $right_scr_cluster,
        $micros_flag );
    if ( $left_t eq "0" ) { $left_scr = 0; $left_scr_cluster = 0; }
    else                  { $left_scr = $n_left_ts[2]; $left_scr_cluster = 1; }
    if ( $right_t eq "0" ) { $right_scr = 0; $right_scr_cluster = 0; }
    else { $right_scr = $n_right_ts[2]; $right_scr_cluster = 1; }
    if   ( $micros eq "0" ) { $micros_flag = 0; }
    else                    { $micros_flag = ( split /[:;]/, $micros )[1]; }

    print OUT1
"$info1[0]\t$info1[1]\t$info1[2]\t$l_repeat_flag\t$r_repeat_flag\t$rp\t$av_dep\t$flag_lefts[-3]\t$flag_lefts[-1]\t$flag_rights[-3]\t$flag_rights[-1]\t$diflag\t$inflag\t$left_scr\t$right_scr\t$sr\t$left_scr_cluster\t$right_scr_cluster\t$micros_flag\t$gcper\t$len\t$snp_num\t$het_num\t$het_prob\t$info1[-1]\n";

}
&showLog("Done!");
&showLog("Prediction...");
if ( $type eq "DEL" ) {
    `Rscript $home/bin/Prediction.R $home/lib/trained_model_del.rds $outdir/$precnvbase.RFM.txt $outdir/$precnvbase.RFM.prediction`;
}
elsif ( $type eq "DUP" ) {
    `Rscript $home/bin/Prediction.R $home/lib/trained_model_dup.rds $outdir/$precnvbase.RFM.txt $outdir/$precnvbase.RFM.prediction`;
}
&showLog("Done!\nALL JOBS RUN SUCCESSFULLY!");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub showLog {
    my @t = localtime();
    printf STDERR "[%04d-%02d-%02d %02d:%02d:%02d]\t%s\n", $t[5] + 1900,
      $t[4] + 1, @t[ 3, 2, 1, 0 ], $_[0];
}

sub pos_cluster {
    my (@info) = @_;
    my %hash;
    $hash{ $info[0] } = 1;
    for ( my $i = 1 ; $i <= $#info ; $i++ ) {
        my ( $chr, $pos ) = ( split /,/, $info[$i] );
        my $flag = 0;
        for my $k ( sort keys %hash ) {
            my ( $chr_o, $pos_o ) = ( split /,/, $k )[ 0, 1 ];
            if ( $chr eq $chr_o and ( abs( $pos - $pos_o ) <= 200 ) ) {
                $hash{$k}++;
                $flag = 1;
            }
        }
        $hash{ $info[$i] } = 1 if ( $flag == 0 );
    }
    my @merge = ();
    my $count = 0;
    for my $k ( sort { $hash{$b} <=> $hash{$a} } keys %hash ) {
        $count++;
        if ( $count <= 3 ) {
            $hash{$k} = $hash{$k} / $normaldep;
            $hash{$k} = sprintf( "%.3f", $hash{$k} );
            push @merge, "$k:$hash{$k}";
        }
    }
    my $merges = join ";", @merge;
    $merges = 0 if ( !$merges );
    return ($merges);
}

sub micro_cluster {
    my (@info) = @_;
    my %hash;
    for my $i ( 0 .. $#info ) {
        $hash{ $info[$i] }++;
    }
    my @merge = ();
    my $count = 0;
    for my $k ( sort { $hash{$b} <=> $hash{$a} } keys %hash ) {
        $count++;
        if ( $count <= 3 ) {
            $hash{$k} = $hash{$k} / $normaldep;
            $hash{$k} = sprintf( "%.3f", $hash{$k} );
            push @merge, "$k:$hash{$k}";
        }
    }
    my $merges = join ";", @merge;
    $merges = 0 if ( !$merges );
    return ($merges);
}

sub orientation_flag {
    my (@info1) = @_;
    open IN2, "$samtools view -F 1280 $bam $info1[0]:$info1[1]-$info1[2] |"
      or die $!;
    my %flag        = ( "LR" => 0, "LL,RR" => 0, "RL" => 0 );
    my $count       = 0;
    my $revert_flag = 0;
    while (<IN2>) {
        chomp;
        my @info2 = split;
        my $line  = $_;
        if ( ( $info2[1] & 0x10 ) == 0 and ( $info2[1] & 0x20 ) == 32 ) {
            if ( $info2[3] < $info2[7] ) {
                $flag{"LR"}++;
            }
            else {
                $flag{"RL"}++;
            }
        }
        elsif ( ( $info2[1] & 0x10 ) == 16 and ( $info2[1] & 0x20 ) == 0 ) {
            if ( $info2[3] > $info2[7] ) {
                $flag{"LR"}++;
            }
            else {
                $flag{"RL"}++;
            }
        }
        elsif (( ( $info2[1] & 0x10 ) == 0 and ( $info2[1] & 0x20 ) == 0 )
            or ( ( $info2[1] & 0x10 ) == 16 and ( $info2[1] & 0x20 ) == 32 ) )
        {
            $flag{"LL,RR"}++;
        }
        else {
            print "$line\n";
        }
        $count++;
    }
    close IN2;
    $count      = $count / $normaldep;
    $count      = sprintf( "%.3f", $count );
    $flag{"LR"} = $flag{"LR"} / $normaldep;
    $flag{"LR"}    = sprintf( "%.3f", $flag{"LR"} );
    $flag{"RL"}    = $flag{"RL"} / $normaldep;
    $flag{"RL"}    = sprintf( "%.3f", $flag{"RL"} );
    $flag{"LL,RR"} = $flag{"LL,RR"} / $normaldep;
    $flag{"LL,RR"} = sprintf( "%.3f", $flag{"LL,RR"} );
    my $flag_join =
      "$count;LR:$flag{\"LR\"};LL,RR:$flag{\"LL,RR\"};RL:$flag{\"RL\"}";
    return ($flag_join);
}

=head1 Function
    Assess the accuracy of putative CNVs from illumina pair-end WGS data

=head1 Usage
    perl CNV-JACG.pl -s example -p example.precnv -b example.bam -r ref.fa -o out

=head1 Options
    -h|-help            help
    -s|-sample   [s]    sample name
    -p|-precnv   [s]    putative CNVs file (chr,start,end,type(DEL|DUP),other. tab-separated)
    -b|-bam      [s]    bam file
    -r|-ref      [s]    reference genome
    -o|-outdir   [s]    output dir

=head1 Author
    sunny@connect.hku.hk

=head1 Version
    v0.1; 2019-03-12

=cut
