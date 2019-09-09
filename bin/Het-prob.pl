#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;
use Bio::DB::HTS;
use Bio::DB::HTS::Tabix;

my ( $Help, $Ref_f, $MAF_f, $Region_f, $debug, $Seed, $inbam, $outfile );

my ( $BaseQuaCutoff, )
=  ( 20,             );

GetOptions(
    'help|?'          => \$Help,
    'fasta=s'         => \$Ref_f,
    'mafFile=s'       => \$MAF_f,
    'regionFile=s'    => \$Region_f,
    'baseQuaCutoff=i' => \$BaseQuaCutoff,
    'seed=i'     => \$Seed,
    'debug'      => \$debug,
    'inbam=s'    => \$inbam,
    'outfile=s'  => \$outfile,
);
die `pod2text $0` if $Help or !$Ref_f or !$MAF_f or !$Region_f or !$outfile or !$inbam;

defined $Seed ? srand($Seed) : srand(time);

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# open bam file and ref
my $hts = Bio::DB::HTS->new(
    -bam  =>"$inbam",
    -fasta=>"$Ref_f",
);

# get prior probability from maf file
my %prior;
my %alt_at; # alt allele at a pos
my $tabix = Bio::DB::HTS::Tabix->new(filename => $MAF_f);

# posterior probability of SNPs
my @P_het;
my ($chr, $start, $end);

# callback func of pileup
my $callback = sub {

    # $pileup_r is a reference to an array of Bio::DB::HTS::Pileup objects
    my ($seqid, $pos, $pileup_r) = @_;

    # return if not have maf info
    return if !defined $prior{$seqid}{$pos};

    # ref base from fasta
    my $refbase = $hts->segment($seqid,$pos,$pos)->dna;
    $refbase =~ tr/acgtn/ACGTN/;

    # hom or het at each position
    my @baseV;
    my @quaV;

    for my $p (@$pileup_r) {
        my $aln     = $p->alignment;
        next if $p->indel or $p->is_refskip;      # don't deal with these ;-)

        # skip head & tail
        next if ( $p->qpos < 5 or length($aln->qseq) - $p->qpos + 1 < 5 );

        my $qbase  = substr($aln->qseq,$p->qpos,1);
        next if $qbase eq "N";

        my $qscore = $aln->qscore->[$p->qpos];
        next if $qscore < $BaseQuaCutoff;

        push @baseV, $qbase;
        push @quaV,  $qscore;
    }

    return if scalar(@baseV) == 0;
    push @P_het, posterior_het(\@baseV, \@quaV, $refbase, $alt_at{$seqid}{$pos}, $prior{$seqid}{$pos});
};

open my $OUT, ">$outfile" or die $!;
open my $RE, $Region_f or die $!;
while (<$RE>) {
    chomp;
    ($chr, $start, $end) = split;
    print "original: ($chr, $start, $end)\n" if $debug;

    %prior = ();
    %alt_at = ();
    @P_het = ();
    my @dots = ();

    # fetch maf of candidate snps
    my $iter = $tabix->query("$chr:$start-$end");
    if ( !defined $iter ) {
        printf $OUT "$_\t0\t0\t-\n";
        next;
    }

    while ( my $line = $iter->next ) {
        push @dots, $line;
    }

    my %to_call = ();
    my $snpN = scalar @dots;
    print "snpN: $snpN\n" if $debug;

    # random select 100 snps to call
    if ( $snpN > 100 ) {
        my $snpN_cnt = 0;
        while (1) {
            my $i = int( rand($snpN) );
            if ( !defined $to_call{ $dots[$i] } ) {
                $to_call{ $dots[$i] }++;
                $snpN_cnt++;
            }

            last if $snpN_cnt == 100;
        }
    }
    elsif ( $snpN > 0 ) {
        map { $to_call{$_}++ } @dots;
    }
    else {
        printf $OUT "$_\t0\t0\t-\n";
        next;
    }

    print "to call: ".scalar(keys %to_call)."\n" if $debug;

    # calculate prior probability for pileup part
    for my $line ( keys %to_call ) {
        chomp;
        my ($c, $p, $r, $a, $maf) = split /\s+/, $line;

        $prior{$c}{$p}{"$r$r"} = (1-$maf)**2;
        $prior{$c}{$p}{"$r$a"} = 2 * (1-$maf) * $maf;
        $prior{$c}{$p}{"$a$a"} = $maf**2;
        $alt_at{$c}{$p} = $a;

        $hts->fast_pileup("$chr:$p-$p", $callback);
        %prior = ();
        %alt_at = ();
    }

    printf $OUT "$_\t%d\t%d\t", scalar(@P_het), hetCount(\@P_het);
    my $mean_P_het = _mean(\@P_het);
    defined $mean_P_het ? printf $OUT "%.2f\n", $mean_P_het : print $OUT "-\n";
}
close $RE;
$tabix->close;
close $OUT;

&showLog("Done!");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub hetCount {
    my ($r, ) = @_;

    my $cnt = 0;
    print "P_het: @$r\n" if $debug;
    map { $cnt++ if $_ > 0.5 } @$r;

    return $cnt;
}

sub _mean {
    my ($r, ) = @_;
    return if scalar(@$r) == 0;

    my $sum = 0;
    map { $sum += $_ } @$r;
    return $sum/scalar(@$r);
}

# ->@baseV, ->@quaV, $refbase, $altbase ->%prior
# return het SNP's posterior probablity
#
sub posterior_het {
    my ($bs, $q, $ref, $alt, $pr) = @_;

    my %cnt = ();
    map { $cnt{$_}++ } @$bs;

    if ($debug) {
        map { print "$cnt{$_}$_ " } sort keys %cnt;
        print " - $ref$alt\n";
    }

    my $depth = scalar @$bs;

    # hom point
    return 0 if !defined $cnt{$ref} || !defined $cnt{$alt};
    return 0 if $cnt{$alt}/$depth < 0.01 || $cnt{$ref}/$depth < 0.01;

    # het point
    return 1 if $cnt{$alt} >= 3 && $cnt{$ref} >= 3 && $cnt{$alt}/$depth > 0.2;

    # likelihood of AA aa Aa
    my $lh_AA = hom_logLikelihood($bs,$q,$ref);
    my $lh_aa = hom_logLikelihood($bs,$q,$alt);
    my $lh_Aa = het_logLikelihood($bs,$q,$ref,$alt);

    print "AA $lh_AA, aa $lh_aa, Aa $lh_Aa, [prior_Aa ". $pr->{"$ref$alt"}."]\n" if $debug;

    my $i = $pr->{"$ref$alt"}; # het
    my $j = $pr->{"$ref$ref"};
    my $k = $pr->{"$alt$alt"};

    my $hood = 1 / (1 + ( exp($lh_AA-$lh_Aa)*$j + exp($lh_aa-$lh_Aa)*$k ) /$i);
    print "het ~~ $hood\n" if $debug;

    return $hood;
}

sub hom_logLikelihood {
    my ($bs, $q, $base) = @_;

    my $sum = 0;
    for my $i ( 0..$#$bs ) {
        my $e = errorRate( $q->[$i] );

        $sum += ( $bs->[$i] eq $base ? log( 1 - $e ) : log( $e/3 ) );
    }

    return $sum;
}

sub het_logLikelihood {
    my ($bs, $q, $ref, $alt) = @_;

    my $sum = 0;
    for my $i ( 0..$#$bs ) {
        my $e = errorRate( $q->[$i] );

        if ( $bs->[$i] eq $ref || $bs->[$i] eq $alt ) {
            $sum += log( 0.5 - $e/3 );
        }
        else {
            $sum += log( $e/3 );
        }
    }

    return $sum;
}

sub errorRate {
    return 10 ** ( - $_[0] / 10 );
}


sub showLog {
    my @t = localtime();
    printf STDERR "[%04d-%02d-%02d %02d:%02d:%02d]\t%s\n", $t[5] + 1900, $t[4] + 1, @t[3,2,1,0], $_[0];
}

__END__




=head1 Function
    Calculate het SNPs probability in given genome regions

=head1 Usage
    perl geno.pl -f ref.fa -m maf_f.gz -r region_f -i in.bam [Options] -o out

=head1 Options
    -h|-help                help
    -i|-inbam         [s]   inbam file
    -f|-fasta         [s]   reference.fa
    -m|-mafFile       [s]   bgziped maf file (chr, pos, ref, alt, maf)
    -r|-regionFile    [s]   region file (chr, start, end)
    -b|-baseQuaCutoff [i]   base quality cutoff [20]
    -d|-debug               debug mode
    -s|-seed          [i]   seed for srand [time]
    -o|-outfile       [s]   outfile

=head1 Author
    yerui@connect.hku.hk

=head1 Version
    v0.1;    2018-09-11

=cut





