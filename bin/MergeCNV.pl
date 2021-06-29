#!/usr/bin/perl -w
use strict;
use File::Spec;
use File::Basename;
use Getopt::Long;

my ( $Help, $sample, $cnvnator, $delly, $lumpy, $seeksv, $type, $outdir );
GetOptions(
        'help|?'     => \$Help,
        'sample=s'   => \$sample,
        'CNVnator=s' => \$cnvnator,
        'Delly=s'    => \$delly,
        'Lumpy=s'    => \$lumpy,
        'Seeksv=s'   => \$seeksv,
        'Type=s'     => \$type,
        'outdir=s'   => \$outdir
);

die `pod2text $0` if ( $Help or !$sample or !$cnvnator or !$delly or !$lumpy or !$seeksv or !$type or !$outdir );

#Configurate file
my $path_curf = File::Spec->rel2abs(__FILE__);
my ( $vol, $dirs, $file ) = File::Spec->splitpath($path_curf);
my $home = dirname($dirs);
my $bedtools = "$home/bin/bedtools";

open OUT, ">$outdir/$sample.$type.overlap";

my %hash = (
    $cnvnator => "CNVnator",
    $delly => "Delly",
    $lumpy => "Lumpy",
    $seeksv => "Seeksv",
);

my @files = ();
for my $f ( $cnvnator, $delly, $lumpy, $seeksv ) {
    if ( ! -z $f ) {
        push @files, $f;
    }
}
for my $fs (0..$#files) {
    $" = ",";
    my @of = ($files[$fs]);
    my @tools = ($hash{$files[$fs]});
    for my $other (0..$#files) {
        if ( $other != $fs ) {
            push @of, $files[$other];
            push @tools, $hash{$files[$other]};
        }
    }
    overlap (@of,@tools);
}

my %dedup;
my $count = 0;
sub overlap {
    my ( @bed ) = @_; $count++;
    my $bfiles = "";
    my %software1;
    for my $n (1..($#bed+1)/2-1) {
        $bfiles = "$bfiles $bed[$n]";
        my $s = $n + ($#bed+1)/2;
        $software1{$n} = $bed[$s];
    }
    open IN1, "$bedtools intersect -wao -f 0.5 -F 0.5 -a $bed[0] -b $bfiles |" or die $!;
    my %sf1; my %hash1;

    while (<IN1>) {
        chomp;
        my @info = split;
        next if ( $info[0] =~ /Un/ );
        next if ( $info[0] =~ /random/ );
        next if ( $info[0] !~ /chr/ );
        my $pos = "$info[0]\t$info[1]\t$info[2]";
        if ( $info[5] == -1 ) {
            $hash1{$pos} = $pos;
        }
        else {
            if ( !$hash1{$pos} ) {
                my @tmp = ();
                if ( $hash{$bed[0]} ne "CNVnator" and $hash{$bed[$info[3]]} ne "CNVnator" ) {
                    @tmp = ( $info[1], $info[2], $info[5], $info[6] );
                }
                elsif ( $hash{$bed[0]} eq "CNVnator" and $hash{$bed[$info[3]]} ne "CNVnator" ) {
                    @tmp = ( $info[5], $info[6] );
                }
                elsif ( $hash{$bed[0]} ne "CNVnator" and $hash{$bed[$info[3]]} eq "CNVnator" ) {
                    @tmp = ( $info[1], $info[2] );
                }
                @tmp = sort {$a<=>$b} @tmp;
                $sf1{$pos} = "$software1{$info[3]}";
                $hash1{$pos} = "$info[0]\t$tmp[0]\t$tmp[-1]";
            }
            else {
                my @tmp = ();
                my @pre = split /\t/, $hash1{$pos};
                if ( $hash{$bed[0]} ne "CNVnator" and $hash{$bed[$info[3]]} ne "CNVnator" ) {
                    @tmp = ( $pre[1], $pre[2], $info[1], $info[2], $info[5], $info[6] );
                }
                elsif ( $hash{$bed[0]} eq "CNVnator" and $hash{$bed[$info[3]]} ne "CNVnator" ) {
                    @tmp = ( $pre[1], $pre[2], $info[5], $info[6] );
                }
                elsif ( $hash{$bed[0]} ne "CNVnator" and $hash{$bed[$info[3]]} eq "CNVnator" ) {
                    @tmp = ( $pre[1], $pre[2], $info[1], $info[2] );
                }
                @tmp = sort {$a<=>$b} @tmp;
                $sf1{$pos} = "$sf1{$pos},$software1{$info[3]}";
                $hash1{$pos} = "$info[0]\t$tmp[0]\t$tmp[-1]";
            }
        }
    }
    my $ft = ($#bed+1)/2;
    for my $k ( sort keys %hash1 ) {
        if ( !$sf1{$k} ) {
            $sf1{$k} = $bed[$ft];
        }
        else {
            $sf1{$k} = "$bed[$ft],$sf1{$k}";
        }
        my @sf1s = split /,/, $sf1{$k};
        @sf1s = sort @sf1s;
        my %dedup1;
        @sf1s = grep { ++$dedup1{$_} < 2 } @sf1s;
        my $join = join ",", @sf1s;
        if ( !$dedup{"$hash1{$k}\t$join"} ) {
            my @detail = split /\t/, $hash1{$k};
            my $len = $detail[2] - $detail[1] + 1;
            print OUT "$hash1{$k}\t$len\t$join\n";
            $dedup{"$hash1{$k}\t$join"} = 1;
        }
    }
    close IN1;
}

close OUT;

`less $outdir/$sample.$type.overlap |sort -k1,1 -k2,2n > $outdir/$sample.$type.overlap.tmp`;

open IN, "$bedtools intersect -wao -f 0.5 -F 0.5 -a $outdir/$sample.$type.overlap.tmp -b $outdir/$sample.$type.overlap.tmp |" or die $!;
open OUT, ">$outdir/$sample.$type.overlap";
my %f_dedup;
while (<IN>) {
    chomp;
    my @info = split;
    if ( "$info[0]\t$info[1]\t$info[2]\t$info[4]" ne "$info[5]\t$info[6]\t$info[7]\t$info[9]" ) {
        $f_dedup{"$info[0]\t$info[1]\t$info[2]\t$info[4]"} = 1;
        $f_dedup{"$info[5]\t$info[6]\t$info[7]\t$info[9]"} = 1;
        next;
    }
    next if ( $f_dedup{"$info[0]\t$info[1]\t$info[2]\t$info[4]"} or $f_dedup{"$info[5]\t$info[6]\t$info[7]\t$info[9]"} );
    print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\n";
}

`rm -f $outdir/$sample.$type.overlap.tmp`;

=head1 Function
    Combine CNVs called from CNVnator, Delly, Lumpy and Seeksv

=head1 Usage
    perl MergeCNV.pl [Options]

=head1 Options
    -h|-help            help
    -s|-sample   [s]    sample name
    -C|-CNVnator [s]    CNVs called by CNVnator (file: chr, start, end)
    -D|-Delly    [s]    CNVs called by Delly (file: chr, start, end)
    -L|-Lumpy    [s]    CNVs called by Lumpy (file: chr, start, end)
    -S|-Seeksv   [s]    CNVs called by Seeksv (file: chr, start, end)
    -t|-Type     [s]    CNV type (DEL|DUP)
    -o|-outdir   [s]    output direction

=head1 Author
    sunnyzxh@connect.hku.hk

=head1 Version
    v0.1;  2019.3.12

=cut
