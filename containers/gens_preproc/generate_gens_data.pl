#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw(sum);
use File::Basename qw(dirname);

if ($ARGV[0] eq "--version") {
    print "generate_gens_data.pl 1.0.2\n";
    exit 0;
}

my $SCRIPT_ROOT = dirname($0);

my @COV_WINDOW_SIZES = (100000, 25000, 5000, 1000, 100);
my @BAF_SKIP_N = (160, 40, 10, 4, 1);
my @PREFIXES = qw( o a b c d );
my $cov_fn = $ARGV[0];
my $gvcf_fn = $ARGV[1];

my $SAMPLE_ID = $ARGV[2];
my $GNOMAD = $ARGV[3];

my $COV_OUTPUT = $SAMPLE_ID.".cov.bed";
my $BAF_OUTPUT = $SAMPLE_ID.".baf.bed";

print STDERR "Calculating coverage data\n";
# Calculate coverage data
open( COVOUT, ">".$COV_OUTPUT );
for my $i (0..$#COV_WINDOW_SIZES) {
    generate_cov_bed($cov_fn, $COV_WINDOW_SIZES[$i], $PREFIXES[$i]);
}
close COVOUT;


print STDERR "Calculating BAFs from gvcf...\n";
# Calculate BAFs
system( 
    $SCRIPT_ROOT."/gvcfvaf.pl ".
    "$gvcf_fn $GNOMAD > baf.tmp"
    );
open( BAFOUT, ">".$BAF_OUTPUT );
for my $i (0..$#BAF_SKIP_N) {
    print STDERR "Outputting BAF $PREFIXES[$i]...\n";
    generate_baf_bed("baf.tmp", $BAF_SKIP_N[$i], $PREFIXES[$i]);
}
close BAFOUT;


system("bgzip -f -\@10 $BAF_OUTPUT");
system("tabix -f -p bed $BAF_OUTPUT.gz");
system("bgzip -f -\@10 $COV_OUTPUT");
system("tabix -f -p bed $COV_OUTPUT.gz");
unlink("baf.tmp");

sub generate_baf_bed {
    my( $fn, $skip, $prefix ) = @_;
    open( my $fh, $fn );
    my $i = 0;
    while(<$fh>) {
	if( $i++ % $skip == 0 ) {
	    chomp;
	    my @a = split /\t/;
	    print BAFOUT $prefix."_".$a[0]."\t".($a[1]-1)."\t".$a[1]."\t".$a[2]."\n";
	}
    }
    close $fh;
}

sub generate_cov_bed {

    my( $fn, $win_size, $prefix ) = @_;
    
    open(my $fh, $fn);
    my( $reg_start, $reg_end, $reg_chr, $force_end );
    my @reg_ratios;
    while(<$fh>) {
	next if /^@/ or /^CONTIG/;
	chomp;
	my ($chr, $start, $end, $ratio ) = split /\t/;
	my $orig_end = $end;
	unless( $reg_start ) {
	    $reg_start = $start;
	    $reg_end = $end;
	    $reg_chr = $chr;
	}

	if( $chr eq $reg_chr ) {
	    if( $start - $reg_end < $win_size ) {
		push @reg_ratios, $ratio;
		$reg_end = $end;
	    }

	    # If there is a large gap to the next region, prematurely end region
	    else {
		$force_end = 1;
		$end = $reg_end;
	    }
	}
	else {
	    $force_end = 1;
	    $end = $reg_end;
	}
	if( $end - $reg_start + 1 >= $win_size or $force_end ) {
	    my $mid_point = $reg_start + int(($end - $reg_start)/2);
	    print COVOUT $prefix."_".$reg_chr."\t".($mid_point-1)."\t".$mid_point."\t".mean(@reg_ratios)."\n";
	    undef $reg_start;
	    undef $reg_end;
	    undef $reg_chr;
	    undef @reg_ratios;
	}

	if( $force_end ) {
	    $reg_start = $start;
	    $reg_end = $orig_end;
	    $reg_chr = $chr;
	    push @reg_ratios, $ratio;
	    undef $force_end;
	}
    }
    close $fh;
}
    
sub mean {
    return sum(@_)/@_;
}
