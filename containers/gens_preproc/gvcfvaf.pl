#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use List::Util qw(max);
use File::Basename qw(dirname);

die "Give GVCF as only argument!" unless -s $ARGV[0];

my $SCRIPT_ROOT = dirname($0);

my @chr_order = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT);
my %chr_order_hash;
foreach my $i (0..$#chr_order ) {
    $chr_order_hash{$chr_order[$i]} = $i;
}

open(GNOMAD, $ARGV[1]);

open(GVCF, "zcat $ARGV[0]|");
my $gvcf_line;
while($gvcf_line = <GVCF>) {
    last unless $gvcf_line =~ /^#/;
}

#my $gvcf = parse_gvcf_entry($gvcf_line);
my $gvcf = gvcf_position($gvcf_line);
my $skipped = 0;
while(<GNOMAD>) {
    chomp;
    my( $gnomad_chr, $gnomad_pos ) = split /\t/;
    print STDERR "Skipping gnomad $gnomad_chr $gnomad_pos ".$gvcf->{chr}." ".$gvcf->{start}."\n" if chr_less($gnomad_chr, $gvcf->{chr});
    next if chr_less($gnomad_chr, $gvcf->{chr});
    if( $gvcf->{start} < $gnomad_pos or chr_less($gvcf->{chr}, $gnomad_chr) ) {
	my $prev_gvcf_chr = $gvcf->{chr};
	do {
	    $gvcf_line = <GVCF>;
	    #$gvcf = parse_gvcf_entry($a);
	    $gvcf = gvcf_position($gvcf_line);
	} until eof(GVCF) or ( $gvcf->{end} >= $gnomad_pos or $gvcf->{start} >= $gnomad_pos or $gvcf->{chr} ne $prev_gvcf_chr );
    }
    if( $gnomad_pos >= $gvcf->{start} and $gnomad_chr eq $gvcf->{chr} ) {
	if( $gnomad_pos <= $gvcf->{end} ) {
	    $gvcf = parse_gvcf_entry($gvcf_line);
	    print $gnomad_chr."\t".$gnomad_pos."\t".($gvcf->{frq} or 0)."\n" if defined $gvcf->{frq};
	}
	else {
	    $skipped ++;
	}
    }
    last if eof(GVCF);
}   
print STDERR "$skipped variants skipped!\n";


sub chr_less {
    my( $chr1, $chr2 ) = @_;
    return $chr_order_hash{$chr1} < $chr_order_hash{$chr2};
}

sub gvcf_position {
    my $str = shift;
    my %data;
    my @a = split /\t/, $str;
    $data{chr} = $a[0];
    $data{start} = $a[1];
    $a[7] =~ /(^|;)END=(.*?)(;|$)/;
    $data{end} = ($2 or $a[1]);
    return \%data;
}
 
sub parse_gvcf_entry {
    my $str = shift;
    chomp $str;
    my %data;
    my @a = split /\t/, $str;
    $data{str} = $str;
    $data{chr} = $a[0];
    $data{start} = $a[1];
    $a[7] =~ /(^|;)END=(.*?)(;|$)/;
    $data{end} = ($2 or $a[1]);

    return \%data if length($a[3]) > 1;


    unless( $a[8] =~ /:AD:/ ) {
	$data{frq} = 0;
    }
    else {
	my @ALT = split /,/, $a[4];
    	my @fmt = split /:/, $a[8];
	my @sam = split /:/, $a[9];
	my( $alt, $alt_cnt, $dp );
	for my $i (0..@fmt) {
	    if( $fmt[$i] eq "GT" ) {
		my( $a, $b ) = (split /\//, $sam[$i]);
		$alt = $b;
		return \%data if $alt != 0 and( !defined $ALT[$alt-1] or length($ALT[$alt-1]) > 1);	
	    }
	}
	for my $i (0..@fmt) {
	    if( $fmt[$i] eq "AD" ) {
		if( $alt != 0 ) {
		    $alt_cnt = (split /,/, $sam[$i])[$alt];
		}
		else {
		    my @cnts = split /,/, $sam[$i];
		    shift @cnts;
		    pop @cnts;
		    $alt_cnt = max(@cnts);
		}
	    }
	}
	for my $i (0..@fmt) {
	    if( $fmt[$i] eq "DP" ) {
		$dp = $sam[$i];
		last;
	    }   
	}
	return \%data if $dp < 10;
	$data{frq} = $alt_cnt/$dp;
	$data{alt} = $ALT[$alt];
	$data{ref} = $a[3];
	$data{all} = $a[9];
    }
    
    
    return \%data;
}
