package MIP::Set::Contigs;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ set_contigs };

}

sub set_contigs {

## set_contigs

## Function : Set contig prefix and contig names depending on reference used.
## Returns  :
## Arguments: $file_info_href, $human_genome_reference
##          : $file_info_href         => File info hash {REF}
##          : $human_genome_reference => Human genome reference

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $human_genome_reference;

    my $tmpl = {
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        human_genome_reference => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$human_genome_reference
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Refseq - prefix and M
    if ( $human_genome_reference =~ / hg\d+ /sxm ) {

        # Contigs for filtering of bam file
        @{ $file_info_href->{bam_contigs} } = qw{
          chr1 chr2 chr3 chr4 chr5 chr6
          chr7 chr8 chr9 chr10 chr11 chr12
          chr13 chr14 chr15 chr16 chr17 chr18
          chr19 chr20 chr21 chr22 chrX chrY
          chrM };

        # Contigs for filtering of bam file
        @{ $file_info_href->{bam_contigs_size_ordered} } = qw{
          chr1 chr2 chr3 chr4 chr5 chr6
          chr7 chrX chr8 chr9 chr10 chr11
          chr12 chr13 chr14 chr15 chr16 chr17
          chr18 chr19 chr20 chr21 chr22 chrY
          chrM };

        # Contigs for filtering of vcf file
        @{ $file_info_href->{contigs} } = qw{
          chr1 chr2 chr3 chr4 chr5 chr6
          chr7 chr8 chr9 chr10 chr11 chr12
          chr13 chr14 chr15 chr16 chr17 chr18
          chr19 chr20 chr21 chr22 chrX chrY
          chrM };

        # Contigs for filtering of vcf file
        @{ $file_info_href->{contigs_size_ordered} } = qw{
          chr1 chr2 chr3 chr4 chr5 chr6
          chr7 chrX chr8 chr9 chr10 chr11
          chr12 chr13 chr14 chr15 chr16 chr17
          chr18 chr19 chr20 chr21 chr22 chrY
          chrM };
    }
    elsif ( $human_genome_reference =~ / grch\d+ /xsm ) {
        ## Ensembl - no prefix and MT

        # Contigs for filtering of bam file
        @{ $file_info_href->{bam_contigs} } = qw{
          1 2 3 4 5 6 7 8 9 10
          11 12 13 14 15 16 17 18 19 20
          21 22 X Y MT };

        # Contigs for filtering of bam file
        @{ $file_info_href->{bam_contigs_size_ordered} } = qw{
          1 2 3 4 5 6 7 X 8 9
          10 11 12 13 14 15 16 17 18 19
          20 21 22 Y MT };

        # Contigs for filtering of vcf file
        @{ $file_info_href->{contigs} } = qw{
          1 2 3 4 5 6 7 8 9 10
          11 12 13 14 15 16 17 18 19 20
          21 22 X Y MT };

        # Contigs for filtering of vcf file
        @{ $file_info_href->{contigs_size_ordered} } = qw{
          1 2 3 4 5 6 7 X 8 9
          10 11 12 13 14 15 16 17 18 19
          20 21 22 Y MT };

    }
    return;
}

1;
