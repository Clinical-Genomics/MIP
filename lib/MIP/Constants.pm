package MIP::Constants;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported

    our @EXPORT_OK = qw{
      $AMPERSAND
      %ANALYSIS
      $ASTERISK
      $BACKTICK
      $BACKWARD_SLASH
      $CLOSE_BRACE
      $CLOSE_BRACKET
      $COLON
      $COMMA
      $DASH
      $DOLLAR_SIGN
      $DOT
      $DOUBLE_QUOTE
      $EMPTY_STR
      $EQUALS
      $ESCAPE
      $FORWARD_SLASH
      $LOG
      $MIP_VERSION
      $MOOSEX_APP_SCEEN_WIDTH
      $NEWLINE
      $OPEN_BRACE
      $OPEN_BRACKET
      $PIPE
      %PRIMARY_CONTIG
      $SEMICOLON
      $SINGLE_QUOTE
      %SO_CONSEQUENCE_SEVERITY
      $SPACE
      $TAB
      $UNDERSCORE
    };
}

## Constants

## Analysis
Readonly our %ANALYSIS => (
    ANNOTATION_DISTANCE    => 5000,
    ANNOTATION_DISTANCE_MT => 0,
    JAVA_GUEST_OS_MEMORY   => 4,
);

## Set MIP version
Readonly our $MIP_VERSION => q{v7.1.0};

## Cli
Readonly our $MOOSEX_APP_SCEEN_WIDTH => 160;

## Log
Readonly our $LOG => q{MIP_ANALYSE};

## SO-terms
Readonly our %SO_CONSEQUENCE_SEVERITY => (
    transcript_ablation               => { rank                      => 1 },
    transcript_ablation               => { genetic_region_annotation => q{exonic} },
    splice_donor_variant              => { rank                      => 2 },
    splice_donor_variant              => { genetic_region_annotation => q{splicing} },
    splice_acceptor_variant           => { rank                      => 2 },
    splice_acceptor_variant           => { genetic_region_annotation => q{splicing} },
    stop_gained                       => { rank                      => 3 },
    stop_gained                       => { genetic_region_annotation => q{exonic} },
    frameshift_variant                => { rank                      => 4 },
    frameshift_variant                => { genetic_region_annotation => q{exonic} },
    stop_lost                         => { rank                      => 5 },
    stop_lost                         => { genetic_region_annotation => q{exonic} },
    start_lost                        => { rank                      => 5 },
    start_lost                        => { genetic_region_annotation => q{exonic} },
    initiator_codon_variant           => { rank                      => 6 },
    initiator_codon_variant           => { genetic_region_annotation => q{exonic} },
    inframe_insertion                 => { rank                      => 6 },
    inframe_insertion                 => { genetic_region_annotation => q{exonic} },
    inframe_deletion                  => { rank                      => 6 },
    inframe_deletion                  => { genetic_region_annotation => q{exonic} },
    missense_variant                  => { rank                      => 6 },
    missense_variant                  => { genetic_region_annotation => q{exonic} },
    protein_altering_variant          => { rank                      => 6 },
    protein_altering_variant          => { genetic_region_annotation => q{exonic} },
    transcript_amplification          => { rank                      => 7 },
    transcript_amplification          => { genetic_region_annotation => q{exonic} },
    splice_region_variant             => { rank                      => 8 },
    splice_region_variant             => { genetic_region_annotation => q{splicing} },
    incomplete_terminal_codon_variant => { rank                      => 9 },
    incomplete_terminal_codon_variant => { genetic_region_annotation => q{exonic} },
    synonymous_variant                => { rank                      => 10 },
    synonymous_variant                => { genetic_region_annotation => q{exonic} },
    stop_retained_variant             => { rank                      => 10 },
    stop_retained_variant             => { genetic_region_annotation => q{exonic} },
    start_retained_variant            => { rank                      => 10 },
    start_retained_variant            => { genetic_region_annotation => q{exonic} },
    coding_sequence_variant           => { rank                      => 11 },
    coding_sequence_variant           => { genetic_region_annotation => q{exonic} },
    mature_miRNA_variant              => { rank                      => 12 },
    mature_miRNA_variant              => { genetic_region_annotation => q{ncRNA_exonic} },
    q{5_prime_UTR_variant}            => { rank                      => 13 },
    q{5_prime_UTR_variant}            => { genetic_region_annotation => q{5UTR} },
    q{3_prime_UTR_variant}            => { rank                      => 14 },
    q{3_prime_UTR_variant}            => { genetic_region_annotation => q{3UTR} },
    non_coding_transcript_exon_variant => { rank => 15 },
    non_coding_transcript_exon_variant =>
      { genetic_region_annotation => q{ncRNA_exonic} },
    non_coding_transcript_variant => { rank                      => 15 },
    non_coding_transcript_variant => { genetic_region_annotation => q{ncRNA} },
    intron_variant                => { rank                      => 16 },
    intron_variant                => { genetic_region_annotation => q{intronic} },
    NMD_transcript_variant        => { rank                      => 17 },
    NMD_transcript_variant        => { genetic_region_annotation => q{ncRNA} },
    upstream_gene_variant         => { rank                      => 18 },
    upstream_gene_variant         => { genetic_region_annotation => q{upstream} },
    downstream_gene_variant       => { rank                      => 19 },
    downstream_gene_variant       => { genetic_region_annotation => q{downstream} },
    TFBS_ablation                 => { rank                      => 20 },
    TFBS_ablation                 => { genetic_region_annotation => q{TFBS} },
    TFBS_amplification            => { rank                      => 21 },
    TFBS_amplification            => { genetic_region_annotation => q{TFBS} },
    TF_binding_site_variant       => { rank                      => 22 },
    TF_binding_site_variant       => { genetic_region_annotation => q{TFBS} },
    regulatory_region_variant     => { rank                      => 22 },
    regulatory_region_variant  => { genetic_region_annotation => q{regulatory_region} },
    regulatory_region_ablation => { rank                      => 23 },
    regulatory_region_ablation => { genetic_region_annotation => q{regulatory_region} },
    regulatory_region_amplification => { rank => 24 },
    regulatory_region_amplification =>
      { genetic_region_annotation => q{regulatory_region} },
    feature_elongation => { rank                      => 25 },
    feature_elongation => { genetic_region_annotation => q{genomic_feature} },
    feature_truncation => { rank                      => 26 },
    feature_truncation => { genetic_region_annotation => q{genomic_feature} },
    intergenic_variant => { rank                      => 27 },
    intergenic_variant => { genetic_region_annotation => q{intergenic} },
);

## Contigs
Readonly our %PRIMARY_CONTIG => (
    38 => {
        contigs => [
            qw{
              chr1 chr2 chr3 chr4 chr5 chr6
              chr7 chr8 chr9 chr10 chr11 chr12
              chr13 chr14 chr15 chr16 chr17 chr18
              chr19 chr20 chr21 chr22 chrX chrY
              chrM },
        ],
        contigs_size_ordered => [
            qw{
              chr1 chr2 chr3 chr4 chr5 chr6
              chr7 chrX chr8 chr9 chr10 chr11
              chr12 chr13 chr14 chr15 chr16 chr17
              chr18 chr19 chr20 chr21 chr22 chrY
              chrM }
        ],
        synonyms_map => {
            chr1  => 1,
            chr2  => 2,
            chr3  => 3,
            chr4  => 4,
            chr5  => 5,
            chr6  => 6,
            chr7  => 7,
            chr8  => 8,
            chr9  => 9,
            chr10 => 10,
            chr11 => 11,
            chr12 => 12,
            chr13 => 13,
            chr14 => 14,
            chr15 => 15,
            chr16 => 16,
            chr17 => 17,
            chr18 => 18,
            chr19 => 19,
            chr20 => 20,
            chr21 => 21,
            chr22 => 22,
            chrX  => q{X},
            chrY  => q{Y},
            chrM  => q{M},
        },
    },
    37 => {
        contigs => [
            qw{
              1 2 3 4 5 6 7 8 9 10
              11 12 13 14 15 16 17 18 19 20
              21 22 X Y MT },
        ],
        contigs_size_ordered => [
            qw{
              1 2 3 4 5 6 7 X 8 9
              10 11 12 13 14 15 16 17 18 19
              20 21 22 Y MT },
        ],
    },
    19 => {
        contigs => [
            qw{
              chr1 chr2 chr3 chr4 chr5 chr6
              chr7 chr8 chr9 chr10 chr11 chr12
              chr13 chr14 chr15 chr16 chr17 chr18
              chr19 chr20 chr21 chr22 chrX chrY
              chrM },
        ],
        contigs_size_ordered => [
            qw{
              chr1 chr2 chr3 chr4 chr5 chr6
              chr7 chrX chr8 chr9 chr10 chr11
              chr12 chr13 chr14 chr15 chr16 chr17
              chr18 chr19 chr20 chr21 chr22 chrY
              chrM }
        ],
    },
);

## Symbols
Readonly our $AMPERSAND      => q{&};
Readonly our $ASTERISK       => q{*};
Readonly our $BACKTICK       => q{`};
Readonly our $BACKWARD_SLASH => q{\\};
Readonly our $CLOSE_BRACE    => q{\}};
Readonly our $CLOSE_BRACKET  => q{]};
Readonly our $COLON          => q{:};
Readonly our $COMMA          => q{,};
Readonly our $DASH           => q{-};
Readonly our $DOLLAR_SIGN    => q{$};
Readonly our $DOT            => q{.};
Readonly our $DOUBLE_QUOTE   => q{"};
Readonly our $EMPTY_STR      => q{};
Readonly our $ESCAPE         => q{\\};
Readonly our $EQUALS         => q{=};
Readonly our $FORWARD_SLASH  => q{/};
Readonly our $NEWLINE        => qq{\n};
Readonly our $OPEN_BRACE     => q{\{};
Readonly our $OPEN_BRACKET   => q{[};
Readonly our $PIPE           => q{|};
Readonly our $SEMICOLON      => q{;};
Readonly our $SINGLE_QUOTE   => q{'};
Readonly our $SPACE          => q{ };
Readonly our $TAB            => qq{\t};
Readonly our $UNDERSCORE     => q{_};

1;
