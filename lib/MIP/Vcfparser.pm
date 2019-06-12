package MIP::Vcfparser;

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

## MIPs lib/
use MIP::Constants qw{ $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ add_vcf_header_info define_select_data };
}

sub add_vcf_header_info {

## Function : Adds arbitrary INFO fields to hash based on supplied headers unless header is already defined.
## Returns  :
## Arguments: $header_ref          => Header from range file {REF}
##          : $meta_data_href      => Hash to store meta_data in {REF}
##          : $position_ref        => Column position in supplied range file {REF}
##          : $range_file_key      => Range file key
##          : $range_file_path_ref => Range file path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $header_ref;
    my $meta_data_href;
    my $position_ref;
    my $range_file_key;
    my $range_file_path_ref;

    my $tmpl = {
      header_ref => {
          defined     => 1,
          default     => \$$,
          required    => 1,
          store       => \$header_ref,
          strict_type => 1,
      },
        meta_data_href => {
          default     => {},
            defined     => 1,
            required    => 1,
            store       => \$meta_data_href,
            strict_type => 1,
        },
        position_ref => {
          default     => \$$,
            defined     => 1,
            required    => 1,
            store       => \$position_ref,
            strict_type => 1,
        },
        range_file_key => {
            defined     => 1,
            required    => 1,
            store       => \$range_file_key,
            strict_type => 1,
        },
        range_file_path_ref => {
            default     => \$$,
            defined     => 1,
            required    => 1,
            store       => \$range_file_path_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add INFO from predefined entries
    if ( defined $meta_data_href->{$range_file_key}{$$header_ref} )
    {

        $meta_data_href->{present}{$$header_ref}{info} =
          $meta_data_href->{$range_file_key}{$$header_ref}{info};

        ## Column position in supplied range input file
        $meta_data_href->{present}{$$header_ref}{column_order} =
          $$position_ref;
    }
    else {
        ## Add arbitrary INFO field using input header

        $meta_data_href->{present}{$$header_ref}{info} =
            q?##INFO=<ID=?
          . $$header_ref
          . q?,Number=.,Type=String,Description="String taken from ?
          . $$range_file_path_ref . q?">?;

        ## Column position in supplied -sf_ac
        $meta_data_href->{present}{$$header_ref}{column_order} =
          $$position_ref;
    }
    return;
}

sub define_select_data {

## Function : Defines arbitrary INFO fields based on headers in select file
## Returns  : %select_data
## Arguments: None

    my %select_data;

    $select_data{select_file}{HGNC_symbol}{info} =
      q?##INFO=<ID=HGNC_symbol,Number=.,Type=String,Description="The HGNC gene symbol">?;
    $select_data{select_file}{Ensembl_gene_id}{info} =
q?##INFO=<ID=Ensembl_gene_id,Number=.,Type=String,Description="Ensembl gene identifier">?;
    $select_data{select_file}{OMIM_morbid}{info} =
q?##INFO=<ID=OMIM_morbid,Number=.,Type=String,Description="OMIM morbid ID associated with gene(s)">?;
    $select_data{select_file}{Phenotypic_disease_model}{info} =
q?##INFO=<ID=Phenotypic_disease_model,Number=.,Type=String,Description="Known disease gene(s) phenotype inheritance model">?;
    $select_data{select_file}{Clinical_db_gene_annotation}{info} =
q?##INFO=<ID=Clinical_db_gene_annotation,Number=.,Type=String,Description="Gene disease group association">?;
    $select_data{select_file}{Reduced_penetrance}{info} =
q?##INFO=<ID=Reduced_penetrance,Number=.,Type=String,Description="Pathogenic gene which can exhibit reduced penetrance">?;
    $select_data{select_file}{Disease_associated_transcript}{info} =
q?##INFO=<ID=Disease_associated_transcript,Number=.,Type=String,Description="Known pathogenic transcript(s) for gene">?;
    $select_data{select_file}{Ensembl_transcript_to_refseq_transcript}{info} =
q?##INFO=<ID=Ensembl_transcript_to_refseq_transcript,Number=.,Type=String,Description="The link between ensembl transcript and refSeq transcript IDs">?;
    $select_data{select_file}{Gene_description}{info} =
q?##INFO=<ID=Gene_description,Number=.,Type=String,Description="The HGNC gene description">?;
    $select_data{select_file}{Genetic_disease_model}{info} =
q?##INFO=<ID=Genetic_disease_model,Number=.,Type=String,Description="Known disease gene(s) inheritance model">?;
    $select_data{select_file}{No_hgnc_symbol}{info} =
q?##INFO=<ID=No_hgnc_symbol,Number=.,Type=String,Description="Clinically relevant genetic regions lacking a HGNC_symbol or Ensembl gene ">?;
    return %select_data;
}

1;
