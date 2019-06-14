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
    our @EXPORT_OK = qw{ define_select_data_headers set_vcf_header_info };
}

sub set_vcf_header_info {

## Function : Adds arbitrary INFO fields to hash based on supplied header key
##            unless header key is already defined
## Returns  :
## Arguments: $feature_file_key  => Feature file key
##          : $feature_file_path => Feature file path
##          : $header_key        => Header key from feature file
##          : $meta_data_href    => Hash to store meta_data in {REF}
##          : $position          => Column position in supplied range file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $feature_file_key;
    my $feature_file_path;
    my $header_key;
    my $meta_data_href;
    my $position;

    my $tmpl = {
        feature_file_key => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_key,
            strict_type => 1,
        },
        feature_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_path,
            strict_type => 1,
        },
        header_key => {
            defined     => 1,
            required    => 1,
            store       => \$header_key,
            strict_type => 1,
        },
        meta_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$meta_data_href,
            strict_type => 1,
        },
        position => {
            defined     => 1,
            required    => 1,
            store       => \$position,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## For not previously defined header keys in feature files definition
    my $arbitrary_info_field =
        q?##INFO=<ID=?
      . $header_key
      . q?,Number=.,Type=String,Description="String taken from ?
      . $feature_file_path . q?">?;

    ## Add INFO from predefined entries
    if ( defined $meta_data_href->{$feature_file_key}{$header_key} ) {

        $meta_data_href->{present}{$header_key}{info} =
          $meta_data_href->{$feature_file_key}{$header_key}{info};
    }
    else {
        ## Add arbitrary INFO field using feature file header key

        $meta_data_href->{present}{$header_key}{info} = $arbitrary_info_field;
    }

    ## Column position in supplied tsv feature file
    $meta_data_href->{present}{$header_key}{column_order} =
      $position;

    return;
}

sub define_select_data_headers {

## Function : Defines arbitrary INFO fields headers based on information in select file
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
