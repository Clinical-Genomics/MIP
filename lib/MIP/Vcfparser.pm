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
use Set::IntervalTree;

## MIPs lib/
use MIP::Constants qw{ $SEMICOLON $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ build_interval_tree define_select_data_headers parse_feature_file_data parse_feature_file_header set_vcf_header_info };
}

sub build_interval_tree {

## Function : Creates the interval tree(s) for feature file types (e.g range and select)
## Returns  :
## Arguments: $feature_columns_ref => Feature annotations columns {REF}
##          : $feature_file_type   => Feature file type e.g. "select_feature"
##          : $line_elements_ref   => Feature file line elements {REF}
##          : $padding             => Nucleotide distance to pad the feature with
##          : $tree_href           => Interval tree hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s
    my $feature_columns_ref;
    my $feature_file_type;
    my $line_elements_ref;
    my $padding;
    my $tree_href;

    my $tmpl = {
        feature_columns_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$feature_columns_ref,
            strict_type => 1,
        },
        feature_file_type => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_type,
            strict_type => 1,
        },
        line_elements_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$line_elements_ref,
            strict_type => 1,
        },
        padding => {
            defined     => 1,
            required    => 1,
            store       => \$padding,
            strict_type => 1,
        },
        tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$tree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Coordinates to set in tree
    my $contig        = $line_elements_ref->[0];
    my $feature_start = $line_elements_ref->[1] - $padding;
    my $feature_stop  = $line_elements_ref->[2] + $padding;

    ## Features to collect (Format: ";" separated elements)
    my $feature_string = join $SEMICOLON,
      @{$line_elements_ref}[ @{$feature_columns_ref} ];

    ## Translate not allowed whitespace from INFO field to underscore
    $feature_string =~ tr/ /_/;

    ## Only create once per key and contig
    if ( not defined $tree_href->{$feature_file_type}{$contig} ) {

        ## Create tree
        $tree_href->{$feature_file_type}{$contig} =
          Set::IntervalTree->new();
    }

    ## Store feature_file_type, contig and ";" sep feature string in tree
    $tree_href->{$feature_file_type}{$contig}
      ->insert( $feature_string, $feature_start, $feature_stop );
    return;
}

sub define_select_data_headers {

## Function : Defines arbitrary INFO fields headers based on information in select file
## Returns  : %select_data
## Arguments: None

    my %select_data;

    $select_data{select_file}{HGNC_symbol}{info} =
      q{##INFO=<ID=HGNC_symbol,Number=.,Type=String,Description="The HGNC gene symbol">};
    $select_data{select_file}{Ensembl_gene_id}{info} =
q{##INFO=<ID=Ensembl_gene_id,Number=.,Type=String,Description="Ensembl gene identifier">};
    $select_data{select_file}{OMIM_morbid}{info} =
q{##INFO=<ID=OMIM_morbid,Number=.,Type=String,Description="OMIM morbid ID associated with gene(s)">};
    $select_data{select_file}{Phenotypic_disease_model}{info} =
q{##INFO=<ID=Phenotypic_disease_model,Number=.,Type=String,Description="Known disease gene(s) phenotype inheritance model">};
    $select_data{select_file}{Clinical_db_gene_annotation}{info} =
q{##INFO=<ID=Clinical_db_gene_annotation,Number=.,Type=String,Description="Gene disease group association">};
    $select_data{select_file}{Reduced_penetrance}{info} =
q{##INFO=<ID=Reduced_penetrance,Number=.,Type=String,Description="Pathogenic gene which can exhibit reduced penetrance">};
    $select_data{select_file}{Disease_associated_transcript}{info} =
q{##INFO=<ID=Disease_associated_transcript,Number=.,Type=String,Description="Known pathogenic transcript(s) for gene">};
    $select_data{select_file}{Ensembl_transcript_to_refseq_transcript}{info} =
q{##INFO=<ID=Ensembl_transcript_to_refseq_transcript,Number=.,Type=String,Description="The link between ensembl transcript and refSeq transcript IDs">};
    $select_data{select_file}{Gene_description}{info} =
q{##INFO=<ID=Gene_description,Number=.,Type=String,Description="The HGNC gene description">};
    $select_data{select_file}{Genetic_disease_model}{info} =
q{##INFO=<ID=Genetic_disease_model,Number=.,Type=String,Description="Known disease gene(s) inheritance model">};
    $select_data{select_file}{No_hgnc_symbol}{info} =
q{##INFO=<ID=No_hgnc_symbol,Number=.,Type=String,Description="Clinically relevant genetic regions lacking a HGNC_symbol or Ensembl gene ">};
    return %select_data;
}

sub set_vcf_header_info {

## Function : Adds arbitrary INFO fields to hash based on supplied header key
##            unless header key is already defined
## Returns  :
## Arguments: $feature_file_type  => Feature file key
##          : $feature_file_path => Feature file path
##          : $header_key        => Header key from feature file
##          : $meta_data_href    => Hash to store meta_data in {REF}
##          : $position          => Column position in supplied range file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $feature_file_type;
    my $feature_file_path;
    my $header_key;
    my $meta_data_href;
    my $position;

    my $tmpl = {
        feature_file_type => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_type,
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
        q{##INFO=<ID=}
      . $header_key
      . q{,Number=.,Type=String,Description="String taken from }
      . $feature_file_path . q{">};

    ## Add INFO from predefined entries
    if ( defined $meta_data_href->{$feature_file_type}{$header_key} ) {

        $meta_data_href->{present}{$header_key}{info} =
          $meta_data_href->{$feature_file_type}{$header_key}{info};
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

sub parse_feature_file_data {

## Function : Parse feature file data and build interval tree from the data
## Returns  :
## Arguments: $data_line                      => Data line
##          : $feature_columns_ref            => Feature columns to include {REF}
##          : $feature_data_href              => Feature file hash {REF}
##          : $feature_file_type               => Feature file key used to distinguish feature file(s) i.e., select or range
##          : $padding                        => Padding distance
##          : $select_feature_matching_column => Column in the select file to match with vcf key annotation {Optional}
##          : $tree_href                      => Interval tree hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $data_line;
    my $feature_columns_ref;
    my $feature_data_href;
    my $feature_file_type;
    my $padding;
    my $select_feature_matching_column;
    my $tree_href;

    my $tmpl = {
        data_line => {
            defined     => 1,
            required    => 1,
            store       => \$data_line,
            strict_type => 1,
        },
        feature_columns_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$feature_columns_ref,
            strict_type => 1,
        },
        feature_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$feature_data_href,
            strict_type => 1,
        },
        feature_file_type => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_type,
            strict_type => 1,
        },
        padding => {
            defined     => 1,
            required    => 1,
            store       => \$padding,
            strict_type => 1,
        },
        select_feature_matching_column =>
          { store => \$select_feature_matching_column, strict_type => 1, },
        tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$tree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Vcfparser qw{ build_interval_tree };

    ## Split data into array elements
    my @data_features = split $TAB, $data_line;

    if ( defined $select_feature_matching_column ) {

        my $data_feature = $data_features[$select_feature_matching_column];

        # Replace whitespace with underscore
        $data_feature =~ s/\s/_/gsxm;

        ## Set matching column data feature to feature data
        $feature_data_href->{$data_feature} = $data_feature;
    }

    ## Create Interval Tree
    if ( @{$feature_columns_ref} ) {

        ## Annotate vcf with features from feature file
        build_interval_tree(
            {
                feature_columns_ref => $feature_columns_ref,
                feature_file_type   => $feature_file_type,
                line_elements_ref   => \@data_features,
                padding             => $padding,
                tree_href           => $tree_href,
            }
        );
    }
    return 1;
}

sub parse_feature_file_header {

## Function : Get feature file header
## Returns  :
## Arguments: $feature_columns_ref => Feature columns to include {REF}
##          : $feature_data_href   => Feature file hash {REF}
##          : $feature_file_type    => Feature file key used to distinguish feature file(s) i.e., select or range
##          : $feature_file_path   => Feature file path
##          : $header_line         => Header line

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $feature_columns_ref;
    my $feature_data_href;
    my $feature_file_type;
    my $feature_file_path;
    my $header_line;

    my $tmpl = {
        feature_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$feature_data_href,
            strict_type => 1,
        },
        feature_columns_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$feature_columns_ref,
            strict_type => 1,
        },
        feature_file_type => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_type,
            strict_type => 1,
        },
        feature_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_path,
            strict_type => 1,
        },
        header_line => {
            defined     => 1,
            required    => 1,
            store       => \$header_line,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Vcfparser qw{ set_vcf_header_info };

    ## Split headers into array elements
    my @headers = split $TAB, $header_line;

    ## Defines what headers to store from feature file
    while ( my ( $feature_index, $feature_position ) = each @{$feature_columns_ref} ) {

        ## Alias
        my $header_key = $headers[$feature_position];

        set_vcf_header_info(
            {
                feature_file_type => $feature_file_type,
                feature_file_path => $feature_file_path,
                header_key        => $header_key,
                meta_data_href    => $feature_data_href,
                position          => $feature_index,
            }
        );
    }
    return 1;
}
1;
