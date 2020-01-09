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
use MIP::Constants
  qw{ $AMPERSAND $COLON $COMMA $EMPTY_STR $EQUALS $NEWLINE $PIPE $SEMICOLON $SINGLE_QUOTE %SO_CONSEQUENCE_SEVERITY $SPACE $TAB $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      add_feature_file_meta_data_to_vcf
      add_most_severe_csq_to_feature
      add_program_to_meta_data_header
      add_transcript_to_feature_file
      build_interval_tree
      check_data_terms
      %CSQ_FIELD_MAP
      define_select_data_headers
      parse_consequence
      parse_vcf_format_line
      parse_vep_csq
      parse_vep_csq_consequence
      parse_vep_csq_schema
      parse_vep_csq_transcripts
      set_most_severe_ann_to_vcf_record
      set_most_severe_pli
      write_feature_file_csq
      write_info_field
      write_info_addition_fields
      write_line_elements
      write_meta_data
    };
}

## Constants
Readonly our %CSQ_FIELD_MAP => (
    Allele      => q{allele},
    Consequence => q{consequence_field},
    Feature     => q{transcript_id},
    HGNC_ID     => q{hgnc_id},
    SYMBOL      => q{hgnc_symbol},
);

sub add_feature_file_meta_data_to_vcf {

## Function : Adds feature file meta data annotation headers to meta data hash.
## Returns  :
## Arguments: $data_href                      => Range file hash {REF}
##          : $feature_annotation_columns_ref => Range columns to include {REF}
##          : $file_key                       => Range file key used to seperate range file(s) i.e., select and range
##          : $meta_data_href                 => Vcf meta data {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $data_href;
    my $feature_annotation_columns_ref;
    my $file_key;
    my $meta_data_href;

    my $tmpl = {
        data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$data_href,
            strict_type => 1,
        },
        feature_annotation_columns_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$feature_annotation_columns_ref,
            strict_type => 1,
        },
        file_key => {
            defined     => 1,
            required    => 1,
            store       => \$file_key,
            strict_type => 1,
        },
        meta_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$meta_data_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Feature file annotations
    return if ( not @{$feature_annotation_columns_ref} );

  ANNOTATION:
    for my $annotation ( keys %{ $data_href->{present} } ) {

        ## Skip if INFO header is already present in VCF file
        next ANNOTATION if ( defined $meta_data_href->{INFO}{$annotation} );

        ## Add specific feature INFO meta data header
        $meta_data_href->{$file_key}{INFO}{$annotation} =
          $data_href->{present}{$annotation}{INFO};
    }
    return;
}

sub add_most_severe_csq_to_feature {

## Function : Set the most severe consequence and transcript for variant genes
## Returns  :
## Arguments: $hgnc_id                  => Hgnc id
##          : $most_severe_consequence  => Most severe consequence
##          : $most_severe_feature_href => Store most severe annotation per feature
##          : $most_severe_transcript   => Most severe transcript
##          : $per_gene                 => Only collect most severe transcript per gene
##          : $select_data_href         => Select file data {REF}
##          : $vcf_record_href          => VCF record {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $hgnc_id;
    my $most_severe_consequence;
    my $most_severe_feature_href;
    my $most_severe_transcript;
    my $select_data_href;
    my $vcf_record_href;

    ## Default(s)
    my $per_gene;

    my $tmpl = {
        hgnc_id => {
            required    => 1,
            store       => \$hgnc_id,
            strict_type => 1,
        },
        most_severe_consequence => {
            required    => 1,
            store       => \$most_severe_consequence,
            strict_type => 1,
        },
        most_severe_feature_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$most_severe_feature_href,
            strict_type => 1,
        },
        most_severe_transcript => {
            required    => 1,
            store       => \$most_severe_transcript,
            strict_type => 1,
        },
        per_gene => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$per_gene,
            strict_type => 1,
        },
        select_data_href => {
            default     => {},
            store       => \$select_data_href,
            strict_type => 1,
        },
        vcf_record_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vcf_record_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Exists in selected features
    if ( $select_data_href->{$hgnc_id} ) {

        push @{ $most_severe_feature_href->{select} }, $most_severe_consequence;
    }
    push @{ $most_severe_feature_href->{range} }, $most_severe_consequence;

    return if ( not $per_gene );

    ## Add to vcf record
    push @{ $vcf_record_href->{range_transcripts} }, $most_severe_transcript;

    return if ( not $select_data_href->{$hgnc_id} );

    push @{ $vcf_record_href->{select_transcripts} }, $most_severe_transcript;
    return;
}

sub add_program_to_meta_data_header {

## Function : Adds the program version and run date to the vcf meta data header hash
## Returns  :
## Arguments: $add_software_tag  => Write software tag to vcf header switch
##          : $meta_data_href    => Vcf meta data {REF}
##          : $vcfparser_version => Vcfparser version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $add_software_tag;
    my $meta_data_href;
    my $vcfparser_version;

    my $tmpl = {
        add_software_tag => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$add_software_tag,
            strict_type => 1,
        },
        meta_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$meta_data_href,
            strict_type => 1,
        },
        vcfparser_version => {
            defined     => 1,
            required    => 1,
            store       => \$vcfparser_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use DateTime;
    use File::Basename qw{ basename };

    return if ( not $add_software_tag );

    ## Get current date
    my $current_date = DateTime->now->ymd;

    ## Get script name
    my $program_name = basename($PROGRAM_NAME);

    ## Add to meta_data_href
    $meta_data_href->{software}{$program_name} =
      qq{##Software=<ID=$program_name,Version=$vcfparser_version,Date=$current_date};
    return;
}

sub add_transcript_to_feature_file {

## Function : Adds INFO key value pairs to record hash
## Returns  :
## Arguments: $hgnc_id          => Hgnc id
##          : $select_data_href => Select file data {REF}
##          : $transcripts_ref  => Transcript {REF}
##          : $vcf_record_href  => Hash for variant line data {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $hgnc_id;
    my $select_data_href;
    my $transcripts_ref;
    my $vcf_record_href;

    my $tmpl = {
        hgnc_id => {
            store       => \$hgnc_id,
            strict_type => 1,
        },
        select_data_href => {
            default     => {},
            store       => \$select_data_href,
            strict_type => 1,
        },
        transcripts_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$transcripts_ref,
            strict_type => 1,
        },
        vcf_record_href => {
            default  => {},
            defined  => 1,
            required => 1,
            store    => \$vcf_record_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add all transcripts to range transcripts
    push @{ $vcf_record_href->{range_transcripts} }, @{$transcripts_ref};

    ## Do not add to select feature
    return if ( not keys %{$select_data_href} or not defined $hgnc_id );

    ## Return if gene is not part of selected features
    return if ( not $select_data_href->{$hgnc_id} );

    ## Add all transcripts to selected transcripts
    push @{ $vcf_record_href->{select_transcripts} }, @{$transcripts_ref};

    return;
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

sub check_data_terms {

## Function : Check the found terms in the vcf correspond to known terms - otherwise croak and exit.
## Returns  :
## Arguments: $data_category_name => Origin of the term i.e SO
##          : $data_href          => Term hash {REF}
##          : $term               => Current term

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $data_category_name;
    my $data_href;
    my $term;

    my $tmpl = {
        data_category_name => {
            defined     => 1,
            required    => 1,
            store       => \$data_category_name,
            strict_type => 1,
        },
        data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$data_href,
            strict_type => 1,
        },
        term => { required => 1, store => \$term, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return 1 if exists $data_href->{$term};

    croak(  q{Could not find }
          . $data_category_name
          . q{ term from vcf in corresponding hash. Update hash to contain term}
          . qq{$COLON $SINGLE_QUOTE $term $SINGLE_QUOTE $NEWLINE} );
}

sub define_select_data_headers {

## Function : Defines arbitrary INFO fields headers based on information in select file
## Returns  : %select_data
## Arguments: None

    my %select_data;

    $select_data{select_file}{HGNC_symbol}{INFO} =
      q{##INFO=<ID=HGNC_symbol,Number=.,Type=String,Description="The HGNC gene symbol">};
    $select_data{select_file}{Ensembl_gene_id}{INFO} =
q{##INFO=<ID=Ensembl_gene_id,Number=.,Type=String,Description="Ensembl gene identifier">};
    $select_data{select_file}{OMIM_morbid}{INFO} =
q{##INFO=<ID=OMIM_morbid,Number=.,Type=String,Description="OMIM morbid ID associated with gene(s)">};
    $select_data{select_file}{Phenotypic_disease_model}{INFO} =
q{##INFO=<ID=Phenotypic_disease_model,Number=.,Type=String,Description="Known disease gene(s) phenotype inheritance model">};
    $select_data{select_file}{Reduced_penetrance}{INFO} =
q{##INFO=<ID=Reduced_penetrance,Number=.,Type=String,Description="Pathogenic gene which can exhibit reduced penetrance">};
    $select_data{select_file}{Disease_associated_transcript}{INFO} =
q{##INFO=<ID=Disease_associated_transcript,Number=.,Type=String,Description="Known pathogenic transcript(s) for gene">};
    $select_data{select_file}{Ensembl_transcript_to_refseq_transcript}{INFO} =
q{##INFO=<ID=Ensembl_transcript_to_refseq_transcript,Number=.,Type=String,Description="The link between ensembl transcript and refSeq transcript IDs">};
    $select_data{select_file}{Gene_description}{INFO} =
q{##INFO=<ID=Gene_description,Number=.,Type=String,Description="The HGNC gene description">};
    $select_data{select_file}{Genetic_disease_model}{INFO} =
q{##INFO=<ID=Genetic_disease_model,Number=.,Type=String,Description="Known disease gene(s) inheritance model">};
    $select_data{select_file}{no_hgnc_symbol}{INFO} =
q{##INFO=<ID=no_hgnc_symbol,Number=.,Type=String,Description="Clinically relevant genetic regions lacking a HGNC_symbol or Ensembl gene ">};
    return %select_data;
}

sub parse_consequence {

## Function : Parse consequence for most severe annotations
## Returns  : 1
## Arguments: $consequence_href         => Variant consequence {REF}
##          : $hgnc_map_href            => Hgnc map {REF}
##          : $most_severe_feature_href => Store most severe annotation per feature file(s)
##          : $most_severe_pli_href     => Store most severe pli per feature file(s)
##          : $per_gene                 => Only collect most severe transcript per gene
##          : $pli_score_href           => Pli score hash
##          : $select_data_href         => Select file data {REF}
##          : $vcf_record_href          => VCF record {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consequence_href;
    my $hgnc_map_href;
    my $most_severe_feature_href;
    my $most_severe_pli_href;
    my $pli_score_href;
    my $select_data_href;
    my $vcf_record_href;

    ## Default(s)
    my $per_gene;

    my $tmpl = {
        consequence_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$consequence_href,
            strict_type => 1,
        },
        hgnc_map_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$hgnc_map_href,
            strict_type => 1,
        },
        most_severe_feature_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$most_severe_feature_href,
            strict_type => 1,
        },
        most_severe_pli_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$most_severe_pli_href,
            strict_type => 1,
        },
        per_gene => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$per_gene,
            strict_type => 1,
        },
        pli_score_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pli_score_href,
            strict_type => 1,
        },
        select_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$select_data_href,
            strict_type => 1,
        },
        vcf_record_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vcf_record_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  HGNC_ID:
    for my $hgnc_id ( keys %{$consequence_href} ) {

        ## Unpack
        my $hgnc_symbol = $hgnc_map_href->{$hgnc_id};
        my $pli_score   = $pli_score_href->{$hgnc_symbol};

        ## For pli value and if current pli is more than stored
        set_most_severe_pli(
            {
                hgnc_id              => $hgnc_id,
                most_severe_pli_href => $most_severe_pli_href,
                pli_score            => $pli_score,
                select_data_href     => $select_data_href,
            }
        );

      ALLEL:
        for my $allele ( keys %{ $consequence_href->{$hgnc_id} } ) {

            ## Unpack
            my $most_severe_consequence =
              $consequence_href->{$hgnc_id}{$allele}{most_severe_consequence};
            my $most_severe_transcript =
              $consequence_href->{$hgnc_id}{$allele}{most_severe_transcript};

            add_most_severe_csq_to_feature(
                {
                    hgnc_id                  => $hgnc_id,
                    most_severe_consequence  => $most_severe_consequence,
                    most_severe_feature_href => $most_severe_feature_href,
                    most_severe_transcript   => $most_severe_transcript,
                    per_gene                 => $per_gene,
                    vcf_record_href          => $vcf_record_href,
                    select_data_href         => $select_data_href,
                }
            );
        }
    }
    return 1;
}

sub parse_vcf_format_line {

## Function : Parse VCF format line (#CHROM) and writes line to filehandle(s)
## Returns  : @vcf_format_columns
## Arguments: $filehandle       => The filehandle to write to
##          : $format_line      => VCF format line
##          : $selectfilehandle => The select filehandle to write to {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $format_line;
    my $selectfilehandle;

    my $tmpl = {
        filehandle  => { defined => 1, required => 1, store => \$filehandle, },
        format_line => {
            defined     => 1,
            required    => 1,
            store       => \$format_line,
            strict_type => 1,
        },
        selectfilehandle => { store => \$selectfilehandle, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Split vcf format/schema line
    my @vcf_format_columns = split $TAB, $format_line;

    ## Write #CHROM header line
    _write_to_file(
        {
            filehandle       => $filehandle,
            meta_data_line   => $format_line,
            selectfilehandle => $selectfilehandle,
        }
    );

    return @vcf_format_columns;
}

sub parse_vep_csq {

## Function : Parse VEP CSQ field
## Returns  :
## Arguments: $consequence_href             => Variant consequence {REF}
##          : $per_gene                     => Only collect most severe transcript per gene
##          : $pli_score_href               => Pli score hash
##          : $record_href                  => VCF record {REF}
##          : $select_data_href             => Select file data {REF}
##          : $vep_format_field_column_href => VEP format columns {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consequence_href;
    my $pli_score_href;
    my $record_href;
    my $select_data_href;
    my $vep_format_field_column_href;

    ## Default(s)
    my $per_gene;

    my $tmpl = {
        consequence_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$consequence_href,
            strict_type => 1,
        },
        per_gene => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$per_gene,
            strict_type => 1,
        },
        pli_score_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pli_score_href,
            strict_type => 1,
        },
        record_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$record_href,
            strict_type => 1,
        },
        select_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$select_data_href,
            strict_type => 1,
        },
        vep_format_field_column_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vep_format_field_column_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Vcfparser qw{
      add_transcript_to_feature_file
      parse_consequence
      parse_vep_csq_transcripts
      set_most_severe_ann_to_vcf_record };

    my @feature_type_keys = qw{ range select};

    ## Convert between hgnc_id and hgnc_symbol
    my %hgnc_map;

    ## Store most severe annotations
    my %most_severe_pli;
    my %most_severe_feature;

    ## Initilize pli score for feature keys
    @most_severe_pli{@feature_type_keys} = 0;

    if ( $record_href->{INFO_key_value}{CSQ} ) {

        ## Split into transcripts
        my @transcripts =
          split $COMMA, $record_href->{INFO_key_value}{CSQ};

        parse_vep_csq_transcripts(
            {
                consequence_href             => $consequence_href,
                hgnc_map_href                => \%hgnc_map,
                per_gene                     => $per_gene,
                select_data_href             => $select_data_href,
                transcripts_ref              => \@transcripts,
                vcf_record_href              => $record_href,
                vep_format_field_column_href => $vep_format_field_column_href,
            }
        );

        ## Parse consequence for most severe annotations
        parse_consequence(
            {
                consequence_href         => $consequence_href,
                hgnc_map_href            => \%hgnc_map,
                most_severe_feature_href => \%most_severe_feature,
                most_severe_pli_href     => \%most_severe_pli,
                per_gene                 => $per_gene,
                pli_score_href           => $pli_score_href,
                select_data_href         => $select_data_href,
                vcf_record_href          => $record_href,
            }
        );
    }

    set_most_severe_ann_to_vcf_record(
        {
            feature_type_keys_ref    => \@feature_type_keys,
            most_severe_feature_href => \%most_severe_feature,
            most_severe_pli_href     => \%most_severe_pli,
            vcf_record_href          => $record_href,
        }
    );
    return 1;
}

sub parse_vep_csq_transcripts {

## Function : Parse veo csq to get annotations per gene, allele and transcripts_id
## Returns  : 1
## Arguments: $consequence_href             => Variant consequence {REF}
##          : $hgnc_map_href                => Hgnc map {REF}
##          : $per_gene                     => Only collect most severe transcript per gene
##          : $transcripts_ref              => CSQ transcripts
##          : $select_data_href             => Select file data {REF}
##          : $vcf_record_href              => VCF record {REF}
##          : $vep_format_field_column_href => VEP format columns {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consequence_href;
    my $hgnc_map_href;
    my $select_data_href;
    my $transcripts_ref;
    my $vcf_record_href;
    my $vep_format_field_column_href;

    ## Default(s)
    my $per_gene;

    my $tmpl = {
        consequence_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$consequence_href,
            strict_type => 1,
        },
        hgnc_map_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$hgnc_map_href,
            strict_type => 1,
        },
        per_gene => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$per_gene,
            strict_type => 1,
        },
        select_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$select_data_href,
            strict_type => 1,
        },
        transcripts_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$transcripts_ref,
            strict_type => 1,
        },
        vcf_record_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vcf_record_href,
            strict_type => 1,
        },
        vep_format_field_column_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vep_format_field_column_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Vcf qw{ get_transcript_effects };

  TRANSCRIPT:
    foreach my $transcript ( @{$transcripts_ref} ) {

        ## Split transcript into VEP CSQ fields
        my @transcript_effects = split /[|]/sxm, $transcript;

        my %transcript_csq = get_transcript_effects(
            {
                transcript_effects_ref       => \@transcript_effects,
                vep_format_field_column_href => $vep_format_field_column_href
            }
        );

        ## If hgnc id
        if ( defined $transcript_csq{hgnc_id}
            and $transcript_csq{hgnc_id} )
        {

            ## VEP CSQ HGNC_ID format changed from grch37 "\d+" to grch38 "HGNC:\d+"
            $transcript_csq{hgnc_id} =~ s/HGNC://sxmg;

            ## Set symbol to hgnc map
            $hgnc_map_href->{ $transcript_csq{hgnc_id} } = $transcript_csq{hgnc_symbol};

            ## Parse the most severe consequence or prediction to gene
            parse_vep_csq_consequence(
                {
                    allele            => $transcript_csq{allele},
                    consequence_field => $transcript_csq{consequence_field},
                    consequence_href  => $consequence_href,
                    hgnc_id           => $transcript_csq{hgnc_id},
                    transcript        => $transcript,
                }
            );

            next TRANSCRIPT if ($per_gene);

            add_transcript_to_feature_file(
                {
                    hgnc_id          => $transcript_csq{hgnc_id},
                    select_data_href => $select_data_href,
                    transcripts_ref  => [$transcript],
                    vcf_record_href  => $vcf_record_href,
                }
            );
            next TRANSCRIPT;
        }

        ## Not part of a coding region
        add_transcript_to_feature_file(
            {
                transcripts_ref => [$transcript],
                vcf_record_href => $vcf_record_href,
            }
        );
    }
    return 1;
}

sub set_most_severe_ann_to_vcf_record {

## Function : Set most severe feature annotations to vcf record
## Returns  :
## Arguments: $feature_type_keys_ref    => Feature type keys
##          : $most_severe_feature_href => Store most severe annotation per feature file(s)
##          : $most_severe_pli_href     => Store most severe pli per feature file(s)
##          : $vcf_record_href          => VCF record {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $feature_type_keys_ref;
    my $most_severe_feature_href;
    my $most_severe_pli_href;
    my $vcf_record_href;

    my $tmpl = {
        feature_type_keys_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$feature_type_keys_ref,
            strict_type => 1,
        },
        most_severe_feature_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$most_severe_feature_href,
            strict_type => 1,
        },
        most_severe_pli_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$most_severe_pli_href,
            strict_type => 1,
        },
        vcf_record_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vcf_record_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  FEATURE_TYPE_KEY:
    foreach my $feature_type_key ( @{$feature_type_keys_ref} ) {

        my $vcf_key = join $UNDERSCORE,
          ( qw{INFO addition}, $feature_type_key, qw{ feature } );

        if ( $most_severe_pli_href->{$feature_type_key} ) {

            $vcf_record_href->{$vcf_key}{most_severe_pli} =
              $most_severe_pli_href->{$feature_type_key};
        }
        if ( exists $most_severe_feature_href->{$feature_type_key}
            and @{ $most_severe_feature_href->{$feature_type_key} } )
        {

            $vcf_record_href->{$vcf_key}{most_severe_consequence} = join $COMMA,
              @{ $most_severe_feature_href->{$feature_type_key} };
        }
    }

    return;
}

sub set_most_severe_pli {

## Function : Set the most severe pli keys for variant genes
## Returns  :
## Arguments: $hgnc_id              => Hgnc id
##          : $most_severe_pli_href => Store most severe pli {REF}
##          : $pli_score            => Pli score for gene
##          : $select_data_href     => Select file data {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $hgnc_id;
    my $most_severe_pli_href;
    my $pli_score;
    my $select_data_href;

    my $tmpl = {
        hgnc_id => {
            required    => 1,
            store       => \$hgnc_id,
            strict_type => 1,
        },
        most_severe_pli_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$most_severe_pli_href,
            strict_type => 1,
        },
        pli_score => {
            required    => 1,
            store       => \$pli_score,
            strict_type => 1,
        },
        select_data_href => {
            default     => {},
            store       => \$select_data_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not defined $pli_score );

    ## For pli value and if current pli is more than stored
    if ( $most_severe_pli_href->{range} < $pli_score ) {

        $most_severe_pli_href->{range} = $pli_score;

        return if ( not exists $select_data_href->{$hgnc_id} );

        $most_severe_pli_href->{select} = $pli_score;

    }
    return;
}

sub parse_vep_csq_consequence {

## Function : Parse the most severe consequence or prediction to gene
## Returns  :
## Arguments: $allele            => Allele
##          : $consequence_field => Consequence field from transcript
##          : $consequence_href  => Consequence hash {REF}
##          : $hgnc_id           => Hgnc id
##          : $transcript        => Transcript

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $allele;
    my $consequence_field;
    my $consequence_href;
    my $hgnc_id;
    my $transcript;

    my $tmpl = {
        allele => {
            defined     => 1,
            required    => 1,
            store       => \$allele,
            strict_type => 1,
        },
        consequence_field => {
            defined     => 1,
            required    => 1,
            store       => \$consequence_field,
            strict_type => 1,
        },
        consequence_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$consequence_href,
            strict_type => 1,
        },
        hgnc_id => {
            defined     => 1,
            required    => 1,
            store       => \$hgnc_id,
            strict_type => 1,
        },
        transcript => {
            defined     => 1,
            required    => 1,
            store       => \$transcript,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Vcf qw{ set_in_consequence_hash };

    ## Constants
    Readonly my $SO_RANK_INIT => 999;

    my %consequence_severity = %SO_CONSEQUENCE_SEVERITY;

    ## Alias set rank
    my $set_so_rank_ref = \$consequence_href->{$hgnc_id}{$allele}{rank};

    ## Initilialize arbitrary rank if unset
    if ( not ${$set_so_rank_ref} ) {
        ${$set_so_rank_ref} = $SO_RANK_INIT;
    }

    # Split consequence field for transcript
    my @consequences = split $AMPERSAND, $consequence_field;

  CONSEQUENCE:
    foreach my $consequence_term (@consequences) {

        check_data_terms(
            {
                data_href          => \%consequence_severity,
                term               => $consequence_term,
                data_category_name => q{SO},
            }
        );

        ## Build most severe consequence format
        my $most_severe_consequence =
          $hgnc_id . $COLON . $allele . $PIPE . $consequence_term;

        ## Unpack
        my $current_so_rank = $consequence_severity{$consequence_term}{rank};

        ## Map of what to set to consequence
        my %set_key = (
            most_severe_consequence => $most_severe_consequence,
            most_severe_transcript  => $transcript,
            rank                    => $current_so_rank,
        );

        ### Compare to previous set so consequence term
        ## Set new rank if $current_so_rank has lower rank than set_so_rank
        next CONSEQUENCE
          if ( $current_so_rank > ${$set_so_rank_ref} );

        ## Set most severe consequence key set in hash
        set_in_consequence_hash(
            {
                allele           => $allele,
                consequence_href => $consequence_href,
                hgnc_id          => $hgnc_id,
                set_key_href     => \%set_key,
            }
        );
    }
    return;
}

sub parse_vep_csq_schema {

## Function : Parse VEP CSQ format field and adds the format field index
## Returns  :
## Arguments: $meta_data_href               => Vcf meta data {REF}
##          : $parse_vep                    => Parse VEP output
##          : $vep_format_field_column_href => Vep format schema {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $meta_data_href;
    my $parse_vep;
    my $vep_format_field_column_href;

    my $tmpl = {
        meta_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$meta_data_href,
            strict_type => 1,
        },
        parse_vep => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$parse_vep,
            strict_type => 1,
        },
        vep_format_field_column_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vep_format_field_column_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Skip if CSQ line was not present in VCF header
    return if ( not exists $meta_data_href->{INFO}{CSQ} );

    ## Get "Format" of VEP feature fields within VEP CSQ header line
    ( my $vep_format ) = $meta_data_href->{INFO}{CSQ} =~ /Format:\s(\S+)"\>/sxm;

    my @vep_format_fields = split /[$PIPE]/sxm, $vep_format;

  FIELD:
    while ( my ( $field_index, $field ) = each @vep_format_fields ) {

        ## Set order of VEP features fields
        $vep_format_field_column_href->{$field} = $field_index;
    }

    ## Check that VEP FORMAT field schema match constant vcfparser CSQ map
  CSQ_FIELD:
    foreach my $csq_format_field_key ( keys %CSQ_FIELD_MAP ) {

        check_data_terms(
            {
                data_category_name => q{VEP_CSQ},
                data_href          => $vep_format_field_column_href,
                term               => $csq_format_field_key,
            }
        );
    }

    return if ( not $parse_vep );

    if (    exists $vep_format_field_column_href->{HGNC_ID}
        and exists $vep_format_field_column_href->{Consequence} )
    {

        $meta_data_href->{INFO}{most_severe_consequence} =
q{##INFO=<ID=most_severe_consequence,Number=.,Type=String,Description="Most severe genomic consequence.">};

    }
    return 1;
}

sub write_feature_file_csq {

## Function : Write vcf record to feature files
## Returns  : $info_field_counter
## Arguments: $filehandle         => The filehandle to write to
##          : $info_field_counter => Count the number of key_value pairs written
##          : $select_fh          => The select filehandle to write to {Optional}
##          : $vcf_record_href    => VCF record {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $select_fh;
    my $vcf_record_href;

    ## Default(s)
    my $info_field_counter;

    my $tmpl = {
        filehandle         => { defined => 1, required => 1, store => \$filehandle, },
        info_field_counter => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 0,
            defined     => 1,
            required    => 1,
            store       => \$info_field_counter,
            strict_type => 1,
        },
        select_fh       => { store => \$select_fh, },
        vcf_record_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vcf_record_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return
      if (  not exists $vcf_record_href->{INFO_key_value}{CSQ}
        and not $vcf_record_href->{INFO_key_value}{CSQ} );

    if ( $vcf_record_href->{range_transcripts} ) {

        print {$filehandle} q{CSQ} . $EQUALS . join $COMMA,
          @{ $vcf_record_href->{range_transcripts} };
    }
    if ( $vcf_record_href->{select_transcripts} ) {

        print {$select_fh} q{CSQ} . $EQUALS . join $COMMA,
          @{ $vcf_record_href->{select_transcripts} };
    }
    delete $vcf_record_href->{INFO_key_value}{CSQ};
    $info_field_counter++;

    return $info_field_counter;
}

sub write_info_addition_fields {

## Function : Write info addition fields of vcf record to feature files
## Returns  :
## Arguments: $filehandle         => The filehandle to write to
##          : $select_fh          => The select filehandle to write to {Optional}
##          : $vcf_record_href    => VCF record {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $select_fh;
    my $vcf_record_href;

    my $tmpl = {
        filehandle      => { defined => 1, required => 1, store => \$filehandle, },
        select_fh       => { store   => \$select_fh, },
        vcf_record_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vcf_record_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @info_addition_keys =
      qw( INFO_addition INFO_addition_range_feature INFO_addition_select_feature );

  INFO_ADDITION_KEY:
    foreach my $info_addition_key (@info_addition_keys) {

      INFO_KEY:
        foreach my $key ( keys %{ $vcf_record_href->{$info_addition_key} } ) {

            my $info_string =
              $SEMICOLON . $key . $EQUALS . $vcf_record_href->{$info_addition_key}{$key};
            print {$filehandle} $info_string;

            ## Never to select file for this info addition key
            next INFO_KEY if ( $info_addition_key eq q{INFO_addition_range_feature} );

            ## Write to select file
            next INFO_KEY if ( not $vcf_record_href->{select_transcripts} );

            print {$select_fh} $info_string;
        }
    }
    return;
}

sub write_info_field {

## Function : Write vcf record to feature files
## Returns  : $info_field_counter
## Arguments: $filehandle         => The filehandle to write to
##          : $info_field_counter => Count the number of key_value pairs written
##          : $select_fh          => The select filehandle to write to {Optional}
##          : $vcf_record_href    => VCF record {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $select_fh;
    my $vcf_record_href;

    ## Default(s)
    my $info_field_counter;

    my $tmpl = {
        filehandle         => { defined => 1, required => 1, store => \$filehandle, },
        info_field_counter => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 0,
            defined     => 1,
            required    => 1,
            store       => \$info_field_counter,
            strict_type => 1,
        },
        select_fh       => { store => \$select_fh, },
        vcf_record_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vcf_record_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %{ $vcf_record_href->{INFO_key_value} } ) {

        ## Build key value pair
        my $info_string = $SEMICOLON . $key;
        if ( not $info_field_counter ) {

            $info_string = $key;
        }

        if ( defined $value ) {

            $info_string .= $EQUALS . $value;
        }

        if ( $vcf_record_href->{select_transcripts} ) {

            print {$select_fh} $info_string;
        }
        print {$filehandle} $info_string;
        $info_field_counter++;
    }
    return $info_field_counter;
}

sub write_line_elements {

## Function : Writes metadata to filehandle specified by order in meta_data_sections.
## Returns  :
## Arguments: $filehandle        => The filehandle to write to
##          : $first_separator   => Separator to write before start index
##          : $last_index        => Last element index to write of line array
##          : $last_separator    => Separator to write after last index
##          : $line_elements_ref => Vcf record line
##          : $start_index       => Start element index to write of line array
##          : $select_fh         => The select filehandle to write to {Optional}
##          : $vcf_record_href   => Hash for vcf data {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $first_separator;
    my $last_index;
    my $last_separator;
    my $line_elements_ref;
    my $select_fh;
    my $start_index;
    my $vcf_record_href;

    my $tmpl = {
        filehandle      => { defined => 1, required => 1, store => \$filehandle, },
        first_separator => {
            allow       => [ $EMPTY_STR, $TAB ],
            default     => $EMPTY_STR,
            store       => \$first_separator,
            strict_type => 1,
        },
        last_index => {
            allow       => qr{ \A \d+ \z }xsm,
            store       => \$last_index,
            strict_type => 1,
        },
        last_separator => {
            allow       => [ $TAB, $NEWLINE ],
            default     => $NEWLINE,
            store       => \$last_separator,
            strict_type => 1,
        },
        line_elements_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$line_elements_ref,
            strict_type => 1,
        },
        select_fh   => { store => \$select_fh, },
        start_index => {
            allow       => qr{ \A \d+ \z }xsm,
            default     => 0,
            store       => \$start_index,
            strict_type => 1,
        },
        vcf_record_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vcf_record_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    print {$filehandle} $first_separator,
      join( $TAB, @{$line_elements_ref}[ $start_index .. $last_index ] ), $last_separator;

    return if ( not $vcf_record_href->{select_transcripts} );

    print {$select_fh} $first_separator,
      join( $TAB, @{$line_elements_ref}[ $start_index .. $last_index ] ),
      $last_separator;

    return;
}

sub write_meta_data {

## Function : Writes metadata to filehandle specified by order in meta_data_sections.
## Returns  :
## Arguments: $filehandle       => The filehandle to write to
##          : $meta_data_href   => Hash for meta_data {REF}
##          : $selectfilehandle => The select filehandle to write to {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $meta_data_href;
    my $selectfilehandle;

    my $tmpl = {
        filehandle     => { defined => 1, required => 1, store => \$filehandle, },
        meta_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$meta_data_href,
            strict_type => 1,
        },
        selectfilehandle => { store => \$selectfilehandle, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Determine order to print for standard vcf schema
    my @meta_data_vcf_schemas =
      qw{ fileformat ALT FILTER FORMAT INFO contig software other };

    ## Dispatch table of how to write meta data
    my %write_record = (
        contig => \&_write_vcf_schema,    # Written "as is"
        other  => \&_write_vcf_schema,
        vcf_id =>
          \&_write_vcf_schema_id_line, # All standard vcf_schema with vcf_id except contig
    );

  VCF_SCHEMA:
    foreach my $vcf_schema (@meta_data_vcf_schemas) {

        ## Meta data record exists
        next VCF_SCHEMA if ( not exists $meta_data_href->{$vcf_schema} );

        if ( exists $write_record{$vcf_schema} ) {

            $write_record{$vcf_schema}->(
                {
                    filehandle       => $filehandle,
                    meta_data_href   => $meta_data_href,
                    selectfilehandle => $selectfilehandle,
                    vcf_schema       => $vcf_schema,
                }
            );
            next VCF_SCHEMA;
        }

        $write_record{vcf_id}->(
            {
                filehandle       => $filehandle,
                meta_data_href   => $meta_data_href,
                selectfilehandle => $selectfilehandle,
                vcf_schema       => $vcf_schema,
            }
        );
    }
    return;
}

sub _write_to_file {

## Function : Writes metadata line to filehandle(s)
## Returns  :
## Arguments: $filehandle       => The filehandle to write to
##          : $selectfilehandle => The select filehandle to write to {Optional}
##          : $meta_data_line   => Meta data line

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $meta_data_line;
    my $selectfilehandle;

    my $tmpl = {
        filehandle     => { defined => 1, required => 1, store => \$filehandle, },
        meta_data_line => { defined => 1, required => 1, store => \$meta_data_line, },
        selectfilehandle => { store => \$selectfilehandle, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {$filehandle} $meta_data_line;

    if ( defined $selectfilehandle ) {

        say {$selectfilehandle} $meta_data_line;
    }
    return;
}

sub _write_vcf_schema {

## Function : Writes vcf schema records metadata to filehandle(s)
## Returns  :
## Arguments: $filehandle       => The filehandle to write to
##          : $meta_data_href   => Hash for meta_data {REF}
##          : $selectfilehandle => The select filehandle to write to {Optional}
##          : $vcf_schema       => Vcf schema

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $meta_data_href;
    my $selectfilehandle;
    my $vcf_schema;

    my $tmpl = {
        filehandle     => { defined => 1, required => 1, store => \$filehandle, },
        meta_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$meta_data_href,
            strict_type => 1,
        },
        selectfilehandle => { store   => \$selectfilehandle, },
        vcf_schema       => { defined => 1, required => 1, store => \$vcf_schema, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  RECORD_LINE:
    foreach my $header_line ( @{ $meta_data_href->{$vcf_schema}{$vcf_schema} } ) {

        _write_to_file(
            {
                filehandle       => $filehandle,
                meta_data_line   => $header_line,
                selectfilehandle => $selectfilehandle,
            }
        );
    }
    return;
}

sub _write_vcf_schema_id_line {

## Function : Writes vcf id metadata to filehandle(s)
## Returns  :
## Arguments: $filehandle       => The filehandle to write to
##          : $meta_data_href   => Hash for meta_data {REF}
##          : $selectfilehandle => The select filehandle to write to {Optional}
##          : $vcf_schema       => Vcf schema

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $meta_data_href;
    my $selectfilehandle;
    my $vcf_schema;

    my $tmpl = {
        filehandle     => { defined => 1, required => 1, store => \$filehandle, },
        meta_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$meta_data_href,
            strict_type => 1,
        },
        selectfilehandle => { store   => \$selectfilehandle, },
        vcf_schema       => { defined => 1, required => 1, store => \$vcf_schema, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  VCF_ID:
    foreach my $vcf_id ( sort keys %{ $meta_data_href->{$vcf_schema} } ) {

        my $header_line = $meta_data_href->{$vcf_schema}{$vcf_id};

        _write_to_file(
            {
                filehandle       => $filehandle,
                meta_data_line   => $header_line,
                selectfilehandle => $selectfilehandle,
            }
        );
    }

    ## Map of feature file type and corresponding filehandle
    my %feature_annotation = (
        range  => $filehandle,
        select => $selectfilehandle,
    );

    ## Add select specific annotations
  FEATURE_FILE_TYPE:
    while ( my ( $feature_file_type, $ANNOTATION_FH ) = each %feature_annotation ) {

      VCF_ID:
        foreach
          my $vcf_id ( sort keys %{ $meta_data_href->{$feature_file_type}{$vcf_schema} } )
        {

            my $header_line = $meta_data_href->{$feature_file_type}{$vcf_schema}{$vcf_id};
            if ( defined $ANNOTATION_FH ) {

                say {$ANNOTATION_FH} $header_line;
            }
        }
    }
    return;
}

1;
