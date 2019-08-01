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
use MIP::Constants qw{ $PIPE $SEMICOLON $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      add_feature_file_meta_data_to_vcf
      add_program_to_meta_data_header
      build_interval_tree
      define_select_data_headers
      parse_vep_csq_schema
      write_meta_data
    };
}

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
        push
          @{ $meta_data_href->{$file_key}{info}{$annotation} },
          $data_href->{present}{$annotation}{info};
    }
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
    push
      @{ $meta_data_href->{software}{$program_name} },
      qq{##Software=<ID=$program_name,Version=$vcfparser_version,Date=$current_date};
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

    ## Default(s)

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
    ( my $vep_format ) = $meta_data_href->{INFO}{CSQ}[0] =~ /Format:\s(\S+)"\>/sxm;

    my @vep_format_fields = split /[$PIPE]/sxm, $vep_format;

  FIELD:
    while ( my ( $field_index, $field ) = each @vep_format_fields ) {

        ## Set order of VEP features fields
        $vep_format_field_column_href->{$field} = $field_index;
    }

    return if ( not $parse_vep );

    if (    exists $vep_format_field_column_href->{HGNC_ID}
        and exists $vep_format_field_column_href->{Consequence} )
    {

        push
          @{ $meta_data_href->{info}{most_severe_consequence} },
q{##INFO=<ID=most_severe_consequence,Number=.,Type=String,Description="Most severe genomic consequence.">};
    }
    return 1;
}

sub write_meta_data {

## Function : Writes metadata to filehandle specified by order in meta_data_sections.
## Returns  :
## Arguments: $FILEHANDLE       => The filehandle to write to
##          : $meta_data_href   => Hash for meta_data {REF}
##          : $SELECTFILEHANDLE => The select filehandle to write to {Optional}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $meta_data_href;
    my $SELECTFILEHANDLE;

    my $tmpl = {
        FILEHANDLE     => { defined => 1, required => 1, store => \$FILEHANDLE, },
        meta_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$meta_data_href,
            strict_type => 1,
        },
        SELECTFILEHANDLE => { store => \$SELECTFILEHANDLE, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Determine order to print for standard records
    my @meta_data_records =
      qw{ fileformat ALT FILTER FORMAT INFO FIX_INFO contig Software };

  RECORD:
    foreach my $record (@meta_data_records) {

        ## Meta data record exists
        if ( $meta_data_href->{$record} ) {

            ## Should not be sorted, but printed "as is"
            if ( $record eq q{contig} ) {

              CONTIG_LINE:
                foreach my $contig_line ( @{ $meta_data_href->{$record}{$record} } ) {

                    say {$FILEHANDLE} $contig_line;

                    if ( defined $SELECTFILEHANDLE ) {

                        say {$SELECTFILEHANDLE} $contig_line;
                    }
                }
                ## Enable print of rest later
                delete $meta_data_href->{$record};
            }
            else {

                foreach my $line ( sort( keys %{ $meta_data_href->{$record} } ) ) {

                    say {$FILEHANDLE} @{ $meta_data_href->{$record}{$line} };

                    if ( defined $SELECTFILEHANDLE ) {

                        say {$SELECTFILEHANDLE} @{ $meta_data_href->{$record}{$line} };
                    }
                }
                if ( defined $SELECTFILEHANDLE ) {

                    foreach
                      my $line ( sort( keys %{ $meta_data_href->{select}{$record} } ) )
                    {

                        say {$SELECTFILEHANDLE}
                          @{ $meta_data_href->{select}{$record}{$line} };
                    }
                }

                ## Range
                foreach my $line ( sort( keys %{ $meta_data_href->{range}{$record} } ) ) {

                    say {$FILEHANDLE} @{ $meta_data_href->{range}{$record}{$line} };
                }

                ## Enable print of rest later
                delete( $meta_data_href->{$record} );

                if ( $meta_data_href->{select}{$record} ) {

                    delete( $meta_data_href->{select}{$record} );
                }
                if ( $meta_data_href->{range}{$record} ) {

                    delete( $meta_data_href->{range}{$record} );
                }
            }
        }
    }
    ## Writing "the rest"
    for my $keys ( keys %{$meta_data_href} ) {

        for my $second_key ( keys %{ $meta_data_href->{$keys} } ) {

            if ( ref $meta_data_href->{$keys}{$second_key} eq q{HASH} ) {

                for my $line ( sort( keys %{ $meta_data_href->{$keys}{$second_key} } ) ) {

                    say {$FILEHANDLE} @{ $meta_data_href->{$keys}{$second_key}{$line} };

                    if ( defined $SELECTFILEHANDLE ) {

                        say {$SELECTFILEHANDLE}
                          @{ $meta_data_href->{$keys}{$second_key}{$line} };
                    }
                    delete $meta_data_href->{$keys}{$second_key}{$line};
                }
            }
            elsif ( ref $meta_data_href->{$keys}{$second_key} eq q{ARRAY} ) {

                foreach my $element ( @{ $meta_data_href->{$keys}{$second_key} } ) {

                    say $element;

                    if ( defined $SELECTFILEHANDLE ) {

                        say {$SELECTFILEHANDLE} $element;
                    }
                }
                delete $meta_data_href->{$keys}{$second_key};
            }
            else {

                say {$FILEHANDLE} @{ $meta_data_href->{$keys}{$second_key} };

                if ( defined $SELECTFILEHANDLE ) {

                    say {$SELECTFILEHANDLE} @{ $meta_data_href->{$keys}{$second_key} };
                }
                delete $meta_data_href->{$keys}{$second_key};
            }
        }
    }
    return;
}

1;
