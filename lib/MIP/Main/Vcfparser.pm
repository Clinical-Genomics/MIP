package MIP::Main::Vcfparser;

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
use MIP::Constants qw{ $COLON $LOG_NAME $NEWLINE %SO_CONSEQUENCE_SEVERITY $TAB };
use MIP::File::Format::Feature_file qw{ read_feature_file };
use MIP::File::Format::Pli qw{ load_pli_file };
use MIP::Log::MIP_log4perl qw{ retrieve_log };
use MIP::Vcfparser qw{ define_select_data_headers };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = q{1.2.18};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ mip_vcfparser };
}

sub mip_vcfparser {

## Function : Split vcf file into selecet (clinical) and research variants. Enables custom annotation via range and feature files
## Returns  :
## Arguments: $padding                               => Padding distance
##          : $parse_vep                             => Parse VEP output
##          : $per_gene                              => Only collect most severe transcript per gene
##          : $pli_values_file_path                  => Pli values file_path
##          : $range_feature_annotation_columns_ref  => Range feature columns {REF}
##          : $range_feature_file                    => Range feature file
##          : $select_feature_annotation_columns_ref => Select feature columns {REF}
##          : $select_feature_file                   => Select feature file
##          : $select_feature_matching_column        => Select feature matching column
##          : $select_outfile_path                   => Select file path
##          : $vcf_in_fh                             => VCF in filehandle
##          : $write_software_tag                    => Write software tag to vcf header switch

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $per_gene;
    my $pli_values_file_path;
    my $range_feature_annotation_columns_ref;
    my $range_feature_file;
    my $select_feature_annotation_columns_ref;
    my $select_feature_matching_column;
    my $select_outfile_path;

    ## Default(s)
    my $padding;
    my $parse_vep;
    my $select_feature_file;
    my $vcf_in_fh;
    my $write_software_tag;

    my $tmpl = {
        padding => {
            defined     => 1,
            required    => 1,
            store       => \$padding,
            strict_type => 1,
        },
        parse_vep => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$parse_vep,
            strict_type => 1,
        },
        per_gene => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$per_gene,
            strict_type => 1,
        },
        pli_values_file_path => {
            store       => \$pli_values_file_path,
            strict_type => 1,
        },
        range_feature_annotation_columns_ref => {
            default     => [],
            store       => \$range_feature_annotation_columns_ref,
            strict_type => 1,
        },
        range_feature_file => {
            default     => 0,
            store       => \$range_feature_file,
            strict_type => 1,
        },
        select_feature_annotation_columns_ref => {
            default     => [],
            store       => \$select_feature_annotation_columns_ref,
            strict_type => 1,
        },
        select_feature_file => {
            default     => 0,
            store       => \$select_feature_file,
            strict_type => 1,
        },
        select_feature_matching_column => {
            store       => \$select_feature_matching_column,
            strict_type => 1,
        },
        select_outfile_path => { store => \$select_outfile_path, strict_type => 1, },
        vcf_in_fh          => { defined => 1, required => 1, store => \$vcf_in_fh, },
        write_software_tag => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$write_software_tag,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Constants
    my %consequence_severity = %SO_CONSEQUENCE_SEVERITY;

    my ( %meta_data, %range_data, %tree, %pli_score );

    ## Retrieve logger object now that log file has been set
    my $log = Log::Log4perl->get_logger(q{Vcfparser});

    if ($pli_values_file_path) {

        $log->info( q{Loading pli value file: } . $pli_values_file_path );
        load_pli_file(
            {
                infile_path    => $pli_values_file_path,
                log            => $log,
                pli_score_href => \%pli_score,
            }
        );
        $log->info(q{Loading pli value file: Done});

        ## Add pli header line to VCF meta data HASH
        $meta_data{INFO}{most_severe_pli} =
q{##INFO=<ID=most_severe_pli,Number=1,Type=Float,Description="Most severe pli score.">};
    }

    my %select_data = define_select_data_headers();

    ## Range feature file
    read_feature_file(
        {
            feature_columns_ref => $range_feature_annotation_columns_ref,
            feature_data_href   => \%range_data,
            log                 => $log,
            feature_file_path   => $range_feature_file,
            padding             => $padding,
            feature_file_type   => q{range_feature},
            tree_href           => \%tree,
        }
    );

    ## Select feature file
    read_feature_file(
        {
            feature_columns_ref     => $select_feature_annotation_columns_ref,
            feature_data_href       => \%select_data,
            log                     => $log,
            feature_file_path       => $select_feature_file,
            padding                 => $padding,
            feature_file_type       => q{select_feature},
            feature_matching_column => $select_feature_matching_column,
            tree_href               => \%tree,
        }
    );

    read_infile_vcf(
        {
            consequence_severity_href            => \%consequence_severity,
            meta_data_href                       => \%meta_data,
            parse_vep                            => $parse_vep,
            per_gene                             => $per_gene,
            pli_score_href                       => \%pli_score,
            range_data_href                      => \%range_data,
            range_feature_annotation_columns_ref => $range_feature_annotation_columns_ref,
            select_data_href                     => \%select_data,
            select_feature_annotation_columns_ref =>
              $select_feature_annotation_columns_ref,
            select_feature_file => $select_feature_file,
            select_outfile_path => $select_outfile_path,
            tree_href           => \%tree,
            vcfparser_version   => $MIP::Main::Vcfparser::VERSION,
            vcf_in_fh           => $vcf_in_fh,
            write_software_tag  => $write_software_tag,
        }
    );

    return;
}

####################
####Sub routines####
####################

sub read_infile_vcf {

## Function : Reads infile in vcf format, adds and parses annotations as well as split transcripts into select subset file
## Returns  :
## Arguments: $consequence_severity_href             => Consequence severity for SO-terms {REF}
##          : $meta_data_href                        => Vcf meta data {REF}
##          : $parse_vep                             => Parse VEP output
##          : $per_gene                              => Only collect most severe transcript per gene
##          : $pli_score_href                        => Pli score hash
##          : $range_data_href                       => Range file data {REF}
##          : $range_feature_annotation_columns_ref  => Range feature columns {REF}
##          : $select_data_href                      => Select file data {REF}
##          : $select_feature_annotation_columns_ref => Select feature columns {REF}
##          : $select_feature_file                   => Select feature file
##          : $select_outfile_path                   => Select file path
##          : $tree_href                             => Interval tree hash {REF}
##          : $vcfparser_version                     => Vcfparser version
##          : $write_software_tag                    => Write software tag to vcf header switch
##          : $vcf_in_fh                             => The filehandle to read from to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consequence_severity_href;
    my $meta_data_href;
    my $per_gene;
    my $pli_score_href;
    my $range_data_href;
    my $range_feature_annotation_columns_ref;
    my $select_data_href;
    my $select_feature_annotation_columns_ref;
    my $select_outfile_path;
    my $tree_href;
    my $vcfparser_version;
    my $vcf_in_fh;

    ## Default(s)
    my $select_feature_file;
    my $parse_vep;
    my $write_software_tag;

    my $tmpl = {
        consequence_severity_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$consequence_severity_href,
            strict_type => 1,
        },
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
        range_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$range_data_href,
            strict_type => 1,
        },
        range_feature_annotation_columns_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$range_feature_annotation_columns_ref,
            strict_type => 1,
        },
        select_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$select_data_href,
            strict_type => 1,
        },
        select_feature_annotation_columns_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$select_feature_annotation_columns_ref,
            strict_type => 1,
        },
        select_feature_file => {
            default     => 0,
            store       => \$select_feature_file,
            strict_type => 1,
        },
        select_outfile_path => { store => \$select_outfile_path, strict_type => 1, },
        tree_href           => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$tree_href,
            strict_type => 1,
        },
        vcfparser_version => {
            defined     => 1,
            required    => 1,
            store       => \$vcfparser_version,
            strict_type => 1,
        },
        vcf_in_fh          => { defined => 1, required => 1, store => \$vcf_in_fh, },
        write_software_tag => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$write_software_tag,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Feature_file qw{ tree_annotations };
    use MIP::File::Format::Vcf qw{
      check_vcf_variant_line parse_vcf_header
      set_info_key_pairs_in_vcf_record
      set_line_elements_in_vcf_record };
    use MIP::Vcfparser qw{
      add_feature_file_meta_data_to_vcf
      add_program_to_meta_data_header
      parse_vcf_format_line
      parse_vep_csq
      parse_vep_csq_schema
      write_feature_file_csq
      write_info_addition_fields
      write_info_field
      write_line_elements
      write_meta_data
    };

    ## Constants
    Readonly my $FORMAT_COLUMN_INDEX => 8;

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger(q{Vcfparser});

    ## Map the VEP CSQ format header line
    my %vep_format_field_column;

    ## Store vcf header "#CHROM" line
    my @vcf_format_columns;

    my $select_fh = _get_select_filehandle(
        {
            select_feature_file => $select_feature_file,
            select_outfile_path => $select_outfile_path,
        }
    );

  LINE:
    while (<$vcf_in_fh>) {

        chomp;

        ## Unpack line
        my $line = $_;

        ## Skip blank lines
        next LINE if ( $line =~ /^\s+$/sxm );

        ## Header meta data
        if ( $line =~ /\A [#]{2}/sxm ) {

            parse_vcf_header(
                {
                    meta_data_href   => $meta_data_href,
                    meta_data_string => $line,
                }
            );

            ## Parse VEP CSQ format field and adds the format field index
            parse_vep_csq_schema(
                {
                    meta_data_href               => $meta_data_href,
                    parse_vep                    => $parse_vep,
                    vep_format_field_column_href => \%vep_format_field_column,
                }
            );
            next;
        }
        if ( $line =~ /\A [#]{1}CHROM/sxm ) {

            add_feature_file_meta_data_to_vcf(
                {
                    data_href => $range_data_href,
                    feature_annotation_columns_ref =>
                      $range_feature_annotation_columns_ref,
                    file_key       => q{Range},
                    meta_data_href => $meta_data_href,
                }
            );
            add_feature_file_meta_data_to_vcf(
                {
                    data_href => $select_data_href,
                    feature_annotation_columns_ref =>
                      $select_feature_annotation_columns_ref,
                    file_key       => q{Select},
                    meta_data_href => $meta_data_href,
                }
            );

            add_program_to_meta_data_header(
                {
                    add_software_tag  => $write_software_tag,
                    meta_data_href    => $meta_data_href,
                    vcfparser_version => $vcfparser_version,
                }
            );

            write_meta_data(
                {
                    filehandle       => *STDOUT,
                    meta_data_href   => $meta_data_href,
                    selectfilehandle => $select_fh,
                }
            );

            @vcf_format_columns = parse_vcf_format_line(
                {
                    filehandle       => *STDOUT,
                    format_line      => $line,
                    selectfilehandle => $select_fh,
                }
            );
            next;
        }

        ## Variant line
        my %consequence;
        my %vcf_record;

        ## Loads vcf file elements
        my @line_elements = split $TAB, $line;

        check_vcf_variant_line(
            {
                input_line_number         => $INPUT_LINE_NUMBER,
                log                       => $log,
                variant_line              => $line,
                variant_line_elements_ref => \@line_elements,
            }
        );

        ## Adds variant line elements to vcf_record hash
        set_line_elements_in_vcf_record(
            {
                line_elements_ref      => \@line_elements,
                vcf_record_href        => \%vcf_record,
                vcf_format_columns_ref => \@vcf_format_columns,
            }
        );

        ## Adds INFO key value pairs to vcf_record hash
        set_info_key_pairs_in_vcf_record( { vcf_record_href => \%vcf_record, } );

        ## Checks if an interval tree exists (per chr) and
        ## collects features from input array and adds annotations to line
        ## noid_region is only for selectfile since all variants are passed to research file
        my %noid_region = tree_annotations(
            {
                alt_allele_field  => $vcf_record{ALT},
                contig            => $vcf_record{q{#CHROM}},
                data_href         => $select_data_href,
                feature_file_type => q{select_feature},
                record_href       => \%vcf_record,
                ref_allele        => $vcf_record{REF},
                start             => $vcf_record{POS},
                tree_href         => $tree_href,
            }
        );

        tree_annotations(
            {
                alt_allele_field  => $vcf_record{ALT},
                contig            => $vcf_record{q{#CHROM}},
                data_href         => $range_data_href,
                feature_file_type => q{range_feature},
                record_href       => \%vcf_record,
                ref_allele        => $vcf_record{REF},
                start             => $vcf_record{POS},
                tree_href         => $tree_href,
            }
        );

        if ($parse_vep) {

            parse_vep_csq(
                {
                    consequence_href             => \%consequence,
                    per_gene                     => $per_gene,
                    pli_score_href               => $pli_score_href,
                    record_href                  => \%vcf_record,
                    select_data_href             => $select_data_href,
                    vep_format_field_column_href => \%vep_format_field_column,
                }
            );
        }

        ### Writing vcf record to files

        ## Get parameters to know how much of line to write
        my ( $last_index, $last_separator ) = _get_line_parameters(
            {
                line_elements_ref => \@line_elements,
                parse_vep         => $parse_vep,
            }
        );

        ## Write until INFO field or end of line
        write_line_elements(
            {
                filehandle        => *STDOUT,
                last_index        => $last_index,
                last_separator    => $last_separator,
                line_elements_ref => \@line_elements,
                select_fh         => $select_fh,
                start_index       => 0,
                vcf_record_href   => \%vcf_record,
            }
        );

        if ($parse_vep) {

            my $info_field_counter = write_feature_file_csq(
                {
                    filehandle         => *STDOUT,
                    info_field_counter => 0,
                    select_fh          => $select_fh,
                    vcf_record_href    => \%vcf_record,
                }
            );

            write_info_field(
                {
                    filehandle         => *STDOUT,
                    info_field_counter => $info_field_counter,
                    select_fh          => $select_fh,
                    vcf_record_href    => \%vcf_record,
                }
            );

            write_info_addition_fields(
                {
                    filehandle      => *STDOUT,
                    select_fh       => $select_fh,
                    vcf_record_href => \%vcf_record,
                }
            );

            ## After INFO to the final column
            write_line_elements(
                {
                    filehandle        => *STDOUT,
                    first_separator   => $TAB,
                    last_index        => $#line_elements,
                    last_separator    => $NEWLINE,
                    line_elements_ref => \@line_elements,
                    select_fh         => $select_fh,
                    start_index       => $FORMAT_COLUMN_INDEX,
                    vcf_record_href   => \%vcf_record,
                }
            );

        }
    }
    if ($select_feature_file) {

        close $select_fh;
    }
    $log->info( q{Finished Processing VCF} . $NEWLINE );
    return;
}

sub _get_line_parameters {

## Function : Get the last index and last seperator value depending on parse vep boolean
## Returns  : $last_index, $last_separator
## Arguments: $line_elements_ref => Line elements array {REF}
##          : $parse_vep         => Scalar description

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $line_elements_ref;

    ## Default(s)
    my $parse_vep;

    my $tmpl = {
        line_elements_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$line_elements_ref,
            strict_type => 1,
        },
        parse_vep => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$parse_vep,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Constants
    Readonly my $FILTER_COLUMN_INDEX => 6;

    ## If we do not need to split the CSQ field we can print the entire vcf
    ## record line
    my $last_index     = $#{$line_elements_ref};
    my $last_separator = $NEWLINE;

    return ( $last_index, $last_separator ) if ( not $parse_vep );

    ## Set parameters to write until INFO field if we need to split CSQ field
    $last_index     = $FILTER_COLUMN_INDEX;
    $last_separator = $TAB;

    return $last_index, $last_separator;
}

sub _get_select_filehandle {

## Function : Returns filehandle for select file
## Returns  : $select_fh or undef
## Arguments: $select_feature_file => Select feature file
##          : $select_outfile_path => Select outfile path (VCF)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $select_feature_file;
    my $select_outfile_path;

    my $tmpl = {
        select_feature_file => {
            required    => 1,
            store       => \$select_feature_file,
            strict_type => 1,
        },
        select_outfile_path => {
            required    => 1,
            store       => \$select_outfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## If we do not have a select file return undef
    return if ( not $select_feature_file );

    my $log = retrieve_log( { log_name => $LOG_NAME, } );

    open my $select_fh, q{>},
      $select_outfile_path
      or $log->logdie( q{Cannot open } . $select_outfile_path . $COLON . $OS_ERROR,
        $NEWLINE );

    return $select_fh;
}

1;
