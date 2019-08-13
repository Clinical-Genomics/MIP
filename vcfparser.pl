#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ basename };
use File::Spec::Functions qw{ catdir catfile devnull };
use FindBin qw{ $Bin };
use Getopt::Long;
use IO::File;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings qw{ FATAL utf8 };

$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

## CPANM
use autodie qw{ open close :all };
use Modern::Perl qw{ 2017 };
use Readonly;
use Set::IntervalTree;

## MIPs lib/
use lib catdir( $Bin, q{lib} );
use MIP::Check::Modules qw{ check_perl_modules };
use MIP::Constants
  qw{ %ANALYSIS $COLON $COMMA $NEWLINE %SO_CONSEQUENCE_SEVERITY $SPACE $TAB };
use MIP::File::Format::Feature_file qw{ read_feature_file };
use MIP::File::Format::Pli qw{ load_pli_file };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Script::Utils qw{ help };
use MIP::Vcfparser qw{ define_select_data_headers };

our $USAGE = build_usage( {} );

BEGIN {

    require MIP::Check::Modules;

    my @modules =
      qw{ Modern::Perl autodie Set::IntervalTree Log::Log4perl MIP::Log::MIP_log4perl };

    ## Evaluate that all modules required are installed
    check_perl_modules(
        {
            modules_ref  => \@modules,
            program_name => $PROGRAM_NAME,
        }
    );
}

## Constants
Readonly my $ANNOTATION_DISTANCE => $ANALYSIS{ANNOTATION_DISTANCE};
my %consequence_severity = %SO_CONSEQUENCE_SEVERITY;

my ( $infile, $pli_values_file_path, $range_feature_file, $select_feature_file,
    $select_feature_matching_column,
    $select_outfile, );

## Scalar parameters with defaults
my ( $write_software_tag, $padding, $log_file ) =
  ( 1, $ANNOTATION_DISTANCE, catfile( cwd(), q{vcfparser.log} ) );

## Boolean
my ( $parse_vep, $per_gene );

my ( @range_feature_annotation_columns, @select_feature_annotation_columns );
my ( %meta_data, %range_data, %tree, %pli_score );

my $VERSION = q{1.2.16};

## Enables cmd "vcfparser.pl" to print usage help
if ( not @ARGV ) {

    help(
        {
            USAGE     => $USAGE,
            exit_code => 0,
        }
    );
}
elsif ( defined $ARGV and $ARGV[0] !~ /^-/sxm ) {
    ## Collect potential infile - otherwise read from STDIN

    $infile = $ARGV[0];
}

### User Options
GetOptions(
    q{pvep|parse_vep}          => \$parse_vep,
    q{rf|range_feature_file:s} => \$range_feature_file,
    q{rf_ac|range_feature_annotation_columns:s} =>
      \@range_feature_annotation_columns,    #Comma separated list
    q{sf|select_feature_file:s}               => \$select_feature_file,
    q{sf_mc|select_feature_matching_column:n} => \$select_feature_matching_column,
    q{sf_ac|select_feature_annotation_columns:s} =>
      \@select_feature_annotation_columns,    #Comma separated list
    q{sof|select_outfile:s}     => \$select_outfile,
    q{wst|write_software_tag:n} => \$write_software_tag,
    q{pad|padding:n}            => \$padding,
    q{peg|per_gene}             => \$per_gene,
    q{pli|pli_values_file:s}    => \$pli_values_file_path,
    q{l|log_file:s}             => \$log_file,
    ## Display help text
    q{h|help} => sub { say {*STDOUT} $USAGE; exit; },
    ## Display version number
    q{v|version} => sub {
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION, $NEWLINE;
        exit;
    },
  )
  or help(
    {
        USAGE     => $USAGE,
        exit_code => 1,
    }
  );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $log_file,
        log_name  => q{Vcfparser},
    }
);

## Basic flag option check
if ( not @range_feature_annotation_columns and $range_feature_file ) {

    $log->info($USAGE);
    $log->fatal(
        q{Need to specify which feature column(s) to use with range feature file: }
          . $range_feature_file
          . q{ when annotating variants by using flag -rf_ac},
        $NEWLINE
    );
    exit 1;
}
if ( not $select_feature_matching_column and $select_feature_file ) {

    $log->info($USAGE);
    $log->fatal(
        q{Need to specify which feature column to use with select feature file: }
          . $select_feature_file
          . q{ when selecting variants by using flag -sf_mc},
        $NEWLINE
    );
    exit 1;
}
if ( not $select_outfile and $select_feature_file ) {

    $log->info($USAGE);
    $log->fatal(
q{Need to specify which a select outfile to use when selecting variants by using flag -sof},
        $NEWLINE
    );
    exit 1;
}

## Enables comma separated annotation columns on cmd
@range_feature_annotation_columns =
  split $COMMA, join $COMMA, @range_feature_annotation_columns;

## Enables comma separated annotation columns on cmd
@select_feature_annotation_columns =
  split $COMMA, join $COMMA, @select_feature_annotation_columns;

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

############
####MAIN####
############

my %select_data = define_select_data_headers();

if ($range_feature_file) {

    read_feature_file(
        {
            feature_columns_ref => \@range_feature_annotation_columns,
            feature_data_href   => \%range_data,
            log                 => $log,
            feature_file_path   => $range_feature_file,
            padding             => $padding,
            feature_file_type   => q{range_feature},
            tree_href           => \%tree,
        }
    );
}

if ($select_feature_file) {

    read_feature_file(
        {
            feature_columns_ref     => \@select_feature_annotation_columns,
            feature_data_href       => \%select_data,
            log                     => $log,
            feature_file_path       => $select_feature_file,
            padding                 => $padding,
            feature_file_type       => q{select_feature},
            feature_matching_column => $select_feature_matching_column,
            tree_href               => \%tree,
        }
    );
}

read_infile_vcf(
    {
        consequence_severity_href             => \%consequence_severity,
        meta_data_href                        => \%meta_data,
        parse_vep                             => $parse_vep,
        per_gene                              => $per_gene,
        pli_score_href                        => \%pli_score,
        range_data_href                       => \%range_data,
        range_feature_annotation_columns_ref  => \@range_feature_annotation_columns,
        select_data_href                      => \%select_data,
        select_feature_annotation_columns_ref => \@select_feature_annotation_columns,
        select_feature_file                   => $select_feature_file,
        select_outfile_path                   => $select_outfile,
        tree_href                             => \%tree,
        vcfparser_version                     => $VERSION,
        write_software_tag                    => $write_software_tag,
    }
);

####################
####Sub routines####
####################

sub build_usage {

## Function : Build the USAGE instructions
## Returns  :
## Arguments: $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options] infile.vcf [OPTIONS] > outfile.vcf
   -pvep/--parse_vep Parse VEP transcript specific entries (Supply flag to enable)
   -rf/--range_feature_file (tsv)
   -rf_ac/--range_feature_annotation_columns
   -sf/--select_feature_file (tsv)
   -sf_mc/--select_feature_matching_column
   -sf_ac/--select_feature_annotation_columns
   -sof/--select_outfile (vcf)
   -pad/--padding (Default: "5000" nucleotides)
   -peg/--per_gene Output most severe consequence transcript (Supply flag to enable)
   -pli/--pli_values_file Pli value file path
   -wst/--write_software_tag (Default: "1")
   -l/--log_file Log file (Default: "vcfparser.log")
   -h/--help Display this help message
   -v/--version Display version
END_USAGE
}

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
      parse_vep_csq_schema
      write_meta_data
    };

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger(q{Vcfparser});

    ## Create anonymous filehandle for select file
    my $SELECT_FH = IO::Handle->new();

    ## Map the VEP CSQ format header line
    my %vep_format_field_column;

    ## Store vcf header "#CHROM" line
    my @vcf_format_columns;

    if ($select_feature_file) {

        open $SELECT_FH, q{>},
          $select_outfile_path
          or $log->logdie( q{Cannot open } . $select_outfile_path . $COLON . $OS_ERROR,
            $NEWLINE );
    }
    else {
        ## If we do not have a select file undef the filehandle
        $SELECT_FH = undef;
    }

  LINE:
    while (<>) {

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
                    FILEHANDLE       => *STDOUT,
                    meta_data_href   => $meta_data_href,
                    SELECTFILEHANDLE => $SELECT_FH,
                }
            );

            @vcf_format_columns = parse_vcf_format_line(
                {
                    FILEHANDLE       => *STDOUT,
                    format_line      => $line,
                    SELECTFILEHANDLE => $SELECT_FH,
                }
            );
            next;
        }

        ## Variant line
        my %consequence;
        my %noid_region;
        my %record;
        my $selected_variant_line;
        my $variant_line;

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

        ## Adds variant line elements to record hash
        set_line_elements_in_vcf_record(
            {
                line_elements_ref      => \@line_elements,
                vcf_record_href        => \%record,
                vcf_format_columns_ref => \@vcf_format_columns,
            }
        );

        ## Adds INFO key value pairs to record hash
        set_info_key_pairs_in_vcf_record( { vcf_record_href => \%record, } );

        ## Checks if an interval tree exists (per chr) and
        ## collects features from input array and adds annotations to line
        ## noid_region is only for selectfile since all variants are passed to research file
        %noid_region = tree_annotations(
            {
                alt_allele_field  => $record{ALT},
                contig            => $record{q{#CHROM}},
                data_href         => $select_data_href,
                feature_file_type => q{select_feature},
                record_href       => \%record,
                ref_allele        => $record{REF},
                start             => $record{POS},
                tree_href         => $tree_href,
            }
        );

        tree_annotations(
            {
                alt_allele_field  => $record{ALT},
                contig            => $record{q{#CHROM}},
                data_href         => $range_data_href,
                feature_file_type => q{range_feature},
                record_href       => \%record,
                ref_allele        => $record{REF},
                start             => $record{POS},
                tree_href         => $tree_href,
            }
        );

        if ($parse_vep) {

            parse_vep_csq(
                {
                    consequence_href             => \%consequence,
                    consequence_severity_href    => $consequence_severity_href,
                    per_gene                     => $per_gene,
                    pli_score_href               => $pli_score_href,
                    record_href                  => \%record,
                    select_data_href             => $select_data_href,
                    vep_format_field_column_href => \%vep_format_field_column,
                }
            );
        }

        for (
            my $line_elements_counter = 0 ;
            $line_elements_counter < scalar(@line_elements) ;
            $line_elements_counter++
          )
        {    #Add until INFO field

            if ( $line_elements_counter < 7 ) {    #Save fields until INFO field

                if ( $record{select_transcripts} ) {

                    print $SELECT_FH $record{ $vcf_format_columns[$line_elements_counter]
                      }
                      . "\t";
                }
                print STDOUT $record{ $vcf_format_columns[$line_elements_counter] }
                  . "\t";
            }

            if ( $line_elements_counter == 7 ) {

                if ( !$parse_vep ) {

                    if ( $record{select_transcripts} ) {

                        print $SELECT_FH
                          $record{ $vcf_format_columns[$line_elements_counter] };
                    }
                    print STDOUT $record{ $vcf_format_columns[$line_elements_counter] };
                }
                else {

                    my $counter = 0;
                    foreach my $key ( keys %{ $record{INFO_key_value} } ) {

                        if ( !$counter ) {

                            if ( defined( $record{INFO_key_value}{$key} ) ) {

                                if ( $key eq "CSQ" ) {

                                    if ( $record{range_transcripts} ) {

                                        print STDOUT $key . "="
                                          . join( ",", @{ $record{range_transcripts} } );
                                    }
                                    if ( $record{select_transcripts} ) {

                                        print $SELECT_FH $key . "="
                                          . join( ",", @{ $record{select_transcripts} } );
                                    }
                                }
                                else {

                                    if ( $record{select_transcripts} ) {

                                        print $SELECT_FH $key . "="
                                          . $record{INFO_key_value}{$key};
                                    }
                                    print STDOUT $key . "="
                                      . $record{INFO_key_value}{$key};
                                }
                            }
                            else {

                                if ( $record{select_transcripts} ) {

                                    print $SELECT_FH $key;
                                }
                                print STDOUT $key;
                            }
                        }
                        else {

                            if ( defined( $record{INFO_key_value}{$key} ) ) {

                                if ( $key eq "CSQ" ) {

                                    if ( $record{range_transcripts} ) {

                                        print STDOUT ";" . $key . "="
                                          . join( ",", @{ $record{range_transcripts} } );
                                    }
                                    if ( $record{select_transcripts} ) {

                                        print $SELECT_FH ";" . $key . "="
                                          . join( ",", @{ $record{select_transcripts} } );
                                    }
                                }
                                else {

                                    if ( $record{select_transcripts} ) {

                                        print $SELECT_FH ";" . $key . "="
                                          . $record{INFO_key_value}{$key};
                                    }
                                    print STDOUT ";" . $key . "="
                                      . $record{INFO_key_value}{$key};
                                }
                            }
                            else {

                                if ( $record{select_transcripts} ) {

                                    print $SELECT_FH ";" . $key;
                                }
                                print STDOUT ";" . $key;
                            }
                        }
                        $counter++;
                    }
                }

                foreach my $key ( keys %{ $record{INFO_addition} } ) {

                    if ( $record{select_transcripts} ) {

                        print $SELECT_FH ";" . $key . "=" . $record{INFO_addition}{$key};
                    }
                    print STDOUT ";" . $key . "=" . $record{INFO_addition}{$key};
                }
                if ( $record{select_transcripts} ) {

                    foreach my $key ( keys %{ $record{INFO_addition_select_feature} } ) {

                        print $SELECT_FH ";" . $key . "="
                          . $record{INFO_addition_select_feature}{$key};
                    }
                }
                foreach my $key ( keys %{ $record{INFO_addition_range_feature} } ) {

                    print STDOUT ";" . $key . "="
                      . $record{INFO_addition_range_feature}{$key};
                }

                if ( $record{select_transcripts} ) {

                    print $SELECT_FH "\t";
                }
                print STDOUT "\t";
            }
            if ( $line_elements_counter > 7 ) {

                if ( $record{select_transcripts} ) {

                    print $SELECT_FH $record{ $vcf_format_columns[$line_elements_counter]
                      }
                      . "\t";
                }
                print STDOUT $record{ $vcf_format_columns[$line_elements_counter] }
                  . "\t";
            }
        }
        if ( $record{select_transcripts} ) {

            print $SELECT_FH "\n";
        }
        print STDOUT "\n";
    }
    if ($select_feature_file) {

        close($SELECT_FH);
    }
    $log->info("Finished Processing VCF\n");
}

sub parse_vep_csq {

## Function : Parse VEP CSQ field
## Returns  :
## Arguments: $consequence_href             => Variant consequence {REF}
##          : $consequence_severity_href    => Consequence severity for SO-terms {REF}
##          : $per_gene                     => Only collect most severe transcript per gene
##          : $pli_score_href               => Pli score hash
##          : $record_href                  => VCF record {REF}
##          : $select_data_href             => Select file data {REF}
##          : $vep_format_field_column_href => VEP format columns {REF}

    my ($arg_href) = @_;

    ## Default(s)

    ## Flatten argument(s)
    my $consequence_href;
    my $consequence_severity_href;
    my $per_gene;
    my $pli_score_href;
    my $record_href;
    my $select_data_href;
    my $vep_format_field_column_href;

    my $tmpl = {
        consequence_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$consequence_href,
            strict_type => 1,
        },
        consequence_severity_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$consequence_severity_href,
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

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    use MIP::File::Format::Vcf qw{ set_in_consequence_hash };

    ## Convert between hgnc_id and hgnc_symbol
    my %hgnc_map;

    if ( $record_href->{INFO_key_value}{CSQ} ) {

        my %most_severe_transcript;    #Collect most severe transcript per gene
        my @transcripts =
          split( /,/, $record_href->{INFO_key_value}{CSQ} );    #Split into transcripts

        foreach my $transcript (@transcripts) {

            my @transcript_effects = split( /\|/, $transcript );    #Split in "|"

            ## Alias hgnc_id column number
            my $hgnc_id_column_ref = \$vep_format_field_column_href->{HGNC_ID};

            ## Alias hgnc_symbol column number
            my $hgnc_symbol_column_ref = \$vep_format_field_column_href->{SYMBOL};

            ## If gene
            if ( ( defined( $transcript_effects[$$hgnc_id_column_ref] ) )
                && $transcript_effects[$$hgnc_id_column_ref] ne "" )
            {

                ## Alias HGNC ID
                my $hgnc_id_ref = \$transcript_effects[$$hgnc_id_column_ref];

                ## Add symbol to hgnc map
                $hgnc_map{$$hgnc_id_ref} = $transcript_effects[$$hgnc_symbol_column_ref];
                my $transcript_id_ref =
                  \$transcript_effects[ $vep_format_field_column_href->{Feature} ]
                  ;    #Alias transcript_id
                my $allele_ref =
                  \$transcript_effects[ $vep_format_field_column_href->{Allele} ]
                  ;    #Alias allele

                my $consequence_field =
                  $transcript_effects[ $vep_format_field_column_href->{Consequence} ];

                ## Parse the most severe consequence or prediction to gene
                parse_vep_csq_consequence(
                    {
                        allele            => $$allele_ref,
                        consequence_field => $consequence_field,
                        consequence_href  => $consequence_href,
                        hgnc_id_ref       => $$hgnc_id_ref,
                        transcript        => $transcript,
                    }
                );

                if ( !$per_gene ) {

                    if ( $select_data_href->{$$hgnc_id_ref} )
                    {    #Exists in selected Features

                        push( @{ $record_href->{select_transcripts} }, $transcript )
                          ;    #Add all transcripts to selected transcripts
                    }
                    push( @{ $record_href->{range_transcripts} }, $transcript )
                      ;        #Add all transcripts to range transcripts
                }
            }
            else {

                ## Intergenic
                push( @{ $record_href->{range_transcripts} }, $transcript )
                  ;            #Add all transcripts to range transcripts
            }
        }
        my @most_severe_range_consequences;
        my @most_severe_select_consequences;
        my $most_severe_range_pli  = 0;
        my $most_severe_select_pli = 0;

      GENE:
        for my $gene ( keys %{$consequence_href} ) {

          ALLEL:
            for my $allele ( keys %{ $consequence_href->{$gene} } ) {

                ## Get hgnc_symbol
                my $hgnc_symbol = $hgnc_map{$gene};

                ## For pli value and if current pli is more than stored
                if ( exists $pli_score_href->{$hgnc_symbol}
                    and $most_severe_range_pli < $pli_score_href->{$hgnc_symbol} )
                {

                    if ( $select_data_href->{$gene} ) {

                        $most_severe_select_pli = $pli_score_href->{$hgnc_symbol};
                    }
                    $most_severe_range_pli = $pli_score_href->{$hgnc_symbol};
                }

                ## Exists in selected features
                if ( $select_data_href->{$gene} ) {

                    push( @most_severe_select_consequences,
                        $consequence_href->{$gene}{$allele}{most_severe_consequence} );

                    if ($per_gene) {

                        push(
                            @{ $record_href->{select_transcripts} },
                            $consequence_href->{$gene}{$allele}{most_severe_transcript}
                        );
                    }
                }
                push( @most_severe_range_consequences,
                    $consequence_href->{$gene}{$allele}{most_severe_consequence} );

                if ($per_gene) {

                    push(
                        @{ $record_href->{range_transcripts} },
                        $consequence_href->{$gene}{$allele}{most_severe_transcript}
                    );
                }
            }
        }
        if (@most_severe_select_consequences) {

            $record_href->{INFO_addition_select_feature}{most_severe_consequence} =
              join( ",", @most_severe_select_consequences );
        }
        if (@most_severe_range_consequences) {

            $record_href->{INFO_addition_range_feature}{most_severe_consequence} =
              join( ",", @most_severe_range_consequences );
        }

        ## Mainly for SV BNDs without consequence and within a gene
        if (    not keys %{$consequence_href}
            and not exists $record_href->{range_transcripts} )
        {
            ## Add all transcripts to range transcripts
            push( @{ $record_href->{range_transcripts} }, @transcripts );
        }
        if ($most_severe_select_pli) {

            $record_href->{INFO_addition_select_feature}{most_severe_pli} =
              $most_severe_select_pli;
        }
        if ($most_severe_range_pli) {

            $record_href->{INFO_addition_range_feature}{most_severe_pli} =
              $most_severe_range_pli;
        }
    }
}

sub add_field_to_element {

##add_field_to_element

##Function : Adds adds a field to an element.
##Returns  : ""
##Arguments: $selected_transcript_tracker_ref, $selected_line_ref, $line_ref, $value_ref, $separator
##         : $selected_transcript_tracker_ref => The selected transcript tracker for belonging to select file {REF}
##         : $selected_line_ref               => Selected line to add annotations to {REF}
##         : $line_ref                        => Variant line to add annotations to {REF}
##         : $value_ref                       => Field value {REF}
##         : $separator                       => Separator for field

    my ($arg_href) = @_;

    ## Default(s)
    my $separator;

    ## Flatten argument(s)
    my $selected_transcript_tracker_ref;
    my $selected_line_ref;
    my $line_ref;
    my $value_ref;

    my $tmpl = {
        selected_transcript_tracker_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$selected_transcript_tracker_ref
        },
        selected_line_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$selected_line_ref
        },
        line_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$line_ref
        },
        separator => { default => ":", strict_type => 1, store => \$separator },
        value_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$value_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    if ($$selected_transcript_tracker_ref) {

        $$selected_line_ref .= $separator . $$value_ref;
    }

    ## Always include all transcripts in orphan list
    $$line_ref .= $separator . $$value_ref;
}

sub collect_consequence_genes {

##collect_consequence_genes

##Function : Collects all consequence and predictors per gene and adds info to line to be written.
##Returns  : ""
##Arguments: $consequence_href, $select_data_href, $fields_ref, $selected_variant_line_ref, $variantVariantLineRef
##         : $consequence_href          => Consequence(s) for each gene {REF}
##         : $select_data_href          => Select file data {REF}
##         : $fields_ref                => Features to be processed as determined by CSQ {REF}
##         : $selected_variant_line_ref => Selected line to add annotations to {REF}
##         : $variantVariantLineRef     => Variant line to add annotations to {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consequence_href;
    my $select_data_href;
    my $fields_ref;
    my $selected_variant_line_ref;
    my $variant_line_ref;

    my $tmpl = {
        consequence_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$consequence_href
        },
        select_data_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$select_data_href
        },
        fields_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$fields_ref
        },
        selected_variant_line_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$selected_variant_line_ref
        },
        variant_line_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$variant_line_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    my %gene_counter;
    my %selected_gene_counter;
    my @temp_fields;
    my @selected_temp_fields;

    for (
        my $field_counter = 0 ;
        $field_counter < scalar(@$fields_ref) ;
        $field_counter++
      )
    {    #Set transcript counter to "0"

        $selected_gene_counter{$field_counter} = 0;
        $gene_counter{$field_counter}          = 0;
    }

  GENES:
    for my $gene ( keys %$consequence_href ) {

      FIELDS:
        for (
            my $field_counter = 0 ;
            $field_counter < scalar(@$fields_ref) ;
            $field_counter++
          )
        {

            if ( $select_data_href->{$gene} ) {    #Exists in selected Features

                collect_consequence_field(
                    {
                        gene_counter_href   => \%selected_gene_counter,
                        consequence_href    => $consequence_href,
                        fields_ref          => \@{$fields_ref},
                        selected_fields_ref => \@selected_temp_fields,
                        field_counter_ref   => \$field_counter,
                        gene_ref            => \$gene,
                    }
                );
            }

            ## Always include all transcripts in research list
            collect_consequence_field(
                {
                    gene_counter_href   => \%gene_counter,
                    consequence_href    => $consequence_href,
                    fields_ref          => \@{$fields_ref},
                    selected_fields_ref => \@temp_fields,
                    field_counter_ref   => \$field_counter,
                    gene_ref            => \$gene,
                }
            );
        }
    }

    ## Adds to present line
    add_to_line(
        {
            fields_ref      => \@{$fields_ref},
            temp_fields_ref => \@selected_temp_fields,
            line_ref        => \$$selected_variant_line_ref,
        }
    );

    ## Adds to present line
    add_to_line(
        {
            fields_ref      => \@{$fields_ref},
            temp_fields_ref => \@temp_fields,
            line_ref        => \$$variant_line_ref,
        }
    );
}

sub collect_consequence_field {

##collect_consequence_field

##Function : Collects consequences for features in @feature_fields to temporary array for adding to line once all information are collected.
##Returns  : ""
##Arguments: $gene_counter_href, $consequence_href, $fields_ref, $selected_fields_ref, $field_counter_ref, $gene_ref
##         : $gene_counter_href   => Counts the number of transcripts per gene {REF}
##         : $consequence_href    => Consequence(s) for each gene {REF}
##         : $fields_ref          => Features to be processed as determined by CSQ {REF}
##         : $selected_fields_ref => Selected array {REF}
##         : $field_counter_ref   => Field number in feature {REF}
##         : $gene_ref            => The gene symbol {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $gene_counter_href;
    my $consequence_href;
    my $fields_ref;
    my $selected_fields_ref;
    my $field_counter_ref;
    my $gene_ref;

    my $tmpl = {
        gene_counter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$gene_counter_href
        },
        consequence_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$consequence_href
        },
        fields_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$fields_ref
        },
        selected_fields_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$selected_fields_ref
        },
        field_counter_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$field_counter_ref
        },
        gene_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$gene_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    my $field_ref = \${$fields_ref}[$$field_counter_ref];    #Alias

    if ( !$$gene_counter_href{$$field_counter_ref} ) {       #First time

        my $allel_counter = 0;

      ALLELS:
        for my $allele ( keys %{ $consequence_href->{$$gene_ref} } ) {    #All alleles

            if ( defined( $$consequence_href{$$gene_ref}{$allele}{$$field_ref} )
                && ( $$consequence_href{$$gene_ref}{$allele}{$$field_ref} ne "" ) )
            {    #If feature exists - else do nothing

                if ( !$allel_counter ) {

                    $$selected_fields_ref[$$field_counter_ref] .= ";"
                      . $$field_ref . "="
                      . $$gene_ref . ":"
                      . $allele . "|"
                      . $$consequence_href{$$gene_ref}{$allele}{$$field_ref};
                }
                else {

                    $$selected_fields_ref[$$field_counter_ref] .= ","
                      . $$gene_ref . ":"
                      . $allele . "|"
                      . $$consequence_href{$$gene_ref}{$allele}{$$field_ref};
                }
                $$gene_counter_href{$$field_counter_ref}++;
                $allel_counter++;
            }
        }
    }
    else {    #Subsequent passes

      ALLELS:
        for my $allele ( keys %{ $consequence_href->{$$gene_ref} } ) {    #All alleles

            if ( defined( $$consequence_href{$$gene_ref}{$allele}{$$field_ref} )
                && ( $$consequence_href{$$gene_ref}{$allele}{$$field_ref} ne "" ) )
            {    #If feature exists - else do nothing

                $$selected_fields_ref[$$field_counter_ref] .= ","
                  . $$gene_ref . ":"
                  . $allele . "|"
                  . $$consequence_href{$$gene_ref}{$allele}{$$field_ref};
            }
        }
    }
}

sub add_to_line {

##add_to_line

##Function : Adds to present line.
##Returns  : ""
##Arguments: $fields_ref, $temp_fields_ref, $line_ref
##         : $fields_ref      => Features to be processed as determined by CSQ {REF}
##         : $temp_fields_ref => Annotations for feature {REF}
##         : $line_ref        => Line to add to {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fields_ref;
    my $temp_fields_ref;
    my $line_ref;

    my $tmpl = {
        fields_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$fields_ref
        },
        temp_fields_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$temp_fields_ref
        },
        line_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$line_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    for (
        my $field_counter = 0 ;
        $field_counter < scalar(@$fields_ref) ;
        $field_counter++
      )
    {

        if ( defined( $$temp_fields_ref[$field_counter] ) ) {

            $$line_ref .= $$temp_fields_ref[$field_counter];
        }
    }
}

sub find_af {

##find_af

##Function : Adds the least alternative allele(s) frequency to each line
##Returns  : ""
##Arguments: $elements_ref, $regexp
##         : $elements_ref => The INFO array {REF}
##         : $regexp       => The regexp to used to locate correct ID field

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $elements_ref;
    my $regexp;

    my $tmpl = {
        elements_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$elements_ref
        },
        regexp => { required => 1, defined => 1, strict_type => 1, store => \$regexp },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    my $temp_maf;

    for my $element (@$elements_ref) {

        if ( $element =~ /$regexp/ ) {    #Find the key=value field

            my @value = split( /=/, $element );    #Split key=value pair

            $temp_maf = $value[1]; #Collect whole string to represent all possible alleles
        }
    }
    return $temp_maf;
}

sub FindLCAF {

##FindLCAF

##Function : Adds the least common alternative allele frequency to each line
##Returns  : ""
##Arguments: $elements_ref, $regexp, $frequencyPosition
##         : $elements_ref      => The INFO array {REF}
##         : $regexp            => The regexp to used to locate correct ID field
##         : $frequencyPosition => The position to extract frequency from

    my ($arg_href) = @_;

    ## Default(s)
    my $frequencyPosition;

    ## Flatten argument(s)
    my $elements_ref;
    my $regexp;

    my $tmpl = {
        elements_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$elements_ref
        },
        regexp => { required => 1, defined => 1, strict_type => 1, store => \$regexp },
        frequencyPosition => {
            default     => 0,
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$frequencyPosition
        },
    };

    check( $tmpl, $arg_href, 1 ) or die qw[Could not parse arguments!];

    my $temp_maf;

    for my $element (@$elements_ref) {

        if ( $element =~ /$regexp/ ) {    #Find the key=value field

            my @value = split( /=/, $element );    #Split key=value pair

            my @temp_mafs =
              sort { $a <=> $b }
              grep { $_ ne "." }
              split( ",", $value[1] )
              ; #Split on ",", remove entries containing only "." and sort remaining entries numerically ascending order

            if ( scalar(@temp_mafs) > 0 ) {

                ## We are interested in the least common allele listed for this position. We cannot connect the frequency position in the list and the multiple alternative alleles. So the best we can do is report the least common allele frequency for multiple alternative allels. Unless the least common frequency is lower than the frequency defined as pathogenic for rare disease (usually 0.01) then this will work. In that case this will be a false positive, but it is better than taking the actual MAF which would be a false negative if the pathogenic variant found in the patient(s) has a lower frequency than the MAF.
                $temp_maf = $temp_mafs[$frequencyPosition];
            }
        }
    }
    return $temp_maf;
}
