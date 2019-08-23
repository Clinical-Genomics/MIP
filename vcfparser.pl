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
use Modern::Perl qw{ 2018 };
use Readonly;
use Set::IntervalTree;

## MIPs lib/
use lib catdir( $Bin, q{lib} );
use MIP::Check::Modules qw{ check_perl_modules };
use MIP::Constants
  qw{ %ANALYSIS $COLON $COMMA $EQUALS $NEWLINE %SO_CONSEQUENCE_SEVERITY $SEMICOLON $SPACE $TAB $UNDERSCORE };
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

    ## Constants
    Readonly my $FILTER_COLUMN_INDEX => 6;
    Readonly my $FORMAT_COLUMN_INDEX => 8;

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
        my %vcf_record;
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
        %noid_region = tree_annotations(
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
                    consequence_severity_href    => $consequence_severity_href,
                    per_gene                     => $per_gene,
                    pli_score_href               => $pli_score_href,
                    record_href                  => \%vcf_record,
                    select_data_href             => $select_data_href,
                    vep_format_field_column_href => \%vep_format_field_column,
                }
            );
        }

	## Writing vcf record to files
        my $last_index     = $FILTER_COLUMN_INDEX;
        my $last_separator = $TAB;

        ## If we do not need to split the CSQ field we can print the entire vcf
        ## record line
        if ( not $parse_vep ) {

            $last_index     = $#line_elements;
            $last_separator = $NEWLINE;
        }

        ## Add until INFO field
        if ( $vcf_record{select_transcripts} ) {

            print {$SELECT_FH} join( $TAB, @line_elements[ 0 .. $last_index ] ),
              $last_separator;
        }
        print {*STDOUT} join( $TAB, @line_elements[ 0 .. $last_index ] ), $last_separator;

        if ($parse_vep) {

            my $counter = 0;
            if ( exists $vcf_record{INFO_key_value}{CSQ} and $vcf_record{INFO_key_value}{CSQ} ) {

                if ( $vcf_record{range_transcripts} ) {

                    print {*STDOUT} q{CSQ} . $EQUALS . join $COMMA,
                      @{ $vcf_record{range_transcripts} };
                }
                if ( $vcf_record{select_transcripts} ) {

                    print {$SELECT_FH} q{CSQ} . $EQUALS . join $COMMA,
                      @{ $vcf_record{select_transcripts} };
                }
                delete $vcf_record{INFO_key_value}{CSQ};
                $counter++;
            }

          KEY_VALUE_PAIR:
            while ( my ( $key, $value ) = each %{ $vcf_record{INFO_key_value} } ) {

                my $info_string = $SEMICOLON . $key;
                if ( not $counter ) {

                    $info_string = $key;
                }

                if ( defined $value ) {

                    $info_string .= $EQUALS . $value;
                }

                if ( $vcf_record{select_transcripts} ) {

                    print {$SELECT_FH} $info_string;
                }
                print {*STDOUT} $info_string;
                $counter++;
            }

            foreach my $key ( keys %{ $vcf_record{INFO_addition} } ) {

                if ( $vcf_record{select_transcripts} ) {

                    print {$SELECT_FH} $SEMICOLON . $key . $EQUALS
                      . $vcf_record{INFO_addition}{$key};
                }
                print {*STDOUT} $SEMICOLON . $key . $EQUALS
                  . $vcf_record{INFO_addition}{$key};
            }
            if ( $vcf_record{select_transcripts} ) {

                foreach my $key ( keys %{ $vcf_record{INFO_addition_select_feature} } ) {

                    print {$SELECT_FH} $SEMICOLON . $key . $EQUALS
                      . $vcf_record{INFO_addition_select_feature}{$key};
                }
            }
            foreach my $key ( keys %{ $vcf_record{INFO_addition_range_feature} } ) {

                print {*STDOUT} $SEMICOLON . $key . $EQUALS
                  . $vcf_record{INFO_addition_range_feature}{$key};
            }

            ## After INFO to the final column
            if ( $vcf_record{select_transcripts} ) {

                say {$SELECT_FH} $TAB, join $TAB,
                  @line_elements[ $FORMAT_COLUMN_INDEX .. $#line_elements ];
            }
            say {*STDOUT} $TAB, join $TAB,
              @line_elements[ $FORMAT_COLUMN_INDEX .. $#line_elements ];
        }
    }
    if ($select_feature_file) {

        close $SELECT_FH;
    }
    $log->info( q{Finished Processing VCF} . $NEWLINE );
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

    ## Flatten argument(s)
    my $consequence_href;
    my $consequence_severity_href;
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

        ## Mainly for SV BNDs without consequence and within a hgnc_id
        if (    not keys %{$consequence_href}
            and not exists $record_href->{range_transcripts} )
        {

            ## Add all transcripts to range transcripts
            add_transcript_to_feature_file(
                {
                    transcripts_ref => \@transcripts,
                    vcf_record_href => $record_href,
                }
            );
        }
    }

    set_most_severe_ann_to_vcf_record(
        {
            feature_type_keys_ref    => \@feature_type_keys,
            most_severe_feature_href => \%most_severe_feature,
            most_severe_pli_href     => \%most_severe_pli,
            vcf_record_href          => $record_href,
        }
    );

    return;
}
