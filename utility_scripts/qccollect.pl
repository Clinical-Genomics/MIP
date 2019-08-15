#!/usr/bin/env perl

#### Collects MPS QC from MIP. Loads information on files to examine and values
#### to extract from in YAML format and outputs exracted metrics in YAML format.

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use Cwd qw{ abs_path };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use English qw{ -no_match_vars };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use Pod::Usage;
use Pod::Text;
use POSIX;
use utf8;
use warnings qw{ FATAL utf8 };

$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

## CPANM
use autodie qw{ open close :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( $Bin, q{lib} );
use MIP::Check::Modules qw{ check_perl_modules };
use MIP::Constants qw{ $COLON $NEWLINE $SPACE $UNDERSCORE };
use MIP::File::Format::Yaml qw{ load_yaml write_yaml };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Qccollect
  qw{ define_evaluate_metric evaluate_case_qc_parameters evaluate_sample_qc_parameters };
use MIP::Qc_data qw{ set_qc_data_recipe_info };
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

BEGIN {

    require MIP::Check::Modules;
    use MIP::Check::Modules qw{ check_perl_modules parse_cpan_file };

    my @modules =
      parse_cpan_file { cpanfile_path => catfile( $Bin, qw{ definitions cpanfile } ), };

    ## Evaluate that all modules required are installed
    check_perl_modules(
        {
            modules_ref  => \@modules,
            program_name => $PROGRAM_NAME,
        }
    );
}

my $VERSION = q{2.1.3};

my ( $evaluate_plink_gender, $print_regexp, $regexp_file, $sample_info_file,
    $skip_evaluation, );

## Scalar parameters with defaults
my ( $log_file, $print_regexp_outfile, $outfile, ) =
  ( catfile( cwd(), q{qccollect.log} ), q{qc_regexp.yaml}, q{qcmetrics.yaml}, );

## Save final output data
my %qc_data;

## Save header(s) in each outfile
my %qc_header;

## Save data in each outfile
my %qc_recipe_data;

### User Options
GetOptions(
    q{si|sample_info_file:s}        => \$sample_info_file,
    q{r|regexp_file:s}              => \$regexp_file,
    q{o|outfile:s}                  => \$outfile,
    q{preg|print_regexp}            => \$print_regexp,
    q{prego|print_regexp_outfile:s} => \$print_regexp_outfile,
    q{ske|skip_evaluation}          => \$skip_evaluation,
    q{epg|evaluate_plink_gender}    => \$evaluate_plink_gender,
    q{l|log_file:s}                 => \$log_file,
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
        log_name  => uc q{mip_qccollect},
    }
);

## Write default regexp to YAML if demanded
regexp_to_yaml(
    {
        is_print_regexp      => $print_regexp,
        log                  => $log,
        print_regexp_outfile => $print_regexp_outfile,
    }
);

if ( not $sample_info_file ) {

    $log->info($USAGE);
    $log->fatal( q{Must supply a '-sample_info_file' (supply whole path)}, $NEWLINE );
    exit;
}
if ( not $regexp_file ) {

    $log->info($USAGE);
    $log->fatal( q{Must supply a '-regexp_file' (supply whole path)}, $NEWLINE );
    exit;
}

###########
####MAIN###
###########

## Loads a YAML file into an arbitrary hash and returns it
$log->info( q{Loading: } . $sample_info_file );
my %sample_info = load_yaml( { yaml_file => $sample_info_file, } );
$log->info( q{Loaded: } . $sample_info_file );

## Loads a reg exp file into an arbitrary hash
$log->info( q{Loading: } . $regexp_file );
my %regexp = load_yaml( { yaml_file => $regexp_file, } );
$log->info( q{Loaded: } . $regexp_file );

## Set qccollect version to qc_data hash
set_qc_data_recipe_info(
    {
        key          => q{version},
        qc_data_href => \%qc_data,
        recipe_name  => q{qccollect},
        value        => $VERSION,
    }
);

## Set supplied regexp file to qc_data hash
set_qc_data_recipe_info(
    {
        key          => q{regexp_file},
        qc_data_href => \%qc_data,
        recipe_name  => q{qccollect},
        value        => $regexp_file,
    }
);

## Extracts all Qc data on sample_id level using information in %sample_info and %regexp
sample_qc(
    {
        qc_data_href        => \%qc_data,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
        regexp_href         => \%regexp,
        sample_info_href    => \%sample_info,
    }
);

## Extracts all Qc data on case level using information in %sample_info_file and %regexp
case_qc(
    {
        qc_data_href        => \%qc_data,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
        regexp_href         => \%regexp,
        sample_info_href    => \%sample_info,
    }
);

## Defines recipes, metrics and thresholds to evaluate
my %evaluate_metric = define_evaluate_metric(
    {
        sample_info_href => \%sample_info,
    }
);

if ( not $skip_evaluation ) {

    ## Evaluate the metrics
    evaluate_case_qc_parameters(
        {
            evaluate_metric_href => \%evaluate_metric,
            qc_data_href         => \%qc_data,
        }
    );

    evaluate_sample_qc_parameters(
        {
            evaluate_metric_href => \%evaluate_metric,
            qc_data_href         => \%qc_data,
        }
    );
}

## Writes a qc data hash to file
write_yaml(
    {
        yaml_file_path => $outfile,
        yaml_href      => \%qc_data,
    }
);
$log->info( q{Wrote: } . $outfile );

######################
####Sub routines######
######################

sub build_usage {

##Function : Build the USAGE instructions
##Returns  :
##Arguments: $program_name => Name of the script

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
 $program_name [options] -si [sample_info.yaml] -r [regexp.yaml] -o [outfile]
   -si/--sample_info_file        Sample info file (YAML format, Supply whole path, mandatory)
   -r/--regexp_file              Regular expression file (YAML format, Supply whole path, mandatory)
   -o/--outfile                  Data file output (Supply whole path, defaults to "qcmetrics.yaml")
   -preg/--print_regexp          Print the regexp used at CMMS switch (defaults to "0" (=no))
   -prego/--print_regexp_outfile Regexp YAML outfile (defaults to "qc_regexp.yaml")
   -ske/--skip_evaluation        Skip evaluation step (boolean)
   -epg/--evaluate_plink_gender  Evaluate plink gender (boolean)
   -l/--log_file                 Log file (Default: "qccollect.log")
   -h/--help                     Display this help message
   -v/--version                  Display version
END_USAGE
}

sub case_qc {

## Function : Extracts all Qc data on case level using information in %sample_info_file and %regexp
## Returns  :
## Arguments: $qc_data_href        => Qc data hash {REF}
##          : $qc_header_href      => Save header(s) in each outfile {REF}
##          : $qc_recipe_data_href => Hash to save data in each outfile {REF}
##          : $regexp_href         => RegExp hash {REF}
##          : $sample_info_href    => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $qc_header_href;
    my $qc_recipe_data_href;
    my $regexp_href;
    my $sample_info_href;

    my $tmpl = {
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        qc_header_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_header_href,
            strict_type => 1,
        },
        qc_recipe_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_recipe_data_href,
            strict_type => 1,
        },
        regexp_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$regexp_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qccollect qw{ plink_gender_check relation_check };
    use MIP::Qc_data
      qw{ parse_qc_recipe_data parse_regexp_hash_and_collect set_qc_data_recipe_info };
    use MIP::Sample_info qw{ get_sample_info_case_recipe_attributes };

  RECIPE:
    for my $recipe ( keys %{ $sample_info_href->{recipe} } ) {

        my %attribute = get_sample_info_case_recipe_attributes(
            {
                recipe_name      => $recipe,
                sample_info_href => \%sample_info,
            }
        );

        my $outdirectory = $attribute{outdirectory};
        my $outfile      = $attribute{outfile};

        if ( exists $attribute{path} ) {

            ( $outfile, $outdirectory ) =
              fileparse( $attribute{path} );
        }

        ## Set package executable version from recipe to metrics hash
        set_qc_data_recipe_info(
            {
                key          => q{version},
                qc_data_href => $qc_data_href,
                recipe_name  => $recipe,
                value        => $attribute{version},
            }
        );

        ## Parses the RegExpHash structure to identify if the info is
        ## 1) Paragraf section(s) (both header and data line(s)
        ## 2) Seperate data line
        parse_regexp_hash_and_collect(
            {
                outdirectory        => $outdirectory,
                outfile             => $outfile,
                qc_recipe_data_href => $qc_recipe_data_href,
                qc_header_href      => $qc_header_href,
                recipe              => $recipe,
                regexp_href         => $regexp_href,
            }
        );

        ## Parse qc_recipe_data and extract information to qc_data
        parse_qc_recipe_data(
            {
                qc_data_href        => $qc_data_href,
                qc_header_href      => $qc_header_href,
                qc_recipe_data_href => $qc_recipe_data_href,
                recipe              => $recipe,
                regexp_href         => $regexp_href,
                sample_info_href    => $sample_info_href,
            }
        );

        ## Check gender for sample_id
        if (    $recipe eq q{plink_sexcheck}
            and $evaluate_plink_gender )
        {

            ## Check that assumed gender is supported by variants on chrX and chrY
            plink_gender_check(
                {
                    qc_data_href     => $qc_data_href,
                    sample_info_href => $sample_info_href,
                }
            );
        }
    }

    if (    defined $qc_data_href->{recipe}{relation_check}{sample_relation_check}
        and defined $qc_data_href->{recipe}{pedigree_check}{sample_order} )
    {

        relation_check(
            {
                qc_data_href => $qc_data_href,
                relationship_values_ref =>
                  \@{ $qc_data_href->{recipe}{relation_check}{sample_relation_check} },
                sample_info_href => $sample_info_href,
                sample_orders_ref =>
                  \@{ $qc_data_href->{recipe}{pedigree_check}{sample_order} },
            }
        );
    }
    return;
}

sub sample_qc {

## Function : Collects all sample qc in files defined by sample_info_file and regular expressions defined by regexp.
## Returns  :
## Arguments: $qc_data_href        => Qc data hash {REF}
##          : $qc_header_href      => Save header(s) in each outfile {REF}
##          : $qc_recipe_data_href => Hash to save data in each outfile {REF}
##          : $regexp_href         => RegExp hash {REF}
##          : $sample_info_href    => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $qc_header_href;
    my $qc_recipe_data_href;
    my $regexp_href;
    my $sample_info_href;

    my $tmpl = {
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        qc_header_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_header_href,
            strict_type => 1,
        },
        qc_recipe_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_recipe_data_href,
            strict_type => 1,
        },
        regexp_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$regexp_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qccollect qw{ chanjo_gender_check };
    use MIP::Qc_data qw{ parse_qc_recipe_data parse_regexp_hash_and_collect };
    use MIP::Sample_info qw{ get_sample_info_sample_recipe_attributes };

  SAMPLE_ID:
    for my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

      RECIPE:
        for my $recipe ( keys %{ $sample_info_href->{sample}{$sample_id}{recipe} } ) {

          INFILE:
            for my $infile (
                keys %{ $sample_info_href->{sample}{$sample_id}{recipe}{$recipe} } )
            {

                my %attribute = get_sample_info_sample_recipe_attributes(
                    {
                        infile           => $infile,
                        recipe_name      => $recipe,
                        sample_id        => $sample_id,
                        sample_info_href => \%sample_info,
                    }
                );

                my $outdirectory = $attribute{outdirectory};
                my $outfile      = $attribute{outfile};

                if ( exists $attribute{path} ) {

                    ( $outfile, $outdirectory ) =
                      fileparse( $attribute{path} );
                }

                ## Parses the RegExpHash structure to identify if the info is
                ## 1) Paragraf section(s) (both header and data line(s)
                ## 2) Seperate data line
                parse_regexp_hash_and_collect(
                    {
                        outdirectory        => $outdirectory,
                        outfile             => $outfile,
                        qc_header_href      => $qc_header_href,
                        qc_recipe_data_href => $qc_recipe_data_href,
                        recipe              => $recipe,
                        regexp_href         => $regexp_href,
                    }
                );

                ## Parse qc_recipe_data and extract information to qc_data
                parse_qc_recipe_data(
                    {
                        infile              => $infile,
                        qc_data_href        => $qc_data_href,
                        qc_header_href      => $qc_header_href,
                        qc_recipe_data_href => $qc_recipe_data_href,
                        recipe              => $recipe,
                        regexp_href         => $regexp_href,
                        sample_id           => $sample_id,
                        sample_info_href    => $sample_info_href,
                    }
                );

                ## Check gender for sample_id
                if ( $recipe eq q{chanjo_sexcheck} ) {

                    ## Check that assumed gender is supported by coverage on chrX and chrY
                    chanjo_gender_check(
                        {
                            infile           => $infile,
                            qc_data_href     => $qc_data_href,
                            recipe_name      => $recipe,
                            sample_id        => $sample_id,
                            sample_info_href => $sample_info_href,
                        }
                    );
                }
            }
        }
    }
    return;
}
