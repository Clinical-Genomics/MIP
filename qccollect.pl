#!/usr/bin/env perl

#### Collects MPS QC from MIP. Loads information on files to examine and values to extract from in YAML format and outputs exracted metrics in YAML format.

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use Cwd qw{ abs_path };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catdir catfile devnull };
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
use Modern::Perl qw{ 2017 };
use Readonly;

##MIPs lib/
use lib catdir( $Bin, q{lib} );
use MIP::Check::Modules qw{ check_perl_modules };
use MIP::File::Format::Yaml qw{ load_yaml write_yaml };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

BEGIN {

    require MIP::Check::Modules;
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

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

my ( $sample_info_file, $regexp_file, $print_regexp, $skip_evaluation,
    $evaluate_plink_gender );

## Scalar parameters with defaults
my ( $outfile, $print_regexp_outfile, $log_file ) =
  ( q{qcmetrics.yaml}, q{qc_regexp.yaml}, catfile( cwd(), q{qccollect.log} ) );

my ( %qc_data, %evaluate_metric );

## Save header(s) in each outfile
my %qc_header;

## Save data in each outfile
my %qc_recipe_data;

my $qccollect_version = q{2.1.0};

###User Options
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
    q{h|help} => sub { say STDOUT $USAGE; exit; },
    ## Display version number
    q{v|version} => sub {
        say STDOUT $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $qccollect_version,
          $NEWLINE;
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
        log_name  => q{Qccollect},
    }
);

if ($print_regexp) {

    ## Write default regexp to YAML
    regexp_to_yaml( { print_regexp_outfile => $print_regexp_outfile, } );
    $log->info( q{Wrote regexp YAML file to: } . $print_regexp_outfile );
    exit;
}

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
my %sample_info = load_yaml( { yaml_file => $sample_info_file, } );
$log->info( q{Loaded: } . $sample_info_file );

## Loads a YAML file into an arbitrary hash and returns it
my %regexp = load_yaml( { yaml_file => $regexp_file, } );
$log->info( q{Loaded: } . $regexp_file );

## Extracts all qcdata on sample_id level using information in %sample_info and %regexp
sample_qc(
    {
        sample_info_href    => \%sample_info,
        regexp_href         => \%regexp,
        qc_data_href        => \%qc_data,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
    }
);

## Extracts all qcdata on case level using information in %sample_info_file and %regexp
case_qc(
    {
        sample_info_href    => \%sample_info,
        regexp_href         => \%regexp,
        qc_data_href        => \%qc_data,
        qc_header_href      => \%qc_header,
        qc_recipe_data_href => \%qc_recipe_data,
    }
);

##Add qcCollect version to qc_data yaml file
$qc_data{recipe}{qccollect}{version}     = $qccollect_version;
$qc_data{recipe}{qccollect}{regexp_file} = $regexp_file;

SAMPLE_ID:
foreach my $sample_id ( keys %{ $sample_info{sample} } ) {

    ## Defines recipes, metrics and thresholds to evaluate
    define_evaluate_metric(
        {
            sample_info_href => \%sample_info,
            sample_id        => $sample_id,
        }
    );
}

if ( not $skip_evaluation ) {

    ## Evaluate the metrics
    evaluate_qc_parameters(
        {
            qc_data_href         => \%qc_data,
            evaluate_metric_href => \%evaluate_metric,
        }
    );
}

## Writes a YAML hash to file
write_yaml(
    {
        yaml_href      => \%qc_data,
        yaml_file_path => $outfile,
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
            strict_type => 1,
            store       => \$program_name,
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

## Function : Extracts all qcdata on case level using information in %sample_info_file and %regexp
## Returns  :
## Arguments: $qc_data_href         => QCData hash {REF}
##          : $qc_header_href       => Save header(s) in each outfile {REF}
##          : $qc_recipe_data_href  => Hash to save data in each outfile {REF}
##          : $regexp_href          => RegExp hash {REF}
##          : $sample_info_href     => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $qc_header_href;
    my $qc_recipe_data_href;
    my $regexp_href;
    my $sample_info_href;

    my $tmpl = {
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        regexp_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$regexp_href,
            strict_type => 1,
        },
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## For every recipe
  RECIPE:
    for my $recipe ( keys %{ $sample_info_href->{recipe} } ) {

        my $outdirectory;
        my $outfile;

        if ( $sample_info_href->{recipe}{$recipe}{version} ) {

            ## Add version to qc_data
            $qc_data_href->{recipe}{$recipe}{version} =
              $sample_info_href->{recipe}{$recipe}{version};
        }
        if ( $sample_info_href->{recipe}{$recipe}{outdirectory} ) {

            ## Extract outdirectory
            $outdirectory =
              $sample_info_href->{recipe}{$recipe}{outdirectory};
        }
        if ( $sample_info_href->{recipe}{$recipe}{outfile} ) {

            ## Extract outfile
            $outfile = $sample_info_href->{recipe}{$recipe}{outfile};
        }
        if ( $sample_info_href->{recipe}{$recipe}{path} ) {

            ( $outfile, $outdirectory ) =
              fileparse( $sample_info_href->{recipe}{$recipe}{path} );
        }

        ## Parses the RegExpHash structure to identify if the info is 1) Paragraf section(s) (both header and data line(s)); 2) Seperate data line.
        parse_regexp_hash_and_collect(
            {
                outdirectory        => $outdirectory,
                outfile             => $outfile,
                recipe              => $recipe,
                qc_recipe_data_href => $qc_recipe_data_href,
                qc_header_href      => $qc_header_href,
                regexp_href         => $regexp_href,
            }
        );

        ## Add extracted information to qc_data
        add_to_qc_data(
            {
                evaluate_plink_gender => $evaluate_plink_gender,
                qc_data_href          => $qc_data_href,
                qc_header_href        => $qc_header_href,
                qc_recipe_data_href   => $qc_recipe_data_href,
                recipe                => $recipe,
                regexp_href           => $regexp_href,
                sample_info_href      => $sample_info_href,
            }
        );
    }
    return;
}

sub sample_qc {

## Function : Collects all sample qc in files defined by sample_info_file and regular expressions defined by regexp.
## Returns  :
## Arguments: $sample_info_href     => Info on samples and case hash {REF}
##          : $regexp_href          => RegExp hash {REF}
##          : $qc_data_href         => QCData hash {REF}
##          : $qc_header_href       => Save header(s) in each outfile {REF}
##          : $qc_recipe_data_href  => Hash to save data in each outfile {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $regexp_href;
    my $qc_data_href;
    my $qc_header_href;
    my $qc_recipe_data_href;

    my $tmpl = {
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        regexp_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$regexp_href
        },
        qc_data_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$qc_data_href
        },
        qc_header_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$qc_header_href
        },
        qc_recipe_data_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$qc_recipe_data_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  SAMPLE_ID:
    for my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

      RECIPE:
        for my $recipe ( keys %{ $sample_info_href->{sample}{$sample_id}{recipe} } ) {

          INFILE:
            for my $infile (
                keys %{ $sample_info_href->{sample}{$sample_id}{recipe}{$recipe} } )
            {

                my $outdirectory;
                my $outfile;
                if ( $sample_info_href->{sample}{$sample_id}{recipe}{$recipe}
                    {$infile}{outdirectory} )
                {

                    ## Extract OutDirectory
                    $outdirectory =
                      $sample_info_href->{sample}{$sample_id}{recipe}
                      {$recipe}{$infile}{outdirectory};
                }
                if ( $sample_info_href->{sample}{$sample_id}{recipe}{$recipe}
                    {$infile}{outfile} )
                {

                    ## Extract OutFile
                    $outfile = $sample_info_href->{sample}{$sample_id}{recipe}
                      {$recipe}{$infile}{outfile};
                }
                if ( $sample_info_href->{sample}{$sample_id}{recipe}{$recipe}
                    {$infile}{path} )
                {

                    ( $outfile, $outdirectory ) =
                      fileparse( $sample_info_href->{sample}{$sample_id}{recipe}
                          {$recipe}{$infile}{path} );
                }

                ## Parses the RegExpHash structure to identify if the info is 1) Paragraf section(s) (both header and data line(s)); 2) Seperate data line.
                parse_regexp_hash_and_collect(
                    {
                        regexp_href         => $regexp_href,
                        qc_recipe_data_href => $qc_recipe_data_href,
                        qc_header_href      => $qc_header_href,
                        recipe              => $recipe,
                        outdirectory        => $outdirectory,
                        outfile             => $outfile,
                    }
                );

                ## Add extracted information to qc_data
                add_to_qc_data(
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
            }
        }
    }
    return;
}

sub parse_regexp_hash_and_collect {

## Function  : Parses the regexp hash structure to identify if the info is 1) Paragraf section(s) (both header and data line(s)); 2) Seperate data line.
## Returns   :
## Arguments : $outdirectory         => Recipes outdirectory
##           : $outfile              => Recipes outfile containing parameter to evaluate
##           : $qc_header_href       => Save header(s) in each outfile {REF}
##           : $qc_recipe_data_href  => Hash to save data in each outfile {REF}
##           : $recipe               => The recipe to examine
##           : $regexp_href          => Regexp hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $outdirectory;
    my $outfile;
    my $qc_header_href;
    my $qc_recipe_data_href;
    my $recipe;
    my $regexp_href;

    my $tmpl = {
        regexp_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$regexp_href,
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
        recipe => {
            defined     => 1,
            required    => 1,
            store       => \$recipe,
            strict_type => 1,
        },
        outdirectory => { store => \$outdirectory, strict_type => 1, },
        outfile      => { store => \$outfile,      strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Holds the current regexp
    my $regexp;

    ## Covers both whitespace and tab. Add other separators if required
    my @separators = ( qw{ \s+ ! }, q{,} );

    ## Find the actual regular expression(s) for each recipe that is used
  REG_EXP:
    for my $regexp_key ( keys %{ $regexp_href->{$recipe} } ) {

        if ( $regexp_key =~ /^header|header$/i ) {
## Detect if the outfile contains paragrafs/header info in the outfile i.e. data is formated as a paragraf with header(s) and line(s). "regexp_key" should either start with or end with "header". This section extracts the header line(s) for the entire outdata file. Necessary to assign correct data entry to header entry later (headers and data are saved in seperate hashes).

            ## Format outfile: Paragraf section
          PARAGRAPH:
            for my $regexp_header_key ( keys %{ $regexp_href->{$recipe}{$regexp_key} } ) {

                ## Regular expression used to collect paragraf header info
                $regexp =
                  $regexp_href->{$recipe}{$regexp_key}{$regexp_header_key};

                ## Loop through possible separators to seperate any eventual header elements
              SEPARATOR:
                for (
                    my $separator_element_counter = 0 ;
                    $separator_element_counter < scalar(@separators) ;
                    $separator_element_counter++
                  )
                {

                    ## Detect if the regexp key is a paragraf header and not paragraf data (header line and data line(s))
                    if ( $regexp_header_key =~ /^header|header$/i ) {

                        ## Collect paragraf header
                        @{ ${$qc_header_href}{$recipe}{$regexp_key}{$regexp_header_key} }
                          = split(
                            /$separators[$separator_element_counter]/,
                            `$regexp $outdirectory/$outfile`
                          );

                        #Then split should have been successful
                        if (
                            defined(
                                ${$qc_header_href}{$recipe}{$regexp_key}
                                  {$regexp_header_key}
                            )
                          )
                        {

                            ## Found correct separator - do not continue
                            last SEPARATOR;
                        }
                    }
                    else {
                        ## For paragraf data line(s)

                        ## Collect paragraf data
                        @{ $qc_recipe_data_href->{$recipe}{$regexp_key}
                              {$regexp_header_key} } = split(
                            /$separators[$separator_element_counter]/,
                            `$regexp $outdirectory/$outfile`
                              );

                        ## Then split should have been successful
                        if (
                            defined(
                                $qc_recipe_data_href->{$recipe}{$regexp_key}
                                  {$regexp_header_key}[1]
                            )
                          )
                        {

                            ## Found correct separator - do not continue
                            last SEPARATOR;
                        }
                    }
                }
            }
        }
        else {
            ##For info contained in Entry --> Value i.e. same line.

            ## The regular expression used to collect info
            $regexp = $regexp_href->{$recipe}{$regexp_key};

            ## Loop through possible separators
          SEPARATOR:
            for (
                my $separator_element_counter = 0 ;
                $separator_element_counter < scalar(@separators) ;
                $separator_element_counter++
              )
            {

                ## Collect data. Use regexp_key as element header
                @{ $qc_recipe_data_href->{$recipe}{$regexp_key} } = split(
                    /$separators[$separator_element_counter]/,
                    `$regexp $outdirectory/$outfile`
                );

                ## Then split should have been successful
                if ( defined( $qc_recipe_data_href->{$recipe}{$regexp_key}[1] ) ) {

                    ## Found correct separator do not continue
                    last SEPARATOR;
                }
            }
        }
    }
    return;
}

sub add_to_qc_data {

## Function  : Add to qc_data hash to enable write to yaml format
## Returns   :
## Arguments : $evaluate_plink_gender => Evaluate plink gender
##           : $infile                => Infile to recipe
##           : $recipe                => Recipe to examine
##           : $qc_data_href          => QCData hash {REF}
##           : $qc_header_href        => Save header(s) in each outfile {REF}
##           : $qc_recipe_data_href   => Hash to save data in each outfile {REF}
##           : $regexp_href           => RegExp hash {REF}
##           : $sample_id             => SampleID
##           : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $evaluate_plink_gender;
    my $infile;
    my $recipe;
    my $qc_data_href;
    my $qc_header_href;
    my $qc_recipe_data_href;
    my $regexp_href;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        evaluate_plink_gender => {
            allow       => [ undef, 0, 1 ],
            store       => \$evaluate_plink_gender,
            strict_type => 1,
        },
        infile => { store => \$infile, strict_type => 1, },
        recipe => {
            defined     => 1,
            required    => 1,
            store       => \$recipe,
            strict_type => 1,
        },
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
        sample_id => { store => \$sample_id, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  REGEXP:
    for my $regexp_key ( keys %{ $regexp_href->{$recipe} } ) {

        ## For info contained in entry --> Value i.e. same line
        if ( $regexp_key !~ /^header|header$/i ) {

            ## Enable seperation of writing array or key-->value in qc_data
            if ( scalar( @{ $qc_recipe_data_href->{$recipe}{$regexp_key} } ) == 1 ) {

                if ( ($sample_id) && ($infile) ) {

                    ## key-->value for sample_id
                    $qc_data_href->{sample}{$sample_id}{$infile}{$recipe}{$regexp_key} =
                      $qc_recipe_data_href->{$recipe}{$regexp_key}[0];
                }
                else {
                    ## Family level

                    ## key-->value for caseID
                    $qc_data_href->{recipe}{$recipe}{$regexp_key} =
                      $qc_recipe_data_href->{$recipe}{$regexp_key}[0];
                }

                ## Check gender for sample_id
                if ( $recipe eq q{chanjo_sexcheck} ) {

                    ## Array_ref
                    my $chanjo_sexcheck =
                      @{ $qc_recipe_data_href->{$recipe}{$regexp_key} }[0];

                    ## Check that assumed gender is supported by coverage on chrX and chrY
                    _chanjo_gender_check(
                        {
                            sample_info_href       => $sample_info_href,
                            qc_data_href           => $qc_data_href,
                            sample_id              => $sample_id,
                            infile                 => $infile,
                            chanjo_sexcheck_gender => $chanjo_sexcheck,
                        }
                    );
                }
            }
            else {
                ## Write array to qc_data

                for (
                    my $regexp_key_counter = 0 ;
                    $regexp_key_counter <
                    scalar( @{ $qc_recipe_data_href->{$recipe}{$regexp_key} } ) ;
                    $regexp_key_counter++
                  )
                {

                    if ( ($sample_id) && ($infile) ) {

                        $qc_data_href->{sample}{$sample_id}{$infile}{$recipe}
                          {$regexp_key}[$regexp_key_counter] =
                          $qc_recipe_data_href->{$recipe}{$regexp_key}
                          [$regexp_key_counter];

                    }
                    else {

                        $qc_data_href->{recipe}{$recipe}{$regexp_key}[$regexp_key_counter]
                          = $qc_recipe_data_href->{$recipe}{$regexp_key}
                          [$regexp_key_counter];
                    }
                    ## Check gender for sample_id
                    if (   $recipe eq q{plink_sexcheck}
                        && $evaluate_plink_gender )
                    {

                        ## Array ref
                        my @sexchecks = split( q{:},
                            @{ $qc_recipe_data_href->{$recipe}{$regexp_key} }
                              [$regexp_key_counter] );

                        ## Check that assumed gender is supported by variants on chrX and chrY
                        plink_gender_check(
                            {
                                sample_info_href          => $sample_info_href,
                                qc_data_href              => $qc_data_href,
                                sample_id_ref             => \$sexchecks[0],
                                plink_sexcheck_gender_ref => \$sexchecks[1],
                            }
                        );
                    }
                }
                if (
                    defined(
                        $qc_data_href->{recipe}{relation_check}{sample_relation_check}
                    )
                    && (
                        defined( $qc_data_href->{recipe}{pedigree_check}{sample_order} ) )
                  )
                {

                    relation_check(
                        {
                            sample_info_href        => $sample_info_href,
                            qc_data_href            => $qc_data_href,
                            relationship_values_ref => \@{
                                $qc_data_href->{recipe}{relation_check}
                                  {sample_relation_check}
                            },
                            sample_orders_ref =>
                              \@{ $qc_data_href->{recipe}{pedigree_check}{sample_order} },
                        }
                    );
                }
            }
        }
        else {
            ## Paragraf data i.e. header and subsequent data lines

          HEADER_INFO:
            for
              my $regexp_header_key ( keys %{ ${$qc_header_href}{$recipe}{$regexp_key} } )
            {

              PARAGRAPH_KEYS:
                for
                  my $regexp_key_header ( keys %{ $regexp_href->{$recipe}{$regexp_key} } )
                {
                    ## All paragraf keys (header and data line(s))

                    ## Detect if the regexp id for headers and not data.
                    if ( $regexp_key_header !~ /^header|header$/i ) {

                        ## For all collected headers
                        for (
                            my $qc_headers_counter = 0 ;
                            $qc_headers_counter < scalar(
                                @{
                                    ${$qc_header_href}{$recipe}{$regexp_key}
                                      {$regexp_header_key}
                                }
                            ) ;
                            $qc_headers_counter++
                          )
                        {

                            if ( ($sample_id) && ($infile) ) {

                                ## Add to qc_data using header element[X] --> data[X] to correctly position elements in qc_data hash
                                $qc_data_href->{sample}{$sample_id}{$infile}
                                  {$recipe}{$regexp_header_key}{$regexp_key_header}
                                  { ${$qc_header_href}{$recipe}{$regexp_key}
                                      {$regexp_header_key}[$qc_headers_counter] } =
                                  $qc_recipe_data_href->{$recipe}
                                  {$regexp_key}{$regexp_key_header}[$qc_headers_counter];
                            }
                            else {

                                ## Add to qc_data using header element[X] --> data[X] to correctly position elements in qc_data hash
                                $qc_data_href->{$recipe}{$regexp_header_key}
                                  {$regexp_key_header}
                                  { ${$qc_header_href}{$recipe}{$regexp_key}
                                      {$regexp_header_key}[$qc_headers_counter] } =
                                  $qc_recipe_data_href->{$recipe}
                                  {$regexp_key}{$regexp_key_header}[$qc_headers_counter];

                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

sub define_evaluate_metric {

## Function  : Sets recipes and recipe metrics and thresholds to be evaluated
## Returns   :
## Arguments : $sample_info_href => Info on samples and case hash {REF}
##           : $sample_id        => Sample ID

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $sample_id;

    my $tmpl = {
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $evaluate_metric{$sample_id}{bamstats}{percentage_mapped_reads}{lt} = 95;
    $evaluate_metric{$sample_id}{collecthsmetrics}{PCT_TARGET_BASES_10X}{lt} =
      0.95;
    $evaluate_metric{$sample_id}{collectmultiplemetrics}{PCT_PF_READS_ALIGNED}{lt} = 0.95;
    $evaluate_metric{$sample_id}{collectmultiplemetrics}{PCT_ADAPTER}{gt} =
      0.0005;
    $evaluate_metric{$sample_id}{markduplicates}{fraction_duplicates}{gt} = 0.2;

    if ( exists $sample_info_href->{sample}{$sample_id}{expected_coverage} ) {

        $evaluate_metric{$sample_id}{collecthsmetrics}{MEAN_TARGET_COVERAGE}{lt} =
          $sample_info_href->{sample}{$sample_id}{expected_coverage};
    }

    $evaluate_metric{mendel}{fraction_of_errors}{gt} = 0.06;
    $evaluate_metric{father}{fraction_of_common_variants}{lt} =
      0.55;

    return;
}

sub evaluate_qc_parameters {

## Function : Evaluate parameters to detect parameters falling below threshold
## Returns  :
## Arguments: $qc_data_href         => QC data hash {REF}
##          : $evaluate_metric_href => Hash for metrics to evaluate

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $evaluate_metric_href;

    my $tmpl = {
        qc_data_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$qc_data_href
        },
        evaluate_metric_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$evaluate_metric_href
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  RECIPE:
    for my $recipe ( keys %{ $qc_data_href->{recipe} } ) {

        ## Recipe to be evaluated
        if ( exists $evaluate_metric_href->{$recipe} ) {

          METRIC:
            for my $metric ( keys %{ $qc_data_href->{recipe}{$recipe} } ) {

              FAMILY_LEVEL:
                if ( exists $evaluate_metric_href->{$recipe}{$metric} ) {

                    check_metric(
                        {
                            qc_data_href => $qc_data_href,
                            reference_metric_href =>
                              $evaluate_metric_href->{$recipe}{$metric},
                            recipe          => $recipe,
                            metric          => $metric,
                            qc_metric_value => $qc_data_href->{recipe}{$recipe}{$metric},
                        }
                    );

                }
            }
        }
    }

  SAMPLE_LEVEL:
    for my $sample_id ( keys %{ $qc_data_href->{sample} } ) {

      INFILE:
        for my $infile ( keys %{ $qc_data_href->{sample}{$sample_id} } ) {

            ## Special case skip evaluation
            next INFILE if ( $infile =~ /evaluation/ );

            ## Special case do not evaluate fastq files with Undetermined in file name
            next INFILE if ( $infile =~ /Undetermined/i );

            ## Special case
            if ( $infile =~ /relation_check/ ) {

                if ( $qc_data_href->{sample}{$sample_id}{$infile} ne q{PASS} ) {

                    my $status =
                        q{Status:}
                      . $infile . q{:}
                      . $qc_data_href->{sample}{$sample_id}{$infile};
                    ## Add to QC data at case level
                    push @{ $qc_data_href->{evaluation}{$infile} }, $status;
                }
                next INFILE;
            }

          RECIPE:
            for my $recipe ( keys %{ $qc_data_href->{sample}{$sample_id}{$infile} } ) {

                ## Recipe to be evaluated
                if ( exists $evaluate_metric_href->{$sample_id}{$recipe} ) {

                  METRIC:
                    for my $metric (
                        keys %{ $evaluate_metric_href->{$sample_id}{$recipe} } )
                    {

                        if (
                            exists $qc_data_href->{sample}{$sample_id}{$infile}
                            {$recipe}{$metric} )
                        {

                            check_metric(
                                {
                                    qc_data_href => $qc_data_href,
                                    reference_metric_href =>
                                      $evaluate_metric_href->{$sample_id}
                                      {$recipe}{$metric},
                                    recipe => $recipe,
                                    metric => $metric,
                                    qc_metric_value =>
                                      $qc_data_href->{sample}{$sample_id}
                                      {$infile}{$recipe}{$metric},
                                }
                            );
                        }
                        else {

                            if (
                                exists $qc_data_href->{sample}{$sample_id}
                                {$infile}{$recipe}{header} )
                            {

                              HEADER:
                                for my $data_header (
                                    keys %{
                                        $qc_data_href->{sample}{$sample_id}
                                          {$infile}{$recipe}{header}
                                    }
                                  )
                                {

                                    if (
                                        exists $qc_data_href->{sample}
                                        {$sample_id}{$infile}{$recipe}{header}
                                        {$data_header}{$metric} )
                                    {

                                        check_metric(
                                            {
                                                qc_data_href => $qc_data_href,
                                                reference_metric_href =>
                                                  $evaluate_metric_href->{$sample_id}
                                                  {$recipe}{$metric},
                                                recipe => $recipe,
                                                metric => $metric,
                                                qc_metric_value =>
                                                  $qc_data_href->{sample}
                                                  {$sample_id}{$infile}{$recipe}{header}
                                                  {$data_header}{$metric},
                                            }
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

sub check_metric {

## Function : Check and add result of check if below threshold
## Returns  :
## Arguments: $qc_data_href           => QCData hash {REF}
##          : $reference_metric_href  => Metrics to evaluate
##          : $recipe                 => The recipe to examine
##          : $metric                 => Metric to evaluate
##          : $qc_metric_value        => Qc metric value

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $reference_metric_href;
    my $recipe;
    my $metric;
    my $qc_metric_value;

    my $tmpl = {
        qc_data_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$qc_data_href
        },
        reference_metric_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$reference_metric_href
        },
        recipe => { required => 1, defined => 1, strict_type => 1, store => \$recipe },
        metric => { required => 1, defined => 1, strict_type => 1, store => \$metric },
        qc_metric_value => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$qc_metric_value
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $status = q{FAILED:};

    if ( exists $reference_metric_href->{lt} ) {

        ## Determine status - if lower than add to hash. otherwise PASS and do not include
        if ( $qc_metric_value < $reference_metric_href->{lt} ) {

            $status .= $recipe . q{_} . $metric . q{:} . $qc_metric_value;
            push @{ $qc_data_href->{evaluation}{$recipe} }, $status;
        }
    }

    if ( exists $reference_metric_href->{gt} ) {

        ## Determine status - if greater than add to hash. otherwise PASS and do not include
        if ( $qc_metric_value > $reference_metric_href->{gt} ) {

            $status .= $recipe . q{_} . $metric . q{:} . $qc_metric_value;
            push @{ $qc_data_href->{evaluation}{$recipe} }, $status;
        }
    }
    return;
}

sub relation_check {

## Function : Uses the .mibs file produced by PLINK to test if case members are indeed related.
## Returns  :
## Arguments: $sample_info_href        => Info on samples and case hash {REF}
##          : $qc_data_href            => QCData hash {REF}
##          : $sample_orders_ref       => The sample order so that correct estimation can be connected to the correct sample_ids {REF}
##          : $relationship_values_ref => All relationship estimations {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $qc_data_href;
    my $relationship_values_ref;
    my $sample_orders_ref;

    my $tmpl = {
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        qc_data_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$qc_data_href
        },
        relationship_values_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$relationship_values_ref
        },
        sample_orders_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_orders_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

#Stores case relations and pairwise comparisons case{$sample_id}{$sample_id}["column"] -> [pairwise]
    my %case;
    my $sample_id_counter  = 0;
    my $incorrect_relation = 0;
    my @pairwise_comparisons;

    ##Copy array to avoid removing actual values in later splice
    my @relationship_values = @{$relationship_values_ref};

    ## Splice all relationship extimations from regexp into pairwise comparisons calculated for each sample_id
    for (
        my $realtionship_counter = 0 ;
        $realtionship_counter < scalar(@relationship_values) ;
        $realtionship_counter++
      )
    {

        ## Splices array into each sample_ids line
        my @pairwise_comparisons =
          splice( @relationship_values, 0, scalar( @{$sample_orders_ref} ) );

        ## All columns in .mibs file
        for ( my $column = 0 ; $column < scalar(@$sample_orders_ref) ; $column++ ) {

            ## Store sample_id, case membersID (including self) and each pairwise comparison. Uses array for to accomodate sibling info.
            push(
                @{
                    $case{ $sample_orders_ref->[$sample_id_counter] }
                      { $sample_orders_ref->[$column] }
                },
                $pairwise_comparisons[$column]
            );
        }
        $sample_id_counter++;
    }
    ## Father_id for the case
    my $father_id = "YYY";

    ## Mother_id for the case
    my $mother_id = "XXX";

    ## Collect father and mother id
  SAMPLE_ID:
    for my $sample_id ( keys %case ) {    #For all sample_ids

        ## Currently only 1 father or Mother per pedigree is supported

        ## Save father_id if not 0
        if ( $sample_info_href->{sample}{$sample_id}{father} ne 0 ) {

            $father_id = $sample_info_href->{sample}{$sample_id}{father};
        }

        ## Save mother_id if not 0
        if ( $sample_info_href->{sample}{$sample_id}{mother} ne 0 ) {

            $mother_id = $sample_info_href->{sample}{$sample_id}{mother};
        }
    }

  SAMPLE_ID:
    for my $sample_id ( keys %case ) {

      MEMBER:
        for my $members ( keys %{ $case{$sample_id} } ) {
            ## For every relation within case (mother/father/child)

          RELATIVES:
            for (
                my $members_count = 0 ;
                $members_count < scalar( @{ $case{$sample_id}{$members} } ) ;
                $members_count++
              )
            {
                ## Necessary for siblings

                ## Should only hit self
                if ( $case{$sample_id}{$members}[$members_count] == 1 ) {

                    if ( $sample_id eq $members ) {

#print "Self: ".$sample_id,"\t", $members, "\t", $case{$sample_id}{$members}[$members_count], "\n";
                    }
                    else {

                        $incorrect_relation++;
                        $qc_data_href->{sample}{$sample_id}{relation_check} =
                          "FAIL: Duplicated sample?;";

#print  "Incorrect should be self: ".$sample_id,"\t", $members, "\t", $case{$sample_id}{$members}[$members_count], "\n";
                    }
                }
                elsif ( $case{$sample_id}{$members}[$members_count] >= 0.70 )
                { #Should include parent to child and child to siblings unless inbreed parents

                    if (
                        ( ( $sample_id ne $father_id ) && ( $sample_id ne $mother_id ) )
                        || (   ( $members ne $father_id )
                            && ( $members ne $mother_id ) )
                      )
                    {    #Correct
                         #print "Parent-to-child or child-to-child: ".$sample_id,"\t", $members, "\t", $case{$sample_id}{$members}[$members_count], "\n";
                    }
                    else {

                        $incorrect_relation++;
                        $qc_data_href->{sample}{$sample_id}{relation_check} =
                          "FAIL: Parents related?;";

#print "Incorrect: ".$sample_id,"\t", $members, "\t", $case{$sample_id}{$members}[$members_count], "\n";
                    }
                }
                elsif ( $case{$sample_id}{$members}[$members_count] < 0.70 )
                {        #Parents unless inbreed

                    if (   ( $sample_id eq $father_id )
                        && ( $members eq $mother_id ) )
                    {

#print "Parents: ".$sample_id,"\t", $members, "\t", $case{$sample_id}{$members}[$members_count], "\n";
                    }
                    elsif (( $sample_id eq $mother_id )
                        && ( $members eq $father_id ) )
                    {

#print "Parents: ".$sample_id,"\t", $members, "\t", $case{$sample_id}{$members}[$members_count], "\n";
                    }
                    else {

                        $incorrect_relation++;
                        $qc_data_href->{sample}{$sample_id}{relation_check} =
                          "FAIL:" . $sample_id . " not related to " . $members . ";";

#print "Incorrect: ".$sample_id,"\t", $members, "\t", $case{$sample_id}{$members}[$members_count], "\n";
                    }
                }
            }
        }
        if ( $incorrect_relation == 0 ) {

            $qc_data_href->{sample}{$sample_id}{relation_check} = "PASS";
        }
    }
    return;
}

sub _chanjo_gender_check {

## Function : Checks that the gender predicted by chanjo_sexcheck is confirmed in the pedigee for the sample
## Returns  :
## Arguments: $chanjo_sexcheck_gender => Chanjo calculated gender
##          : $infile                 => Infile {REF}
##          : $qc_data_href           => QCData hash {REF}
##          : $sample_id              => Sample ID
##          : $sample_info_href       => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chanjo_sexcheck_gender;
    my $infile;
    my $qc_data_href;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        chanjo_sexcheck_gender => {
            defined     => 1,
            required    => 1,
            store       => \$chanjo_sexcheck_gender,
            strict_type => 1,
        },
        infile => {
            defined     => 1,
            required    => 1,
            store       => \$infile,
            strict_type => 1,
        },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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

    ## Get sample id sex
    my $sample_id_sex = $sample_info_href->{sample}{$sample_id}{sex};

    ## Female
    if (   $chanjo_sexcheck_gender eq q{female}
        && $sample_id_sex =~ /2|female/ )
    {

        $qc_data_href->{sample}{$sample_id}{$infile}{gender_check} =
          q{PASS};
    }
    elsif ($chanjo_sexcheck_gender eq q{male}
        && $sample_id_sex =~ /1|^male/ )
    {
        ## Male

        $qc_data_href->{sample}{$sample_id}{$infile}{gender_check} =
          q{PASS};
    }
    elsif ( $sample_id_sex =~ /other|unknown/ ) {
        ## Other|Unknown

        $qc_data_href->{sample}{$sample_id}{$infile}{gender_check} =
          q{PASS};
    }
    else {

        $qc_data_href->{sample}{$sample_id}{$infile}{gender_check} =
          q{FAIL};
    }
    return;
}

sub plink_gender_check {

##Function : Checks that the gender predicted by Plink sexcheck is confirmed in the pedigee for the sample
##Returns  :
##Arguments: $sample_info_href          => Info on samples and case hash {REF}
##         : $qc_data_href              => QCData hash {REF}
##         : $sample_id_ref             => SampleID {REF}
##         : $plink_sexcheck_gender_ref => Plink calculated gender {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;
    my $qc_data_href;
    my $sample_id_ref;
    my $plink_sexcheck_gender_ref;

    my $tmpl = {
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        qc_data_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$qc_data_href
        },
        sample_id_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$sample_id_ref
        },
        plink_sexcheck_gender_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$plink_sexcheck_gender_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Alias
    my $sample_id_sex_ref = \$sample_info_href->{sample}{$$sample_id_ref}{sex};

    ## Female
    if (   ( $$plink_sexcheck_gender_ref eq q{2} )
        && ( $$sample_id_sex_ref =~ /2|female/ ) )
    {

        push @{ $qc_data_href->{recipe}{plink_gender_check} }, $$sample_id_ref . q{:PASS};
    }
    elsif (( $$plink_sexcheck_gender_ref eq q{1} )
        && ( $$sample_id_sex_ref =~ /1|^male/ ) )
    {
        ## Male

        push @{ $qc_data_href->{recipe}{plink_gender_check} }, $$sample_id_ref . q{:PASS};
    }
    elsif ( $$sample_id_sex_ref =~ /other|unknown/ ) {
        ## Other|Unknown

        push @{ $qc_data_href->{recipe}{plink_gender_check} }, $$sample_id_ref . q{:PASS};
    }
    else {

        push @{ $qc_data_href->{recipe}{plink_gender_check} }, $$sample_id_ref . q{:FAIL};
    }
    return;
}

sub regexp_to_yaml {

##Function : Write default regexp to YAML
##Returns  :
##Arguments: $print_regexp_outfile => File to print regexp to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $print_regexp_outfile;

    my $tmpl = {
        print_regexp_outfile => {
            defined     => 1,
            required    => 1,
            store       => \$print_regexp_outfile,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %regexp;

    #Add to %regexp to enable print in YAML
    $regexp{fastqc}{version} =
      q?perl -nae' if ($_=~/##FastQC\\s+(\\S+)/) {print $1;last;}' ?
      ;    #Collect FastQC version

    $regexp{fastqc}{encoding} =
q?perl -nae' if ($_=~/Encoding\s+(\S+\s\S+\s\S+\s\S+|\S+\s\S+)/) { my $encoding = $1;$encoding=~s/\s/\_/g; print $encoding;last;}' ?
      ;    #Collect Encoding

    $regexp{fastqc}{sequence_length} =
      q?perl -nae' if ($_=~/Sequence length\s(\d+)/) {print $1;last;}' ?
      ;    #Collect Sequence length

    $regexp{fastqc}{total_number_of_reads} =
      q?perl -nae' if ($_=~/Total Sequences\s(\d+)/) {print $1;last;}' ?
      ;    #Collect Total sequences

    $regexp{fastqc}{gc} =
      q?perl -nae' if ($_=~/%GC\s(\d+)/) {print $1;last;}' ?;    #Collect GC content

    $regexp{fastqc}{sequence_duplication} =
      q?perl -nae' if ($_=~/#Total Duplicate Percentage\s+(\d+.\d)/) {print $1;last;}' ?
      ;    #Collect Sequence duplication level

    $regexp{fastqc}{basic_statistics} =
      q?perl -nae' if ($_=~/>>Basic Statistics\s+(\S+)/) {print $1;last;}' ?
      ;    #Collect Basic Statistics

    $regexp{fastqc}{per_base_sequence_quality} =
      q?perl -nae' if ($_=~/>>Per base sequence quality\s+(\S+)/) {print $1;last;}' ?
      ;    #Collect Per base sequence quality

    $regexp{fastqc}{per_sequence_quality_scores} =
      q?perl -nae' if ($_=~/>>Per sequence quality scores\s+(\S+)/) {print $1;last;}' ?
      ;    #Collect Per sequence quality scores

    $regexp{fastqc}{per_base_sequence_content} =
      q?perl -nae' if ($_=~/>>Per base sequence content\s+(\S+)/) {print $1;last;}' ?
      ;    #Collect Per base sequence content

    $regexp{fastqc}{per_base_gc_content} =
      q?perl -nae' if ($_=~/>>Per base GC content\s+(\S+)/) {print $1;last;}' ?
      ;    #Collect Per base GC content

    $regexp{fastqc}{per_sequence_gc_content} =
      q?perl -nae' if ($_=~/>>Per sequence GC content\s+(\S+)/) {print $1;last;}' ?
      ;    #Collect Per sequence GC content

    $regexp{fastqc}{per_base_n_content} =
      q?perl -nae' if ($_=~/>>Per base N content\s+(\S+)/) {print $1;last;}' ?
      ;    #Collect Per base N content

    $regexp{fastqc}{sequence_duplication_levels} =
      q?perl -nae' if ($_=~/>>Sequence Duplication Levels\s+(\S+)/) {print $1;last;}' ?
      ;    #Collect Sequence Duplication Levels

    $regexp{fastqc}{overrepresented_sequences} =
      q?perl -nae' if ($_=~/>>Overrepresented sequences\s+(\S+)/) {print $1;last;}' ?
      ;    #Collect Overrepresented sequences

    $regexp{fastqc}{kmer_content} =
      q?perl -nae' if ($_=~/>>Kmer Content\s+(\S+)/) {print $1;last;}' ?
      ;    #Collect Kmer Content

    $regexp{bamstats}{percentage_mapped_reads} =
      q?perl -nae 'if($_=~/percentage mapped reads:\s+(\S+)/) {print $1;last}' ?
      ;    #Collect % mapped reads from BAM alignment

    $regexp{bamstats}{raw_total_sequences} =
      q?perl -nae 'if($_=~/raw total sequences:\s+(\S+)/) {print $1;last}' ?
      ;    #Collect raw total sequences from BAM alignment

    $regexp{bamstats}{reads_mapped} =
      q?perl -nae 'if($_=~/reads mapped:\s+(\S+)/) {print $1;last}' ?
      ;    #Collect reads mapped from BAM alignment

    $regexp{chanjo_sexcheck}{gender} =
      q?perl -nae 'if( ($F[0]!~/^#/) && ($F[2] =~/\S+/) ) {print $F[2];}' ?
      ;    #Collect gender from chanjo_sexcheck

    $regexp{pedigree_check}{sample_order} =
q?perl -nae 'if ($_=~/^#CHROM/) {chomp $_; my @line = split(/\t/,$_); for (my $sample=9;$sample<scalar(@line);$sample++) { print $line[$sample], "\t";}last;}' ?
      ; #Collect sample order from vcf file used to create ".ped", ".map" and hence ".mibs".

    $regexp{inbreeding_factor}{sample_inbreeding_factor} =
q?perl -nae 'my @inbreedingFactor; if ($. > 1) {my @temp = split(/\s/,$_);push(@inbreedingFactor, $F[0].":".$F[5]); print $inbreedingFactor[0], "\t"; }' ?;

    $regexp{plink_sexcheck}{sample_sexcheck} =
q?perl -nae 'my @sexCheckFactor; if ($. > 1) {my @temp = split(/\s+/,$_);push(@sexCheckFactor,$temp[2].":".$temp[4]); print $sexCheckFactor[0], "\t"; }' ?;

    $regexp{relation_check}{sample_relation_check} =
      q?perl -nae 'print $_;' ?;    #Note will return whole file

    $regexp{markduplicates}{fraction_duplicates} =
      q?perl -nae 'if($_=~/Fraction Duplicates\: (\S+)/) {print $1;}' ?
      ;                             #Collect fraction duplicates

    $regexp{collecthsmetrics}{header_info}{header} =
      q?perl -nae' if ($_ =~/^BAIT_SET/ ) {print $_;last;}' ?
      ;                             #Note return whole line (header)

    $regexp{collecthsmetrics}{header_info}{data} =
      q?perl -nae' if ( ($. ==8) && ($_ =~/(\S+)/) ) {print $_;last;}' ?
      ;    #Note return whole line and only look at line 8, where the data action is

    $regexp{collectmultiplemetrics}{header_info}{header} =
      q?perl -nae' if ($_ =~/^CATEGORY/ ) {print $_;last;}' ?
      ;    #Note return whole line (header)

    $regexp{collectmultiplemetrics}{header_info}{first_of_pair} =
      q?perl -nae' if ($_ =~/^FIRST_OF_PAIR/ ) {print $_;last;}' ?
      ;    #Note return whole line (FIRST_OF_PAIR)

    $regexp{collectmultiplemetrics}{header_info}{second_of_pair} =
      q?perl -nae' if ($_ =~/^SECOND_OF_PAIR/ ) {print $_;last;}' ?
      ;    #Note return whole line (SECOND_OF_PAIR)

    $regexp{collectmultiplemetrics}{header_info}{pair} =
      q?perl -nae' if ($_ =~/^PAIR/ ) {print $_;last;}'  ?; #Note return whole line (PAIR)

    $regexp{collectmultiplemetricsinsertsize}{header_info}{header} =
      q?perl -nae' if ($_ =~/^MEDIAN_INSERT_SIZE/ ) {print $_;last;}' ?
      ;    #Note return whole line (header)

    $regexp{collectmultiplemetricsinsertsize}{header_info}{data} =
      q?perl -nae' if ( ($. ==8) && ($_ =~/(\S+)/) ) {print $_;last;}' ?
      ;    #Note return whole line and only look at line 8, where the data action is

    $regexp{variantevalall}{comp_overlap_header}{comp_overlap_header} =
      q?perl -nae' if ($_ =~/^CompOverlap\s+CompRod/ ) {print $_;last;}' ?
      ;    #Note return whole line (header)

    $regexp{variantevalall}{comp_overlap_header}{comp_overlap_data_all} =
q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/all/) && ($_ =~/none/)) {print $_;last;}' ?
      ;    #Note return whole line

    $regexp{variantevalall}{comp_overlap_header}{comp_overlap_data_known} =
      q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/known\s/) ) {print $_;last;}' ?
      ;    #Note return whole line

    $regexp{variantevalall}{comp_overlap_header}{comp_overlap_data_novel} =
      q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/novel\s/) ) {print $_;last;}' ?
      ;    #Note return whole line

    $regexp{variantevalall}{count_variants_header}{count_variants_header} =
      q?perl -nae' if ($_ =~/^CountVariants\s+CompRod/ ) {print $_;last;}' ?
      ;    #Note return whole line (header)
    $regexp{variantevalall}{count_variants_header}{count_variants_data_all} =
      q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/all\s/) ) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{count_variants_header}{count_variants_data_known} =
      q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/known\s/) ) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{count_variants_header}{count_variants_data_novel} =
      q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/novel\s/) ) {print $_;last;}' ?
      ;    #Note return whole line

    $regexp{variantevalall}{indel_summary_header}{indel_summary_header} =
      q?perl -nae' if ($_ =~/^IndelSummary\s+CompRod/ ) {print $_;last;}' ?
      ;    #Note return whole line (header)
    $regexp{variantevalall}{indel_summary_header}{indel_summary_data_all} =
      q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{indel_summary_header}{indel_summary_data_known} =
      q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{indel_summary_header}{indel_summary_data_novel} =
      q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?
      ;    #Note return whole line

    $regexp{variantevalall}{multiallelic_summary_header}{multiallelic_summary_header} =
      q?perl -nae' if ($_ =~/^MultiallelicSummary\s+CompRod/ ) {print $_;last;}' ?
      ;    #Note return whole line (header)
    $regexp{variantevalall}{multiallelic_summary_header}{multiallelic_summary_data_all} =
q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{multiallelic_summary_header}{multiallelic_summary_data_known}
      = q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{multiallelic_summary_header}{multiallelic_summary_data_novel}
      = q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?
      ;    #Note return whole line

    $regexp{variantevalall}{titv_variant_evaluator_header}{titv_variant_evaluator_header}
      = q?perl -nae' if ($_ =~/^TiTvVariantEvaluator\s+CompRod/ ) {print $_;last;}' ?
      ;    #Note return whole line (header)
    $regexp{variantevalall}{titv_variant_evaluator_header}
      {titv_variant_evaluator_data_all} =
q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/all\s/) ) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{titv_variant_evaluator_header}
      {titv_variant_evaluator_data_known} =
q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/known\s/) ) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{titv_variant_evaluator_header}
      {titv_variant_evaluator_data_novel} =
q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/novel\s/) ) {print $_;last;}' ?
      ;    #Note return whole line

    $regexp{variantevalall}{validation_report_header}{validation_report_header} =
      q?perl -nae' if ($_ =~/^ValidationReport\s+CompRod/ ) {print $_;last;}' ?
      ;    #Note return whole line (header)
    $regexp{variantevalall}{validation_report_header}{validation_report_data_all} =
q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/all\s/) && ($_ =~/none\s/)) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{validation_report_header}{validation_report_data_known} =
q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/known\s/) ) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{validation_report_header}{validation_report_data_novel} =
q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/novel\s/) ) {print $_;last;}' ?
      ;    #Note return whole line

    $regexp{variantevalall}{variant_summary_header}{variant_summary_header} =
      q?perl -nae' if ($_ =~/^VariantSummary\s+CompRod/ ) {print $_;last;}' ?
      ;    #Note return whole line (header)
    $regexp{variantevalall}{variant_summary_header}{variant_summary_data_all} =
      q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{variant_summary_header}{variant_summary_data_known} =
      q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?
      ;    #Note return whole line
    $regexp{variantevalall}{variant_summary_header}{variant_summary_data_novel} =
      q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?
      ;    #Note return whole line

    $regexp{variantevalexome} = $regexp{variantevalall};

    $regexp{genmod}{version} =
q?perl -nae 'if($_=~/##Software=<ID=genmod,Version=(\d+.\d+.\d+|\d+.\d+)/) {print $1;last;}' ?
      ;    #Collect Genmod version

    $regexp{snpeff}{version} =
q?perl -nae 'if($_=~/##SnpSiftVersion=\"(.+),/) {my $ret=$1; $ret=~s/\s/_/g;print $ret;last;}' ?
      ;    #Collect SnpEff version

    $regexp{varianteffectpredictor}{version} =
      q?perl -nae 'if($_=~/##VEP="(\w+)"/) {print $1;last;}' ?
      ;    #Collect varianteffectpredictor version

    $regexp{varianteffectpredictor}{cache} =
      q?perl -nae 'if($_=~/##VEP=\w+\s+cache=(\S+)/) {print $1;last;}' ?
      ;    #Collect varianteffectpredictor cache directory

    $regexp{varianteffectpredictor}{polyphen} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/polyphen=(\S+)/) {print $1;last;}' ?
      ;    #Collect varianteffectpredictor polyPhen version

    $regexp{varianteffectpredictor}{sift} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/sift=sift(\S+)/) {print $1;last;}' ?
      ;    #Collect varianteffectpredictor sift version

    $regexp{varianteffectpredictor}{gene_build} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/genebuild=(\S+)/) {print $1;last;}' ?
      ;    #Collect varianteffectpredictor geneBuild

    $regexp{varianteffectpredictor}{assembly} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/assembly=(\S+)/) {print $1;last;}' ?
      ;    #Collect varianteffectpredictor assembly

    $regexp{varianteffectpredictor}{hgmd_public} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/HGMD-PUBLIC=(\S+)/) {print $1;last;}' ?
      ;    #Collect varianteffectpredictor HGMD-PUBLIC version

    $regexp{varianteffectpredictor}{reg_build} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/regbuild=(\S+)/) {print $1;last;}' ?
      ;    #Collect varianteffectpredictor regbuild version

    $regexp{varianteffectpredictor}{gencode} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/gencode=\S+\s+(\d+)/) {print $1;last;}' ?
      ;    #Collect varianteffectpredictor gencode version

    $regexp{vcfparser}{version} =
q?perl -nae 'if($_=~/##Software=<ID=vcfParser.pl,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?
      ;    #Collect vcfparser version

    $regexp{bwa}{version} =
      q?perl -nae 'if($_=~/\[main\]\sVersion:\s(\S+)/) {print $1;last;}' ?
      ;    #Collect Bwa version

    $regexp{chanjo}{version} =
      q?perl -nae 'if($_=~/version\s(\d+.\d+.\d+)/) {print $1;last;}' ?
      ;    #Collect Chanjo version

    $regexp{vt}{version} =
      q?perl -nae 'if($_=~/decompose\sv(\S+)/) {print $1;last;}' ?;    #Collect vt version

    $regexp{samtools}{version} =
      q?perl -nae 'if($_=~/samtoolsVersion=(\S+)/) {print $1;last;}' ?
      ;    #Collect Samtools version

    $regexp{bcftools}{version} =
      q?perl -nae 'if($_=~/bcftools_\w+Version=(\S+)/) {print $1;last;}' ?
      ;    #Collect Bcftools version

    $regexp{freebayes}{version} =
      q?perl -nae 'if($_=~/source=freeBayes\s(\S+)/) {print $1;last;}' ?
      ;    #Collect Freebayes version

    $regexp{delly}{version} =
      q?perl -nae 'if($_=~/SVMETHOD=EMBL\.DELLY(v\d+\.\d+\.\d+)/) {print $1;last }' ?
      ;    #Collect Delly version

    $regexp{manta}{version} =
      q?perl -nae 'if($_=~/GenerateSVCandidates\s+(\S+)/) {print $1;last}' ?
      ;    #Collect Manta version

    $regexp{sv_combinevariantcallsets}{vcfanno} =
      q?perl -nae 'if($_=~/vcfanno\sversion\s(\S+)/) {print $1;last;}' ?
      ;    #Collect SVVCFAnno version

    $regexp{sv_varianteffectpredictor}{version} =
      q?perl -nae 'if($_=~/##VEP="(\w+)"/) {print $1;last;}' ?
      ;    #Collect sv_varianteffectpredictor version

    $regexp{sv_varianteffectpredictor}{cache} =
      q?perl -nae 'if($_=~/##VEP=\w+\s+cache=(\S+)/) {print $1;last;}' ?
      ;    #Collect sv_varianteffectpredictor cache directory

    $regexp{sv_varianteffectpredictor}{polyphen} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/polyphen=(\S+)/) {print $1;last;}' ?
      ;    #Collect sv_varianteffectpredictor polyPhen version

    $regexp{sv_varianteffectpredictor}{sift} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/sift=sift(\S+)/) {print $1;last;}' ?
      ;    #Collect sv_varianteffectpredictor sift version

    $regexp{sv_varianteffectpredictor}{gene_build} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/genebuild=(\S+)/) {print $1;last;}' ?
      ;    #Collect sv_varianteffectpredictor geneBuild

    $regexp{sv_varianteffectpredictor}{assembly} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/assembly=(\S+)/) {print $1;last;}' ?
      ;    #Collect sv_varianteffectpredictor assembly

    $regexp{sv_varianteffectpredictor}{hgmd_public} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/HGMD-PUBLIC=(\S+)/) {print $1;last;}' ?
      ;    #Collect sv_varianteffectpredictor HGMD-PUBLIC version

    $regexp{sv_varianteffectpredictor}{reg_build} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/regbuild=(\S+)/) {print $1;last;}' ?
      ;    #Collect sv_varianteffectpredictor regbuild version

    $regexp{sv_varianteffectpredictor}{gencode} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/gencode=\S+\s+(\d+)/) {print $1;last;}' ?
      ;    #Collect sv_varianteffectpredictor gencode version

    $regexp{sv_vcfparser}{version} =
q?perl -nae 'if($_=~/##Software=<ID=vcfParser.pl,Version=(\d+.\d+.\d+)/) {print $1;last;} else { if($_=~/#CHROM/) {last;} }' ?
      ;    #Collect sv_vcfparser version

    $regexp{sv_genmod}{version} =
q?perl -nae 'if($_=~/##Software=<ID=genmod,Version=(\d+.\d+.\d+|\d+.\d+)/) {print $1;last;} else { if($_=~/#CHROM/) {last;} } ' ?
      ;    #Collect SVGenmod version

    $regexp{plink2}{version} =
q?perl -nae 'if($_=~/PLINK\s(\S+\s\S+\s\S+\s\S+\s\S+)/) {my $ret = $1;$ret =~s/\s/_/g;print $ret;last;}' ?
      ;    #Collect Plink2 version

    $regexp{mendel}{fraction_of_errors} =
      q?perl -nae 'unless ($_=~/^#/) {print $F[1];last;}' ?;

    $regexp{mendel}{mendelian_errors} =
      q?perl -nae 'unless ($_=~/^#/) {print $F[2];last;}' ?;

    $regexp{father}{fraction_of_common_variants} =
      q?perl -nae 'unless ($_=~/^#/) {print $F[1];last;}' ?;

    $regexp{father}{common_variants} =
      q?perl -nae 'unless ($_=~/^#/) {print $F[2];last;}' ?;

    $regexp{tiddit}{version} =
q?perl -nae 'if($_=~/^##source=TIDDIT-(\S+)/) { print $1; last; } else { if($_=~/#CHROM/) { last;} }' ?;

    $regexp{svdb}{version} =
q?perl -nae 'if($_=~/^##SVDB_version=(\S+)/) { print $1; last; } else { if($_=~/#CHROM/) { last;} }' ?;

    $regexp{vcf2cytosure_version}{version} =
      q?perl -nae 'if($_=~/cytosure\s+(\d+[.]\d+[.]\d+)/xsm) { print $1;last; }' ?;

    ## Writes a YAML hash to file
    write_yaml(
        {
            yaml_href      => \%regexp,
            yaml_file_path => $print_regexp_outfile,
        }
    );
    return;
}
