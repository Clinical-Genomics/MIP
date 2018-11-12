package MIP::Recipes::Analysis::Analysisrunstatus;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use FindBin qw{ $Bin };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_analysisrunstatus };

}

## Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $TAB        => qq{\t};
Readonly my $UNDERSCORE => q{_};

sub analysis_analysisrunstatus {

## Function : Execute last in MAIN chain, tests that all recorded files exists, have a file sixe greater than zero, checks QC-metrics for PASS or FAIL and sets analysis run status flag to finished.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_path_entries };
    use MIP::Get::Parameter qw{ get_recipe_parameters };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set recipe mode
    my $recipe_mode = $active_parameter_href->{$recipe_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$recipe_name}{chain};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ($recipe_file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $case_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            log                   => $log,
            process_time          => $time,
            recipe_directory      => $recipe_name,
            recipe_name           => $recipe_name,
        }
    );

    ## Set status flagg so that perl not_finished remains in sample_info_file
    say {$FILEHANDLE} q?STATUS="0"?;

    my @paths;

    ## Collects all recipes file path(s) created by MIP located in %sample_info
    get_path_entries(
        {
            paths_ref        => \@paths,
            sample_info_href => $sample_info_href,
        }
    );

    ### Test all file that are supposed to exists as they are present in the sample_info file
    _check_mip_analysis_files(
        {
            FILEHANDLE => $FILEHANDLE,
            paths_ref  => \@paths
        }
    );

    ## Test varianteffectpredictor fork status. If varianteffectpredictor is unable to fork it will prematurely end the analysis and we will lose variants.
    my $variant_effect_predictor_file =
      $sample_info_href->{recipe}{varianteffectpredictor}{stderrfile}{path};

    ## Test peddy warnings
    my $peddy_file =
      $sample_info_href->{recipe}{peddy_ar}{stderr}{path};

    ## Test if FAIL exists in qccollect file i.e. issues with samples e.g. Sex and seq data correlation, relationship etc
    my $qccollect_file;
    if ( not $active_parameter_href->{qccollect_skip_evaluation} ) {

        $qccollect_file = $sample_info_href->{recipe}{qccollect_ar}{path};
    }

    my %files_to_check = (
        q{FAIL}                   => $qccollect_file,
        q{pedigree warning:}      => $peddy_file,
        q{WARNING Unable to fork} => $variant_effect_predictor_file,
    );

  CHECK_FILE:
    while ( my ( $file_string_to_match, $file ) = each %files_to_check ) {

        _check_string_within_file(
            {
                file            => $file,
                FILEHANDLE      => $FILEHANDLE,
                string_to_match => $file_string_to_match,
            }
        );
    }

    ## Test integrity of vcf data keys in header and body
    my %vcf_file = (
        sv_vcf_file => [qw{ clinical research }],
        vcf_file    => [qw{ clinical research }],
    );

    _check_vcf_header_and_keys(
        {
            analysis_config_file => $active_parameter_href->{config_file_analysis},
            FILEHANDLE           => $FILEHANDLE,
            sample_info_href     => $sample_info_href,
            vcf_file_href        => \%vcf_file,
        }
    );

    ## Eval status flag
    _eval_status_flag(
        {
            FILEHANDLE       => $FILEHANDLE,
            sample_info_file => $active_parameter_href->{sample_info_file},
        }
    );

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    if ( $recipe_mode == 1 ) {

        submit_recipe(
            {
                dependency_method   => q{add_to_all},
                job_dependency_type => q{afterok},
                job_id_href         => $job_id_href,
                log                 => $log,
                job_id_chain        => $job_id_chain,
                recipe_file_path    => $recipe_file_path,
                submission_profile  => $active_parameter_href->{submission_profile},
            }
        );
    }
    return;
}

sub _eval_status_flag {

## Function : Eval status flag
## Returns  :
## Arguments: $FILEHANDLE       => Filehandle to write to
##          : $sample_info_file => Sample info file for the analysis

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $sample_info_file;

    my $tmpl = {
        FILEHANDLE => {
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
        },
        sample_info_file => {
            defined  => 1,
            required => 1,
            store    => \$sample_info_file,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Eval status value
    say {$FILEHANDLE} q?if [ $STATUS -ne 1 ]; then?;

    ## Execute perl
    print {$FILEHANDLE} $TAB . q?perl -i -p -e '?;

    ## Find analysisrunstatus line
    print {$FILEHANDLE} q?if($_=~/analysisrunstatus\:/) { ?;

    ## All ok - set runstatus mode to finished
    say {$FILEHANDLE} q?s/not_finished/finished/g }' ? . $sample_info_file . q? ?;

    ## Found discrepancies - exit
    say {$FILEHANDLE} q?else?;
    say {$FILEHANDLE} $TAB . q?exit 1?;
    say {$FILEHANDLE} q?fi?, $NEWLINE;

    return;
}

sub _check_mip_analysis_files {

## Function : Test all file that are supposed to exists as they are present in the sample_info file
## Returns  :
## Arguments: $FILEHANDLE => Filehandle to write to
##          : $paths_ref  => Paths to files to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $paths_ref;

    my $tmpl = {
        FILEHANDLE => {
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
        },
        paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$paths_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Create bash array
    print {$FILEHANDLE} q?readonly FILES=(?;

  PATH:
    foreach my $path ( @{$paths_ref} ) {

        ## First analysis and dry run will otherwise cause try to print uninitialized values
        next PATH if ( not defined $path );

        ## Add to array
        print {$FILEHANDLE} q?"? . $path . q?" ?;
    }

    ## Close bash array
    say {$FILEHANDLE} q?)?;

    ## Loop over files
    say {$FILEHANDLE} q?for file in "${FILES[@]}"?;

    ## For each element in array do
    say {$FILEHANDLE} q?do? . $SPACE;

    ## File exists and is larger than zero
    say {$FILEHANDLE} $TAB . q?if [ -s "$file" ]; then?;

    ## Echo
    say {$FILEHANDLE} $TAB x 2 . q?echo "Found file $file"?;
    say {$FILEHANDLE} $TAB . q?else?;

    ## Redirect to STDERR
    say {$FILEHANDLE} $TAB x 2 . q?echo "Could not find $file" >&2?;

    ## Set status flagg so that perl notFinished remains in sample_info_file
    say {$FILEHANDLE} $TAB x 2 . q?STATUS="1"?;
    say {$FILEHANDLE} $TAB . q?fi?;
    say {$FILEHANDLE} q?done ?, $NEWLINE;

    return;
}

sub _check_string_within_file {

## Function : Test presence of string within file
## Returns  :
## Arguments: $file            => Files to check
##          : $FILEHANDLE      => Filehandle to write to
##          : $string_to_match => String to match within file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file;
    my $FILEHANDLE;
    my $string_to_match;

    my $tmpl = {
        file => {
            required    => 1,
            store       => \$file,
            strict_type => 1,
        },
        FILEHANDLE => {
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
        },
        string_to_match => {
            defined     => 1,
            required    => 1,
            store       => \$string_to_match,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( defined $file ) {

        ## Not output the matched text only return the exit status code
        print {$FILEHANDLE} q?if grep -q "? . $string_to_match . q?" ?;

        ## Infile
        say {$FILEHANDLE} $file . q?; then?;

        ## Found pattern
        say {$FILEHANDLE} $TAB . q?STATUS="1"?;

        ## Echo FAILED
        say {$FILEHANDLE} $TAB
          . q?echo "String match status=FAILED for file: ?
          . $file
          . q?" >&2?;

        ## Infile is clean
        say {$FILEHANDLE} q?else?;

        ## Echo PASSED
        say {$FILEHANDLE} $TAB
          . q?echo "String match status=PASSED for file: ?
          . $file
          . q?" >&2?;
        say {$FILEHANDLE} q?fi?, $NEWLINE;
    }
    return;
}

sub _check_vcf_header_and_keys {

## Function : Test integrity of vcf data keys in header and body
## Returns  :
## Arguments: $analysis_config_file => Config file for the analysis
##          : $FILEHANDLE           => Filehandle to write to
##          : $vcf_file_href        => Files to check
##          : $sample_info_href     => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_config_file;
    my $FILEHANDLE;
    my $vcf_file_href;
    my $sample_info_href;

    my $tmpl = {
        analysis_config_file => {
            defined  => 1,
            required => 1,
            store    => \$analysis_config_file,
        },
        FILEHANDLE => {
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
        },
        vcf_file_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vcf_file_href,
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

  FILE:
    foreach my $file ( keys %{$vcf_file_href} ) {

      MODE:
        foreach my $mode ( @{ $vcf_file_href->{$file} } ) {

            next MODE
              if ( not defined $sample_info_href->{$file}{$mode}{path} );

            ## Execute on cmd
            print {$FILEHANDLE} q?perl -MTest::Harness -e ' ?;

            ## Adjust arguments to harness object
            print {$FILEHANDLE} q?my %args = (?;

            ## Print individual test results to STDOUT
            print {$FILEHANDLE} q?verbosity => 1, ?;

            ##Argument to test script
            print {$FILEHANDLE} q?test_args => { ?;

            ## Add test for select file using alias
            print {$FILEHANDLE} q?"test ? . $mode . $SPACE . $file . q?" => [ ?;

            ## Infile
            print {$FILEHANDLE} q?"? . $sample_info_href->{$file}{$mode}{path} . q?", ?;

            ##ConfigFile
            print {$FILEHANDLE} q?"? . $analysis_config_file . q?", ?;
            print {$FILEHANDLE} q?], ?;

            print {$FILEHANDLE} q?}); ?;

            ## Create harness using arguments provided
            print {$FILEHANDLE} q?my $harness = TAP::Harness->new( \%args ); ?;

            ## Execute test(s)
            print {$FILEHANDLE} q?$harness->runtests( ?;

            print {$FILEHANDLE} q?["?
              . catfile( $Bin, qw{ t mip_analysis.test } )
              . q?", "test ?
              . $mode
              . $SPACE
              . $file . q?"], ?;

            print {$FILEHANDLE} q?)'?;
            say   {$FILEHANDLE} $NEWLINE;
        }
    }
    return;
}

1;
