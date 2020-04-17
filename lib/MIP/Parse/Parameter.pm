package MIP::Parse::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.15;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      parse_conda_env_name
      parse_download_reference_parameter
      parse_infiles
    };

}

sub parse_conda_env_name {

## Function : Build conda environment names depending on input parameters
## Returns  : $conda_environment_name
## Arguments: $base_name     => Degfault base environment name
##          : $date          => Date
##          : $environment   => Installation environment
##          : parameter_href => Parmeter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $base_name;
    my $date;
    my $environment;
    my $parameter_href;

    my $tmpl = {
        base_name => {
            defined     => 1,
            required    => 1,
            store       => \$base_name,
            strict_type => 1,
        },
        date => {
            defined     => 1,
            required    => 1,
            store       => \$date,
            strict_type => 1,
        },
        environment => {
            defined     => 1,
            required    => 1,
            store       => \$environment,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $environment_name = $parameter_href->{environment_name}{$environment};

    ## Give the env a default name if not given
    if ( not $environment_name ) {

        ## Strip the first character, i.e. 'e' from the environment string
        my $env_postfix = substr $environment, 1;
        $environment_name = $base_name . $UNDERSCORE . $env_postfix;
    }

    ## Prepend environemnt prefix
    if ( $parameter_href->{environment_prefix} ) {

        $environment_name =
          $parameter_href->{environment_prefix} . $UNDERSCORE . $environment_name;
    }

    ## Add environment date
    if ( $parameter_href->{add_environment_date} ) {

        $environment_name = $environment_name . $UNDERSCORE . $date;
    }

    ## Append environment suffix
    if ( $parameter_href->{environment_suffix} ) {

        $environment_name =
          $environment_name . $UNDERSCORE . $parameter_href->{environment_suffix};
    }

    return $environment_name;
}

sub parse_download_reference_parameter {

## Function : Remodel depending on if "--reference" was used or not as the user info is stored as a scalar per reference_id while yaml is stored as arrays per reference_id
## Returns  :
## Arguments: $reference_href => Reference hash {REF}

    my ($arg_href) = @_;

## Flatten argument(s)
    my $reference_href;

    my $tmpl = {
        reference_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  VERSION_REF:
    foreach my $versions_ref ( values %{$reference_href} ) {

        if ( ref $versions_ref ne q{ARRAY} ) {

            ## Make scalar from CLI '--ref key=value' option into array
            $versions_ref = [$versions_ref];
        }
    }

    return;
}

sub parse_infiles {

## Function : Collects the ".fastq(.gz)" files from the supplied infiles directory. Checks if any files exist.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_info_href        => File info hash {REF}
##          : $log                   => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $log;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_matching_values_key };
    use MIP::Check::Parameter qw{ check_infiles };
    use MIP::File::Path qw{ get_file_names };
    use MIP::Set::File qw{ set_infiles };

    ## Collect input files governed by sample_ids
  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Return the key if the hash value exists
        my $infile_directory = get_matching_values_key(
            {
                active_parameter_href => $active_parameter_href,
                parameter_name        => q{infile_dirs},
                query_value           => $sample_id,
            }
        );

        _check_infile_directory(
            {
                infile_directory => $infile_directory,
                sample_id        => $sample_id,
            }
        );

        ## Get the file(s) from filesystem
        my @infiles = get_file_names(
            {
                file_directory   => $infile_directory,
                rule_name        => q{*.fastq*},
                rule_skip_subdir => q{original_fastq_files},
            }
        );

        ## Check infiles found and that they contain sample_id
        check_infiles(
            {
                infiles_ref      => \@infiles,
                infile_directory => $infile_directory,
                log              => $log,
                sample_id        => $sample_id,
            }
        );

        ## Set the infile features i.e. dir and files
        set_infiles(
            {
                file_info_href   => $file_info_href,
                infiles_ref      => \@infiles,
                infile_directory => $infile_directory,
                sample_id        => $sample_id,
            }
        );

        ## Broadcast to user
        $log->info(q{Reads from platform:});
        $log->info( q{Sample id: } . $sample_id );
        $log->info(qq{\tInputfiles:});

        ## Log each file from platform
      FILE:
        foreach my $file (@infiles) {

            # Indent for visability
            $log->info( qq{\t\t}, $file );
        }
    }
    return 1;
}

sub _check_infile_directory {

## Function : Check if infile directory exists per sample id
## Returns  :
## Arguments: $infile_directory => Infile directory
##          : $sample_id        => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_directory;
    my $sample_id;

    my $tmpl = {
        infile_directory => {
            required    => 1,
            store       => \$infile_directory,
            strict_type => 1,
        },
        sample_id => {
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Log::MIP_log4perl qw{ retrieve_log };

    my $log = retrieve_log( { log_name => $LOG_NAME, } );

    return if ( defined $infile_directory );

    $log->fatal(
        q{Could not detect any supplied '--infile_dirs' for sample: } . $sample_id );
    exit 1;
}

1;
