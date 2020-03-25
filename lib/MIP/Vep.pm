package MIP::Vep;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { any };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_vep_api_cache_versions };
}

sub check_vep_api_cache_versions {

## Function : Compare VEP API and VEP chache versions. Exit if non-match
## Returns  :
## Arguments: $vep_binary_path     => VEP binary path {REF}
##          : $vep_directory_cache => VEP cache directory path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_directory_cache;

    ## Default(s)
    my $vep_binary_path;

    my $tmpl = {
        vep_binary_path => {
            default     => q{vep},
            defined     => 1,
            store       => \$vep_binary_path,
            strict_type => 1,
        },
        vep_directory_cache => {
            defined     => 1,
            required    => 1,
            store       => \$vep_directory_cache,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use File::Find::Rule;
    use MIP::Environment::Executable qw{ get_binary_version };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my $vep_version = get_binary_version(
        {
            binary      => q{vep},
            binary_path => $vep_binary_path,
        }
    );

    ## Check that a version number was picked up
    if ( not $vep_version ) {
        $log->warn(
q{Could not retrieve VEP version. Skipping checking that VEP api and cache matches.}
        );
        return;
    }

    ## Get absolute path to VEP species cache as it is commonly linked
    my $vep_cache_species_dir_path =
      _get_vep_cache_species_dir_path( { vep_directory_cache => $vep_directory_cache, } );

    ## Get folders in cache directory
    # Build rule
    my $rule = File::Find::Rule->new();

    # Find directories
    $rule->directory;

    # Only get directories in the current folder
    $rule->maxdepth(1);

    # Get relative paths
    $rule->relative;

    # Apply rule to find directories
    my @vep_cache_version_folders = $rule->in($vep_cache_species_dir_path);

    ## Check that
    if ( not @vep_cache_version_folders ) {
        $log->warn(
q{Could not retrieve VEP cache version. Skipping checking that VEP api and cache matches.}
        );
        return;
    }

    ## Check if the VEP api version and cache versions matches
    if ( not any { /$vep_version/xms } @vep_cache_version_folders ) {

        $log->fatal( q{Differing versions between 'VEP API':}
              . $SPACE
              . $vep_version
              . $SPACE
              . q{and '--vep_directory_cache':}
              . $SPACE
              . $vep_directory_cache );
        exit 1;
    }
    return 1;
}

sub _get_vep_cache_species_dir_path {

## Function : Get the vep cache species dir path
## Returns  : $vep_cache_species_dir_path
## Arguments: $vep_directory_cache => VEP cache directory path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_directory_cache;

    my $tmpl = {
        vep_directory_cache => {
            defined     => 1,
            required    => 1,
            store       => \$vep_directory_cache,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Cwd qw{ abs_path };

    my @vep_species_cache = qw{ homo_sapiens homo_sapiens_merged };
    my $vep_cache_species_dir_path;

  SPECIES_CACHE:
    foreach my $species_cache (@vep_species_cache) {

        ## Get species specific cache dir
        $vep_cache_species_dir_path =
          abs_path( catdir( $vep_directory_cache, $species_cache ) );
        last if ( -e $vep_cache_species_dir_path );
    }
    return $vep_cache_species_dir_path;
}

1;
