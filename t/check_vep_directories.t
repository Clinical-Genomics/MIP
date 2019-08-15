#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Check::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Check::Parameter qw{ check_vep_directories };

diag(   q{Test check_vep_directories from Parameter.pm v}
      . $MIP::Check::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir      = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{test},
    }
);

## Given matching vep API and cache version
my $vep_directory_path =
  catdir( $Bin, qw{ data modules miniconda envs test_env ensembl-vep-91 } );
my $vep_directory_cache =
  catdir( $Bin, qw{ data modules miniconda envs test_env ensembl-tools-91 cache } );
## When comparing API and cache version
my $match = check_vep_directories(
    {
        log                 => $log,
        vep_directory_cache => $vep_directory_cache,
        vep_directory_path  => $vep_directory_path,
    }
);
## Then return true
ok( $match, q{Return on matching versions} );

## Given non matching API and cache
$vep_directory_cache =
  catdir( $Bin, qw{ data modules miniconda envs test_env ensembl-tools-92 cache } );
## When comparing API and cache version
trap {
    check_vep_directories(
        {
            log                 => $log,
            vep_directory_cache => $vep_directory_cache,
            vep_directory_path  => $vep_directory_path,
        }
    )
};
## Then exit and print fatal error message
ok( $trap->exit, q{Exit on non matching versions} );
like( $trap->stderr, qr/FATAL/xms, q{Throw FATAL log message} );

## Given a cache folder lacking the homo_sapiens subdirectory
$vep_directory_cache =
  catdir( $Bin, qw{ data modules miniconda envs test_env ensembl-vep } );
## When trying to retireve the cache versions
trap {
    check_vep_directories(
        {
            log                 => $log,
            vep_directory_cache => $vep_directory_cache,
            vep_directory_path  => $vep_directory_path,
        }
    )
};
## Then return and print warning message
ok( $trap->return, q{Return on unknown cache version} );
like( $trap->stderr, qr/WARN/xms, q{Warn for unknown VEP cache version} );

## Given a directory lacking the version number in the .version/ensembl file
$vep_directory_path =
  catdir( $Bin, qw{ data modules miniconda envs test_env_1 ensembl-vep } );
$vep_directory_cache =
  catdir( $Bin, qw{ data modules miniconda envs test_env ensembl-tools-91 cache } );
## When trying to retrieve API version
trap {
    check_vep_directories(
        {
            log                 => $log,
            vep_directory_cache => $vep_directory_cache,
            vep_directory_path  => $vep_directory_path,
        }
    )
};
## Then return and print warning message
ok( $trap->return, q{Return on unknown API version} );
like( $trap->stderr, qr/WARN/xms, q{Warn for unknown VEP api version} );

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

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
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help     Display this help message
    -v/--version  Display version
END_USAGE
}
