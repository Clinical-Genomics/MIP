#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.06;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter}    => [qw{ parse_program_executables }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ parse_program_executables };

diag(   q{Test parse_program_executables from Active_parameter.pm v}
      . $MIP::Active_parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );
$active_parameter{conda_path} = catfile( $Bin, qw{ data modules miniconda } );
$active_parameter{bwa_mem}    = 0;

my %parameter = test_mip_hashes( { mip_hash_name => q{define_parameter}, } );

## Given switched off active recipe parameter
my $return = parse_program_executables(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
    }
);

## Then return undef
is( $return, undef, q{Skip check for inactive recipes} );

## Given switched on active parameter, defined program and program_executables parameter, when program is in path and executable
$active_parameter{bwa_mem} = 1;

trap {
    parse_program_executables(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    )
};

## Then INFO message should broadcast
like( $trap->stderr, qr/samtools/xms,
    q{Found bin and executable: Throw INFO log message} );

## Given switched on active parameter, defined recipe and program_executables parameter, but no executable
$parameter{bwa_mem}{program_executables} = [qw{ not_in_path }];

trap {
    parse_program_executables(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    );
};

## Then exit and throw FATAL message
is( $trap->leaveby, q{die}, q{Exit if binary cannot be found} );
like( $trap->stderr, qr/not_in_path/xms,
    q{No bin and executable - Throw FATAL log message} );

done_testing();
