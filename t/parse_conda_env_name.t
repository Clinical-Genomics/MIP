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
use MIP::Constants qw{ $COMMA $NEWLINE $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli test_log };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Parse::Parameter} => [qw{ parse_conda_env_name }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Parameter qw{ parse_conda_env_name };

diag(   q{Test parse_conda_env_name from Parameter.pm v}
      . $MIP::Parse::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

Readonly my $DATE => 190624;

## Some test variables
my %parameter = (
    environment_name => {
        emip      => q{Test},
        etest_env => undef,
    },
);
my $base_name = q{base};

## Given an already defined name wihout any prefix, suffix or env date
my $environment_name = parse_conda_env_name(
    {
        base_name      => $base_name,
        date           => $DATE,
        environment    => q{emip},
        parameter_href => \%parameter,
    }
);

## Then return an unmodified enviornment name
is( $environment_name, q{Test}, q{Return environment name} );

## Given a request to add date, prefix and suffix to an undefined environment name
$parameter{add_environment_date} = 1;
$parameter{environment_prefix}   = q{D};
$parameter{environment_suffix}   = q{JD};

$environment_name = parse_conda_env_name(
    {
        base_name      => $base_name,
        date           => $DATE,
        environment    => q{etest_env},
        parameter_href => \%parameter,
    }
);

## Then return a environment name with all the bells and whistles
is( $environment_name, q{D_base_test_env_190624_JD},
    q{Set undefined environment name with prefix, date and suffix} );

done_testing();
