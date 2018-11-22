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
use Modern::Perl qw{ 2014 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Check::Path}    => [qw{ check_vcfanno_toml }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Path qw{ check_vcfanno_toml };

diag(   q{Test check_vcfanno_toml from Path.pm v}
      . $MIP::Check::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log();

## Given a toml config file
my $fqf_vcfanno_config =
  catfile( $Bin,
    qw{ data references GRCh37_frequency_vcfanno_filter_config_-v1.0-.toml  } );

my $is_ok = check_vcfanno_toml(
    {
        log               => $log,
        parameter_name    => q{fqf_vcfanno_config},
        vcfanno_file_toml => $fqf_vcfanno_config,
    }
);

## Then return true
ok( $is_ok, q{Passed check for toml file} );

## Given a tomlconfig file, when mandatory features are absent
my $faulty_fqf_vcfanno_config_file = catfile( $Bin,
    qw{ data references GRCh37_frequency_vcfanno_filter_config_bad_data_-v1.0-.toml } );

trap {
    check_vcfanno_toml(
        {
            log               => $log,
            parameter_name    => q{fqf_vcfanno_config},
            vcfanno_file_toml => $faulty_fqf_vcfanno_config_file,
        }
      )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the record does not match} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message for non matching reference} );

done_testing();
