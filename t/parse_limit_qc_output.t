#! /usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Main::Qccollect} => [qw{ parse_limit_qc_output }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Main::Qccollect qw{ parse_limit_qc_output };

diag(   q{Test parse_limit_qc_output from Qccollect.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a qc data hash
my %qc_data = (
    sample => {
        sample_id => {
            recipe_output_file => {
                variantevalall => {
                    header => q{data},
                },
            },
            another_file => {
                collecthsmetrics => {
                    header => q{data},
                },
            },
        },
    },
);

## When executing sub
parse_limit_qc_output(
    {
        limit_qc_output => 1,
        qc_href         => \%qc_data,
    }
);

## Then remove key variantevalall
my %expected_qc_data = (
    sample => {
        sample_id => {
            recipe_output_file => {},
            another_file       => {
                collecthsmetrics => {
                    header => q{data},
                },
            },
        },
    },
);

is_deeply( \%qc_data, \%expected_qc_data, q{Delete qc data from hash} );

done_testing();
