#!/usr/bin/env perl

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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Pedigree}       => [qw{ check_pedigree_sample_allowed_values }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ check_pedigree_sample_allowed_values };

diag(   q{Test check_pedigree_sample_allowed_values from Pedigree.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

my %pedigree = (
    samples => [
        {
            analysis_type => q{wes},
            phenotype     => q{affected},
            sex           => q{female},
        },
        {
            analysis_type => q{wgs},
            phenotype     => q{unaffected},
            sex           => q{male},
        },
        {
            analysis_type => q{wts},
            phenotype     => q{unknown},
            sex           => q{other},
        },
        {
            analysis_type => q{wgs},
            phenotype     => q{unknown},
            sex           => q{unknown},
        },
    ],
);

## Given allowed values in pedigree file
my $is_ok = check_pedigree_sample_allowed_values(
    {
        pedigree_file_path => $Bin,
        pedigree_href      => \%pedigree,
    }
);

## Then return true
ok( $is_ok, q{Only allowed values in pedigree yaml file} );

## Given a not allowed value in pedigree file
$pedigree{samples}[0]{analysis_type} = q{How rude};

trap {
    check_pedigree_sample_allowed_values(
        {
            pedigree_file_path => $Bin,
            pedigree_href      => \%pedigree,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if not sample allowed value is found} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if not sample allowed value is found} );

done_testing();
