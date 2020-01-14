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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Pedigree}       => [qw{ check_pedigree_vs_user_input_sample_ids }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ check_pedigree_vs_user_input_sample_ids };

diag(   q{Test check_pedigree_vs_user_input_sample_ids from Pedigree.pm v}
      . $MIP::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given matching sample ids
my @pedigree_sample_ids   = qw{ sample_1 sample_2 sample_3 sample_4 };
my @user_input_sample_ids = qw{ sample_1 sample_2 sample_3 sample_4 };

my $is_ok = check_pedigree_vs_user_input_sample_ids(
    {
        pedigree_file_path        => $Bin,
        pedigree_sample_ids_ref   => \@pedigree_sample_ids,
        user_input_sample_ids_ref => \@user_input_sample_ids,
    }
);

## Then return true
ok( $is_ok, q{Sample ids matched} );

## Given a sample id that is not present in the pedigree
push @user_input_sample_ids, q{Han Solo};

trap {
    check_pedigree_vs_user_input_sample_ids(
        {
            pedigree_file_path        => $Bin,
            pedigree_sample_ids_ref   => \@pedigree_sample_ids,
            user_input_sample_ids_ref => \@user_input_sample_ids,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the path cannot be found} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message} );

done_testing();
