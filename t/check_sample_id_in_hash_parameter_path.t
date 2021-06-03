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
        q{MIP::Active_parameter} => [qw{ check_sample_id_in_hash_parameter_path }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ check_sample_id_in_hash_parameter_path };

diag(   q{Test check_sample_id_in_hash_parameter_path from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given parameter and sample_ids, when sample id is duplicated in parameter
my %active_parameter = (
    infile_dirs => {
        infile_dir   => qw{sample-1},
        infile_dir_1 => qw{sample-1},
        infile_dir_2 => qw{sample-2},
    },
    sample_ids => [qw{ sample-1 sample-2 }],
);

trap {
    check_sample_id_in_hash_parameter_path(
        {
            active_parameter_href => \%active_parameter,
            parameter_names_ref   => [qw{ infile_dirs }],
            sample_ids_ref        => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the sample id is duplicated in parameter} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message for duplicated sample id} );

## Given parameter and sample_ids, when sample id is missing from parameter
%active_parameter = (
    infile_dirs => {
        infile_dir   => qw{not_included_sample},
        infile_dir_1 => qw{sample-1},
    },
    sample_ids => [qw{ sample-1 sample-2 }],
);

trap {
    check_sample_id_in_hash_parameter_path(
        {
            active_parameter_href => \%active_parameter,
            parameter_names_ref   => [qw{ infile_dirs }],
            sample_ids_ref        => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the sample id for parameter is missing} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message for missing sample id} );

## Given parameter and sample_ids, when all is ok
%active_parameter = (
    infile_dirs => {
        infile_dir_1 => qw{sample-1},
        infile_dir_2 => qw{sample-2},
    },
    sample_ids => [qw{ sample-1 sample-2 }],
);

my $is_ok = check_sample_id_in_hash_parameter_path(
    {
        active_parameter_href => \%active_parameter,
        parameter_names_ref   => [qw{ infile_dirs }],
        sample_ids_ref        => \@{ $active_parameter{sample_ids} },
    }
);
ok( $is_ok, q{Passed check} );

done_testing();
