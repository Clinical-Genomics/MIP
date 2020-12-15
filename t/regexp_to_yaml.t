#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Qcc_regexp}     => [qw{ regexp_to_yaml }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qcc_regexp qw{ regexp_to_yaml };

diag(   q{Test regexp_to_yaml from Qcc_regexp.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

my $test_dir       = File::Temp->newdir();
my $test_file_path = catfile( $test_dir, q{test.yaml} );

## Given a file path when not to print
my $is_ok = regexp_to_yaml(
    {
        log                  => $log,
        print_regexp_outfile => undef,
    }
);

## Then return true
ok( $is_ok, q{Return if not writing to file} );

## Given a file path when to print
trap {
    regexp_to_yaml(
        {
            log                  => $log,
            print_regexp_outfile => $test_file_path,
        }
    )
};

## Then write file, broadcast log message and exit
ok( -e $test_file_path, q{Wrote yaml file} );
is( $trap->leaveby, q{exit}, q{Exit if after writing} );
like( $trap->stderr, qr/Wrote\s+regexp\s+YAML/xms, q{Broadcast log message} );

done_testing();
