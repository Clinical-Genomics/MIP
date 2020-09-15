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
        q{MIP::Reference}      => [qw{ check_nist_file_name }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ check_nist_file_name };

diag(   q{Test check_nist_file_name from Reference.pm v}
      . $MIP::Reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given nist info
my $file_name      = q{grch37_nist_hg001_-na12878_v3.3.2-.vcf};
my $nist_id        = q{NA12878};
my $nist_parameter = q{nist_call_set_vcf};
my $nist_version   = q{3.3.2};

my $is_ok = check_nist_file_name(
    {
        file_name      => $file_name,
        nist_id        => $nist_id,
        nist_parameter => $nist_parameter,
        nist_version   => $nist_version,
    }
);

## Then return true
ok( $is_ok, q{Checked nist file name is defined in nist hash parameters} );

## Given an undefined file name

trap {
    check_nist_file_name(
        {
            file_name      => undef,
            nist_id        => $nist_id,
            nist_parameter => $nist_parameter,
            nist_version   => $nist_version,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if file name is undef} );
like(
    $trap->stderr,
    qr/Please\s+supply\s+a\s+file\s+name/xms,
    q{Throw fatal log message}
);

done_testing();
