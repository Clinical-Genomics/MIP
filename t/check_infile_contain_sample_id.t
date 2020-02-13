#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = '1.0.1';

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
        q{MIP::Check::Parameter} => [qw{ check_infile_contain_sample_id }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Parameter qw{ check_infile_contain_sample_id };

diag(   q{Test check_infile_contain_sample_id from Parameter.pm v}
      . $MIP::Check::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given proper input data
my $infile_sample_id = q{sample-1};
my $file_index       = 0;
my $sample_id        = q{sample-1};

my %active_parameter = ( sample_ids  => [ qw{ sample-1 sample-2 }, ], );
my %infile           = ( q{sample-1} => { mip_infiles => [qw{ sample-1.fastq }], }, );

my $is_ok = check_infile_contain_sample_id(
    {
        infile_name      => $infile{q{sample-1}}{mip_infiles}[0],
        infile_sample_id => $infile_sample_id,
        log              => $log,
        sample_id        => $sample_id,
        sample_ids_ref   => \@{ $active_parameter{sample_ids} },
    }
);
## Then return true
ok( $is_ok, q{Sample ids match} );

## Given non matching sample id
$infile_sample_id = q{a.k.a gone missing};

trap {
    check_infile_contain_sample_id(
        {
            infile_name      => $infile{q{sample-1}}{mip_infiles}[0],
            infile_sample_id => $infile_sample_id,
            log              => $log,
            sample_id        => $sample_id,
            sample_ids_ref   => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if sample_id do not match} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message if sample_id do not match} );

done_testing();
