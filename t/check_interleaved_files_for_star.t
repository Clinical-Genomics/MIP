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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Star}           => [qw{ check_interleaved_files_for_star }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Star qw{ check_interleaved_files_for_star };

diag(   q{Test check_interleaved_files_for_star from Star.pm v}
      . $MIP::Star::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given
my %file_info = test_mip_hashes( { mip_hash_name => q{file_info}, } );

my %active_parameter = ( sample_ids => [qw{ADM1059A1 ADM1059A2 ADM1059A3 }], );

trap {
    check_interleaved_files_for_star(
        {
            file_info_href => \%file_info,
            sample_ids_ref => \@{ $active_parameter{sample_ids} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if interleaved fastq files are found} );
like(
    $trap->stderr,
    qr/MIP \s+ rd_rna \s+  does \s+  not \s+  support \s+ interleaved/xms,
    q{Throw fatal log message if interleaved fastq files are found}
);

done_testing();
