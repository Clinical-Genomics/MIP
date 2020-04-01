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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE $TAB };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_log test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Gatk}           => [qw{ check_gatk_sample_map_paths }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Gatk qw{ check_gatk_sample_map_paths };

diag(   q{Test check_gatk_sample_map_paths from Gatk.pm v}
      . $MIP::Gatk::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log      = test_log( {} );
my $test_dir = File::Temp->newdir();

## Given a sample file with paths that exists
my $gatk_genotypegvcfs_ref_gvcf = catfile( $test_dir, q{sample_map_file.txt} );

## Create anonymous filehandle
my $filehandle = IO::Handle->new();

open $filehandle, q{>}, $gatk_genotypegvcfs_ref_gvcf
  or $log->logdie( q{Cannot open '} . $gatk_genotypegvcfs_ref_gvcf . q{': } . $OS_ERROR );

## Write to file
say {$filehandle} q{100-1-2A-REF}
  . $TAB
  . catfile( $Bin,
    qw{ data references grch37_merged_reference_infiles_-2014-.g.100-1-2A-REF.vcf } );
close $filehandle;

my $is_ok = check_gatk_sample_map_paths(
    {
        sample_map_path => $gatk_genotypegvcfs_ref_gvcf,
    }
);
## Then return true
ok( $is_ok, q{Found all paths in sample_map_path} );

## Given a non existing file path with trailing garbage
open $filehandle, q{>>}, $gatk_genotypegvcfs_ref_gvcf
  or $log->logdie( q{Cannot open '} . $gatk_genotypegvcfs_ref_gvcf . q{': } . $OS_ERROR );

## Write to file
say {$filehandle} q{sample_1} . $TAB
  . catfile( $Bin, q{not_a_file_path} . $TAB . q{garbage} );
say {$filehandle} q{sample_1} . $TAB . catfile( $Bin, q{not_a_file_path_again} );

close $filehandle;

trap {
    check_gatk_sample_map_paths(
        {
            sample_map_path => $gatk_genotypegvcfs_ref_gvcf,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the path cannot be found} );
like( $trap->stderr, qr{FATAL}xms, q{Throw fatal log message} );
like( $trap->stderr, qr{Unexpected\s+trailing\s+garbage}xms, q{Found trailing garbage} );

done_testing();
