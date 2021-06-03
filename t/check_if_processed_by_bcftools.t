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
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Check::Reference} => [qw{ check_if_processed_by_bcftools }],
        q{MIP::Test::Fixtures}   => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Reference qw{ check_if_processed_by_bcftools };

diag(   q{Test check_if_processed_by_bcftools from Reference.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Given no reference path
#my $bcftools_binary_path = catfile( $Bin, qw{ data modules miniconda envs mip_ci bin bcftools } );
my $return = check_if_processed_by_bcftools(
    {
        reference_file_path => q{file_does_not_exists},
    }
);
is( $return, undef, q{No reference file to check} );

## Given a reference file path, when needing bcftools processing
my $reference_file_path_no_processing =
  catfile( $Bin, qw{ data references grch37_all_wgs_-phase3_v5b.2013-05-02-.vcf.gz } );

## When checking if bcftools has processed references using regexp
my @checked_references = check_if_processed_by_bcftools(
    {
        reference_file_path => $reference_file_path_no_processing,
    }
);

## Then return references to be processed
is( 1, scalar @checked_references, q{Detected Bcftools norm processing is needed} );

## Given a reference file path, when not needing bcftools processing
my $reference_file_path_bcftools =
  catfile( $Bin, qw{ data references grch37_gnomad.genomes_-r2.0.1-.vcf.gz } );

## When checking if bcftools has processed references using regexp
@checked_references = check_if_processed_by_bcftools(
    {
        reference_file_path => $reference_file_path_bcftools,
    }
);

## Then no references should be returned
is( 0, scalar @checked_references, q{Detected Bcftools norm processing} );

done_testing();
