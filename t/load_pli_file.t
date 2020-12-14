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
        q{MIP::File::Format::Pli} => [qw{ load_pli_file }],
        q{MIP::Test::Fixtures}    => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Pli qw{ load_pli_file };

diag(   q{Test load_pli_file from Pli.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $ADK_PLI => 0.91;

my $log = test_log( {} );

## Given a pli file
my $pli_values_file_path =
  catfile( $Bin, qw{ data references gnomad_pli_per_gene_-_r2.1.1-.txt } );
my %pli_score;

my $is_ok = load_pli_file(
    {
        infile_path    => $pli_values_file_path,
        log            => $log,
        pli_score_href => \%pli_score,
    }
);

my %expected_pli_score = ( ADK => $ADK_PLI );

## Then load plI file should return true
ok( $is_ok, q{Loaded plI file} );

## Then pli_score hash should be loaded with hgnc_symbol and corresponding pli score
is_deeply( \%pli_score, \%expected_pli_score, q{Loaded plI data into hash} );

done_testing();
