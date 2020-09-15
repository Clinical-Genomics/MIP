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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::File_info}      => [qw{ parse_select_file_contigs }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ parse_select_file_contigs };

diag(   q{Test parse_select_file_contigs from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 0, } );

## Given a vcfparser select file path and contigs
my $consensus_analysis_type = q{wes};
my %file_info               = (
    contigs => [qw{ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT }],
    contigs_size_ordered =>
      [qw{ 1 2 3 4 5 6 7 X 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 Y MT }]
);

my $vcfparser_select_file_path =
  catfile( $Bin, qw{ data 643594-miptest aggregated_gene_panel_test.txt } );

my $is_ok = parse_select_file_contigs(
    {
        consensus_analysis_type => $consensus_analysis_type,
        file_info_href          => \%file_info,
        select_file_path        => $vcfparser_select_file_path,
    }
);

## Then
ok( $is_ok, q{Parsed select file contigs} );

done_testing();
