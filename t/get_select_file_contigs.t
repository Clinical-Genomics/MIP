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
        q{MIP::Reference}      => [qw{ get_select_file_contigs }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ get_select_file_contigs };

diag(   q{Test get_select_file_contigs from File.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $AUTOSOMAL_CONTIG_NR => 22;

## Creates log object
my $log = test_log( {} );

## Given proper input data
my %file_info;
my $select_file_path =
  catfile( $Bin, qw{ data 643594-miptest aggregated_gene_panel_test.txt } );

@{ $file_info{select_file_contigs} } = get_select_file_contigs(
    {
        select_file_path => $select_file_path,
    }
);
my @expected_contigs = ( 1 .. $AUTOSOMAL_CONTIG_NR, qw{ X Y MT} );
@expected_contigs = sort @expected_contigs;

## Then return the expected contigs
is_deeply( \@{ $file_info{select_file_contigs} },
    \@expected_contigs, q{Got select file contigs} );

## Given bed file with no contigs
my $wrong_file =
  catfile( $Bin, qw{ data 643594-miptest aggregated_gene_panel_test_no_contigs.txt } );

trap {
    get_select_file_contigs(
        {
            select_file_path => $wrong_file,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if contigs cannot be found} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if contigs cannot be found} );

done_testing();
