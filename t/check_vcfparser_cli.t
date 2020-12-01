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

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Vcfparser}      => [qw{ check_vcfparser_cli }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ check_vcfparser_cli };

diag(   q{Test check_vcfparser_cli from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $SELECT_FEATURE_MATCHING_COLUMN => 3;

## Initiate global log name
test_log( {} );

my $range_feature_file               = catfile(qw{ a range file });
my @range_feature_annotation_columns = ( 1, 2 );
my $select_feature_file              = catfile(qw{ a test select_file.bed });
my $select_outfile                   = catfile(qw{ a test select_outfile.vcf });

## Given a select feature file and no select feature matching column
trap {
    check_vcfparser_cli(
        {
            range_feature_file              => $range_feature_file,
            range_feature_annotation_column => $range_feature_annotation_columns[0],
            select_feature_file             => $select_feature_file,
            select_feature_matching_column  => undef,
            select_outfile                  => $select_outfile,
        }
    );
};

## Then exit and throw FATAL log message if supplied select feature file and no matching column
ok( $trap->exit, q{Exit if supplied select feature file and no matching column} );
like( $trap->stderr, qr/feature\s+column/xms,
    q{Throw fatal log message if supplied select feature file and no matching column} );

## Given a select feature file and not select feature outfile
trap {
    check_vcfparser_cli(
        {
            range_feature_file              => $range_feature_file,
            range_feature_annotation_column => $range_feature_annotation_columns[0],
            select_feature_file             => $select_feature_file,
            select_feature_matching_column  => $SELECT_FEATURE_MATCHING_COLUMN,
            select_outfile                  => undef,
        }
    )
};

## Then exit and throw FATAL log message if supplied select feature file and no matching column
ok( $trap->exit, q{Exit if supplied select feature file and no select outfile} );
like( $trap->stderr, qr/select\s+outfile/xms,
    q{Throw fatal log message if supplied select feature file and no select outfile} );

## Given a range feature file and not range feature file columns
trap {
    check_vcfparser_cli(
        {
            range_feature_annotation_column => undef,
            range_feature_file              => $range_feature_file,
            select_feature_file             => $select_feature_file,
            select_feature_matching_column  => $SELECT_FEATURE_MATCHING_COLUMN,
            select_outfile                  => $select_outfile,
        }
    )
};

## Then exit and throw FATAL log message if supplied select feature file and no matching column
ok( $trap->exit, q{Exit if supplied range feature file and no range anno columns} );
like( $trap->stderr, qr/feature\s+column[(]s[)]/xms,
    q{Throw fatal log message if supplied select feature file and no range anno columns}
);

done_testing();
