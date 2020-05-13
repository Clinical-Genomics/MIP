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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::File_info}      => [qw{ add_sample_fastq_file_lanes }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ add_sample_fastq_file_lanes };

diag(   q{Test add_sample_fastq_file_lanes from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample_id and lane
my %file_info;
my $lane      = 1;
my $sample_id = q{sample_id};

## When direction is not one
my $direction = 2;

add_sample_fastq_file_lanes(
    {
        direction      => $direction,
        file_info_href => \%file_info,
        lane           => $lane,
        sample_id      => $sample_id,
    }
);

## Then do not add lanes
is( @{ $file_info{$sample_id}{lanes} },
    0, q{Did not add lane for direction other than 1} );

## When direction is one
$direction = 1;

add_sample_fastq_file_lanes(
    {
        direction      => $direction,
        file_info_href => \%file_info,
        lane           => $lane,
        sample_id      => $sample_id,
    }
);

## Then add lanes
is( @{ $file_info{$sample_id}{lanes} }, 1, q{Added lane for direction equals 1} );

done_testing();
