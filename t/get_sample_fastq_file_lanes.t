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

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::File_info} => [qw{ get_sample_fastq_file_lanes }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ get_sample_fastq_file_lanes };

diag(   q{Test get_sample_fastq_file_lanes from File_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample_id with multiple lanes
my $sample_id = q{a_sample_id};
my %file_info = (
    $sample_id => {
        lanes => [ 1, 2 ],
    },
);
my @sample_fastq_lanes = get_sample_fastq_file_lanes(
    {
        file_info_href => \%file_info,
        sample_id      => $sample_id,
    }
);

## Then return lanes
my @expected = ( 1, 2 );
is_deeply( \@sample_fastq_lanes, \@expected, q{Returned lanes for sample id } );

done_testing();
