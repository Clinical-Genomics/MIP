#! /usr/bin/env perl

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
        q{MIP::Reference}      => [qw{ get_transcript_annotation_file_path }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ get_transcript_annotation_file_path };

diag(   q{Test get_transcript_annotation_file_path from Reference.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

my %active_parameter = ( transcript_annotation => catfile(qw{ path to annotation.gtf }), );
my %file_info =
  ( transcript_annotation_file_endings => [qw{ .bed .rrna.interval_list .refflat }], );

## Given a bed file annotation request
my $transcript_annotation_file_path = get_transcript_annotation_file_path(
    {
        active_parameter_href => \%active_parameter,
        file_format           => q{bed},
        file_info_href        => \%file_info,
    }
);

## Then return bed file
my $expected = catfile(qw{ path to annotation.gtf.bed });
is( $transcript_annotation_file_path, $expected, q{Return transcript annotation} );

## Given a gtf file annotation request, but the endings list isn't matching
trap {
    $transcript_annotation_file_path = get_transcript_annotation_file_path(
        {
            active_parameter_href => \%active_parameter,
            file_format           => q{refflat},
            file_info_href        => \%file_info,
        }
    )
};

## Then croak with error message
is( $trap->leaveby, q{die}, q{Exit if annotation format is wrong} );
like( $trap->die, qr/The\sretrieved\sannotation/xms, q{Throw error if annotation format is wrong} );

done_testing();
