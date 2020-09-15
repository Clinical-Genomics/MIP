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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Reference}      => [qw{ parse_exome_target_bed }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ parse_exome_target_bed };

diag(   q{Test parse_exome_target_bed from Reference.pm v}
      . $MIP::Reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log({no_screen => 1,});

## Test hashes
my %file_info = (
    human_genome_reference_version => q{37},
    human_genome_reference_source  => q{grch},
);

my %active_parameter = (
    exome_target_bed => {
        q{genome_reference_source_version_agilent_sureselect_targets_cre_-v1-.bed} =>
          q{sample_1},
        q{test_capture.bed} => q{sample_2},
    },
);

parse_exome_target_bed(
    {
        exome_target_bed_file_href     => $active_parameter{exome_target_bed},
        human_genome_reference_source  => $file_info{human_genome_reference_source},
        human_genome_reference_version => $file_info{human_genome_reference_version},
    }
);

CAPTURE_FILE:
foreach my $capture_file ( keys %{ $active_parameter{exome_target_bed} } ) {

    if ( $capture_file ne q{test_capture.bed} ) {

        like( $capture_file, qr/grch/xsm, q{Updated genome source} );

        like( $capture_file, qr/37/xsm, q{Updated genome version} );
    }
    else {

        unlike( $capture_file, qr/grch/xsm, q{Did not updated genome source} );

        unlike( $capture_file, qr/37/xsm, q{Did not updated genome version} );
    }
}

done_testing();
