#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
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
        q{MIP::Reference}      => [qw{ check_human_genome_file_endings }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ check_human_genome_file_endings };
use MIP::File_info qw{ set_human_genome_reference_features };

diag(   q{Test check_human_genome_file_endings from Reference.pm v}
      . $MIP::Reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( { no_screen => 1, } );

my %parameter;

my %active_parameter = ( human_genome_reference =>
      catfile( $Bin, qw{ data references grch37_homo_sapiens_-d5-.fasta } ), );

## File info hash
my %file_info = (

    # Human genome meta files
    human_genome_reference_file_endings => [qw{ .dict .fai }],
);

## Detect version and source of the human_genome_reference: Source (hg19 or GRCh).
set_human_genome_reference_features(
    {
        file_info_href         => \%file_info,
        human_genome_reference => basename( $active_parameter{human_genome_reference} ),
        parameter_href         => \%parameter,
    }
);

check_human_genome_file_endings(
    {
        human_genome_reference_file_endings_ref =>
          $file_info{human_genome_reference_file_endings},
        human_genome_reference_path => $active_parameter{human_genome_reference},
        parameter_href              => \%parameter,
        parameter_name              => q{human_genome_reference_file_endings},
    }
);

is( $parameter{human_genome_reference_file_endings}{build_file},
    0, q{Set build file switch for human genome reference to 0} );

$active_parameter{human_genome_reference} = q{not an existing reference};

check_human_genome_file_endings(
    {
        human_genome_reference_file_endings_ref =>
          $file_info{human_genome_reference_file_endings},
        human_genome_reference_path => $active_parameter{human_genome_reference},
        parameter_href              => \%parameter,
        parameter_name              => q{human_genome_reference_file_endings},
    }
);

is( $parameter{human_genome_reference_file_endings}{build_file},
    1, q{Set build file switch for human genome reference file endings to 1} );

done_testing();
