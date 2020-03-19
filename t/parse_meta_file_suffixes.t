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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

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
        q{MIP::Reference}      => [qw{ parse_meta_file_suffixes }],
        q{MIP::Io::Read}       => [qw{ read_from_file }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ parse_meta_file_suffixes };
use MIP::Io::Read qw{ read_from_file };

diag(   q{Test parse_meta_file_suffixes from Reference.pm v}
      . $MIP::Reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %parameter = read_from_file(
    {
        format => q{yaml},
        path   => catfile( dirname($Bin), qw{ definitions rd_dna_parameters.yaml} ),
    }
);

my %active_parameter = (
    exome_target_bed => {
        catfile( $Bin,
            qw{ data references grch37_agilent_sureselect_targets_cre_-v1-.bed } ) =>
          q{sample1},
    },
);

## File info hash
my %file_info = (

    # BWA human genome reference file endings
    bwa_build_reference => [qw{ .bwt .ann .amb .pac .sa }],

    exome_target_bed => [qw{ .interval_list .pad100.interval_list }],

    # Human genome meta files
    human_genome_reference_file_endings => [qw{ .dict .fai }],

    # RTG human genome reference file endings
    rtg_vcfeval_reference_genome => [qw{ _sdf_dir }],
);

## Active parameter
my $parameter_name = q{exome_target_bed};
my $parameter      = $active_parameter{$parameter_name};

## Given Hash entries when files exists
PATH:
for my $path ( keys %{$parameter} ) {

    parse_meta_file_suffixes(
        {
            active_parameter_href  => \%active_parameter,
            file_name              => $path,
            meta_file_suffixes_ref => \@{ $file_info{$parameter_name} },
            parameter_href         => \%parameter,
            parameter_name         => $parameter_name,
        }
    );
}

## Then set build switch to false
is( $parameter{$parameter_name}{build_file},
    0, q{Set build file switch for hash parameter reference to 0} );

## Given Hash entries where files do not exist
%active_parameter = (
    exome_target_bed => {
        catfile( $Bin, qw{ data references does_not_exists.bed } ) => q{sample1},
    },
);
$parameter = $active_parameter{$parameter_name};

PATH:
for my $path ( keys %{$parameter} ) {

    parse_meta_file_suffixes(
        {
            active_parameter_href  => \%active_parameter,
            file_name              => $path,
            meta_file_suffixes_ref => \@{ $file_info{$parameter_name} },
            parameter_href         => \%parameter,
            parameter_name         => $parameter_name,
        }
    );
}

## Then set build file to true
is( $parameter{exome_target_bed}{build_file},
    1, q{Set build file switch for hash parameter reference to 1} );

done_testing();
