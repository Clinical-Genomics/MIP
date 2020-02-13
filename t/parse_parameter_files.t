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
        q{MIP::Parameter}      => [qw{ parse_parameter_files }],
        q{MIP::Io::Read}       => [qw{ read_from_file }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Io::Read qw{ read_from_file };
use MIP::Parameter qw{ parse_parameter_files };

diag(   q{Test parse_parameter_files from Parameter.pm v}
      . $MIP::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given
my %active_parameter = (
    gatk_baserecalibration_known_sites =>
      [ catfile( $Bin, qw{ data references grch37_dbsnp_-138-.vcf} ), ],
    mip                    => 1,
    gatk_baserecalibration => 1,
);
my $consensus_analysis_type = q{wgs};
my %parameter               = read_from_file(
    {
        format => q{yaml},
        path   => catfile( dirname($Bin), qw{ definitions rd_dna_parameters.yaml} ),
    }
);

my $is_ok = parse_parameter_files(
    {
        active_parameter_href   => \%active_parameter,
        consensus_analysis_type => $consensus_analysis_type,
        parameter_href          => \%parameter,
    }
);

## Then
ok( $is_ok, q{Parsed parameter files for existence} );

done_testing();
