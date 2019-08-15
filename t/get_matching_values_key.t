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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Get::File}      => [qw{ get_matching_values_key }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::File qw{ get_matching_values_key };

diag(   q{Test get_matching_values_key from File.pm v}
      . $MIP::Get::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a parameter name when not existing
my $sample_id = q{sample-1};
my %active_parameter;

my $infile_directory = get_matching_values_key(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => q{infile_dirs},
        query_value           => $sample_id,
    }
);

## Then return undef
is( $infile_directory, undef, q{Return for non-existing key} );

## Given matching sample_id
%active_parameter = ( infile_dirs => { q{for_sample-1} => q{sample-1}, }, );

$infile_directory = get_matching_values_key(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => q{infile_dirs},
        query_value           => $sample_id,
    }
);

## Then return the infile directory
is( $infile_directory, q{for_sample-1}, q{Returned infile_directory} );

## Given non-matching sample_id
$active_parameter{infile_dirs}{q{for_sample-1}} = q{no_match_sample_id};

$infile_directory = get_matching_values_key(
    {
        active_parameter_href => \%active_parameter,
        parameter_name        => q{infile_dirs},
        query_value           => $sample_id,
    }
);

## Then return undef
is( $infile_directory, undef, q{Returned undef with no matching sample_id} );

done_testing();
