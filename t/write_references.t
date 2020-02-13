#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
        q{MIP::File::Format::Reference} => [qw{ write_references }],
        q{MIP::Io::Read}                => [qw{ read_from_file }],
        q{MIP::Test::Fixtures}          => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Reference qw{ write_references };
use MIP::Io::Read qw{ read_from_file };

diag(   q{Test write_references from Reference.pm v}
      . $MIP::File::Format::Reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

my $test_dir     = File::Temp->newdir();
my $outfile_path = catfile( $test_dir, q{reference.yaml} );

## Given
my %active_parameter = (
    array_parameter => [qw{ ref_1 ref_2 }],
    hash_parameter  => {
        ref_1 => q{ref_3},
        ref_2 => q{ref_4}
    },
    not_a_reference  => q{lack_the_key},
    scalar_parameter => q{a_ref},
);
my %parameter = (
    array_parameter  => { is_reference => 1, },
    hash_parameter   => { is_reference => 1, },
    scalar_parameter => { is_reference => 1, },
);

write_references(
    {
        active_parameter_href => \%active_parameter,
        outfile_path          => $outfile_path,
        parameter_href        => \%parameter,
    }
);

## Then file should be created
ok( -e $outfile_path, q{Created reference file path} );

my %expected_reference = (
    array_parameter => [qw{ ref_1 ref_2 }],
    hash_parameter  => {
        ref_1 => q{ref_3},
        ref_2 => q{ref_4},
    },
    scalar_parameter => q{a_ref},
);
my %reference = read_from_file(
    {
        format => q{yaml},
        path   => $outfile_path,
    }
);

## Then only parameters with is_reference should be printed
is_deeply( \%reference, \%expected_reference, q{Wrote reference hash correctly} );

done_testing();
