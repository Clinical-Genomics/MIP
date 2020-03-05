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
        q{MIP::Set::File}      => [qw{ set_file_compression_features }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::File qw{ set_file_compression_features };

diag(   q{Test set_file_compression_features from File.pm v}
      . $MIP::Set::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given uncompressed file
my $file_name = q{a_file.fastq};

my ( $is_file_compressed, $read_file_command ) =
  set_file_compression_features( { file_name => $file_name, } );

## Then return false
is( $is_file_compressed, 0, q{File was not compressed} );

## Then set read command to handle uncompressed file
is( $read_file_command, q{cat}, q{Set read file command for uncompressed file} );

## Given compressed file
$file_name .= q{.gz};

( $is_file_compressed, $read_file_command ) =
  set_file_compression_features( { file_name => $file_name, } );

## Then return true
ok( $is_file_compressed, q{File was compressed} );

## Then set read command to handle uncompressed file
is( $read_file_command, q{gzip -d -c}, q{Set read file command for compressed file} );

done_testing();
