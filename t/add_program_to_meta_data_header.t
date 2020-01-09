#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use DateTime;
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
        q{MIP::Vcfparser}      => [qw{ add_program_to_meta_data_header }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ add_program_to_meta_data_header };

diag(   q{Test add_program_to_meta_data_header from Vcfparser.pm v}
      . $MIP::Vcfparser::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a vcfparser version and no addition of software tag
my %meta_data;
my $current_date      = DateTime->now->ymd;
my $program_name      = basename($PROGRAM_NAME);
my $vcfparser_version = 1.0;
my $header_line =
  qq{##Software=<ID=$program_name,Version=$vcfparser_version,Date=$current_date};

add_program_to_meta_data_header(
    {
        meta_data_href    => \%meta_data,
        vcfparser_version => $vcfparser_version,
    }
);

## Then software key should not exists
is( $meta_data{software}, undef, q{Skip addition of software tag} );

## Given a vcfparser version and addition of software tag
add_program_to_meta_data_header(
    {
        add_software_tag  => 1,
        meta_data_href    => \%meta_data,
        vcfparser_version => $vcfparser_version,
    }
);

## Then software key should exists
is( $meta_data{software}{$program_name}, $header_line, q{Added software tag to hash} );

done_testing();
