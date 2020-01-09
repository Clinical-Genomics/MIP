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
use MIP::Constants qw{ $COMMA $SPACE $TAB };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };
use Test::Trap;

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::File::Format::Vcf} => [qw{ check_vcf_variant_line }],
        q{MIP::Test::Fixtures}    => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Vcf qw{ check_vcf_variant_line };

diag(   q{Test check_vcf_variant_line from Vcf.pm v}
      . $MIP::File::Format::Vcf::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $INFO_COL_NR => 7;

my $log = test_log( {} );

## Given a correct variant line
my $input_line_number = 1;
my @line_elements     = qw{ 1 1 rsid1 A G . PASS AF=1 GT  0/1 };
my $variant_line      = join $TAB, @line_elements;
my $is_ok             = check_vcf_variant_line(
    {
        input_line_number         => $input_line_number,
        log                       => $log,
        variant_line              => $variant_line,
        variant_line_elements_ref => \@line_elements,
    }
);

## Then return true
ok( $is_ok, q{Correct line} );

## Given a variant line when no INFO
$line_elements[$INFO_COL_NR] = undef;
$variant_line = join $TAB, grep { defined } @line_elements;

trap {
    check_vcf_variant_line(
        {
            input_line_number         => $input_line_number,
            log                       => $log,
            variant_line              => $variant_line,
            variant_line_elements_ref => \@line_elements,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the INFO field cannot be found} );
like( $trap->stderr, qr/No\s+INFO\s+field\s+at/xms, q{Throw fatal log message} );

done_testing();
