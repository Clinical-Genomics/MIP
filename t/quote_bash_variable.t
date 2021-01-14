#! /usr/bin/env perl

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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Language::Shell} => [qw{ quote_bash_variable }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Language::Shell qw{ quote_bash_variable };

diag(   q{Test quote_bash_variable from Shell.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a string with a variable to quote
my $return_string =
  quote_bash_variable( { string_with_variable_to_quote => q{$TEST_VARIABLE} } );

## Then return quoted variable
is( $return_string, q{"$TEST_VARIABLE"}, q{quote_bash_variable} );

done_testing();
