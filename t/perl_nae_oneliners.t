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
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Language::Perl} => [qw{ perl_nae_oneliners }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Language::Perl qw{ perl_nae_oneliners };

diag(   q{Test perl_nae_oneliners from Perl.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ perl };

my %base_argument = (
    filehandle => {
        input           => undef,
        expected_output => \@function_base_commands,
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    stderrfile_path_append => {
        input           => q{stderrfile.test},
        expected_output => q{2>> stderrfile.test},
    },
    stdoutfile_path => {
        input           => q{stdoutfile.test},
        expected_output => q{1> stdoutfile.test},
    },
    stdoutfile_path_append => {
        input           => q{stdoutfile.test},
        expected_output => q{1>> stdoutfile.test},
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = ();

my %specific_argument = (
    oneliner_cmd => {
        input           => q?'if($_=~/ensembl-vep\s+:\s(\d+)/xms) {print $1;}'?,
        expected_output => q?'if($_=~/ensembl-vep\s+:\s(\d+)/xms) {print $1;}'?,
    },
    oneliner_name => {
        input => q{synonyms_grch37_to_grch38},
        expected_output =>
          q?'if($_=~s/^MT/chrM/g) {} elsif ($_=~s/^([^#])/chr$1/g) {} print $_'?,
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&perl_nae_oneliners;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {

    test_function(
        {
            argument_href              => $argument_href,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
        }
    );
}

## Given more oneliner_names
my %oneliner_map = ( synonyms_grch38_to_grch37 =>
      q?'if($_=~s/^chrM/MT/g) {} elsif ($_=~s/^chr(.+)/$1/g) {} print $_'?, );

ONELINER:
while ( my ( $name, $oneliner ) = each %oneliner_map ) {
    $specific_argument{oneliner_name}{input}           = $name;
    $specific_argument{oneliner_name}{expected_output} = $oneliner;

    test_function(
        {
            argument_href              => \%specific_argument,
            do_test_base_command       => 0,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
        }
    );
}

## Given escape argument
my @perl_cmds = perl_nae_oneliners(
    {
        escape_oneliner => 1,
        oneliner_name   => q{synonyms_grch38_to_grch37},
    }
);

## Then escape the perl expression
my @expected_oneliner = (
    qw{ perl -n -a -e },
    q?\'if($_=~s/^chrM/MT/g) {} elsif ($_=~s/^chr(.+)/$1/g) {} print $_\'?
);
is_deeply( \@perl_cmds, \@expected_oneliner, q{Escape oneliner} );

done_testing();
