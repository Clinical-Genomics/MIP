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
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Commands qw{ test_function };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Docker} => [qw{ docker_run }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Docker qw{ docker_run };

diag(   q{Test docker_run from Docker.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my @function_base_commands = qw{ docker run };

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
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = (
    image => {
        input           => q{docker.io/clinicalgenomics/bwa:0.7.17},
        expected_output => q{docker.io/clinicalgenomics/bwa:0.7.17},
    },
);

my %specific_argument = (
    bind_paths_ref => {
        inputs_ref => [ catdir(qw{ path host }) . $COLON . catdir(qw{ path container }) ],
        expected_output => q{--volume}
          . $SPACE
          . catdir(qw{ path host })
          . $COLON
          . catdir(qw{ path container }),
    },
    container_cmds_ref => {
        inputs_ref      => [q{Hello_world.py}],
        expected_output => qw{ Hello_world.py },
    },
    entrypoint => {
        input           => q{/bin/bash},
        expected_output => q{--entrypoint} . $SPACE . q{/bin/bash},
    },
    image => {
        input           => q{docker.io/clinicalgenomics/bwa:0.7.17},
        expected_output => q{docker.io/clinicalgenomics/bwa:0.7.17},
    },
    remove => {
        input           => 1,
        expected_output => q{--rm},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&docker_run;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {

    my @commands = test_function(
        {
            argument_href              => $argument_href,
            do_test_base_command       => 1,
            function_base_commands_ref => \@function_base_commands,
            module_function_cref       => $module_function_cref,
            required_argument_href     => \%required_argument,
        }
    );
}

done_testing();
