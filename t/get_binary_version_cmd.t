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
use MIP::Constants qw{ $COMMA $EMPTY_STR $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Constants}               => [qw{ set_container_cmd }],
        q{MIP::Environment::Executable} => [qw{ get_binary_version_cmd }],
        q{MIP::Test::Fixtures}          => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Constants qw{ set_container_cmd };
use MIP::Environment::Executable qw{ get_binary_version_cmd };

diag(   q{Test get_binary_version_cmd from Executable.pm v}
      . $MIP::Environment::Executable::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an executable name
my @version_cmd = get_binary_version_cmd(
    {
        binary     => q{vep},
        binary_cmd => q{vep}
    }
);

my @expected_version_cmd = (
    qw{ vep | perl -n -a -e },
    q?'my ($version) = /ensembl-vep\s+:\s(\d+)/xms; if($version) {print $version;last;}'?,
);

## Then return version command
is_deeply( \@version_cmd, \@expected_version_cmd, q{Got executable version command} );

## Given an existing perl command in CONTAINER_CMD
my $container_base_command =
  q{singularity exec docker.io/clinicalgenomics/perl:5.26 perl};
my %container_cmd = ( perl => $container_base_command, );

set_container_cmd( { container_cmd_href => \%container_cmd, } );

## When calling using container
@version_cmd = get_binary_version_cmd(
    {
        binary        => q{vep},
        binary_cmd    => q{vep},
        use_container => 1,
    }
);

my @expected_version_cmd_container = (
    qw{ vep | },
    q{singularity exec docker.io/clinicalgenomics/perl:5.26 perl},
    qw{ -n -a -e },
    q?'my ($version) = /ensembl-vep\s+:\s(\d+)/xms; if($version) {print $version;last;}'?,
);

## Then return version command
is_deeply(
    \@version_cmd,
    \@expected_version_cmd_container,
    q{Got executable version command with perl container}
);

## Given an executable name and stdoutfile path to append to
my $stdoutfile_path = q{an_outfile_path};
@version_cmd = get_binary_version_cmd(
    {
        binary                 => q{vep},
        binary_cmd             => q{vep},
        stdoutfile_path_append => $stdoutfile_path,
    }
);

my @expected_version_cmd_stdout = (
    qw{ vep | perl -n -a -e },
    q?'my ($version) = /ensembl-vep\s+:\s(\d+)/xms; if($version) {print $version;last;}'?,
    q{1>>} . $SPACE . $stdoutfile_path,
);

## Then return version command
is_deeply(
    \@version_cmd,
    \@expected_version_cmd_stdout,
    q{Got executable version command appending to stdout}
);
done_testing();
