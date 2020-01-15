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
use MIP::Constants
  qw{ $COMMA $DOUBLE_QUOTE $EQUALS $SEMICOLON set_singularity_constants @SINGULARITY_BIND_PATHS $SPACE };
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
        q{MIP::Parse::Singularity} => [qw{ parse_sing_bind_paths }],
        q{MIP::Test::Fixtures}     => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Singularity qw{ parse_sing_bind_paths };

diag(   q{Test parse_sing_bind_paths from Singularity.pm v}
      . $MIP::Parse::Singularity::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my %active_parameter = (
    outdata_dir                  => q{an_outdata_dir},
    reference_dir                => q{a_ref_dir},
    singularity_recipe_bind_path => { bwa_mem => [qw{ infiles_dir }], },
    temp_directory               => q{a_temp_dir},
    with_singularity             => 1,
);

my @source_environment_cmds = qw{ conda activate test };
set_singularity_constants( { active_parameter_href => \%active_parameter, } );

parse_sing_bind_paths(
    {
        active_parameter_href       => \%active_parameter,
        package_name                => q{bwa_mem},
        source_environment_cmds_ref => \@source_environment_cmds,
    }
);

## Constant dir paths and specific one(s)
my $singularity_bind = join $COMMA, ( @SINGULARITY_BIND_PATHS, q{infiles_dir} );

my @expected_source_environment_cmds = (
    qw{ conda activate test },
    q{export SINGULARITY_BIND}
      . $EQUALS
      . $DOUBLE_QUOTE
      . $singularity_bind
      . $DOUBLE_QUOTE
      . $SEMICOLON
);

## Then the recipe specific list of paths should ahve been added to the source env comd
is_deeply(
    \@source_environment_cmds,
    \@expected_source_environment_cmds,
    q{Bound extra paths}
);

done_testing();
