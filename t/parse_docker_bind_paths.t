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
  qw{ $COLON $COMMA $DOUBLE_QUOTE $EQUALS $SEMICOLON set_container_constants @CONTAINER_BIND_PATHS $SPACE };
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
        q{MIP::Environment::Container} => [qw{ parse_container_bind_paths }],
        q{MIP::Test::Fixtures}         => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Container qw{ parse_container_bind_paths };

diag(   q{Test parse_container_bind_paths from Container.pm v}
      . $MIP::Environment::Container::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my %active_parameter = (
    outdata_dir       => q{an_outdata_dir},
    reference_dir     => q{a_ref_dir},
    recipe_bind_path  => { bwa_mem => [qw{ infiles_dir }], },
    temp_directory    => q{a_temp_dir},
    container_manager => q{docker},
);

my @source_environment_cmds = qw{ conda activate test };
set_container_constants( { active_parameter_href => \%active_parameter, } );

parse_container_bind_paths(
    {
        active_parameter_href       => \%active_parameter,
        package_name                => q{bwa_mem},
        source_environment_cmds_ref => \@source_environment_cmds,
    }
);

my @docker_binds = map { $_ . $COLON . $_ } ( @CONTAINER_BIND_PATHS, q{infiles_dir} );
my $docker_bind  = join $SPACE . q{--volume} . $SPACE, @docker_binds;

my @expected_source_environment_cmds = (
    qw{ conda activate test },
    q{export MIP_BIND}
      . $EQUALS
      . $DOUBLE_QUOTE
      . $docker_bind
      . $DOUBLE_QUOTE
      . $SEMICOLON
);

## Then set MIP_BIND for docker
is_deeply(
    \@source_environment_cmds,
    \@expected_source_environment_cmds,
    q{Got docker bind paths}
);
done_testing();
