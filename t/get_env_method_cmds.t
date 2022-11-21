#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
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
    my %perl_module = ( q{MIP::Environment::Manager} => [qw{ get_env_method_cmds }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Manager qw{ get_env_method_cmds };

diag(   q{Test get_env_method_cmds from Parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given env method, name and load action
my $env_method = q{conda};
my $env_name   = q{test_env};

my @env_cmds = get_env_method_cmds(
    {
        action          => q{load},
        conda_init_path => catfile(qw{path to conda etc profile.d conda.sh}),
        env_name        => $env_name,
        env_method      => $env_method,
    }
);
my @expected_load_cmds = (
    q{source},
    catfile(qw{path to conda etc profile.d conda.sh}),
    qw{ ; conda activate }, $env_name,
);

## Then return load env commands
is_deeply( \@expected_load_cmds, \@env_cmds, q{Got env load cmds} );

## Given env method, name and unload action
@env_cmds = get_env_method_cmds(
    {
        action     => q{unload},
        env_name   => $env_name,
        env_method => $env_method,
    }
);

my @expected_unload_cmds = (qw{ conda deactivate });

## Then return unload env commands
is_deeply( \@expected_unload_cmds, \@env_cmds, q{Got env unload cmds} );

done_testing();
