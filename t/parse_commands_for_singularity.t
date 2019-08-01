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
use Modern::Perl qw{ 2017 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
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
        q{MIP::Parse::Parameter} => [qw{ parse_commands_for_singularity }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Parameter qw{ parse_commands_for_singularity };
use MIP::Constants qw{ set_analysis_constants };

diag(   q{Test parse_commands_for_singularity from Parse::Parameter.pm v}
      . $MIP::Parse::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %active_parameter = (
    with_singularity      => 1,
    singularity_container => {
        q{TIDDIT.py} => q{TIDDIT.simg},
    },
);

## Make constants
set_analysis_constants( { active_parameter_href => \%active_parameter, } );

## Given some commands were some are executables within a singularity container.
my @commands = (
    qw{ TIDDIT.py --sv },
    q{-p 6},
    q{-o /output/path},
    q{--ref }
      . catfile( dirname($Bin), qw{ t data references grch37_homo_sapiens_-d5-.fasta } ),
    q{--bam } . catdir( dirname($Bin), qw{ t data input.bam } ),
);

## Then construct and prepend singularity command
my @expected = (
    qw{ singularity exec },
    q{--bind } . catfile( dirname($Bin), qw{ t data } ),
    q{TIDDIT.simg}, @commands,
);

parse_commands_for_singularity( { commands_ref => \@commands, } );

is_deeply( \@commands, \@expected, q{Return paths} );

done_testing();
