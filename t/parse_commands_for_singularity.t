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
        q{MIP::Parse::Singularity} => [qw{ parse_commands_for_singularity }],
        q{MIP::Test::Fixtures}     => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Singularity qw{ parse_commands_for_singularity };
use MIP::Constants qw{ set_analysis_constants };

diag(   q{Test parse_commands_for_singularity from Parse::Singularity.pm v}
      . $MIP::Parse::Singularity::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %active_parameter = (
    infile_dirs => {
        catdir(qw{ infile dir }) => q{sample_id},
    },
    outdata_dir           => catdir( dirname($Bin), qw{ t data } ),
    pedigree_file         => catfile( dirname($Bin), qw{ t data pedigree } ),
    reference_dir         => catdir( dirname($Bin), qw{ t data references } ),
    singularity_container => {
        q{test.py} => {
            container_path   => catfile(qw{ path to test.sif }),
            extra_bind_paths => [ catfile( dirname($Bin), qw{ t data input.bam } ) ],
        },
    },
    temp_directory   => catdir(qw{ temp folder }),
    with_singularity => 1,
);

## Make constants
set_analysis_constants( { active_parameter_href => \%active_parameter, } );

## Given some commands were some are executables within a singularity container.
my @commands = (qw{ test.py });

## Then construct and prepend singularity command
my @expected = (
    qw{ singularity exec },
    q{--bind }
      . catdir(qw{ infile dir })
      . $COMMA
      . catdir(qw{ temp folder })
      . $COMMA
      . catdir( dirname($Bin), qw{ t data } ),
    catfile(qw{ path to test.sif }),
    @commands,
);

parse_commands_for_singularity( { commands_ref => \@commands, } );

is_deeply( \@commands, \@expected, q{Return parsed singularity command} );

done_testing();
