#!/usr/bin/env perl

use 5.018;
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
use Modern::Perl qw{ 2014 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = '1.0.0';

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COLON => q{:};
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Build::Bwa_prerequisites} =>
          [qw{ build_bwa_prerequisites }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Build::Bwa_prerequisites qw{ build_bwa_prerequisites };

diag(   q{Test build_bwa_prerequisites from Bwa_prerequisites.pm v}
      . $MIP::Recipes::Build::Bwa_prerequisites::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log();

## Create temp logger
my $test_dir = File::Temp->newdir();

## Given build parameters
my $call_type            = q{both};
my $parameter_build_name = q{bwa_build_reference};
my $program_name         = q{bwa_mem};
my $temp_dir             = catfile($test_dir);
my $test_program_name    = q{bwa_mem};

my %active_parameter = (
    bash_set_errexit                => 1,
    bash_set_nounset                => 1,
    bash_set_pipefail               => 1,
    email_types                     => [qw{ BEGIN FAIL }],
    family_id                       => q{family_1},
    human_genome_reference          => q{human_genome.fasta},
    outaligner_dir                  => q{bwa},
    outdata_dir                     => catfile( $test_dir, q{test_data_dir} ),
    outscript_dir                   => catfile( $test_dir, q{test_script_dir} ),
    $program_name                   => 2,
    project_id                      => q{wamdu},
    sample_ids                      => [qw{ sample-1 sample-2 }],
    slurm_quality_of_service        => q{low},
    source_environment_commands_ref => [qw{ source activate test }],
    temp_directory                  => $temp_dir,
);
my %file_info = (
    human_genome_compressed => q{not_compressed},
    $parameter_build_name   => [qw{ .bwt .ann .amb .pac .sa }],
);
my %infile_lane_prefix;
my %job_id;
my %parameter = (
    human_genome_reference => { build_file => 0, },
    $parameter_build_name  => { build_file => 1, },
    $program_name          => { chain      => q{TEST}, },
);
my %sample_info;

trap {
    build_bwa_prerequisites(
        {
            active_parameter_href => \%active_parameter,
            bwa_build_reference_file_endings_ref =>
              \@{ $file_info{$parameter_build_name} },
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            log                     => $log,
            parameter_href          => \%parameter,
            program_name            => $program_name,
            sample_info_href        => \%sample_info,
        }
      )
};

## Then broadcast info log message
my $log_msg =
  q{Will\s+try\s+to\s+create\s+required\s+human_genome.fasta\s+index\s+files};
like( $trap->stderr, qr/$log_msg/msx, q{Broadcast bwa build log message} );

done_testing();
