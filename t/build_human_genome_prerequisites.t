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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Recipes::Build::Human_genome_prerequisites} => [qw{ build_human_genome_prerequisites }],
        q{MIP::Test::Fixtures} =>
          [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Build::Human_genome_prerequisites qw{ build_human_genome_prerequisites };

diag(   q{Test build_human_genome_prerequisites from Human_genome_prerequisites.pm v}
      . $MIP::Recipes::Build::Human_genome_prerequisites::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log();

## Given build parameters
my $parameter_build_name = q{human_genome_reference_file_endings};
my $program_name = q{bwa_mem};

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        program_name  => $program_name,
    }
);
my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        program_name  => $program_name,
    }
);
$file_info{human_genome_reference_name_prefix} = q{human_genome};

my %infile_lane_prefix;
my %job_id;
my %parameter = test_mip_hashes( { mip_hash_name => q{parameter}, } );

my %sample_info;

trap {
    build_human_genome_prerequisites(
        {
            active_parameter_href   => \%active_parameter,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            log                     => $log,
	 parameter_build_suffixes_ref =>
	 \@{ $file_info{$parameter_build_name} },
            parameter_href          => \%parameter,
            program_name            => $program_name,
            sample_info_href        => \%sample_info,
        }
      )
};

## Then broadcast info log message
my $log_msg = q{Will\s+try\s+to\s+create\s+.dict};
like( $trap->stderr, qr/$log_msg/msx, q{Broadcast human genome build log message} );

done_testing();
