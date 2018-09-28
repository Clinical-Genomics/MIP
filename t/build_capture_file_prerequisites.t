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
        q{MIP::Recipes::Build::Capture_file_prerequisites} =>
          [qw{ build_capture_file_prerequisites }],
        q{MIP::Test::Fixtures} =>
          [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Build::Capture_file_prerequisites
  qw{ build_capture_file_prerequisites };

diag(
q{Test build_capture_file_prerequisites from Capture_file_prerequisites.pm v}
      . $MIP::Recipes::Build::Capture_file_prerequisites::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log();

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to}
  . $SPACE
  . $file_content
  . $COLON
  . $SPACE
  . $OS_ERROR;

## Given build parameters
my $parameter_build_name = q{exome_target_bed};
my $program_name         = q{picardtools_collecthsmetrics};

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
my %infile_lane_prefix;
my %job_id;
my %parameter = test_mip_hashes( { mip_hash_name => q{recipe_parameter}, } );

my %sample_info;

my $infile_list_suffix          = $file_info{exome_target_bed}[0];
my $padded_infile_list_suffix   = $file_info{exome_target_bed}[1];
my $padded_interval_list_suffix = $file_info{exome_target_bed}[2];

trap {
    build_capture_file_prerequisites(
        {
            active_parameter_href   => \%active_parameter,
            FILEHANDLE              => $FILEHANDLE,
            file_info_href          => \%file_info,
            infile_lane_prefix_href => \%infile_lane_prefix,
            job_id_href             => \%job_id,
            log                     => $log,
            parameter_build_suffixes_ref =>
              \@{ $file_info{$parameter_build_name} },
            parameter_href   => \%parameter,
            program_name     => $program_name,
            sample_info_href => \%sample_info,
        }
      )
};

close $FILEHANDLE;

## Then broadcast info log message
my $log_msg = q{Will\s+try\s+to\s+create\s+required};
like( $trap->stderr, qr/$log_msg/msx,
    q{Broadcast exome_target_bed build log message} );

done_testing();
