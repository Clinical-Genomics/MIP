#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
use MIP::Constants qw{ $COLON $COMMA $SPACE };
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
        q{MIP::Script::Setup_script} => [qw{ create_script_error_trap }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Setup_script qw{ create_script_error_trap };

diag(   q{Test create_script_error_trap from Setup_script.pm v}
      . $MIP::Script::Setup_script::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir = File::Temp->newdir();
my $temp_dir = catfile($test_dir);

# For storing info to write
my $file_content;

# Store file content in memory by using referenced variable
open my $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

my %active_parameter = (
    log_file            => catfile( $temp_dir, q{log.status} ),
    sacct_format_fields => [qw{ jobid }],
    submission_profile  => q{slurm},
);
my %job_id = ( ALL => { ALL => [ qw{ jobid_1 }, ], }, );

## When error trap is zero
my $return = create_script_error_trap(
    {
        error_trap              => 0,
        filehandle              => $filehandle,
        job_ids_ref             => $job_id{ALL}{ALL},
        log_file_path           => $active_parameter{log_file},
        sacct_format_fields_ref => $active_parameter{sacct_format_fields},
    }
);

## Then the recipe zero should be returned
is( $return, 0, q{Skip sub - error_trap is false} );

## When error_trap is true
create_script_error_trap(
    {
        error_trap              => 1,
        filehandle              => $filehandle,
        job_ids_ref             => $job_id{ALL}{ALL},
        log_file_path           => $active_parameter{log_file},
        sacct_format_fields_ref => $active_parameter{sacct_format_fields},
    }
);

close $filehandle;

## Then instructions should be written to file to create error trap
ok( $file_content =~ / local \s+ program /xms, q{Create error trap} );

done_testing();
