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
        q{MIP::Script::Setup_script} => [qw{ create_script_temp_dir }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Setup_script qw{ create_script_temp_dir };

diag(   q{Test create_script_temp_dir from Setup_script.pm v}
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
my %job_id = ( PAN => { PAN => [ qw{ jobid_1 }, ], }, );

## When no temp directory supplied
my $return = create_script_temp_dir(
    {
        filehandle              => $filehandle,
        job_ids_ref             => $job_id{PAN}{PAN},
        log_file_path           => $active_parameter{log_file},
        sacct_format_fields_ref => $active_parameter{sacct_format_fields},
    }
);

## Then the recipe file path should be returned
is( $return, 0, q{Skip sub - no temp directory supplied} );

## Given a temp directory
create_script_temp_dir(
    {
        filehandle              => $filehandle,
        job_ids_ref             => $job_id{PAN}{PAN},
        log_file_path           => $active_parameter{log_file},
        sacct_format_fields_ref => $active_parameter{sacct_format_fields},
        temp_directory          => $temp_dir,
    }
);

close $filehandle;

## Then instructions should be written to file
ok( $file_content =~ / \A [#]{2} \s+ Create \s+ temporary \s+ directory /xms,
    q{Create temp dir} );

done_testing();
