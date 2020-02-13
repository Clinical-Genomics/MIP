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
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Script::Setup_script} => [qw{ write_return_to_conda_environment }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Setup_script qw{ write_return_to_conda_environment };

diag(   q{Test write_return_to_conda_environment from Setup_script.pm v}
      . $MIP::Script::Setup_script::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given main environment array
my $source_main_environment_commands_ref = [q{conda activate test}];

# Create anonymous filehandle
my $filehandle = IO::Handle->new();

# For storing info to write
my $file_content;

# Store file content in memory by using referenced variable
open $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

write_return_to_conda_environment(
    {
        source_main_environment_commands_ref => $source_main_environment_commands_ref,
        filehandle                           => $filehandle,
    }
);

# Close the filehandle
close $filehandle;

## Then env load command shoudl be written to file
my ($load_command) = $file_content =~ /conda\sactivate\stest/xms;
ok( $load_command, q{Wrote env load command} );

## Given no main environment
$source_main_environment_commands_ref = [];

# Create anonymous filehandle
$filehandle = IO::Handle->new();

# For storing info to write
my $file_content_2;

# Store file content in memory by using referenced variable
open $filehandle, q{>}, \$file_content_2
  or croak q{Cannot write to} . $SPACE . $file_content_2 . $COLON . $SPACE . $OS_ERROR;

write_return_to_conda_environment(
    {
        source_main_environment_commands_ref => $source_main_environment_commands_ref,
        filehandle                           => $filehandle,
    }
);

# Close the filehandle
close $filehandle;

## Then env load command shoudl be written to file
my ($load_command_2) = $file_content_2 =~ /conda\sdeactivate/xms;
ok( $load_command_2, q{Wrote conda deactivate command} );

done_testing();
