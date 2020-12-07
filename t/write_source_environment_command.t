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
        q{MIP::Environment::Manager} => [qw{ write_source_environment_command }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Manager qw{ write_source_environment_command };

diag(   q{Test write_source_environment_command from Manager.pm v}
      . $MIP::Environment::Manager::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a filehandle
my $filehandle = IO::Handle->new();

## When no environment command
my $return = write_source_environment_command(
    {
        filehandle => $filehandle,
    }
);

## Then return zero
is( $return, 0, q{Skipped writing of command} );

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given a filehandle and a environment command
write_source_environment_command(
    {
        filehandle                      => $filehandle,
        source_environment_commands_ref => [qw{ conda activate test }],
    }
);

# Close the filehandle
close $filehandle;

## Then env comment and env should be written to file
my ($title) = $file_content =~ /\A ## Activate environment \z/mxs;

ok( $title, q{Wrote environment title} );

my ($write_source_environment_command) = $file_content =~ /^(conda\s+activate\s+test)/mxs;

ok( $write_source_environment_command, q{Wrote environment command} );

done_testing();
