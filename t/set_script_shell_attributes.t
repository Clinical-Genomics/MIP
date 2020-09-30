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
        q{MIP::Script::Setup_script} => [qw{ set_script_shell_attributes }],
        q{MIP::Test::Fixtures}       => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Setup_script qw{ set_script_shell_attributes };

diag(   q{Test set_script_shell_attributes from Setup_script.pm v}
      . $MIP::Script::Setup_script::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# For storing info to write
my $file_content;

# Store file content in memory by using referenced variable
open my $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## When calling sub
set_script_shell_attributes(
    {
        filehandle   => $filehandle,
        set_errexit  => 1,
        set_nounset  => 1,
        set_pipefail => 1,
        sleep        => 1,
        ulimit_n     => 2,
    }
);

close $filehandle;

## Then instructions should be written to file
ok( $file_content =~ / readonly \s+ PROGNAME /xms, q{Wrote shell attributes} );

done_testing();
