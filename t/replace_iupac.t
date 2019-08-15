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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Program::Variantcalling::Perl} => [qw{ replace_iupac }],
        q{MIP::Test::Fixtures}                => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Variantcalling::Perl qw{ replace_iupac };

diag(   q{Test replace_iupac from Perl.pm v}
      . $MIP::Program::Variantcalling::Perl::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given xargs request
replace_iupac(
    {
        FILEHANDLE => $FILEHANDLE,
        xargs      => 1,
    }
);

## Close the filehandle
close $FILEHANDLE;

## Then write with escape chars
my ($returned_command_xargs) = $file_content =~ /(\'if)/xms;
ok( $returned_command_xargs, q{Wrote xargs perl command} );

## Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given xargs request
replace_iupac(
    {
        FILEHANDLE => $FILEHANDLE,
        xargs      => 0,
    }
);

## Close the filehandle
close $FILEHANDLE;

## Then write without escape chars
my ($returned_command) = $file_content =~ /('if)/xms;

ok( $returned_command, q{Wrote perl command} );

done_testing();
