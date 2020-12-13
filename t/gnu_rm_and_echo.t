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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Gnu::Coreutils} => [qw{ gnu_rm_and_echo }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gnu::Coreutils qw{ gnu_rm_and_echo };

diag(   q{Test gnu_rm_and_echo from Coreutils.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open my $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given a hash with files and new content
my $file_string = q{Test};
my %file        = (
    file_1 => q{Test},
    file_2 => q{Test},
);

## Then write rm and echo
my $write_command = gnu_rm_and_echo(
    {
        file_href  => \%file,
        filehandle => $filehandle,
        force      => 1,
    }
);
close $filehandle;

ok( $write_command, q{Write rm and echo} );

done_testing();
