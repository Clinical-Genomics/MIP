#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use File::Temp;
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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Language::Shell} => [qw{ check_mip_process_paths }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Language::Shell qw{ check_mip_process_paths };
use MIP::Environment::Child_process qw{ child_process };

diag(   q{Test check_mip_process_paths from Shell.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir       = File::Temp->newdir();
my $test_file_path = catfile( $test_dir, q{test.sh} );

open my $filehandle, q{>}, $test_file_path
  or croak q{Cannot write to} . $SPACE . $test_file_path . $COLON . $SPACE . $OS_ERROR;

## Given path that exists and one that does not
my @paths = ( $test_file_path, catfile( cwd(), q{does_not_exist.test} ) );

my $is_ok = check_mip_process_paths(
    {
        filehandle => $filehandle,
        paths_ref  => \@paths,
    }
);

close $filehandle;

## Then bash function should be written
ok( $is_ok, q{Wrote bash test} );

my %process_return = child_process(
    {
        commands_ref => [ q{bash } . $test_file_path, ],
        process_type => q{open3},
    }
);

## Then "test.sh" file should be found
like( $process_return{stdouts_ref}[0], qr/Found\s+file/sxm, q{Found file test} );

## Then does_not_exist.test file should not be found
like( $process_return{stderrs_ref}[0],
    qr/Could\s+not\s+find/sxm, q{Could not find file test} );

done_testing();
