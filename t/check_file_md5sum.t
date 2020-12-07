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
        q{MIP::Validate::File} => [qw{ check_file_md5sum }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Validate::File qw{ check_file_md5sum };

diag(   q{Test check_file_md5sum from File.pm v}
      . $MIP::Validate::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create anonymous filehandle
my $filehandle = IO::Handle->new();

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## Given an undefined method
my $method;
my $outfile_path  = catfile(qw{ a dir file.fastq.gz });
my $md5_file_path = catfile(qw{ a dir file.fastq.gz.md5 });

my $return = check_file_md5sum(
    {
        filehandle    => $filehandle,
        md5_file_path => $md5_file_path,
        check_method  => $method,
    }
);

## Then return undef
is( $return, undef, q{Undef check method - skip} );

## Given a defined method, when no md5 suffix
$method = q{md5sum};
$return = check_file_md5sum(
    {
        filehandle    => $filehandle,
        md5_file_path => $outfile_path,
        check_method  => $method,
    }
);

## Then return undef
is( $return, undef, q{Wrong file suffix - skip} );

## Given a defined method, when a md5 suffix exists
my $is_ok = check_file_md5sum(
    {
        filehandle    => $filehandle,
        md5_file_path => $md5_file_path,
        check_method  => $method,
    }
);

## Close the filehandle
close $filehandle;

## Then return TRUE
ok( $is_ok, q{Checked file} );

## Then perl should be in file_content
my $base_command = q{perl};
my ($returned_base_command) = $file_content =~ /^($base_command)/xms;
is( $returned_base_command, $base_command, q{Wrote perl regexp} );

## Then md5sum should be in file content
$base_command = q{md5sum};
($returned_base_command) = $file_content =~ /^($base_command)/xms;
is( $returned_base_command, $base_command, q{Wrote md5sum check} );

done_testing();
