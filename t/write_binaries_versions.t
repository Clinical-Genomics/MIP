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
use MIP::Constants qw{ $COLON $COMMA $EMPTY_STR $SPACE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Environment::Executable} => [qw{ write_binaries_versions }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Executable qw{ write_binaries_versions };

diag(   q{Test write_binaries_versions from Executable.pm}
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

## Given an existing perl command in CONTAINER_CMD
my $container_perl_command = q{singularity exec docker.io/clinicalgenomics/perl:5.26 perl};
my $container_vep_command  = q{singularity exec docker.io/ensemblorg/ensembl-vep:release_103.1 vep};

my %container_cmd = (
    perl => $container_perl_command,
    vep  => $container_vep_command,
);

## Given an executable name
my $outfile_path = q{an_outfile};

write_binaries_versions(
    {
        binary_info_href => \%container_cmd,
        filehandle       => $filehandle,
        outfile_path     => $outfile_path,
    }
);

close $filehandle;

## Then return version command
my ($returned_base_command) = $file_content =~ /($outfile_path)/mxs;
is( $returned_base_command, $outfile_path, q{Wrote binary version command} );

done_testing();
