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
use MIP::Constants qw{ $COLON $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Gatk}  => [qw{ gatk_concatenate_variants }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gatk qw{ gatk_concatenate_variants };

diag(   q{Test gatk_concatenate_variants from Gatk.pm}
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

## Given files and input
my %active_parameter = (
    java_use_large_pages => 1,
    gatk_logging_level   => q{INFO},
    temp_directory       => catdir(qw{ a test dir}),
);
my @contigs        = qw{ 1 2 3 };
my $infile_prefix  = q{infile};
my $outfile_suffix = q{.vcf};

gatk_concatenate_variants(
    {
        active_parameter_href => \%active_parameter,
        continue              => 1,
        filehandle            => $filehandle,
        elements_ref          => \@contigs,
        infile_prefix         => $infile_prefix,
        outfile_suffix        => $outfile_suffix,
    }
);

## Close the filehandle
close $filehandle;

## Then write concatenate instructions
my ($returned_command) = $file_content =~ /(GatherVcfsCloud)/xms;
ok( $returned_command, q{Wrote concatenate instructions} );

## Then also write ampersand
my ($wrote_ampersand) = $file_content =~ /(&)/xms;
ok( $wrote_ampersand, q{Wrote continue instructions} );

done_testing();
