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
use MIP::Constants qw{ $COMMA $SPACE };
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
        q{MIP::File_info}      => [qw{ set_sample_file_attribute }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_sample_file_attribute };

diag(   q{Test set_sample_file_attribute from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a file to read
my %file_info;
my $file_name         = q{a_file.gz};
my $read_file_command = q{gzip -d -c};
my $sample_id         = q{a_sample};

## When file is compressed with gzip
set_sample_file_attribute(
    {
        attribute       => q{read_file_command},
        attribute_value => $read_file_command,
        file_info_href  => \%file_info,
        file_name       => $file_name,
        sample_id       => $sample_id,
    }
);

my $expected_read_file_command = q{gzip -d -c};

## Then set the read file command for file
is( $file_info{$sample_id}{$file_name}{read_file_command},
    $expected_read_file_command, q{Set sample file attribute for file} );

done_testing();
