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
use Modern::Perl qw{ 2014 };
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
Readonly my $COMMA                => q{,};
Readonly my $GENOME_BUILD_VERSION => 37;
Readonly my $NEWLINE              => qq{\n};
Readonly my $SPACE                => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::QC::Record}     => [qw{ add_to_sample_info }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::QC::Record qw{ add_to_sample_info };

diag(   q{Test add_to_sample_info from Record.pm v}
      . $MIP::QC::Record::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given parameter paths
my %active_parameter = (
    human_genome_reference => catfile(qw{ a test path genome build}),
    log_file               => catfile(qw{ a test dir and log_path}),
    pedigree_file          => catfile(qw{ a test pedigree path }),
);

my %file_info = (
    human_genome_reference_source  => q{GRCh},
    human_genome_reference_version => $GENOME_BUILD_VERSION,
);
my %sample_info;

add_to_sample_info(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        sample_info_href      => \%sample_info,
    }
);

## Then these entries should be set in sample info
is(
    $sample_info{human_genome_build}{path},
    $active_parameter{human_genome_reference},
    q{Added genome build path}
);
is(
    $sample_info{human_genome_build}{source},
    $file_info{human_genome_reference_source},
    q{Added genome build source}
);
is(
    $sample_info{human_genome_build}{version},
    $file_info{human_genome_reference_version},
    q{Added genome build version}
);
is(
    $sample_info{pedigree_file}{path},
    $active_parameter{pedigree_file},
    q{Added pedigree path}
);
is(
    $sample_info{log_file_dir},
    dirname( dirname( $active_parameter{log_file} ) ),
    q{Added log dir path}
);
is(
    $sample_info{last_log_file_path},
    $active_parameter{log_file},
    q{Added log file path}
);

done_testing();
