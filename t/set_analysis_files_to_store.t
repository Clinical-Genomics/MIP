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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Store} =>
          [qw{ define_analysis_files_to_store set_analysis_files_to_store }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Store qw{ define_analysis_files_to_store set_analysis_files_to_store };

diag(   q{Test set_analysis_files_to_store from Store.pm v}
      . $MIP::Store::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a path
my %active_parameter = (
    case_id              => q{a_case_id},
    config_file          => catfile(qw{ a config_file }),
    config_file_analysis => catfile(qw{ a config_analysis_file }),
    log_file             => catfile(qw{ a log_file }),
    pedigree_file        => catfile(qw{ a pedigree_file }),
    pedigree_fam_file    => catfile(qw{ a pedigree_fam_file }),
    reference_info_file  => catfile(qw{ a reference_info_file }),
    sample_info_file     => catfile(qw{ a sample_info_file }),
);

my %sample_info;

## Then path should be recorded under store in sample_info
set_analysis_files_to_store(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
);

my %expected_sample_info = (
    files => [
        {
            format     => q{meta},
            id         => $active_parameter{case_id},
            path       => $active_parameter{config_file},
            path_index => undef,
            step       => q{mip_analyse},
            tag        => q{config},
        },
    ],
);

while ( my ( $file_href_index, $file_href ) = each @{ $sample_info{files} } ) {

    if ( $expected_sample_info{files}[0]{tag} eq $file_href->{tag} ) {

        is_deeply(
            $file_href,
            \%{ $expected_sample_info{files}[0] },
            q{Set analysis files}
        );
    }
}
done_testing();
