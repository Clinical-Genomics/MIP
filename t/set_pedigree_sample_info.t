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
use MIP::Constants qw{ $COMMA $SPACE };
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
        q{MIP::Pedigree}       => [qw{ set_pedigree_sample_info }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ set_pedigree_sample_info };

diag(   q{Test set_pedigree_sample_info from Pedigree.pm v}
      . $MIP::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $TOTAL_SAMPLE_IDS => 4;
Readonly my $SAMPLE_3_INDEX   => 3;

my %active_parameter;

my %pedigree = (
    case    => q{case_1},
    samples => [
        {
            analysis_type => q{wes},
            capture_kit   => q{agilent_sureselect.v5},
            father        => 0,
            mother        => 0,
            phenotype     => q{affected},
            sample_id     => q{sample_1},
            sex           => q{female},
        },
        {
            analysis_type => q{wgs},
            capture_kit   => q{agilent_sureselect.v4},
            father        => 0,
            mother        => 0,
            phenotype     => q{unaffected},
            sample_id     => q{sample_2},
            sex           => q{male},
        },
        {
            analysis_type => q{wts},
            capture_kit   => q{agilent_sureselect.v5},
            father        => 0,
            mother        => 0,
            phenotype     => q{unknown},
            sample_id     => q{sample_3},
            sex           => q{other},
        },
        {
            analysis_type => q{wgs},
            father        => q{sample_1},
            mother        => q{sample_2},
            phenotype     => q{unknown},
            sample_id     => q{sample_4},
            sex           => q{unknown},
        },
    ],
);

my %sample_info = (
    sample => {
        sample_1 => {
            analysis_type     => q{wes},
            expected_coverage => 30,
        },
    },
);

my @user_input_sample_ids;
my %is_user_supplied = ( sample_ids => 0 );

my @pedigree_sample_ids = set_pedigree_sample_info(
    {
        active_parameter_href     => \%active_parameter,
        is_user_supplied_href     => \%is_user_supplied,
        pedigree_href             => \%pedigree,
        sample_info_href          => \%sample_info,
        user_input_sample_ids_ref => \@user_input_sample_ids,
    }
);

## Return of sample_ids
is( scalar @pedigree_sample_ids, $TOTAL_SAMPLE_IDS, q{Found all sample_ids} );

## Transfer to active parameters
is( scalar @{ $active_parameter{sample_ids} },
    $TOTAL_SAMPLE_IDS, q{Set all sample_ids to active parameter} );

## Sample level addition to sample info
SAMPLE_HREF:
foreach my $pedigree_sample_href ( @{ $pedigree{samples} } ) {

    # Alias
    my $sample_id = $pedigree_sample_href->{sample_id};

    foreach my $key ( keys %{$pedigree_sample_href} ) {

        is( $sample_info{sample}{$sample_id}{$key}, $pedigree_sample_href->{$key},
                q{Set key: }
              . $key
              . q{ to '}
              . $pedigree_sample_href->{$key}
              . q{' for '}
              . $sample_id
              . q{'} );
    }
}

## User input
@user_input_sample_ids = qw{ sample_1 sample_2 };
%is_user_supplied      = ( sample_ids => 1 );
%active_parameter      = ();
%sample_info           = ();

@pedigree_sample_ids = set_pedigree_sample_info(
    {
        active_parameter_href     => \%active_parameter,
        is_user_supplied_href     => \%is_user_supplied,
        pedigree_href             => \%pedigree,
        sample_info_href          => \%sample_info,
        user_input_sample_ids_ref => \@user_input_sample_ids,
    }
);

## Sample level addition to sample info for user supplied sample ids
foreach my $key ( keys %{ $sample_info{sample}{sample_1} } ) {

    is( $sample_info{sample}{sample_1}{$key}, $pedigree{samples}[0]{$key},
            q{Set key: }
          . $key
          . q{ to '}
          . $pedigree{samples}[0]{$key}
          . q{' for ' sample_1 '} );
}

## Sample level addition to sample info for not supplied sample ids
is( $sample_info{sample}{sample_3}{analysis_type}, undef,
        q{Did not set key: analysis_type to '}
      . $pedigree{samples}[$SAMPLE_3_INDEX]{analysis_type}
      . q{' for ' sample_3 '} );

is( $active_parameter{sample_ids}, undef, q{Did not set sample_ids to active parameter} );

done_testing();
