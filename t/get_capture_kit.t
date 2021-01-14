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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Parameter}      => [qw{ get_capture_kit }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parameter qw{ get_capture_kit };

diag(   q{Test get_capture_kit from Parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Test defined switch but no user info - default ie. "latest" capture kit
my $capture_kit = q{latest};

my %parameter = (
    supported_capture_kit => {
        q{agilent_sureselect.v5} =>
          q{genome_reference_source_version_agilent_sureselect_targets_-v5-.bed},
        latest =>
          q{genome_reference_source_version_agilent_sureselect_targets_cre_-v1-.bed},
    },
);

my %is_user_supplied = ( exome_target_bed => 0, );

my $exome_target_bed_file = get_capture_kit(
    {
        capture_kit                => $capture_kit,
        is_set_by_user             => $is_user_supplied{exome_target_bed},
        supported_capture_kit_href => $parameter{supported_capture_kit},
    }
);
is(
    $exome_target_bed_file,
    q{genome_reference_source_version_agilent_sureselect_targets_cre_-v1-.bed},
    q{Got latest default capture kit}
);

## Test defined switch but no user info - no default capture kit
$capture_kit           = q{agilent_sureselect.v4};
$exome_target_bed_file = get_capture_kit(
    {
        capture_kit                => $capture_kit,
        is_set_by_user             => $is_user_supplied{exome_target_bed},
        supported_capture_kit_href => $parameter{supported_capture_kit},
    }
);
is( $exome_target_bed_file, q{agilent_sureselect.v4}, q{Got supplied capture kit} );

## Test undefined switch but no user info - default capture kit
$capture_kit      = q{agilent_sureselect.v5};
%is_user_supplied = ();

$exome_target_bed_file = get_capture_kit(
    {
        capture_kit                => $capture_kit,
        is_set_by_user             => $is_user_supplied{exome_target_bed},
        supported_capture_kit_href => $parameter{supported_capture_kit},
    }
);
is(
    $exome_target_bed_file,
    q{genome_reference_source_version_agilent_sureselect_targets_-v5-.bed},
    q{Got capture kit with undef user switch}
);

## ## Test undefined switch but no user info - no default capture kit
$capture_kit           = q{agilent_sureselect.v4};
$exome_target_bed_file = get_capture_kit(
    {
        capture_kit                => $capture_kit,
        is_set_by_user             => $is_user_supplied{exome_target_bed},
        supported_capture_kit_href => $parameter{supported_capture_kit},
    }
);
is( $exome_target_bed_file, q{agilent_sureselect.v4},
    q{Got supplied capture kit with undef user switch} );

## Test defined switch and user info - Get no cpature kit
%is_user_supplied = ( exome_target_bed => 1, );

$exome_target_bed_file = get_capture_kit(
    {
        capture_kit                => $capture_kit,
        is_set_by_user             => $is_user_supplied{exome_target_bed},
        supported_capture_kit_href => $parameter{supported_capture_kit},
    }
);
is( $exome_target_bed_file, undef, q{Did not get any capture kit} );

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
