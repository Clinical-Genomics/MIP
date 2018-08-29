#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use Getopt::Long;
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
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE
          . basename($PROGRAM_NAME)
          . $SPACE
          . $VERSION
          . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Get::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Get::Parameter qw{ get_dynamic_conda_path get_program_version  };

diag(   q{Test get_program_version from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given parameter paths to use for regexps
my %active_parameter = (
    gatk_path        => catfile(qw{ a test path gatk-3.8}),
    picardtools_path => catfile(qw{ a test path picard-tools-2.14.1}),
);
my %sample_info;

## Define program features
my %program_feature = (
    gatk_path => {
        cmd          => q{To be disclosed},
        regexp       => q?gatk-([^,]+)?,
        program_name => q{gatk},
    },
    picardtools_path => {
        cmd          => q{To be disclosed},
        regexp       => q?picard-tools-([^,]+)?,
        program_name => q{picard.jar},
    },
    sambamba_depth => {
        cmd          => q{To be disclosed},
        regexp       => q?Not relevant?,
        program_name => q{sambamba},
    },
);
my %return;

## Given parameters to find version for, when valid regexps
PARAMETER:
foreach my $parameter_name ( keys %program_feature ) {

    my $version = get_program_version(
        {
            active_parameter_href => \%active_parameter,
            cmd                   => $program_feature{$parameter_name}{cmd},
            parameter_name        => $parameter_name,
            regexp                => $program_feature{$parameter_name}{regexp},
            sample_info_href      => \%sample_info,
        }
    );
    $return{$parameter_name}{regexp_return} = $version;
}

## Then version should be returned based on the parameter paths
is( $return{gatk_path}{regexp_return}, q{3.8},
    q{Added gatk version by regexp} );
is( $return{picardtools_path}{regexp_return},
    q{2.14.1}, q{Added picard version by regexp} );

## Given parameters to find version for, when invalid regexps
# Get the environment binary path
PARAMETER:
foreach my $parameter_name ( keys %program_feature ) {

    $active_parameter{$parameter_name} = get_dynamic_conda_path(
        {
            active_parameter_href => \%active_parameter,
            bin_file => $program_feature{$parameter_name}{program_name},
        }
    );
}

## Build cmd for system call to get version from program output
my $gatk_cmd =
    q{java -jar }
  . catfile( $active_parameter{gatk_path}, q{GenomeAnalysisTK.jar} )
  . $SPACE
  . q{--version 2>&1};
my $picard_cmd =
    q{java -jar }
  . catfile( $active_parameter{picardtools_path}, q{picard.jar} )
  . $SPACE
  . q{CreateSequenceDictionary --version 2>&1};
my $sambamba_cmd =
  q?sambamba 2>&1 | perl -nae 'if($_=~/sambamba\s(\S+)/) {print $1;last;}'?;

## Set in program features hash
$program_feature{gatk_path}{cmd}        = $gatk_cmd;
$program_feature{picardtools_path}{cmd} = $picard_cmd;
$program_feature{sambamba_depth}{cmd}   = $sambamba_cmd;

## Scramble the regexps, that need scrambling
$program_feature{gatk_path}{regexp}        = q{Not valid};
$program_feature{picardtools_path}{regexp} = q{Not valid};

## Activate parameter
$active_parameter{sambamba_depth} = 2;

PARAMETER:
foreach my $parameter_name ( keys %program_feature ) {

    my $version = get_program_version(
        {
            active_parameter_href => \%active_parameter,
            cmd                   => $program_feature{$parameter_name}{cmd},
            parameter_name        => $parameter_name,
            regexp                => $program_feature{$parameter_name}{regexp},
            sample_info_href      => \%sample_info,
        }
    );
    $return{$parameter_name}{cmd_return} = $version;
}

## Then version should be true
ok( $return{gatk_path}{cmd_return},        q{Added gatk version by cmd} );
ok( $return{picardtools_path}{cmd_return}, q{Added picard version by cmd} );
is( $return{sambamba_depth}{cmd_return}, undef, q{Skipped sambamba version} );

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
    -h/--help     Display this help message
    -v/--version  Display version
END_USAGE
}
