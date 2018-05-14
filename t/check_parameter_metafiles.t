#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use Getopt::Long;
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
    my %perl_module = (
        q{MIP::File::Format::Yaml} => [qw{ load_yaml }],
        q{MIP::Script::Utils}      => [qw{ help }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Check::Reference});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Reference qw{ check_parameter_metafiles };

diag(   q{Test check_parameter_metafiles from Reference.pm v}
      . $MIP::Check::Reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %parameter = load_yaml(
    {
        yaml_file => catfile(
            dirname($Bin), qw{ definitions rare_disease_parameters.yaml}
        ),
    }
);

## File info hash
my %file_info = (

    # BWA human genome reference file endings
    bwa_build_reference => [qw{ .bwt .ann .amb .pac .sa }],

    exome_target_bed =>
      [qw{ .infile_list .pad100.infile_list .pad100.interval_list }],

    # Human genome meta files
    human_genome_reference_file_endings => [qw{ .dict .fai }],

    # RTG human genome reference file endings
    rtg_vcfeval_reference_genome => [qw{ _sdf_dir }],
);

my $parameter_name = q{exome_target_bed};

## Given hash entries with no active parameter and existing files
my %active_parameter = (
    not_correct_key => {
        catfile( $Bin,
            qw{ data references GRCh37_agilent_sureselect_targets_cre_-v1-.bed }
        ) => q{sample1},
    },
);

check_parameter_metafiles(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
    }
);

## Then set build switch should not be set to zero i.e. default of "1" is kept
## as this will not be evaluated in sub or downstream
is( $parameter{$parameter_name}{build_file},
    1, q{No active parameter and existing file} );

## Given hash entries with active parameter, existing files, and no active associated programs
%active_parameter = (
    exome_target_bed => {
        catfile( $Bin,
            qw{ data references GRCh37_agilent_sureselect_targets_cre_-v1-.bed }
        ) => q{sample1},
    },
);

check_parameter_metafiles(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
    }
);

## Then set build switch should not be set to zero i.e. default of "1" is kept
## as this will not be evaluated in sub or downstream
is( $parameter{$parameter_name}{build_file},
    1,
    q{Active parameter and existing file, and no active associated program} );

## Given hash entries with active parameter, files exists, and active associated programs
%active_parameter = (
    exome_target_bed => {
        catfile( $Bin,
            qw{ data references GRCh37_agilent_sureselect_targets_cre_-v1-.bed }
        ) => q{sample1},
    },
    ppicardtools_collecthsmetrics => 1,
);

check_parameter_metafiles(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
    }
);

## Then set build switch should be set to zero to not build meta files as
## they are all ready and do not need to be built
is( $parameter{$parameter_name}{build_file},
    0, q{Active parameter, existing file and active associated program} );

## Given hash entries where one path exists and one does not
$active_parameter{exome_target_bed}
  { catfile( $Bin, qw{ data references does_not_exists.bed } ) } = q{sample1};

check_parameter_metafiles(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
    }
);

## Then set build file to true to rebuild for all parameter keys
is( $parameter{exome_target_bed}{build_file}, 1,
q{Set build file switch for hash parameter reference with mixed existence to 1}
);

## Given scalar entries with active parameter, files exists, and active associated programs
$active_parameter{bwa_build_reference} = 1;
$active_parameter{human_genome_reference} =
  catfile( $Bin, qw{ data references GRCh37_homo_sapiens_-d5-.fasta } );
$active_parameter{pbwa_mem} = 1;

check_parameter_metafiles(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
    }
);

## Then set build file to false
is( $parameter{bwa_build_reference}{build_file}, 0,
q{Active parameter, existing file and active associated program for scalar parameter reference}
);

## Given scalar entries with active parameter, files do not exists, and active associated programs
%file_info = (

    # BWA human genome reference file endings
    bwa_build_reference => [qw{ .does_not_exist .ann .amb .pac .sa }],
);

check_parameter_metafiles(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
    }
);

## Then set build file to true
is( $parameter{bwa_build_reference}{build_file}, 1,
q{Active parameter, no existing file and active associated program for scalar parameter reference}
);

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
