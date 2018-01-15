#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
Readonly my $COMMA      => q{,};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

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
        q{MIP::File::Format::Yaml} => [qw{load_yaml}],
        q{MIP::Get::Parameter}     => [qw{ get_capture_kit }],
        q{MIP::Log::MIP_log4perl}  => [qw{ initiate_logger }],
        q{MIP::Script::Utils}      => [qw{ help }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Set::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Get::Parameter qw{ get_capture_kit };
use MIP::Set::Parameter qw{ set_custom_default_to_active_parameter };

diag(   q{Test set_custom_default_to_active_parameter from Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $test_dir = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);

my %active_parameter = (
    cluster_constant_path  => catfile(qw{constant path}),
    family_id              => 1,
    human_genome_reference => q{human_genom_reference.fasta},
    outdata_dir            => catfile(qw{ a outdata dir }),
    sample_ids             => [qw{ sample_1 }],
);

my %parameter = load_yaml(
    {
        yaml_file =>
          catfile( dirname($Bin), qw{ definitions define_parameters.yaml } ),
    }
);

my @custom_default_parameters =
  qw{ analysis_type bwa_build_reference exome_target_bed sample_info_file };

foreach my $parameter_name (@custom_default_parameters) {

    set_custom_default_to_active_parameter(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
            parameter_name        => $parameter_name,
        }
    );
}

is( $active_parameter{analysis_type}{sample_1},
    q{wgs}, q{Set analysis_type default} );

is(
    $active_parameter{human_genome_reference},
    q{human_genom_reference.fasta},
    q{Set human_genome_reference default for bwa}
);

foreach my $capture_kit ( keys %{ $active_parameter{exome_target_bed} } ) {

    is(
        $capture_kit,
        $parameter{supported_capture_kit}{default}{latest},
        q{Set default capture kit}
    );

    is( $active_parameter{exome_target_bed}{$capture_kit},
        q{sample_1}, q{Set default capture kit for sample_1} );

}

my $sample_info_file = catfile(
    $active_parameter{outdata_dir},
    $active_parameter{family_id},
    $active_parameter{family_id} . $UNDERSCORE . q{qc_sample_info.yaml}
);

is( $parameter{sample_info_file}{default},
    $sample_info_file, q{Set default sample_info_file} );

is( $parameter{qccollect_sampleinfo_file}{default},
    $sample_info_file, q{Set default qccollect sample_info_file} );

## Reset to test infile_dirs whose default is dependent on analysis_type
set_custom_default_to_active_parameter(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
        parameter_name        => q{infile_dirs},
    }
);

my $path = catfile(
    $active_parameter{cluster_constant_path},
    $active_parameter{family_id},
    $active_parameter{analysis_type}{sample_1},
    q{sample_1}, q{fastq}
);
is( $active_parameter{infile_dirs}{$path},
    q{sample_1}, q{Set default infile_dirs} );

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
