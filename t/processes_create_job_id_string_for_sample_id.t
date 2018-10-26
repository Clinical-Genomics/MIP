#!/usr/bin/env perl

#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
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
use Test::More;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

## Constants
Readonly my $COMMA      => q{,};
Readonly my $EMPTY_STR  => q{};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

my $VERBOSE = 1;
our $VERSION = q{1.0.1};

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
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION . $NEWLINE;
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
    my @modules = (q{MIP::Processmanagement::Processes});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Processmanagement::Processes qw{ create_job_id_string_for_sample_id };

diag(   q{Test create_job_id_string_for_sample_id from Processes.pm v}
      . $MIP::Processmanagement::Processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $case_id               = q{case1};
my $sample_id             = q{sample1};
my $path                  = q{MAIN};
my $sbatch_script_tracker = 0;
my $case_id_chain_key     = $case_id . $UNDERSCORE . $path;
my $sample_id_chain_key   = $sample_id . $UNDERSCORE . $path;
my $sample_id_parallel_chain_key =
  $sample_id . $UNDERSCORE . q{parallel} . $UNDERSCORE . $path . $sbatch_script_tracker;

my %job_id = (
    $case_id_chain_key => {
        q{sample1} . $UNDERSCORE . $path => [qw{job_id_1 job_id_2}],
        q{sample2} . $UNDERSCORE . $path => [qw{job_id_3}],
        q{sample3} . $UNDERSCORE . $path => [qw{job_id_4 job_id_5 job_id_8}],
        q{sample4} . $UNDERSCORE . $path => [undef],
        $case_id_chain_key               => [qw{job_id_6}],
    },
);

## Given job ids from MAIN chain, when using sample id

## Add job_ids from MAIN chain to job_id_string
my $job_ids_string = create_job_id_string_for_sample_id(
    {
        job_id_href         => \%job_id,
        case_id             => $case_id,
        sample_id           => $sample_id,
        case_id_chain_key   => $case_id_chain_key,
        sample_id_chain_key => $sample_id_chain_key,
        path                => $path,
    }
);

## Then add job_ids for sample1 from MAIN
my $expected_job_id_string = q{:job_id_1:job_id_2};
is( $job_ids_string, $expected_job_id_string, q{Added job_id from MAIN job_id chain} );

## Given job id string using job ids from other chain with no previous job ids for sample id
my $path_other                = q{other};
my $case_id_chain_key_other   = $case_id . $UNDERSCORE . $path_other;
my $sample_id_chain_key_other = $sample_id . $UNDERSCORE . $path_other;

## Add job_ids from MAIN chain to job_id_string
$job_ids_string = create_job_id_string_for_sample_id(
    {
        job_id_href         => \%job_id,
        case_id             => $case_id,
        sample_id           => $sample_id,
        case_id_chain_key   => $case_id_chain_key_other,
        sample_id_chain_key => $sample_id_chain_key_other,
        path                => $path_other,
    }
);

## Then add job_ids for sample1 inherited from MAIN chain job ids
my $expected_job_id_string_empty_other = q{:job_id_1:job_id_2};
is(
    $job_ids_string,
    $expected_job_id_string_empty_other,
    q{Added job_id from MAIN job_id chain initializing other chain}
);

### Inherit job_ids from other chain

%job_id = (
    $case_id_chain_key_other => {
        q{sample1} . $UNDERSCORE . $path_other => [qw{job_id_9 job_id_10}],
    },
);

## Given job id string using job ids from other chain for sample id
$job_ids_string = create_job_id_string_for_sample_id(
    {
        job_id_href         => \%job_id,
        case_id             => $case_id,
        sample_id           => $sample_id,
        case_id_chain_key   => $case_id_chain_key_other,
        sample_id_chain_key => $sample_id_chain_key_other,
        path                => $path_other,
    }
);

my $expected_job_id_string_other = q{:job_id_9:job_id_10};

## Then add job_ids for other chain
is( $job_ids_string, $expected_job_id_string_other,
    q{Added job_id from other job_id chain} );

## Given job id string using parallel job ids from other chain with no previous job ids for sample id

## Clean-up for new test
%job_id = ();

$job_id{$case_id_chain_key}{$sample_id_parallel_chain_key} =
  [qw{ job_id_11 job_id_12 }];

## Add job_ids from MAIN chain to job_id_string
$job_ids_string = create_job_id_string_for_sample_id(
    {
        job_id_href           => \%job_id,
        case_id               => $case_id,
        sample_id             => $sample_id,
        case_id_chain_key     => $case_id_chain_key_other,
        sample_id_chain_key   => $sample_id_chain_key_other,
        path                  => $path_other,
        sbatch_script_tracker => $sbatch_script_tracker,
    }
);

## Then add job_ids for sample1 inherited from MAIN chain parallel job ids
my $expected_job_id_string_parallel_other = q{:job_id_11:job_id_12};
is(
    $job_ids_string,
    $expected_job_id_string_parallel_other,
    q{Added parallel job_id from MAIN job_id chain initializing other chain}
);

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

##build_usage

##Function : Build the USAGE instructions
##Returns  : ""
##Arguments: $program_name
##         : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$program_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
