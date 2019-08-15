#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Getopt::Long;
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
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.1';

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

##Modules
    my @modules = ('MIP::Processmanagement::Processes');

  MODULES:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load } . $module;
    }
}

use MIP::Processmanagement::Processes
  qw{add_parallel_job_id_to_sample_id_dependency_tree};

diag(   q{Test add_parallel_job_id_to_sample_id_dependency_tree from Processes.pm v}
      . $MIP::Processmanagement::Processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $sample_id           = q{sample1};
my $path                = q{MAIN};
my $case_id_chain_key   = q{case1} . $UNDERSCORE . $path;
my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;
my $chain_key_parallel  = q{parallel};

my %infile_lane_prefix = (
    sample1 => [qw{1_lane1 1_lane2}],
    sample2 => [qw{2_lane1}],
    sample3 => [qw{3_lane4 3_lane5}],
);

my %job_id = (
    $case_id_chain_key => {
        $sample_id_chain_key => [qw{job_id_1 job_id_2}],
        q{sample2_MAIN}      => [qw{job_id_3}],
        q{sample3_MAIN}      => [qw{job_id_4 job_id_5}],
        $case_id_chain_key   => [qw{job_id_6}],
    },
);

### No previous parallel job_ids
add_parallel_job_id_to_sample_id_dependency_tree(
    {
        infile_lane_prefix_href => \%infile_lane_prefix,
        job_id_href             => \%job_id,
        case_id_chain_key       => $case_id_chain_key,
        sample_id_chain_key     => $sample_id_chain_key,
        path                    => $path,
        sample_id               => $sample_id,
    }
);
my $no_parallel_push_result = join $SPACE,
  @{ $job_id{$case_id_chain_key}{$sample_id_chain_key} };
is( $no_parallel_push_result, q{job_id_1 job_id_2}, q{No parallel job_id} );

### Previous parallel jobs
## Set-up previous parallel job
INFILES:
while ( my ($infile_index) = each @{ $infile_lane_prefix{$sample_id} } ) {

    # Set key
    my $chain_key_parallel_job =
        $sample_id
      . $UNDERSCORE
      . $chain_key_parallel
      . $UNDERSCORE
      . $path
      . $infile_index;

    push @{ $job_id{$case_id_chain_key}{$chain_key_parallel_job} },
      q{job_id_} . $infile_index;
}

add_parallel_job_id_to_sample_id_dependency_tree(
    {
        infile_lane_prefix_href => \%infile_lane_prefix,
        job_id_href             => \%job_id,
        case_id_chain_key       => $case_id_chain_key,
        sample_id_chain_key     => $sample_id_chain_key,
        path                    => $path,
        sample_id               => $sample_id,
    }
);

my $parallel_push_result = join $SPACE,
  @{ $job_id{$case_id_chain_key}{$sample_id_chain_key} };
is(
    $parallel_push_result,
    q{job_id_1 job_id_2 job_id_0 job_id_1},
    q{Push parallel job_id}
);

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

## Function : Build the USAGE instructions
## Returns  :
## Arguments: $program_name => Name of the script

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

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
