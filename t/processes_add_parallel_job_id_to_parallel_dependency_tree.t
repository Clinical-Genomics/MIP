#!/usr/bin/env perl

#### Copyright 2017 Henrik Stranneheim

use Modern::Perl qw{ 2018 };
use warnings qw{FATAL utf8};
use autodie;
use 5.026;    #Require at least perl 5.18
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{-no_match_vars};
use Params::Check qw{check allow last_error};

use FindBin qw{$Bin};    #Find directory of script
use File::Basename qw{dirname basename};
use File::Spec::Functions qw{catdir};
use Getopt::Long;
use Test::More;
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Script::Utils qw{help};

our $USAGE = build_usage( {} );

##Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

my $VERBOSE = 1;
our $VERSION = q{1.0.0};

###User Options
GetOptions(
    'h|help' => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },    #Display help text
    'v|version' => sub {
        done_testing();
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION, $NEWLINE;
        exit;
    },    #Display version number
    'vb|verbose' => $VERBOSE,
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
    my %perl_module;

    $perl_module{'MIP::Script::Utils'} = [qw{help}];

  MODULES:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load } . $module;
    }

    ## Modules
    my @modules = ('MIP::Processmanagement::Processes');

  MODULES:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load } . $module;
    }
}

use MIP::Processmanagement::Processes qw{add_parallel_job_id_to_parallel_dependency_tree};

diag(
"Test add_parallel_job_id_to_parallel_dependency_tree $MIP::Processmanagement::Processes::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

## Base arguments
my $case_id               = q{case1};
my $sample_id             = q{sample2};
my $path                  = q{MAIN};
my $case_id_chain_key     = $case_id . $UNDERSCORE . $path;
my $sample_id_chain_key   = $sample_id . $UNDERSCORE . $path;
my $sbatch_script_tracker = 0;
my $case_id_parallel_chain_key =
  $case_id . $UNDERSCORE . q{parallel} . $UNDERSCORE . $path . $sbatch_script_tracker;

my %job_id = (
    $case_id_chain_key => {
        $sample_id_chain_key        => [qw{job_id_1 job_id_2}],
        q{sample2_MAIN}             => [qw{job_id_3}],
        q{sample3_MAIN}             => [qw{job_id_4 job_id_5}],
        $case_id_parallel_chain_key => [qw{job_id_10 job_id_11}],
        $case_id_chain_key          => [qw{job_id_6}],
    },
);

### Parallel jobs

add_parallel_job_id_to_parallel_dependency_tree(
    {
        job_id_href           => \%job_id,
        case_id_chain_key     => $case_id_chain_key . q{no_parallel},
        id                    => $case_id,
        path                  => $path,
        sbatch_script_tracker => $sbatch_script_tracker,
        job_id_returned       => q{job_id_12},
    }
);

my $no_push_result = join $SPACE,
  @{ $job_id{$case_id_chain_key}{$case_id_parallel_chain_key} };
is(
    $no_push_result,
    q{job_id_10 job_id_11},
    q{No pushed parallel job_ids to parallel id}
);

add_parallel_job_id_to_parallel_dependency_tree(
    {
        job_id_href           => \%job_id,
        case_id_chain_key     => $case_id_chain_key,
        id                    => $case_id,
        path                  => $path,
        sbatch_script_tracker => $sbatch_script_tracker,
        job_id_returned       => q{job_id_12},
    }
);

my $push_result = join $SPACE,
  @{ $job_id{$case_id_chain_key}{$case_id_parallel_chain_key} };
is(
    $push_result,
    q{job_id_10 job_id_11 job_id_12},
    q{Pushed parallel job_ids to parallel id}
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
