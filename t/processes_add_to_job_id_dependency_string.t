#!/usr/bin/env perl

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
Readonly my $EMPTY_STR  => q{};
Readonly my $UNDERSCORE => q{_};

my $VERBOSE = 1;

###User Options
GetOptions(
    'h|help' => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },    #Display help text
    'v|version' => sub {
        done_testing();
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $NEWLINE;
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

use MIP::Processmanagement::Processes qw{add_to_job_id_dependency_string};

diag(
"Test add_to_job_id_dependency_string $MIP::Processmanagement::Processes::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

## Base arguments
my $sample_id         = q{sample2};
my $path              = q{MAIN};
my $case_id_chain_key = q{case1} . $UNDERSCORE . $path;

my %job_id = (
    $case_id_chain_key => {
        q{sample1} . $UNDERSCORE . $path => [qw{job_id_1 job_id_2}],
        q{sample2} . $UNDERSCORE . $path => [qw{job_id_3}],
        q{sample3} . $UNDERSCORE . $path => [qw{job_id_4 job_id_5 job_id_8}],
        q{sample4} . $UNDERSCORE . $path => [undef],
        $case_id_chain_key               => [qw{job_id_6}],
    },
);

### Sample job

## Add 1 job_id to job_id_string
my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

my $job_ids_string = add_to_job_id_dependency_string(
    {
        job_id_href       => \%job_id,
        case_id_chain_key => $case_id_chain_key,
        chain_key         => $sample_id_chain_key,
    }
);
my $expected_job_id_string = q{:job_id_3};
is( $job_ids_string, $expected_job_id_string, q{Added 1 job_id to job_id_string} );

## Add 2 job_ids to job_id_string
$sample_id           = q{sample1};
$sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

## Add to job_id string
$job_ids_string = add_to_job_id_dependency_string(
    {
        job_id_href       => \%job_id,
        case_id_chain_key => $case_id_chain_key,
        chain_key         => $sample_id_chain_key,
    }
);

$expected_job_id_string = q{:job_id_1:job_id_2};

is( $job_ids_string, $expected_job_id_string, q{Added 2 job_ids to job_id_string} );

## Add 3 job_ids to job_id_string
$sample_id           = q{sample3};
$sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

## Add to job_id string
$job_ids_string = add_to_job_id_dependency_string(
    {
        job_id_href       => \%job_id,
        case_id_chain_key => $case_id_chain_key,
        chain_key         => $sample_id_chain_key,
    }
);

$expected_job_id_string = q{:job_id_4:job_id_5:job_id_8};

is( $job_ids_string, $expected_job_id_string, q{Added 3 job_ids to job_id_string} );

## Do not add undef job_ids to job_id_string
$sample_id           = q{sample4};
$sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

## Add to job_id string
$job_ids_string = add_to_job_id_dependency_string(
    {
        job_id_href       => \%job_id,
        case_id_chain_key => $case_id_chain_key,
        chain_key         => $sample_id_chain_key,
    }
);

$expected_job_id_string = $EMPTY_STR;

is( $job_ids_string, $expected_job_id_string, q{Nothing was added to job_id_string} );

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
