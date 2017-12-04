#!/usr/bin/env perl

#### Copyright 2017 Henrik Stranneheim

use Modern::Perl qw{2014};
use warnings qw{FATAL utf8};
use autodie;
use 5.018;    #Require at least perl 5.18
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
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $VERSION,
          $NEWLINE;
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

  PERL_MODULES:
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

use MIP::Processmanagement::Processes qw{create_job_id_string_for_sample_id};

diag(
"Test create_job_id_string_for_sample_id $MIP::Processmanagement::Processes::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

## Base arguments
my $family_id           = q{family1};
my $sample_id           = q{sample1};
my $path                = q{MAIN};
my $family_id_chain_key = $family_id . $UNDERSCORE . $path;
my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

my %job_id = (
    $family_id_chain_key => {
        q{sample1} . $UNDERSCORE . $path => [qw{job_id_1 job_id_2}],
        q{sample2} . $UNDERSCORE . $path => [qw{job_id_3}],
        q{sample3} . $UNDERSCORE . $path => [qw{job_id_4 job_id_5 job_id_8}],
        q{sample4} . $UNDERSCORE . $path => [undef],
        $family_id_chain_key             => [qw{job_id_6}],
    },
);

### Creation of job id string using job ids from MAIN chain for sample id

## Add job_ids from MAIN chain to job_id_string
my $job_ids_string = create_job_id_string_for_sample_id(
    {
        job_id_href         => \%job_id,
        family_id           => $family_id,
        sample_id           => $sample_id,
        family_id_chain_key => $family_id_chain_key,
        sample_id_chain_key => $sample_id_chain_key,
        path                => $path,
    }
);

my $expected_job_id_string = q{:job_id_1:job_id_2};
is( $job_ids_string, $expected_job_id_string,
    q{Added job_id from MAIN job_id chain} );

### Creation of job id string using job ids from other chain with no previous job ids for sample id
my $path_other                = q{other};
my $family_id_chain_key_other = $family_id . $UNDERSCORE . $path_other;
my $sample_id_chain_key_other = $sample_id . $UNDERSCORE . $path_other;

## Add job_ids from MAIN chain to job_id_string
$job_ids_string = create_job_id_string_for_sample_id(
    {
        job_id_href         => \%job_id,
        family_id           => $family_id,
        sample_id           => $sample_id,
        family_id_chain_key => $family_id_chain_key_other,
        sample_id_chain_key => $sample_id_chain_key_other,
        path                => $path_other,
    }
);

my $expected_job_id_string_empty_other = q{:job_id_1:job_id_2};
is(
    $job_ids_string,
    $expected_job_id_string_empty_other,
    q{Added job_id from MAIN job_id chain initializing other chain}
);

### Inherit job_ids from other chain

%job_id = (
    $family_id_chain_key_other => {
        q{sample1} . $UNDERSCORE . $path_other => [qw{job_id_9 job_id_10}],
    },
);

## Creation of job id string using job ids from other chain for sample id
$job_ids_string = create_job_id_string_for_sample_id(
    {
        job_id_href         => \%job_id,
        family_id           => $family_id,
        sample_id           => $sample_id,
        family_id_chain_key => $family_id_chain_key_other,
        sample_id_chain_key => $sample_id_chain_key_other,
        path                => $path_other,
    }
);

my $expected_job_id_string_other = q{:job_id_9:job_id_10};

is( $job_ids_string, $expected_job_id_string_other,
    q{Added job_id from other job_id chain} );

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
