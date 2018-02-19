#!/usr/bin/env perl

use 5.018;
use autodie;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use Modern::Perl qw{ 2014 };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = q{1.0.1};

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
        q{MIP::Log::MIP_log4perl} => [qw{ initiate_logger }],
        q{MIP::Script::Utils}     => [qw{ help }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Update::Programs});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Update::Programs qw{ update_program_mode_for_analysis_type };

diag(   q{Test update_program_mode_for_analysis_type from Programs.pm v}
      . $MIP::Update::Programs::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $test_dir = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Create log object
my $log = initiate_logger(
    {
        file_path      => $test_log_path,
        log_name       => q{TEST},
        categories_ref => [qw{ TRACE LogFile }],
    }
);

my @programs =
  qw{ cnvnator delly_call delly_reformat samtools_subsample_mt tiddit };
my %active_parameter = (
    pcnvnator              => 1,
    pdelly_call            => 1,
    pdelly_reformat        => 1,
    pmanta                 => 1,
    psamtools_subsample_mt => 1,
    ptiddit                => 1,
);

my @warning_msgs = update_program_mode_for_analysis_type(
    {
        active_parameter_href   => \%active_parameter,
        consensus_analysis_type => q{wgs},
        log                     => $log,
        programs_ref            => \@programs,
    }
);

is( @warning_msgs, 0, q{No updates to programs mode} );

@warning_msgs = update_program_mode_for_analysis_type(
    {
        active_parameter_href   => \%active_parameter,
        consensus_analysis_type => q{wes},
        log                     => $log,
        programs_ref            => \@programs,
    }
);
## Alias
my $cnvnator_mode              = $active_parameter{pcnvnator};
my $delly_call_mode            = $active_parameter{pdelly_call};
my $delly_reformat_mode        = $active_parameter{pdelly_reformat};
my $manta_mode                 = $active_parameter{pmanta};
my $samtools_subsample_mt_mode = $active_parameter{psamtools_subsample_mt};
my $tiddit_mode                = $active_parameter{ptiddit};

## Test program mode updates and warnings
is( $cnvnator_mode,       0, q{Updated programs mode for cnvnator} );
is( $delly_call_mode,     0, q{Updated programs mode for delly_call} );
is( $delly_reformat_mode, 0, q{Updated programs mode for delly_reformat} );
is( $manta_mode,          1, q{Updated programs mode for manta} );
is( $samtools_subsample_mt_mode, 0,
    q{Updated programs mode for samtools_subsample_mt} );
is( $tiddit_mode, 0, q{Updated programs mode for tiddit} );
isnt( @warning_msgs, 0, q{Generated warning message} );

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
            strict_type => 1,
            store       => \$program_name,
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
