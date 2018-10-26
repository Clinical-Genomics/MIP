#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
use MIP::Log::MIP_log4perl qw{ initiate_logger };
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
        q{MIP::Log::MIP_log4perl} => [qw{ initiate_logger }],
        q{MIP::Script::Utils}     => [qw{ help }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Parse::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Parse::Parameter qw{ parse_prioritize_variant_callers };

diag(   q{Test parse_prioritize_variant_callers from Parameter.pm v}
      . $MIP::Parse::Parameter::VERSION
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

## Given no structural active callers, when priority string is ok
my %active_parameter = (
    gatk_combinevariants_prioritize_caller => q{gatk,bcftools,freebayes},
    bcftools_mpileup                       => 1,
    freebayes                              => 1,
    gatk_variantrecalibration              => 1,
);

my %parameter = (
    dynamic_parameter => {
        variant_callers =>
          [qw{ freebayes bcftools_mpileup gatk_variantrecalibration}],
        structural_variant_callers => [qw{ pdelly }],
    },
    bcftools_mpileup          => { outdir_name => q{bcftools}, },
    freebayes                 => { outdir_name => q{freebayes}, },
    gatk_variantrecalibration => { outdir_name => q{gatk}, },
);

my $return = parse_prioritize_variant_callers(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_href        => \%parameter,
    }
);

## Then return undef
is( $return, undef, q{No active structural variant callers} );

## Given structural active callers, when priority string is ok
$active_parameter{sv_svdb_merge_prioritize} = q{delly};
$active_parameter{pdelly}                   = 1;
$parameter{pdelly}{outdir_name}             = q{delly};

my $is_ok = parse_prioritize_variant_callers(
    {
        active_parameter_href => \%active_parameter,
        log                   => $log,
        parameter_href        => \%parameter,
    }
);

ok( $is_ok, q{Passed parsing} );

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
