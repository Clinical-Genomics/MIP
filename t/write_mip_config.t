#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Path qw{ remove_tree };
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
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::File::Format::Yaml qw{ load_yaml };
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
    my %perl_module = (
        q{MIP::File::Format::Yaml} => [qw{ load_yaml }],
        q{MIP::Log::MIP_log4perl}  => [qw{ initiate_logger }],
        q{MIP::Script::Utils}      => [qw{ help }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::File::Format::Config});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::File::Format::Config qw{ write_mip_config };

diag(   q{Test write_mip_config from Config.pm v}
      . $MIP::File::Format::Config::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $test_dir      = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);

## Given a config file path
my %active_parameter = (
    config_file_analysis => catfile( $Bin, qw{ data test_data config_file } ),
    associated_program   => q{some_program},
    test_key             => q{This is a value},
);
my %sample_info;
trap {
    write_mip_config(
        {
            active_parameter_href => \%active_parameter,
            log                   => $log,
            remove_keys_ref       => [qw{ associated_program }],
            sample_info_href      => \%sample_info,
        }
    )
};

## Then we should broadcast an INFO message
like( $trap->stderr, qr/INFO/xms, q{Send info log message} );

## Then some keys should be removed
is( $active_parameter{associated_program}, undef, q{Removed keys from active parameter} );

## Then the file should have been writen to disc and be reloaded
my %config_parameter =
  load_yaml( { yaml_file => $active_parameter{config_file_analysis}, } );

is_deeply( \%config_parameter, \%active_parameter, q{Wrote config file} );

## Then the path of the config should be transfered to sample_info
is(
    $sample_info{config_file_analysis},
    $active_parameter{config_file_analysis},
    q{Transfered config file analysis path to sample_info}
);

## Clean-up
remove_tree( $active_parameter{config_file_analysis} );

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
