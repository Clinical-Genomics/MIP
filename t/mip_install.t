#!/usr/bin/env perl

use Carp;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

BEGIN {

    ## Special case to initiate testing
    my @modules = qw{ Test::More };

    ## Evaluate that all modules required are installed
    check_modules( \@modules );

    sub check_modules {

        ## Function : Evaluate that all modules required are installed
        ## Returns  :
        ## Arguments: $modules_ref
        ##          : $modules_ref => Array of module names

        my ($modules_ref) = @_;

      MODULE:
        foreach my $module ( @{$modules_ref} ) {

# Replace "::" with "/" since the automatic replacement magic only occurs for barewords.
            $module =~ s/::/\//gxms;

            # Add perl module ending for the same reason
            $module .= q{.pm};

            if ( eval { require $module; 1 } ) {
            }
            if ($EVAL_ERROR) {

                warn q{NOTE: }
                  . $module
                  . q{ not installed - Please install to run install.t.}
                  . qq{\n};
                warn q{NOTE: Aborting!} . qq{\n};
                exit 1;
            }
        }
        return;
    }

    ## Initate tests
    say {*STDOUT} q{Initiate tests:} . qq{\n};
    say {*STDOUT} q{Testing default perl modules and selected functions}
      . qq{\n};

    ## More proper testing
    use Test::More;

    ## Modules with import
    my %perl_module;

    $perl_module{charnames}                = [qw{ :full :short }];
    $perl_module{English}                  = [qw{ -no_match_vars }];
    $perl_module{q{File::Basename}}        = [qw{ dirname basename }];
    $perl_module{q{File::Spec::Functions}} = [qw{ catfile catdir devnull }];
    $perl_module{FindBin}                  = [qw{ $Bin }];
    $perl_module{q{IPC::Cmd}}              = [qw{ can_run run }];
    $perl_module{open}                     = [qw{ :encoding(UTF-8) :std }];
    $perl_module{q{Params::Check}}         = [qw{ check allow last_error }];
    $perl_module{warnings}                 = [qw{ FATAL utf8 }];

  PERL_MODULE:
    while ( my ( $module, $module_imports_ref ) = each %perl_module ) {

        use_ok( $module, @{$module_imports_ref} )
          or BAIL_OUT q{Can't load } . $module;
    }

    ## Modules
    @modules = qw{ Carp Cwd Getopt::Long utf8 strict warnings };

  MODULE:
    for my $module (@modules) {

        require_ok($module) or BAIL_OUT q{Can't load } . $module;
    }
}

use Cwd;
use FindBin qw{ $Bin };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catfile catdir devnull };
use Getopt::Long;
use IPC::Cmd qw{ can_run run };
use Params::Check qw{ check allow last_error };
use Test::More;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $conda_dir_path;
my $verbose = 1;
our $VERSION = q{1.0.1};

###User Options
GetOptions(
    q{cdp|conda_dir_path:s} => \$conda_dir_path,
    q{vb|verbose}           => $verbose,
    q{h|help}               => sub {
        done_testing();
        print {*STDOUT} $USAGE, "\n";
        exit;
    },    #Display help text
    q{v|version} => sub {
        done_testing();
        print {*STDOUT} "\n" . basename($PROGRAM_NAME) . q{  } . $VERSION,
          "\n\n";
        exit;
    },    #Display version number
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

ok( can_run(q{conda}), q{Checking can run of conda binary} );

ok(
    catfile( dirname($Bin), q{mip_install.pl} ),
    q{Locating install script in MIP dir}
);

my $install_script = catfile( dirname($Bin), q{mip_install.pl} );

## Test execution of install.pl
# Create array ref for cmd
my $cmds_ref = [ qw{ perl }, $install_script, qw{ -sp mip_scripts } ];
if ($conda_dir_path) {

    push @{$cmds_ref}, q{--conda_dir_path}, $conda_dir_path;
}

my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) = run(
    command => $cmds_ref,
    verbose => $verbose
);

ok( $success, q{Executed mip_install.pl} );

is( -e catfile( getcwd(), q{mip.sh} ),
    1, q{Locating created mip.sh in MIP dir} );

## Clean-up mip_install.pl output
my @outfiles = qw{ mip.sh mip_install_*.log };

foreach my $outfile (@outfiles) {

    my $cmd = q{rm } . catfile( getcwd(), $outfile );
    if (
        scalar run(
            command => $cmd,
            timeout => 20,
            verbose => $verbose,
        )
      )
    {
        say {STDOUT} q{Removed outfile successfully: } . $outfile;
    }
}

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

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    return <<"END_USAGE";
 $program_name [options]
    -cdp/--conda_dir_path The conda directory path (Default: "HOME/miniconda")
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
