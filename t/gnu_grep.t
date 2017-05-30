#!/usr/bin/env perl

###Copyright 2017 Henrik Stranneheim

use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie;
use 5.018;    #Require at least perl 5.18
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);
use Params::Check qw[check allow last_error];

use FindBin qw( $Bin );    #Find directory of script
use File::Basename qw( dirname basename );
use File::Spec::Functions qw( catdir catfile devnull );
use Getopt::Long;
use Test::More;

##MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use Script::Utils qw( help );

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '0.0.0';

###User Options
GetOptions(
    'h|help' => sub { done_testing(); print {*STDOUT} $USAGE, "\n"; exit; }
    ,    #Display help text
    'v|version' => sub {
        done_testing();
        print {*STDOUT} "\n" . basename($PROGRAM_NAME) . q{  } . $VERSION,
          "\n\n";
        exit;
    },    #Display version number
    'vb|verbose' => $VERBOSE,
  )
  or (
    done_testing(),
    Script::Utils::help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
##Modules with import
    my %perl_module;

    $perl_module{'Script::Utils'} = [qw(help)];

    while ( my ( $module, $module_import ) = each %perl_module ) {

        use_ok( $module, @{$module_import} ) or BAIL_OUT "Can't load $module";
    }

##Modules
    my @modules = ('MIP::Gnu::Software::Gnu_grep');

    for my $module (@modules) {

        require_ok($module) or BAIL_OUT "Can't load $module";
    }
}

use MIP::Gnu::Software::Gnu_grep;

diag(
"Test Gnu_grep $MIP::Gnu::Software::Gnu_grep::VERSION, Perl $^V, $EXECUTABLE_NAME"
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

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
