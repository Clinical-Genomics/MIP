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

use MIP::Gnu::Software::Gnu_grep qw(gnu_grep);
use MIP::Test::Commands qw(generate_call);

diag(
"Test Gnu_grep $MIP::Gnu::Software::Gnu_grep::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

## Base parameters
my %base_parameters = (
    outfile_path => {
        input           => 'outfile.test',
        expected_output => '> outfile.test',
    },
    stderrfile_path => {
        input           => 'stderrfile.test',
        expected_output => '2> stderrfile.test',
    },
);

my %required_parameters = (
    infile_path => {
        input           => 'infile.test',
        expected_output => 'infile.test',
    },
);

## Specific parameters
my %specific_parameter = (
    invert_match => {
        input           => 1,
        expected_output => '--invert-match',
    },
    filter_file_path => {
        input           => 'test_file',
        expected_output => '--file=test_file',
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&gnu_grep;

## Test both base and function specific parameters
my @parameters = ( \%base_parameters, \%specific_parameter );

foreach my $parameter_href (@parameters) {

    my @commands = generate_call(
        {
            parameter_href           => $parameter_href,
            required_parameters_href => \%required_parameters,
            module_function_cref     => $module_function_cref,
            function_base_command    => 'grep',
        }
    );
}

## Test writing to file
my $FILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

## Store file in memory
my $variable;

open $FILEHANDLE, '>', \$variable or croak q{Can't open STDOUT: } . $OS_ERROR;
gnu_grep(
    {
        infile_path => 'infile.test',
        FILEHANDLE  => $FILEHANDLE,
    }
);
close $FILEHANDLE;

ok( $variable =~ /infile.test/x, 'Test write commands to file' );

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
