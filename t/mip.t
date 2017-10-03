#!/usr/bin/env perl

###Copyright 2016 Henrik Stranneheim

use Modern::Perl qw(2014);    #CPAN
use warnings qw( FATAL utf8 );
use autodie qw(open close :all);    #CPAN
use 5.018;                          #Require at least perl 5.18
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);

use Test::More;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);
use File::Path qw(remove_tree);
use File::Spec::Functions qw(catfile catdir devnull);
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case
use IPC::Cmd qw[can_run run];
use Getopt::Long;
use Cwd;

##MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Script::Utils qw(help);

our $USAGE = build_usage( {} );

my $verbose = 1;
my $VERSION = '0.0.1';

my $config_file = catfile( dirname($Bin), qw(templates mip_config.yaml) );

###User Options
GetOptions(
    'c|config_file:s' => \$config_file,
    'h|help'    => sub { print STDOUT $USAGE, "\n"; exit; },  #Display help text
    'v|version' => sub {
        print STDOUT "\n" . basename($PROGRAM_NAME) . q{ } . $VERSION, "\n\n";
        exit;
    },    #Display version number
    'vb|verbose' => $verbose,
  )
  or help(
    {
        USAGE     => $USAGE,
        exit_code => 1,
    }
  );

ok( can_run('mip'), 'Checking can run mip' );

## Test execution of mip
# Create array ref for cmd
my $cmds_ref = [
    'mip',
    qw(-f 643594-miptest),
    '-c',
    $config_file,
    '-ifd',
    catfile(qw(data 643594-miptest test_data ADM1059A1 fastq=ADM1059A1)),
    '-ifd',
    catfile(qw(data 643594-miptest test_data ADM1059A2 fastq=ADM1059A2)),
    '-ifd',
    catfile(qw(data 643594-miptest test_data ADM1059A3 fastq=ADM1059A3)),
    qw(-rio 1),
    qw(-dra 2),
    qw(-pvep 0),
    qw(-psvv 0),
];

my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
  run( command => $cmds_ref, verbose => $verbose );
ok( $success, 'Executed mip' );

my $qc_pedigree = catfile( getcwd(),
    qw(data 643594-miptest analysis 643594-miptest qc_pedigree.yaml) );
ok( -f $qc_pedigree, 'Checking for qc_pedigree.yaml' );

##Clean-up
remove_tree( catfile( getcwd(), qw(data 643594-miptest analysis) ) );

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
    -c/--config_file YAML config file for analysis parameters (defaults to ../templates/mip_config.yaml")
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
