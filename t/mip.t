#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw( :full :short );
use English qw(-no_match_vars);
use File::Basename qw(dirname basename);
use File::Path qw(remove_tree);
use File::Spec::Functions qw(catfile catdir devnull);
use FindBin qw($Bin);
use Getopt::Long;
use open qw( :encoding(UTF-8) :std );
use Params::Check qw[check allow last_error];
use Test::More;
use utf8;
use warnings qw( FATAL utf8 );

## Cpanm
use autodie qw(open close :all);
use IPC::Cmd qw[can_run run];
use Modern::Perl qw(2014);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Script::Utils qw(help);

our $USAGE = build_usage( {} );

my $verbose = 1;
our $VERSION = '0.0.2';

## Set default
my $config_file = catfile( dirname($Bin), qw(templates mip_config.yaml) );
my $cluster_constant_path = catfile( dirname($Bin), qw{ t data} );

### User Options
GetOptions(
    q{c|config_file:s} => \$config_file,
    ## Display help text
    q{h|help} => sub { print {*STDOUT} $USAGE, "\n"; exit; },
    ## Display version number
    q{v|version} => sub {
        print {*STDOUT} "\n" . basename($PROGRAM_NAME) . q{ } . $VERSION,
          "\n\n";
        exit;
    },
    q{vb|verbose} => $verbose,
  )
  or help(
    {
        USAGE     => $USAGE,
        exit_code => 1,
    }
  );

ok( can_run(q{mip}), q{Checking can run mip} );

## Test execution of mip
# Create array ref for cmd
my $cmds_ref = [
    q{mip},
    qw(-f 643594-miptest),
    q{-c},
    $config_file,
    q{-ccp},
    $cluster_constant_path,
    q{-ifd},
    catfile(
        $cluster_constant_path,
        qw( 643594-miptest test_data ADM1059A1 fastq=ADM1059A1)
    ),
    q{-ifd},
    catfile(
        $cluster_constant_path,
        qw( 643594-miptest test_data ADM1059A2 fastq=ADM1059A2)
    ),
    q{-ifd},
    catfile(
        $cluster_constant_path,
        qw( 643594-miptest test_data ADM1059A3 fastq=ADM1059A3)
    ),
    qw(--rio),
    qw(--dra),
    qw(--psvv 0),
    qw(--pfreebayes 0),
    qw(--ptiddit 0),
	qw(--sv_svdb_merge_prioritize manta,delly,cnvnator),
    qw(--pvcf2cytosure 0),
	qw(--gatk_combinevariants_prioritize_caller gatk,bcftools),
];

my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
  run( command => $cmds_ref, verbose => $verbose );
ok( $success, q{Executed mip} );

my $qc_sample_info_file = catfile( $cluster_constant_path,
    qw{ 643594-miptest analysis 643594-miptest_qc_sample_info.yaml} );
ok( -f $qc_sample_info_file,
    q{Checking for 643594-miptest_qc_sample_info.yaml} );

## Clean-up
remove_tree( catfile( $cluster_constant_path, qw( 643594-miptest analysis) ) );

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

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -c/--config_file YAML config file for analysis parameters (defaults to ../templates/mip_config.yaml")
    -vb/--verbose    Verbose
    -h/--help        Display this help message
    -v/--version     Display version
END_USAGE
}
