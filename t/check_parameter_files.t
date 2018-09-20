#!/usr/bin/env perl

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use Params::Check qw{ check allow last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };
use 5.018;

## CPANM
use autodie;
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};
Readonly my $COMMA   => q{,};

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
    my @modules = (qw{ MIP::Check::Path });

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Path qw{ check_parameter_files };

diag(   q{Test check_parameter_files from Path.pm v}
      . $MIP::Check::Path::VERSION
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

my @order_parameters =
  qw{ gatk_baserecalibration_known_sites gatk_genotypegvcfs_ref_gvcf human_genome_reference snpsift_annotation_files sv_vcfparser_select_file vcfparser_select_file };

my %active_parameter = (

    # To test array parameter
    gatk_baserecalibration_known_sites =>
      [ catfile( $Bin, qw{ data references GRCh37_dbsnp_-138-.vcf} ), ],
    gatk_genotypegvcfs_ref_gvcf => q{test_file},

    # To test scalar parameter
    human_genome_reference =>
      catfile( $Bin, qw{data references GRCh37_homo_sapiens_-d5-.fasta} ),
    mip                    => 1,
    gatk_baserecalibration => 1,
    gatk_genotypegvcfs     => 1,
    snpeff                 => 1,
    sv_vcfparser           => 0,

    # To test hash parameter
    snpsift_annotation_files => {
        catfile( $Bin,
            qw{data references GRCh37_anon-swegen_snp_-1000samples-.vcf.gz} )
          => q{AF},
    },
    sv_vcfparser_select_file => q{test_file},
);

my %parameter = load_yaml(
    {
        yaml_file => catfile(
            dirname($Bin), qw{ definitions rare_disease_parameters.yaml}
        ),
    }
);

$parameter{dynamic_parameter}{consensus_analysis_type} = q{wgs};

PARAMETER:
foreach my $parameter_name (@order_parameters) {

    check_parameter_files(
        {
            active_parameter_href => \%active_parameter,
            associated_programs_ref =>
              \@{ $parameter{$parameter_name}{associated_program} },
            log                    => $log,
            parameter_exists_check => $parameter{$parameter_name}{exists_check},
            parameter_href         => \%parameter,
            parameter_name         => $parameter_name,
        }
    );
}

is( $active_parameter{gatk_genotypegvcfs_ref_gvcf},
    q{test_file}, q{Returned for not required exome mode parameter} );

is( $active_parameter{sv_vcfparser_select_file},
    q{test_file}, q{Returned for not associated program} );

is( $active_parameter{vcfparser_select_file},
    undef, q{Returned for not active parameter} );

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
