#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use MIP::Log::MIP_log4perl qw{ initiate_logger };
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
    my @modules = (q{MIP::Check::Reference});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Reference qw{ check_references_for_vt };

diag(   q{Test check_references_for_vt from Reference.pm v}
      . $MIP::Check::Reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create a temp logger
my $test_dir = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

# Create log object
## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);

my %active_parameter_test = (
    pgatk_realigner                  => 1,
    pgatk_baserecalibration          => 1,
    pgatk_variantevalall             => 1,
    pgatk_variantevalexome           => 1,
    psnpeff                          => 1,
    gatk_realigner_indel_known_sites => [
        qw{ GRCh37_1000g_indels_-phase1-.vcf GRCh37_mills_and_1000g_indels_-gold_standard-.vcf },
    ],
    gatk_baserecalibration_known_sites => [
        qw{ GRCh37_dbsnp_-138-.vcf GRCh37_1000g_indels_-phase1-.vcf GRCh37_mills_and_1000g_indels_-gold_standard-.vcf },
    ],
    gatk_varianteval_dbsnp   => q{GRCh37_dbsnp_-138_esa_129-.vcf},
    snpsift_annotation_files => {
        SWEREF => q{GRCh37_anon-swegen_snp_-1000samples-.vcf.gz},
        EXAC   => q{GRCh37_exac_reheader_-r0.3.1-.vcf.gz},
    },
);

my %parameter_test = (
    gatk_realigner_indel_known_sites =>
      { associated_program => [qw{ pgatk_realigner }], data_type => q{ARRAY}, },
    gatk_baserecalibration_known_sites => {
        associated_program => [qw{ pgatk_baserecalibration }],
        data_type          => q{ARRAY},
    },
    gatk_varianteval_dbsnp => {
        associated_program =>
          [qw{ pgatk_variantevalall pgatk_variantevalexome }],
        data_type => q{SCALAR},
    },
    snpsift_annotation_files =>
      { associated_program => [qw{ psnpeff }], data_type => q{HASH}, },
);

my @vt_references_test =
  qw{ gatk_realigner_indel_known_sites gatk_baserecalibration_known_sites gatk_varianteval_dbsnp snpsift_annotation_files };

my @refs_to_process = check_references_for_vt(
    {
        active_parameter_href => \%active_parameter_test,
        log                   => Log::Log4perl->get_logger(q{TEST}),
        parameter_href        => \%parameter_test,
        vt_references_ref     => \@vt_references_test,
    }
);

is( @refs_to_process, 0, q{Test passed, 0 references to process.} );

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
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
