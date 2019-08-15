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
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

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
        q{MIP::Log::MIP_log4perl} => [qw{ initiate_logger }],
        q{MIP::Script::Utils}     => [qw{ help }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Check::Parameter});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Parameter qw{ check_snpsift_keys };

diag(   q{Test check_snpsift_keys from Parameter.pm v}
      . $MIP::Check::Parameter::VERSION
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

## Given matching snpsift annotation and outinfo keys
my %active_parameter = (
    snpsift_annotation_files => {
        q{grch37_clinvar_reformat_-2017-10-29-.vcf.gz} => q{CLNSIG,CLNVID,CLNREVSTAT},
        q{grch37_gnomad.genomes_-r2.0.1-.vcf.gz}       => q{AF,AF_POPMAX},
        q{grch37_anon-swegen_str_nsphs_-1000samples-.vcf.gz} =>
          q{AF,AC_Hom,AC_Het,AC_Hemi},
        q{grch37_loqusdb_-2017-05-22-.vcf.gz}            => q{Obs,Hom},
        q{grch37_genbank_haplogroup_-2015-08-01-.vcf.gz} => q{MTAF},
    },
    snpsift_annotation_outinfo_key => {
        q{grch37_gnomad.genomes_-r2.0.1-.vcf.gz}             => q{GNOMAD},
        q{grch37_anon-swegen_str_nsphs_-1000samples-.vcf.gz} => q{SWEGEN},
    },
);

my $is_ok = check_snpsift_keys(
    {
        log => $log,
        snpsift_annotation_files_href =>
          \%{ $active_parameter{snpsift_annotation_files} },
        snpsift_annotation_outinfo_key_href =>
          \%{ $active_parameter{snpsift_annotation_outinfo_key} },
    }
);

## Then all is ok
ok( $is_ok, q{All snpsift files matched} );

## Given non-matching snpsift annotation and outinfo keys
$active_parameter{snpsift_annotation_outinfo_key}{not_matching} = q{You shall not pass};

trap {
    check_snpsift_keys(
        {
            log => $log,
            snpsift_annotation_files_href =>
              \%{ $active_parameter{snpsift_annotation_files} },
            snpsift_annotation_outinfo_key_href =>
              \%{ $active_parameter{snpsift_annotation_outinfo_key} },
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if no matching snpsift annotation file} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if no matching snpsift annotation file} );

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
