#!/usr/bin/env perl

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp qw{ tempdir tempfile };
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
    my @modules = (q{MIP::Check::Path});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Check::Path qw{ check_filesystem_objects_and_index_existance };

diag(   q{Test check_filesystem_objects_and_index_existance from Path.pm v}
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

## Create a temp dir
my $dir = tempdir( CLEANUP => 1 );

## Create a temp file using newly created temp dir
my ( $fh, $file_name ) = tempfile( DIR => $dir );

my %parameter = (
    dir       => $dir,
    file_name => { build_file => 1 },
);

my %active_parameter = (
    gatk_baserecalibration_known_sites =>
      catfile( $Bin, qw{ data references GRCh37_dbsnp_-138-.vcf} ),
    human_genome_reference =>
      catfile( $Bin, qw{data references GRCh37_homo_sapiens_-d5-.fasta} ),
    snpsift_annotation_files =>
      catfile( $Bin, qw{data references GRCh37_clinvar_-2017-01-04-.vcf.gz} ),
);

### TEST
my ($exist) = check_filesystem_objects_and_index_existance(
    {
        log         => $log,
        object_name => $active_parameter{human_genome_reference},
        ,
        object_type    => q{file},
        parameter_href => \%parameter,
        parameter_name => q{human_genome_reference},
        path           => $active_parameter{human_genome_reference},
    }
);

is( $exist, undef, q{Return for build file} );

($exist) = check_filesystem_objects_and_index_existance(
    {
        log            => $log,
        object_name    => $active_parameter{gatk_baserecalibration_known_sites},
        object_type    => q{file},
        parameter_href => \%parameter,
        parameter_name => q{gatk_baserecalibration_known_sites},
        path           => $active_parameter{gatk_baserecalibration_known_sites},
    }
);

is( $exist, undef, q{Return for file} );

($exist) = check_filesystem_objects_and_index_existance(
    {
        log            => $log,
        object_name    => $active_parameter{snpsift_annotation_files},
        object_type    => q{file},
        parameter_href => \%parameter,
        parameter_name => q{snpsift_annotation_files},
        path           => $active_parameter{snpsift_annotation_files},
    }
);

is( $exist, undef, q{Return for index file} );

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
