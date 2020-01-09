#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp qw{ tempdir tempfile };
use FindBin qw{ $Bin };
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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Check::Path} => [qw{ check_filesystem_objects_and_index_existance }],
        q{MIP::File::Format::Yaml} => [qw{ load_yaml }],
        q{MIP::Test::Fixtures}     => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
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
my $test_dir      = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = test_log( { no_screen => 0, } );

my %active_parameter = (
    gatk_baserecalibration_known_sites =>
      catfile( $Bin, qw{ data references grch37_dbsnp_-138-.vcf} ),
    human_genome_reference_file_endings => catfile( $Bin, qw{ Should_not_exist_file } ),
    gatk_variantrecalibration_resource_indel =>
      catfile( $Bin, qw{ data references grch37_dbsnp_-138-.vcf } ),
);

my %parameter = load_yaml(
    {
        yaml_file => catfile( dirname($Bin), qw{ definitions rd_dna_parameters.yaml} ),
    }
);

### Given a reference which can be built
my ($exist) = check_filesystem_objects_and_index_existance(
    {
        log         => $log,
        object_name => $active_parameter{human_genome_reference_file_endings},
        ,
        object_type    => q{file},
        parameter_href => \%parameter,
        parameter_name => q{human_genome_reference_file_endings},
        path           => $active_parameter{human_genome_reference_file_endings},
    }
);

## Then do not check anything and return undef
is( $exist, undef, q{Return for build file} );

## Given a file that exists
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

## Then return true
ok( $exist, q{File exists} );

## Given a file with index
($exist) = check_filesystem_objects_and_index_existance(
    {
        log            => $log,
        object_name    => $active_parameter{gatk_variantrecalibration_resource_indel},
        object_type    => q{file},
        parameter_href => \%parameter,
        parameter_name => q{gatk_variantrecalibration_resource_indel},
        path           => $active_parameter{gatk_variantrecalibration_resource_indel},
    }
);

## Then return true
ok( $exist, q{Index file exists} );

## Given a file that do not exist
$active_parameter{gatk_variantrecalibration_resource_indel} = q{file_do_not_exist};

trap {
    check_filesystem_objects_and_index_existance(
        {
            log            => $log,
            object_name    => $active_parameter{gatk_variantrecalibration_resource_indel},
            object_type    => q{file},
            parameter_href => \%parameter,
            parameter_name => q{gatk_variantrecalibration_resource_indel},
            path           => $active_parameter{gatk_variantrecalibration_resource_indel},
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if file cannot be found} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message if file cannot be found} );

## Given a file with an index that do not exist
$active_parameter{gatk_variantrecalibration_resource_indel} =
  catfile( $Bin, qw{ data references index_do_not_exist.gz} );

trap {
    check_filesystem_objects_and_index_existance(
        {
            log            => $log,
            object_name    => $active_parameter{gatk_variantrecalibration_resource_indel},
            object_type    => q{file},
            parameter_href => \%parameter,
            parameter_name => q{gatk_variantrecalibration_resource_indel},
            path           => $active_parameter{gatk_variantrecalibration_resource_indel},
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if file index cannot be found} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if file index cannot be found} );

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
