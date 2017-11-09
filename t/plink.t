#!/usr/bin/env perl

use Modern::Perl qw{ 2014 };
use warnings qw{ FATAL utf8 };
use autodie;
use 5.018;
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

use FindBin qw{ $Bin };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile };
use Getopt::Long;
use Test::More;

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = 1.0.0;

## Constants
Readonly my $SPACE               => q{ };
Readonly my $NEWLINE             => qq{\n};
Readonly my $COMMA               => q{,};
Readonly my $INDEP_WINDOW_SIZE   => 50;
Readonly my $INDEP_STEP_SIZE     => 5;
Readonly my $INDEP_VIF_THRESHOLD => 2;

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
    my %perl_module;

    $perl_module{q{MIP::Script::Utils}} = [qw{ help }];

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Program::Variantcalling::Plink});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Program::Variantcalling::Plink qw{ plink };
use MIP::Test::Commands qw{ test_function };

diag(   q{Test plink from Plink v}
      . $MIP::Program::Variantcalling::Plink::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $function_base_command = q{plink2};

my %base_argument = (
    stdoutfile_path => {
        input           => q{stdoutfile.test},
        expected_output => q{1> stdoutfile.test},
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{2> stderrfile.test},
    },
    FILEHANDLE => {
        input           => undef,
        expected_output => $function_base_command,
    },
);

## Can be duplicated with %base_argument and/or %specific_argument
## to enable testing of each individual argument
my %required_argument = ();

my %specific_argument = (
    regions_ref => {
        inputs_ref      => [qw{ 1 2 3 }],
        expected_output => q{--chr 1,2,3},
    },
    vcffile_path => {
        input           => catfile(qw{ dir infile.vcf }),
        expected_output => q{--vcf} . $SPACE . catfile(qw{ dir infile.vcf }),
    },
    outfile_prefix => {
        input           => catfile(qw{ temp_directory $family_id _data }),
        expected_output => q{--out}
          . $SPACE
          . catfile(qw{ temp_directory $family_id _data }),
    },
    make_bed => {
        input           => 1,
        expected_output => q{--make-bed},
    },
    binary_fileset_prefix => {
        input           => catfile(qw{ temp_directory $family_id _data }),
        expected_output => q{--bfile}
          . $SPACE
          . catfile(qw{ temp_directory $family_id _data }),
    },
    split_x => {
        input           => q{hg19},
        expected_output => q{--split-x} . $SPACE . q{hg19},
    },
    ## subcondition of $split_x
    #no_fail => {
    #    input           => 1,
    #    expected_output => q{no-fail},
    #},
    set_missing_var_ids => {
        input           => q?@:#[hg19]\$1,\$2?,
        expected_output => q{--set-missing-var-ids}
          . $SPACE
          . q?@:#[hg19]\$1,\$2?,
    },
    indep => {
        input           => 1,
        expected_output => q{--indep},
    },
    ## subconditions of INDEP
    #indep_window_size => {
    #    input           => $INDEP_WINDOW_SIZE,
    #    expected_output => $INDEP_WINDOW_SIZE,
    #},
    #indep_step_size => {
    #    input           => $INDEP_STEP_SIZE,
    #    expected_output => $INDEP_STEP_SIZE,
    #},
    #indep_vif_threshold => {
    #    input           => $INDEP_VIF_THRESHOLD,
    #    expected_output => $INDEP_VIF_THRESHOLD,
    #},
    extract_file => {
        input => catfile(qw{ temp_directory, family_id _data.prune.in }),
        expected_output => q{--extract}
          . $SPACE
          . catfile(qw{ temp_directory, family_id _data.prune.in }),
    },
    fam_file_path => {
        input => catfile(qw{ outfamily_file_directory family_id .fam }),
        expected_output => q{--fam}
          . $SPACE
          . catfile(qw{ outfamily_file_directory family_id .fam }),
    },
    ped_file_path => {
        input => catfile(qw{ temp_directory family_id_ref _data.ped }),
        expected_output => q{--ped}
          . $SPACE
          . catfile(qw{ temp_directory family_id_ref _data.ped }),
    },
    map_file_path => {
        input => catfile(qw{ temp_directory family_id_ref _data.map }),
        expected_output => q{--map}
          . $SPACE
          . catfile(qw{ temp_directory family_id_ref _data.map }),
    },
    keep_file_path => {
        input           => catfile(qw{ temp_directory family_id_ref file }),
        expected_output => q{--keep}
          . $SPACE
          . catfile(qw{ temp_directory family_id_ref file }),
    },
    read_freqfile_path => {
        input           => catfile(qw{ temp_directory family_id _data.frqx }),
        expected_output => q{--read-freq}
          . $SPACE
          . catfile(qw{ temp_directory family_id _data.frqx }),
    },
    const_fid => {
        input           => q{family_id},
        expected_output => q{--const-fid} . $SPACE . q{family_id},
    },
    check_sex => {
        input           => 1,
        expected_output => q{--check-sex},
    },
    ## subcondition of sex check
    #sex_check_min_f => {
    #    input           => q{family_id},
    #    expected_output => q{family_id},
    #},
    #y_only => {
    #    input           => 1,
    #    expected_output => q{y-only},
    #},
    vcf_require_gt => {
        input           => 1,
        expected_output => q{--vcf-require-gt},
    },
    het => {
        input           => 1,
        expected_output => q{--het},
    },
    ## subcondition of $het
    #small_sample => {
    #    input           => 1,
    #    expected_output => q{small-sample},
    #},
    inbreeding_coefficients => {
        input           => 1,
        expected_output => q{--ibc},
    },
    recode => {
        input           => 1,
        expected_output => q{--recode},
    },
    make_just_fam => {
        input           => 1,
        expected_output => q{--make-just-fam},
    },
    cluster => {
        input           => 1,
        expected_output => q{--cluster},
    },
    matrix => {
        input           => 1,
        expected_output => q{--matrix},
    },
    freqx => {
        input           => 1,
        expected_output => q{--freqx},
    },
    allow_no_sex => {
        input           => 1,
        expected_output => q{--allow-no-sex},
    },
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&plink;

## Test both base and function specific arguments
my @arguments = ( \%base_argument, \%specific_argument );

ARGUMENT_HASH_REF:
foreach my $argument_href (@arguments) {
    my @commands = test_function(
        {
            argument_href          => $argument_href,
            required_argument_href => \%required_argument,
            module_function_cref   => $module_function_cref,
            function_base_command  => $function_base_command,
            do_test_base_command   => 1,
        }
    );
}

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
