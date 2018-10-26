#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Getopt::Long;
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
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Set::Pedigree});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::Set::Pedigree qw{ set_pedigree_phenotype_info };

diag(   q{Test set_pedigree_phenotype_info from Pedigree.pm v}
      . $MIP::Set::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %pedigree = (
    case    => q{case_1},
    samples => [
        {
            analysis_type => q{wes},
            father        => 0,
            mother        => 0,
            phenotype     => q{affected},
            sample_id     => q{sample_1},
            sex           => q{female},
        },
        {
            analysis_type => q{wgs},
            father        => 0,
            mother        => 0,
            phenotype     => q{unaffected},
            sample_id     => q{sample_2},
            sex           => q{male},
        },
        {
            analysis_type => q{wts},
            father        => 0,
            mother        => 0,
            phenotype     => q{unknown},
            sample_id     => q{sample_3},
            sex           => q{other},
        },
        {
            analysis_type => q{wgs},
            father        => q{sample_1},
            mother        => q{sample_2},
            phenotype     => q{unknown},
            sample_id     => q{sample_4},
            sex           => q{unknown},
        },
    ],
);
my %parameter;
set_pedigree_phenotype_info(
    {
        pedigree_href  => \%pedigree,
        parameter_href => \%parameter,
    }
);

my @got_unaffected_samples = @{ $parameter{dynamic_parameter}{unaffected} };
my @got_affected_samples   = @{ $parameter{dynamic_parameter}{affected} };
my @got_unknown_samples    = @{ $parameter{dynamic_parameter}{unknown} };

is( scalar @got_unaffected_samples, 1, q{Got all samples with unaffected phenotype} );

is( scalar @got_affected_samples, 1, q{Got all samples with affected phenotype} );

is( scalar @got_unknown_samples, 2, q{Got all samples with unknown phenotype} );

my $affected_plink_phenotype =
  $parameter{dynamic_parameter}{sample_1}{plink_phenotype};
my $unaffected_plink_phenotype =
  $parameter{dynamic_parameter}{sample_2}{plink_phenotype};
my $unknown_plink_phenotype =
  $parameter{dynamic_parameter}{sample_3}{plink_phenotype};

is( $affected_plink_phenotype, 2, q{Reformated to plink affected phenotype} );

is( $unaffected_plink_phenotype, 1, q{Reformated to plink unaffected phenotype} );

is( $unknown_plink_phenotype, 0, q{Reformated to plink unknown phenotype} );

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
