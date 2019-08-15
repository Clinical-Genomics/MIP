#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
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
    my @modules = (q{MIP::File::Format::Pedigree});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::File::Format::Pedigree qw{ detect_sample_id_gender };

diag(   q{Test detect_sample_id_gender from Pedigree.pm v}
      . $MIP::File::Format::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given sample ids and genders
my %active_parameter = ( sample_ids => [qw{ sample-1 sample-2 sample-3 }], );

my %sample_info = (
    sample => {
        q{sample-1} => { sex => q{male}, },
        q{sample-2} => { sex => q{female}, },
        q{sample-3} => { sex => q{other}, },
    },
);
(

    $active_parameter{found_male},
    $active_parameter{found_female},
    $active_parameter{found_other},
    $active_parameter{found_other_count},
  )
  = detect_sample_id_gender(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
  );

my %expected_result = (
    found_male        => 1,
    found_female      => 1,
    found_other       => 1,
    found_other_count => 1,
);

GENDER:
foreach my $found_gender ( keys %expected_result ) {

## Then all genders count should be one
    is( $active_parameter{$found_gender}, $expected_result{$found_gender},
        $found_gender );
}

## Given no males or females
%sample_info = (
    sample => {
        q{sample-1} => { sex => q{xyz}, },
        q{sample-2} => { sex => q{xyz}, },
        q{sample-3} => { sex => q{other}, },
    },
);

%expected_result = (
    found_male        => 1,
    found_female      => 0,
    found_other       => 1,
    found_other_count => 3,
);
(

    $active_parameter{found_male},
    $active_parameter{found_female},
    $active_parameter{found_other},
    $active_parameter{found_other_count},
  )
  = detect_sample_id_gender(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
  );

GENDER:
foreach my $found_gender ( keys %expected_result ) {

## Then one male should be found and a other count of three
    is( $active_parameter{$found_gender}, $expected_result{$found_gender},
        $found_gender );
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
