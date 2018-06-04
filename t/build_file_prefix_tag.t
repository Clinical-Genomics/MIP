#!/usr/bin/env perl

use 5.018;
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
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::File::Format::Mip});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::File::Format::Mip qw{ build_file_prefix_tag };

diag(   q{Test build_file_prefix_tag from Mip.pm v}
      . $MIP::File::Format::Mip::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $sample_id     = q{homer};
my $current_chain = q{MAIN};
my $other_chain   = q{SV};

my @order_parameters =
  qw{ pbwa_mem pmerge pmark pmanta random_test_parameter random_test_parameter-2 };

my %active_parameter = (
    family_id             => q{simpsons},
    pbwa_mem              => 1,
    pmanta                => 1,
    pmark                 => 1,
    pmerge                => 0,
    random_test_parameter => undef,
    random_test_parameter => q{not_a_program},
    sample_ids            => [$sample_id],
);
my %file_info;
my %parameter = (
    dynamic_parameter => { program => [qw{ pbwa_mem pmark pmanta pmerge }], },
    pbwa_mem          => {
        chain    => $current_chain,
        file_tag => q{mem},
    },
    pmark => {
        chain    => $current_chain,
        file_tag => q{md},
    },
    pmanta => {
        chain    => $other_chain,
        file_tag => q{manta},
    },
    pmerge => {
        chain    => $current_chain,
        file_tag => q{nofile_tag},
    },
);

## Given file tags, when multiple chains
build_file_prefix_tag(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        order_parameters_ref  => \@order_parameters,
        parameter_href        => \%parameter,
    }
);

my %expected_file_tag = (
    $sample_id => {
        pbwa_mem => { file_tag => q{mem}, },
        pmark    => { file_tag => q{memmd}, },
        pmanta   => { file_tag => q{memmdmanta}, },
    },
    q{simpsons} => {
        pbwa_mem => { file_tag => q{mem}, },
        pmark    => { file_tag => q{memmd}, },
        pmanta   => { file_tag => q{memmdmanta}, },
    },
);

## Then these file tags should be set according to %expected_file_tag
is_deeply( \%file_info, \%expected_file_tag, q{Built file endings} );

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
