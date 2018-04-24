#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use Getopt::Long;
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $DOT        => q{.};
Readonly my $COMMA      => q{,};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

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
    my @modules = (q{MIP::QC::Record});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::QC::Record qw{ add_most_complete_vcf };

diag(   q{Test add_most_complete_vcf from Record.pm v}
      . $MIP::QC::Record::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $path              = catfile(qw{ a test vcf_file });
my $program_name_test = q{sv_rankvariants};
my $file_suffix       = $DOT . q{vcf};
my $vcf_file_key =
  q{sv} . $UNDERSCORE . substr( $file_suffix, 1 ) . $UNDERSCORE . q{file};

my %active_parameter = ( q{p} . $program_name_test => 1, );

my %sample_info;
add_most_complete_vcf(
    {
        active_parameter_href     => \%active_parameter,
        path                      => $path,
        program_name              => $program_name_test,
        sample_info_href          => \%sample_info,
        vcfparser_outfile_counter => 0,
        vcf_file_key              => $vcf_file_key,
    }
);
is( $sample_info{$vcf_file_key}{research}{path}, $path,
    q{Added research path} );

add_most_complete_vcf(
    {
        active_parameter_href     => \%active_parameter,
        path                      => $path,
        program_name              => $program_name_test,
        sample_info_href          => \%sample_info,
        vcfparser_outfile_counter => 1,
        vcf_file_key              => $vcf_file_key,
    }
);
is( $sample_info{$vcf_file_key}{clinical}{path}, $path,
    q{Added clinical path} );

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
