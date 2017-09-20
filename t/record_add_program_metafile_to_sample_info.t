#!/usr/bin/env perl

#### Copyright 2017 Henrik Stranneheim

use Modern::Perl qw(2014);
use warnings qw(FATAL utf8);
use autodie;
use 5.018;    # Require at least perl 5.18
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);
use Params::Check qw(check allow last_error);

use FindBin qw($Bin);    # Find directory of script
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catfile catdir devnull);
use Getopt::Long;
use Test::More;

use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use Script::Utils qw(help);

our $USAGE = build_usage( {} );

my $VERBOSE = 0;
our $VERSION = '1.0.0';

###User Options
GetOptions(
    'h|help' => sub {
        done_testing();
        print {*STDOUT} $USAGE, "\n";
        exit;
    },    #Display help text
    'v|version' => sub {
        done_testing();
        print {*STDOUT} "\n" . basename($PROGRAM_NAME) . q{  } . $VERSION,
          "\n\n";
        exit;
    },    #Display version number
    'vb|verbose' => $VERBOSE,
  )
  or (
    done_testing(),
    Script::Utils::help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
    ## Modules with import
    my %perl_module = ( 'Script::Utils' => [qw{help}], );

    while ( my ( $module, $module_import ) = each %perl_module ) {

        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load } . $module;
    }

    ## Modules
    my @modules = (qw{MIP::QC::Record});

    for my $module (@modules) {

        require_ok($module) or BAIL_OUT q{Cannot load } . $module;
    }
}

use MIP::QC::Record qw(add_program_metafile_to_sample_info);

diag(
"Test add_program_metafile_to_sample_info $MIP::QC::Record::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

# Init hash
my %sample_info;

# Test variables
my $test_program_name = q{test_program};
my $metafile          = q{test_metafile};
my $directory         = q{test_directory};
my $file              = q{test.yaml};
my $path              = catfile( $directory, $file );
my $version           = q{1.0.1};

## Family level
add_program_metafile_to_sample_info(
    {
        sample_info_href => \%sample_info,
        program_name     => $test_program_name,
        metafile_tag     => $metafile,
        directory        => $directory,
        file             => $file,
        path             => $path,
        version          => $version,
    }
);

## Test
is( exists $sample_info{program}{$test_program_name}{$metafile},
    1, q{Created family level hash key} );

is( $sample_info{program}{$test_program_name}{$metafile}{directory},
    $directory, q{Assigned correct value to family level directory} );

is( $sample_info{program}{$test_program_name}{$metafile}{file},
    $file, q{Assigned correct value to family level file} );

is( $sample_info{program}{$test_program_name}{$metafile}{path},
    $path, q{Assigned correct value to family level path} );

is( $sample_info{program}{$test_program_name}{$metafile}{version},
    $version, q{Assigned correct value to family level version} );

## Sample level
my $sample_id = q{test_sample_id};
my $infile    = q{test_infile};

add_program_metafile_to_sample_info(
    {
        sample_info_href => \%sample_info,
        sample_id        => $sample_id,
        infile           => $infile,
        program_name     => $test_program_name,
        metafile_tag     => $metafile,
        directory        => $directory,
        file             => $file,
        path             => $path,
        version          => $version,
    }
);

## Test
is(
    exists $sample_info{sample}{$sample_id}{program}{$test_program_name}
      {$infile}{$metafile},
    1,
    q{Created sample level hash key}
);

is(
    $sample_info{sample}{$sample_id}{program}{$test_program_name}
      {$infile}{$metafile}{directory},
    $directory, q{Assigned correct value to sample level directory}
);

is(
    $sample_info{sample}{$sample_id}{program}{$test_program_name}
      {$infile}{$metafile}{file},
    $file, q{Assigned correct value to sample level file}
);

is(
    $sample_info{sample}{$sample_id}{program}{$test_program_name}
      {$infile}{$metafile}{path},
    $path, q{Assigned correct value to sample level path}
);

is(
    $sample_info{sample}{$sample_id}{program}{$test_program_name}
      {$infile}{$metafile}{version},
    $version, q{Assigned correct value to sample level version}
);

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

##build_usage

##Function : Build the USAGE instructions
##Returns  : ""
##Arguments: $program_name
##         : $program_name => Name of the script

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

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}
