package MIP::Program::Interval::Picardtools;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };
use MIP::Language::Java qw{java_core};
use MIP::Program::Base::Picardtools qw{ picardtools_base};

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ picardtools_intervallisttools };
}

## Constants
Readonly my $SPACE => q{ };

sub picardtools_intervallisttools {

## Function : Perl wrapper for writing picardtools intervallisttools recipe to $FILEHANDLE. Based on picardtools 2.5.0.
##Returns  : @commands
## Arguments: $infile_paths_ref     => Infile paths {REF}
##          : $outfile_path         => Outfile path
##          : $referencefile_path   => Genome reference file
##          : $FILEHANDLE           => Sbatch filehandle to write to
##          : $stderrfile_path      => Stderrfile path
##          : $memory_allocation    => Memory allocation for java
##          : $temp_directory       => Redirect tmp files to java temp
##          : $java_use_large_pages => Use java large pages
##          : $java_jar             => Java jar
##          : $create_index         => Create index
##          : $padding              => The amount to pad each end of the intervals by before other operations are undertaken

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $referencefile_path;
    my $FILEHANDLE;
    my $stderrfile_path;
    my $memory_allocation;
    my $temp_directory;
    my $java_jar;

    ## Default(s)
    my $java_use_large_pages;
    my $create_index;
    my $padding;

    my $tmpl = {
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        FILEHANDLE        => { store       => \$FILEHANDLE },
        stderrfile_path   => { strict_type => 1, store => \$stderrfile_path },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        temp_directory    => { strict_type => 1, store => \$temp_directory },
        java_jar          => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        create_index => {
            default     => q{false},
            allow       => [qw{ true false }],
            strict_type => 1,
            store       => \$create_index
        },
        padding => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$padding
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands;

    ## Return java core commands
    if ($java_jar) {

        @commands = java_core(
            {
                memory_allocation    => $memory_allocation,
                java_use_large_pages => $java_use_large_pages,
                temp_directory       => $temp_directory,
                java_jar             => $java_jar,
            }
        );
    }

    ## Picardtools intervallisttools
    push @commands, q{IntervalListTools};

    ## Picardtools base args
    @commands = picardtools_base(
        {
            commands_ref       => \@commands,
            referencefile_path => $referencefile_path,
            create_index       => $create_index,
        }
    );

    ##Options
    if ($padding) {

        push @commands, q{PADDING=} . $padding;
    }

    ## Infile
    push @commands, q{INPUT=} . join $SPACE . q{INPUT=}, @{$infile_paths_ref};

    ## Output
    push @commands, q{OUTPUT=} . $outfile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

1;

