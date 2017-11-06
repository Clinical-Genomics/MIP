package MIP::Program::Variantcalling::Freebayes;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ freebayes_calling };
}

## Constants
Readonly my $SPACE => q{ };

sub freebayes_calling {

## Function : Perl wrapper for generic commands module.
## Returns  : @commands

## Arguments: $infile_paths_ref           => Infile paths {REF}
##          : $referencefile_path         => Reference sequence file
##          : $outfile_path               => Outfile path
##          : $stderrfile_path            => Stderrfile path
##          : $FILEHANDLE                 => Filehandle to write to
##          : $apply_standard_filter      => Use stringent input base and mapping quality filters. Equivalent to -m 30 -q 20 -R 0 -S 0
##          : $calculate_genotype_quality => Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $referencefile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $apply_standard_filter;
    my $calculate_genotype_quality;

    my $tmpl = {
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        FILEHANDLE            => { store => \$FILEHANDLE },
        apply_standard_filter => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$apply_standard_filter
        },
        calculate_genotype_quality => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$calculate_genotype_quality
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{freebayes};

    ## Options
    if ($apply_standard_filter) {

        #Equivalent to -m 30 -q 20 -R 0 -S 0
        push @commands, q{--standard-filters};
    }

    if ($calculate_genotype_quality) {

        push @commands, q{--genotype-qualities};
    }

    if ($referencefile_path) {

        #Reference sequence file
        push @commands, q{--fasta-reference} . $SPACE . $referencefile_path;
    }

    ## Infile
    push @commands, join $SPACE, @{$infile_paths_ref};

    #Specify output filename
    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
    }

    push @commands, unix_standard_streams(
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
