package MIP::Program::Base::Bcftools;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ bcftools_base };
}

sub bcftools_base {

## Function : Perl wrapper for bcftools base. Based on Bcftools 1.9
## Returns  : @commands
## Arguments: $commands_ref      => List of commands added earlier
##          : $filehandle        => Filehandle to write to
##          : $outfile_path      => Outfile path
##          : $output_type       => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref       => Regions to process {REF}
##          : $samples_file_path => File of samples to annotate
##          : $samples_ref       => Samples to include or exclude if prefixed with "^"
##          : $threads           => Extra compression threds in addition to main thread

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $filehandle;
    my $outfile_path;
    my $output_type;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $threads;

    my $tmpl = {
        commands_ref => {
            default     => [],
            store       => \$commands_ref,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        output_type => {
            allow       => [ undef, qw{ b u z v} ],
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => {
            default     => [],
            store       => \$regions_ref,
            strict_type => 1,
        },
        samples_file_path => {
            store       => \$samples_file_path,
            strict_type => 1,
        },
        samples_ref => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        threads => {
            allow       => qr/\A \d+ \z | undef /xms,
            store       => \$threads,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = @{$commands_ref};

    if ($samples_file_path) {

        push @commands, q{--samples-file} . $SPACE . $samples_file_path;
    }
    if ( @{$samples_ref} ) {

        push @commands, q{--samples} . $SPACE . join $COMMA, @{$samples_ref};
    }
    if ( @{$regions_ref} ) {

        # Limit output to regions
        push @commands, q{--regions} . $SPACE . join $COMMA, @{$regions_ref};
    }
    if ($outfile_path) {

        # Specify output filename
        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    if ($output_type) {

        #Specify output type
        push @commands, q{--output-type} . $SPACE . $output_type;
    }

    if ($threads) {

        push @commands, q{--threads} . $SPACE . $threads;
    }

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            filehandle   => $filehandle,
        }
    );
    return @commands;
}

1;
