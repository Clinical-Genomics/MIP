package MIP::Program::Alignment::Samtools;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use Carp;
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };

use FindBin qw{$Bin};    # Find directory of script
use File::Basename qw{dirname};
use File::Spec::Functions qw{catdir};

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Inherit from Exporter to export functions and variables
    use base qw {Exporter};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{samtools_view samtools_index samtools_stats samtools_mpileup samtools_faidx};

}

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw{unix_standard_streams};
use MIP::Unix::Write_to_file qw{unix_write_to_file};

use Params::Check qw{check allow last_error};
use Readonly;

## Constants
Readonly my $SPACE => q{ };
Readonly my $COMMA => q{,};

sub samtools_view {

## samtools_view

## Function : Perl wrapper for writing samtools view recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
## Returns  : "@commands"
## Arguments: $regions_ref, $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $thread_number, $with_header, $output_format, $auto_detect_input_format, $uncompressed_bam_output
##          : $regions_ref              => The regions to process {REF}
##          : $infile_path              => Infile path
##          : $outfile_path             => Outfile path
##          : $stderrfile_path          => Stderrfile path
##          : $FILEHANDLE               => Sbatch filehandle to write to
##          : $thread_number            => Number of BAM/CRAM compression threads
##          : $with_header              => Include header
##          : $output_format            => Output format
##          : $auto_detect_input_format => Ignored (input format is auto-detected)
##          : $uncompressed_bam_output  => Uncompressed bam output
##          : $stderrfile_path_append   => Stderrfile path append

    my ($arg_href) = @_;

    ## Default(s)
    my $with_header;
    my $output_format;
    my $auto_detect_input_format;
    my $uncompressed_bam_output;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $thread_number;
    my $stderrfile_path_append;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        thread_number => {
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$thread_number
        },
        FILEHANDLE  => { store => \$FILEHANDLE },
        with_header => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$with_header
        },
        output_format => {
            default     => q{bam},
            allow       => [qw{sam bam cram json}],
            strict_type => 1,
            store       => \$output_format
        },
        auto_detect_input_format => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$auto_detect_input_format
        },
        uncompressed_bam_output => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$uncompressed_bam_output
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{samtools view};

    ## Options
    if ($thread_number) {

        #Number of threads
        push @commands, q{--threads} . $SPACE . $thread_number;
    }

    if ($with_header) {

        #Include header
        push @commands, q{-h};
    }

    if ($output_format) {

        #Output format
        push @commands, q{--output-fmt} . $SPACE . uc $output_format;
    }

    if ($auto_detect_input_format) {

        push @commands, q{-S};
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{-o} . $SPACE . $outfile_path;
    }

    if ($uncompressed_bam_output) {

        push @commands, q{-u};
    }

    ## Infile
    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        #Limit output to regions
        push @commands, join $SPACE, @{$regions_ref};
    }

    # Redirect stderr output to program specific stderr file
    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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

sub samtools_index {

## samtools_index

## Function : Perl wrapper for writing samtools index recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
## Returns  : "@commands"
## Arguments: $infile_path, $stderrfile_path, $FILEHANDLE, $bai_format
##          : $infile_path              => Infile path
##          : $stderrfile_path          => Stderrfile path
##          : $FILEHANDLE               => Sbatch filehandle to write to
##          : $bai_format               => Generate BAI-format index for BAM files
##          : $stderrfile_path_append   => Stderrfile path append

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $stdout_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $bai_format;
    my $stderrfile_path_append;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        stdout_path     => { strict_type => 1, store => \$stdout_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE => { store       => \$FILEHANDLE },
        bai_format => { strict_type => 1, store => \$bai_format },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{samtools index};

    ## Options
    if ($bai_format) {

        #Generate BAI-format index for BAM files
        push @commands, q{-b};
    }

    ## Infile
    push @commands, $infile_path;

    # Redirect stderr output to program specific stderr file
    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdout_path            => $stdout_path
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

sub samtools_stats {

## samtools_stats

## Function : Perl wrapper for writing samtools stats recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
## Returns  : "@commands"
## Arguments: $regions_ref, $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $auto_detect_input_format
##          : $regions_ref              => The regions to process {REF}
##          : $infile_path              => Infile path
##          : $outfile_path             => Outfile path
##          : $stderrfile_path          => Stderrfile path
##          : $FILEHANDLE               => Sbatch filehandle to write to
##          : $auto_detect_input_format => Ignored (input format is auto-detected)
##          : $stderrfile_path_append   => Stderrfile path append

    my ($arg_href) = @_;

    ## Default(s)
    my $auto_detect_input_format;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $stderrfile_path_append;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE               => { store => \$FILEHANDLE },
        auto_detect_input_format => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$auto_detect_input_format
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{samtools stats};

    if ($auto_detect_input_format) {

        push @commands, q{-s};
    }

    ## Infile
    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        #Limit output to regions
        push @commands, join $SPACE, @{$regions_ref};
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{>} . $SPACE . $outfile_path;
    }

    # Redirect stderr output to program specific stderr file
    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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

sub samtools_mpileup {

## samtools_mpileup

## Function : Perl wrapper for writing samtools mpileup recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
## Returns  : "@commands"
## Arguments: $infile_paths_ref, $output_tags_ref, $outfile_path, $referencefile_path, $stderrfile_path, $FILEHANDLE, $output_bcf, $adjust_mq
##          : $infile_paths_ref                 => Infile paths {REF}
##          : $output_tags_ref                  => Optional tags to output {REF}
##          : $outfile_path                     => Outfile path
##          : $referencefile_path               => Reference sequence file
##          : $stderrfile_path                  => Stderrfile path
##          : $FILEHANDLE                       => Sbatch filehandle to write to
##          : $region                           => The regions to process {REF}
##          : $output_bcf                       => Generate genotype likelihoods in BCF format
##          : $per_sample_increased_sensitivity => Apply -m and -F per-sample for increased sensitivity
##          : $adjust_mq                        => Adjust mapping quality
##          : $stderrfile_path_append           => Stderrfile path append

    my ($arg_href) = @_;

    ## Default(s)
    my $per_sample_increased_sensitivity;
    my $adjust_mq;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $output_tags_ref;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $region;
    my $output_bcf;
    my $stderrfile_path_append;

    my $tmpl = {
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        output_tags_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$output_tags_ref
        },
        outfile_path       => { strict_type => 1, store => \$outfile_path },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE => { store       => \$FILEHANDLE },
        region     => { strict_type => 1, store => \$region },
        output_bcf => { strict_type => 1, store => \$output_bcf },
        per_sample_increased_sensitivity => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$per_sample_increased_sensitivity
        },
        adjust_mq => {
            default     => 50,
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$adjust_mq
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{samtools mpileup};

    ## Options
    push @commands, q{--adjust-MQ} . $SPACE . $adjust_mq;

    if ($per_sample_increased_sensitivity) {

        push @commands, q{--per-sample-mF};
    }

    if ( @{$output_tags_ref} ) {

        push @commands, q{--output-tags} . $SPACE . join $COMMA,
          @{$output_tags_ref};
    }

    if ($region) {

        #Limit output to region
        push @commands, q{--region} . $SPACE . $region;
    }

    if ($referencefile_path) {

        #Reference sequence file
        push @commands, q{--fasta-ref} . $SPACE . $referencefile_path;
    }

    if ($output_bcf) {

        push @commands, q{--BCF};
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    ## Infile
    push @commands, join $SPACE, @{$infile_paths_ref};

    # Redirect stderr output to program specific stderr file
    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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

sub samtools_faidx {

## samtools_faidx

## Function : Perl wrapper for writing samtools faidx recipe to $FILEHANDLE. Based on samtools 1.3.1 (using htslib 1.3.1).
## Returns  : "@commands"
## Arguments: $regions_ref, $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE
##          : $regions_ref                      => The regions to process {REF}
##          : $infile_path                      => Infile path
##          : $outfile_path                     => Outfile path
##          : $stderrfile_path                  => Stderrfile path
##          : $FILEHANDLE                       => Sbatch filehandle to write to
##          : $stderrfile_path_append           => Stderrfile path append

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $FILEHANDLE;
    my $stderrfile_path_append;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{samtools faidx};

    ## Infile
    push @commands, $infile_path;

    if ( @{$regions_ref} ) {

        #Limit output to regions
        push @commands, join $SPACE, @{$regions_ref};
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{>} . $SPACE . $outfile_path;
    }

    # Redirect stderr output to program specific stderr file
    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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
