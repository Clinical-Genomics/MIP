package MIP::Language::Perl;

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
use MIP::Constants qw{ $BACKWARD_SLASH $DASH $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.11;

    our @EXPORT_OK = qw{ perl_base perl_nae_oneliners };
}

Readonly my $MINUS_ONE => -1;

sub perl_base {

## Function : Perl base and switches
## Returns  :
## Arguments: $autosplit    => Turns on autosplit mode when used with a -n or -p
##          : $command_line => Enter one line of program
##          : $n            => Iterate over filename arguments

    my ($arg_href) = @_;

    ## Flatten argument(s)

    ## Default(s)
    my $autosplit;
    my $command_line;
    my $n;

    my $tmpl = {
        autosplit => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$autosplit,
            strict_type => 1,
        },
        command_line => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$command_line,
            strict_type => 1,
        },
        n => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$n,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ perl };

    if ($n) {

        push @commands, q{-n};
    }
    if ($autosplit) {

        push @commands, q{-a};
    }
    if ($command_line) {

        push @commands, q{-e};
    }

    return @commands;
}

sub perl_nae_oneliners {

## Function : Return predifined one liners
## Returns  : @commands
## Arguments: $autosplit              => Turns on autosplit mode when used with a -n or -p
##          : $command_line           => Enter one line of program
##          : $escape_oneliner        => Escape perl oneliner program
##          : $filehandle             => Filehandle to write to
##          : $n                      => Iterate over filename arguments
##          : $oneliner_cmd           => Command to execute
##          : $oneliner_name          => Perl oneliner name
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $escape_oneliner;
    my $filehandle;
    my $oneliner_cmd;
    my $oneliner_name;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $autosplit;
    my $command_line;
    my $n;

    my $tmpl = {
        autosplit => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$autosplit,
            strict_type => 1,
        },
        command_line => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$command_line,
            strict_type => 1,
        },
        escape_oneliner => {
            allow       => [ undef, 0, 1 ],
            store       => \$escape_oneliner,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        n => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$n,
            strict_type => 1,
        },
        oneliner_cmd => {
            store       => \$oneliner_cmd,
            strict_type => 1,
        },
        oneliner_name => {
            store       => \$oneliner_name,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdinfile_path  => { store => \$stdinfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Oneliner dispatch table
    my %oneliner = (
        genepred_to_refflat                  => \&_genepred_to_refflat,
        get_dict_contigs                     => \&_get_dict_contigs,
        q{get_fastq_header_v1.4}             => \&_get_fastq_header_v1_4,
        q{get_fastq_header_v1.4_interleaved} => \&_get_fastq_header_v1_4_interleaved,
        q{get_fastq_header_v1.8}             => \&_get_fastq_header_v1_8,
        q{get_fastq_header_v1.8_interleaved} => \&_get_fastq_header_v1_8_interleaved,
        get_fastq_read_length                => \&_get_fastq_read_length,
        get_rrna_transcripts                 => \&_get_rrna_transcripts,
        get_select_contigs_by_col            => \&_get_select_contigs_by_col,
        remove_decomposed_asterisk_records   => \&_remove_decomposed_asterisk_records,
        synonyms_grch37_to_grch38            => \&_synonyms_grch37_to_grch38,
        synonyms_grch38_to_grch37            => \&_synonyms_grch38_to_grch37,
        write_contigs_size_file              => \&_write_contigs_size_file,
    );

    ## Stores commands depending on input parameters
    my @commands = perl_base(
        {
            autosplit    => $autosplit,
            command_line => $command_line,
            n            => $n,
        }
    );

    if (    defined $oneliner_name
        and exists $oneliner{$oneliner_name}
        and not $oneliner_cmd )
    {

        $oneliner_cmd = $oneliner{$oneliner_name}->();
    }

    if ( $oneliner_cmd and $escape_oneliner ) {

        $oneliner_cmd = $BACKWARD_SLASH . $oneliner_cmd;
        substr $oneliner_cmd, $MINUS_ONE, 0, $BACKWARD_SLASH;

    }

    if ($oneliner_cmd) {

        push @commands, $oneliner_cmd;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdinfile_path         => $stdinfile_path,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub _genepred_to_refflat {

## Function : Convert extended genePred format to refFlat format
## Returns  : $genepred_to_refflat
## Arguments:

    my ($arg_href) = @_;

    ## Put gene name first followed by the original first 10 fields
    my $genepred_to_refflat = q?'say STDOUT join qq{\t}, ($F[11], @F[0..9])'?;

    return $genepred_to_refflat;
}

sub _get_dict_contigs {

## Function : Return predifined one liners to get contig names from dict file
## Returns  : $get_dict_contigs
## Arguments:

    my ($arg_href) = @_;

    # Find contig line
    my $get_dict_contigs = q?'if($F[0]=~/^\@SQ/) { ?;

    # Collect contig name
    $get_dict_contigs .= q? if($F[1]=~/SN\:(\S+)/) { ?;

    # Alias capture
    $get_dict_contigs .= q?my $contig_name = $1; ?;

    # Write to STDOUT
    $get_dict_contigs .= q?print $contig_name, q{,};} }'?;

    return $get_dict_contigs;
}

sub _get_fastq_header_v1_4 {

## Function : Return regexp for header elements for fastq file format version 1.4
## Returns  : $get_fastq_header_regexp
## Arguments:

    my ($arg_href) = @_;

    ## Remove newline
    my $get_fastq_header_regexp = q?'chomp; ?;

    ## Get header features
    $get_fastq_header_regexp .=
q?my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)[\/](\d+)/; ?;

    ## If  we found correct header version
    $get_fastq_header_regexp .= q?if($instrument_id) { ?;

    # Print header features
    $get_fastq_header_regexp .=
q?print join " ", ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction);} ?;

    # Process line one and then exit
    $get_fastq_header_regexp .= q?if($.=1) {last;}' ?;

    return $get_fastq_header_regexp;
}

sub _get_fastq_header_v1_4_interleaved {

## Function : Return reg exp to get read direction for interleaved fastq file format version 1.4
## Returns  : $get_fastq_header_regexp
## Arguments:

    my ($arg_href) = @_;

    ## Remove newline
    my $get_fastq_header_regexp = q?'chomp; ?;

    ## If line one or five
    $get_fastq_header_regexp .= q?if($.==1 or $.==5) { ?;

    ## Get header features
    $get_fastq_header_regexp .=
q?my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)[\/](\d+)/; ?;

    ## If  we found correct header version
    $get_fastq_header_regexp .= q?if($instrument_id) { ?;

    # Print direction feature
    $get_fastq_header_regexp .= q?print $direction;} } ?;

    # Process to line six and then exit
    $get_fastq_header_regexp .= q?elsif ($.==6) {last;}' ?;

    return $get_fastq_header_regexp;
}

sub _get_fastq_header_v1_8 {

## Function : Return regexp for header elements for fastq interleaved file format version 1.8
## Returns  : $get_fastq_header_regexp
## Arguments:

    my ($arg_href) = @_;

    ## Remove newline
    my $get_fastq_header_regexp = q?'chomp; ?;

    ## Get header features
    $get_fastq_header_regexp .=
q?my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction, $filtered, $control_bit, $index,) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)\s(\d+):(\w+):(\d+):(\w+)/; ?;

    ## If  we found correct header version
    $get_fastq_header_regexp .= q?if($instrument_id) { ?;

    # Print header features
    $get_fastq_header_regexp .=
q?print join " ", ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction, $filtered, $control_bit, $index); } ?;

    # Process line one and then exit
    $get_fastq_header_regexp .= q?if($.=1) {last;}' ?;

    return $get_fastq_header_regexp;
}

sub _get_fastq_header_v1_8_interleaved {

## Function : Return read direction for interleaved fastq file format version 1.8
## Returns  : $get_fastq_header_regexp
## Arguments:

    my ($arg_href) = @_;

    ## Remove newline
    my $get_fastq_header_regexp = q?'chomp; ?;

    ## If line one or five
    $get_fastq_header_regexp .= q?if($.==1 or $.==5) { ?;

    ## Get header features
    $get_fastq_header_regexp .=
q?my ($instrument_id, $run_number, $flowcell, $lane, $tile, $x_pos, $y_pos, $direction, $filtered, $control_bit, $index,) = /^(@[^:]*):(\d+):(\w+):(\d+):(\d+):(\d+):(\d+)\s(\d+):(\w+):(\d+):(\w+)/; ?;

    ## If  we found correct header version
    $get_fastq_header_regexp .= q?if($instrument_id) { ?;

    # Print direction feature
    $get_fastq_header_regexp .= q?print $direction;} } ?;

    # Process to line six and then exit
    $get_fastq_header_regexp .= q?elsif ($.==6) {last;}' ?;

    return $get_fastq_header_regexp;
}

sub _get_fastq_read_length {

## Function : Return read length from a fastq infile
## Returns  : $get_fastq_read_length
## Arguments:

    my ($arg_href) = @_;

    ## Prints sequence length and exits

    # Skip header line
    my $read_length_regexp = q?'if ($_!~/@/) {?;

    # Remove newline
    $read_length_regexp .= q?chomp;?;

    # Count chars
    $read_length_regexp .= q?my $seq_length = length;?;

    # Print and exit
    $read_length_regexp .= q?print $seq_length;last;}' ?;

    return $read_length_regexp;
}

sub _get_rrna_transcripts {

## Function : Return rRNA transcripts from gtf file
## Returns  : $rrna_transcripts
## Arguments:

    my ($arg_href) = @_;

    # Print header line
    my $rrna_transcripts = q?'if (/^#/) {print $_} ?;

    # For rRNA, rRNA_pseudogenes or Mt_rRNA
    $rrna_transcripts .=
      q?elsif ($_ =~ / gene_type \s \"(rRNA|rRNA_pseudogene|Mt_rRNA)\" /nxms)?;

    # Print
    $rrna_transcripts .= q?{print $_}'?;

    return $rrna_transcripts;
}

sub _get_select_contigs_by_col {

## Function : Return contig names from column one of bed file
## Returns  : $get_select_contigs
## Arguments:

    my ($arg_href) = @_;

    # Initilize hash
    my $get_select_contigs = q?'my %contig; ?;

    # Loop per line
    $get_select_contigs .= q?while (<>) { ?;

    # Get contig name
    $get_select_contigs .= q?my ($contig_name) = $_=~/(\w+)/sxm; ?;

    # Skip if on black list
    $get_select_contigs .=
      q?next if($contig_name =~/browser|contig|chromosome|gene_panel|track/); ?;

    # Set contig name in hash
    $get_select_contigs .= q?if($contig_name) { $contig{$contig_name}=undef } } ?;

    # Print contig names to string
    $get_select_contigs .= q?print join ",", sort keys %contig; ?;

    # Exit perl program
    $get_select_contigs .= q?last;' ?;

    return $get_select_contigs;
}

sub _remove_decomposed_asterisk_records {

## Function : Remove decomposed '*' records
## Returns  : $remove_star_regexp
## Arguments:

    my ($arg_href) = @_;

    ## As long as the fourth column isn't an asterisk
    my $remove_star_regexp = q?'unless( $F[4] eq q{*} ) { ?;

    ## Print record
    $remove_star_regexp .= q?print $_ }'?;

    return $remove_star_regexp;
}

sub _synonyms_grch37_to_grch38 {

## Function : Return predifined one liners to modify chr prefix from genome version 37 to 38
## Returns  : $modify_chr_prefix
## Arguments:

    my ($arg_href) = @_;

    ## Add "chr" prefix to chromosome name and rename "M" to "MT"
    my $modify_chr_prefix = q?'if($_=~s/^M/chrMT/g) {} ?;

    ## Add "chr" prefix to chromosome name
    $modify_chr_prefix .= q?elsif ($_=~s/^(.+)/chr$1/g) {} ?;

## Print line
    $modify_chr_prefix .= q?print $_'?;

    return $modify_chr_prefix;
}

sub _synonyms_grch38_to_grch37 {

## Function : Return predifined one liners to modify chr prefix from genome version 37 to 38
## Returns  : $modify_chr_prefix
## Arguments:

    my ($arg_href) = @_;

## Remove "chr" prefix from chromosome name and rename "MT" to "M"
    my $modify_chr_prefix = q?'if($_=~s/^chrMT/M/g) {} ?;

## Remove "chr" prefix from chromosome name
    $modify_chr_prefix .= q?elsif ($_=~s/^chr(.+)/$1/g) {} ?;

## Print line
    $modify_chr_prefix .= q?print $_'?;

    return $modify_chr_prefix;
}

sub _write_contigs_size_file {

## Function : Return predifined one liners to write contig names and length from fai file
## Returns  : $write_contigs_size
## Arguments:

    my ($arg_href) = @_;

    ## Contig name ($F[0]), contig length ($F[1])
    my $write_contigs_size = q?'say STDOUT $F[0] . "\t" . $F[1] '?;

    return $write_contigs_size;
}

1;
