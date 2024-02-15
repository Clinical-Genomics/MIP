package MIP::Program::Fastp;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ fastp };
}

Readonly my $BASE_COMMAND => q{fastp};

sub fastp {

## Function : Perl wrapper for Fastp version 0.23.4
## Returns  : @commands
## Arguments: $detect_pe_adapter           => Detect adapter for parired end reads
##          : $dont_eval_duplication       => Don't evaluate duplication rate to saves time and memory
##          : $filehandle                  => Filehandle to write to
##          : $first_infile_path           => Path to input read 1 fastq
##          : $first_outfile_path          => Path to output read 1 fastq
##          : $interleaved_in              => Interleaved fastq file
##          : $length_required             => Minimum required length after trimming
##          : $low_complexity_filter       => Filter low complexity regions
##          : $overrepresentation_analysis => Enable overrepresentation analysis
##          : $report_html                 => Trimming report in html format
##          : $report_json                 => Trimming report in json format
##          : $second_infile_path          => Path to input read 2 fastq
##          : $second_outfile_path         => Path to output read 2 fastq
##          : $stderrfile_path             => Stderrfile path
##          : $stderrfile_path_append      => Append stderr info to file path
##          : $stdoutfile_path             => Stdoutfile path
##          : $thread                      => Number of threads
##          : $trim_poly_g                 => Force poly g trimming

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $detect_pe_adapter;
    my $dont_eval_duplication;
    my $filehandle;
    my $first_infile_path;
    my $first_outfile_path;
    my $interleaved_in;
    my $length_required;
    my $low_complexity_filter;
    my $overrepresentation_analysis;
    my $report_html;
    my $report_json;
    my $second_infile_path;
    my $second_outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $threads;
    my $trim_poly_g;

    ## Default(s)

    my $tmpl = {
        detect_pe_adapter => {
            allow       => [ 0, 1 ],
            store       => \$detect_pe_adapter,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        first_infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$first_infile_path,
            strict_type => 1,
        },
        first_outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$first_outfile_path,
            strict_type => 1,
        },
        interleaved_in => {
            allow       => [ 0, 1 ],
            store       => \$interleaved_in,
            strict_type => 1,
        },
        length_required => {
            allow       => qr/\A \d+ \z/xms,
            store       => \$length_required,
            strict_type => 1,
        },
        low_complexity_filter => {
            allow       => [ 0, 1 ],
            store       => \$low_complexity_filter,
            strict_type => 1,
        },
        overrepresentation_analysis => {
            allow       => [ 0, 1 ],
            store       => \$overrepresentation_analysis,
            strict_type => 1,
        },
        report_html => {
            allow       => qr/ \S.html \z/xms,
            required    => 1,
            store       => \$report_html,
            strict_type => 1,
        },
        report_json => {
            allow       => qr/ \S.json \z/xms,
            required    => 1,
            store       => \$report_json,
            strict_type => 1,
        },
        second_infile_path => {
            store       => \$second_infile_path,
            strict_type => 1,
        },
        second_outfile_path => {
            store       => \$second_outfile_path,
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
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        threads => {
            allow       => qr/\A \d+ \z/xms,
            defined     => 1,
            store       => \$threads,
            strict_type => 1,
        },
        trim_poly_g => {
            allow       => [ 0, 1 ],
            store       => \$trim_poly_g,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), );

    push @commands, q{--in1} . $SPACE . $first_infile_path;

    push @commands, q{--out1} . $SPACE . $first_outfile_path;

    push @commands, q{--html} . $SPACE . $report_html;

    push @commands, q{--json} . $SPACE . $report_json;

    if ($detect_pe_adapter) {

        push @commands, q{--detect_adapter_for_pe};
    }

    if ($interleaved_in) {

        push @commands, q{--interleaved_in};
    }

    if ($length_required) {

        push @commands, q{--length_required} . $SPACE . $length_required;
    }

    if ($low_complexity_filter) {

        push @commands, q{--low_complexity_filter};
    }

    if ($overrepresentation_analysis) {

        push @commands, q{--overrepresentation_analysis};
    }

    if ($second_infile_path) {

        push @commands, q{--in2} . $SPACE . $second_infile_path;
    }

    if ($second_outfile_path) {

        push @commands, q{--out2} . $SPACE . $second_outfile_path;
    }

    if ($threads) {

        push @commands, q{--thread} . $SPACE . $threads;
    }

    if ($trim_poly_g) {

        push @commands, q{--trim_poly_g};
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
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

1;

