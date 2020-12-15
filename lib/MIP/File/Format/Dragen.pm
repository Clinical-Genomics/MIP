package MIP::File::Format::Dragen;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ make_path };
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
use MIP::Constants qw{ $COMMA $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ create_dragen_fastq_list_sample_id };
}

sub create_dragen_fastq_list_sample_id {

## Function : Create dragen fastq list sample id file for supplying input fastq files to dragen
## Returns  :
## Arguments: $fastq_list_lines_ref => Infile info per fastq (pair)
##          : $fastq_list_file_path => Where to write the dragen fastq list file
##          : $include_header       => Include header ("1") or not ("0")
##          : $log                  => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fastq_list_lines_ref;
    my $fastq_list_file_path;
    my $log;

    ## Default(s)
    my $include_header;

    my $tmpl = {
        fastq_list_lines_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$fastq_list_lines_ref,
            strict_type => 1,
        },
        fastq_list_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$fastq_list_file_path,
            strict_type => 1,
        },
        include_header => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$include_header,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @fastq_list_headers = qw{ RGID RGSM RGLB Lane Read1File Read2File };

    # Create anonymous filehandle
    my $filehandle_SYS = IO::Handle->new();

    ## Create dir if it does not exists
    make_path( dirname($fastq_list_file_path) );

    open $filehandle_SYS, q{>}, $fastq_list_file_path
      or $log->logdie(qq{Can't open $fastq_list_file_path: $ERRNO });

    ## Adds the information from the samples in pedigree_lines, separated by \n
    ## Add @fastq_list_headers
    if ($include_header) {

        say {$filehandle_SYS} join $COMMA, @fastq_list_headers;
    }
  LINE:
    foreach my $line ( @{$fastq_list_lines_ref} ) {

        say {$filehandle_SYS} $line;
    }
    $log->info( q{Wrote: } . $fastq_list_file_path, $NEWLINE );
    close $filehandle_SYS;
    return;
}

1;
