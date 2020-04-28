package MIP::Fastq;

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

## MIPs lib/
use MIP::Constants qw{ $COMMA $LOG_NAME $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_fastq_infiles_format };
}

sub parse_fastq_infiles_format {

## Function : Parse infile according to MIP filename convention
## Returns  : %infile_info or undef
## Arguments: $file_name => File name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;

    my $tmpl = {
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    # Store fastq file name features
    my %fastq_file_name;
    my @missing_feature;

    ## Define MIP fastq file name formats matching regexp
    my %fastq_file_name_regexp = _fastq_file_name_regexp();

    # Parse fastq file name
    my @file_features = $file_name =~ /$fastq_file_name_regexp{regexp}/sxm;

  FEATURE:
    while ( my ( $index, $feature ) = each @{ $fastq_file_name_regexp{features} } ) {

        ## Return undef if not all expected features found
        if ( not $file_features[$index] ) {
            push @missing_feature, $feature;
        }

        # Store feature
        $fastq_file_name{$feature} = $file_features[$index];
    }
    if (@missing_feature) {
        $log->warn(qq{Could not detect MIP file name convention for file: $file_name });
        $log->warn( q{Missing file name feature: } . join $COMMA . $SPACE,
            @missing_feature );
        return;
    }
    return %fastq_file_name;
}

sub _fastq_file_name_regexp {

## Function : Define MIP fastq file name formats matching regexp
## Returns  : %fastq_file_name_regexp
## Arguments:

    my ($arg_href) = @_;

    my %fastq_file_name_regexp = (
        features => [qw{ lane date flowcell infile_sample_id index direction }],
        regexp   => q?(\d+)_(\d+)_([^_]+)_([^_]+)_([^_]+)_(\d).fastq?,
    );

    return %fastq_file_name_regexp;
}

1;