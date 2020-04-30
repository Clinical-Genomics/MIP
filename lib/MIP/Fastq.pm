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
use MIP::Constants qw{ $COMMA $LOG_NAME $SPACE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ define_mip_fastq_file_features get_read_length parse_fastq_infiles_format };
}

sub define_mip_fastq_file_features {

## Function : Define MIP internal fastq file features derived from supplied fastq infile name
## Returns  : $mip_file_format, $mip_file_format_with_direction, $original_file_name_prefix, $run_barcode
## Arguments: $date               => Flow-cell sequencing date
##          : $direction          => Sequencing read direction
##          : $flowcell           => Flow-cell id
##          : $index              => The DNA library preparation molecular barcode
##          : $lane               => Flow-cell lane
##          : $original_file_name => Original file name
##          : $sample_id          => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $date;
    my $direction;
    my $flowcell;
    my $index;
    my $lane;
    my $original_file_name;
    my $sample_id;

    my $tmpl = {
        date => {
            defined     => 1,
            required    => 1,
            store       => \$date,
            strict_type => 1,
        },
        direction => {
            allow       => [ 1, 2 ],
            defined     => 1,
            required    => 1,
            store       => \$direction,
            strict_type => 1,
        },
        flowcell => {
            defined     => 1,
            required    => 1,
            store       => \$flowcell,
            strict_type => 1,
        },
        index => { defined => 1, required => 1, store => \$index, strict_type => 1, },
        lane  => {
            allow       => qr{ \A\d+\z }xsm,
            defined     => 1,
            required    => 1,
            store       => \$lane,
            strict_type => 1,
        },
        original_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$original_file_name,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $mip_file_format = join $UNDERSCORE,
      ( $sample_id, $date, $flowcell, $index, q{lane} . $lane );

    my $mip_file_format_with_direction = join $UNDERSCORE,
      ( $mip_file_format, $direction );

    my $original_file_name_prefix = substr $original_file_name, 0,
      index $original_file_name, q{.fastq};

    my $run_barcode = join $UNDERSCORE, ( $date, $flowcell, $lane, $index );

    return $mip_file_format, $mip_file_format_with_direction,
      $original_file_name_prefix, $run_barcode;
}

sub get_read_length {

## Function : Collect read length from a fastq infile
## Returns  : $read_length
## Arguments: $file_path => File to parse
##          : $read_file => Command used to read file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $read_file_command;

    my $tmpl = {
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        read_file_command => {
            defined     => 1,
            required    => 1,
            store       => \$read_file_command,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Language::Perl qw{ perl_nae_oneliners };

    ## Build regexp to find read length
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_name => q{get_fastq_read_length},
        }
    );

    my $read_length_cmd = qq{$read_file_command $file_path | @perl_commands};

    my %process_return = child_process(
        {
            commands_ref => [$read_length_cmd],
            process_type => q{ipc_cmd_run},
        }
    );

    ## Return read length
    return $process_return{stdouts_ref}[0];
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
