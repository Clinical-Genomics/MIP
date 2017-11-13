package MIP::Get::File;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use Readonly;
use List::MoreUtils qw{ any };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ get_exom_target_bed_file get_file_suffix get_merged_infile_prefix get_select_file_contigs get_seq_dict_contigs };
}

## Constants
Readonly my $COMMA   => q{,};
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub get_exom_target_bed_file {

## Function : Get exome_target_bed file for specfic sample_id and add file_ending from file_info hash if supplied
## Returns  : $exome_target_bed_file
## Arguments: $exome_target_bed_href => Exome target bed files lnked to sample ids
##          : $sample_id             => Sample id
##          : $log                   => Log object
##          : $file_ending           => File ending to add to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exome_target_bed_href;
    my $sample_id;
    my $log;
    my $file_ending;

    my $tmpl = {
        exome_target_bed_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$exome_target_bed_href,
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
        file_ending => { store => \$file_ending },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## To track sample_ids with capture kits
    my %seen;

    ## Make local hash copy to keep in scope
    my %exome_target_bed = %{$exome_target_bed_href};

  BED_FILE:
    while ( my ( $exome_target_bed_file, $sample_id_string ) =
        each %exome_target_bed )
    {

        my @capture_kit_samples = split $COMMA, $sample_id_string;

        ## Count number of times sample_id has been seen
      SAMPLE:
        foreach my $samples (@capture_kit_samples) {

            $seen{$samples}++;
        }

        ## If capture_kit sample_id is associated with exome_target_bed file
        if ( any { $_ eq $sample_id } @capture_kit_samples ) {

            if ( defined $file_ending ) {

                $exome_target_bed_file .= $file_ending;
            }
            return $exome_target_bed_file;
        }
    }
    if ( not defined $seen{$sample_id} ) {

        $log->fatal(
            q{Could not detect }
              . $sample_id
              . q{ in '-exome_target_bed' associated files in sub routine get_exom_target_bed_file},
            $NEWLINE
        );
        exit 1;
    }
    return;
}

sub get_file_suffix {

## Function : Return the current file suffix for this jobid chain or program
## Returns  : $file_suffix
## Arguments: $parameter_href => Holds all parameters
##          : $suffix_key     => Suffix key
##          : $job_id_chain   => Job id chain for program
##          : $program_name   => Program name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $suffix_key;
    my $job_id_chain;
    my $program_name;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        suffix_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$suffix_key
        },
        jobid_chain  => { strict_type => 1, store => \$job_id_chain },
        program_name => { strict_type => 1, store => \$program_name },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $file_suffix;

    ## Jobid chain specific suffix
    if ( defined $job_id_chain ) {

        $file_suffix = $parameter_href->{$suffix_key}{$job_id_chain};
    }
    elsif ( defined $program_name ) {
        ## Program  specific

        $file_suffix = $parameter_href->{$program_name}{$suffix_key};
    }

    ## If suffix was found
    if ( defined $file_suffix && $file_suffix ) {

        return $file_suffix;
    }
    else {
        ## Broadcast no suffix was found
        if ( defined $job_id_chain ) {

            say $log->fatal(
                q{Could not get requested infile_suffix for jobid_chain:}
                  . $job_id_chain );
        }
        elsif ( defined $program_name ) {

            say $log->fatal(
                q{Could not get requested infile_suffix for program:}
                  . $program_name );
        }
        exit 1;
    }
    return;
}

sub get_merged_infile_prefix {

## Function : Get the merged infile prefix for sample id
## Returns  : $merged_infile_prefix
## Arguments: $file_info_href => File info hash {REF}
##          : $sample_id      => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return $file_info_href->{$sample_id}{merged_infile};
}

sub get_select_file_contigs {

## Function : Collects sequences contigs used in select file
## Returns  :
## Arguments: $select_file_path => Select file path
##          : $log              => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $select_file_path;
    my $log;

    my $tmpl = {
        select_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$select_file_path
        },
        log => { store => \$log },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Execute perl
    my $find_contig_name = q?perl -nae ?;

    ## Get contig name
    $find_contig_name .= q?'if ($_=~/ contig=(\w+) /xsm) { ?;

    ## Alias capture
    $find_contig_name .= q?my $contig_name = $1; ?;

    ## Write contig name and comma
    $find_contig_name .= q?print $contig_name, q{,};} ?;

    ## Quit if #CHROM found in line
    $find_contig_name .= q?if($_=~/ [#]CHROM /xsm) {last;}' ?;

    ## Returns a comma seperated string of sequence contigs from file
    my @contigs = `$find_contig_name $select_file_path `;

    @contigs = split $COMMA, join $COMMA, @contigs;

    my $error_msg;

    if ( not @contigs ) {

        $error_msg =
q{Could not detect any '##contig' in meta data header in select file: }
          . $select_file_path;

        return $error_msg;
    }
    return $error_msg, @contigs;
}

sub get_seq_dict_contigs {

## Function : Collects sequences contigs used in analysis from human genome sequence dictionnary associated with $human_genome_reference
## Returns  :
## Arguments: $dict_file_path => Dict file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dict_file_path;

    my $tmpl = {
        dict_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$dict_file_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### Build regexp to find contig names
    ## Execute perl
    my $find_contig_name = q?perl -nae '?;

    ## Find contig line
    $find_contig_name .= q?if($F[0]=~/^\@SQ/) { ?;

    ## Collect contig name
    $find_contig_name .= q? if($F[1]=~/SN\:(\S+)/) { ?;

    ## Alias capture
    $find_contig_name .= q?my $contig_name = $1; ?;

    ## Write to STDOUT
    $find_contig_name .= q?print $contig_name, q{,};} }' ?;

    ## Returns a comma seperated string of sequence contigs from dict file
    my @contigs = `$find_contig_name $dict_file_path `;

    ## Save contigs
    @contigs = split $COMMA, join $COMMA, @contigs;

    my $error_msg;

    if ( not @contigs ) {

        $error_msg = q{Could not detect any 'SN:contig_names' in dict file: }
          . $dict_file_path;
        return $error_msg;
    }
    return $error_msg, @contigs;
}

1;
