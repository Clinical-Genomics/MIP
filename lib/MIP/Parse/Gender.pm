package MIP::Parse::Gender;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $BACKWARD_SLASH
  $DASH
  $DOT
  $DOUBLE_QUOTE
  $EQUALS
  $LOG_NAME
  $PIPE
};

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      parse_fastq_for_gender
      update_gender_info
    };
}

sub parse_fastq_for_gender {

## Function : Parse fastq infiles for gender. Update contigs depending on results.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $consensus_analysis_type => Consensus analysis_type
##          : $file_info_href          => File info hash {REF}
##          : $sample_info_href        => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $consensus_analysis_type;
    my $file_info_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_active_parameter_attribute };
    use MIP::Recipes::Analysis::Estimate_gender qw{ get_number_of_male_reads };
    use MIP::Fastq qw{ get_stream_fastq_file_cmd };
    use MIP::File_info qw{ get_sample_file_attribute get_sampling_fastq_files };
    use MIP::Program::Gnu::Coreutils qw{ gnu_cut };
    use MIP::Program::Gnu::Software::Gnu_grep qw{ gnu_grep };
    use MIP::Program::Bwa qw{ bwa_mem2_mem };

    my $other_gender_count = get_active_parameter_attribute(
        {
            active_parameter_href => $active_parameter_href,
            attribute             => q{others},
            parameter_name        => q{gender},
        }
    );
    ## All sample ids have a gender - no need to continue
    return if ( not $other_gender_count );

    ## Unpack
    my $log                = Log::Log4perl->get_logger($LOG_NAME);
    my $referencefile_path = $active_parameter_href->{human_genome_reference};

  SAMPLE_ID:
    for my $sample_id ( @{ $active_parameter_href->{gender}{others} } ) {

        next SAMPLE_ID
          if ( not $active_parameter_href->{analysis_type}{$sample_id} eq q{wgs} );

        $log->warn(qq{Detected gender "other/unknown" for sample_id: $sample_id});
        $log->warn(q{Sampling reads from fastq file to estimate gender});

        ### Estimate gender from reads

        my %file_info_sample = get_sample_file_attribute(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );

        my ( $is_interleaved_fastq, @fastq_files ) = get_sampling_fastq_files(
            {
                file_info_sample_href => \%file_info_sample,
            }
        );

        ## Add infile dir to fastq files
        @fastq_files = map { catfile( $file_info_sample{mip_infiles_dir}, $_ ) } @fastq_files;

        ## Build command for streaming of chunk from fastq file(s)
        my @bwa_infiles = get_stream_fastq_file_cmd( { fastq_files_ref => \@fastq_files, } );

        ## Make reference dir available
        my @commands = ( q{SINGULARITY_BIND} . $EQUALS . $active_parameter_href->{reference_dir} );

        ## Build bwa mem command
        push @commands,
          bwa_mem2_mem(
            {
                idxbase                 => $referencefile_path,
                infile_path             => $bwa_infiles[0],
                interleaved_fastq_file  => $is_interleaved_fastq,
                mark_split_as_secondary => 1,
                second_infile_path      => $bwa_infiles[1],
                thread_number           => 2,
            }
          );
        push @commands, $PIPE;

        ## Cut column with contig name
        push @commands,
          gnu_cut(
            {
                infile_path => $DASH,
                list        => 3,
            }
          );
        push @commands, $PIPE;

        ## Count contig name for male contig
        push @commands,
          gnu_grep(
            {
                count   => 1,
                pattern => $DOUBLE_QUOTE . q{chrY} . $BACKWARD_SLASH . $PIPE . q{Y} . $DOUBLE_QUOTE,
            }
          );

        ## Get the number of male reads by aligning fastq read chunk and counting "chrY" or "Y" aligned reads
        my $y_read_count = get_number_of_male_reads( { commands_ref => \@commands, } );

        ## Update gender info in active_parameter and update contigs depending on results
        update_gender_info(
            {
                active_parameter_href   => $active_parameter_href,
                consensus_analysis_type => $consensus_analysis_type,
                file_info_href          => $file_info_href,
                sample_id               => $sample_id,
                sample_info_href        => $sample_info_href,
                y_read_count            => $y_read_count,
            }
        );
    }
    return 1;
}

sub update_gender_info {

## Function : Update gender info in active_parameter and update contigs depending on results.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $consensus_analysis_type => Consensus analysis_type
##          : $file_info_href          => File info hash {REF}
##          : $sample_id               => Sample id
##          : $sample_info_href        => File info hash {REF}
##          : $y_read_count            => Y read count

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $consensus_analysis_type;
    my $file_info_href;
    my $sample_id;
    my $sample_info_href;
    my $y_read_count;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        y_read_count => {
            required => 1,
            store    => \$y_read_count,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter
      qw{ add_gender remove_sample_id_from_gender set_gender_estimation set_include_y };
    use MIP::Contigs qw{ update_contigs_for_run };
    use MIP::Sample_info qw{ set_sample_gender };

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Constants
    Readonly my $MALE_THRESHOLD => 36;

    my $gender  = $y_read_count > $MALE_THRESHOLD ? q{male} : q{female};
    my $genders = $gender . q{s};
    $log->info(qq{Found $gender according to fastq reads});

    ## Update in active parameter hash
    add_gender(
        {
            active_parameter_href => $active_parameter_href,
            sample_id             => $sample_id,
            gender                => $genders,
        }
    );

    ## For tracability
    set_gender_estimation(
        {
            active_parameter_href => $active_parameter_href,
            gender                => $gender,
            sample_id             => $sample_id,
        }
    );

    remove_sample_id_from_gender(
        {
            active_parameter_href => $active_parameter_href,
            gender                => q{others},
            sample_id             => $sample_id,
        }
    );

    set_include_y(
        {
            active_parameter_href => $active_parameter_href,
        }
    );

    ## Update gender in sample info hash
    set_sample_gender(
        {
            gender           => $gender,
            sample_id        => $sample_id,
            sample_info_href => $sample_info_href,
        }
    );

    ## Update cache

    ## Update contigs depending on settings in run (wes or if only male samples)
    update_contigs_for_run(
        {
            consensus_analysis_type => $consensus_analysis_type,
            exclude_contigs_ref     => \@{ $active_parameter_href->{exclude_contigs} },
            file_info_href          => $file_info_href,
            include_y               => $active_parameter_href->{include_y},
        }
    );
    return 1;
}

1;
