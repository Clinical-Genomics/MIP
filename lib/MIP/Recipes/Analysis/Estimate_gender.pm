package MIP::Recipes::Analysis::Estimate_gender;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Path qw{ make_path };
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
use MIP::Constants qw{ $COLON $LOG_NAME $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_number_of_male_reads update_gender_info };
}

sub get_number_of_male_reads {

## Function : Get the number of male reads by aligning fastq read chunk and counting "chrY" or "Y" aligned reads
## Returns  : $y_read_count
## Arguments: $commands_ref  => Command array for estimating number of male reads {REF}
##          : $outscript_dir => Directory to write the bash script to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $outscript_dir;

    my $tmpl = {
        commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$commands_ref,
            strict_type => 1,
        },
        outscript_dir => => {
            defined     => 1,
            required    => 1,
            store       => \$outscript_dir,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    make_path($outscript_dir);

    ## Recipe bash file for commands to estimate gender
    my $bash_file = catfile( $outscript_dir, q{estimate_gender_from_reads} . q{.sh} );

    open my $filehandle, q{>}, $bash_file
      or croak q{Cannot write to} . $SPACE . $bash_file . $COLON . $SPACE . $OS_ERROR;

    $log->info(qq{Writing estimation of gender recipe to: $bash_file});

    ## Write to file
    say {$filehandle} join $SPACE, @{$commands_ref};

    my $cmds_ref       = [ q{bash}, $bash_file ];
    my %process_return = child_process(
        {
            commands_ref => $cmds_ref,
            process_type => q{ipc_cmd_run},
        }
    );

    my $y_read_count = $process_return{stdouts_ref}->[0];

    ## Clean-up
    close $filehandle;

    return $y_read_count;
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

    $log->info(qq{Requiring greater than $MALE_THRESHOLD reads from chromosome Y to set as male});
    $log->info(qq{Number of reads from chromosome Y: $y_read_count});
    $log->info(qq{Found $gender according to fastq reads});
    $log->info(qq{Updating sample: $sample_id to gender: $gender for analysis});

    ## Update in active parameter hash
    add_gender(
        {
            active_parameter_href => $active_parameter_href,
            gender                => $genders,
            sample_id             => $sample_id,
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
