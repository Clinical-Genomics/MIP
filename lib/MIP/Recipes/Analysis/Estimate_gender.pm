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
    our @EXPORT_OK = qw{ get_number_of_male_reads };
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
    $log->info(qq{Number of read from chromosome Y: $y_read_count});

    ## Clean-up
    close $filehandle;

    return $y_read_count;
}

1;
