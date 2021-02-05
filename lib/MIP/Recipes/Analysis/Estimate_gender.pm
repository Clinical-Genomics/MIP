package MIP::Recipes::Analysis::Estimate_gender;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use FindBin qw{ $Bin };
use File::Path qw{ remove_tree };
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
use MIP::Constants qw{ $COLON $LOG_NAME $SPACE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_number_of_male_reads };
}

sub get_number_of_male_reads {

## Function : Get the number of male reads by aligning fastq read chunk and counting "chrY" or "Y" aligned reads
## Returns  : $y_read_count
## Arguments: $commands_ref => Command array for estimating number of male reads {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;

    my $tmpl = {
        commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$commands_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };

    ## Constants
    Readonly my $MAX_RANDOM_NUMBER => 10_000;

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Generate a random integer between 0-10,000.
    my $random_integer = int rand $MAX_RANDOM_NUMBER;

    ## Temporary bash file for commands
    my $bash_temp_file =
      catfile( $Bin, q{estimate_gender_from_reads} . $UNDERSCORE . $random_integer . q{.sh} );

    open my $filehandle, q{>}, $bash_temp_file
      or croak q{Cannot write to} . $SPACE . $bash_temp_file . $COLON . $SPACE . $OS_ERROR;

    ## Write to file
    say {$filehandle} join $SPACE, @{$commands_ref};

    my $cmds_ref       = [ q{bash}, $bash_temp_file ];
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
    remove_tree($bash_temp_file);

    return $y_read_count;
}

1;
