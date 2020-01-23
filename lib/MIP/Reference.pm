package MIP::Reference;

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
use MIP::Constants qw{ $COMMA $LOG_NAME $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_seq_dict_contigs };
}

sub get_seq_dict_contigs {

## Function : Collects sequence contigs used in analysis from human genome sequence
##          : dictionnary (.dict file)
## Returns  : @contigs
## Arguments: $dict_file_path => Dict file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dict_file_path;

    my $tmpl = {
        dict_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$dict_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use IPC::Cmd qw{ run };
    use MIP::Language::Perl qw{ perl_nae_oneliners };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Build regexp to find contig names
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_name => q{get_seq_dict_contigs},
        }
    );

    my @get_seq_dict_contigs_cmds = join $SPACE, ( @perl_commands, $dict_file_path );

    # System call
    my (
        $success_ref,    $error_message_ref, $full_buf_ref,
        $stdout_buf_ref, $stderr_buf_ref
    ) = run( command => \@get_seq_dict_contigs_cmds, verbose => 0 );

    # Save contigs
    my @contigs = split $COMMA, join $COMMA, @{$stdout_buf_ref};

    return @contigs if (@contigs);

    $log->fatal(
        q{Could not detect any 'SN:contig_names' in dict file: } . $dict_file_path );
    exit 1;
}

1;
