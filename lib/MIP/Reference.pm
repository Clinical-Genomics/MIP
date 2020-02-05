package MIP::Reference;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catfile };
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
use MIP::Constants qw{ $COMMA $LOG_NAME $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_human_genome_file_endings
      get_dict_contigs
      write_contigs_size_file
    };
}

sub check_human_genome_file_endings {

## Function : Check the existance of associated human genome files
## Returns  :
## Arguments: $human_genome_reference_file_endings_ref => Human genome reference file endings
##          : $human_genome_reference_path             => Human genome reference file path
##          : $parameter_href                          => Parameter hash {REF}
##          : $parameter_name                          => The parameter under evaluation

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $human_genome_reference_file_endings_ref;
    my $human_genome_reference_path;
    my $parameter_href;
    my $parameter_name;

    my $tmpl = {
        human_genome_reference_file_endings_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_file_endings_ref,
            strict_type => 1,
        },
        human_genome_reference_path => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_path,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => { store => \$parameter_name, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_filesystem_objects_existance };
    use MIP::Parameter qw{ set_parameter_build_file_status };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Count the number of files that exists
    my $existence_check_counter = 0;

  FILE_ENDING:
    foreach my $file_ending ( @{$human_genome_reference_file_endings_ref} ) {

        ## Duplicate to keep fresh with each loop
        my $path = $human_genome_reference_path;

        ## Dict requires no fastq(.gz) ending
        if ( $file_ending eq q{.dict} ) {

            ## Removes ".file_ending" in filename.FILENDING(.gz)
            my ( $file, $dir_path ) =
              fileparse( $path, qr/ [.]fasta | [.]fasta[.]gz /sxm );
            $path = catfile( $dir_path, $file );
        }

        ## Add current ending
        $path = $path . $file_ending;

        my ($does_exist) = check_filesystem_objects_existance(
            {
                object_name    => $path,
                object_type    => q{file},
                parameter_name => $parameter_name,
            }
        );

        ## Sum up the number of file that exists
        $existence_check_counter = $existence_check_counter + $does_exist;
    }
    ## Files need to be built
    if ( $existence_check_counter != scalar @{$human_genome_reference_file_endings_ref} )
    {

        set_parameter_build_file_status {
            (
                parameter_href => $parameter_href,
                parameter_name => $parameter_name,
                status         => 1,
            )
        };
        return;
    }

    # All files exist in this check
    set_parameter_build_file_status {
        (
            parameter_href => $parameter_href,
            parameter_name => $parameter_name,
            status         => 0,
        )
    };
    return;
}

sub get_dict_contigs {

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
            oneliner_name  => q{get_dict_contigs},
            stdinfile_path => $dict_file_path,
        }
    );

    my @get_dict_contigs_cmds = join $SPACE, ( @perl_commands, );

    # System call
    my (
        $success_ref,    $error_message_ref, $full_buf_ref,
        $stdout_buf_ref, $stderr_buf_ref
    ) = run( command => \@get_dict_contigs_cmds, verbose => 0 );

    # Save contigs
    my @contigs = split $COMMA, join $COMMA, @{$stdout_buf_ref};

    return @contigs if (@contigs);

    $log->fatal(
        q{Could not detect any 'SN:contig_names' in dict file: } . $dict_file_path );
    exit 1;
}

sub write_contigs_size_file {

## Function : Write contig size file from human genome sequence fai (.fai) file
## Returns  :
## Arguments: $fai_file_path => Fai file path
##          : $outfile_path  => Chromosome size file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fai_file_path;
    my $outfile_path;

    my $tmpl = {
        fai_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$fai_file_path,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
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
            oneliner_name   => q{write_contigs_size_file},
            stdinfile_path  => $fai_file_path,
            stdoutfile_path => $outfile_path,
        }
    );

    my @write_contigs_size_cmd = join $SPACE, ( @perl_commands, );

    # System call
    my (
        $success_ref,    $error_message_ref, $full_buf_ref,
        $stdout_buf_ref, $stderr_buf_ref
    ) = run( command => \@write_contigs_size_cmd, verbose => 0 );

    return if ( not @{$stderr_buf_ref} );

    $log->fatal(q{Could not write contigs size file});
    $log->fatal( q{Error: } . join $NEWLINE, @{$stderr_buf_ref} );
    exit 1;
}

1;
