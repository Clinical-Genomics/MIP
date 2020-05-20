package MIP::File::Format::Star_fusion;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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
use MIP::Constants qw{ $NEWLINE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ create_star_fusion_sample_file };
}

## Constants
Readonly my $TAB => q{\t};

sub create_star_fusion_sample_file {

## Function : Create the samples file for STAR-fusion.
## Returns  :
## Arguments: $filehandle        => Filehandle to write to
##          : $file_info_href    => File info hash {REF}
##          : $infile_paths_ref  => Infile paths for sample {REF}
##          : $sample_id         => Sample id
##          : $samples_file_path => Path to STAR-Fuison sample file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $file_info_href;
    my $infile_paths_ref;
    my $sample_id;
    my $samples_file_path;

    my $tmpl = {
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        samples_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$samples_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_sample_file_attribute };
    use MIP::Program::Gnu::Coreutils qw{ gnu_echo };

    say {$filehandle} q{# Generating STAR-fusion 'samples_file'};

    my %sample_line;
    my @strings;

    # Too adjust infile_index for paired-ends
    my $paired_end_tracker = 0;

    my %file_info_sample = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            sample_id      => $sample_id,
        }
    );

    ## Perform per single-end or read pair
  INFILE_PREFIX:
    while ( my ( $infile_prefix, $sequence_run_type ) =
        each %{ $file_info_sample{file_prefix_no_direction} } )
    {

        ## Add sample id file index array
        push @{ $sample_line{$infile_prefix} }, $sample_id;

        ## Add read one to file index array
        push @{ $sample_line{$infile_prefix} }, $infile_paths_ref->[$paired_end_tracker];

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            # Increment to collect correct read 2 from infiles
            $paired_end_tracker++;

            ## Add read two file index array
            push @{ $sample_line{$infile_prefix} },
              $infile_paths_ref->[$paired_end_tracker];
        }

        ## Increment paired end tracker
        $paired_end_tracker++;

        ## Add tab to each element and add as string to array
        push @strings, join $TAB, @{ $sample_line{$infile_prefix} };
    }

    ## Add newline to each string (line)
    @strings = map { $_ . $NEWLINE } @strings;

    gnu_echo(
        {
            enable_interpretation => 1,
            filehandle            => $filehandle,
            no_trailing_newline   => 1,
            outfile_path          => $samples_file_path,
            strings_ref           => \@strings,
        }
    );
    return 1;
}

1;
