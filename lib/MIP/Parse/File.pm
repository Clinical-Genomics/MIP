package MIP::Parse::File;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $EMPTY_STR $LOG_NAME $SPACE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_io_outfiles };

}

sub parse_io_outfiles {

## Function : Set and get the io files per chain, id and stream
## Returns  : %io
## Arguments: $chain_id               => Chain of recipe
##          : $file_info_href         => File info hash {REF}
##          : $file_name_prefixes     => Build outfile using file name prefix
##          : $file_name_prefixes_ref => Build outfile using file name prefixes {REF}
##          : $file_paths_ref         => File paths {REF}
##          : $id                     => Id (sample or case)
##          : $iterators_ref          => Build outfile using iterator (e.g contigs) {REF}
##          : $outdata_dir            => Outdata directory
##          : $parameter_href         => Parameter hash {REF}
##          : $recipe_name            => Recipe name
##          : $stream                 => Stream (out)
##          : $temp_directory         => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_id;
    my $file_info_href;
    my $file_name_prefix;
    my $file_name_prefixes_ref;
    my $file_paths_ref;
    my $id;
    my $iterators_ref;
    my $outdata_dir;
    my $parameter_href;
    my $recipe_name;
    my $temp_directory;

    ## Default(s)
    my $stream;

    my $tmpl = {
        chain_id => {
            defined     => 1,
            required    => 1,
            store       => \$chain_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_name_prefix => {
            store       => \$file_name_prefix,
            strict_type => 1,
        },
        file_name_prefixes_ref => {
            default     => [],
            store       => \$file_name_prefixes_ref,
            strict_type => 1,
        },
        file_paths_ref => {
            default     => [],
            store       => \$file_paths_ref,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        iterators_ref => {
            default     => [],
            store       => \$iterators_ref,
            strict_type => 1,
        },
        outdata_dir => {
            store       => \$outdata_dir,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        stream => {
            allow       => [qw{ out }],
            default     => q{out},
            store       => \$stream,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Parameter qw{ get_recipe_attributes };
    use MIP::Set::File qw{ set_io_files };

    my @file_paths = @{$file_paths_ref};

    ## Build default @file_paths
    if ( not @file_paths and $outdata_dir ) {

        my %rec_atr = get_recipe_attributes(
            {
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
            }
        );
        my $outfile_tag    = $rec_atr{file_tag}       //= $EMPTY_STR;
        my $outfile_suffix = $rec_atr{outfile_suffix} //= $EMPTY_STR;
        my $directory      = catdir( $outdata_dir, $id, $recipe_name );

        ## Default paths with iterators
        if ( @{$iterators_ref} and $file_name_prefix ) {

            ## Localize as we will mutate elements
            my @iterators = @{$iterators_ref};
            foreach my $iterator (@iterators) {

                ## Add "." if not empty string
                if ( $iterator ne $EMPTY_STR ) {

                    $iterator = $DOT . $iterator;
                }
            }
            @file_paths =
              map { catfile( $directory, $file_name_prefix . $outfile_tag . $_ . $outfile_suffix ) }
              @iterators;
        }
        ## Default paths without iterators
        else {

            ## $file_name_prefixes_needs to be set
            croak q{Missing argument!} if not @{$file_name_prefixes_ref};
            @file_paths =
              map { catfile( $directory, $_ . $outfile_tag . $outfile_suffix ) }
              @{$file_name_prefixes_ref};
        }
    }

    ## Set the io files per chain and stream
    set_io_files(
        {
            chain_id       => $chain_id,
            id             => $id,
            file_paths_ref => \@file_paths,
            file_info_href => $file_info_href,
            recipe_name    => $recipe_name,
            stream         => $stream,
            temp_directory => $temp_directory,
        }
    );

    my %io = get_io_files(
        {
            id             => $id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => $stream,
        }
    );

    return %io;
}

1;
