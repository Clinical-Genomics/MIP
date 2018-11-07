package MIP::Main::Download;

use 5.026;
use Carp;
use charnames qw( :full :short );
use Cwd;
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catfile catdir devnull };
use FindBin qw{ $Bin };
use Getopt::Long;
use IO::Handle;
use warnings qw{ FATAL utf8 };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;

## CPANM
use autodie qw{ open close :all };
use List::Util qw{ any };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Gnu::Coreutils qw{ gnu_mkdir };
use MIP::Language::Shell qw{ create_bash_file };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Script::Utils qw{ create_temp_dir };

## Constants
Readonly my $COLON        => q{:};
Readonly my $DOT          => q{.};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $SPACE        => q{ };
Readonly my $UNDERSCORE   => q{_};

BEGIN {
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = q{0.0.4};

    # Functions and variables that can be optionally exported
    our @EXPORT_OK = qw{ mip_download };

}

sub mip_download {

## Function : Main script for generating MIP download scripts
## Returns  :
## Arguments: $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments};

    ## Transfer to lexical variables
    my %parameter = %{$parameter_href};

    ## Create default log name
    if ( not $parameter{log_file} ) {

        ## Get local time
        my $date_time       = localtime;
        my $date_time_stamp = $date_time->datetime;

        $parameter{log_file} =
          catfile( q{mip_donwload} . $UNDERSCORE . $date_time_stamp . $DOT . q{log} );
    }

    ## Initiate logger
    my $log = initiate_logger(
        {
            file_path => $parameter{log_file},
            log_name  => q{mip_download},
        }
    );
    $log->info( q{Writing log messages to} . $COLON . $SPACE . $parameter{log_file} );

## Set parameter default
    if ( not $parameter{reference_dir} ) {

        $parameter{reference_dir} = cwd();
    }

    our $VERSION = q{0.0.4};

    check_user_reference(
        {
            cmd_reference_ref => \%{ $parameter{cmd_reference} },
            reference_ref     => \%{ $parameter{reference} },
        }
    );

## Change relative path to absolute path for certain parameters
    update_to_absolute_path( { parameter_href => \%parameter, } );

    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    # Downloads instruction file
    my $bash_file_path = catfile( cwd(), q{download_reference} . $DOT . q{sh} );

    open $FILEHANDLE, '>', $bash_file_path
      or
      $log->logdie( q{Cannot write to '} . $bash_file_path . q{' :} . $OS_ERROR . "\n" );

    # Install temp directory
    my $temp_dir = create_temp_dir( {} );

## Create bash file for writing install instructions
    create_bash_file(
        {
            file_name   => $bash_file_path,
            FILEHANDLE  => $FILEHANDLE,
            remove_dir  => $temp_dir,
            log         => $log,
            set_errexit => 1,
            set_nounset => 1,
        }
    );

    say {$FILEHANDLE} q{## Create temp dir};

    gnu_mkdir(
        {
            FILEHANDLE       => $FILEHANDLE,
            indirectory_path => $temp_dir,
            parents          => 1,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    $log->info( q{Will write install instructions to '} . $bash_file_path,
        $SINGLE_QUOTE );

## Build install references recipe in bash file
    build_reference_install_recipe(
        {
            parameter_href => \%parameter,
            FILEHANDLE     => $FILEHANDLE,
        }
    );

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    return;
}

#################
###SubRoutines###
#################

sub build_reference_install_recipe {

## Function : Build install references recipe in bash file
## Returns  :
## Arguments: $parameter_href => Holds all parameters
##          : $FILEHANDLE     => Filehandle to write to
##          : $quiet          => Be quiet
##          : $verbose        => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    ## Default(s)
    my $quiet;
    my $verbose;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        quiet      => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        verbose => {
            default     => $arg_href->{parameter_href}{verbose},
            strict_type => 1,
            store       => \$verbose
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Bash qw{ gnu_cd };
    use MIP::Gnu::Coreutils qw{ gnu_mkdir };

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger(q{Download_reference});

    my $pwd = cwd();

    say {$FILEHANDLE} q{## Create reference directory};
    gnu_mkdir(
        {
            indirectory_path => $parameter_href->{reference_dir},
            parents          => 1,
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Since all commands should assume working directory to be the reference directory
    gnu_cd(
        {
            directory_path => $parameter_href->{reference_dir},
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

  REFERENCE:
    while ( my ( $reference_id, $versions_ref ) = each %{ $parameter_href->{reference} } )
    {

        ## Remodel depending on if "--reference" was used or not as the user info is stored as a scalar per reference_id while yaml is stored as arrays per reference_id

        my @reference_versions;
        if ( ref $versions_ref eq q{ARRAY} ) {

            @reference_versions = @{$versions_ref};
        }
        else {

            push @reference_versions, $versions_ref;
        }

      REFERENCE_VERSION:
        foreach my $reference_version (@reference_versions) {

          GENOME_VERSION:
            foreach
              my $genome_version ( @{ $parameter_href->{reference_genome_versions} } )
            {

                ## Standardize case
                $genome_version = lc $genome_version;

                my $reference_href =
                  $parameter_href->{$reference_id}{$genome_version}{$reference_version};

                next GENOME_VERSION
                  if ( not exists $parameter_href->{$reference_id}{$genome_version} );

                next GENOME_VERSION
                  if (
                    not exists $parameter_href->{$reference_id}{$genome_version}
                    {$reference_version} );

                ## Build file name and path
                my $outfile_name = $reference_href->{outfile};
                my $outfile_path =
                  catfile( $parameter_href->{reference_dir}, $outfile_name );

                ## Check if reference already exists in reference directory
                next GENOME_VERSION if ( -f $outfile_path );

                $log->warn( q{Cannot find reference file:} . $outfile_path );
                $log->warn(
                        q{Will try to download: }
                      . $reference_id
                      . q{ version: }
                      . $reference_version,
                );

                reference_install_recipe(
                    {
                        parameter_href => $parameter_href,
                        reference_href => $reference_href,
                        reference_id   => $reference_id,
                        FILEHANDLE     => $FILEHANDLE,
                        quiet          => $quiet,
                        verbose        => $verbose,
                    }
                );

                ## Writes post processing commands associated with reference e.g. tabix
                write_post_processing_command(
                    {
                        reference_href => $reference_href,
                        FILEHANDLE     => $FILEHANDLE,
                    }
                );
            }
        }
    }

    ## Move back to original
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;
    return;
}

sub reference_install_recipe {

## Function : Write reference install recipe
## Returns  :
## Arguments: $parameter_href => Holds all parameters {REF}
##          : $reference_href => Reference hash {REF}
##          : $reference_id  => Reference id
##          : $FILEHANDLE     => Filehandle to write to
##          : $quiet          => Quiet (no output)
##          : $verbose        => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $reference_href;
    my $reference_id;
    my $FILEHANDLE;

    ## Default(s)
    my $quiet;
    my $verbose;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        reference_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$reference_href,
        },
        reference_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$reference_id
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        quiet      => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        verbose => {
            default     => $arg_href->{parameter_href}{verbose},
            strict_type => 1,
            store       => \$verbose
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Potential download files
    my @file_keys = qw{ file file_check
      file_index file_index_check };

  REFERENCE_FILE:
    foreach my $key (@file_keys) {

        next REFERENCE_FILE
          if ( not exists $reference_href->{$key} );

        ## Install reference
        my $file         = $reference_href->{$key};
        my $outfile      = $reference_href->{ q{out} . $key };
        my $outfile_path = catfile( $parameter_href->{reference_dir}, $outfile );

        download(
            {
                parameter_href => $parameter_href,
                FILEHANDLE     => $FILEHANDLE,
                url            => $reference_href->{url_prefix} . $file,
                outfile_path   => $outfile_path,
                file_id        => $reference_id,
            }
        );

        ## Check if file needs to be decompress and write decompression if so
        decompress_file(
            {
                parameter_href => $parameter_href,
                FILEHANDLE     => $FILEHANDLE,
                outfile_path   => $outfile_path,
                file_decompress =>
                  $reference_href->{ q{out} . $key . $UNDERSCORE . q{decompress} },
            }
        );

        ## Check file integrity of file
        check_file(
            {
                FILEHANDLE         => $FILEHANDLE,
                outfile_path       => $outfile_path,
                outfile_path_check => $outfile_path,
                check_method =>
                  $reference_href->{ q{out} . $key . $UNDERSCORE . q{method} },
            }
        );
    }
    return;
}

sub download {

## Function : Downloads files
## Returns  :
## Arguments: $parameter_href => Holds all parameters
##          : $FILEHANDLE     => Filehandle to write to
##          : $url            => Url to use for download
##          : $outfile_path   => Outfile path
##          : $program        => Program to use for download
##          : $file_id        => File id
##          : $quiet          => Quiet (no output)
##          : $verbose        => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $url;
    my $outfile_path;
    my $file_id;

    ## Default(s)
    my $program;
    my $quiet;
    my $verbose;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        FILEHANDLE => { required => 1, defined => 1, store       => \$FILEHANDLE },
        url        => { required => 1, defined => 1, strict_type => 1, store => \$url },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        file_id => { required => 1, defined => 1, strict_type => 1, store => \$file_id },
        program => {
            default     => q{wget},
            allow       => [qw{ wget }],
            strict_type => 1,
            store       => \$program
        },
        quiet => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        verbose => {
            default     => $arg_href->{parameter_href}{verbose},
            strict_type => 1,
            store       => \$verbose
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Download::Wget qw{ wget };

    ## Download
    say {$FILEHANDLE} q{## Download } . $file_id . $NEWLINE;

    if ( $program eq q{wget} ) {

        wget(
            {
                url          => $url,
                FILEHANDLE   => $FILEHANDLE,
                quiet        => $quiet,
                verbose      => $verbose,
                outfile_path => $outfile_path,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return;
}

sub remove_file_ending {

## Function : Removes ".file_ending" in filename.file_ending(.gz)
## Returns  : File name with supplied $file_ending or $file_ending(.gz) removed
## Arguments: $file_name   => File name
##          : $file_ending => File ending to be removed

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;
    my $file_ending;

    my $tmpl = {
        file_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_name
        },
        file_ending => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_ending
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $file_name_noending;

    if ( defined $file_name
        && $file_name =~ / (\S+)($file_ending$ | $file_ending.gz$) /x )
    {

        $file_name_noending = $1;
    }
    return $file_name_noending;
}

sub decompress_file {

## Function : Check if file needs to be decompress and write decompression if so
## Returns  :
## Arguments: $parameter_href  => Holds all parameters
##          : $FILEHANDLE      => Filehandle to write to
##          : $outfile_path    => Outfile path
##          : $file_decompress => Decompress the downloaded file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $outfile_path;
    my $file_decompress;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        FILEHANDLE   => { required => 1, defined => 1, store => \$FILEHANDLE, },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path,
        },
        file_decompress => { strict_type => 1, store => \$file_decompress, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Compression::Gzip qw{ gzip };
    use MIP::Program::Compression::Tar qw{ tar };
    use MIP::Program::Compression::Zip qw{ unzip };

    return if ( not defined $outfile_path );

    if ( defined $file_decompress && $file_decompress eq q{gzip} ) {

        ## Removes ".file_ending" in filename.FILENDING(.gz)
        my $outfile_path_no_ending = remove_file_ending(
            {
                file_name   => $outfile_path,
                file_ending => $DOT . q{gz},
            }
        );

        gzip(
            {
                infile_path  => $outfile_path,
                outfile_path => $outfile_path_no_ending,
                force        => 1,
                quiet        => 1,
                decompress   => 1,
                stdout       => 1,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    if ( defined $file_decompress && $file_decompress eq q{unzip} ) {

        unzip(
            {
                infile_path => $outfile_path,
                outdir_path => $parameter_href->{reference_dir},
                FILEHANDLE  => $FILEHANDLE,
            }
        );
    }

    if ( defined $file_decompress && $file_decompress eq q{tar} ) {

        tar(
            {
                extract           => 1,
                filter_gzip       => 1,
                file_path         => $outfile_path,
                outdirectory_path => $parameter_href->{reference_dir},
                FILEHANDLE        => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return;
}

sub check_file {

## Function : Check file integrity of file
## Returns  :
## Arguments: $FILEHANDLE         => Filehandle to write to
##          : $outfile_path       => Outfile path
##          : $outfile_path_check => File to check
##          : $check_method       => Method to perform file check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $outfile_path;
    my $outfile_path_check;
    my $check_method;

    my $tmpl = {
        FILEHANDLE   => { required => 1, defined => 1, store => \$FILEHANDLE, },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path,
        },
        outfile_path_check => { strict_type => 1, store => \$outfile_path_check, },
        check_method       => { strict_type => 1, store => \$check_method, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_md5sum gnu_rm };

    return if ( not defined $check_method );

    if ( $check_method eq q{md5sum} ) {

        ## Removes ".file_ending" in filename.FILENDING(.gz)
        my $outfile_path_no_ending = remove_file_ending(
            {
                file_name   => $outfile_path,
                file_ending => $DOT . q{md5},
            }
        );

        return if ( not defined $outfile_path_no_ending );

        my $perl_regexp =
            q?perl -nae 'print $F[0]."  ?
          . $outfile_path_no_ending . q?" ' ?
          . $outfile_path_check;
        print {$FILEHANDLE} $perl_regexp . q{ > md5sum_check.txt};
        say   {$FILEHANDLE} $NEWLINE;

        gnu_md5sum(
            {
                check       => 1,
                infile_path => q{md5sum_check.txt},
                FILEHANDLE  => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Clean-up
        gnu_rm(
            {
                infile_path => q{md5sum_check.txt},
                force       => 1,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return;
}

sub update_to_absolute_path {

## Function : Change relative path to absolute path for certain parameter_names
## Returns  :
## Arguments: $parameter_href => The parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Set::File qw{ set_absolute_path };

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger(q{mip_download});

  PARAMETER:
    foreach my $parameter_name ( @{ $parameter_href->{absolute_paths} } ) {

        next PARAMETER if ( not $parameter_href->{$parameter_name} );

        ## Alias
        my $parameter = \$parameter_href->{$parameter_name};

        ## Array reference
        if ( ref($parameter) eq q{ARRAY} ) {

            foreach my $parameter_value ( @{$parameter} ) {
                ## Replace original input with abolute path for supplied path or croaks and exists if path does not exists
                $parameter_value = set_absolute_path(
                    {
                        path           => $parameter_value,
                        parameter_name => $parameter_name,
                        log            => $log,
                    }
                );
            }
        }
        elsif ( ref($parameter) eq 'HASH' ) {
            ## Hash reference

            foreach my $key ( keys %{ $parameter_href->{$parameter_name} } )
            {    #Cannot use each since we are updating key

                ## Find aboslute path for supplied path or croaks and exists if path does not exists
                my $updated_key = set_absolute_path(
                    {
                        path           => $key,
                        parameter_name => $parameter_name,
                        log            => $log,
                    }
                );
                $parameter_href->{$parameter_name}{$updated_key} =
                  delete( $parameter_href->{$parameter_name}{$key} );
            }
        }
        else {
            ## Scalar - not a reference

            ## Find aboslute path for supplied path or croaks and exists if path does not exists
            $parameter_href->{$parameter_name} = set_absolute_path(
                {
                    path           => $parameter_href->{$parameter_name},
                    parameter_name => $parameter_name,
                }
            );
        }
    }
    return;
}

sub check_user_reference {

## Function : Check that the user supplied reference id and version
## Returns  :
## Arguments: $cmd_reference_ref => User supplied reference id and version
##          : $reference_ref     => Defined reference id and version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cmd_reference_ref;
    my $reference_ref;

    my $tmpl = {
        cmd_reference_ref => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$cmd_reference_ref
        },
        reference_ref => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$reference_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger(q{mip_download});

  REFERENCE:
    while ( my ( $reference_id, $version ) = each %{$cmd_reference_ref} ) {

        if ( not exists $reference_ref->{$reference_id} ) {

            $log->fatal( q{Cannot find reference key:} . $reference_id );
            exit 1;
        }
        elsif (
            not any { $_ eq $version }
            @{ $reference_ref->{$reference_id} }
          )
        {
            ## If element is part of array

            $log->fatal(
                    q{Cannot find version key: }
                  . $version
                  . q{ reference key:}
                  . $reference_id,
            );
            exit 1;
        }
    }
    return;
}

sub write_post_processing_command {

## Function : Writes post processing commands associated with reference e.g. tabic
## Returns  :
## Arguments: $reference_href => Reference hash {REF}
##          : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

## Flatten argument(s)
    my $reference_href;
    my $FILEHANDLE;

    my $tmpl = {
        reference_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$reference_href,
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Process reference with commands
    my @processing_commands =
      qw{ outfile_reformat_command outfile_bgzip_command outfile_tabix_command};

    if ( scalar @{$reference_href}{@processing_commands} ) {

        ## We have a command to process add ;
        print {$FILEHANDLE} q{;};
    }
  COMMAND:
    foreach my $command (@processing_commands) {

        next COMMAND if ( not exists $reference_href->{$command} );

        ## Command
        say {$FILEHANDLE} $reference_href->{$command}, $NEWLINE;
    }
    return;
}

1;
