package MIP::Main::Download;

use 5.026;
use Carp;
use charnames qw( :full :short );
use Cwd;
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catfile catdir devnull };
use FindBin qw{ $Bin };
use Getopt::Long;
use IO::Handle;
use warnings qw{ FATAL utf8 };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use Time::Piece;
use utf8;

## CPANM
use autodie qw{ open close :all };
use List::Util qw{ any };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $COLON $DOT $MIP_VERSION $NEWLINE $SINGLE_QUOTE $SPACE $UNDERSCORE };
use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Gnu::Coreutils qw{ gnu_mkdir };
use MIP::Language::Shell qw{ create_bash_file };
use MIP::Log::MIP_log4perl qw{ initiate_logger set_default_log4perl_file };
use MIP::Script::Utils qw{ create_temp_dir };

BEGIN {
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

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

    ## Get local time
    my $date_time       = localtime;
    my $date_time_stamp = $date_time->datetime;
    my $date            = $date_time->ymd;

    # Catches script name and removes ending
    my $script = fileparse( basename( $PROGRAM_NAME, $DOT . q{pl} ) );

    # Add specific MIP process
    $script .= $UNDERSCORE . q{download};

    ## Set the default Log4perl file using supplied dynamic parameters.
    $parameter{log_file} = set_default_log4perl_file(
        {
            cmd_input       => $parameter{log_file},
            date            => $date,
            date_time_stamp => $date_time_stamp,
            script          => $script,
        }
    );

    ## Initiate logger
    my $log = initiate_logger(
        {
            file_path => $parameter{log_file},
            log_name  => uc q{mip_download},
        }
    );
    $log->info( q{MIP Version: } . $MIP_VERSION );
    $log->info( q{Writing log messages to} . $COLON . $SPACE . $parameter{log_file} );

## Set parameter default
    if ( not $parameter{reference_dir} ) {

        $parameter{reference_dir} = cwd();
    }

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

    # Install temp directory
    my $temp_dir = create_temp_dir( {} );

    if ( $parameter{sbatch_mode} ) {

        if ( not $parameter{project_id} ) {

            $log->fatal(q{Please provide a sbatch project id with option '--project_id'});
            exit 1;
        }
        $log->info(
q{Will write sbatch install instructions to for sbatch enabled references to individual sbatch scripts}
        );

    }

    # Downloads instruction file
    my $bash_file_path = catfile( cwd(), q{download_reference} . $DOT . q{sh} );

    open $FILEHANDLE, '>', $bash_file_path
      or
      $log->logdie( q{Cannot write to '} . $bash_file_path . q{' :} . $OS_ERROR . "\n" );

## Create bash file for writing install instructions
    create_bash_file(
        {
            FILEHANDLE  => $FILEHANDLE,
            file_name   => $bash_file_path,
            log         => $log,
            remove_dir  => $temp_dir,
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

    if ( $parameter{sbatch_mode} ) {

## Build install references recipe in bash file
        sbatch_build_reference_install_recipe(
            {
                FILEHANDLE     => $FILEHANDLE,
                parameter_href => \%parameter,
                temp_directory => $temp_dir,
            }
        );
    }
    else {

        ## Build install references recipe in bash file
        build_reference_install_recipe(
            {
                FILEHANDLE     => $FILEHANDLE,
                parameter_href => \%parameter,
                temp_directory => $temp_dir,
            }
        );
    }
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
##          : $temp_directory => Temporary directory
##          : $verbose        => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    ## Default(s)
    my $quiet;
    my $temp_directory;
    my $verbose;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        FILEHANDLE => { defined => 1, required => 1, store => \$FILEHANDLE, },
        quiet      => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbose => {
            default     => $arg_href->{parameter_href}{verbose},
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Bash qw{ gnu_cd };
    use MIP::Gnu::Coreutils qw{ gnu_mkdir };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Recipes::Download::Human_reference qw{ download_human_reference };

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger(q{Download_reference});

    my $pwd = cwd();

    ### Download recipes
    ## Create code reference table for download recipes
    my %download_recipe = ( human_reference => \&download_human_reference, );

    # Storing job_ids from SLURM
    my %job_id;

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

                get_reference(
                    {
                        FILEHANDLE     => $FILEHANDLE,
                        parameter_href => $parameter_href,
                        recipe_name    => $reference_id,
                        reference_href => $reference_href,
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

sub sbatch_build_reference_install_recipe {

## Function : Build install references recipe in bash file
## Returns  :
## Arguments: $parameter_href => Holds all parameters
##          : $FILEHANDLE     => Filehandle to write to
##          : $quiet          => Be quiet
##          : $temp_directory => Temporary directory
##          : $verbose        => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;

    ## Default(s)
    my $quiet;
    my $temp_directory;
    my $verbose;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        FILEHANDLE => { defined => 1, required => 1, store => \$FILEHANDLE, },
        quiet      => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbose => {
            default     => $arg_href->{parameter_href}{verbose},
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Bash qw{ gnu_cd };
    use MIP::Gnu::Coreutils qw{ gnu_mkdir };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Recipes::Download::Human_reference qw{ download_human_reference };

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger(q{Download_reference});

    my $pwd = cwd();

    ### Download recipes
    ## Create code reference table for download recipes
    my %download_recipe = ( human_reference => \&download_human_reference, );

    # Storing job_ids from SLURM
    my %job_id;

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

                if ( exists $download_recipe{$reference_id} ) {

                    $download_recipe{$reference_id}->(
                        {
                            job_id_href       => \%job_id,
                            parameter_href    => $parameter_href,
                            recipe_name       => $reference_id,
                            reference_href    => $reference_href,
                            reference_version => $reference_version,
                            quiet             => $quiet,
                            temp_directory    => $temp_directory,
                        }
                    );
                }
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

  COMMAND:
    foreach my $command (@processing_commands) {

        next COMMAND if ( not exists $reference_href->{$command} );

        ## Command
        say {$FILEHANDLE} $reference_href->{$command}, $NEWLINE;
    }
    return;
}

1;
