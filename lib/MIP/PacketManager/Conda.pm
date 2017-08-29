package MIP::PacketManager::Conda;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

use Getopt::Long;
use Cwd;
use Cwd qw{ abs_path };
use FindBin qw{ $Bin };               #Find directory of script
use IO::Handle;
use File::Basename qw{ dirname basename fileparse };
use File::Spec::Functions qw{ catfile catdir devnull };
use Readonly;
use List::Util qw{ none first };
use IPC::Cmd qw{ can_run run };

## MIPs lib/
use MIP::Unix::Write_to_file qw{ unix_write_to_file };
use Program::Download::Wget qw(wget);
use MIP::Gnu::Bash qw(gnu_cd);
use MIP::Gnu::Coreutils qw(gnu_cp gnu_rm gnu_mv gnu_mkdir gnu_link );

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};

BEGIN {

    use base qw{Exporter};
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{conda_create conda_source_activate conda_source_deactivate conda_update conda_check conda_install };
}

sub conda_create {

##conda_create

##Function : Create Conda environment
##Returns  : @commands
##Arguments: $packages_ref, $env_name, $python_version, $quiet, $no_confirmation,
##           $FILEHANDLE
##         : $packages_ref    => Packages to be installed
##         : $FILEHANDLE      => Filehandle to write to
##         : $env_name        => Name of environment to create
##         : $quiet           => Do not display progress bar
##         : $no_confirmation => Do not ask for confirmation

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $env_name;
    my $quiet;
    my $no_confirmation;
    my $packages_ref;
    my $FILEHANDLE;

    my $tmpl = {
        packages_ref => {
            default     => [],
            strict_type => 1,
            store       => \$packages_ref
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE
        },
        env_name => {
            default     => q{},
            strict_type => 1,
            store       => \$env_name
        },
        quiet => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        no_confirmation => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$no_confirmation
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ##Stores commands depending on input parameters
    # Basic command
    my @commands = q{conda create};

    if ($env_name) {
        push @commands, q{--name} . $SPACE . $env_name;
    }

    if ($quiet) {
        push @commands, q{--quiet};
    }

    if ($no_confirmation) {
        push @commands, q{--yes};
    }

    if ( @{$packages_ref} ) {
        push @commands, join $SPACE, @{$packages_ref};
    }

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub conda_source_activate {

##conda_source_activate

##Function : Activate conda environment
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $env_name,
##         : $FILEHANDLE => Filehandle to write to
##         : $env_name   => Name of conda environment

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $env_name;

    my $tmpl = {
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE
        },
        env_name => {
            required    => 1,
            strict_type => 1,
            store       => \$env_name
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters

    # Basic command
    my @commands = q{source activate};

    # Activates env, default root
    push @commands, $env_name;

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub conda_source_deactivate {

## conda_source_deactivate

##Function : Deactivate conda environment
##Returns  : "@commands"
##Arguments: $FILEHANDLE
##         : $FILEHANDLE => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = {
        FILEHANDLE => {
            required => 1,
            defined  => 1,
            store    => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = q{source deactivate};

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
        },
    );

    return @commands;
}

sub conda_update {

## conda_update

## Function  : Update conda
## Returns   : @commands
## Arguments : $FILEHANDLE, $no_confirmation
##           : $FILEHANDLE      => Filehandle to write to
##           : $no_confirmation => Do not ask for confirmation

    my ($arg_href) = @_;

    ## Flatten arguments
    my $FILEHANDLE;
    my $no_confirmation;

    my $tmpl = {
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE
        },
        no_confirmation => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$no_confirmation
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = q{conda update};

    if ($no_confirmation) {
        push @commands, q{--yes};
    }

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub conda_check {

## conda_check

## Function  : Check if the path to conda is correct and if a conda
##           : environment is active (exit if true). Returns path to conda.
## Returns   : $conda_dir_path
## Arguments : $conda_dir_path_ref
##           : $conda_dir_path_ref => Path to conda directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_dir_path_ref;

    my $tmpl = {
        conda_dir_path_ref => {
            default     => [],
            required    => 1,
            strict_type => 1,
            store       => \$conda_dir_path_ref,
        },
    };
    
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};


    my $conda_dir_path;

    # Search for conda directory in supplied path to conda
    if ( none { -d $_ } @{$conda_dir_path_ref} ) {
        say STDERR q{Could not find miniconda directory in:} . $SPACE
          . join $SPACE, @{$conda_dir_path_ref};
        exit 1;
    }
    else {
        $conda_dir_path = first { -d $_ } @{$conda_dir_path_ref};
    }

    ## Deactivate any activate env prior to installation
    #   Perl options:
    #   n : loop over input
    #   a : automatically split input and store in array @F
    #   e : execute code
    #
    # Unless the active environment is root the expression will return true
    #   and print the environment name
    my $detect_active_conda_env =
      q?perl -nae 'if( ($_!~/^root/) && ($_=~/\*/) ) {print $F[0]}'?;

    # Pipes the output from the shell command "conda info --envs"
    #   to $detect_active_conda_env.
    #   Output is captured in $run_output.
    my $run_output;
    run(
        command => qq{conda info --envs | $detect_active_conda_env},
        buffer  => \$run_output
    );

    if ($run_output) {
        say STDOUT q{Found activated conda env:} . $SPACE . $run_output;
        say STDOUT q{Please exit conda env:} . $SPACE . $run_output . $SPACE
          . q{with 'source deactivate' before executing install script};
        exit 1;
    }
    
    return $conda_dir_path;
}

sub conda_install {

## conda_install

##Function : Install packages into conda environment
##Returns  : @commands
##Arguments: 
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $packages_ref;
    my $conda_channel;
    my $env_name;
    my $FILEHANDLE;
    my $quiet;
    my $no_confirmation;

    my $tmpl = {
        parameter_href => {
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        packages_ref => {
            required => 1,
            defined => 1,
            default => [],
            strict_type => 1,
            store => \$packages_ref
        },
        conda_channel => {
            defined => 1,
            strict_type => 1,
            store => \$conda_channel
        },
        env_name => {
            strict_type => 1,
            store => \$env_name
        },
        FILEHANDLE => { 
            required => 1, 
            defined => 1, 
            store => \$FILEHANDLE 
        },
        quiet => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        no_confirmation => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$no_confirmation
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = q{conda install};

    if ( $env_name ) {
        push @commands, q{--name} . $SPACE . $env_name;
    }

    if ( $quiet ) {
        #Do not display progress bar
        push @commands, q{--quiet};   
    }

    if ( $no_confirmation ) {
        push @commands, q{--yes};
    }

    if ( $conda_channel ) {
        push @commands, q{--channel} . $SPACE . $conda_channel; 
    }

    push @commands, join $SPACE, @{$packages_ref};

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    
    ## Install all bioconda packages
    foreach my $program ( keys %{ $parameter_href->{bioconda} } ) {

        print $FILEHANDLE $program . '='
          . $parameter_href->{bioconda}{$program} . q{ };
    }

    print $FILEHANDLE "\n\n";

    ## Custom
    foreach my $program ( keys %{ $parameter_href->{bioconda} } ) {

        if ( $program eq 'bwakit' ) {

            ## Define binaries
            my @bwakit_binaries = (
                'k8',                'seqtk',
                'bwa-postalt.js',    'run-HLA',
                'typeHLA.sh',        'fermi2',
                'fermi2.pl',         'ropebwt2',
                'typeHLA-selctg.js', 'typeHLA.js'
            );

            foreach my $binary (@bwakit_binaries) {
                # Specifying target and link paths
                my $target_path = catfile(
                    $parameter_href->{conda_prefix_path}, q{share},
                    q{bwakit-} . $parameter_href->{bioconda}{bwakit}
                    . $parameter_href->{bioconda_bwakit_patch}, $binary
                );
                my $link_path = catfile(
                    $parameter_href->{conda_prefix_path}, q{bin}, $binary
                );
                gnu_link(
                    {
                        FILEHANDLE  => $FILEHANDLE,
                        target_path => $target_path,
                        link_path   => $link_path,
                        symbolic    => 1,
                        force       => 1,
                    }
                );
                print $FILEHANDLE $NEWLINE;
            }
            print $FILEHANDLE $NEWLINE;

            gnu_cp(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    recursive   => 1,
                    force       => 1,
                    infile_path => catdir(
                        $parameter_href->{conda_prefix_path},
                        'share',
                        'bwakit-'
                          . $parameter_href->{bioconda}{bwakit}
                          . $parameter_href->{bioconda_bwakit_patch},
                        'resource-human-HLA'
                    ),
                    outfile_path =>
                      catdir( $parameter_href->{conda_prefix_path}, 'bin' ),
                }
            );
            print $FILEHANDLE "\n\n";
        }

        if ( $program eq 'picard' ) {
            # Specifying target and link paths
            my $target_path = catfile(
                $parameter_href->{conda_prefix_path}, q{share}, q{picard-}
                . $parameter_href->{bioconda}{picard}
                . $parameter_href->{bioconda_picard_patch}, q{picard.jar}
            );
            my $link_path = catfile(
                $parameter_href->{conda_prefix_path}, q{picard.jar}
            );
            gnu_link(
                {
                    FILEHANDLE  => $FILEHANDLE,
                    target_path => $target_path,
                    link_path   => $link_path,
                    symbolic    => 1,
                    force       => 1,
                }
            );
            print $FILEHANDLE $NEWLINE;
        }

        if ( $program eq 'snpeff' ) {
            ## Define binaries
            my @snpeff_binaries = qw(snpEff.jar snpEff.config);

            foreach my $binary (@snpeff_binaries) {
                # Specifying target and link paths
                my $target_path = catfile(
                    $parameter_href->{conda_prefix_path}, q{share},
                    q{snpeff-} . $parameter_href->{bioconda}{snpeff}
                    . $parameter_href->{bioconda_snpeff_patch}, $binary
                );
                my $link_path = catfile(
                    $parameter_href->{conda_prefix_path}, $binary
                );
                gnu_link(
                    {
                        FILEHANDLE  => $FILEHANDLE,
                        target_path => $target_path,
                        link_path   => $link_path,
                        symbolic    => 1,
                        force       => 1,
                    }
                );
                print $FILEHANDLE $NEWLINE;
            }
            print $FILEHANDLE $NEWLINE;

            foreach my $genome_version (
                @{ $parameter_href->{snpeff_genome_versions} } )
            {

                ## Check and if required add the vertebrate mitochondrial codon table to snpeff config
                check_mt_codon_table(
                    {
                        parameter_href => $parameter_href,
                        FILEHANDLE     => $FILEHANDLE,
                        share_dir      => catdir(
                            $parameter_href->{conda_prefix_path},
                            'share',
                            'snpeff-'
                              . $parameter_href->{bioconda}{snpeff}
                              . $parameter_href->{bioconda_snpeff_patch}
                        ),
                        config_file        => 'snpEff.config',
                        genome_version_ref => \$genome_version,
                    }
                );

                unless (
                    -d catdir(
                        $parameter_href->{conda_prefix_path},
                        'share',
                        'snpeff-'
                          . $parameter_href->{bioconda}{snpeff}
                          . $parameter_href->{bioconda_snpeff_patch},
                        'data',
                        $genome_version
                    )
                  )
                {

                    ## Write instructions to download snpeff database.
                    ## This is done by install script to avoid race conditin when doing first analysis run in MIP
                    snpeff_download(
                        {
                            parameter_href     => $parameter_href,
                            FILEHANDLE         => $FILEHANDLE,
                            genome_version_ref => \$genome_version,
                        }
                    );
                }
            }
        }

        if ( $program eq 'snpsift' ) {
            ## Define binaries
            my @snpsift_binaries = qw(SnpSift.jar);

            foreach my $binary (@snpsift_binaries) {
                ## Specifying target and link paths
                my $target_path = catfile(
                    $parameter_href->{conda_prefix_path}, q{share},
                    q{snpsift-} . $parameter_href->{bioconda}{snpsift}
                    . $parameter_href->{bioconda_snpsift_patch}, $binary
                );
                my $link_path = catfile(
                    $parameter_href->{conda_prefix_path}, $binary
                );
                gnu_link(
                    {
                        FILEHANDLE  => $FILEHANDLE,
                        target_path => $target_path,
                        link_path   => $link_path,
                        symbolic    => 1,
                        force       => 1,
                    }
                );
                print $FILEHANDLE $NEWLINE;
            }
            print $FILEHANDLE $NEWLINE;
        }

        if ( $program eq 'manta' ) {
            ## Define binaries
            my @manta_binaries = qw(configManta.py configManta.py.ini);

            foreach my $binary (@manta_binaries) {
                ## Specifying target and link paths
                my $target_path = catfile(
                    $parameter_href->{conda_prefix_path}, q{share}, q{manta-}
                    . $parameter_href->{bioconda}{manta}
                    . $parameter_href->{bioconda_manta_patch}, q{bin}, $binary
                );
                my $link_path = catfile(
                    $parameter_href->{conda_prefix_path}, $binary
                );
                gnu_link(
                    {
                        FILEHANDLE  => $FILEHANDLE,
                        target_path => $target_path,
                        link_path   => $link_path,
                        symbolic    => 1,
                        force       => 1,
                    }
                );
                print $FILEHANDLE $NEWLINE;
            }
            print $FILEHANDLE $NEWLINE;

            ## Make file executable
            enable_executable(
                {
                    parameter_href => $parameter_href,
                    FILEHANDLE     => $FILEHANDLE,
                    binary         => q?configManta.py?,
                }
            );
        }
    }
    return;
}

sub enable_executable {

##enable_executable

##Function : Make file executable
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $binary
##         : $parameter_href => Holds all parameters
##         : $FILEHANDLE     => FILEHANDLE to write to
##         : $binary         => The binary file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $binary;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        binary =>
          { required => 1, defined => 1, strict_type => 1, store => \$binary },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();

    ## Enable executable
    print $FILEHANDLE '## Enable executable', "\n";
    gnu_cd(
        {
            directory_path =>
              catdir( $parameter_href->{conda_prefix_path}, 'bin' ),
            FILEHANDLE => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n";

    print $FILEHANDLE 'chmod a+x ';
    print $FILEHANDLE $binary . q{ };
    print $FILEHANDLE "\n\n";

    ## Move to back
    print $FILEHANDLE '## Move to original working directory', "\n";
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    print $FILEHANDLE "\n\n";

    return;
}

sub check_mt_codon_table {

##check_mt_codon_table

##Function : Check and if required add the vertebrate mitochondrial codon table to snpeff config
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $share_dir, $config_file, $genome_version_ref
##         : $parameter_href     => Holds all parameters
##         : $FILEHANDLE         => FILEHANDLE to write to
##         : $share_dir          => The conda env shared directory
##         : $config_file        => The config config_file
##         : $genome_version_ref => snpeff genome version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $share_dir;
    my $config_file;
    my $genome_version_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        share_dir  => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$share_dir
        },
        config_file => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$config_file
        },
        genome_version_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$genome_version_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my $pwd = cwd();
    my $detect_regexp =
        q?perl -nae 'if($_=~/?
      . $$genome_version_ref
      . q?.MT.codonTable/) {print 1}' ?;
    my $add_regexp =
        q?perl -nae 'if($_=~/?
      . $$genome_version_ref
      . q?.reference/) {print $_; print "?
      . $$genome_version_ref
      . q?.MT.codonTable : Vertebrate_Mitochondrial\n"} else {print $_;}' ?;
    my $ret;

    if ( -f catfile( $share_dir, $config_file ) ) {

        $ret = `$detect_regexp $share_dir/$config_file`;
    }
    if ( !$ret ) {    #No MT.codonTable in config

        print $FILEHANDLE '## Adding '
          . $$genome_version_ref
          . '.MT.codonTable : Vertebrate_Mitochondrial to '
          . $share_dir
          . $config_file, "\n";

        ## Add MT.codon Table to config
        print $FILEHANDLE $add_regexp . q{ }
          . catfile( $share_dir, $config_file ) . ' > '
          . catfile( $share_dir, $config_file . '.tmp' ), "\n";
        gnu_mv(
            {
                infile_path  => catfile( $share_dir, $config_file . '.tmp' ),
                outfile_path => catfile( $share_dir, $config_file ),
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        print $FILEHANDLE "\n\n";

    }
    else {

        print STDERR 'Found MT.codonTable in '
          . catfile( $share_dir, 'snpEff.config' )
          . '. Skipping addition to snpEff config', "\n";
    }
    return;
}


sub snpeff_download {

##snpeff_download

##Function : Write instructions to download snpeff database. This is done by install script to avoid race conditin when doing first analysis run in MIP
##Returns  : ""
##Arguments: $parameter_href, $FILEHANDLE, $genome_version_ref
##         : $parameter_href     => Holds all parameters
##         : $FILEHANDLE         => FILEHANDLE to write to
##         : $genome_version_ref => snpeff genome version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $FILEHANDLE;
    my $genome_version_ref;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        genome_version_ref => {
            required    => 1,
            defined     => 1,
            default     => \$$,
            strict_type => 1,
            store       => \$genome_version_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    ## Only activate conda environment if supplied by user
    if ( $parameter_href->{conda_environment} ) {
        ## Activate conda environment
        say $FILEHANDLE q{## Activate conda environment};
        conda_source_activate(
            {
                FILEHANDLE => $FILEHANDLE,
                env_name   => $parameter_href->{conda_environment},
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    print $FILEHANDLE 'java -Xmx2g ';
    print $FILEHANDLE '-jar '
      . catfile( $parameter_href->{conda_prefix_path}, 'bin', 'snpEff.jar' )
      . q{ };
    print $FILEHANDLE 'download ';
    print $FILEHANDLE ' -v ';
    print $FILEHANDLE $$genome_version_ref . q{ };
    print $FILEHANDLE '-c '
      . catfile( $parameter_href->{conda_prefix_path}, 'bin', 'snpEff.config' )
      . q{ };
    print $FILEHANDLE "\n\n";

    ## Deactivate conda environment if conda_environment exists
    if ( $parameter_href->{conda_environment} ) {
        say $FILEHANDLE q{## Deactivate conda environment};
        conda_source_deactivate(
            {
                FILEHANDLE => $FILEHANDLE,
            }
        );
        say $FILEHANDLE $NEWLINE;
    }

    return;
}


1;
