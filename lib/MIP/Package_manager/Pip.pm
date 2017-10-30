package MIP::Package_manager::Pip;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use Readonly;
use IPC::Cmd qw{ run };

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pip_install pip_check_package};
}

## Constants
Readonly my $SPACE => q{ };

sub pip_install {

## Function : Perl wrapper for writing pip install command. Based on pip v9.0.1.
## Returns  : "@commands"
## Arguments: $packages_ref           => Array of packages to be installed {REF}
##          : $quiet                  => Quiet output
##          : $verbose                => Verbose output
##          : $requirement            => Install from a requirement file
##          : $editable               => Install in editable mode from a local project path or a VCS url.
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $packages_ref;
    my $quiet;
    my $verbose;
    my $requirement;
    my $editable;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        packages_ref => {
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$packages_ref
        },
        quiet => {
            allow => [ undef, 0, 1 ],
            store => \$quiet
        },
        verbose => {
            allow => [ undef, 0, 1 ],
            store => \$verbose
        },
        requirement => {
            strict_type => 1,
            store       => \$requirement
        },
        editable => {
            strict_type => 1,
            store       => \$editable
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{pip install};

    if ($quiet) {
        push @commands, q{--quiet};
    }

    if ($verbose) {
        push @commands, q{--verbose};
    }

    if ($requirement) {
        push @commands, q{--requirement} . $SPACE . $requirement;
    }
    elsif ($editable) {
        push @commands, q{--editable} . $SPACE . $editable;
    }
    else {
        push @commands, join $SPACE, @{$packages_ref};
    }

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub pip_check_package {

## Function : Check if the package has been installed via pip
## Returns  : $status
## Arguments: $package           => Pip package to check
##          : $version           => optional version check
##          : $conda_environment => Name of conda environment
##          : $conda_prefix_path => Path to conda environment

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $package;
    my $version;
    my $conda_environment;
    my $conda_prefix_path;

    my $tmpl = {
        package => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$package,
        },
        version => {
            strict_type => 1,
            store       => \$version,
        },
        conda_environment => {
            strict_type => 1,
            store       => \$conda_environment,
        },
        conda_prefix_path => {
            strict_type => 1,
            store       => \$conda_prefix_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $detect_regexp;

    if ($version) {
        $detect_regexp =

          # Execute perl, loop over input and split on whitespace
          q?perl -nae? . $SPACE

          # Check for a matching package
          . q?'if( ($_=~/? . $package . q?/) && ($_=~/? . $version . q?/) )?

          # Print 1 in case of match
          . q?{print 1}'?;
    }
    else {
        $detect_regexp =

          # Execute perl, loop over input and split on whitespace
          q?perl -nae? . $SPACE

          # Check for a matching package
          . q?'if($_=~/? . $package . q?/)?

          # Print 1 in case of match
          . q?{print 1}' ?;
    }

    my $status;

    # shell command to launch
    my $command;

    # Check if the program is to be installed into a conda env
    if ($conda_environment) {

        # Check if the environemnt already exists
        if ( -d $conda_prefix_path ) {

            # Test if the program already exists in that environment
            $command =
qq{conda list -n $conda_environment | grep 'pip' | $detect_regexp};
            run(
                command => $command,
                buffer  => \$status
            );
        }
    }
    else {
        #Test if the program is already installed in the root env
        $command = qq{pip list --format columns | $detect_regexp};
        run(
            command => $command,
            buffer  => \$status
        );

    }

    return $status;

}
1;
