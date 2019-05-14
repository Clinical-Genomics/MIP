package MIP::Package_manager::Pip;

use 5.026;
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
use MIP::Language::Python qw{ python_core };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_pip_package pip_install };
}

## Constants
Readonly my $SPACE => q{ };

sub check_pip_package {

## Function : Check if the package has been installed via pip
## Returns  : $status
## Arguments: $conda_environment => Name of conda environment
##          : $conda_prefix_path => Path to conda environment
##          : $package           => Pip package to check
##          : $version           => Optional version check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_environment;
    my $conda_prefix_path;
    my $package;
    my $version;

    my $tmpl = {
        conda_environment => {
            store       => \$conda_environment,
            strict_type => 1,
        },
        conda_prefix_path => {
            store       => \$conda_prefix_path,
            strict_type => 1,
        },
        package => {
            defined     => 1,
            required    => 1,
            store       => \$package,
            strict_type => 1,
        },
        version => {
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $check_pip_package_regexp = _build_package_check_regexp(
        {
            package => $package,
            version => $version,
        }
    );

    # Variable for storing pip package status
    my $status;

    # Shell command to launch
    my $command;

    # Check if the program is to be installed into a conda env
    if ($conda_environment) {

        # Check if the environemnt already exists
        if ( $conda_prefix_path and -d $conda_prefix_path ) {

            # Test if the program already exists in that environment
            $command =
qq{conda list -n $conda_environment | grep 'pip' | $check_pip_package_regexp};
            run(
                command => $command,
                buffer  => \$status
            );
        }
        return $status;
    }

    # Test if the program is already installed in the root env
    $command = qq{pip list --format columns | $check_pip_package_regexp};
    run(
        command => $command,
        buffer  => \$status
    );

    return $status;
}

sub pip_install {

## Function : Perl wrapper for writing pip install command. Based on pip v18.1
## Returns  : @commands
##          : $editable               => Install in editable mode from a local project path or a VCS url.
##          : $FILEHANDLE             => Filehandle to write to
##          : $packages_ref           => Array of packages to be installed {REF}
##          : $python_module          => Execute pip as python module
##          : $quiet                  => Quiet output
##          : $requirement            => Install from a requirement file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $verbose                => Verbose output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $editable;
    my $FILEHANDLE;
    my $packages_ref;
    my $no_cache_dir;
    my $python_module;
    my $quiet;
    my $requirement;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $verbose;

    my $tmpl = {
        editable => {
            store       => \$editable,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        packages_ref => {
            default     => [],
            defined     => 1,
            store       => \$packages_ref,
            strict_type => 1,
        },
        python_module => {
            allow => [ undef, 0, 1 ],
            store => \$python_module,
        },
        quiet => {
            allow => [ undef, 0, 1 ],
            store => \$quiet,
        },
        requirement => {
            store       => \$requirement,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        verbose => {
            allow => [ undef, 0, 1 ],
            store => \$verbose,
        },
        no_cache_dir => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$no_cache_dir,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands;

    ## Prepend 'python -m' to base command
    if ($python_module) {

        push @commands, join $SPACE, python_core( { module_mode => 1, } );
    }

    ## Push base command
    push @commands, qw{ pip install };

    if ($no_cache_dir) {
        push @commands, q{--no-cache-dir};
    }

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

sub _build_package_check_regexp {

## Function : Build regexp for package check
## Returns  : $check_pip_package_regexp
## Arguments: $package => Pip package to check
##          : $version => Optional version check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $package;
    my $version;

    my $tmpl = {
        package => {
            defined     => 1,
            required    => 1,
            store       => \$package,
            strict_type => 1,
        },
        version => {
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $check_pip_package_regexp;

    if ($version) {

        $check_pip_package_regexp =

          # Execute perl, loop over input and split on whitespace
          q?perl -nae? . $SPACE

          # Check for a matching package
          . q?'if( ($_=~/? . $package . q?/) && ($_=~/? . $version . q?/) )?

          # Print 1 in case of match
          . q?{print 1}'?;
    }
    else {

        $check_pip_package_regexp =

          # Execute perl, loop over input and split on whitespace
          q?perl -nae? . $SPACE

          # Check for a matching package
          . q?'if($_=~/? . $package . q?/)?

          # Print 1 in case of match
          . q?{print 1}' ?;
    }

    return $check_pip_package_regexp;

}

1;
