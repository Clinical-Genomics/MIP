package MIP::Workloadmanager::Slurm;

#### Copyright 2017 Henrik Stranneheim

use strict;
use warnings;
use warnings qw(FATAL utf8);
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

use FindBin qw($Bin);                 #Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw(unix_standard_streams);
use MIP::Unix::Write_to_file qw(unix_write_to_file);

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(slurm_sacct slurm_sbatch);
}

sub slurm_sacct {

##slurm_sacct

##Function : Perl wrapper for writing SLURM sacct recipe to already open $FILEHANDLE or return commands array. Based on SLURM sacct 2.6.0.
##Returns  : "@commands"
##Arguments: $fields_format_ref, $job_ids_ref, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $stderrfile_path_append,
##         : $job_ids_ref        => Slurm job id
##         : $fields_format_ref  => List of format fields
##         : $stderrfile_path    => Stderrfile path
##         : $stdoutfile_path    => Stdoutfile path
##         : $FILEHANDLE         => Filehandle to write to
##         : $stderrfile_path_append => Append stderr info to file

    my ($arg_href) = @_;

    ## Default(s)
    my $stderrfile_path_append;

    ## Flatten argument(s)
    my $fields_format_ref;
    my $job_ids_ref;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        fields_format_ref =>
          { default => [], strict_type => 1, store => \$fields_format_ref },
        job_ids_ref =>
          { default => [], strict_type => 1, store => \$job_ids_ref },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    my $SPACE = q{ };
    my $COMMA = q{,};

    ## sacct
    my @commands = qw(sacct);    #Stores commands depending on input parameters

    ## Options
    if ( @{$fields_format_ref} ) {

        push @commands, '--format=' . join $COMMA, @{$fields_format_ref};
    }
    if ( @{$job_ids_ref} ) {

        push @commands, '--jobs ' . join $COMMA, @{$job_ids_ref};
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

sub slurm_sbatch {

##slurm_sbatch

##Function : Perl wrapper for writing SLURM sbatch recipe to already open $FILEHANDLE or return commands array. Based on SLURM sbatch 2.6.0.
##Returns  : "@commands"
##Arguments: infile_path, $dependency_type, $job_ids_string, $stderrfile_path, $stdoutfile_path, $FILEHANDLE, $stderrfile_path_append,
##         : $infile_path            => Infile_path
##         : $dependency_type        => Type of slurm job dependency
##         : $job_ids_string         => Slurm job ids string (:job_id:job_id...n)
##         : $stderrfile_path        => Stderrfile path
##         : $stdoutfile_path        => Stdoutfile path
##         : $FILEHANDLE             => Filehandle to write to
##         : $stderrfile_path_append => Append stderr info to file

    my ($arg_href) = @_;

    ## Default(s)
    my $stderrfile_path_append;

    ## Flatten argument(s)
    my $infile_path;
    my $dependency_type;
    my $job_ids_string;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        dependency_type => { strict_type => 1, store => \$dependency_type },
        job_ids_string  => { strict_type => 1, store => \$job_ids_string },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        FILEHANDLE => { store => \$FILEHANDLE },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    my $SPACE = q{ };

    ## sbatch
    my @commands = qw(sbatch);    #Stores commands depending on input parameters

    ## Options
    if ( ($dependency_type) && ($job_ids_string) ) {

        push @commands, '--dependency=' . $dependency_type . $job_ids_string;
    }

    # Infile
    push @commands, $infile_path;

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

1;
