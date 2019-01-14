package MIP::Program::Variantcalling::Svdb;

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
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ svdb_merge svdb_query };
}

## Constants
Readonly my $SPACE => q{ };

sub svdb_merge {

## Function : Perl wrapper for writing svdb merge recipe to $FILEHANDLE or return commands array. Based on svdb 1.0.7.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_paths_ref       => Infile path {REF}
##          : $notag                  => Do not add the the VARID and set entries to the info field
##          : $outfile_path           => Outfile path
##          : $priority               => Priority order of structural variant calls
##          : $same_order             => Across all input vcf files, the order of the sample columns are the same
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $notag;
    my $outfile_path;
    my $priority;
    my $same_order;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE       => { store => \$FILEHANDLE, },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        notag           => { store => \$notag,        strict_type => 1, },
        outfile_path    => { store => \$outfile_path, strict_type => 1, },
        priority        => { store => \$priority,     strict_type => 1, },
        same_order      => { store => \$same_order,   strict_type => 1, },
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## svdb
    my @commands = qw{ svdb --merge };

    ## Options
    if ($priority) {

        ## Priority order of structural variant calls
        push @commands, q{--priority} . $SPACE . $priority;
    }
    if ($notag) {

        ## Do not tag variant with origin file
        push @commands, q{--notag};
    }
    if ($same_order) {

        ## Same sample order across vcf files
        push @commands, q{--same_order};
    }

    ## Infile
    push @commands, q{--vcf} . $SPACE . join $SPACE, @{$infile_paths_ref};

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );

    return @commands;
}

sub svdb_query {

## Function : Perl wrapper for writing svdb query recipe to $FILEHANDLE or return commands array. Based on svdb 0.1.2.
## Returns  : @commands
## Arguments: $bnd_distance           => Maximum distance between two similar precise breakpoints
##          : $dbfile_path            => Svdb database file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $frequency_tag          => Tag used to describe the frequency of the variant
##          : $hit_tag                => The tag used to describe the number of hits within the info field of the output vcf
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $overlap                => Overlap required to merge two events
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bnd_distance;
    my $dbfile_path;
    my $FILEHANDLE;
    my $frequency_tag;
    my $hit_tag;
    my $infile_path;
    my $outfile_path;
    my $overlap;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        bnd_distance => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$bnd_distance
        },
        dbfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$dbfile_path
        },
        FILEHANDLE    => { store       => \$FILEHANDLE },
        frequency_tag => { strict_type => 1, store => \$frequency_tag },
        hit_tag       => { strict_type => 1, store => \$hit_tag },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => { strict_type => 1, store => \$outfile_path },
        overlap      => {
            allow       => qr/ ^\d+ | d+[.]d+$ /sxm,
            strict_type => 1,
            store       => \$overlap
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## svdb
    my @commands = qw{ svdb --query };

    ## Options
    if ($bnd_distance) {

        push @commands, q{--bnd_distance} . $SPACE . $bnd_distance;
    }
    if ($overlap) {

        push @commands, q{--overlap} . $SPACE . $overlap;
    }
    if ($hit_tag) {

        push @commands, q{--hit_tag} . $SPACE . $hit_tag;
    }
    if ($frequency_tag) {

        push @commands, q{--frequency_tag} . $SPACE . $frequency_tag;
    }

    push @commands, q{--db} . $SPACE . $dbfile_path;

    ## Infile
    push @commands, q{--query_vcf} . $SPACE . $infile_path;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
