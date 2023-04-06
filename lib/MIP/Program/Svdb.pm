package MIP::Program::Svdb;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ svdb_merge svdb_query };
}

Readonly my $BASE_COMMAND => q{svdb};

sub svdb_merge {

## Function : Perl wrapper for writing svdb merge recipe to $filehandle or return commands array. Based on svdb 1.0.7.
## Returns  : @commands
## Arguments: $bnd_distance           => Maximum distance between two similar precise breakpoints
##          : $filehandle             => Filehandle to write to
##          : $infile_paths_ref       => Infile path {REF}
##          : $notag                  => Do not add the the VARID and set entries to the info field
##          : $outfile_path           => Outfile path
##          : $overlap                => Overlap required to merge two events
##          : $pass_only              => Only merge variants labeled PASS
##          : $priority               => Priority order of structural variant calls
##          : $same_order             => Across all input vcf files, the order of the sample columns are the same
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bnd_distance;
    my $filehandle;
    my $infile_paths_ref;
    my $notag;
    my $outfile_path;
    my $overlap;
    my $pass_only;
    my $priority;
    my $same_order;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        bnd_distance => {
            allow       => [ undef, qr/ \A \d+ \z /sxm ],
            strict_type => 1,
            store       => \$bnd_distance
        },
        filehandle       => { store => \$filehandle, },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        notag        => { store => \$notag,        strict_type => 1, },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        overlap      => {
            allow       => [ undef, qr/ \A [-]? \d+ | \d+[.]\d+ \z /sxm ],
            strict_type => 1,
            store       => \$overlap
        },
        pass_only => {
            allow       => [ undef, 0, 1 ],
            store       => \$pass_only,
            strict_type => 1,
        },
        priority        => { store => \$priority,   strict_type => 1, },
        same_order      => { store => \$same_order, strict_type => 1, },
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

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ --merge } );

    if ($bnd_distance) {

        push @commands, q{--bnd_distance} . $SPACE . $bnd_distance;
    }
    if ($notag) {

        push @commands, q{--notag};
    }
    if ($overlap) {

        push @commands, q{--overlap} . $SPACE . $overlap;
    }
    if ($pass_only) {

        push @commands, q{--pass_only},;
    }
    if ($priority) {

        push @commands, q{--priority} . $SPACE . $priority;
    }
    if ($same_order) {

        push @commands, q{--same_order};
    }

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
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );

    return @commands;
}

sub svdb_query {

## Function : Perl wrapper for writing svdb query recipe to $filehandle or return commands array. Based on svdb 2.0.0.
## Returns  : @commands
## Arguments: $bedpedb_path           => File path to bedpe file
##          : $bnd_distance           => Maximum distance between two similar precise breakpoints
##          : $dbfile_path            => Svdb database file path
##          : $filehandle             => Filehandle to write to
##          : $in_frequency_tag       => The frequency count tag, if used, this tag must be present in the INFO column of the input DB (usually AF or FRQ)
##          : $in_allele_count_tag    => The allele count tag, if used, this tag must be present in the INFO column of the input DB (usually AC or OCC)
##          : $infile_path            => Infile path
##          : $no_var                 => Run without taking variant type into consideration
##          : $outfile_path           => Outfile path
##          : $out_allele_count_tag   => The allele count tag output name
##          : $out_frequency_tag      => The frequency count tag output name
##          : $overlap                => Overlap required to merge two events
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bedpedb_path;
    my $bnd_distance;
    my $dbfile_path;
    my $filehandle;
    my $in_frequency_tag;
    my $in_allele_count_tag;
    my $infile_path;
    my $no_var;
    my $outfile_path;
    my $out_allele_count_tag;
    my $out_frequency_tag;
    my $overlap;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        bedpedb_path => {
            store       => \$bedpedb_path,
            strict_type => 1,
        },
        bnd_distance => {
            allow       => qr/ \A \d+ \z /sxm,
            strict_type => 1,
            store       => \$bnd_distance
        },
        dbfile_path => {
            strict_type => 1,
            store       => \$dbfile_path
        },
        filehandle          => { store       => \$filehandle },
        in_frequency_tag    => { strict_type => 1, store => \$in_frequency_tag },
        in_allele_count_tag => { strict_type => 1, store => \$in_allele_count_tag },
        infile_path         => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        no_var => {
            allow       => [ undef, 0, 1, ],
            store       => \$no_var,
            strict_type => 1,
        },
        outfile_path         => { strict_type => 1, store => \$outfile_path },
        out_allele_count_tag => { strict_type => 1, store => \$out_allele_count_tag },
        out_frequency_tag    => { strict_type => 1, store => \$out_frequency_tag },
        overlap              => {
            allow       => qr/ \A [-]? \d+ | \d+[.]\d+ \z /sxm,
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

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ --query } );

    if ($bnd_distance) {

        push @commands, q{--bnd_distance} . $SPACE . $bnd_distance;
    }
    if ($overlap) {

        push @commands, q{--overlap} . $SPACE . $overlap;
    }
    if ($no_var) {

        push @commands, q{--no_var};
    }
    if ($in_allele_count_tag) {

        push @commands, q{--in_occ} . $SPACE . $in_allele_count_tag;
    }
    if ($out_allele_count_tag) {

        push @commands, q{--out_occ} . $SPACE . $out_allele_count_tag;
    }
    if ($in_frequency_tag) {

        push @commands, q{--in_frq} . $SPACE . $in_frequency_tag;
    }
    if ($out_frequency_tag) {

        push @commands, q{--out_frq} . $SPACE . $out_frequency_tag;
    }
    if ($bedpedb_path) {

        push @commands, q{--bedpedb} . $SPACE . $bedpedb_path;
    }
    if ($dbfile_path) {

        push @commands, q{--db} . $SPACE . $dbfile_path;
    }

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
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
