package MIP::Program::Mip;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ mip_qccollect mip_vcfparser mip_vercollect };
}

sub mip_qccollect {

## Function : Perl wrapper for qcCollect. Collects metrics information from each analysis run.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $log_file_path          => Log file path
##          : $outfile_path           => Outfile path
##          : $regexp_file_path       => Regular expression file
##          : $skip_evaluation        => Skip evaluation step
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $log_file_path;
    my $outfile_path;
    my $regexp_file_path;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $skip_evaluation;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        log_file_path => { store => \$log_file_path, strict_type => 1, },
        outfile_path  => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        regexp_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$regexp_file_path,
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
        skip_evaluation => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$skip_evaluation,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ mip qccollect};

    ## Options
    if ($log_file_path) {

        push @commands, q{--log_file} . $SPACE . $log_file_path;
    }

    if ($regexp_file_path) {

        push @commands, q{--regexp_file} . $SPACE . $regexp_file_path;
    }

    if ($skip_evaluation) {

        push @commands, q{--skip_evaluation};
    }

    ## Infile
    if ($infile_path) {

        push @commands, q{--sample_info_file} . $SPACE . $infile_path;
    }

    ## Outfile
    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

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

sub mip_vcfparser {

## Function : Perl wrapper for MIPs vcfparser to separate clinical variants from research
## Returns  : @commands
## Arguments: $filehandle                            => Filehandle to write to
##          : $infile_path                           => Infile path
##          : $log_file_path                         => Log file path
##          : $padding                               => Pad each gene with X number of nucleotides
##          : $parse_vep                             => Parse VEP transcript specific entries
##          : $per_gene                              => Output most severe consequence transcript
##          : $pli_values_file_path                  => Pli value file path
##          : $range_feature_annotation_columns_ref  => Range file annotation columns
##          : $range_feature_file_path               => Path to range feature file
##          : $select_feature_annotation_columns_ref => Range file annotation columns
##          : $select_feature_file_path              => Path to range feature file
##          : $select_feature_matching_column        => Select feature matching column
##          : $select_outfile                        => Select outfile
##          : $stderrfile_path                       => Stderrfile path
##          : $stderrfile_path_append                => Append stderr info to file path
##          : $stdoutfile_path                       => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $log_file_path;
    my $pli_values_file_path;
    my $range_feature_annotation_columns_ref;
    my $range_feature_file_path;
    my $select_feature_annotation_columns_ref;
    my $select_feature_file_path;
    my $select_feature_matching_column;
    my $select_outfile;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $padding;
    my $parse_vep;
    my $per_gene;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        log_file_path => {
            store       => \$log_file_path,
            strict_type => 1,
        },
        padding => {
            allow       => [ undef, qr{ \A\d+\z }sxm, ],
            store       => \$padding,
            strict_type => 1,
        },
        parse_vep => {
            allow       => [ undef, 0, 1, 2 ],
            default     => 0,
            store       => \$parse_vep,
            strict_type => 1,
        },
        per_gene => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$per_gene,
            strict_type => 1,
        },
        pli_values_file_path => {
            store       => \$pli_values_file_path,
            strict_type => 1,
        },
        range_feature_annotation_columns_ref => {
            default     => [],
            store       => \$range_feature_annotation_columns_ref,
            strict_type => 1,
        },
        range_feature_file_path =>
          { store => \$range_feature_file_path, strict_type => 1, },
        select_feature_annotation_columns_ref => {
            default     => [],
            store       => \$select_feature_annotation_columns_ref,
            strict_type => 1,
        },
        select_feature_file_path =>
          { store => \$select_feature_file_path, strict_type => 1, },
        select_feature_matching_column => {
            allow       => [ undef, qr{ \A\d+\z }sxm, ],
            store       => \$select_feature_matching_column,
            strict_type => 1,
        },
        select_outfile  => { store => \$select_outfile, strict_type => 1, },
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

    ## Stores commands depending on input parameters
    my @commands = qw{ mip vcfparser };

    ## Infile
    push @commands, $infile_path;

    ## Options
    if ($log_file_path) {

        push @commands, q{--log_file} . $SPACE . $log_file_path;
    }

    if ($parse_vep) {

        push @commands, q{--parse_vep};
    }

    if ($per_gene) {

        push @commands, q{--per_gene};
    }

    if ( defined $padding ) {

        push @commands, q{--padding} . $SPACE . $padding;
    }
    if ( defined $pli_values_file_path ) {

        push @commands, q{--pli_values_file} . $SPACE . $pli_values_file_path;
    }

    if ($range_feature_file_path) {

        push @commands, q{--range_feature_file} . $SPACE . $range_feature_file_path;
    }

    # Limit output to regions
    if ( @{$range_feature_annotation_columns_ref} ) {

        push @commands,
          q{--range_feature_annotation_columns} . $SPACE . join $COMMA,
          @{$range_feature_annotation_columns_ref};
    }

    if ($select_feature_file_path) {

        push @commands, q{--select_feature_file} . $SPACE . $select_feature_file_path;
    }

    if ($select_feature_matching_column) {

        push @commands,
          q{--select_feature_matching_column} . $SPACE . $select_feature_matching_column;
    }

    # Limit output to regions
    if ( @{$select_feature_annotation_columns_ref} ) {

        push @commands,
          q{--select_feature_annotation_columns} . $SPACE . join $COMMA,
          @{$select_feature_annotation_columns_ref};
    }

    if ($select_outfile) {

        push @commands, q{--select_outfile} . $SPACE . $select_outfile;
    }

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

sub mip_vercollect {

## Function : Perl wrapper for vercollect
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $log_file_path          => Log file path
##          : $outfile_path           => Outfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $log_file_path;
    my $outfile_path;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        log_file_path => { store => \$log_file_path, strict_type => 1, },
        outfile_path  => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ mip vercollect};

    ## Options
    if ($log_file_path) {

        push @commands, q{--log_file} . $SPACE . $log_file_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, q{--infile} . $SPACE . $infile_path;
    }

    ## Outfile
    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

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
