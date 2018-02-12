package MIP::Program::Variantcalling::Mip_vcfparser;

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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ mip_vcfparser };
}

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

sub mip_vcfparser {

## Function : Perl wrapper for generic commands module.
## Returns  : @commands
## Arguments: $FILEHANDLE                            => Filehandle to write to
##          : $infile_path                           => Infile path
##          : $padding                               => Pad each gene with X number of nucleotides
##          : $parse_vep                             => Parse VEP transcript specific entries
##          : $per_gene                              => Output most severe consequence transcript
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
    my $FILEHANDLE;
    my $infile_path;
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
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path,
        },
        padding => {
            allow       => [ undef, qr/^\d+$/ ],
            strict_type => 1,
            store       => \$padding,
        },
        parse_vep => {
            default     => 0,
            allow       => [ undef, 0, 1, 2 ],
            strict_type => 1,
            store       => \$parse_vep,
        },
        per_gene => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$per_gene,
        },
        range_feature_annotation_columns_ref => {
            default     => [],
            strict_type => 1,
            store       => \$range_feature_annotation_columns_ref,
        },
        range_feature_file_path =>
          { strict_type => 1, store => \$range_feature_file_path, },

        select_feature_annotation_columns_ref => {
            default     => [],
            strict_type => 1,
            store       => \$select_feature_annotation_columns_ref,
        },
        select_feature_file_path =>
          { strict_type => 1, store => \$select_feature_file_path, },
        select_feature_matching_column => {
            allow       => [ undef, qr/^\d+$/ ],
            strict_type => 1,
            store       => \$select_feature_matching_column,
        },
        select_outfile  => { strict_type => 1, store => \$select_outfile, },
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

    ## Stores commands depending on input parameters
    my @commands = q{vcfparser};

    ## Infile
    push @commands, $infile_path;

    ## Options
    if ($parse_vep) {

        push @commands, q{--parse_vep};
    }

    if ($per_gene) {

        push @commands, q{--per_gene};
    }

    if ( $padding ) {

        push @commands, q{--padding} . $SPACE . $padding;
    }

    if ($range_feature_file_path) {

        push @commands,
          q{--range_feature_file} . $SPACE . $range_feature_file_path;
    }

    # Limit output to regions
    if ( @{$range_feature_annotation_columns_ref} ) {

        push @commands,
          q{--range_feature_annotation_columns} . $SPACE . join $COMMA,
          @{$range_feature_annotation_columns_ref};
    }

    if ($select_feature_file_path) {

        push @commands,
          q{--select_feature_file} . $SPACE . $select_feature_file_path;
    }

    if ($select_feature_matching_column) {

        push @commands,
            q{--select_feature_matching_column}
          . $SPACE
          . $select_feature_matching_column;
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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
