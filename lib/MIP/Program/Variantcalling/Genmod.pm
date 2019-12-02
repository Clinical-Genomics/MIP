package MIP::Program::Variantcalling::Genmod;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DASH $NEWLINE $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ genmod_annotate genmod_compound genmod_filter genmod_models genmod_score };
}

sub genmod_annotate {

## Function : Perl wrapper for writing Genmod annotate recipe to $filehandle or return commands array. Based on genmod 3.7.0.
## Returns  : @commands
## Arguments: $annotate_region        => Annotate what regions a variant belongs to
##          : $append_stderr_info     => Append stderr info to file
##          : $cadd_file_paths_ref    => Specify the path to a bgzipped cadd file (with index) with variant scores
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $max_af                 => If the MAX AF should be annotated
##          : $outfile_path           => Outfile path to write to
##          : $spidex_file_path       => Specify the path to a bgzipped tsv file (with index) with spidex information
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $temp_directory_path    => Directory for storing intermediate files
##          : $thousand_g_file_path   => Specify the path to a bgzipped vcf file (with index) with 1000g variants
##          : $verbosity              => Increase output verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cadd_file_paths_ref;
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $spidex_file_path;
    my $temp_directory_path;
    my $thousand_g_file_path;

    ## Default(s)
    my $annotate_region;
    my $append_stderr_info;
    my $max_af;
    my $verbosity;

    my $tmpl = {
        annotate_region => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$annotate_region,
            strict_type => 1,
        },
        append_stderr_info => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$append_stderr_info,
            strict_type => 1,
        },
        cadd_file_paths_ref =>
          { default => [], store => \$cadd_file_paths_ref, strict_type => 1, },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        max_af => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$max_af,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path, strict_type => 1, },
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
        spidex_file_path     => { store => \$spidex_file_path,     strict_type => 1, },
        temp_directory_path  => { store => \$temp_directory_path,  strict_type => 1, },
        thousand_g_file_path => { store => \$thousand_g_file_path, strict_type => 1, },
        verbosity            => {
            allow       => qr/ \A \w+ \z /sxm,
            store       => \$verbosity,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Genmod annotate
    my @commands = qw{ genmod };

    ## Options
    if ($verbosity) {

        push @commands, $DASH . $verbosity;
    }

    ## Genmod sub command
    push @commands, q{annotate};

    if ($temp_directory_path) {

        push @commands, q{--temp_dir} . $SPACE . $temp_directory_path;
    }
    if ($annotate_region) {

        push @commands, q{--annotate_regions};
    }
    if ($thousand_g_file_path) {

        push @commands, q{--thousand-g} . $SPACE . $thousand_g_file_path;
    }
    if ($spidex_file_path) {

        push @commands, q{--spidex} . $SPACE . $spidex_file_path;
    }
    if ( @{$cadd_file_paths_ref} ) {

        push @commands,
          q{--cadd_file} . $SPACE . join $SPACE . q{--cadd_file} . $SPACE,
          @{$cadd_file_paths_ref};
    }
    if ($max_af) {

        push @commands, q{--max_af};
    }
    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

    ## Infile
    push @commands, $infile_path;

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

sub genmod_compound {

## Function : Perl wrapper for writing Genmod compound recipe to $filehandle or return commands array. Based on genmod 3.7.0.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $temp_directory_path    => Directory for storing intermediate files
##          : $thread_number          => Define how many processes that should be use for annotation
##          : $verbosity              => Increase output verbosity
##          : $vep                    => If variants are annotated with the Variant Effect Predictor

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $temp_directory_path;

    ## Default(s)
    my $thread_number;
    my $verbosity;
    my $vep;

    my $tmpl = {
        filehandle  => { store => \$filehandle },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path, strict_type => 1, },
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
        temp_directory_path => { store => \$temp_directory_path, strict_type => 1, },
        thread_number       => {
            allow       => qr/ ^\d+$ /sxm,
            default     => 0,
            store       => \$thread_number,
            strict_type => 1,
        },
        verbosity => {
            allow       => qr/ ^\w+$ /sxm,
            store       => \$verbosity,
            strict_type => 1,
        },
        vep => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$vep,
            strict_type => 1,
        },
        verbosity => {
            allow       => qr/ \A \w+ \z /sxm,
            store       => \$verbosity,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Genmod compound
    my @commands = qw{ genmod };

    ## Options
    if ($verbosity) {

        push @commands, q{-} . $verbosity;
    }
    push @commands, q{compound};

    if ($temp_directory_path) {

        push @commands, q{--temp_dir} . $SPACE . $temp_directory_path;
    }
    if ($thread_number) {

        push @commands, q{--processes} . $SPACE . $thread_number;
    }
    if ($vep) {

        push @commands, q{--vep};
    }
    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

    ## Infile
    push @commands, $infile_path;

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

sub genmod_filter {

## Function : Perl wrapper for writing Genmod filter recipe to $filehandle or return commands array. Based on genmod 3.7.0.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $threshold              => Directory for storing intermediate files
##          : $verbosity              => Increase output verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $threshold;

    ## Default(s)
    my $append_stderr_info;
    my $verbosity;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path, strict_type => 1, },
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
        threshold => {
            allow       => qr/ ^\d+$ | ^\d+.\d+$ /xsm,
            store       => \$threshold,
            strict_type => 1,
        },
        verbosity => {
            allow       => qr/ \A \w+ \z /sxm,
            store       => \$verbosity,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Genmod filter
    my @commands = qw{ genmod };

    ## Options
    if ($verbosity) {

        push @commands, q{-} . $verbosity;
    }
    push @commands, q{filter};

    if ($threshold) {

        push @commands, q{--threshold} . $SPACE . $threshold;
    }
    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

    ## Infile
    push @commands, $infile_path;

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

sub genmod_models {

## Function : Perl wrapper for writing Genmod models recipe to $filehandle or return commands array. Based on genmod 3.7.0.
## Returns  : @commands
## Arguments: $case_file                    => Family file
##          : $case_type                    => Setup of family file
##          : $filehandle                   => Filehandle to write to
##          : $infile_path                  => Infile path to read from
##          : $outfile_path                 => Outfile path to write to
##          : $reduced_penetrance_file_path => Reduced penetrance file path
##          : $stderrfile_path              => Stderrfile path
##          : $stderrfile_path_append       => Append stderr info to file path
##          : $stdoutfile_path              => Stdoutfile path
##          : $temp_directory_path          => Directory for storing intermediate files
##          : $thread_number                => Define how many processes that should be use for annotation
##          : $vep                          => If variants are annotated with the Variant Effect Predictor
##          : $verbosity                    => Increase output verbosity
##          : $whole_gene                   => If compounds should be checked for all variants in a gene

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_file;
    my $case_type;
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $reduced_penetrance_file_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $temp_directory_path;

    ## Default(s)
    my $thread_number;
    my $verbosity;
    my $vep;
    my $whole_gene;

    my $tmpl = {
        case_file => {
            defined     => 1,
            required    => 1,
            store       => \$case_file,
            strict_type => 1,
        },
        case_type => {
            allow       => [qw{ped alt cmms mip }],
            store       => \$case_type,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        reduced_penetrance_file_path =>
          { store => \$reduced_penetrance_file_path, strict_type => 1, },
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
        temp_directory_path => { store => \$temp_directory_path, strict_type => 1, },
        thread_number       => {
            allow       => qr/ ^\d+$ /sxm,
            default     => 1,
            store       => \$thread_number,
            strict_type => 1,
        },
        verbosity => {
            allow       => qr/ ^\w+$ /sxm,
            store       => \$verbosity,
            strict_type => 1,
        },
        vep => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$vep,
            strict_type => 1,
        },
        verbosity => {
            allow       => qr/ \A \w+ \z /sxm,
            store       => \$verbosity,
            strict_type => 1,
        },
        whole_gene => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$whole_gene,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Genmod annotate
    my @commands = qw{ genmod };

    ## Options
    if ($verbosity) {

        push @commands, q{-} . $verbosity;
    }
    push @commands, q{models};

    if ($temp_directory_path) {

        push @commands, q{--temp_dir} . $SPACE . $temp_directory_path;
    }

    ## Case file
    push @commands, q{--family_file} . $SPACE . $case_file;

    if ($case_type) {

        push @commands, q{--family_type} . $SPACE . $case_type;
    }
    if ($reduced_penetrance_file_path) {

        push @commands, q{--reduced_penetrance} . $SPACE . $reduced_penetrance_file_path;
    }
    if ($thread_number) {

        push @commands, q{--processes} . $SPACE . $thread_number;
    }
    if ($vep) {

        push @commands, q{--vep};
    }
    if ($whole_gene) {

        push @commands, q{--whole_gene};
    }

    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

    ## Infile
    push @commands, $infile_path;

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

sub genmod_score {

## Function : Perl wrapper for writing Genmod score recipe to $filehandle or return commands array. Based on genmod 3.7.0.
## Returns  : @commands
## Arguments: $case_file               => Family file
##          : $case_type               => Setup of family file
##          : $filehandle              => Filehandle to write to
##          : $infile_path             => Infile path to read from
##          : $outfile_path            => Outfile path to write to
##          : $rank_model_file_path    => The plug-in config file
##          : $rank_result             => Add a info field that shows how the different categories contribute to the rank score
##          : $stderrfile_path         => Stderrfile path
##          : $stderrfile_path_append  => Append stderr info to file path
##          : $stdoutfile_path         => Stdoutfile path
##          : $temp_directory_path     => Directory for storing intermediate files
##          : $verbosity               => Increase output verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_file;
    my $case_type;
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $rank_model_file_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $temp_directory_path;

    ## Default(s)
    my $rank_result;
    my $verbosity;

    my $tmpl = {
        case_file => {
            defined     => 1,
            required    => 1,
            store       => \$case_file,
            strict_type => 1,
        },
        case_type => {
            allow       => [qw{ ped alt cmms mip }],
            store       => \$case_type,
            strict_type => 1,
        },
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path         => { store => \$outfile_path, strict_type => 1, },
        rank_model_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$rank_model_file_path,
            strict_type => 1,
        },
        rank_result => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$rank_result,
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
        temp_directory_path => { store => \$temp_directory_path, strict_type => 1, },
        verbosity           => {
            allow       => qr/ \A \w+ \z /sxm,
            store       => \$verbosity,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Genmod score
    my @commands = qw{ genmod };

    ## Options
    if ($verbosity) {

        push @commands, q{-} . $verbosity;
    }

    push @commands, q{score};

    if ($temp_directory_path) {

        push @commands, q{--temp_dir} . $SPACE . $temp_directory_path;
    }

    ## Case file
    push @commands, q{--family_file} . $SPACE . $case_file;

    if ($case_type) {

        push @commands, q{--family_type} . $SPACE . $case_type;
    }
    if ($rank_result) {

        push @commands, q{--rank_results};
    }

    ## Rank model file
    push @commands, q{--score_config} . $SPACE . $rank_model_file_path;

    if ($outfile_path) {

        push @commands, q{--outfile} . $SPACE . $outfile_path;
    }

    ## Infile
    push @commands, $infile_path;

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
            filehandle   => $filehandle,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
