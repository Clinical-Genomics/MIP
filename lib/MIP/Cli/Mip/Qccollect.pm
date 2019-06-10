package MIP::Cli::Mip::Qccollect;

use 5.026;
use Carp;
use Cwd;
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Moose::Util::TypeConstraints;
use MooseX::App::Command;
use MooseX::Types::Moose qw{ ArrayRef Bool HashRef Int Str };

## MIPs lib/
use MIP::Constants qw{ $NEWLINE };

our $VERSION = 1.00;

command_short_description(q{MIP qccollect command});

command_long_description(q{Entry point for collecting MIP QC metrics});

command_usage(
    q{qccollect <options> -si [sample_info.yaml] -r [regexp.yaml] -o [outfile]});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {

    ## Input from Cli
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    # Flatten argument(s)
    my $evaluate_plink_gender = $arg_href->{evaluate_plink_gender};
    my $print_regexp          = $arg_href->{print_regexp};
    my $regexp_file           = $arg_href->{regexp_file};
    my $sample_info_file      = $arg_href->{sample_info_file};
    my $skip_evaluation       = $arg_href->{skip_evaluation};
    my $log_file              = $arg_href->{log_file};
    my $print_regexp_outfile  = $arg_href->{print_regexp_outfile};
    my $outfile               = $arg_href->{outfile};

    use MIP::Log::MIP_log4perl qw{ initiate_logger };
    use MIP::Main::Qccollect qw{ mip_qccollect };
    use MIP::Qcc_regexp qw{ regexp_to_yaml };

    ## Creates log object
    my $log = initiate_logger(
        {
            file_path => $log_file,
            log_name  => uc q{mip_qccollect},
        }
    );

    ## Write default regexp to YAML if demanded
    regexp_to_yaml(
        {
            log                  => $log,
            print_regexp_outfile => $print_regexp_outfile,
        }
    );

    mip_qccollect(
        {
            evaluate_plink_gender => $evaluate_plink_gender,
            regexp_file           => $regexp_file,
            sample_info_file      => $sample_info_file,
            skip_evaluation       => $skip_evaluation,
            log                   => $log,
            outfile               => $outfile,
        }
    );
    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    my ($arg_href) = @_;

    option(
        q{evaluate_plink_gender} => (
            cmd_aliases   => [qw{ epg }],
            documentation => q{Evaluate plink gender},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{log_file} => (
            cmd_aliases   => [qw{ l }],
            default       => catfile( cwd(), q{qccollect.log} ),
            documentation => q{Log file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{outfile} => (
            cmd_aliases   => [qw{ o }],
            cmd_tags      => [q{YAML}],
            default       => q{qcmetrics.yaml},
            documentation => q{Data file output},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{print_regexp_outfile} => (
            cmd_aliases   => [qw{ prego }],
            cmd_flag      => q{regexp_outfile},
            cmd_tags      => [q{YAML}],
            documentation => q{Regexp YAML outfile},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{regexp_file} => (
            cmd_aliases   => [qw{ rxf }],
            cmd_tags      => [q{YAML}],
            documentation => q{Regular expression file path},
            is            => q{rw},
            isa           => Str,
            required      => 1,
        )
    );

    option(
        q{sample_info_file} => (
            cmd_aliases   => [qw{ sif }],
            cmd_tags      => [q{YAML}],
            documentation => q{File for sample info used in the analysis},
            is            => q{rw},
            isa           => Str,
            required      => 1,
        )
    );

    option(
        q{skip_evaluation} => (
            cmd_aliases   => [qw{ ske }],
            documentation => q{Skip evaluation step},
            is            => q{rw},
            isa           => Bool,
        )
    );
    return;
}

1;
