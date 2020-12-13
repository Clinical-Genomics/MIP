package MIP::Validate::Case;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $DOT $LOG_NAME $SINGLE_QUOTE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_infiles check_infile_contain_sample_id check_sample_ids };
}

sub check_infiles {

## Function : Check infiles found and that they contain sample_id
## Returns  :
## Arguments: $infiles_ref      => Infiles to check {REF}
##          : $infile_directory => Infile directory
##          : $sample_id        => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infiles_ref;
    my $infile_directory;
    my $sample_id;

    my $tmpl = {
        infile_directory => {
            defined     => 1,
            required    => 1,
            store       => \$infile_directory,
            strict_type => 1,
        },
        infiles_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infiles_ref,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## No fastq infiles
    if ( not @{$infiles_ref} ) {

        $log->fatal(
qq{Could not find any fastq files in supplied infiles directory: $infile_directory}
        );
        exit 1;
    }

    ## Check that infiledirs/infile contains sample_id in filename
  INFILE:
    foreach my $infile ( @{$infiles_ref} ) {

        next INFILE if ( $infile =~ /$sample_id/sxm );

        $log->fatal( q{Could not detect sample_id: }
              . $sample_id
              . q{ in supplied infile: }
              . catfile( $infile_directory, $infile ) );
        $log->fatal(
q{Check that: '--sample_ids' and '--infile_dirs' contain the same sample_id and that the filename of the infile contains the sample_id.}
        );
        exit 1;
    }
    return 1;
}

sub check_infile_contain_sample_id {

## Function : Check that the sample_id provided and sample_id in infile name match
## Returns  : 1
## Arguments: $infile_name      => Infile name
##          : $infile_sample_id => Sample_id collect with regexp from infile
##          : $sample_id        => Sample id from user

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_name;
    my $infile_sample_id;
    my $sample_id;

    my $tmpl = {
        infile_name => {
            defined     => 1,
            required    => 1,
            store       => \$infile_name,
            strict_type => 1,
        },
        infile_sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$infile_sample_id,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    return 1 if ( $sample_id eq $infile_sample_id );

    $log->fatal(
qq{$sample_id supplied and sample_id $infile_sample_id found in file : $infile_name does not match}
    );
    $log->fatal(qq{Please rename file to match sample_id: $sample_id});
    exit 1;
}

sub check_sample_ids {

## Function : Check that the case_id and the sample_id(s) exists and are unique. Check if id sample_id contains "_".
## Returns  : 1
## Arguments: $case_id        => Case id
##          : $sample_ids_ref => Sample ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $sample_ids_ref;

    my $tmpl = {
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Hash to test duplicate sample_ids later
    my %seen;

    if ( not @{$sample_ids_ref} ) {

        $log->fatal(q{Please provide sample_id(s)});
        exit 1;
    }

  SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        ## Increment instance to check duplicates later
        $seen{$sample_id}++;

        ## Family_id cannot be the same as sample_id
        if ( $case_id eq $sample_id ) {

            $log->fatal( q{Case_id: '}
                  . $case_id
                  . q{' equals sample_id: '}
                  . $sample_id
                  . $SINGLE_QUOTE );
            $log->fatal(q{Please make sure that the case_id and sample_id(s) are unique});
            exit 1;
        }
        ## Check for unique sample_ids
        if ( $seen{$sample_id} > 1 ) {

            $log->fatal( q{Sample_id: '} . $sample_id . q{' is not uniqe.} );
            exit 1;
        }
        ## Sample_id contains "_", not allowed in filename convention
        if ( $sample_id =~ /$UNDERSCORE/sxm ) {

            $log->fatal( q{Sample_id: '} . $sample_id . q{' contains '_'} );
            $log->fatal(
q{Please rename sample_id according to MIP's filename convention, removing the '_'.}
            );
            exit 1;
        }
    }
    return 1;
}

1;
