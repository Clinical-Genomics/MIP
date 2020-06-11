package MIP::Check::File;

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
use List::Util qw{ any };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $SPACE $TAB $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.08;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_file_md5sum check_mip_process_files };
}

## Constants
Readonly my $MAX_RANDOM_NUMBER => 10_000;

sub check_file_md5sum {

## Function : Check file integrity of file using md5sum
## Returns  :
## Arguments: $check_method  => Method to perform file check (undef or md5sum)
##          : $filehandle    => Filehandle to write to
##          : $md5_file_path => File to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $check_method;
    my $filehandle;
    my $md5_file_path;

    my $tmpl = {
        check_method => {
            allow       => [ undef, qw{ md5sum } ],
            store       => \$check_method,
            strict_type => 1,
        },
        filehandle => { defined => 1, required => 1, store => \$filehandle, },
        md5_file_path =>
          { defined => 1, required => 1, store => \$md5_file_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Gnu::Coreutils qw{ gnu_md5sum gnu_rm };
    use MIP::Parse::File qw{ parse_file_suffix };

    ## Skip file
    return if ( not defined $check_method );

    my $random_integer = int rand $MAX_RANDOM_NUMBER;
    ## Parse file suffix in filename.suffix(.gz).
    ## Removes suffix if matching else return undef
    my $file_path = parse_file_suffix(
        {
            file_name   => $md5_file_path,
            file_suffix => $DOT . q{md5},
        }
    );

    #md5_file_path did not have a ".md5" suffix - skip
    return if ( not defined $file_path );

    my $md5sum_check_file = q{md5sum_check} . $UNDERSCORE . $random_integer . q{.txt};
    _write_md5sum_check_file(
        {
            filehandle        => $filehandle,
            file_path         => $file_path,
            md5sum_check_file => $md5sum_check_file,
            md5_file_path     => $md5_file_path,
        }
    );

    ## Perform md5sum check
    gnu_md5sum(
        {
            check       => 1,
            filehandle  => $filehandle,
            infile_path => $md5sum_check_file,
        }
    );
    say {$filehandle} $NEWLINE;

    ## Clean-up
    gnu_rm(
        {
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $md5sum_check_file,
        }
    );
    say {$filehandle} $NEWLINE;
    return 1;
}

sub check_mip_process_files {

## Function : Test all file that are supposed to exists after process
## Returns  :
## Arguments: $filehandle => Filehandle to write to
##          : $paths_ref  => Paths to files to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $paths_ref;

    my $tmpl = {
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
        paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$paths_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Create bash array
    print {$filehandle} q?readonly FILES=(?;

  PATH:
    foreach my $path ( @{$paths_ref} ) {

        ## First analysis and dry run will otherwise cause try to print uninitialized values
        next PATH if ( not defined $path );

        ## Add to array
        print {$filehandle} q?"? . $path . q?" ?;
    }

    ## Close bash array
    say {$filehandle} q?)?;

    ## Loop over files
    say {$filehandle} q?for file in "${FILES[@]}"?;

    ## For each element in array do
    say {$filehandle} q?do? . $SPACE;

    ## File exists and is larger than zero
    say {$filehandle} $TAB . q?if [ -s "$file" ]; then?;

    ## Echo
    say {$filehandle} $TAB x 2 . q?echo "Found file $file"?;
    say {$filehandle} $TAB . q?else?;

    ## Redirect to STDERR
    say {$filehandle} $TAB x 2 . q?echo "Could not find $file" >&2?;

    ## Set status flagg so that perl notFinished remains in sample_info_file
    say {$filehandle} $TAB x 2 . q?STATUS="1"?;
    say {$filehandle} $TAB . q?fi?;
    say {$filehandle} q?done ?, $NEWLINE;

    return 1;
}

sub _write_md5sum_check_file {

## Function : Write md5sum file
## Returns  :
## Arguments: $filehandle        => Filehandle to write to
##          : $file_path         => Outfile path
##          : $md5sum_check_file => Md5sum check file
##          : $md5_file_path     => File to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $file_path;
    my $md5sum_check_file;
    my $md5_file_path;

    my $tmpl = {
        filehandle => { defined => 1, required => 1, store => \$filehandle, },
        file_path  => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        md5sum_check_file => {
            defined     => 1,
            required    => 1,
            store       => \$md5sum_check_file,
            strict_type => 1,
        },
        md5_file_path =>
          { defined => 1, required => 1, store => \$md5_file_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Build perl regexp
    # Execute perl
    my $perl_regexp = q?perl -nae '?;

    # Print first md5sum from $md5_file_path with 2 white spaces
    $perl_regexp .= q?print $F[0].q{  ?;

    # Print file name in the same line that correspond to md5sum hash
    $perl_regexp .= $file_path . q?} ' ?;

    ## Write perl command to create md5sum check file
    print {$filehandle} $perl_regexp . $md5_file_path . q{ > } . $md5sum_check_file;
    say   {$filehandle} $NEWLINE;
    return;
}

1;
