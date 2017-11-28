package MIP::Test::Writefile;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use Test::More;

## CPANM
use Readonly;
use autodie;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ test_write_to_file };
}

## Constants
Readonly my $SPACE => q{ };

sub test_write_to_file {

##Function : Test of writing to file using anonymous FILEHANDLE
##Returns  :
##Arguments: $module_function_cref => Module method to test
##         : $args_ref             => Arguments to function call
##         : $base_command         => First word in command line usually name of executable
##         : $separator            => Separator to use when writing

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $module_function_cref;
    my $args_ref;
    my $base_command;
    my $separator;

    my $tmpl = {
        module_function_cref =>
          { required => 1, defined => 1, store => \$module_function_cref },
        args_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$args_ref
        },
        base_command => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$base_command
        },
        separator => {
            default     => q{ },
            strict_type => 1,
            store       => \$separator
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    # Add new FILEHANDLE to args
    push @{$args_ref}, q{FILEHANDLE}, $FILEHANDLE;

    # For storing info to write
    my $file_content;

    ## Store file content in memory by using referenced variable
    open $FILEHANDLE, q{>}, \$file_content
      or croak q{Cannot write to}
      . $SPACE
      . $file_content . q{:}
      . $SPACE
      . $OS_ERROR;

    $module_function_cref->( { @{$args_ref} } );

    close $FILEHANDLE;

    ## Perform test
    my $return_base_command;
    if ( $file_content =~ /^($base_command)/ ) {

        $return_base_command = $1;
    }
    is( $return_base_command, $base_command,
        q{Write commands with '} . $separator . q{' to file} );
    return;
}

1;
