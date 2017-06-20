package MIP::Test::Writefile;

#### Copyright 2017 Henrik Stranneheim

use strict;
use warnings;
use warnings qw(FATAL utf8);
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);
use autodie;
use Test::More;

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(test_write_to_file);
}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

sub test_write_to_file {

## test_write_to_file

##Function : Test of writing to file using anonymous FILEHANDLE
##Returns  : ""
##Arguments: $module_function_cref, $args_ref, $base_command, $separator
##         : $module_function_cref => Module method to test
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

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    # Add new FILEHANDLE to args
    push @{$args_ref}, 'FILEHANDLE', $FILEHANDLE;

    # For storing info to write
    my $file_content;

    ## Store file content in memory by using referenced variable
    open $FILEHANDLE, '>', \$file_content
      or croak 'Cannot write to ' . $file_content . ': ' . $OS_ERROR;

    $module_function_cref->( { @{$args_ref} } );

    close $FILEHANDLE;

    ## Perform test
    ok( $file_content =~ /^$base_command/,
        q{Write commands with '} . $separator . q{' to file} );

    return;
}

1;
