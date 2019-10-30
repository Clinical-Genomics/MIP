package MIP::Test::Writefile;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use Test::More;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COLON $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ test_write_to_file write_toml_config };
}

sub test_write_to_file {

## Function : Test of writing to file using anonymous filehandle
## Returns  :
## Arguments: $args_ref             => Arguments to function call
##          : $base_commands_ref    => Base commands {REF}
##          : $module_function_cref => Module method to test
##          : $separator            => Separator to use when writing

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $args_ref;
    my $base_commands_ref;
    my $module_function_cref;
    my $separator;

    my $tmpl = {
        args_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$args_ref,
            strict_type => 1,
        },
        base_commands_ref => {
            defined     => 1,
            default     => [],
            required    => 1,
            store       => \$base_commands_ref,
            strict_type => 1,
        },
        module_function_cref =>
          { defined => 1, required => 1, store => \$module_function_cref, },
        separator => {
            default     => $SPACE,
            store       => \$separator,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    # Add new filehandle to args
    push @{$args_ref}, q{filehandle}, $filehandle;

    # For storing info to write
    my $file_content;

    ## Store file content in memory by using referenced variable
    open $filehandle, q{>}, \$file_content
      or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

    $module_function_cref->( { @{$args_ref} } );

    close $filehandle;

    my $base_command = join $separator, @{$base_commands_ref};

    ## Perform test
    my ($returned_base_command) = $file_content =~ /^($base_command)/ms;

    is( $returned_base_command, $base_command,
        q{Write commands with '} . $separator . q{' to file} );
    return;
}

sub write_toml_config {

## Function : Copy toml template and update to system specific path
## Returns  :
## Arguments: $test_reference_path => Test reference path
##          : $toml_config_path    => Path to new toml config
##          : $toml_template_path  => Path to toml template

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $test_reference_path;
    my $toml_config_path;
    my $toml_template_path;

    my $tmpl = {
        test_reference_path => {
            store       => \$test_reference_path,
            strict_type => 1,
        },
        toml_config_path => {
            store       => \$toml_config_path,
            strict_type => 1,
        },
        toml_template_path => {
            store       => \$toml_template_path,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use File::Copy qw{ cp };
    use Path::Tiny qw{ path };

    ## Create copy of template
    cp( $toml_template_path, $toml_config_path );

    ## Create Path::Tiny object
    $toml_config_path = path($toml_config_path);

    ## Read into memory
    my $toml_config = $toml_config_path->slurp_utf8;

    ## Replace with system specific path
    $toml_config =~ s/TEST_REFERENCES!/$test_reference_path/xms;

    ## Write to file
    $toml_config_path->spew_utf8($toml_config);

    return;
}

1;
