package MIP::Qc_data;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE $NEWLINE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      add_qc_data_evaluation_info
      add_qc_data_recipe_info
      add_qc_data_regexp_return
      get_qc_data_case_recipe_attributes
      get_qc_data_sample_recipe_attributes
      get_regexp_qc_data
      set_qc_data_recipe_info
    };
}

sub add_qc_data_evaluation_info {

## Function : Add recipe evaluation info in qc_data hash
## Returns  :
## Arguments: $qc_data_href => Qc_data hash {REF}
##          : $recipe_name  => Recipe to set attributes for
##          : $value        => Value to store

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $recipe_name;
    my $value;

    my $tmpl = {
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        value => {
            defined     => 1,
            required    => 1,
            store       => \$value,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add recipe key value pair for arbitrary info on case level
    push @{ $qc_data_href->{evaluation}{$recipe_name} }, $value;
    return;
}

sub add_qc_data_recipe_info {

## Function : Add recipe arbitrary info in qc_data hash
## Returns  :
## Arguments: $infile       => Infile key
##          : $key          => Metafile tag
##          : $qc_data_href => Qc_data hash {REF}
##          : $recipe_name  => Recipe to set attributes for
##          : $sample_id    => Sample ID
##          : $value        => Value to store

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile;
    my $key;
    my $qc_data_href;
    my $recipe_name;
    my $sample_id;
    my $value;

    my $tmpl = {
        infile => {
            store       => \$infile,
            strict_type => 1,
        },
        key => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$key,
        },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        value => {
            defined     => 1,
            required    => 1,
            store       => \$value,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set recipe key value pair for arbitrary info on sample and infile level
    if ( $sample_id and $infile ) {

        ## Add to array
        push @{ $qc_data_href->{sample}{$sample_id}{$infile}{$recipe_name}{$key} },
          $value;
        return;
    }

    ## Add recipe key value pair for arbitrary info on case level
    push @{ $qc_data_href->{recipe}{$recipe_name}{$key} }, $value;
    return;
}

sub add_qc_data_regexp_return {

## Function  : Use reg exp to collect data via system call and split potential return
##             by seperator
## Returns   : 1 or undef
## Arguments : $data_file_path => Path to data file from which to collect
##           : $qc_href        => Save header(s) or data from each outfile {REF}
##           : $recipe_name    => The recipe to examine
##           : $regexp         => Regular expression to collect data
##           : $regexp_key     => Regexp key to store data in

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $data_file_path;
    my $qc_href;
    my $recipe_name;
    my $regexp;
    my $regexp_key;

    my $tmpl = {
        data_file_path => { store => \$data_file_path, strict_type => 1, },
        qc_href        => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        regexp => {
            defined     => 1,
            required    => 1,
            store       => \$regexp,
            strict_type => 1,
        },
        regexp_key => {
            defined     => 1,
            required    => 1,
            store       => \$regexp_key,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @regexp_returns = get_regexp_qc_data(
        {
            data_file_path => $data_file_path,
            regexp         => $regexp,
        }
    );

    ## Covers both whitespace and tab. Add other separators if required
    my @separators = ( qw{ \s+ ! }, q{,} );

    ## Loop through possible separators to seperate any eventual header/data elements
  SEPARATOR:
    foreach my $separator (@separators) {

        ## Add to qc_data
        @{ $qc_href->{$recipe_name}{$regexp_key} } = split /$separator/sxm,
          join $NEWLINE, @regexp_returns;

        ## Return true if seperation of data was successful
        return 1
          if ( defined $qc_href->{$recipe_name}{$regexp_key} );
    }
    return;
}

sub get_qc_data_case_recipe_attributes {

## Function : Get case recipe attributes from qc_data hash
## Returns  : "$attribute" or "attributes_ref" or "$attribute_href"
## Arguments: $attribute    => Attribute key
##          : $qc_data_href => Sample info hash {REF}
##          : $recipe_name  => Recipe to get attributes from

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $qc_data_href;
    my $recipe_name;

    my $tmpl = {
        attribute => {
            store       => \$attribute,
            strict_type => 1,
        },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get and return attribute value
    if ( defined $attribute && $attribute ) {

        return $qc_data_href->{recipe}{$recipe_name}{$attribute};
    }

    ## Return recipe attribute hash
    return %{ $qc_data_href->{recipe}{$recipe_name} };
}

sub get_qc_data_sample_recipe_attributes {

## Function : Get sample recipe attributes from qc_data hash
## Returns  : "$attribute" or "$attribute_href"
## Arguments: $attribute    => Attribute key
##          : $infile       => Infile key
##          : $qc_data_href => Sample info hash {REF}
##          : $recipe_name  => Recipe to get attributes from
##          : $sample_id    => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $infile;
    my $qc_data_href;
    my $recipe_name;
    my $sample_id;

    my $tmpl = {
        attribute => {
            store       => \$attribute,
            strict_type => 1,
        },
        infile => {
            defined     => 1,
            required    => 1,
            store       => \$infile,
            strict_type => 1,
        },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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

    ## Get and return attribute value
    if ( defined $attribute && $attribute ) {

        return $qc_data_href->{sample}{$sample_id}{$infile}{$recipe_name}{$attribute};
    }

    ## Get recipe attribute hash
    return %{ $qc_data_href->{sample}{$sample_id}{$infile}{$recipe_name} };
}

sub get_regexp_qc_data {

    ## Function  : Use reg exp to collect qc data via system call
    ## Returns   : @{ $chld_handler{output} } or @{ $chld_handler{error} }
    ## Arguments : $data_file_path => Path to data file from which to collect
    ##           : $regexp         => Regular expression to collect data

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $data_file_path;
    my $regexp;

    my $tmpl = {
        data_file_path => { store => \$data_file_path, strict_type => 1, },
        regexp         => {
            defined     => 1,
            required    => 1,
            store       => \$regexp,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Unix::System qw{ system_cmd_call };

    ## Get return from reg exp system call
    my %chld_handler =
      system_cmd_call( { command_string => qq{$regexp $data_file_path}, } );

    ## Print stderr if returned from regexp
    if ( @{ $chld_handler{error} } ) {

        ## Be verbose that something went wrong
        say {*STDERR} join $NEWLINE, @{ $chld_handler{error} };
        return @{ $chld_handler{error} };
    }
    return @{ $chld_handler{output} };
}

sub set_qc_data_recipe_info {

## Function : Set recipe arbitrary info in qc_data hash
## Returns  :
## Arguments: $infile       => Infile key
##          : $key          => Metafile tag
##          : $qc_data_href => Qc_data hash {REF}
##          : $recipe_name  => Recipe to set attributes for
##          : $sample_id    => Sample ID
##          : $value        => Value to store

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile;
    my $key;
    my $qc_data_href;
    my $recipe_name;
    my $sample_id;
    my $value;

    my $tmpl = {
        infile => {
            store       => \$infile,
            strict_type => 1,
        },
        key => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$key,
        },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
        recipe_name => {
            store       => \$recipe_name,
            strict_type => 1,
        },
        value => {
            store       => \$value,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not defined $value );

    ## Set recipe key value pair for arbitrary info on sample and infile level
    if ( $sample_id and $infile ) {

        $qc_data_href->{sample}{$sample_id}{$infile}{$recipe_name}{$key} = $value;
        return;
    }
    if ( $sample_id and $recipe_name ) {

        $qc_data_href->{sample}{$sample_id}{$recipe_name}{$key} = $value;
        return;
    }
    if ($sample_id) {

        $qc_data_href->{sample}{$sample_id}{$key} = $value;
        return;
    }

    $qc_data_href->{recipe}{$recipe_name}{$key} = $value;

    return;
}

1;
