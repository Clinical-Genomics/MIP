package MIP::Qc_data;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
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

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      add_qc_data_evaluation_info
      add_qc_data_recipe_info
      add_qc_data_regexp_return
      add_to_qc_data
      get_qc_data_case_recipe_attributes
      get_qc_data_sample_recipe_attributes
      get_qc_metric
      get_qc_recipe_data
      get_regexp_qc_data
      parse_qc_recipe_data
      parse_qc_recipe_table_data
      parse_regexp_hash_and_collect
      set_header_metrics_to_qc_data
      set_metrics_to_store
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
        push @{ $qc_data_href->{sample}{$sample_id}{$infile}{$recipe_name}{$key} }, $value;
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
          if ( @{ $qc_href->{$recipe_name}{$regexp_key} } );
    }
    return;
}

sub add_to_qc_data {

## Function  : Parse qc_recipe_data for data metric(s) and add metric(s) to qc_data hash to enable write to yaml format
## Returns   :
## Arguments : $attribute           => Attribute to collect for
##           : $infile              => Infile to recipe
##           : $qc_data_href        => Qc data hash {REF}
##           : $qc_header_href      => Save header(s) in each outfile {REF}
##           : $qc_recipe_data_href => Hash to save data in each outfile {REF}
##           : $recipe              => Recipe to examine
##           : $sample_id           => Sample ID
##           : $sample_info_href    => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $infile;
    my $qc_data_href;
    my $qc_header_href;
    my $qc_recipe_data_href;
    my $recipe;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        attribute => {
            defined     => 1,
            required    => 1,
            store       => \$attribute,
            strict_type => 1,
        },
        infile       => { store => \$infile, strict_type => 1, },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        qc_header_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_header_href,
            strict_type => 1,
        },
        qc_recipe_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_recipe_data_href,
            strict_type => 1,
        },
        recipe => {
            defined     => 1,
            required    => 1,
            store       => \$recipe,
            strict_type => 1,
        },
        sample_id        => { store => \$sample_id, strict_type => 1, },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    use MIP::Qc_data qw{ add_qc_data_recipe_info set_qc_data_recipe_info };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get corresponding data metrics for header array
    my @data_metrics = get_qc_recipe_data(
        {
            attribute           => $attribute,
            recipe_name         => $recipe,
            qc_recipe_data_href => $qc_recipe_data_href,
        }
    );

    ## Enable seperation of writing array or key-->value in qc_data
    if ( scalar @data_metrics == 1 ) {

        my $data_metric = $data_metrics[0];

        set_qc_data_recipe_info(
            {
                key          => $attribute,
                infile       => $infile,
                qc_data_href => $qc_data_href,
                recipe_name  => $recipe,
                sample_id    => $sample_id,
                value        => $data_metric,
            }
        );
        set_metrics_to_store(
            {
                id           => $sample_id,
                input        => $infile,
                metric_name  => $attribute,
                metric_value => $data_metric,
                qc_data_href => $qc_data_href,
                recipe_name  => $recipe,
            }
        );
    }
    elsif ( not exists $qc_header_href->{$recipe} ) {
        ## Write array to qc_data for metrics without header

      DATA_METRIC:
        foreach my $data_metric (@data_metrics) {

            add_qc_data_recipe_info(
                {
                    key          => $attribute,
                    infile       => $infile,
                    qc_data_href => $qc_data_href,
                    recipe_name  => $recipe,
                    sample_id    => $sample_id,
                    value        => $data_metric,
                }
            );
            set_metrics_to_store(
                {
                    id           => $sample_id,
                    input        => $infile,
                    metric_name  => $attribute,
                    metric_value => $data_metric,
                    qc_data_href => $qc_data_href,
                    recipe_name  => $recipe,
                }
            );
        }
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

sub get_qc_metric {

## Function : Get metric and meta data in qc_data hash
## Returns  : @metrics
## Arguments: $header       => Metrics table header
##          : $id           => Id associated with metric (sample_id|case_id)
##          : $input        => Input source used to generate metric from
##          : $metric_name  => Name of metric
##          : $metric_value => Value to store
##          : $qc_data_href => Qc_data hash {REF}
##          : $recipe_name  => Recipe to set attributes for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $header;
    my $id;
    my $input;
    my $metric_name;
    my $metric_value;
    my $qc_data_href;
    my $recipe_name;

    my $tmpl = {
        header => {
            store       => \$header,
            strict_type => 1,
        },
        id => {
            store       => \$id,
            strict_type => 1,
        },
        input => {
            store       => \$input,
            strict_type => 1,
        },
        metric_name => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$metric_name,
        },
        metric_value => {
            store       => \$metric_value,
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
            store       => \$recipe_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @metrics;

    foreach my $metric_href ( @{ $qc_data_href->{metrics} } ) {

        if (    $metric_href->{step} eq $recipe_name
            and $metric_href->{id} eq $id
            and $metric_href->{name} eq $metric_name )
        {
            push @metrics, $metric_href;
        }
    }
    return @metrics;
}

sub get_qc_recipe_data {

## Function : Get sample recipe attributes from qc_data hash
## Returns  : "$data_metric",  "$data_metrics" or "data_metrics_hash"
## Arguments: $attribute           => Attribute key
##          : $qc_recipe_data_href => Hash to save data in each outfile {REF}
##          : $qc_header_index     => Qc header element index
##          : $recipe_name         => Recipe to examine

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $qc_recipe_data_href;
    my $qc_header_index;
    my $recipe_name;

    my $tmpl = {
        attribute => {
            store       => \$attribute,
            strict_type => 1,
        },
        qc_recipe_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_recipe_data_href,
            strict_type => 1,
        },
        qc_header_index => {
            store       => \$qc_header_index,
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
    if ( $attribute && defined $qc_header_index ) {

        return $qc_recipe_data_href->{$recipe_name}{$attribute}[$qc_header_index];
    }

    ## Get and return data metrics array
    if ($attribute) {

        return @{ $qc_recipe_data_href->{$recipe_name}{$attribute} };
    }

    ## Get qc data recipe data metric hash
    return %{ $qc_recipe_data_href->{$recipe_name} };
}

sub get_regexp_qc_data {

## Function  : Use reg exp to collect qc data via system call
## Returns   : @{ $process_return{stdouts_ref} } or @{ $process_return{stderrs_ref} }
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

    use MIP::Environment::Child_process qw{ child_process };

    ## Get return from reg exp system call
    my %process_return = child_process(
        {
            commands_ref => [ qq{$regexp $data_file_path}, ],
            process_type => q{open3},
        }
    );

    ## Print stderr if returned from regexp
    if ( @{ $process_return{stderrs_ref} } ) {

        ## Be verbose that something went wrong
        say {*STDERR} join $NEWLINE, @{ $process_return{stderrs_ref} };
        return;
    }
    return @{ $process_return{stdouts_ref} };
}

sub parse_qc_recipe_data {

## Function  : Parse qc_recipe_data and add to qc_data hash to enable write to yaml format
## Returns   :
## Arguments : $infile              => Infile to recipe
##           : $qc_data_href        => Qc data hash {REF}
##           : $qc_header_href      => Save header(s) in each outfile {REF}
##           : $qc_recipe_data_href => Hash to save data in each outfile {REF}
##           : $recipe              => Recipe to examine
##           : $regexp_href         => RegExp hash {REF}
##           : $sample_id           => SampleID
##           : $sample_info_href    => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile;
    my $qc_data_href;
    my $qc_header_href;
    my $qc_recipe_data_href;
    my $recipe;
    my $regexp_href;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        infile       => { store => \$infile, strict_type => 1, },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        qc_header_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_header_href,
            strict_type => 1,
        },
        qc_recipe_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_recipe_data_href,
            strict_type => 1,
        },
        recipe => {
            defined     => 1,
            required    => 1,
            store       => \$recipe,
            strict_type => 1,
        },
        regexp_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$regexp_href,
            strict_type => 1,
        },
        sample_id        => { store => \$sample_id, strict_type => 1, },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qc_data qw{ add_to_qc_data parse_qc_recipe_table_data set_header_metrics_to_qc_data };

  REG_EXP_ATTRIBUTE:
    for my $attribute ( keys %{ $regexp_href->{$recipe} } ) {

        ## For info contained in entry --> Value i.e. same line without header
        if ( $attribute !~ /^header|header$/ixsm ) {

            add_to_qc_data(
                {
                    attribute           => $attribute,
                    infile              => $infile,
                    qc_data_href        => $qc_data_href,
                    qc_header_href      => $qc_header_href,
                    qc_recipe_data_href => $qc_recipe_data_href,
                    recipe              => $recipe,
                    sample_id           => $sample_id,
                    sample_info_href    => $sample_info_href,
                }
            );
        }
        else {
            ## Table data i.e. header and subsequent data line(s).
            ## Can be multiple tables per file

            parse_qc_recipe_table_data(
                {
                    infile              => $infile,
                    qc_data_href        => $qc_data_href,
                    qc_header_href      => $qc_header_href,
                    qc_recipe_data_href => $qc_recipe_data_href,
                    regexp_href         => $regexp_href,
                    recipe              => $recipe,
                    sample_id           => $sample_id,
                }
            );
        }
    }
    return 1;
}

sub parse_qc_recipe_table_data {

## Function  : Parse qc_recipe_data in table form and set metrics to qc_data hash to enable write to yaml format
## Returns   :
## Arguments : $infile              => Infile to recipe
##           : $qc_data_href        => Qc data hash {REF}
##           : $qc_header_href      => Save header(s) in each outfile {REF}
##           : $qc_recipe_data_href => Hash to save data in each outfile {REF}
##           : $recipe              => Recipe to examine
##           : $regexp_href         => RegExp hash {REF}
##           : $sample_id           => SampleID

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile;
    my $qc_data_href;
    my $qc_header_href;
    my $qc_recipe_data_href;
    my $recipe;
    my $regexp_href;
    my $sample_id;

    my $tmpl = {
        infile       => { store => \$infile, strict_type => 1, },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        qc_header_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_header_href,
            strict_type => 1,
        },
        qc_recipe_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_recipe_data_href,
            strict_type => 1,
        },
        recipe => {
            defined     => 1,
            required    => 1,
            store       => \$recipe,
            strict_type => 1,
        },
        regexp_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$regexp_href,
            strict_type => 1,
        },
        sample_id => { store => \$sample_id, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qc_data qw{ get_qc_recipe_data };

  TABLE_HEADER_KEY:
    for my $regexp_header_key ( keys %{ $qc_header_href->{$recipe} } ) {

      PARAGRAPH_KEY:
        for my $regexp_key ( keys %{ $regexp_href->{$recipe} } ) {

            ## Detect if the regexp id is for data and not header
            next PARAGRAPH_KEY if ( $regexp_key =~ /^header|header$/isxm );

            ## For all collected headers for this paragraph
          HEADER_VALUE:
            while ( my ( $qc_header_index, $qc_header ) =
                each @{ $qc_header_href->{$recipe}{$regexp_header_key} } )
            {
                ## Get corresponding data metric for header element index
                my $data_metric = get_qc_recipe_data(
                    {
                        attribute           => $regexp_key,
                        recipe_name         => $recipe,
                        qc_recipe_data_href => $qc_recipe_data_href,
                        qc_header_index     => $qc_header_index,
                    }
                );

                ## Set table metric data to qc_data hash
                set_header_metrics_to_qc_data(
                    {
                        infile            => $infile,
                        key               => $qc_header,
                        qc_data_href      => $qc_data_href,
                        regexp_header_key => $regexp_header_key,
                        regexp_key        => $regexp_key,
                        recipe_name       => $recipe,
                        sample_id         => $sample_id,
                        value             => $data_metric,
                    }
                );
                set_metrics_to_store(
                    {
                        header       => $regexp_key,
                        id           => $sample_id,
                        input        => $infile,
                        metric_name  => $qc_header,
                        metric_value => $data_metric,
                        qc_data_href => $qc_data_href,
                        recipe_name  => $recipe,
                    }
                );
            }
        }
    }
    return 1;
}

sub parse_regexp_hash_and_collect {

## Function  : Parses the regexp hash structure to identify if the info is
##             1) Table section(s) (both header(s) and data line(s) 2) Seperate data line.
## Returns   :
## Arguments : $outdirectory        => Recipes outdirectory
##           : $outfile             => Recipes outfile containing parameter to evaluate
##           : $qc_header_href      => Save header(s) in each outfile {REF}
##           : $qc_recipe_data_href => Hash to save data in each outfile {REF}
##           : $recipe              => The recipe to examine
##           : $regexp_href         => Regexp hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $outdirectory;
    my $outfile;
    my $qc_header_href;
    my $qc_recipe_data_href;
    my $recipe;
    my $regexp_href;

    my $tmpl = {
        outdirectory   => { store => \$outdirectory, strict_type => 1, },
        outfile        => { store => \$outfile,      strict_type => 1, },
        qc_header_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_header_href,
            strict_type => 1,
        },
        qc_recipe_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_recipe_data_href,
            strict_type => 1,
        },
        recipe => {
            defined     => 1,
            required    => 1,
            store       => \$recipe,
            strict_type => 1,
        },
        regexp_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$regexp_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qcc_regexp qw{ get_qcc_regexp_recipe_attribute };
    use MIP::Qc_data qw{ add_qc_data_regexp_return };

    <<"FUNCTION";
        ## Detect if the outfile contains table info in the outfile
        ## i.e. data is formated as a paragraf with header(s) and line(s).
        ## "regexp_key" should either start with or end with "header". This
        ## section extracts the header/data line(s) for the entire outdata file.
        ## Necessary to assign correct data entry to header entry later
        ## (headers and data are saved in seperate hashes).
FUNCTION

    ## Find the regular expression(s) for each recipe that is used
  REG_EXP_KEY:
    for my $regexp_key ( keys %{ $regexp_href->{$recipe} } ) {

        ## Regular expression used to collect paragraf header info
        my $regexp = get_qcc_regexp_recipe_attribute(
            {
                attribute       => $regexp_key,
                qcc_regexp_href => $regexp_href,
                recipe_name     => $recipe,
            }
        );

        ## Detect if the regexp key is a paragraf header line
        if ( $regexp_key =~ /^header|header$/isxm ) {

            ## Add qc data from data outfile using regexp
            add_qc_data_regexp_return(
                {
                    data_file_path => catfile( $outdirectory, $outfile ),
                    qc_href        => $qc_header_href,
                    recipe_name    => $recipe,
                    regexp         => $regexp,
                    regexp_key     => $regexp_key,
                }
            );
            next REG_EXP_KEY;
        }

        ### For metrics contained in data line.
        add_qc_data_regexp_return(
            {
                data_file_path => catfile( $outdirectory, $outfile ),
                qc_href        => $qc_recipe_data_href,
                recipe_name    => $recipe,
                regexp         => $regexp,
                regexp_key     => $regexp_key,
            }
        );
    }
    return 1;
}

sub set_header_metrics_to_qc_data {

## Function : Set table metric data to qc_data hash
## Returns  :
## Arguments: $infile            => Infile key
##          : $key               => Metafile tag
##          : $qc_data_href      => Qc_data hash {REF}
##          : $regexp_header_key => Regexp header key
##          : $regexp_key        => Regexp key
##          : $recipe_name       => Recipe to set attributes for
##          : $sample_id         => Sample ID
##          : $value             => Value to store

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile;
    my $key;
    my $qc_data_href;
    my $regexp_header_key;
    my $regexp_key;
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
        regexp_header_key => {
            defined     => 1,
            required    => 1,
            store       => \$regexp_header_key,
            strict_type => 1,
        },
        regexp_key => {
            defined     => 1,
            required    => 1,
            store       => \$regexp_key,
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

    ## Set recipe key value pair for paragraph header data on sample, infile
    ## and recipe level
    if ( $sample_id and $infile ) {

        $qc_data_href->{sample}{$sample_id}{$infile}{$recipe_name}{$regexp_header_key}
          {$regexp_key}{$key} = $value;
        return;
    }

    $qc_data_href->{$recipe_name}{$regexp_header_key}{$regexp_key}{$key} = $value;

    return;
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

sub set_metrics_to_store {

## Function : Set metric and meta data in qc_data hash
## Returns  :
## Arguments: $header       => Metrics table header
##          : $id           => Id associated with metric (sample_id|case_id)
##          : $input        => Input source used to generate metric from
##          : $metric_name  => Name of metric
##          : $metric_value => Value to store
##          : $qc_data_href => Qc_data hash {REF}
##          : $recipe_name  => Recipe to set attributes for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $header;
    my $id;
    my $input;
    my $metric_name;
    my $metric_value;
    my $qc_data_href;
    my $recipe_name;

    my $tmpl = {
        header => {
            store       => \$header,
            strict_type => 1,
        },
        id => {
            store       => \$id,
            strict_type => 1,
        },
        input => {
            store       => \$input,
            strict_type => 1,
        },
        metric_name => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$metric_name,
        },
        metric_value => {
            store       => \$metric_value,
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
            store       => \$recipe_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not defined $metric_value );

    ## Build metric meta data hash
    my %metric_info = (
        header => $header,
        id     => $id,
        input  => $input,
        name   => $metric_name,
        step   => $recipe_name,
        value  => $metric_value,
    );

    ## Set metric according to metric info
    push @{ $qc_data_href->{metrics} }, {%metric_info};

    return;
}
1;
