package MIP::Test::Fixtures;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use File::Temp;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $NEWLINE $SPACE };
use MIP::Script::Utils qw{ help };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ test_add_io_for_recipe test_constants test_import test_log test_mip_hashes };
}

sub test_add_io_for_recipe {

## Function : Add io from upstream recipe
## Returns  :
## Arguments: $chain_id          => Chain id of recipe
##          : $file_info_href    => File info hash {REF}
##          : $id                => Sample or case
##          : $parameter_href    => Parameter hash {REF}
##          : $order_recipes_ref => Order of recipes {REF}
##          : $outfile_suffix    => Set outfile suffix in parameters
##          : $recipe_name       => Recipe name
##          : $step              => Level to inherit from

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $parameter_href;
    my $order_recipes_ref;
    my $outfile_suffix;
    my $recipe_name;

    ## Default(s)
    my $id;
    my $chain_id;
    my $step;

    my $tmpl = {
        chain_id => {
            default     => q{TEST},
            store       => \$chain_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        id => {
            default     => q{ADM1059A1},
            store       => \$id,
            strict_type => 1,
        },
        order_recipes_ref => {
            default     => [],
            store       => \$order_recipes_ref,
            strict_type => 1,
        },
        outfile_suffix => {
            store       => \$outfile_suffix,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        step => {
            allow       => [qw{ fastq bam vcf vcf.gz hdf5 tsv }],
            default     => q{fastq},
            store       => \$step,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @order_recipes =
        @{$order_recipes_ref}
      ? @{$order_recipes_ref}
      : ( qw{ first_recipe }, $recipe_name );

    ## Add a chain to the recipe
    $parameter_href->{$recipe_name}{chain} = $chain_id;

    if ( $step eq q{fastq} ) {

        @{ $parameter_href->{cache}{order_recipes_ref} } = $recipe_name;
        $parameter_href->{$recipe_name}{outfile_suffix} =
          $outfile_suffix ? $outfile_suffix : q{.bam};
    }
    if ( $step eq q{bam} ) {

        %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} } = test_mip_hashes(
            {
                mip_hash_name => q{io_bam},
            }
        );
        @{ $parameter_href->{cache}{order_recipes_ref} } = @order_recipes;
        $parameter_href->{$recipe_name}{outfile_suffix} =
          $outfile_suffix ? $outfile_suffix : q{.bam};
    }
    if ( $step eq q{vcf} ) {

        %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} } = test_mip_hashes(
            {
                mip_hash_name => q{io_vcf},
            }
        );
        @{ $parameter_href->{cache}{order_recipes_ref} } = @order_recipes;
        $parameter_href->{$recipe_name}{outfile_suffix} =
          $outfile_suffix ? $outfile_suffix : q{.vcf};
    }
    if ( $step eq q{vcf.gz} ) {

        %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} } = test_mip_hashes(
            {
                mip_hash_name => q{io_vcf_gz},
            }
        );
        @{ $parameter_href->{cache}{order_recipes_ref} } = @order_recipes;
        $parameter_href->{$recipe_name}{outfile_suffix} =
          $outfile_suffix ? $outfile_suffix : q{.vcf.gz};
    }
    if ( $step eq q{hdf5} ) {

        %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} } = test_mip_hashes(
            {
                mip_hash_name => q{io_hdf5},
            }
        );
        @{ $parameter_href->{cache}{order_recipes_ref} } = @order_recipes;
        $parameter_href->{$recipe_name}{outfile_suffix} =
          $outfile_suffix ? $outfile_suffix : q{.hdf5};
    }
    if ( $step eq q{tsv} ) {

        %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} } = test_mip_hashes(
            {
                mip_hash_name => q{io_tsv},
            }
        );
        @{ $parameter_href->{cache}{order_recipes_ref} } = @order_recipes;
        $parameter_href->{$recipe_name}{outfile_suffix} =
          $outfile_suffix ? $outfile_suffix : q{.tsv};
    }
    return;
}

sub test_import {

## Function : Test modules and imports
## Returns  :
## Arguments: $perl_module_href => Modules with imports to test

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $perl_module_href;

    my $tmpl = {
        perl_module_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$perl_module_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Test::More;

    ## Modules with import
  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %{$perl_module_href} ) {

        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;

        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

    return;
}

sub test_log {

## Function : Generate a log object and a temporary log file
## Returns  : $log
## Arguments: $log_level => Log level
##          : $log_name  => Name of logger
##          : $no_screen => Don't log to screen

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log_name;
    my $no_screen;

    ## Default(s)
    my $log_level;

    my $tmpl = {
        log_level => {
            allow       => [qw{ DEBUG ERROR FATAL INFO TRACE WARN }],
            default     => q{TRACE},
            store       => \$log_level,
            strict_type => 1,
        },
        log_name => {
            default     => q{TEST},
            store       => \$log_name,
            strict_type => 1,
        },
        no_screen => {
            allow       => [ undef, 0, 1 ],
            store       => \$no_screen,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Log::MIP_log4perl qw{ initiate_logger };

    ## Create temp logger
    my $test_dir      = File::Temp->newdir();
    my $test_log_path = catfile( $test_dir, q{test.log} );

    ## Set default log categories
    my @categories = ( $log_level, qw{ LogFile ScreenApp } );

    ## Disable print to screen
    if ($no_screen) {
        @categories = ( $log_level, q{LogFile} );
    }

    ## Creates log object
    my $log = initiate_logger(
        {
            categories_ref => \@categories,
            file_path      => $test_log_path,
            log_name       => $log_name,
        }
    );

    return $log;
}

sub test_mip_hashes {

## Function : Loads test MIP hashes with core parameters set e.g. active_parameter
## Returns  : MIP core hash
## Arguments: $mip_hash_name  => MIP core hash to return
##          : $recipe_name    => Recipe name
##          : $temp_directory => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $mip_hash_name;
    my $recipe_name;
    my $temp_directory;

    my $tmpl = {
        mip_hash_name => {
            allow => [
                qw{ active_parameter
                  define_parameter
                  dependency_tree_dna
                  dependency_tree_rna
                  download_active_parameter
                  file_info
                  install_active_parameter
                  io
                  io_bam
                  io_vcf
                  io_vcf_gz
                  io_hdf5
                  io_tsv
                  job_id
                  pedigree
                  primary_contig
                  recipe_parameter
                  qc_sample_info }
            ],
            defined     => 1,
            required    => 1,
            store       => \$mip_hash_name,
            strict_type => 1,
        },
        recipe_name => {
            default     => q{bwa_mem},
            store       => \$recipe_name,
            strict_type => 1,
        },
        temp_directory => {
            default     => catfile( File::Temp->newdir() ),
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Read qw{ read_from_file };

    my %test_hash = (
        active_parameter    => catfile( $Bin, qw{ data test_data recipe_active_parameter.yaml } ),
        define_parameter    => catfile( $Bin, qw{ data test_data define_parameters.yaml } ),
        dependency_tree_dna => catfile( $Bin, qw{ data test_data rd_dna_initiation_map.yaml } ),
        dependency_tree_rna => catfile( $Bin, qw{ data test_data rd_rna_initiation_map.yaml } ),
        download_active_parameter =>
          catfile( $Bin, qw{ data test_data download_active_parameters.yaml } ),
        file_info                => catfile( $Bin, qw{ data test_data recipe_file_info.yaml } ),
        install_active_parameter =>
          catfile( $Bin, qw{ data test_data install_active_parameters.yaml } ),
        io               => catfile( $Bin, qw{ data test_data io.yaml } ),
        io_bam           => catfile( $Bin, qw{ data test_data io_bam.yaml } ),
        io_vcf           => catfile( $Bin, qw{ data test_data io_vcf.yaml } ),
        io_vcf_gz        => catfile( $Bin, qw{ data test_data io_vcf_gz.yaml } ),
        io_hdf5          => catfile( $Bin, qw{ data test_data io_hdf5.yaml } ),
        io_tsv           => catfile( $Bin, qw{ data test_data io_tsv.yaml } ),
        job_id           => catfile( $Bin, qw{ data test_data job_id.yaml } ),
        pedigree         => catfile( $Bin, qw{ data test_data pedigree_wes.yaml } ),
        primary_contig   => catfile( $Bin, qw{ data test_data primary_contig.yaml } ),
        recipe_parameter => catfile( $Bin, qw{ data test_data recipe_parameter.yaml } ),
        qc_sample_info => catfile( $Bin, qw{ data test_data 643594-miptest_qc_sample_info.yaml } ),
    );

    my %hash_to_return = read_from_file(
        {
            format => q{yaml},
            path   => $test_hash{$mip_hash_name},
        }
    );

    ## Add dynamic parameters
    if ( $mip_hash_name eq q{active_parameter} ) {

        ## Adds the recipe name
        $hash_to_return{$recipe_name} = 2;

        ## Adds reference dir
        $hash_to_return{reference_dir} = catfile( $Bin, qw{ data test_data references } );

        ## Adds parameters with temp directory
        $hash_to_return{outdata_dir}    = catfile( $temp_directory, q{test_data_dir} );
        $hash_to_return{outscript_dir}  = catfile( $temp_directory, q{test_script_dir} );
        $hash_to_return{temp_directory} = $temp_directory;
    }
    if ( $mip_hash_name eq q{recipe_parameter} ) {

        ## Adds a recipe chain
        $hash_to_return{$recipe_name}{chain} = q{TEST};
    }
    return %hash_to_return;
}

sub test_constants {

## Function : Generate standard command line interface for test scripts
## Returns  : $verbose
## Arguments: $container_manager        => Set container manager
##          : $test_mode                => Version of test
##          : $test_process_return_href => Process return hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $container_manager;
    my $test_mode;
    my $test_process_return_href;

    my $tmpl = {
        container_manager => {
            allow       => [qw{docker singularity}],
            default     => q{docker},
            store       => \$container_manager,
            strict_type => 1,
        },
        test_mode => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$test_mode,
            strict_type => 1,
        },
        test_process_return_href => {
            default     => {},
            store       => \$test_process_return_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Constants qw{ set_test_constants };

    if ( not keys %{$test_process_return_href} ) {

        %{$test_process_return_href} = (
            buffers_ref   => [],
            error_message => undef,
            stderrs_ref   => [],
            stdouts_ref   => [qw{ one }],
            success       => 1,
        );
    }

    set_test_constants(
        {
            container_manager        => $container_manager,
            test_mode                => $test_mode,
            test_process_return_href => $test_process_return_href,
        }
    );

    return;
}
1;
