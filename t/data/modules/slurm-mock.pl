#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

my ($no_job_id);

###User Options
GetOptions( q{nj|no_job_id} => \$no_job_id, );

## Print without job id
if ($no_job_id) {

    print {*STDOUT} q{Submitted batch job};
    exit;
}

## Correct return
print {*STDOUT} q{Submitted batch job 1234};
