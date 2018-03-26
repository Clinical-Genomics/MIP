# Testing Data Structures
When testing complex data structures use Test::More and the is_deeply sub.

```Perl
my @covariets = qw{ ReadGroupCovariate ContextCovariate CycleCovariate QualityScoreCovariate };
my @expected_covariets = qw{ ReadGroupCovariate ContextCovariate CycleCovariate QualityScoreCovariate };

is_deeply( \@covariets, \@expected_covariets, q{Identical arrays} );

my %colors = (red  => q{car},
	      blue => q{sky},
             );
my %expected_colors = (red  => q{car},
                       blue => q{sky},
                      );
is_deeply(\%colors, \%expected_colors, q{Identical hashes});
```
