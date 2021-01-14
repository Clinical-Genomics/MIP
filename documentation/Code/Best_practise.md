# Perl Best Practices

Following best practices is a good way to create maintainable and readable code and should always be encouraged. However, learning what these best practices are and when they apply in the context of your code can be hard to determine. Luckily, there are several tools to help guide you on your way.

## Code standards
All code that are present with MIP has to be processed by both [Perltidy], [Perlcritic] (Perl code) and [yamllint] (yaml code). MIP uses [pre-commit] to automatically run Perlcritic, Perltidy and yamllint on each commit. Perltidy modifies the files in place and the files are then analyzed by Perl critic with a level 1 severity with a few exceptions as specified in the .perlcriticrc file. These programs can also be used as standalone tools or be integrated into your IDE.

### Perl Critic

[Perlcritic] is a Perl source code analyzer. It is the executable front-end to the Perl::Critic engine, which attempts to identify awkward, hard to read, error-prone, or unconventional constructs in your code. Most of the rules are based on Damian Conway's book Perl Best Practices.

Perl critic allows different degrees of severity. Severity values are integers ranging from 1 (least severe) to 5 (most severe) when analyzing your code. You can also use the severity names if you think it hard to remember the meaning of the integers (1 = brutal and 5 = gentle). The level is controlled with the '--severity' flag. There is also a verbose flag, to print more information about the identified deviations from the perl critic best practises. There are 11 levels of verbosity.

#### Examples

```
perlcritic --severity 4 --verbose 11 my_perl_script.pl
```

Perl critic also has [web interface] to instantly analyze your code.


### Perl Tidy

[Perltidy] is a Perl script which indents and reformats Perl scripts to make them easier to read. If you write Perl scripts, or spend much time reading them, you will probably find it useful. Perltidy is an excellent way to automate the code standardisation with minimum of effort.  

#### Examples

```
perltidy somefile.pl
```

This will produce a file somefile.pl.tdy containing the script reformatted using the default options, which approximate the style suggested in perlstyle(1). The source file somefile.pl is unchanged.

```
perltidy -b -bext='/' file1.pl file2.pl
```

Create backups of files and modify files in place. The backup files file1.pl.bak and file2.pl.bak will be deleted if there are no errors.

### Yamllint
[yamllint] is a linter for YAML files. Yamllint does not only check for syntax validity, but for weirdnesses like key repetition and cosmetic problems such as lines length, trailing spaces, indentation, etc.

[pre-commit]: https://pre-commit.com
[Perlcritic]: http://search.cpan.org/~petdance/Perl-Critic/bin/perlcritic
[Perltidy]: http://perltidy.sourceforge.net/
[web interface]: http://perlcritic.com/
[yamllint]: https://github.com/adrienverge/yamllint
