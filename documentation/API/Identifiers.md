# Identifiers

The single most important practice when creating names is to devise a set of grammar rules to which all names must conform. A grammar rules specifies one or more templates (e.g., Noun :: Adjective :: Adjective) that describe how to form the entity on the left of the arrow (e.g., namespace).

A suitable grammar rule for naming packages and classes is:
```
namespace ➝  Noun :: Adjective :: Adjective 
           | Noun :: Adjective
           | Noun
```

This rule might produce package names such as:

     package Disk;
     package Disk::Audio;
     package Disk::DVD;
     package Disk::DVD::Rewritable;

In this scheme, specialized versions of an existing namespace are named by adding adjectives to the name of the more general namespace.

Variables should be named according to the data they will store, and as specifically as possible. Variables that are used in more than one block should always have a two-part (or longer) name.

A variable is named with a noun, preceded by zero or more adjectives:

```
variable ➝ [adjective _ ]* noun
```

The choice of nouns and adjectives is critical. The nouns in particular should indicate what the variable does in terms of the problem domain, not in terms of the implementation. Use adjectives whenever you can.

There is one extra grammatical variation that applies only to hashes and arrays that are used as look-up tables:

```
lookup_variable ➝ [adjective _ ]* noun preposition
```

Adding a preposition to the end of the name makes hash and array accesses much more readable.

```
my %title_of;
my @sales_from;
```

For subroutines and methods, a suitable grammatical rule for forming names is:

```
routine ➝  imperative_verb [ _ adjective]? _ noun _ preposition
         | imperative_verb [ _ adjective]? _ noun _ participle
         | imperative_verb [ _ adjective]? _ noun
```

This rule results in subroutine names such as:

```
sub get_record;                         # imperative_verb noun
sub get_record_for;                     # imperative_verb noun preposition
sub eat_cookie;                         # imperative_verb noun
sub eat_previous_cookie;                # imperative_verb adjective noun
sub build_profile;                      # imperative_verb noun
sub build_execution_profile;            # imperative_verb adjective noun
sub build_execution_profile_using;      # imperative_verb adjective noun participle
```
These naming rules particularly the two that put participles or prepositions at the ends of names create identifiers that read far more naturally, often eliminating the need for any additional comments.

Recommended imperative verbs usage:

*add_*: Is use to add to an already existing variable, data container or part of a data container

*get_*: Is used to fetch a scalar, array, hash or object from a data container  (e.g. array or hash)

*parse_*: Is used to decode data containers using iteration and/or conditions and usually calls other subroutines internally.

*set_*: Is used to replace a variable, data container (e.g. array or hash) or part of a data container
