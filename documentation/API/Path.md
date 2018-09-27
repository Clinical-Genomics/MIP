# File objects

_"On heights all paths are paved with dagger -Old Seanchan saying"_

1. "a/dir/file_1.1.txt"
2. "a/dir/file_1.txt"

## Definitions:

### Array file features

#### file_paths:
Complete paths
- [0] "a/dir/file_1.1.txt"
- [1] "a/dir/file_1.txt"

#### file_suffixes:
Including first "." and everything after
- [0] ".1.txt"
- [1] ".txt"

#### file_path_prefixes:
Paths without suffix
- [0] "a/dir/file_1"
- [1] "a/dir/file_1"

#### file_names:
File names
- [0] "file_1.1.txt"
- [1] "file_1.txt"

#### file_name_prefixes:
File names without suffix
- [0] "file_1"
- [1] "file_1"

### Scalar file features

#### dir_path:
Directory path
- "a/dir/"

#### dir_path_prefix:
Directory name without trailing slash
- "a/dir"

##### file_suffix:
Everything after last "."
- ".txt"

### Constant file features

#### file_path_prefixes:
Paths without suffix
- "a/dir/file_1"

#### file_name_prefix:
File name without suffix
- "file_1"

##### file_constant_suffix:
Only for identical file_suffixes otherwise undef

For paths:
- [0] "a/dir/file.1.txt"
- [1] "a/dir/file.1.txt"

*file_constant_suffix*:
- ".1.txt"

### Hash file features

iterator: string between dots in file_suffixes, e.g. ".1.txt" => iterator = "1".

#### file_path_href:
Complete paths per iterator.
{ 1 => "a/dir/file_1.1.txt",
}

#### file_name_href:
File name per iterator.
{ 1 => "file_1.1.txt",
}

