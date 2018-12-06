# IO

**Version: 1.0.0**

MIP stores the information on file object IO in the `io` sub-level hash located in `file_info`.
The io hash has three hash reference keys: `in`, `out`, and `temp` with identical sub-level keys. Setting either `in` or `out` also sets `temp`.

## File objects

1. "a/dir/file_1.1.txt"
2. "a/dir/file_1.txt"

```
io:
{
  in {file_paths: [  "a/dir/file_1.1.txt"
                      "a/dir/file_1.txt"
                      ],
      file_suffixes: [ ".1.txt"
                        ".txt"
                        ],
      file_path_prefixes: [ "a/dir/file_1"
                            "a/dir/file_1"
                            ],
      file_names: [ "file_1.1.txt"
                    "file_1.txt"
                    ],
      file_name_prefixes: [ "file_1"
                            "file_1"
                            ],
      dir_path: "a/dir/",
      dir_path_prefix: "a/dir",
      file_suffix: ".txt",
      file_path_prefix: "a/dir/file_1",
      file_name_prefix: "file_1"
      file_constant_suffix: undef,
      file_path_href: { 1 => "a/dir/file_1.1.txt",
			},
      file_name_href: { 1 => "file_1.1.txt",
                        },
	}
  out {"exactly the same keys and types"}                  
  temp {"exactly the same keys and types, but with a temp directory base instead of 'in' or 'out' base"}                  
}
```
