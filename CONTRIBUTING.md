# Contributing

Contributing to Finch is easy: just open a [pull
request](https://help.github.com/articles/using-pull-requests/). Make
`main` the destination branch on the [Finch
repository](https://code.ornl.gov/8s2/finch) and allow edits from
maintainers in the pull request.

Your pull request must pass use the coding style from `.clang-format`
(enforced with clang-format-14) and adding doxygen
documentation, and be reviewed by at least one Finch developer.

`pre-commit` is a useful tool for ensuring feature branches are ready for
review by running automatic checks locally before a commit is made.
[Installation details](https://pre-commit.com/#install) (once per system) and
[activation details](https://pre-commit.com/#usage) (once per repo) are
available.

Other coding style includes:
* Camel case template parameters (`NewTemplateType`)
* Camel case class names (`NewClassName`)
* Lower camel case function names (`newFunctionName`)
* Lower case, underscore separated variables (`new_variable_name`)
* Class members which are `private` are followed by an underscore (`private_class_variable_`)
* Class/struct member type aliases use lower case, underscore separated names (`using integer_type = int;`)

Use the following naming convections for folders, files, libraries, and executables:
* Camel case folders and files (`NewFileName`)
    * the execption being a main source file (`main.cpp`)
* Camel case library names, same as file name (`NewLibraryName`)
* Lower case, underscore executable names (`new_executable_name`)
