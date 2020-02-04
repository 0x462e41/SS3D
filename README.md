# SS3D
> Sequence Similarity 3D

SS3D is a bioinformatics tool for combining sequence and structure alignment data,
evaluating the Sequence Similarity in 3D space between two proteins. Such comparisons
are useful for a plethora of purposes such as determination of conserved positions,
and explaining different functionalities of proteins that evolved from a common ancestor.

For more information about this tool, check the paper at <upcoming link>.

## Installation

OS X & Linux:

The only dependence is the sdt14 library itself, and the gcc compiler. You can also use
others compilers, just change the 'CC' variable present at the top of the ``makefile``.

After downloading/cloning the repository, just enter the folder and type the command:

```sh
make
```

It will create an executable called 'ss3d' inside this folder.

Windows:

It should be easier to import the repository as a new project and compile it using Code::Blocks,
Visual Studio or any c++14 IDE/Compiler.

## Usage example

To get a quick overview of the required and optional parameters, type:

```sh
./ss3d -h
```

or

```sh
./ss3d help
```

_For examples and a tutorial, please refer to the [Wiki][wiki]._

## Meta

Igor Daniel M. Lima – igor.sj13@gmail.com

Distributed under the GNU GPL-3.0 license. See ``LICENSE`` for more information.

[https://github.com/IgorLima/](https://github.com/0x462e41/)

## Contributing

1. Fork it (<https://github.com/0x462e41/SS3D/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

<!-- Markdown link & img dfn's -->
[wiki]: https://github.com/0x462e41/SS3D/wiki/Tutorial
