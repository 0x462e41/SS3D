# SS3D
> Sequence Similarity 3D

SS3D is a bioinformatics tool for combining sequence and structure alignment data,
evaluating the Sequence Similarity in 3D space between two proteins. Such comparisons
are useful for a plethora of purposes such as determination of conserved positions,
and explaining different functionalities of proteins that evolved from a common ancestor.

For more information about this tool, check the paper at ~~upcoming link~~.

## Installation

### OS X & Linux:

The only dependences are the sdt14 library itself, and the gcc compiler.
You can install it using the following command:

```sh
sudo apt install build-essential
```

You can also use other compilers, just change the 'CC' variable present at the top
of the ``makefile``.
After downloading/cloning the repository, just enter the folder and type the command:

```sh
make
```

It will create an executable called 'ss3d' inside the folder.

### Windows:

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

>   SYNOPSIS
>
>   ss3d [-a [.pdb]] [-b [.pdb]] [-o [.xvg]] [-dist <number>] [-skip <number>] [-rad <number>] [-min <number>] [-mat <matrix type>] [-dup]
>
>   OPTIONS:
>
>   Obrigatory options to specify input & output files:
>
>    -a
>           Structural file in .pdb format of the first protein.
>    -b
>           Structural file in .pdb format of the second protein.
>    -o
>           Output text file in .xvg format.
>
>   Others options:
>
>    -dist
>           The maximum distance to be considered a contact between alpha carbons, in Ångströms.
>           (default dinstance is 10 Ångströms)
>    -skip
>           The minimum distance between residues in the primary structure to be considered,
>           in absolute number (default skip is 3).
>    -rad
>           The maximum radius to search for interaction between alpha carbons, in Ångströms.
>           (default radius is 10 Ångströms)
>    -min
>           Minimum number of common residues around the contact to be considered.
>           (default is 1).
>    -mat
>           Select the substitution matrix to be used to evaluate the conscore. The default matrix is blosum62.
>           The options are:
>           pam250 - The PAM250 matrix.
>           blosum62 - The BLOSUM62 matrix.
>    -dup
>           Duplicates the matrix in a mirrored fashion.


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
