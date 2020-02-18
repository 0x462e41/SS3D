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

## Usage

To get a quick overview of the required and optional parameters, type:

```sh
./ss3d help
or
./ss3d -h
```

>   SYNOPSIS
>
>   ss3d [-a [.pdb]] [-b [.pdb]] [-o [.xvg/.pdb]] [-dist <number>] [-skip <number>] [-rad <number>]
>   [-min <number>] [-mat <matrix type>] [-dup] [-norm] [-bfac] [-raw]
>
>   OPTIONS:
>
>   Obligatory options to specify input & output files:
>
>    -a
>
>           Structural file in .pdb format of the first protein.
>
>    -b
>
>           Structural file in .pdb format of the second protein.
>
>    -o
>
>           Output text file in .xvg or .pdb format.
>
>   Other options:
>
>    -dist
>
>           The maximum distance to be considered a contact between two alpha carbons, in Ångströms.
>           (default distance is 10 Ångströms)
>
>    -rad
>
>           The maximum search radius for identifying alpha carbons around a contact, in Ångströms.
>           (default radius is 10 Ångströms)
>
>    -skip
>
>           Exclude this number of residues within the search radius upstream and downstream in primary
>           sequence of the alpha carbon involved in the contact. Used to diminish the contribution of
>           sequentially proximal residues in the scoring. (default skip is 3, which will exclude
>           6 residues total: the three previous and three subsequent to the one involved in the contact).
>
>    -min
>
>           Minimum (absolute) number of common residues between the two proteins around a particular contact
>           for it to be accepted. (default is 1 residue[s]).
>
>    -mat
>
>           Select the substitution matrix to be used to evaluate the score. The default matrix is blosum62.
>           The options are:
>           pam250 - The PAM250 matrix.
>           blosum62 - The BLOSUM62 matrix.
>
>    -dup
>
>           Duplicates the matrix in a mirrored fashion.
>
>    -norm
>
>           Normalizes the result to be between 0 and 1.
>
>    -bfac
>
>           Maps the score value to the B-factor of the first protein PDB.
>
>    -raw
>
>           Produces a raw text file containing the local alignment for each contact, ignoring the
>           parameters ``-bfac``, ``-dup`` and ``-norm``, if present.


_For examples and a step-by-step tutorial, please refer to the [Wiki][wiki]._

## Quick test the SS3D

To quickly test the SS3D tool, we will assume that you have [gnuplot](http://www.gnuplot.info/)
and [vmd](http://www.ks.uiuc.edu/Research/vmd/) previously installed on your computer.

Access the "tutorial" folder and type the following into the terminal:

### Generate matrix

```sh
../ss3d -a 1a5r_processed_renumbered.pdb -b 1ubq_processed_renumbered.pdb -o text-output -norm -dup
gnuplot gen_matrix.gnu
```

### Visualize pdb with SS3D coloring:

```sh
../ss3d -a 1a5r_processed_renumbered.pdb -b 1ubq_processed_renumbered.pdb -o pdb-output -norm -bfac
vmd -e view_bfactor.vmd
```

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
