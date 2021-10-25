# pdbman

Analyzes and modifies PDB files for QM/MM calculations.

## Description

This program takes a PDB file as input and analyzes or changes it as needed to
serve as input for Orca QM/MM calculations. 

Note that for the time being two version of this software is provided: a Python and a Rust version. They are roughly equal in terms of features with the Rust version being faster and more actively developed. The Python version will be deprecated at some point.

### Capabilities

- Query specific atoms and residues by ID or name
- Find atoms or residues within a user-defined sphere of a given atom
- Analyze atoms and/or residues present in QM or active region
- Find clashes and atomic contacts in given PDB file
- Add or remove atoms or residues to QM or active region by ID or name
- Add or remove atoms and residues to QM or active region by calculating a sphere of given radius around a given atom
- Add or remove only the sidechain or backbone of residues by ID or name

The several options are provided via command line flags (the ordering of the flags does not matter). Additional information can be obtained from `help`

---

### Options

Four different subcommands are available (equivalent aliases given in parentheses):

- `Analyze` (`analyze`, `ana`, `Y`, `y`), analyzes the QM and active region definitions currently in place
- `Remove` (`remove`, `rem`, `R`, `r`), removes atoms or residues from QM or active region
- `Add` (`add`, `A`, `a`), adds atoms or residues from QM or active regions
- `Query` (`query`, `que`, `Q`, `q`), queries the PDB file for information on atoms, residues, atom names etc.

Each subcommand can be called by any of its aliases like so:

`analyze`

`rem`

`A -rql 1,5,18`

&nbsp;

Each subcommand has a set of options they accept which are in turn divided into groups:

#### Region

Defines which region pdbman acts on. Relevant for `Add`, `Remove` and `Analyze` subcommands.

- QM1 (`-q`)
- QM2 (`-o`), this is only relevant for ONIOM calculations.
- Active (`-a`)

#### Target

Defines whether atoms or residues are targeted. Relevant for all subcommands.

- atoms (`-t`)
- residues (`-r`)

#### Source

The source of user-input. Relevant for `Add`, `Remove` and `Query` subcommands.

- list (`-l`), takes a comma-separated list of IDs **or** names, ranges are also supported.
- sphere (`-s`), takes an atom ID and a radius in Angström as arguments.

#### Partials

For cases when only a part of one or more residues should be manipulated. Relevant for `Add` and `Remove` subcommands. Requires the `residues` option.

- backbone (`-b`)
- sidechain (`-d`)

Help messages are available for all subcommands like so:

`pdbman Add --help`

## Installation

### Binary

Linux and Windows executables are provided (`pdbman` and `pdbman.exe`, respectively). Just put them somewhere in your $PATH and you're good to go.

### From Source

In case you're on MacOS or the provided executables do not work (e.g. because of outdated versions of underlying libraries), you need to compile fom source:

1. Download `rustup` from [here](https://rustup.rs/) and follow installation instructions
2. Download the `src` folder, the `cargo.toml` file and put them in a directory as they are.
3. Open a terminal, navigate to the folder containing `src` and `cargo.toml`
4. Type `cargo install --path .` (Do not forget the dot)

Rust will now compile the binary and place it in your cargo root folder (location differs based on operating system). It can now be accessed from the command line without any further steps.

## Usage

As illustration what `pdbman` is capable of a typical workflow will be shown.

It is always called with a single PDB file as command line option like so:

![](/home/floeserb/Nextcloud/Code/rust/pdbman/example/pdbman_init.png)

As you can see from the changed prompt, a new shell environment has been started. It comes with a few convenience features, like a command history (accessible via the ' up' arrow key or by pressing `Ctrl + n`/`Ctrl + p` for the next or previous history item, respectively) or suggestions generated from previous commands.

In this shell all of the remaining commands can be entered. Help can be obtained by typing `help`/`--help` whereas `e`/`exit` will exit the shell. If any changes have been made to the file, they will be saved in a separate file.

First, the file of interest (here called `myfile.pdb`) will be analyzed:

`Y`

![](/home/floeserb/Nextcloud/Code/rust/pdbman/example/pdbman_y_initial.png)

This is typical for files that were preprocessed with Ambertools and needs to be changed.
Next, we look for potential clashes between atoms. This can happen if, e.g., a water shell has been added around the molecule and results in spurious high contributions from the involved Lennard-Jones-terms.

`Y -c`

![](/home/floeserb/Nextcloud/Code/rust/pdbman/example/pdbman_clash.png)

(Different input file shown in figure because `myfile.pdb` has no clashes)

This will look for any atoms not belonging to the same residue whose distance is smaller than 1.0 Angströms. If distances are significantly smaller than that they will be highlighted in yellow or red. For more granular control you can also search the local environment of a specific atom like so:

`Q -ts 2589 2`

which returns all atoms within 2 Angströms of the atom with ID 2589. Keep in mind that atoms within the same residue are not included here.

In order to make use of the PDB file with ORCA it's usually useful to clean all atoms from QM and active regions and start fresh:

`R`

It is also possible to only reset the QM or active region by giving the respective flags:

`R -q`/`R -a`

The changes made will be saved in a separate file upon exiting the `pdbman` shell so the input file will never be overwritten.

Next, we want to build an active region of residues around a metal ion in the center of the region of interest to us. In order to find the ID of a suitable atom you can do:

`Q -tl Cu`

![](/home/floeserb/Nextcloud/Code/rust/pdbman/example/pdbman_query.png)

This finds all atoms with the name `Cu` in the PDB file. The search is case-insensitive. You can also search for residues with the `-r` flag but be aware that the residue names can differ from the usually more intuitive atom names.

After we found the atom ID we need, we can build a spherical active space with a radius of our choosing (in this case 8 Angströms):

`A -ras 2589 8`

The radius argument can be any floating point number (i.e. decimal points are allowed). If the `-r` flag is given, all residues that contain an atom within the given radius to the given atom will be included. We can look at the results:

`Y`

or with more detail:

`Y -ra`

In general it is a good idea to keep checking your progress in case of unanticipated behaviour.
Next, we build a QM region. We have to decide on residues to include here which is usually done via inspection of the active site with some graphical program. Once you decide what to include you can use pdbman to apply it to your PDB file in a few commands.

In this case, we have a mononuclear Cu center, surrounded by several amino acids, a hydrogen peroxide and a polysaccharide. We determined by hand the IDs of the residues to include but in case we forget a number, we can use `pdbman` to query for it like so:

`Q -rl HIS`

This searches for all occurences of the residue name 'HIS' (input is case-insensitive) including the atoms in the residues.

In order to finally build our QM region we first add a couple of residues with all atoms:

`A -rql 1,171,460,203,204`

Since we only want the sidechains of some of the amino acids, we can choose to only add those:

`A -rqdl 85,87,160`

Note that in the case of GLY no atoms will be added to the chosen region. If residue insertion codes are present, these need to be included when residue list input is given like so:

`A -rql 85A,87A,9999A-4B`

This is mostly relevant if more than 9999 residues are present.

If we now want to make some more granular changes to the QM region we can add or remove specific atoms. First we query for their IDs:

`Q -rl 1`

The we remove the carbon and oxygen atom:

`R -tql 17,18`

The `-l` option supports ranges of atoms or residues which can be given with a dash or colon as separator like so:

`A -tql 1-8,19:27,5`

A final look at the region declarations we made:

`Y`

`Y -rq`

![](/home/floeserb/Nextcloud/Code/rust/pdbman/example/pdbman_ana_final.png)

Looks good! Now we're good to go for QM/MM jobs!

## Requirements

`Rust` > 1.55 since it uses the 2021 edition.

## Known Issues

If the number of residues is over 9999, automatically created PDB files will usually wrap around and start counting at 1 again. `pdbman` can read these but uses an internal numbering scheme that keeps on counting after 9999. 

This does **not** affect the printing of the file but has an effect on the querying and analysis functions as the residues will be displayed with their internal numbering, not with what is present in the input or output files.

## Contributor

Contributed by Benedikt Flöser
