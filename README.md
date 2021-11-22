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
- Write changes to stdout or file
- Lists of QM or active atoms can be saved to or loaded from files

The several options are provided via command line flags (the ordering of the flags does not matter). Additional information can be obtained from `help`

---

### Modes of operation

There are three ways to use `pdbman`:

1. shell mode

2. directly from the command line

3. with input loaded from a file

Depending on what you want to use the program for, each of these may come in handy. 

#### Shell mode

For workflows where you want to enter several commands in quick succession the shell mode is recommended. It is started simply by calling `pdbman` with the file you want to work with as argument and the `--interactive` option (or `-i` in short):

```
pdbman myfile.pdb -i
```

This will drop you in a shell where you can then enter commands:

```
pdbman>
```

This is convenient because the PDB file has only to be loaded once and will be kept in memory as long as you're working with it. Furthermore, a few convenience functions are present, such as a command history of the session, comman suggestions based theron and tab completion for files.

For quitting the shell type `exit`, `e`, `quit` or press `Ctrl-C` or `Ctrl-D`.

#### Command line

For single queries or scripting purposes pdbman can be called from command line, even with multiple commands. For this, call the program with the PDB file you want to use as first argument and any commands after that. 

If several commands are to be run in succession, they can be chained by using `/` as a separator like so:

```
pdbman myfile.pdb command1 / command2 / command3
```

This makes it easier to use pdbman in scripts but will also make the output less legible when calling several queries in succession.

#### Load commands from file

If many commands are to be executed in an automated fashion, they can be saved in a file and called from there. In this case the `--File` (short: `-f`) option needs to be given, followed by the name of the file which holds the commands to be executed:

```
pdbman myfile.pdb -f mycommands.txt
```

`pdbman` expects one command per line in the given file:

```
command1
command2
command2
```

It is also possible to utilize this mode of operation for scripting in which case a file holding the commands needs to be present.

---

### Commands

Five different subcommands are available (equivalent aliases given in parentheses):

- `Analyze` (`analyze`, `ana`, `Y`, `y`), analyzes the QM and active region definitions currently in place
- `Remove` (`remove`, `rem`, `R`, `r`), removes atoms or residues from QM or active region
- `Add` (`add`, `A`, `a`), adds atoms or residues from QM or active regions
- `Query` (`query`, `que`, `Q`, `q`), queries the PDB file for information on atoms, residues, atom names etc.
- `Write` (`write`, `wri`, `W`, `w`), writes information about the current state of the PDB structure in memory to stdout or file

Each subcommand can be called by any of its aliases like so:

`analyze`

`rem`

`A -rql 1,5,18`

&nbsp;

#### Analyze

This mode will show the number of atoms and residues in the QM1, QM2 and active regions. If only a part of a residue is present in one region, the residue will still be counted.

If the `--residues`/`-r` or `--atoms`/`-t` flag is given in conjunction with a `--qm1`/`-q`, `--qm2`/`-o` or `--active`/`-a` flag the residues/atoms of the respective region will be listed.

Examples:

```
y -rq
y -ta
```

If the `--clashes`/`-c` or `--contacts`/`-n` flags are given, van-der-Waals clashes or contacts will be listed, respectively, if present. Especially close contacts will be colorized.

#### Remove

This mode will remove atoms and/or residues from the specified region.

If no other flags are given, all atoms will be remove from all regions. If only a region flag (``--qm1`/`-q`, `--qm2`/`-o` or `--active`/`-a) is given, all atoms will be removed from the specified region.

Example:

```
r
r -q
r -o
```

If an `--atoms`/`-t` or `--residues`/`-r` flag is also present, specific atoms/residues will be remove from the specified region. In this case a source of atoms/residues must be given. This can be a command line list given after a `--list`/`-l` flag, a list read from a file specified with `--file`/`-f` or a calculated sphere around a given atom specified with the `--sphere`/`-s` flag.

In the latter case an atom ID and a radius must be given. All atoms within the given radius around the given atom will be processed. If the `--residues`/`-r` flag is active, all residues that have at least one atom within the given sphere will be removed.

A `--list`/`-l` flag must be followed by a comma-separated list of items that are to be processed. This can consist of atom/residue IDs or names thereof, but not both.

`pdbman` recognizes ranges of lists separated with a colon (`:`) or dash (`-`).

Examples:

```
r -tql 16,17
r -ral 1-12,49,128
r -tof atomlist.txt
r -tas 3230 10
```

#### Add

The syntax is exactly the same as for the `Remove` subcommand, except that is does not accept blanket additions of whole regions. 

Examples:

```
a -tql 16,17
a -ral 1-12,49,128
a -tof atomlist.txt
a -tas 3230 10
```

#### Query

This command allows getting information on atoms and residues in the PDB file. It has to be called with a target (`--atoms`/`-t` or `--residues`/`-r`) and a source (`--list`/`l`, `--file`/`-f` or `--sphere`/`-s`) which specify what is to be queried.

This way, atoms/residues with a specific name or ID can be targeted or the surrounding atoms/residues of a specific atom.

Examples:

```
q -tl cu
q -tl 23:30
q -rl 1,2
q -ts 2589 2.5
```

If the number of residues queried is one or the atoms queried all belong to the same residue and the respective residue is an amino acid, an ASCII representation of it will be shown:

```
pdbman> q -tl 1

                    |                 HE1
                  H-N             __  /
                    |   HB1    ND1--CE1
                    |   |     /      |
                 HA-CA--CB--CG       |
                    |   |     \\     |
                    |   HB2    CD2--NE2
                  O=C           |     \
                    |          HD2    HE2

+---------+-----------+------------+--------------+----+--------+
| Atom ID | Atom name | Residue ID | Residue Name | QM | Active |
+---------+-----------+------------+--------------+----+--------+
| 1       | N         | 1          | HIE          | 2  | 1      |
+---------+-----------+------------+--------------+----+--------+
```

#### Write

This command will write the current state of the PDB structure held in memory to stdout or a file. If no further options are given, the whole PDB structure will be written to stdout. The output target can be modified by providing a `--overwrite`/`-w` or `--file`/`-f` flag. These will write the output to the input PDB file or a specified file. Note that in the first case, the input PDB will be overwritten!

If a region flag (`--qm1`/`-q`, `--qm2`/`-o` or `--active`/`-a`) is given, the atom IDs currently making up that region will be printed to stdout or a file, depending on the user-given flag. This is useful to quickly transfer the state of one PDB file to another.

Examples:

```
w -w
w -q
w -af activeatoms.txt
```

---

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

## Example Usage

As illustration what `pdbman` is capable of a typical workflow will be shown.

It is always called with a single PDB file as command line option like so:

```
pdbman myfile.pdb

pdbman>
```

In this shell all of the remaining commands can be entered. Help can be obtained by typing `help`/`--help` whereas `e`/`exit` will exit the shell. If any changes have been made to the file, they will be saved in a separate file.

First, the file of interest (here called `myfile.pdb`) will be analyzed:

`Y`

```
pdbman> y
+--------+------------+---------------+
|        | # of Atoms | # of Residues |
+--------+------------+---------------+
| QM1    | 116        | 9             |
+--------+------------+---------------+
| QM2    | 1          | 1             |
+--------+------------+---------------+
| Active | 294        | 21            |
+--------+------------+---------------+
```

This is typical for files that were preprocessed with Ambertools and needs to be changed.
Next, we look for potential clashes between atoms. This can happen if, e.g., a water shell has been added around the molecule and results in spurious high contributions from the involved Lennard-Jones-terms.

`Y -c`

```
pdbman> y -c
+--------+------------+---------------+
|        | # of Atoms | # of Residues |
+--------+------------+---------------+
| QM1    | 0          | 0             |
+--------+------------+---------------+
| QM2    | 0          | 0             |
+--------+------------+---------------+
| Active | 0          | 0             |
+--------+------------+---------------+

Clash Analysis
+-----------+-------------+----------------+-----------+-------------+----------------+----------+
| Atom ID 1 | Atom Name 1 | Residue Name 1 | Atom ID 2 | Atom Name 2 | Residue Name 2 | Distance |
+-----------+-------------+----------------+-----------+-------------+----------------+----------+
| 14        | HE2         | HIE            | 81        | O           | WAT            | 0.66     |
+-----------+-------------+----------------+-----------+-------------+----------------+----------+
```

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

```
pdbman> q -tl cu
+---------+-----------+------------+--------------+----+--------+
| Atom ID | Atom name | Residue ID | Residue Name | QM | Active |
+---------+-----------+------------+--------------+----+--------+
| 2589    | CU        | 171        | CU           | 1  | 1      |
+---------+-----------+------------+--------------+----+--------+
```

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

```
pdbman> y -rq
+--------+------------+---------------+
|        | # of Atoms | # of Residues |
+--------+------------+---------------+
| QM1    | 104        | 8             |
+--------+------------+---------------+
| QM2    | 0          | 0             |
+--------+------------+---------------+
| Active | 441        | 37            |
+--------+------------+---------------+

QM1 Residues
+------------+--------------+------------+----------------+
| Residue ID | Residue Name | # of Atoms | # of QM1 Atoms |
+------------+--------------+------------+----------------+
| 1          | HIE          | 18         | 16             |
+------------+--------------+------------+----------------+
| 85         | ALA          | 10         | 4              |
+------------+--------------+------------+----------------+
| 87         | HID          | 17         | 11             |
+------------+--------------+------------+----------------+
| 160        | PHE          | 20         | 14             |
+------------+--------------+------------+----------------+
| 171        | CU           | 1          | 1              |
+------------+--------------+------------+----------------+
| 203        | 4YB          | 27         | 27             |
+------------+--------------+------------+----------------+
| 204        | 4YB          | 27         | 27             |
+------------+--------------+------------+----------------+
| 460        | PER          | 4          | 4              |
+------------+--------------+------------+----------------+
```

Looks good! Now we're good to go for QM/MM jobs!

## Requirements

`Rust` > 1.55 since it uses the 2021 edition.

## Known Issues

If the number of residues is over 9999 or the number of atoms is over 99999, automatically created PDB files will usually wrap around and start counting at 1 again. `pdbman` can read these but uses an internal numbering scheme that keeps on counting after reaching the limit. 

This discrepancy exists because the number of digits for atom and residues IDs is limited by the file format but these are used by `pdbman` to distiguish between unique atoms and residues. 

This does **not** affect the reading or printing of the file but has an effect on the querying and analysis functions as the residues will be displayed with their internal numbering, not with what is present in the input or output files.

## Contributor

Contributed by Benedikt Flöser
