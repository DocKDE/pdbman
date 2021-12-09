# pdbman

Analyzes and modifies PDB files for QM/MM calculations.

## Description

This program takes a PDB file as input and analyzes or changes it as needed to
serve as input for Orca QM/MM calculations. 

### Capabilities

- Query specific atoms or residues by ID or name
- Find atoms or residues within a user-defined sphere of a given atom
- Analyze atoms and/or residues present in QM or active region
- Find clashes and atomic contacts in given PDB file
- Add or remove atoms or residues to QM or active region by ID or name
- Add or remove atoms and residues to QM or active region by calculating a sphere of given radius around a given atom
- Add or remove only the sidechain or backbone of residues
- Write PDB structure, atom or residue ID lists to stdout or file
- Lists of QM or active atoms can be saved to or loaded from files

The several options are provided via command line flags (the ordering of the flags does not matter). Additional information can be obtained by giving the `--help`/'`-h` option anywhere in the program.

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

This is convenient because the PDB file has only to be loaded once and will be kept in memory as long as you're working with it. Furthermore, a few convenience functions are present:

- Command history
  
  - Access last command with "up" arrow key
  
  - Reverse-search in command history with `Ctrl-R`
  
  - Cycle through command history with `Ctrl-p` and `Ctrl-n` (previous and next item, respectively)

- Undo/redo functionality
  
  - Type `undo` or `redo` to perform the respective actions
  
  - Only relevant for commands that induce changes in the PDB structure

- Command suggestions based on history
  
  - To accept a suggested command press the "right" arrow key

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

Five different subcommands are available:

- `Analyze` analyzes the QM and active region definitions currently in place
- `Remove` removes atoms or residues from QM or active region
- `Add` adds atoms or residues from QM or active regions
- `Query` queries the PDB file for information on atoms, residues, atom names etc.
- `Write` writes information about the current state of the PDB structure in memory to stdout or file

Each subcommand can be called by various aliases (list not exhaustive):

Analyze: `ana`, `Y`, `y`

Query: `Que`, `Q`, `q`

Remove: `rem`, `R`, `r`

Add: `Add`, `A`, `a`

Write: `Wri`, `W`, `w`

In general subcommand inputs are not case-sensitive and `pdbman` will try to infer the desired subcommand from any given abbreviation (which has to start at the beginning of the word). 

&nbsp;

#### Analyze

This mode will show the number of atoms and residues in the QM1, QM2 and active regions. If only a part of a residue is present in one region, the residue will still be counted.

If the `--residues`/`-r` or `--atoms`/`-t` flag is given in conjunction with a `--qm1`/`-q`, `--qm2`/`-o` or `--active`/`-a` flag the residues/atoms of the respective region will be listed.

Examples:

```Python
# Analyze QM1 residues
y -rq 
# Analyze active atoms
y -ta
```

If the `--clashes`/`-c` or `--contacts`/`-n` flag is given, van-der-Waals clashes or contacts will be listed, respectively, if present. Especially close contacts will be colorized.

#### Remove

This mode will remove atoms and/or residues from the specified region.

If no other flags are given, all atoms will be remove from all regions. If only a region flag (`--qm1`/`-q`, `--qm2`/`-o` or `--active`/`-a`) is given, all atoms will be removed from the specified region.

Examples:

```Python
# Remove all atoms from all regions
r
# Remove all atoms from QM1 region
r -q
# Remove all atoms from QM2 region
r -o
```

If an `--atoms`/`-t` or `--residues`/`-r` flag is also present, specific atoms/residues will be remove from the specified region. In this case a source of atoms/residues must be given. This can be a command line list given after a `--list`/`-l` flag, a list read from a file specified with `--file`/`-f` or a calculated sphere around a given atom specified with the `--sphere`/`-s` flag.

In the latter case an atom ID and a radius must be given. All atoms within the given radius around the given atom will be processed. If the `--residues`/`-r` flag is active, all residues that have at least one atom within the given sphere will be removed.

A `--list`/`-l` flag must be followed by a comma-separated list of items that are to be processed. This can consist of atom/residue IDs or names thereof, but not both.

`pdbman` recognizes ranges of lists separated with a colon (`:`) or dash (`-`).

If a list from a file is to be used with the `--file`/`-f` flag, a filename holding the list must be specified. It follows the same syntax rules as the command line list mentioned previously.

Examples:

```Python
# Remove atoms 16 and 17 from QM1 region
r -tql 16,17
# Remove given list of atoms from active region
r -ral 1-12,49,128
# Remove list of atoms given in 'atomlist.txt' from QM2 region
r -tof atomlist.txt
# Remove all atoms within 10 Å of the atom with ID 3230 from the active region
r -tas 3230 10
# Remove all residues that have an atom within 6 Å of the atom with ID 3230 from the QM1 region
r -rqs 3230 6
```

#### Add

The syntax is exactly the same as for the `Remove` subcommand, except that is does not accept blanket additions of whole regions. 

Examples:

```Python
# Same syntax as for 'Remove' (see above)
a -tql 16,17
a -ral 1-12,49,128
a -tof atomlist.txt
a -tas 3230 10
```

#### Query

This command allows getting information on atoms and residues in the PDB file. It has to be called with a target (`--atoms`/`-t` or `--residues`/`-r`) and a source (`--list`/`l`, `--file`/`-f` or `--sphere`/`-s`) which specify what is to be queried.

This way, atoms/residues with a specific name or ID can be targeted or the surrounding atoms/residues of a specific atom.

Examples:

```Python
# Find atom(s) with name 'cu' (case-insensitive)
q -tl cu
# Show atoms with IDs 23-30
q -tl 23:30
# Show residues with IDs 1 and 2
q -rl 1,2
# Show all atoms within 2.5 Å of atom with ID 2589
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

This command will write information about the current state of the PDB structure held in memory to stdout or a file. This can either be the PDB structure itself or a list of the atom or residue IDs in a given region.

 If no further options are given, the whole PDB structure will be written to stdout. The output target can be modified by providing a `--overwrite`/`-w` or `--file`/`-f` flag. These will write the output to the input PDB file or a specified file. Note that in the first case, the input PDB will be overwritten!

If a region flag (`--qm1`/`-q`, `--qm2`/`-o` or `--active`/`-a`) is given, the atom or residue IDs currently making up that region will be printed to stdout or a file, depending on the user-given flags. Accordingly, either an `--atoms`/`-t` or `--residues`/`-r` flag must be given. 

This is useful to quickly transfer the state of one PDB file to another. Note that if the `--residues` flag is given, all residue IDs containing at least one atom in the given region will be printed so for this use case, using a list of atoms is the safe option.

Examples:

```Python
# Write PDB structure to input file, overwriting it
w -w
# Write list of QM1 atoms to stdout
w -qt
# Write list of active residues to file 'activeatoms.txt'
w -arf activeatoms.txt
```

---

Help messages are available for all subcommands like so:

`pdbman [Subcommand] --help`

Or from within the shell:

```
pdbman> [Subcommand] --help
```



## Installation

### Binary

Linux and Windows executables are provided (`pdbman` and `pdbman.exe`, respectively). Just put them somewhere in your $PATH and you're good to go.

### From Source

In case you're on MacOS or the provided executables do not work (e.g. because of outdated versions of underlying libraries), you need to compile fom source:

1. Download `rustup` from [here](https://rustup.rs/) and follow installation instructions.
2. Download or clone the repository.
3. Open a terminal, navigate to the folder containing the `src` folder and the `Cargo.toml` file.
4. Type `cargo install --path .` (Do not forget the dot).

Rust will now compile the binary and place it in your cargo root folder (location differs based on operating system). It can now be accessed from the command line without any further steps. 

It is also possible to just build the binary without installing it anywhere on your system so you can handle it as you see fit. In this case do:

```
cargo build --release
```

After a successful build a release binary will be placed in the `./target/release` folder (relative to the root folder of the repository).

## Example Usage

As illustration of what `pdbman` is capable of typical workflows will be shown.

### Shell mode

First the shell needs to be started like so:

```
pdbman myfile.pdb -i

pdbman>
```

In this shell all of the remaining commands can be entered. Help can be obtained by typing `h`/`help`/`--help` whereas `e`/`exit` will exit the shell. Note that no changes made will be saved unless explicitly requested by the user.

First, the file of interest (here called `myfile.pdb`) will be analyzed:

```
pdbman> y
+--------+------------+---------------+
|        | # of Atoms | # of Residues |
+--------+------------+---------------+
| QM1    | 87449      | 28301         |
+--------+------------+---------------+
| QM2    | 0          | 0             |
+--------+------------+---------------+
| Active | 0          | 0             |
+--------+------------+---------------+
```

This is typical for files that were preprocessed with Ambertools and needs to be changed.
Next, we look for potential clashes between atoms. This can happen if, e.g., a water shell has been added around the molecule and results in spurious high contributions from the involved Lennard-Jones-terms.

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

(Different input file was used here because `myfile.pdb` has no clashes)

This will look for any atoms not belonging to the same residue whose distance is smaller than 1.0 Angströms. If distances are significantly smaller than that they will be highlighted in yellow or red. For more granular control you can also search the local environment of a specific atom like so:

```
pdbman> Q -ts 2589 2
```

which returns all atoms within 2 Angströms of the atom with ID 2589. Keep in mind that atoms within the same residue are not included here.

In order to make use of the PDB file with ORCA it's usually useful to clean all atoms from QM and active regions and start fresh:

```
pdbman> R
```

It is also possible to only reset the QM or active region by giving the respective flags:

```
pdbman> R -q/`R -a`
```

We will save the changes made here at the end of our editing session.

Next, we want to build an active region of residues around a metal ion in the center of the region of interest to us. In order to find the ID of a suitable atom you can do, e.g.:

```
pdbman> Q -tl cu
+---------+-----------+------------+--------------+----+--------+
| Atom ID | Atom name | Residue ID | Residue Name | QM | Active |
+---------+-----------+------------+--------------+----+--------+
| 2589    | CU        | 171        | CU           | 1  | 1      |
+---------+-----------+------------+--------------+----+--------+
```

This finds all atoms with the name `Cu` in the PDB file. The search is case-insensitive. You can also search for residues with the `-r` flag but be aware that the residue names can differ from the usually more intuitive atom names.

After we found the atom ID we need, we can build a spherical active space with a radius of our choosing (in this case 8 Angströms):

```
pdbman> A -ras 2589 8
```

The radius argument can be any floating point number (i.e. decimal points are allowed). If the `-r` flag is given, all residues that contain an atom within the given radius to the given atom will be included. We can look at the results:

```
pdbman> y
+--------+------------+---------------+
|        | # of Atoms | # of Residues |
+--------+------------+---------------+
| QM1    | 0          | 0             |
+--------+------------+---------------+
| QM2    | 0          | 0             |
+--------+------------+---------------+
| Active | 441        | 37            |
+--------+------------+---------------+
```

or with more detail (output not shown here):

```
pdbman> Y -ra
```

In general it is a good idea to keep checking your progress in case of unanticipated behaviour.
Next, we build a QM region. We have to decide on residues to include here which is usually done via inspection of the active site with some graphical program. Once you decide what to include you can use `pdbman` to apply it to your PDB file in a few commands.

In this case, we have a mononuclear Cu center, surrounded by several amino acids, a hydrogen peroxide and a polysaccharide. We determined by hand the IDs of the residues to include but in case we forget a number, we can use `pdbman` to query for it like so:

```
pdbman> Q -rl HIS
```

This searches for all occurences of the residue name 'HIS' (input is case-insensitive) including the atoms in the residues.

In order to finally build our QM region we first add a couple of residues with all atoms:

```
pdbman> A -rql 1,171,460,203,204
```

Since we only want the sidechains of some of the amino acids, we can choose to only add those:

```
pdbman> A -rqdl 85,87,160
```

Note that in the case of GLY no atoms will be added to the chosen region, if the `--sidechain`/`-d` flag is given. 

If we now want to make some more granular changes to the QM region we can add or remove specific atoms. First we query for their IDs:

```
pdbman> Q -rl 1

                    |                 HE1
                  H-N             __  /
                    |   HB1    ND1--CE1
                    |   |     /      |
                 HA-CA--CB--CG       |
                    |   |     \\     |
                    |   HB2    CD2--NEM
                  O=C           |     \
                    |          HD2     CME--HM3
                                      /  \
                                    HM1  HM2

+---------+-----------+------------+--------------+----+--------+
| Atom ID | Atom name | Residue ID | Residue Name | QM | Active |
+---------+-----------+------------+--------------+----+--------+
| 1       | N         | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 2       | H1        | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 3       | H2        | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 4       | CA        | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 5       | HA        | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 6       | CB        | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 7       | HB2       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 8       | HB3       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 9       | CG        | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 10      | ND1       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 11      | CE1       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 12      | HE1       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 13      | NEM       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 14      | CD2       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 15      | HD2       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 16      | C         | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 17      | O         | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 18      | CME       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 19      | HM1       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 20      | HM2       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
| 21      | HM3       | 1          | HIM          | 1  | 0      |
+---------+-----------+------------+--------------+----+--------+
```

The we remove the carbon and oxygen atom:

```
pdbman> R -tql 16,17
```

The `-l` option supports ranges of atoms or residues which can be given with a dash or colon as separator like so:

```
pdbman> A -tql 1-8,19:27,5
```

A final look at the region declarations we made:

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

Looks good! Now we only have to save our progress.

```
pdbman> w -w
```

This will overwrite the input file, alternatively we can do:

```
pdbman> w -f output.pdb
```

which will save the edited PDB structure to `output.pdb` leaving the input file untouched.

Now we're finally ready for some QM/MM calculations with ORCA!

### Command line mode

All commands shown in the previous section can also be used directly from command line. Typically this would be desirable if no analysis or querying is necessary and several commands can be chained. One use case would be to transfer the state of one PDB file to another:

First, the atoms in the respective sections of the source PDB file need to be saved to a file each:

```
pdbman myfile.pdb w -qtf qmatoms.txt / w -atf activeatoms.txt
```

Chaining the commands like this has the advantage that the PDB file is only read and parsed once. Next, the files can be read by `pdbman` to serve as input for another PDB file:

```
pdbman otherfile.pdb r / a -qtf qmatoms.txt / a -atf activeatoms.txt / w -w
```

Make sure to clean the QM and active regions beforehand to remove potential leftovers. Also note that no changes will be written to the file unless specifically requested so a call to the `Write` subcommand is necessary.

Now the state of `myfile.pdb` has been transferred to `otherfile.pdb` successfully. The commands for `pdbman` could also have been read from a file if desired and both modes of operation are useful for scripting and automation workflows.

Finally, for quick queries and analyses `pdbman` may be called from command line directly as above:

```
pdbman myfile.pdb y
```

which will output something like this (just as in shell mode, see above):

```
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

However, if several queries or analyses are to be run in succession, it is recommended to use shell mode to avoid parsing the PDB file with every call of `pdbman`. Especially for larger files this quickly adds up.

## Requirements

`Rust` > 1.55 since it uses the 2021 edition.

## Known Issues

If the number of residues is over 9999 or the number of atoms is over 99999, automatically created PDB files will usually wrap around and start counting at 1 again. `pdbman` can read these but uses an internal numbering scheme that keeps on counting after reaching the wraparound point. 

This discrepancy exists because the number of digits for atom and residues IDs is limited by the file format but these are used by `pdbman` to distiguish between unique atoms and residues. 

This does **not** affect the reading or printing of the file but has an effect on the querying and analysis functions as the residues will be displayed with their internal numbering, not with what is present in the input or output files.

## Contributor

Contributed by Benedikt Flöser
