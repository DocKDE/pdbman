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
- Write PDB structure and commands to recreate it from scratch to stdout or file
- Measure distances, angles or dihedrals between selected atoms
- Measure distances in a sphere around an atom
- Lists of QM or active atoms can be saved to or loaded from files

The several options are provided via command line flags (the ordering of the flags does not matter). 
Additional information can be obtained by giving the `--help`/`-h` option anywhere in the program.

---

### Modes of operation

There are three ways to use `pdbman`:

1. shell mode

2. directly from command line with input given as arguments

3. directly from command line with input loaded from a file

Depending on what you want to use the program for, each of these may come in handy. 

#### Shell mode

For workflows where you want to enter several commands in quick succession the shell mode is recommended. 
It is started simply by calling `pdbman` with the file you want to work with as argument and the 
`--interactive` option (or `-i` in short):

```
pdbman myfile.pdb -i
```

This will drop you in a shell where you can then enter commands:

```
pdbman>
```

This is convenient because the PDB file has only to be loaded once and will be kept in memory as 
long as you're working with it. Furthermore, a few convenience functions are present:

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

For single queries or scripting purposes pdbman can be called from command line, even with multiple commands. 
For this, call the program with the PDB file you want to use as first argument and any commands after that. 

If several commands are to be run in succession, they can be chained by using `/` as a separator like so:

```
pdbman myfile.pdb command1 / command2 / command3
```

This makes it easier to use pdbman in scripts but will also make the output less legible when calling several 
queries in succession.

#### Load commands from file

If many commands are to be executed in an automated fashion, they can be saved in a file and called from there. 
In this case the `--file` (short: `-f`) option needs to be given, followed by the name of the file which holds the commands to be executed:

```
pdbman myfile.pdb -f commands.txt
```

`pdbman` expects one command per line in the given file:

```
command1
command2
command2
```

It is also possible to utilize this mode of operation for scripting in which case a file holding 
the commands needs to be present.
The primary use case for this is to transfer the state of one PDB structure to another in a
different file. To achieve this the commands to recreate the PDB structure need to be written to a
file (see below) and read in to be applied to a different file with the command given above.

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

```
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

```
# Remove all atoms from all regions
r
# Remove all atoms from QM1 region
r -q
# Remove all atoms from QM2 region
r -o
```
If specific atoms or residues are to be removed, these need to be selected after giving the region
flag appropriate to the region from which the atoms/residues are to be removed.

The selection syntax for atoms or residues is keyword-based and accepts six of them (case-insensitive). Each of them selects something different:
- id -> atoms by ID 
- resid -> residues by ID
- name -> atoms by name
- resn(ame) -> residues by name
- s(phere) -> atoms in a sphere around a central atom
- ressphere/rs -> whole residues in a sphere around a central atom

Each keyword needs to be followed by appropriate input to select for as given in the following examples.

```
# Remove atoms 16 and 17 from QM1 region
r -q id 16,17
# Remove residues with the name 'HIS' from QM2 region
r -o resn his
# Remove given list of residues from active region. Note that ranges are supported
r -a resid 1-12,49,128,3:5
# Remove all atoms within 10 Å of the atom with ID 3230 from the active region
r -a sphere 3230 10
# Remove all residues that have an atom within 6 Å of the atom with ID 3230 from the QM1 region
r -q rs 3230 6
```

The selections can be chained for more finegrained control:
```
# Remove the C atoms of GLY residues from QM1 region
r -q name c and resn gly
# Remove all waters combining different naming schemes for it
r -a resn wat or resn hoh or resn h2o
```
For chaining either the `and`/`or` keywords or the equivalent `&`/`|` operators can be used.

Inverting a selection is also possible:
```
# Remove all non-Water atoms from QM1 region
r -q not resname wat and not resname hoh
```
Instead of the `not` keyword an exclamation mark (`!`) can also be given.
As many selections as desired can be chained, note that they are left-associate, i.e. they will be
evaluated from left to right.

#### Add

The syntax is exactly the same as for the `Remove` subcommand, except that is does not accept blanket additions of whole regions. 

Examples:

```
# Same syntax as for 'Remove' (see above)
a -q id 16,17
a -a resid 1-12,49,128
a -o name cu or resname glu
a -a s 3230 10
```

#### Query

This command allows getting information on atoms and residues in the PDB file. As the `Add` and `Remove` commands 
the selection of which atoms/residues is keyword-driven and follows the same syntax.

This way, atoms/residues with a specific name or ID can be targeted or the surrounding atoms/residues of a specific atom.

Examples:

```
# Find atom(s) with name 'cu' (case-insensitive)
q name cu
# Show atoms with IDs 23-30
q id 23:30
# Show residues with IDs 1 and 2
q resid 1,2
# Show all atoms within 2.5 Å of atom with ID 2589
q s 2589 2.5
```

If the number of residues queried is one or the atoms queried all belong to the same residue 
and the respective residue is an amino acid, an ASCII representation of it will be shown:

```
pdbman> q id 1

                    |                 HE1
                  H-N             __  /
                    |   HB1    ND1--CE1
                    |   |     /      |
                 HA-CA--CB--CG       |
                    |   |     \\     |
                    |   HB2    CD2--NE2
                  O=C           |     \
                    |          HD2    HE2

╭─────────┬───────────┬────────────┬──────────────┬────┬────────╮
│ Atom ID │ Atom name │ Residue ID │ Residue Name │ QM │ Active │
╞═════════╪═══════════╪════════════╪══════════════╪════╪════════╡
│ 1       │ N         │ 1          │ HIE          │ 2  │ 1      │
╰─────────┴───────────┴────────────┴──────────────┴────┴────────╯
```

#### Measure

The measure command takes either an `--atom`/`-a` or a `--sphere`/`-s` option to indicate
what is to be measured. The `--atom` option expects two, three or four atom indices as
arguments which will prompt `pdbman` to calculate the distance, angle or dihedral of the
given atoms, depending on how many indices were given.

Note that the order of the atoms is of importance for angles and dihedrals.

If the `--sphere` option is given, an atom ID and a radius are expected as arguments.
`pdbman` then print a table of all atoms within the given radius around the given atom
with the respective distances.

Examples:

```
# Get distance between atom ID 2 and 8
m -a 2 8
# Get dihedral given by atom IDs 1,2,3 and 4
m -a 1 2 3 4
# Show atoms within 3 Å of atom ID 2589
m -s 2589 3
```

#### Write

This command will write information about the current state of the PDB structure held in 
memory to stdout or a file. This can either be the PDB structure itself or a list of 
commands that reproduce the current state.

 If no further options are given, the whole PDB structure will be written to stdout. 
 The output target can be modified by providing a `--overwrite`/`-w` or `--file`/`-f` flag. 
 These will write the output to the input PDB file or a specified file. Note that in the 
 first case, the input PDB will be overwritten!

If the state flag (`--state`/`-s`) is given, a list of commands to recreate the current 
state of the PDB file will be written to stdout.
If, additionally, the `--file`/`-f` option is given followed by a file path, the output 
will be written to the given file.

This is useful to quickly transfer the state of one PDB file to another. 

Examples:

```
# Write PDB structure to input file, overwriting it
w -w
# Write commands to recreate state to stdout
w -s
# Write commands to recreate state to file
w -sf commands.txt
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

In this shell all of the remaining commands can be entered. Help can be obtained by typing 
`h`/`help`/`--help` whereas `e`/`exit` will exit the shell. Note that no changes made will 
be saved unless explicitly requested by the user.

First, the file of interest (here called `myfile.pdb`) will be analyzed:

```
pdbman> y
╭────────┬────────────┬───────────────╮
│        │ # of Atoms │ # of Residues │
╞════════╪════════════╪═══════════════╡
│ QM1    │ 29178      │ 5413          │
├────────┼────────────┼───────────────┤
│ QM2    │ 0          │ 0             │
├────────┼────────────┼───────────────┤
│ Active │ 0          │ 0             │
├────────┼────────────┼───────────────┤
│ Total  │ 29178      │ 5413          │
╰────────┴────────────┴───────────────╯
```

This is typical for files that were preprocessed with Ambertools and needs to be changed.
Next, we look for potential clashes between atoms. This can happen if, e.g., a water shell 
has been added around the molecule and results in spurious high contributions from the 
involved Lennard-Jones-terms.

```
pdbman> y -c
╭────────┬────────────┬───────────────╮
│        │ # of Atoms │ # of Residues │
╞════════╪════════════╪═══════════════╡
│ QM1    │ 0          │ 0             │
├────────┼────────────┼───────────────┤
│ QM2    │ 0          │ 0             │
├────────┼────────────┼───────────────┤
│ Active │ 0          │ 0             │
├────────┼────────────┼───────────────┤
│ Total  │ 83         │ 7             │
╰────────┴────────────┴───────────────╯

Clash Analysis
╭────────────────────────────────────────────────────────────────────────────────────────────────╮
│ Atom ID 1   Atom Name 1   Residue Name 1   Atom ID 2   Atom Name 2   Residue Name 2   Distance │
╞════════════════════════════════════════════════════════════════════════════════════════════════╡
│ 14          HE2           HIE              81          O             WAT              0.66     │
╰────────────────────────────────────────────────────────────────────────────────────────────────╯
```

(Different input file was used here because `myfile.pdb` has no clashes)

This will look for any atoms not belonging to the same residue whose distance is smaller 
than 1.0 Angströms and displays the resulting atom pairs ordered by distance. 
For more granular control you can also search the local 
environment of a specific atom like so:

```
pdbman> m -s 2589 2
```

which returns all atoms within 2 Angströms of the atom with ID 2589 ordered by distance. 

In order to make use of the PDB file with ORCA it's usually useful to clean all atoms from 
QM and active regions and start fresh:

```
pdbman> R
```

It is also possible to only reset the QM or active region by giving the respective flags:

```
pdbman> R -q/R -a
```

We will save the changes made here at the end of our editing session.

Next, we want to build an active region of residues around a metal ion in the center of 
the region of interest to us. In order to find the ID of a suitable atom you can do, e.g.:

```
pdbman> q name cu
╭─────────┬───────────┬────────────┬──────────────┬────┬────────╮
│ Atom ID │ Atom name │ Residue ID │ Residue Name │ QM │ Active │
╞═════════╪═══════════╪════════════╪══════════════╪════╪════════╡
│ 2589    │ CU        │ 171        │ CU           │ 1  │ 1      │
╰─────────┴───────────┴────────────┴──────────────┴────┴────────╯
```

This finds all atoms with the name `Cu` in the PDB file. The search is case-insensitive. 
You can also search for residues with the `-r` flag but be aware that the residue names 
can differ from the usually more intuitive atom names.

After we found the atom ID we need, we can build a spherical active space with a radius 
of our choosing (in this case 8 Angströms):

```
pdbman> A -a s 2589 8
```

The radius argument can be any floating point number (i.e. decimal points are allowed). If 
the `-r` flag is given, all residues that contain an atom within the given radius to the 
given atom will be included. We can look at the results:

```
pdbman> y
╭────────┬────────────┬───────────────╮
│        │ # of Atoms │ # of Residues │
╞════════╪════════════╪═══════════════╡
│ QM1    │ 0          │ 0             │
├────────┼────────────┼───────────────┤
│ QM2    │ 0          │ 0             │
├────────┼────────────┼───────────────┤
│ Active │ 243        │ 36            │
├────────┼────────────┼───────────────┤
│ Total  │ 29178      │ 5413          │
╰────────┴────────────┴───────────────╯
```

or with more detail (output not shown here):

```
pdbman> Y -ra
```

In general it is a good idea to keep checking your progress in case of unanticipated behaviour.
Next, we build a QM region. We have to decide on residues to include here which is usually 
done via inspection of the active site with some graphical program. Once you decide what 
to include you can use `pdbman` to apply it to your PDB file in a few commands.

In this case, we have a mononuclear Cu center, surrounded by several amino acids, a 
hydrogen peroxide and a polysaccharide. We determined by hand the IDs of the residues 
to include but in case we forget a number, we can use `pdbman` to query for it like so:

```
pdbman> Q resn HIS
```

This searches for all occurences of the residue name 'HIS' (input is case-insensitive) 
including the atoms in the residues.

In order to finally build our QM region we first add a couple of residues with all atoms:

```
pdbman> A -q resid 1,171,460,203,204
```

Since we only want the sidechains of some of the amino acids, we can choose to only 
add those:

```
pdbman> A -qd resid 85,87,160
```

Note that in the case of GLY no atoms will be added to the chosen region, if the 
`--sidechain`/`-d` flag is given. 

If we now want to make some more granular changes to the QM region we can add or 
remove specific atoms. First we query for their IDs:

```
pdbman> q resid 1

                    |                 HE1
                  H-N             __  /
                    |   HB1    ND1--CE1
                    |   |     /      |
                 HA-CA--CB--CG       |
                    |   |     \\     |
                    |   HB2    CD2--NE2
                  O=C           |     \
                    |          HD2    HE2

╭─────────┬───────────┬────────────┬──────────────┬────┬────────╮
│ Atom ID │ Atom name │ Residue ID │ Residue Name │ QM │ Active │
╞═════════╪═══════════╪════════════╪══════════════╪════╪════════╡
│ 1       │ N         │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 2       │ H1        │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 3       │ H2        │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 4       │ CA        │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 5       │ HA        │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 6       │ CB        │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 7       │ HB2       │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 8       │ HB3       │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 9       │ CG        │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 10      │ ND1       │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 11      │ CE1       │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 12      │ HE1       │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 13      │ NE2       │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 14      │ HE2       │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 15      │ CD2       │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 16      │ HD2       │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 17      │ C         │ 1          │ HIE          │ 0  │ 1      │
├─────────┼───────────┼────────────┼──────────────┼────┼────────┤
│ 18      │ O         │ 1          │ HIE          │ 0  │ 1      │
╰─────────┴───────────┴────────────┴──────────────┴────┴────────╯
```

The we remove the carbon and oxygen atom:

```
pdbman> R -q id 16,17
```

The `id` and `resid` keywords support ranges of atoms or residues which can be given with a dash or colon 
as separator like so:

```
pdbman> A -q id 1-8,19:27,5
```

A final look at the region declarations we made:

```
pdbman> y -rq
╭────────┬────────────┬───────────────╮
│        │ # of Atoms │ # of Residues │
╞════════╪════════════╪═══════════════╡
│ QM1    │ 122        │ 8             │
├────────┼────────────┼───────────────┤
│ QM2    │ 0          │ 0             │
├────────┼────────────┼───────────────┤
│ Active │ 243        │ 36            │
├────────┼────────────┼───────────────┤
│ Total  │ 29178      │ 5413          │
╰────────┴────────────┴───────────────╯

QM1 Residues
╭────────────┬──────────────┬────────────┬────────────────╮
│ Residue ID │ Residue Name │ # of Atoms │ # of QM1 Atoms │
╞════════════╪══════════════╪════════════╪════════════════╡
│ 1          │ HIE          │ 18         │ 16             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 85         │ ALA          │ 10         │ 10             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 87         │ HID          │ 17         │ 17             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 160        │ PHE          │ 20         │ 20             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 171        │ CU           │ 1          │ 1              │
├────────────┼──────────────┼────────────┼────────────────┤
│ 203        │ 4YB          │ 27         │ 27             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 204        │ 4YB          │ 27         │ 27             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 460        │ PER          │ 4          │ 4              │
╰────────────┴──────────────┴────────────┴────────────────╯
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

All commands shown in the previous section can also be used directly from command line. 
Typically this would be desirable if no analysis or querying is necessary and several 
commands can be chained. One use case would be to transfer the state of one PDB file to another:

First, the state needs to be saved to a file:

```
pdbman myfile.pdb w -sf state.txt
```

Next, the state can be imported into another file like this:

```
pdbman otherfile.pdb -f state.txt
```

Note that this will delete the previous state of the chosen PDB file and overwrite it with the
imported one.
Now the state of `myfile.pdb` has been transferred to `otherfile.pdb` successfully. 

Finally, for quick queries and analyses `pdbman` may be called from command line directly as above:

```
pdbman myfile.pdb y
```

which will output something like this (just as in shell mode, see above):

```
pdbman> y -rq
╭────────┬────────────┬───────────────╮
│        │ # of Atoms │ # of Residues │
╞════════╪════════════╪═══════════════╡
│ QM1    │ 122        │ 8             │
├────────┼────────────┼───────────────┤
│ QM2    │ 0          │ 0             │
├────────┼────────────┼───────────────┤
│ Active │ 243        │ 36            │
├────────┼────────────┼───────────────┤
│ Total  │ 29178      │ 5413          │
╰────────┴────────────┴───────────────╯

QM1 Residues
╭────────────┬──────────────┬────────────┬────────────────╮
│ Residue ID │ Residue Name │ # of Atoms │ # of QM1 Atoms │
╞════════════╪══════════════╪════════════╪════════════════╡
│ 1          │ HIE          │ 18         │ 16             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 85         │ ALA          │ 10         │ 10             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 87         │ HID          │ 17         │ 17             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 160        │ PHE          │ 20         │ 20             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 171        │ CU           │ 1          │ 1              │
├────────────┼──────────────┼────────────┼────────────────┤
│ 203        │ 4YB          │ 27         │ 27             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 204        │ 4YB          │ 27         │ 27             │
├────────────┼──────────────┼────────────┼────────────────┤
│ 460        │ PER          │ 4          │ 4              │
╰────────────┴──────────────┴────────────┴────────────────╯
```

However, if several queries or analyses are to be run in succession, it is recommended to use shell mode to avoid 
parsing the PDB file with every call of `pdbman`. Especially for larger files this quickly adds up.

## Requirements

`Rust` > 1.55 since it uses the 2021 edition.

## Known Issues

If the number of residues is over 9999 or the number of atoms is over 99999, automatically created PDB files 
will usually wrap around and start counting at 1 again. `pdbman` can read these but uses an internal 
numbering scheme that keeps on counting after reaching the wraparound point. 

This discrepancy exists because the number of digits for atom and residues IDs is limited by the 
file format but these are used by `pdbman` to distiguish between unique atoms and residues. 

This does **not** affect the reading or printing of the file but has an effect on the querying 
and analysis functions as the residues will be displayed with their internal numbering, not 
with what is present in the input or output files.

## Contributor

Contributed by Benedikt Flöser
