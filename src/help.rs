pub const HELP_LONG: &str = "
pdbman 0.9.2
    
Benedikt M. Fl\u{f6}ser <benedikt.floeser@cec.mpg.de>
    
Analyzes and edits PDB files for usage in QM/MM calculations with the ORCA Quantum Chemistry package

USAGE:
    pdbman <PDBFILE> <[OPTIONS]|[SUBCOMMAND]>
    
ARGS:
    <PDBFILE>    Path to PDB file
    
OPTIONS:
    -f, --file <File>    Read commands from file
    -h, --help           Display help message
    -i, --interactive    Enter interactive mode

SUBCOMMANDS:
    Add                  Add atoms or residues to QM1/QM2/Active region
    Remove               Remove atoms or residues from QM1/QM2/Active region
    Query                Query atoms or residues
    Analyse              Analyze PDB file and QM1/QM2/Active region
    Write                Write PDB structure information to file or stdout
    Measure              Measure distances, angles and dihedrals between atoms

Calling a subcommand with the '--help/-h' flag will display a help message for it";

pub const HELP_SHORT: &str = "
USAGE:
    pdbman <PDBFILE> <[OPTIONS]|[SUBCOMMAND]>
    
ARGS:
    <PDBFILE>    Path to PDB file
    
OPTIONS:
    -f, --file <File>    Read commands from file
    -h, --help           Display help message
    -i, --interactive    Enter interactive mode

SUBCOMMANDS:
    Add                  Add atoms or residues to QM1/QM2/Active region
    Remove               Remove atoms or residues from QM1/QM2/Active region
    Query                Query atoms or residues
    Analyse              Analyze PDB file and QM1/QM2/Active region
    Write                Write PDB structure information to file or stdout
    Measure              Measure distances, angles and dihedrals between atoms

Calling a subcommand with the '--help/-h' flag will display a help message for it";

pub const HELP_INTER: &str = "
USAGE:
    [SUBCOMMAND] 
    
SUBCOMMANDS:
    Add                  Add atoms or residues to QM1/QM2/Active region
    Remove               Remove atoms or residues from QM1/QM2/Active region
    Query                Query atoms or residues
    Analyse              Analyze PDB file and QM1/QM2/Active region
    Write                Write PDB structure information to file or stdout
    Measure              Measure distances, angles and dihedrals between atoms

Calling a subcommand with the '--help/-h' flag will display a help message for it";
