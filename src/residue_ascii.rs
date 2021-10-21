use std::collections::HashMap;

lazy_static! {
    pub static ref RESIDUE_ASCII: HashMap<&'static str, &'static str> = {
        let mut m = HashMap::new();
        m.insert(
            "ALA",
            r"
                    |
                  H-N
                    |     HB1
                    |    /
                 HA-CA--CB-HB2
                    |    \
                    |     HB3
                  O=C
                    |  
        ",
        );
        m.insert(
            "ARG",
            r"
                    |                      HH11
                  H-N                       |
                    |   HB1 HG1 HD1 HE     NH1-HH12
                    |   |   |   |   |    //(+)
                 HA-CA--CB--CG--CD--NE--CZ
                    |   |   |   |         \
                    |   HB2 HG2 HD2        NH2-HH22
                  O=C                       |
                    |                      HH21
        ",
        );
        m.insert(
            "ASN",
            r"
                    |
                  H-N
                    |   HB1 OD1    HD21 (cis to OD1)
                    |   |   ||    /
                 HA-CA--CB--CG--ND2
                    |   |         \
                    |   HB2        HD22 (trans to OD1)
                  O=C
                    | 
        ",
        );
        m.insert(
            "ASP",
            r"
                    |
                  H-N
                    |   HB1   OD1
                    |   |    //
                 HA-CA--CB--CG
                    |   |    \
                    |   HB2   OD2(-)
                  O=C
                    |
        ",
        );
        let ash = r"
                    |
                  H-N
                    |   HB1   OD1
                    |   |    //
                 HA-CA--CB--CG
                    |   |    \
                    |   HB2   OD2-HD2
                  O=C
                    |
        ";
        m.insert("ASH", ash);
        m.insert("ASPP", ash);
        m.insert(
            "CYM",
            r"
                    |
                  H-N
                    |   HB1
                    |   |     -
                 HA-CA--CB--SG   (thiolate)
                    |   |     
                    |   HB2    
                  O=C
                    |
        ",
        );
        m.insert(
            "CYS",
            r"
                    |
                  H-N
                    |   HB1
                    |   |
                 HA-CA--CB--SG
                    |   |     \
                    |   HB2    HG1
                  O=C
                    |
        ",
        );
        m.insert(
            "CYX",
            r"
                    |
                  H-N
                    |   HB1
                    |   |
                 HA-CA--CB--SG
                    |   |     \
                    |   HB2  (other residue)
                  O=C
                    |
        ",
        );
        m.insert(
            "GLN",
            r"
                    |
                  H-N
                    |   HB1 HG1 OE1   HE21 (cis to OE1)
                    |   |   |   ||    /
                 HA-CA--CB--CG--CD--NE2
                    |   |   |         \
                    |   HB2 HG2       HE22 (trans to OE1)
                  O=C
                    |
        ",
        );
        m.insert(
            "GLU",
            r"
                    |
                  H-N
                    |   HB1 HG1   OE1
                    |   |   |    //
                 HA-CA--CB--CG--CD
                    |   |   |    \
                    |   HB2 HG2   OE2(-)
                  O=C
                    |
        ",
        );
        let glh = r"
                    |
                  H-N
                    |   HB1 HG1   OE1
                    |   |   |    //
                 HA-CA--CB--CG--CD
                    |   |   |    \
                    |   HB2 HG2   OE2-HE2
                  O=C
                    |
        ";
        m.insert("GLH", glh);
        m.insert("GLUP", glh);
        m.insert(
            "GLY",
            "
                    |
                  H-N
                    |
                    |
                HA2-CA--HA3
                    |
                    |
                  O=C
                    |
        ",
        );
        m.insert(
            "HIM",
            r"
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
        ",
        );
        m.insert(
            "HIS",
            r"
                    |          (HD1)   HE1
                  H-N           | __  /
                    |   HB1    ND1--CE1
                    |   |     /      ||
                 HA-CA--CB--CG       ||
                    |   |     \\     ||
                    |   HB2    CD2--NE2
                  O=C           |     \
                    |          HD2   (HE2)
        ",
        );
        let hid = r"
                    |          HD1    HE1
                  H-N           |     /
                    |   HB1    ND1--CE1
                    |   |     /      ||
                 HA-CA--CB--CG       ||
                    |   |     \\     ||
                    |   HB2    CD2--NE2
                  O=C           |
                    |          HD2
        ";
        m.insert("HID", hid);
        m.insert("HSD", hid);
        let hie = r"
                    |                 HE1
                  H-N             __  /
                    |   HB1    ND1--CE1
                    |   |     /      |
                 HA-CA--CB--CG       |
                    |   |     \\     |
                    |   HB2    CD2--NE2
                  O=C           |     \
                    |          HD2    HE2
        ";
        m.insert("HIE", hie);
        m.insert("HSE", hie);
        let hip = r"
                    |          HD1    HE1
                  H-N           |     /
                    |   HB1    ND1--CE1
                    |   |     /      ||
                 HA-CA--CB--CG       ||
                    |   |     \\     ||
                    |   HB2    CD2--NE2(+)
                  O=C           |     \
                    |          HD2    HE2
        ";
        m.insert("HIP", hip);
        m.insert("HSP", hip);
        m.insert(
            "ILE",
            r"
                    |    HG21 HG22
                  H-N      | /
                    |     CG2--HG23
                    |    /
                 HA-CA--CB-HB    HD1
                    |    \       /
                    |     CG1--CD--HD2
                  O=C    / \     \
                    | HG11 HG12  HD3
        ",
        );
        m.insert(
            "LEU",
            r"
                    |        HD11 HD12
                  H-N          | /
                    |   HB1   CD1--HD13
                    |   |    /
                 HA-CA--CB--CG-HG
                    |   |    \
                    |   HB2   CD2--HD23
                  O=C          | \
                    |        HD21 HD22
        ",
        );
        m.insert(
            "LYS",
            r"
                    |
                  H-N
                    |   HB1 HG1 HD1 HE1    HZ1
                    |   |   |   |   |     /
                 HA-CA--CB--CG--CD--CE--NZ--HZ2
                    |   |   |   |   |     \
                    |   HB2 HG2 HD2 HE2    HZ3
                  O=C
                    |
        ",
        );
        let lyn = r"
                    |
                  H-N
                    |   HB1 HG1 HD1 HE1    HZ1
                    |   |   |   |   |     /
                 HA-CA--CB--CG--CD--CE--NZ
                    |   |   |   |   |     \
                    |   HB2 HG2 HD2 HE2    HZ3
                  O=C
                    |
        ";
        m.insert("LYN", lyn);
        m.insert("LSN", lyn);
        m.insert(
            "MET",
            r"
                    |
                  H-N
                    |   HB1 HG1     HE1
                    |   |   |       |
                 HA-CA--CB--CG--SD--CE--HE3
                    |   |   |       |
                    |   HB2 HG2     HE2
                  O=C
                    |
        ",
        );
        m.insert(
            "PHE",
            r"
                    |        HD1  HE1
                  H-N         |    |
                    |   HB1  CD1--CE1
                    |   |    //     \\
                 HA-CA--CB--CG      CZ--HZ
                    |   |    \  __  /
                    |   HB2  CD2--CE2
                  O=C         |    |
                    |        HD2  HE2
        ",
        );
        m.insert(
            "PRO",
            r"
                    | HD1 HD2
                    |   \ /
                    N---CD   HG1
                    |     \  /
                    |      CG
                    |     /  \
                 HA-CA--CB   HG2
                    |   / \
                    | HB1 HB2
                  O=C
                    |
        ",
        );
        m.insert(
            "SER",
            r"
                    |
                  H-N
                    |   HB1
                    |   |
                 HA-CA--CB--OG
                    |   |     \
                    |   HB2    HG1
                  O=C
                    |
        ",
        );
        m.insert(
            "THR",
            r"
                    |
                  H-N
                    |     OG1--HG1
                    |    /
                 HA-CA--CB-HB
                    |    \
                    |     CG2--HG21
                  O=C    / \
                    | HG21 HG22
        ",
        );
        m.insert(
            "TRP",
            r"
                    |                  HE3
                  H-N                   |
                    |   HB1            CE3
                    |   |             /  \\
                 HA-CA--CB---CG-----CD2   CZ3-HZ3
                    |   |    ||     ||     |
                    |   HB2  CD1    CE2   CH2-HH2
                  O=C       /   \   / \  //
                    |     HD1    NE1   CZ2
                                  |     |
                                 HE1   HZ2
        ",
        );
        m.insert(
            "TYR",
            r"
                    |        HD1  HE1
                  H-N         |    |
                    |   HB1  CD1--CE1
                    |   |   //      \\
                 HA-CA--CB--CG      CZ--OH
                    |   |    \  __  /     \
                    |   HB2  CD2--CE2     HH
                  O=C         |    |
                    |        HD2  HE2
        ",
        );
        m.insert(
            "VAL",
            r"
                    |    HG11 HG12
                  H-N      | /
                    |     CG1--HG13
                    |    /
                 HA-CA--CB-HB
                    |    \
                    |     CG2--HG21
                  O=C    / \
                    | HG21 HG22
        ",
        );
        m
    };
}
