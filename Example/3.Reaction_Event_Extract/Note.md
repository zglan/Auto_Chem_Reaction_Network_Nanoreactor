dynReacExtr.py is needed.

Input:
sampling      # Directory
├─control.sh(Or dynReacExtr.py)
└─[number]     # Sub-directory of sampling, eg: 1, 2, ...
    └─xtb.trj    # trajectory file of mtd.

Usage:
```
dynReacExtr.py -i trj --xtbopt --ts --refine
```

Output:
eg: 
```
*  1 times  |  [H]C1C([H])C(N([H])C(O)N([H])C2C([H])C(Cl)C(Cl)C([H])C2[H])C([H])C([H])C1Cl+[H]OCl->Cl+[H]C1C([H])C(NC(O)N([H])C2C([H])C(Cl)C(Cl)C([H])C2[H])C([H])C([H])C1Cl+[H]O[H]
*  2 times  |  [H]OCl->Cl+[H]O
*  1 times  |  [H]C1C([H])C(NC(O)N([H])C2C([H])C([H])C(Cl)C(Cl)C2[H])C([H])C([H])C1Cl+[H]O->[H]C1C([H])C(NC(O)NC2C([H])C(Cl)C(Cl)C([H])C2[H])C([H])C([H])C1Cl+[H]O[H]
*  1 times  |  [H]OOC(NC1C([H])C([H])C(Cl)C([H])C1[H])NC1C([H])C(Cl)C(Cl)C([H])C1[H]+[H]O[H]->[H]C1C([H])C(NC(O)NC2C([H])C(Cl)C(Cl)C([H])C2[H])C([H])C([H])C1Cl+[H]O+[H]O[H]
*  1 times  |  [H]O+[H]OCl->OCl+[H]O[H]
*  1 times  |  Cl+[H]C1C([H])C(NC(O)NC2C([H])C(Cl)C(Cl)C([H])C2[H])C([H])C([H])C1Cl->[H]C1C([H])C(NC(O)NC2C([H])C([H])C(Cl)C(Cl)C2([H])Cl)C([H])C([H])C1Cl
*  1 times  |  [H]C1C([H])C(NC(O)NC2C([H])C([H])C(Cl)C(Cl)C2([H])Cl)C([H])C([H])C1Cl->Cl+[H]C1C([H])C(NC(O)NC2C([H])C(Cl)C(Cl)C([H])C2[H])C([H])C([H])C1Cl
*  1 times  |  [H]C1C([H])C(NC(O)NC2C([H])C(Cl)C(Cl)C([H])C2[H])C([H])C([H])C1Cl->[H]C1C([H])C(NC2NC3C([H])C([H])C(Cl)C([H])C3([H])O2)C([H])C(Cl)C1Cl
*  1 times  |  Cl+[H]C1C(Cl)C(Cl)C([H])C(NC2NC3C([H])C([H])C(Cl)C([H])C3([H])O2)C1[H]->[H]C1C([H])C2C(OC(NC3C([H])C([H])C(Cl)C(Cl)C3[H])N2)C([H])C1Cl+[H]Cl
*  1 times  |  OCl+[H]Cl->Cl+[H]OCl
*  1 times  |  [H]OCl+[H]O[H]->Cl+[H]O+[H]O[H]
*  1 times  |  [H]OC12OC(NC3C([H])C(Cl)C(Cl)C([H])C3[H])NC1C([H])(Cl)C([H])C(Cl)C2[H]->Cl+[H]OC12OC(NC3C([H])C(Cl)C(Cl)C([H])C3[H])NC1C([H])C([H])C(Cl)C2[H]
*  2 times  |  Cl+Cl->ClCl
```

[number]     # Sub-directory of sampling, eg: 1, 2, ...
├─sp            # sp directory, for single point calculation. The command "xtb [id].xyz --input [id].inp" have been used.
|  └─[number]   # Sub-directory of sp
|    ├─[id].charges     # charges data
|    ├─[id].xyz         # xyz file
|    ├─[id].inp         # job file
|    ├─[id].log         
|    └─[id].wbo
├─ts            # ts directory, for ts optimization by gaussian.
|  └─[number]   # Sub-directory of ts
|    └─[id].gjf         # TS initial guess structure file.
├─xtbopt        # xtbopt directory, for optimization by xtb. The command "xtb [id].xyz --input [id].inp --opt" have been used.
|  └─[number]   # Sub-directory of xtbopt
|    ├─[id]-opt.xyz     # xyz file after optimization.
|    ├─[id].inp         # job file.
|    ├─[id].log         
|    └─[id].xyz         # xyz file.
├─_job.csv
├─_opt_xtb_data.csv
├─_ref_relation.csv
├─_relation.csv
├─.ini_data.npy
└─tsjob.info


```
# Initial Analysis
.ini_data.npy
_relation.csv
_job.csv

# --refine: Reactant, Products(In SMILES).
_ref_relation.csv

# --xtbopt: Optimized Structure in xtbopt file.
xtbopt
_opt_xtb_data.csv

# --ts: Folder ts has initial guess structure for TS.
tsjob.info
ts
```