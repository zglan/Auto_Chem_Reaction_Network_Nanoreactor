trjSplit.py     # (or python2 trjSplit.py)

Input:
```xtb.trj```

Output:
```
sampling      # Directory
└─[number]     # Sub-directory of sampling, eg: 1, 2, ...
    └─stru.xyz    # File of Sub-directory
```

Note: 
    1. The script trjSplit.py needs to be run via python2.
    2. If there are too much structures in xtb.trj, try to change the interval of trjSplit.py, or directly use scoord file by xtb 