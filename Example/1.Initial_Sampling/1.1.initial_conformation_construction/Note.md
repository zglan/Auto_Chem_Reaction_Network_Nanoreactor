# Construct initial xyz file by sphereMaker.py
```sphereMaker.py -i H2O.xyz HClO.xyz TCC2.xyz -n 10 4 1 -d 1.5 --dis 3 --center TCC2.xyz --name TCC_Cl4```

```
*** Sphere Radii is 6.90 A.
*** Sphere Radii is 13.04 Bohr.
=== TCC_Cl4.xyz is built
```

# Optimize md_init.xyz file by xtb if needed.(Optional)
```xtb md_init.xyz --opt```
```mv xtbopt.xyz md_init.xyz```