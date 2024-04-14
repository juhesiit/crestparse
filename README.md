# Crestparse 0.2.1

Crestparse is a command line tool for analyzing and manipulating multi-xyz files produced by crest.

## Usage

Analyzing a conformation file is carried out by simply running

```
crestparse crest_conformers.xyz
```

which will output a default

```
#       E (Hartree)     dE (Hartree)    dE (kcal/mol)   Boltzmann T=298.15 K
 1      -21.044884      0.000000        0.000000        72.390641%
 2      -21.043020      0.001864        1.169905        10.049309%
 3      -21.042881      0.002003        1.257052        8.674748%
 4      -21.042213      0.002671        1.675864        4.278241%
 5      -21.041661      0.003223        2.022288        2.384156%
 6      -21.040953      0.003931        2.466671        1.126159%
 7      -21.040082      0.004802        3.013393        0.447558%
 8      -21.039917      0.004967        3.117012        0.375748%
 9      -21.039617      0.005267        3.305325        0.273440%
```

Conformations can be extracted to separate files (`conf1.xyz`) with `-e` or `--extract` followed by the desired conformer index (or multiple indeces). By providing just `--extract`, all conformers within the energy window will be exported as `.xyz` files. If a cutoff limit is invoked, only conformers until the energy cutoff limit will be extracted.

```
crestparse crest_conformers.xyz -e 1,3,4
```

Boltzmann temperature can be changed with `-t` or `--temperature` followed by the desired temperature in Kelvin.

```
crestparse crest_conformers.xyz -t 100
```

High energy conformers can be excluded with `-c` or `--cutoff` followed by the cutoff energy in kcal/mol. For example, to list all conformers with relative energies lower than 2.0 kcal/mol.

```
crestparse crest_conformers.xyz -c 2.0
```

Atomic distances can be tabulated with `-d` or `--distance` followed by two atom indeces (indexing starts from 0). For example, to list distance between atoms 4 and 7 in all conformers.

```
crestparse crest_conformers.xyz -d 4 7
```

will output

```
#       E (Hartree)     dE (Hartree)    dE (kcal/mol)   Boltzmann T=298.15 K    Distance 4-7
 1      -21.044884      0.000000        0.000000        72.390641%              2.961969
 2      -21.043020      0.001864        1.169905        10.049309%              3.646953
 3      -21.042881      0.002003        1.257052        8.674748%               3.293627
 4      -21.042213      0.002671        1.675864        4.278241%               2.948717
 5      -21.041661      0.003223        2.022288        2.384156%               3.737814
 6      -21.040953      0.003931        2.466671        1.126159%               3.849287
 7      -21.040082      0.004802        3.013393        0.447558%               4.164731
 8      -21.039917      0.004967        3.117012        0.375748%               3.706732
 9      -21.039617      0.005267        3.305325        0.273440%               3.328362
```

Angles between three atoms can be tabulated with `-a` or `--angle` followed by three atom indeces (indexing starts from 0). For example, to list angles between atoms 0, 1 and 2 in all conformers.

```
crestparse crest_conformers.xyz -a 0 1 2
```

will output

```
#       E (Hartree)     dE (Hartree)    dE (kcal/mol)   Boltzmann T=298.15 K    Angle 0-1-2
 1      -21.044884      0.000000        0.000000        72.390641%              112.652563
 2      -21.043020      0.001864        1.169905        10.049309%              115.229758
 3      -21.042881      0.002003        1.257052        8.674748%               108.829629
 4      -21.042213      0.002671        1.675864        4.278241%               111.963991
 5      -21.041661      0.003223        2.022288        2.384156%               110.297835
 6      -21.040953      0.003931        2.466671        1.126159%               107.924930
 7      -21.040082      0.004802        3.013393        0.447558%               107.876325
 8      -21.039917      0.004967        3.117012        0.375748%               110.353769
 9      -21.039617      0.005267        3.305325        0.273440%               108.053664
```

Dihedrals between four atoms can be tabulated with `-r` or `--dihedral` followed by four atom indeces (indexing starts from 0). For example, to list dihedrals between atoms 0, 1, 2 ja 3 in all conformers.

```
crestparse crest_conformers.xyz -a 0 1 2
```

will output

```
#       E (Hartree)     dE (Hartree)    dE (kcal/mol)   Boltzmann T=298.15 K    Dihedral 0-1-2-3
 1      -21.044884      0.000000        0.000000        72.390641%              -12.457067
 2      -21.043020      0.001864        1.169905        10.049309%              167.512960
 3      -21.042881      0.002003        1.257052        8.674748%               -33.705215
 6      -21.040953      0.003931        2.466671        1.126159%               -163.834904
 7      -21.040082      0.004802        3.013393        0.447558%               161.605868
 8      -21.039917      0.004967        3.117012        0.375748%               -144.963430
 9      -21.039617      0.005267        3.305325        0.273440%               175.822122
```
