# Crestparse 0.1

Crestparse is a command line tool for analyzing and manipulating multi-xyz files produced by crest.

## Usage

Analyzing a conformation file is carried out by simply doing

```
crestparse crest_conformers.xyz
```

which will output a default

```
Successfully read in 57 structures
Calculating the Boltzmann distribution at 298.15 K for total of 57 conformers...
#       E (Hartree)     dE (Hartree)    dE (kcal/mol)   Boltzmann (%) T = 298.15 K
 1      -33.880227      0.000000        0.000000        17.483862%
 2      -33.880161      0.000066        0.041333        16.305712%
 3      -33.879471      0.000756        0.474277        7.852173%
 4      -33.879375      0.000853        0.535044        7.086752%
 5      -33.878869      0.001359        0.852490        4.147226%
 6      -33.878752      0.001476        0.926027        3.663152%
 7      -33.878735      0.001493        0.936663        3.597979%
 8      -33.878699      0.001528        0.958946        3.465177%
 9      -33.878604      0.001624        1.018903        3.131671%
10      -33.878602      0.001625        1.019656        3.127693%
11      -33.878544      0.001683        1.055982        2.941691%
12      -33.878506      0.001722        1.080404        2.822899%
13      -33.878437      0.001791        1.123620        2.624326%
14      -33.878429      0.001799        1.128659        2.602102%
15      -33.878413      0.001814        1.138592        2.558840%
...
```

Conformations can be extracted to spearate files (`conf1.xyz`) with `-e` or `--extract` followed by the desired conformer index (or multiple indeces)

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
