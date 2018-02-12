#! /bin/bash
## Philippes tests, on the edge of locality.
## Explicit Efield.
#mpirun -np 8 C7 -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -S0 1.0 -E0 1.0 -a0 5e31
## Mimiced Efield.
#mpirun -np 8 C7 -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 5e31

### Pascal's setting for nonlocal test
## First, diffusive, case sets sigma 1e5x higher, which assures SH solution.
mpirun -np 8 C7 -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -sigma 1e17 -Tgrad 180 -Z 2 -ne 1e20 -L 0.1 -minG 2000
## Nonlocal solution very well corresponding to Pascal's solution with Aladin.
## P1 closure.
mpirun -np 8 C7 -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -sigma 1e12 -Tgrad 180 -Z 2 -ne 1e20 -L 0.1 -minG 2000
## M1 closure.
mpirun -np 8 C7 -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -sigma 1e12 -Tgrad 180 -Z 2 -ne 1e20 -L 0.1 -M1 -minG 200
 
## Diffusive asymptotic test. The SH, AWBS (original and corrected), 
## and C7 (proper and mimic Efield) calculations are compared. 
## In reality, the resulting Knudsen number is just on the diffusive limit 
## in the following case.
#mpirun -np 8 C7 -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 2 -ot 1 -no-vis -fa -print -Tmax 800.5 -Tmin 799.5 -Tgrad 2.3 -S0 1.0 -E0 1.0 -sigma 1e11 -Z 4 -ne 5e20 -M1 -minG 2000
## A pure diffusion case. 
## Converged numerical flux from -minG 200. Err 1e-5 -minG 50.
mpirun -np 8 C7 -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 2 -ot 1 -no-vis -fa -print -Tmax 800.5 -Tmin 799.5 -Tgrad 2.3 -S0 1.0 -E0 1.0 -sigma 1e15 -Z 4 -ne 5e20 -M1 -minG 20
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/fe_point_Ecorrect.txt

mpirun -np 8 C7 -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 2 -ot 1 -no-vis -fa -print -Tmax 800.5 -Tmin 799.5 -Tgrad 2.3 -sigma 1e15 -Z 4 -ne 5e20 -minG 15
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/fe_point_Emimic.txt

cd results/fe_analysis
python AWBSf1_integrate.py -Z 4 -T 8e2 -g -2.3 -s 1e15 -n 5e20
