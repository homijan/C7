#! /bin/bash
## Philippes tests, on the edge of locality.
## Explicit Efield.
#mpirun -np 8 C7 -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -S0 1.0 -E0 1.0 -a0 5e31
## Mimiced Efield.
#mpirun -np 8 C7 -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 5e31

MINGSAFE=500
NPROC=8

SIGMA=1e12
SIGMASAFE=1e17
NE=1e20
TMAX=1000
TMIN=100
TGRAD=180
ZBAR=2
L=0.1
XPOINT=0.046 # in cm
PROBLEM=5
### Pascal's setting for nonlocal test
## First, diffusive, case sets sigma 1e5x higher, which assures SH solution.
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMASAFE -Tgrad $TGRAD -Z $ZBAR -ne $NE -L $L -xp $XPOINT -minG $MINGSAFE
cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -T 8e2 -g -2.3 -s $SIGMA -n $NE -xp $XPOINT
cd ../..
## Nonlocal solution very well corresponding to Pascal's solution with Aladin.
## P1 closure.
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMA -Tgrad $TGRAD -Z $ZBAR -ne $NE -L $L -xp $XPOINT -minG $MINGSAFE
cd results/fe_analysis 
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -T 8e2 -g -2.3 -s $SIGMA -n $NE -xp $XPOINT
cd ../..
## M1 closure.
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMA -Tgrad $TGRAD -Z $ZBAR -ne $NE -L $L -xp $XPOINT -minG $MINGSAFE -M1
cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -T 8e2 -g -2.3 -s $SIGMA -n $NE -xp $XPOINT
cd ../..
 
## Diffusive asymptotic test. The SH, AWBS (original and corrected), 
## and C7 (proper and mimic Efield) calculations are compared. 
## In reality, the resulting Knudsen number is just on the diffusive limit 
## in the following case.
#mpirun -np 8 C7 -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 2 -ot 1 -no-vis -fa -print -Tmax 800.5 -Tmin 799.5 -Tgrad 2.3 -S0 1.0 -E0 1.0 -sigma 1e11 -Z 4 -ne 5e20 -M1 -minG 2000
## A pure diffusion case. 
## Converged numerical flux from -minG 200. Err 1e-5 -minG 50.
NE=5e20
ZBAR=4
MING=20
SIGMA=8.45e14    ## Kn 1.0e-10
#SIGMA=8.45e10    ## Kn 1.0e-6
#SIGMA=1.69e10    ## Kn 5.0-6
#SIGMA=0.845e10   ## Kn 1.0-5
#SIGMA=0.535e10   ## Kn 1.6-5 qC7E zero
#SIGMA=0.845e9    ## Kn 1.0-4 qC7E negative
#ZBAR=100
#MING=$MINGSAFE
XPOINT=0.5

mpirun -np $NPROC C7 -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 2 -ot 1 -no-vis -fa -print -Tmax 800.5 -Tmin 799.5 -Tgrad 2.3 -S0 1.0 -E0 1.0 -sigma $SIGMA -Z $ZBAR -ne $NE -M1 -xp $XPOINT -minG $MING
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/fe_point_Ecorrect.txt

mpirun -np $NPROC C7 -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 2 -ot 1 -no-vis -fa -print -Tmax 800.5 -Tmin 799.5 -Tgrad 2.3 -sigma $SIGMA -Z $ZBAR -ne $NE -xp $XPOINT -minG $MING
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/fe_point_Emimic.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -T 8e2 -g -2.3 -s $SIGMA -n $NE -xp $XPOINT
