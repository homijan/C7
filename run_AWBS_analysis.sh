#! /bin/bash
## Philippes tests, on the edge of locality.
## Explicit Efield.
#mpirun -np 8 C7 -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -S0 1.0 -E0 1.0 -a0 5e31
## Mimiced Efield.
#mpirun -np 8 C7 -p 5 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 4 -ot 3 -vis -fa -print -Tmax 1000 -Tmin 100 -a0 5e31

NPROC=8
## Working Efield general diffusive regime.
#RS=7 # Increasing spatial resolution improves the stability the best!
#MINGSAFE=2000
#F1ORDER=3
#F0ORDER=2
#E0=1.0
#SIGMA=1e15 # Kn 1e-5, safe for Z=2
##SIGMA=1e16 # Kn 1e-6, safe for Z=4
#SIGMASAFE=1e18 # Kn 1e-8, safe for any Z

RS=7
MINGSAFE=2000
F1ORDER=3
F0ORDER=2
E0=1.0
#SIGMA=1.05e11 # Kn 1e-1
#SIGMA=0.53e12 # Kn 2e-2
SIGMA=1.06e12 # Kn 1e-2
SIGMASAFE=1.06e18 # Kn 1e-8
XPOINT=0.046775 # in cm

NE=1e20
TMAX=1000
TMIN=100
TGRAD=180
ZBAR=4
L=0.1
PROBLEM=5
### Pascal's setting for nonlocal test
## First, diffusive, case sets sigma 1e5x higher, which assures SH solution.
## C7E
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMASAFE -Tgrad $TGRAD -Z $ZBAR -ne $NE -L $L -xp $XPOINT -minG $MINGSAFE -S0 1.0 -E0 $E0
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/fe_point_Ecorrect.txt
## C7*
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMASAFE -Tgrad $TGRAD -Z $ZBAR -ne $NE -L $L -xp $XPOINT -minG $MINGSAFE
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/fe_point_Emimic.txt
cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMASAFE -n $NE -xp $XPOINT --Emimic --Ecorresults
cd ../..
"""
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMA -Tgrad $TGRAD -Z $ZBAR -ne $NE -L $L -xp $XPOINT -minG $MINGSAFE -S0 1.0 -E0 $E0
cd results/fe_analysis
cp ../tmp/C7_1_fe_point.txt fe_point_Emimic.txt
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMA -n $NE -xp $XPOINT --Ecorresults
cd ../..
"""
## Nonlocal solution very well corresponding to Pascal's solution with Aladin.
## P1 closure.
## C7*
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMA -Tgrad $TGRAD -Z $ZBAR -ne $NE -L $L -xp $XPOINT -minG $MINGSAFE
cd results/fe_analysis
cp ../tmp/C7_1_fe_point.txt fe_point_Emimic.txt
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMA -n $NE -xp $XPOINT --Emimic
cd ../..
## M1 closure.
## C7*
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMA -Tgrad $TGRAD -Z $ZBAR -ne $NE -L $L -xp $XPOINT -minG $MINGSAFE -M1
cd results/fe_analysis
cp ../tmp/C7_1_fe_point.txt fe_point_Emimic.txt
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMA -n $NE -xp $XPOINT --Emimic
cd ../..
"""
## M1 closure.
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMA -Tgrad $TGRAD -Z $ZBAR -ne $NE -L $L -xp $XPOINT -minG $MINGSAFE -M1
cd results/fe_analysis
cp ../tmp/C7_1_fe_point.txt fe_point_Emimic.txt
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMA -n $NE -xp $XPOINT --Emimic
cd ../..
"""
 
## Diffusive asymptotic test. The SH, AWBS (original and corrected), 
## and C7 (proper and mimic Efield) calculations are compared. 
## In reality, the resulting Knudsen number is just on the diffusive limit 
## in the following case.
#mpirun -np 8 C7 -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 2 -ot 1 -no-vis -fa -print -Tmax 800.5 -Tmin 799.5 -Tgrad 2.3 -S0 1.0 -E0 1.0 -sigma 1e11 -Z 4 -ne 5e20 -M1 -minG 2000
## A pure diffusion case. 
## Converged numerical flux from -minG 200. Err 1e-5 -minG 50.
NE=5e20
ZBAR=4
MING=30
TMAX=800.5
TMIN=799.5
TGRAD=2.3
SIGMA=8.45e15    ## Kn 1.0e-10
#SIGMA=8.45e11    ## Kn 1.0e-6
#SIGMA=1.69e11    ## Kn 5.0-6
#SIGMA=0.845e11   ## Kn 1.0-5
#SIGMA=0.535e11   ## Kn 1.6-5 qC7E zero
#SIGMA=0.845e10    ## Kn 1.0-4 qC7E negative
#ZBAR=100
#MING=$MINGSAFE
L=0.1
XPOINT=0.05
PROBLEM=8

mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -S0 1.0 -E0 1.0 -sigma $SIGMA -Z $ZBAR -ne $NE -L $L -M1 -xp $XPOINT -minG $MING
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/fe_point_Ecorrect.txt

mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -sigma $SIGMA -Z $ZBAR -ne $NE -L $L -xp $XPOINT -minG $MING
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/fe_point_Emimic.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMA -n $NE -xp $XPOINT --Ecorresults --Emimic
