#! /bin/bash
NPROC=8

RS=8
MINGSAFE=2000
F1ORDER=3
F0ORDER=2

NI=1e20
TMAX=1000
TMIN=100
TGRAD=180
XPOINT=0.046775 # in cm

#ZBAR=1
#ZBAR=2
#ZBAR=5
#ZBAR=10
#ZBAR=20
#ZBAR=50

MINGSAFE=1000
MINGIMPL=35
MINGIMPLSAFE=1000

L=0.1
PROBLEM=5

# The highest Z.
ZBAR=100

if false; then
SIGMA=1.2e10 # Zbar = 100 -> Kn 1e-4
## Nonlocal solution very well corresponding to Pascal's solution with Aladin.
## P1 closure.
## C7*
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMA -Tgrad $TGRAD -Z $ZBAR -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e5 -S0 1.0 -E0 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMA -n $NI -xp $XPOINT --Emimic --labelEmimic C7implicit
cd ../..
fi

SIGMASAFE=1.2e15 # Zbar = 100 -> Kn 1e-9
if false ; then
## Nonlocal solution very well corresponding to Pascal's solution with Aladin.
## P1 closure.
## C7*
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMASAFE -Tgrad $TGRAD -Z $ZBAR -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e5 -S0 1.0 -E0 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMASAFE -Tgrad $TGRAD -Z $ZBAR -ni $NI -L $L -xp $XPOINT -minG $MINGSAFE -s 2 -cfl 0.25 -S0 1.0 -E0 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Ecorrect_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Ecorrect_data/fe_point_Ecorrect.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMASAFE -n $NI -xp $XPOINT --Emimic --Ecorrect --labelEmimic C7implicit --labelEcorrect C7explicit
cd ../..

SIGMA=1.2e10 # Zbar = 100 -> Kn 1e-4
## Nonlocal solution very well corresponding to Pascal's solution with Aladin.
## P1 closure.
## C7*
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMA -Tgrad $TGRAD -Z $ZBAR -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e5 -S0 1.0 -E0 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMA -n $NI -xp $XPOINT --Emimic --labelEmimic C7implicit
cd ../..

SIGMA=1.2e9 # Zbar = 100 -> Kn 1e-3
#if false ; then
## Nonlocal solution very well corresponding to Pascal's solution with Aladin.
## P1 closure.
## C7*
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMA -Tgrad $TGRAD -Z $ZBAR -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e5 -S0 1.0 -E0 0.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMA -n $NI -xp $XPOINT --Emimic --labelEmimic C7implicit
cd ../..
fi


# The lowest Z.
ZBAR=1

SIGMASAFE=6.1e18 # Zbar = 1 -> Kn 1e-9
if false ; then
## Nonlocal solution very well corresponding to Pascal's solution with Aladin.
## P1 closure.
## C7*
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMASAFE -Tgrad $TGRAD -Z $ZBAR -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e5 -S0 1.0 -E0 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMASAFE -Tgrad $TGRAD -Z $ZBAR -ni $NI -L $L -xp $XPOINT -minG $MINGSAFE -s 2 -cfl 0.25 -S0 1.0 -E0 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Ecorrect_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Ecorrect_data/fe_point_Ecorrect.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMASAFE -n $NI -xp $XPOINT --Emimic --Ecorrect --labelEmimic C7implicit --labelEcorrect C7explicit
cd ../..

SIGMA=6.1e13 # Zbar = 1 -> Kn 1e-4
## Nonlocal solution very well corresponding to Pascal's solution with Aladin.
## P1 closure.
## C7*
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMA -Tgrad $TGRAD -Z $ZBAR -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e5 -S0 1.0 -E0 1.47
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMA -n $NI -xp $XPOINT --Emimic --labelEmimic C7implicit
cd ../..

SIGMA=6.1e12 # Zbar = 1 -> Kn 1e-3
#if false ; then
## Nonlocal solution very well corresponding to Pascal's solution with Aladin.
## P1 closure.
## C7*
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -sigma $SIGMA -Tgrad $TGRAD -Z $ZBAR -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e5 -S0 1.0 -E0 1.45
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMA -n $NI -xp $XPOINT --Emimic --labelEmimic C7implicit
cd ../..
fi



## Diffusive asymptotic test. The SH, AWBS (original and corrected), 
## and C7 (proper and mimic Efield) calculations are compared. 
## In reality, the resulting Knudsen number is just on the diffusive limit 
## in the following case.
#mpirun -np 8 C7 -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 2 -ot 1 -no-vis -fa -print -Tmax 800.5 -Tmin 799.5 -Tgrad 2.3 -S0 1.0 -E0 1.0 -sigma 1e11 -Z 4 -ne 5e20 -M1 -minG 2000
## A pure diffusion case. 
## Converged numerical flux from -minG 200. Err 1e-5 -minG 50.
RS=7
F1ORDER=3
F0ORDER=2
NI=5e20
ZBAR=100
MING=25
if [ $ZBAR -gt 10 ] ; then 
MING=2000 
fi
TMAX=800.5
TMIN=799.5
#TGRAD=2.3
#TGRAD=80
TGRAD=300
SIGMA=8.7e18 ## Kn=1.0e-4 nonlocality limit.
#SIGMA=8.7e9 ## Kn=1.0e-4 nonlocality limit.
#SIGMA=8.7e8 ## Kn=1.0e-3 q=0.8*qSH E0=0.91
#SIGMA=1.75e8 ## Zbar = 100, Kn=5.0e-3 q=0.01*qSH

#SIGMA=8.45e15    ## Kn 1.0e-10
#SIGMA=8.45e11    ## Kn 1.0e-6
#SIGMA=1.69e11    ## Kn 5.0-6
#SIGMA=0.845e11   ## Kn 1.0-5
#SIGMA=0.535e11   ## Kn 1.6-5 qC7E zero
#SIGMA=0.845e10    ## Kn 1.0-4 qC7E negative
#ZBAR=100
#MING=$MINGSAFE
L=0.1
XPOINT=0.05
PROBLEM=5

#if false ; then
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -sigma $SIGMA -Z $ZBAR -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e5 -S0 1.0 -E0 0.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -sigma $SIGMA -Z $ZBAR -ni $NI -L $L -M1 -xp $XPOINT -minG $MING -s 2 -cfl 0.25 -S0 1.0 -E0 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Ecorrect_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Ecorrect_data/fe_point_Ecorrect.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -s $SIGMA -n $NI -xp $XPOINT --Ecorrect --Emimic --AWBSoriginal
#fi
