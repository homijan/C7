#! /bin/bash
SIGMA=8.1027575e17 ## Matching the SH diffusive flux.
CL=10.0 # Coulomb logarithm.

NPROC=8

RS=8
MINGSAFE=2000
F1ORDER=3
F0ORDER=2

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
ZBAR=100

MINGSAFE=1000
MINGIMPL=35

L=0.1
PROBLEM=5

#if false; then
NI=2e25 # Zbar = 100 -> Kn 1e-10
if [ $ZBAR -le 10 ] ; then
   NI=1e29 # Zbar = 1 -> Kn 1e-10
fi
## Nonlocal solution very well corresponding to Pascal's solution with Aladin.
## P1 closure.
## C7*
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -Z $ZBAR -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e10 -S0 1.0 -E0 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -Z $ZBAR -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MINGSAFE -s 2 -cfl 0.25 -S0 1.0 -E0 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Ecorrect_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Ecorrect_data/fe_point_Ecorrect.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -cl $CL -n $NI -xp $XPOINT --Emimic --Ecorrect --labelEmimic C7impl --labelEcorrect C7expl --AWBSoriginal
cd ../..
#fi

#if false; then
NI=2e18 # Zbar = 100 -> Kn 1e-3
if [ $ZBAR -le 10 ] ; then
   NI=1e21 # Zbar = 1 -> Kn 1e-2
fi
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -Z $ZBAR -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e10
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -cl $CL -n $NI -xp $XPOINT --Emimic --labelEmimic 'C7impl*' --AWBSoriginal
cd ../..
#fi

#if false; then
## Highest Kn limit to compute.
NI=2e20 # Zbar = 100 -> Kn 1e-5
if [ $ZBAR -le 10 ] ; then
   NI=1e22 # Zbar = 1 -> Kn 1e-3
fi
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -Z $ZBAR -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e10 -S0 1.0 -E0 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -cl $CL -n $NI -xp $XPOINT --Emimic --labelEmimic 'C7impl' --AWBSoriginal
cd ../..
#fi

## Diffusive asymptotic test. The SH, AWBS (original and corrected), 
## and C7 (proper and mimic Efield) calculations are compared. 
## In reality, the resulting Knudsen number is just on the diffusive limit 
## in the following case.
#mpirun -np 8 C7 -p 8 -m data/segment01.mesh -rs 6 -tf 0.0 -ok 2 -ot 1 -no-vis -fa -print -Tmax 800.5 -Tmin 799.5 -Tgrad 2.3 -S0 1.0 -E0 1.0 -sigma 1e11 -Z 4 -ne 5e20 -M1 -minG 2000
## A pure diffusion case. 
## Converged numerical flux from -minG 200. Err 1e-5 -minG 50.

RS=8
F1ORDER=3
F0ORDER=2
MING=25
if [ $ZBAR -gt 10 ] ; then 
MING=2000 
fi
TMAX=800.5
TMIN=799.5
TGRAD=300

L=0.1
XPOINT=0.05
PROBLEM=5

#if false ; then
#NI=7.044e25 # Zbar = 100, Kn = 1e-10
#NI=7.044e20 # Zbar = 100, Kn = 1e-5
NI=7.044e18 # Zbar = 100, Kn = 1e-3
if [ $ZBAR -le 10 ] ; then
   #NI=3.522e22 # Zbar = 1, Kn = 1e-3
   NI=3.522e21 # Zbar = 1, Kn = 1e-2
fi
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -Z $ZBAR -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MINGIMPL -s 4 -cfl 1e10 -S0 1.0 -E0 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -cl $CL -n $NI -xp $XPOINT --Emimic --labelEmimic 'C7impl' --AWBSoriginal
#fi
