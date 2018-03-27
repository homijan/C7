#! /bin/bash
SIGMA=8.1027575e17 ## Matching the SH diffusive flux.
CL=10.0 # Coulomb logarithm.

NPROC=8

RS=6
F1ORDER=4
F0ORDER=3

TMAX=1000
TMIN=100
TGRAD=180
XPOINT=0.051722 # in cm
#XPOINT=0.046775 # in cm

ZBAR=1
#ZBAR=5
#ZBAR=10
#ZBAR=20
#ZBAR=100

#MINGSAFE=1000
MING=100
#MING=2000

# Challenge SNB ;)
#MING=25
#F1ORDER=2
#F0ORDER=1

L=0.1
PROBLEM=5
#PROBLEM=9

#if false; then
## Highest Kn limit to compute.
if [ $ZBAR -eq 100 ] ; then 
   #NI=2e25 # Zbar = 100 -> Kn 1e-10
   #NI=2e20 # Zbar = 100 -> Kn 1e-5
   #NI=2e19 # Zbar = 100 -> Kn 1e-4
   #NI=4e18 # Zbar = 100 -> Kn 5e-4
   #NI=2e18 # Zbar = 100 -> Kn 1e-3
   NI=2e17 # Zbar = 100 -> Kn 1e-2
   #NI=2e16 # Zbar = 100 -> Kn 1e-1
elif [ $ZBAR -eq 20 ] ; then 
   NI=4.8e21 # Zbar = 20 -> Kn 1e-5
   #NI=4.8e20 # Zbar = 20 -> Kn 1e-4
   #NI=4.8e19 # Zbar = 20 -> Kn 1e-3 
   #NI=4.8e18 # Zbar = 20 -> Kn 1e-2
   #NI=0.963e18 # Zbar = 20 -> Kn 5e-2
elif [ $ZBAR -eq 10 ] ; then 
   NI=1.84e22 # Zbar = 10 -> Kn 1e-5
   #NI=1.84e21 # Zbar = 10 -> Kn 1e-4
   #NI=1.84e20 # Zbar = 10 -> Kn 1e-3
   #NI=0.368e20 # Zbar = 10 -> Kn 5e-3
   #NI=1.84e19 # Zbar = 10 -> Kn 1e-2
elif [ $ZBAR -eq 5 ] ; then 
   NI=6.74e22 # Zbar = 5 -> Kn 1e-5
   #NI=6.74e21 # Zbar = 5 -> Kn 1e-4
   #NI=6.74e20 # Zbar = 5 -> Kn 1e-3
   #NI=6.74e19 # Zbar = 5 -> Kn 1e-2
   #NI=1.348e19 # Zbar = 5 -> Kn 5e-2
elif [ $ZBAR -eq 1 ] ; then
   #NI=1e29 # Zbar = 1 -> Kn 1e-10
   #NI=1e24 # Zbar = 1 -> Kn 1e-5
   #NI=1e23 # Zbar = 1 -> Kn 1e-4
   #NI=1e22 # Zbar = 1 -> Kn 1e-3
   #NI=1e21 # Zbar = 1 -> Kn 1e-2
   #NI=2.02e20 # Zbar = 1 -> Kn 5e-2
   NI=1e20 # Zbar = 1 -> Kn 1e-1, dE 0.1, MING 35
fi
#mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -Z $ZBAR -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -dE 0.001  | tee C7E.out
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -Z $ZBAR -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -dE 0.05  | tee C7E.out
cp results/tmp/C7_1_profiles.* results/fe_analysis/Ecorrect_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Ecorrect_data/fe_point_Ecorrect.txt
cp results/tmp/C7_1_fe_pointmax.txt results/fe_analysis/Ecorrect_data/fe_pointmax_Ecorrect.txt
# Store the output file.
#cp C7E.out results/fe_analysis/C7E/P5_Z1_Kn1e-5.out
#cp C7E.out results/fe_analysis/C7E/P5_Z1_Kn1e-4.out
#cp C7E.out results/fe_analysis/C7E/P5_Z1_Kn1e-3.out
#cp C7E.out results/fe_analysis/C7E/P5_Z1_Kn1e-2.out

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -cl $CL -n $NI --Ecorrect --labelEcorrect 'C7E' #--AWBSoriginal
# Safe figs.
#cp heatflux.png C7E/P5_heatflux_Z1_Kn1e-5.png
#cp kinetics.png C7E/P5_kinetics_Z1_Kn1e-5.png
#cp heatflux.png C7E/P5_heatflux_Z1_Kn1e-4.png
#cp kinetics.png C7E/P5_kinetics_Z1_Kn1e-4.png
#cp heatflux.png C7E/P5_heatflux_Z1_Kn1e-3.png
#cp kinetics.png C7E/P5_kinetics_Z1_Kn1e-3.png
#cp heatflux.png C7E/P5_heatflux_Z1_Kn1e-2.png
#cp kinetics.png C7E/P5_kinetics_Z1_Kn1e-2.png
cd ../..
#fi

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

if false ; then
#NI=7.044e25 # Zbar = 100, Kn = 1e-10
#NI=7.044e19 # Zbar = 100, Kn = 1e-4
#NI=7.0442e18 # Zbar = 100, Kn = 1e-3
NI=1.4442e18 # Zbar = 100, Kn = 5e-3
if [ $ZBAR -le 10 ] ; then
   #NI=3.522e29 # Zbar = 1, Kn = 1e-10
   NI=3.522e22 # Zbar = 1, Kn = 1e-3
   #NI=3.522e21 # Zbar = 1, Kn = 1e-2
   #NI=0.704e21 # Zbar = 1, Kn = 5e-2
fi
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -Z $ZBAR -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -E0 1.0 -dE 0.001
cp results/tmp/C7_1_profiles.* results/fe_analysis/Ecorrect_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Ecorrect_data/fe_point_Ecorrect.txt

mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -Z $ZBAR -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -dE 1.0
cp results/tmp/C7_1_profiles.* results/fe_analysis/Emimic_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Emimic_data/fe_point_Emimic.txt

cd results/fe_analysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -cl $CL -n $NI -xp $XPOINT --Ecorrect --labelEcorrect 'C7' #--Emimic --labelEmimic 'C7*' --AWBSoriginal
fi
