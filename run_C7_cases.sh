#! /bin/bash
SIGMA=8.1027575e17 ## Matching the SH diffusive flux.
CL=7.09 # Coulomb logarithm.

NPROC=8

RS=6
F1ORDER=4
F0ORDER=3

TMAX=1000
TMIN=100
TGRAD=180

XPOINT=0.046775 # in cm qSH maximum

#MING=2000
MING=100

# Challenge SNB ;)
#MING=25
#F1ORDER=1
#F0ORDER=0

L=0.1
PROBLEM=5
#PROBLEM=9

ZBAR=1
#ZBAR=5
#ZBAR=10
#ZBAR=20
#ZBAR=100


## ZBAR 100
#NI=2e25 # Zbar = 100 -> Kn 1e-10
#NI=2e20 # Zbar = 100 -> Kn 1e-5
#NI=2e19 # Zbar = 100 -> Kn 1e-4
#NI=4e18 # Zbar = 100 -> Kn 5e-4
#NI=2e18 # Zbar = 100 -> Kn 1e-3
#NI=2e17 # Zbar = 100 -> Kn 1e-2
#NI=2e16 # Zbar = 100 -> Kn 1e-1
## ZBAR 20 
#NI=4.8e21 # Zbar = 20 -> Kn 1e-5
#NI=4.8e20 # Zbar = 20 -> Kn 1e-4
#NI=4.8e19 # Zbar = 20 -> Kn 1e-3 
#NI=4.8e18 # Zbar = 20 -> Kn 1e-2
#NI=0.963e18 # Zbar = 20 -> Kn 5e-2
## ZBAR 10
#NI=1.84e22 # Zbar = 10 -> Kn 1e-5
#NI=1.84e21 # Zbar = 10 -> Kn 1e-4
#NI=1.84e20 # Zbar = 10 -> Kn 1e-3
#NI=0.368e20 # Zbar = 10 -> Kn 5e-3
#NI=1.84e19 # Zbar = 10 -> Kn 1e-2
#NI=0.73e19 # Zbar = 10 -> Kn 2.2e-2
#NI=0.368e19 # Zbar = 10 -> Kn 5e-2
## ZBAR 5 
#NI=6.74e22 # Zbar = 5 -> Kn 1e-5
#NI=6.74e21 # Zbar = 5 -> Kn 1e-4
#NI=6.74e20 # Zbar = 5 -> Kn 1e-3
#NI=6.74e19 # Zbar = 5 -> Kn 1e-2
#NI=1.348e19 # Zbar = 5 -> Kn 5e-2
## ZBAR 1
#NI=1e29 # Zbar = 1 -> Kn 1e-10
#NI=1e24 # Zbar = 1 -> Kn 1e-5
#NI=1e23 # Zbar = 1 -> Kn 1e-4
#NI=1e22 # Zbar = 1 -> Kn 1e-3
#NI=1e21 # Zbar = 1 -> Kn 1e-2
#NI=2.02e20 # Zbar = 1 -> Kn 5e-2
#NI=1e20 # Zbar = 1 -> Kn 1e-1

DIRroot=$PWD
DIRanalysis="results/fe_analysis/"
DIRout="C7E/"
#declare -a Zarray=("1" "1" "1")
## full flux, half flux, fifth flux
#declare -a NIarray=("1.428e29" "6.25e20" "1.65e20")
#declare -a KNarray=("Kn1e-10" "Kn25e-3" "Kn72e-3")
#declare -a Zarray=("3" "3" "3")
## full flux, half flux, fifth flux
#declare -a NIarray=("2.38e28" "1.6e20" "4.1e19")
#declare -a KNarray=("Kn1e-10" "Kn16e-3" "Kn46e-3")
#declare -a Zarray=("10" "10" "10")
## full flux, half flux, fifth flux
#declare -a NIarray=("2.6e27" "2.65e19" "7.7e18")
#declare -a KNarray=("Kn1e-10" "Kn1e-2" "Kn27e-3")
#declare -a Zarray=("20" "20" "20")
## full flux, half flux, fifth flux
#declare -a NIarray=("6.8e26" "8.1e18" "2.5e18")
#declare -a KNarray=("Kn1e-10" "Kn84e-4" "Kn21e-3")

## Complete set of test cases.
declare -a Zarray=("1" "1" "1" "3" "3" "3" "10" "10" "10")
declare -a NIarray=("1.428e29" "6.25e20" "1.65e20" "2.38e28" "1.6e20" "4.1e19"  "2.6e27" "2.65e19" "7.7e18")
declare -a KNarray=("Kn1e-10" "Kn25e-3" "Kn72e-3" "Kn1e-10" "Kn16e-3" "Kn46e-3"  "Kn1e-10" "Kn1e-2" "Kn27e-3")

# get length of an array
length=${#Zarray[@]}

# use for loop to read all values and indexes
for (( i=0; i<${length}; i++ ));
do
# assign iterated values
ZBAR=${Zarray[$i]}
NI=${NIarray[$i]}
KN=${KNarray[$i]}
echo "ZBAR: " $ZBAR " NI: " $NI
# Run C7.
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -Z $ZBAR -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -dE 0.005 -Em 100 | tee C7E.out
cp results/tmp/C7_1_profiles.* results/fe_analysis/Ecorrect_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Ecorrect_data/fe_point_Ecorrect.txt
cp results/tmp/C7_1_fe_pointmax.txt results/fe_analysis/Ecorrect_data/fe_pointmax_Ecorrect.txt
# Store the output file.
cp C7E.out $DIRanalysis$DIRout"P5_Z"$ZBAR"_"$KN".output"
# Perform analysis.
cd $DIRanalysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -cl $CL -n $NI --Ecorrect --labelEcorrect 'C7E' --vlimshow #--pltshow #-xp #--AWBSstar #--AWBSoriginal
# Safe figs.
cp heatflux.png $DIRout"P5_heatflux_Z"$ZBAR"_"$KN".png"
cp kinetics.png $DIRout"P5_kinetics_Z"$ZBAR"_"$KN".png"
cd $DIRroot
done
