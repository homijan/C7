#! /bin/bash
SIGMA=8.1027575e17 ## Matching the SH diffusive flux.
CL=7.09 # Coulomb logarithm.

NPROC=8

RS=6
F1ORDER=4
F0ORDER=3
#F1ORDER=2
#F0ORDER=1

XPOINT=0.046775 # in cm qSH maximum

#MING=2000
MING=250

MAXITER=100

#F1ORDER=3
#F0ORDER=2
#MING=100
#MAXITER=6

# Challenge SNB ;)
#MING=25
#F1ORDER=1
#F0ORDER=0


DIRroot=$PWD
DIRanalysis="results/fe_analysis/"
DIRout="C7E/"

# Philippe's test!
PROBLEM=9
#L=0.07
#declare -a Zarray=("2")
#declare -a NIarray=("2.5e20")
#NUS0=0.5035
L=0.014
declare -a Zarray=("10")
declare -a NIarray=("5e19")
NUS0=0.5485

declare -a NAMEarray=("case")
## Prepare data profiles.
cd VFPdata
cd Pascal
#python ../loadPhilippegrace.py -f Te_Aladin_tanh_50mic_Zeq2_B0_1.20e-11.txt -o Te_ -x -m 1e3 -s
#python ../loadPhilippegrace.py -f FluxX_Aladin_tanh_50mic_Zeq2_B0_1.20e-11.txt -o Q_ -s
python ../loadPhilippegrace.py -f Te_Aladin_tanh_10mic_Zeq10_B0_1.20e-11.txt -o Te_ -x -m 1e3 -s
python ../loadPhilippegrace.py -f FluxX_Aladin_tanh_10mic_Zeq10_B0_1.20e-11.txt -o Q_ -s
cd ..
#cp Pascal/Te_Te_Aladin_tanh_50mic_Zeq2_B0_1.20e-11.txt.txt temperature.dat
#cp Pascal/Q_FluxX_Aladin_tanh_50mic_Zeq2_B0_1.20e-11.txt.txt flux1.dat
cp Pascal/Te_Te_Aladin_tanh_10mic_Zeq10_B0_1.20e-11.txt.txt temperature.dat
cp Pascal/Q_FluxX_Aladin_tanh_10mic_Zeq10_B0_1.20e-11.txt.txt flux1.dat
cd ..

# get length of an array
length=${#Zarray[@]}

# use for loop to read all values and indexes
for (( i=0; i<${length}; i++ ));
do
# assign iterated values
ZBAR=${Zarray[$i]}
NI=${NIarray[$i]}
NAME=${NAMEarray[$i]}
echo "ZBAR: " $ZBAR " NI: " $NI
# Run C7.
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Z $ZBAR -n0 $NUS0 -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -dE 0.01 -Em $MAXITER | tee C7E.out

cp results/tmp/C7_1_profiles.* results/fe_analysis/C7_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/C7_data/fe_point_C7.txt
cp results/tmp/C7_1_fe_pointmax.txt results/fe_analysis/C7_data/fe_pointmax_C7.txt

# Store the output file.
#cp C7E.out $DIRanalysis$DIRout"P9_Z"$ZBAR"_"$NAME".output"
# Perform analysis.
cd $DIRanalysis
python C7_analysis.py -N $NPROC -s $SIGMA -cl $CL -C7 --labelFluxExt1 Aladin --pltshow -SH --pltTe

# Safe figs.
## 12ps.
#cp heatflux.png $DIRroot/VFPdata/C7_heatflux_12ps.png
#cp kinetics.png $DIRroot/VFPdata/C7_kinetics_12ps.png
cd $DIRroot
done
