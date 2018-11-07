#! /bin/bash
SIGMA=8.1027575e17 ## Matching the SH diffusive flux.
CL=7.09 # Coulomb logarithm.

NPROC=2

RS=6
F1ORDER=4
F0ORDER=3
#F1ORDER=1
#F0ORDER=0

XPOINT=0.0480 # in cm qSH maximum

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

# Universal value.
NUS0=0.5

declare -a NAMEarray=("case")
## Prepare data profiles.

### CASE 1 ###
#XPOINT=0.06441 # in cm qCalder maximum
XPOINT=0.075 # in cm q nonlocal
L=0.094
ZBAR=2
NE=5e20
declare -a Zarray=("2")
declare -a NIarray=("2.5e20") # ne = 5e20
#declare -a NIarray=("1e19")
cd VFPdata/Calder/z2/z2tanh50microns_Bfield0/Te
python $DIRroot/VFPdata/loadPhilippegrace.py -f Te_00110000.txt -o _Te_ -mx 1e-4 -my 1e3 -cf 940 #700 -s
cp _Te_Te_00110000.txt.txt $DIRroot/VFPdata/temperature.dat
cd $DIRroot
cd VFPdata/Calder/z2/z2tanh50microns_Bfield0/qx
python $DIRroot/VFPdata/loadPhilippegrace.py -f qx_00110000.txt -o _Q_ -cf 940 #700 -s
cp _Q_qx_00110000.txt.txt $DIRroot/VFPdata/flux1.dat
cd $DIRroot
#cd VFPdata/Calder/z2/z2tanh50microns_Bfield0/Ex
#python $DIRroot/VFPdata/loadPhilippegrace.py -f Ex_00110000.txt -o _E_ -cf 940 -s
#cp _E_Ex_00110000.txt.txt $DIRroot/VFPdata/Efield1.dat
#cd $DIRroot
cd VFPdata/Calder/z2/z2tanh50microns_Bfield0/f
## Heat flux maximum.
#python $DIRroot/VFPdata/loadPhilippegrace.py -f F0F1x_Calder_Zeq2_tanh_50mic_f0f1_1.10e-11_644.1mic.txt -o _F0_ #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f F0F1x_Calder_Zeq2_tanh_50mic_f0f1_1.10e-11_644.1mic.txt -o _F1_ --column 2 #-s
#cp _F0_F0F1x_Calder_Zeq2_tanh_50mic_f0f1_1.10e-11_644.1mic.txt.txt $DIRroot/VFPdata/F0distribution1.dat
#cp _F1_F0F1x_Calder_Zeq2_tanh_50mic_f0f1_1.10e-11_644.1mic.txt.txt $DIRroot/VFPdata/F1distribution1.dat
## Preheat position
python $DIRroot/VFPdata/loadPhilippegrace.py -f F0F1x_Calder_Zeq2_tanh_50mic_f0f1_1.10e-11_749.4mic.txt -o _F0_ #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f F0F1x_Calder_Zeq2_tanh_50mic_f0f1_1.10e-11_749.4mic.txt -o _F1_ --column 2 #-s
cp _F0_F0F1x_Calder_Zeq2_tanh_50mic_f0f1_1.10e-11_749.4mic.txt.txt $DIRroot/VFPdata/F0distribution1.dat
cp _F1_F0F1x_Calder_Zeq2_tanh_50mic_f0f1_1.10e-11_749.4mic.txt.txt $DIRroot/VFPdata/F1distribution1.dat

## SNBE output.
python $DIRroot/SNBE/SNBE.py -ne $NE -Z $ZBAR -Tinf $DIRroot/VFPdata/temperature.dat -Qo $DIRroot/VFPdata/flux2.dat -F1o $DIRroot/VFPdata/F1distribution2.dat -pt $XPOINT
cd $DIRroot

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
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Z $ZBAR -n0 $NUS0 -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -dE 0.01 -Em $MAXITER -xn 0.0 | tee C7E.out

cp results/tmp/C7_1_profiles.* results/fe_analysis/C7_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/C7_data/fe_point_C7.txt
cp results/tmp/C7_1_fe_pointmax.txt results/fe_analysis/C7_data/fe_pointmax_C7.txt

# Store the output file.
#cp C7E.out $DIRanalysis$DIRout"P9_Z"$ZBAR"_"$NAME".output"
# Perform analysis.
cd $DIRanalysis
# FLUX MAXIMUM DISTRIBUTION
#python C7_analysis.py -N $NPROC -s $SIGMA -cl $CL -fs 18 --labelUseC7 AP1 --labelFluxExt1 Calder -lF2 SNB --pltshow --pltTe -SH -xp -lD1 Calder -lD2 SNB --plotmultvTh 7 -Tpts 0.075 0.06441 #--Efield --labelEfieldExt1 Calder
# NONLOCAL DISTRIBUTION
python C7_analysis.py -N $NPROC -s $SIGMA -cl $CL -fs 18 --labelUseC7 AP1 --labelFluxExt1 Calder -lF2 SNB --pltshow --pltTe -SH -xp -lD1 Calder -lD2 SNB --plotmultvTh 14 -Tpts 0.075 0.06441 #--Efield --labelEfieldExt1 Calder

# Safe figs.
### CASE 1 ###
cp heatflux.png $DIRroot/VFPdata/C7_Calder_case1_heatflux.png
#cp kinetics.png $DIRroot/VFPdata/C7_Calder_case1_kinetics.png
cp kinetics.png $DIRroot/VFPdata/C7_Calder_case1_nonlocal_kinetics.png

cd $DIRroot
done
