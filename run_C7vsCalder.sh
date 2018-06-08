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
XPOINT=0.0442 # in cm qSH maximum
L=0.095 #0.07
declare -a Zarray=("2")
declare -a NIarray=("2.5e20") # ne = 5e20
cd VFPdata/Calder/z2/z2tanh50microns_Bfield0/Te
python $DIRroot/VFPdata/loadPhilippegrace.py -f Te_00110000.txt -o _Te_ -mx 1e-4 -my 1e3 -cf 950 #700 -s
cp _Te_Te_00110000.txt.txt $DIRroot/VFPdata/temperature.dat
cd $DIRroot
cd VFPdata/Calder/z2/z2tanh50microns_Bfield0/qx
python $DIRroot/VFPdata/loadPhilippegrace.py -f qx_00110000.txt -o _Q_ -cf 950 #700 -s
cp _Q_qx_00110000.txt.txt $DIRroot/VFPdata/flux1.dat
cd $DIRroot
#cd VFPdata/Calder/z2/z2tanh50microns_Bfield0/Ex
#python $DIRroot/VFPdata/loadPhilippegrace.py -f Ex_00110000.txt -o _E_ -cf 950 -s
#cp _E_Ex_00110000.txt.txt $DIRroot/VFPdata/Efield1.dat
#cd $DIRroot
#cd VFPdata/Calder/z2/z2tanh50microns_Bfield0/f
#python $DIRroot/VFPdata/loadPhilippegrace.py -f f_00110000.txt -o _F0_ -cf 950 -s
#cp _F0_f_00110000.txt.txt $DIRroot/VFPdata/F0distribution1.dat
#cd $DIRroot

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
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Z $ZBAR -n0 $NUS0 -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -dE 0.01 -Em $MAXITER -xn 1 | tee C7E.out

cp results/tmp/C7_1_profiles.* results/fe_analysis/C7_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/C7_data/fe_point_C7.txt
cp results/tmp/C7_1_fe_pointmax.txt results/fe_analysis/C7_data/fe_pointmax_C7.txt

# Store the output file.
#cp C7E.out $DIRanalysis$DIRout"P9_Z"$ZBAR"_"$NAME".output"
# Perform analysis.
cd $DIRanalysis
python C7_analysis.py -N $NPROC -s $SIGMA -cl $CL -fs 20 --labelUseC7 C7 --labelFluxExt1 Calder --pltshow --pltTe #-SH #-lD1 Calder --Efield --labelEfieldExt1 Calder

# Safe figs.
### CASE 1 ###
cp heatflux.png $DIRroot/VFPdata/C7_Calder_case1_heatflux.png
cp kinetics.png $DIRroot/VFPdata/C7_Calder_case1_kinetics.png

cd $DIRroot
done
