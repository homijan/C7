#! /bin/bash
SIGMA=8.1027575e17 ## Matching the SH diffusive flux.
CL=7.09 # Coulomb logarithm.

NPROC=2

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

# Universal value.
NUS0=0.5

declare -a NAMEarray=("case")
## Prepare data profiles.

### CASE 1 ###
XPOINT=2.1115 # in cm qSH maximum
L=3.495
declare -a Zarray=("2")
declare -a NIarray=("2.5e20") # ne = 5e20
cd VFPdata/Impact/tanh/z2/NoB/2500mic_ramp/241.3ps
python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_TkeV -o _Te_ -mx 1e-4 -my 1e3 #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_QxWcm2 -o _Q_ #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_ExkVpermm -o _E_ -my 1e6 #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_v_f0_near_qmax_21125.0mic -o _F0_ --row 1 -mx 1e2 -mom 2 -my 1e-12 #-s # s^-3 / cm^-6
python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_v_f1x_at_qmax_21100.0mic -o _F1_ --row 1 -mx 1e2 -mom 5 -my 4.5547e-40 #-s # me/2 * s^-3 / cm^-6
## Copy the input data files to dedicated directory. 
cp _Te_IMPACT_xmic_TkeV.txt $DIRroot/VFPdata/temperature.dat
cp _Q_IMPACT_xmic_QxWcm2.txt $DIRroot/VFPdata/flux1.dat
cp _E_IMPACT_xmic_ExkVpermm.txt $DIRroot/VFPdata/Efield1.dat
cp _F0_IMPACT_v_f0_near_qmax_21125.0mic.txt $DIRroot/VFPdata/F0distribution1.dat
cp _F1_IMPACT_v_f1x_at_qmax_21100.0mic.txt $DIRroot/VFPdata/F1distribution1.dat
cd $DIRroot

### CASE 2 ###
#XPOINT=0.1081 # in cm qSH maximum
#L=0.17475
#declare -a Zarray=("2")
#declare -a NIarray=("2.5e20") # ne = 5e20
#cd VFPdata/Impact/tanh/z2/NoB/125mic_ramp/24.1ps
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_TkeV -o _Te_ -mx 1e-4 -my 1e3 #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_QxWcm2 -o _Q_ #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_ExkVpermm -o _E_ -my 1e6 #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_v_f0_near_qmax_1081.2mic -o _F0_ --row 1 -mx 1e2 -mom 2 -my 1e-12 #-s # s^-3 / cm^-6
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_v_f1x_at_qmax_1080.0mic -o _F1_ --row 1 -mx 1e2 -mom 5 -my 4.5547e-40 #-s # me/2 * s^-3 / cm^-6
### Copy the input data files to dedicated directory. 
#cp _Te_IMPACT_xmic_TkeV.txt $DIRroot/VFPdata/temperature.dat
#cp _Q_IMPACT_xmic_QxWcm2.txt $DIRroot/VFPdata/flux1.dat
#cp _E_IMPACT_xmic_ExkVpermm.txt $DIRroot/VFPdata/Efield1.dat
#cp _F0_IMPACT_v_f0_near_qmax_1081.2mic.txt $DIRroot/VFPdata/F0distribution1.dat
#cp _F1_IMPACT_v_f1x_at_qmax_1080.0mic.txt $DIRroot/VFPdata/F1distribution1.dat
#cd $DIRroot

### CASE 3 ###
#XPOINT=0.0437 # in cm qSH maximum
#L=0.07
#declare -a Zarray=("2")
#declare -a NIarray=("2.5e20") # ne = 5e20
#cd VFPdata/Impact/tanh/z2/NoB/50mic_ramp/12.1ps
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_TkeV -o _Te_ -mx 1e-4 -my 1e3 #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_QxWcm2 -o _Q_ #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_ExkVpermm -o _E_ -my 1e6 #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_v_f0_near_qmax_437.5mic -o _F0_ --row 1 -mx 1e2 -mom 2 -my 1e-12 #-s # s^-3 / cm^-6
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_v_f1x_at_qmax_437.0mic -o _F1_ --row 1 -mx 1e2 -mom 5 -my 4.5547e-40 #-s # me/2 * s^-3 / cm^-6
### Copy the input data files to dedicated directory. 
#cp _Te_IMPACT_xmic_TkeV.txt $DIRroot/VFPdata/temperature.dat
#cp _Q_IMPACT_xmic_QxWcm2.txt $DIRroot/VFPdata/flux1.dat
#cp _E_IMPACT_xmic_ExkVpermm.txt $DIRroot/VFPdata/Efield1.dat
#cp _F0_IMPACT_v_f0_near_qmax_437.5mic.txt $DIRroot/VFPdata/F0distribution1.dat
#cp _F1_IMPACT_v_f1x_at_qmax_437.0mic.txt $DIRroot/VFPdata/F1distribution1.dat
#cd $DIRroot

### CASE 4 ###
#XPOINT=0.04457 # in cm qSH maximum
#L=0.07
#declare -a Zarray=("2")
#declare -a NIarray=("2.5e20") # ne = 5e20
#cd VFPdata/Impact/tanh/z2/NoB/25mic_ramp/12.1ps
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_TkeV -o _Te_ -mx 1e-4 -my 1e3 #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_QxWcm2 -o _Q_ #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_ExkVpermm -o _E_ -my 1e6 #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_v_f0_near_qmax_445.8mic -o _F0_ --row 1 -mx 1e2 -mom 2 -my 1e-12 #-s # s^-3 / cm^-6
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_v_f1x_at_qmax_445.5mic -o _F1_ --row 1 -mx 1e2 -mom 5 -my 4.5547e-40 #-s # me/2 * s^-3 / cm^-6
### Copy the input data files to dedicated directory. 
#cp _Te_IMPACT_xmic_TkeV.txt $DIRroot/VFPdata/temperature.dat
#cp _Q_IMPACT_xmic_QxWcm2.txt $DIRroot/VFPdata/flux1.dat
#cp _E_IMPACT_xmic_ExkVpermm.txt $DIRroot/VFPdata/Efield1.dat
#cp _F0_IMPACT_v_f0_near_qmax_445.8mic.txt $DIRroot/VFPdata/F0distribution1.dat
#cp _F1_IMPACT_v_f1x_at_qmax_445.5mic.txt $DIRroot/VFPdata/F1distribution1.dat
#cd $DIRroot

### CASE 5 ###
#XPOINT=0.01347 # in cm qSH maximum
#L=0.021
#declare -a Zarray=("2")
#declare -a NIarray=("2.5e20") # ne = 5e20
#cd VFPdata/Impact/tanh/z2/NoB/7.5mic_ramp/6.0ps
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_TkeV -o _Te_ -mx 1e-4 -my 1e3 #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_QxWcm2 -o _Q_ #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_xmic_ExkVpermm -o _E_ -my 1e6 #-s
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_v_f0_near_qmax_134.8mic -o _F0_ --row 1 -mx 1e2 -mom 2 -my 1e-12 #-s # s^-3 / cm^-6
#python $DIRroot/VFPdata/loadPhilippegrace.py -f IMPACT_v_f1x_at_qmax_134.7mic -o _F1_ --row 1 -mx 1e2 -mom 5 -my 4.5547e-40 #-s # me/2 * s^-3 / cm^-6
### Copy the input data files to dedicated directory. 
#cp _Te_IMPACT_xmic_TkeV.txt $DIRroot/VFPdata/temperature.dat
#cp _Q_IMPACT_xmic_QxWcm2.txt $DIRroot/VFPdata/flux1.dat
#cp _E_IMPACT_xmic_ExkVpermm.txt $DIRroot/VFPdata/Efield1.dat
#cp _F0_IMPACT_v_f0_near_qmax_134.8mic.txt $DIRroot/VFPdata/F0distribution1.dat
#cp _F1_IMPACT_v_f1x_at_qmax_134.7mic.txt $DIRroot/VFPdata/F1distribution1.dat
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
python C7_analysis.py -N $NPROC -s $SIGMA -cl $CL --labelUseC7 C7 --labelFluxExt1 Impact --pltshow -SH --pltTe --Efield --labelEfieldExt1 Impact -lD1 Impact -xp

# Safe figs.
### CASE 1 ###
cp heatflux.png $DIRroot/VFPdata/C7_Impact_case1_heatflux.png
cp kinetics.png $DIRroot/VFPdata/C7_Impact_case1_kinetics.png
### CASE 2 ###
#cp heatflux.png $DIRroot/VFPdata/C7_Impact_case2_heatflux.png
#cp kinetics.png $DIRroot/VFPdata/C7_Impact_case2_kinetics.png
### CASE 3 ###
#cp heatflux.png $DIRroot/VFPdata/C7_Impact_case3_heatflux.png
#cp kinetics.png $DIRroot/VFPdata/C7_Impact_case3_kinetics.png
### CASE 4 ###
#cp heatflux.png $DIRroot/VFPdata/C7_Impact_case4_heatflux.png
#cp kinetics.png $DIRroot/VFPdata/C7_Impact_case4_kinetics.png
### CASE 5 ###
#cp heatflux.png $DIRroot/VFPdata/C7_Impact_case5_heatflux.png
#cp kinetics.png $DIRroot/VFPdata/C7_Impact_case5_kinetics.png

cd $DIRroot
done
