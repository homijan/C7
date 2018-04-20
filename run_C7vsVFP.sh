#! /bin/bash
SIGMA=8.1027575e17 ## Matching the SH diffusive flux.
CL=7.09 # Coulomb logarithm.

NPROC=8

RS=6
F1ORDER=4
F0ORDER=3

XPOINT=0.046775 # in cm qSH maximum

MING=2000
#MING=250

# Challenge SNB ;)
#MING=25
#F1ORDER=1
#F0ORDER=0

DIRroot=$PWD
DIRanalysis="results/fe_analysis/"
DIRout="C7E/"

# Philippe's test!
PROBLEM=9
L=0.035
declare -a Zarray=("4")
## full flux, half flux, fifth flux
declare -a NIarray=("1.25e20")
declare -a NAMEarray=("case")
## Prepare data profiles.
cd VFPdata
python loadgrace.py
## case 6ps
#cp Te_VFP_6ps.dat temperature.dat
#cp flux_VFP_6_ps.dat flux_VFP.dat
#cp flux_SNB_6_ps_r=2.dat flux_SNB.dat
## case 12ps
cp Te_VFP_12_ps.dat temperature.dat
cp flux_VFP_12_ps.dat flux_VFP.dat
cp flux_SNB_12_ps_r=2.dat flux_SNB.dat
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
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Z $ZBAR -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -dE 0.01 -Em 100 | tee C7E.out

cp results/tmp/C7_1_profiles.* results/fe_analysis/Ecorrect_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Ecorrect_data/fe_point_Ecorrect.txt
cp results/tmp/C7_1_fe_pointmax.txt results/fe_analysis/Ecorrect_data/fe_pointmax_Ecorrect.txt
# Store the output file.
cp C7E.out $DIRanalysis$DIRout"P9_Z"$ZBAR"_"$NAME".output"
# Perform analysis.
cd $DIRanalysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -cl $CL -n $NI --Ecorrect --labelEcorrect 'C7' --pltshow #--vlimshow #-xp #--AWBSstar #--AWBSoriginal
#python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -cl $CL -n $NI --Ecorrect --pltshow --AWBSstar --AWBSoriginal
# Safe figs.
## Philippe's test.
## 6ps.
#cp heatflux.png $DIRroot/VFPdata/C7_heatflux_6ps.png
#cp kinetics.png $DIRroot/VFPdata/C7_kinetics_6ps.png
## 12ps.
cp heatflux.png $DIRroot/VFPdata/C7_heatflux_12ps.png
cp kinetics.png $DIRroot/VFPdata/C7_kinetics_12ps.png
cd $DIRroot
done
