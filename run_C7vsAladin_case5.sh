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
#MING=35

MAXITER=100

DIRroot=$PWD
DIRanalysis="results/fe_analysis/"
DIRout="C7E/"

# Philippe's test!
PROBLEM=9

# Universal value of electron collision frequency scaling for AWBS and BGK.
NUS0=0.5

## Prepare data profiles.

### CASE 5 ###  Z = 1 for AWBS paper.
#XPOINT=0.046 # in cm qSH maximum
#XPOINT=0.058 # in cm q nonlocal
L=0.07
ZBAR=1
NI=5e20
NE=5e20
declare -a Pointsarray=("0.046" "0.058")
declare -a Namesarray=("maximum" "nonlocal")

# get length of an array
length=${#Pointsarray[@]}

# use for loop to read all values and indexes
for (( i=0; i<${length}; i++ ));
do
# assign iterated values
XPOINT=${Pointsarray[$i]}
NAME=${Namesarray[$i]}

echo "ZBAR: " $ZBAR " NI: " $NI "Point: " $XPOINT

cd $DIRroot/VFPdata/Aladin/Aladin_cases4AWBSpaper/case5
python $DIRroot/VFPdata/loadPhilippegrace.py -f Te_Aladin_Zeq1_tanh_50mic_2.00e-11.txt -o _Te_ -mx 1e-4 -my 1e3 #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f FluxX_Aladin_Zeq1_tanh_50mic_2.00e-11.txt -o _Q_ #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f ElecX_Aladin_Zeq1_tanh_50mic_2.00e-11.txt -o _E_ #-m 0.33333334e-4 -s
## Copy 4C7 analysis generated data.
cp _Te_Te_Aladin_Zeq1_tanh_50mic_2.00e-11.txt.txt $DIRroot/VFPdata/temperature.dat
cp _Q_FluxX_Aladin_Zeq1_tanh_50mic_2.00e-11.txt.txt $DIRroot/VFPdata/flux1.dat
cp _E_ElecX_Aladin_Zeq1_tanh_50mic_2.00e-11.txt.txt $DIRroot/VFPdata/Efield1.dat
## Heat flux maximum.
if [ "$NAME" == "maximum" ]
then
echo $NAME
python $DIRroot/VFPdata/loadPhilippegrace.py -f F0F1x_Aladin_Zeq1_tanh_50mic_5e20_B0_20ps_f0f1b_2.00e-11_460mic.txt -o _F0_ #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f F0F1x_Aladin_Zeq1_tanh_50mic_5e20_B0_20ps_f0f1b_2.00e-11_460mic.txt -o _F1_ --column 2 #-s
cp _F0_F0F1x_Aladin_Zeq1_tanh_50mic_5e20_B0_20ps_f0f1b_2.00e-11_460mic.txt.txt $DIRroot/VFPdata/F0distribution1.dat
cp _F1_F0F1x_Aladin_Zeq1_tanh_50mic_5e20_B0_20ps_f0f1b_2.00e-11_460mic.txt.txt $DIRroot/VFPdata/F1distribution1.dat
fi
## Nonlocal at preheat point.
if [ "$NAME" == "nonlocal" ]
then
echo $NAME
python $DIRroot/VFPdata/loadPhilippegrace.py -f F0F1x_Aladin_Zeq1_tanh_50mic_5e20_B0_20ps_f0f1b_2.00e-11_580mic.txt -o _F0_ #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f F0F1x_Aladin_Zeq1_tanh_50mic_5e20_B0_20ps_f0f1b_2.00e-11_580mic.txt -o _F1_ --column 2 #-s
cp _F0_F0F1x_Aladin_Zeq1_tanh_50mic_5e20_B0_20ps_f0f1b_2.00e-11_580mic.txt.txt $DIRroot/VFPdata/F0distribution1.dat
cp _F1_F0F1x_Aladin_Zeq1_tanh_50mic_5e20_B0_20ps_f0f1b_2.00e-11_580mic.txt.txt $DIRroot/VFPdata/F1distribution1.dat
fi

## SNBE output.
python $DIRroot/SNBE/SNBE.py -ne $NE -Z $ZBAR -Tinf $DIRroot/VFPdata/temperature.dat -Qo $DIRroot/VFPdata/flux2.dat -Jo $DIRroot/VFPdata/jSNB.dat -F1o $DIRroot/VFPdata/F1distribution2.dat -pt $XPOINT
cd $DIRroot

# Run C7.
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Z $ZBAR -n0 $NUS0 -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -dE 0.01 -Em $MAXITER -xn 0.0 | tee C7E.out

cp results/tmp/C7_1_profiles.* results/fe_analysis/C7_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/C7_data/fe_point_C7.txt
cp results/tmp/C7_1_fe_pointmax.txt results/fe_analysis/C7_data/fe_pointmax_C7.txt

# Perform analysis.
cd $DIRanalysis
## FLUX MAXIMUM DISTRIBUTION
if [ "$NAME" == "maximum" ]
then
python C7_analysis.py -N $NPROC -s $SIGMA -cl $CL -fs 18 --labelUseC7 AP1 --labelFluxExt1 Aladin -lF2 SNB -SH -lD1 Aladin -hftit 'Heat flux Z = 1' -ktit 'Kinetics at maximum point' --pltshow --pltTe -xp --plotmultvTh 7 --Efield --labelEfieldExt1 Aladin -lD2 SNB -Tpts 0.058 0.046
fi
## NONLOCAL DISTRIBUTION
if [ "$NAME" == "nonlocal" ]
then
python C7_analysis.py -N $NPROC -s $SIGMA -cl $CL -fs 18 --labelUseC7 AP1 -lF1 Aladin -lF2 SNB -hftit 'Heat flux Z = 1' -ktit 'Kinetics at preheat point' -lD1 Aladin -lD2 SNB --pltshow --pltTe -xp -SH --plotmultvTh 14 -Tpts 0.058 0.046 --Efield --labelEfieldExt1 Aladin
fi

# Safe figs.
### CASE 5 ### Z = 1 for AWBS paper.
cp heatflux.png $DIRroot/VFPdata/C7_Aladin_case5_heatflux.png
cp Efield.png $DIRroot/VFPdata/C7_Aladin_case5_Efield.png
if [ "$NAME" == "maximum" ]
then
cp kinetics.png $DIRroot/VFPdata/C7_Aladin_case5_kinetics.png
fi
if [ "$NAME" == "nonlocal" ]
then
cp kinetics.png $DIRroot/VFPdata/C7_Aladin_case5_nonlocal_kinetics.png
fi

cd $DIRroot
done
