#! /bin/bash
SIGMA=8.1027575e17 ## Matching the SH diffusive flux.
CL=7.09 # Coulomb logarithm.

NPROC=7

RS=6
F1ORDER=4
F0ORDER=3

TMAX=1000
TMIN=100
TGRAD=180

XPOINT=0.046775 # in cm qSH maximum

#MING=1000
MING=250

# Challenge SNB ;)
#MING=25
#F1ORDER=1
#F0ORDER=0

L=0.1
PROBLEM=5

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
#declare -a Zarray=("100" "100" "100")
## full flux, half flux, fifth flux
#declare -a NIarray=("2.828e25" "6.0e17" "1.1e17")
#declare -a KNarray=("Kn1e-10" "Kn53e-4" "Kn28e-3")


## Complete set of test cases.
#declare -a Zarray=("1" "1" "1" "3" "3" "3" "10" "10" "10" "20" "20" "20" "100" "100" "100")
#declare -a NIarray=("1.428e29" "6.25e20" "1.65e20" "2.38e28" "1.6e20" "4.1e19"  "2.6e27" "2.65e19" "7.7e18" "6.8e26" "8.1e18" "2.5e18" "2.828e25" "6.0e17" "1.1e17")
#declare -a NAMEarray=("fullF" "halfF" "fifthF" "fullF" "halfF" "fifthF"  "fullF" "halfF" "fifthF" "fullF" "halfF" "fifthF" "fullF" "halfF" "fifthF")

#declare -a Zarray=("100")
#declare -a nuS0array=("0.55")
#declare -a Zarray=("20")
#declare -a nuS0array=("0.5535")
#declare -a Zarray=("10")
#declare -a nuS0array=("0.5485")
#declare -a Zarray=("4")
#declare -a nuS0array=("0.53225")
#declare -a Zarray=("2")
#declare -a nuS0array=("0.5035")
declare -a Zarray=("1")
#declare -a nuS0array=("0.455")
## full flux, half flux, fifth flux
declare -a NIarray=("1.428e29")
#declare -a NIarray=("1.428e20")
#declare -a NIarray=("1.428e21")
declare -a NAMEarray=("test")

declare -a nuS0array=("0.5")

# get length of an array
length=${#Zarray[@]}

# use for loop to read all values and indexes
for (( i=0; i<${length}; i++ ));
do
# assign iterated values
ZBAR=${Zarray[$i]}
NUS0=${nuS0array[$i]}
NI=${NIarray[$i]}
NAME=${NAMEarray[$i]}
echo "ZBAR: " $ZBAR " NI: " $NI
# Run C7.
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -Tmax $TMAX -Tmin $TMIN -Tgrad $TGRAD -Z $ZBAR -n0 $NUS0 -cl $CL -ni $NI -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -dE 0.01 -Em 100 | tee C7E.out

cp results/tmp/C7_1_profiles.* results/fe_analysis/Ecorrect_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/Ecorrect_data/fe_point_Ecorrect.txt
cp results/tmp/C7_1_fe_pointmax.txt results/fe_analysis/Ecorrect_data/fe_pointmax_Ecorrect.txt
# Store the output file.
cp C7E.out $DIRanalysis$DIRout"P5_Z"$ZBAR"_"$NAME".output"
# Perform analysis.
cd $DIRanalysis
python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -cl $CL -n $NI --Ecorrect --labelEcorrect 'C7' --pltshow #-xp #--vlimshow #-xp #--AWBSstar #--AWBSoriginal
#python C7_AWBS_SH_analysis.py -N $NPROC -Z $ZBAR -cl $CL -n $NI --Ecorrect --pltshow --AWBSstar --AWBSoriginal
# Safe figs.
cp heatflux.png $DIRout"P5_heatflux_Z"$ZBAR"_"$NAME".png"
cp kinetics.png $DIRout"P5_kinetics_Z"$ZBAR"_"$NAME".png"
cd $DIRroot
done
