#! /bin/bash
#SIGMA=8.1027575e17 ## Matching the SH diffusive flux.
SIGMA=2.425e17 ## Matching the CHIC diffusive flux.
CL=7.09 # Coulomb logarithm.

NPROC=8

RS=9
F1ORDER=4
F0ORDER=3

XPOINT=0.1605 # in cm qSH maximum

#MING=2000
MING=100

# Challenge SNB ;)
#MING=25
#F1ORDER=1
#F0ORDER=0

DIRroot=$PWD
DIRanalysis="results/fe_analysis/"
DIRout="C7E/"

# Philippe's test!
PROBLEM=10
L=0.188175

## Prepare data profiles.
cd VFPdata
cd GD_Hohlraum
#python ../loadPhilippegrace.py -f gdhohlraum_xmic_ne1e20cm3_interp -o ne_ -x -s -m 2e18 #-m 2e24
python ../loadPhilippegrace.py -f gdhohlraum_xmic_ne1e20cm3_interp -o ne_ -x -m 1e20 #-s
python ../loadPhilippegrace.py -f gdhohlraum_xmic_Z_interp -o Zbar_ -x #-s
python ../loadPhilippegrace.py -f gdhohlraum_xmic_10ps_TekeV_interp -o Te_ -x -m 1e3 #-s
python ../loadPhilippegrace.py -f gdhohlraum_xmic_10ps_IMPACTWcm2 -o Q_ #-s
python ../loadPhilippegrace.py -f gdhohlraum_xmic_10ps_separatedsnbWcm2 -o Q_ #-s
python ../loadPhilippegrace.py -f gdhohlraum_xmic_10ps_LocalWcm2 -o Q_ #-s
cd ..
cp GD_Hohlraum/Te_gdhohlraum_xmic_10ps_TekeV_interp.txt temperature.dat
cp GD_Hohlraum/ne_gdhohlraum_xmic_ne1e20cm3_interp.txt ne.dat
cp GD_Hohlraum/Zbar_gdhohlraum_xmic_Z_interp.txt zbar.dat
cp GD_Hohlraum/Q_gdhohlraum_xmic_10ps_IMPACTWcm2.txt flux1.dat
cp GD_Hohlraum/Q_gdhohlraum_xmic_10ps_separatedsnbWcm2.txt flux2.dat
cp GD_Hohlraum/Q_gdhohlraum_xmic_10ps_LocalWcm2.txt flux3.dat
cd ..


# Run C7.
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -sigma $SIGMA -cl $CL -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -dE 0.00000001 -Em 200 | tee C7E.out
## Run C7 10ps MING=3000 converged.
##mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -cl $CL -L $L -xp $XPOINT -minG $MING -s 3 -cfl 1e10 -S0 1.0 -dE 0.0000001 -Em 150 | tee C7E.out

cp results/tmp/C7_1_profiles.* results/fe_analysis/C7_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/C7_data/fe_point_C7.txt
cp results/tmp/C7_1_fe_pointmax.txt results/fe_analysis/C7_data/fe_pointmax_C7.txt
# Store the output file.
#cp C7E.out $DIRanalysis$DIRout"P10_Z"$ZBAR"_"$NAME".output"
# Perform analysis.
cd $DIRanalysis
python C7_analysis.py -N $NPROC -s $SIGMA -cl $CL -C7 --labelFluxExt1 IMPACT --labelFluxExt2 SNBs --labelFluxExt3 local --pltshow -SH #-xp #--Efield #-xp #--vlimshow #--AWBSstar #--AWBSoriginal

# Safe figs.

cd $DIRroot
