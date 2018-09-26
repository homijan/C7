#! /bin/bash
#SIGMA=8.1027575e17 ## Matching the SH diffusive flux.
SIGMA=2.425e17 ## Matching the CHIC diffusive flux.
CL=7.09 # Coulomb logarithm.

NPROC=2

RS=9
F1ORDER=2
F0ORDER=1

XPOINT=0.1605 # in cm qSH maximum

MING=1000
#MING=250

MAXITER=100
#MAXITER=40

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
python $DIRroot/VFPdata/loadPhilippegrace.py -f gdhohlraum_xmic_10ps_TekeV_interp -o _Te_ -mx 1e-4 -my 1e3 #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f gdhohlraum_xmic_ne1e20cm3_interp -o _ne_ -mx 1e-4 -my 1e20 #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f gdhohlraum_xmic_Z_interp -o _Zbar_ -mx 1e-4 #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f gdhohlraum_xmic_10ps_IMPACTWcm2 -o _Q_ #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f gdhohlraum_xmic_10ps_separatedsnbWcm2 -o _Q_ #-s
python $DIRroot/VFPdata/loadPhilippegrace.py -f gdhohlraum_xmic_10ps_LocalWcm2 -o _Q_ #-s
cd ..
cp GD_Hohlraum/_Te_gdhohlraum_xmic_10ps_TekeV_interp.txt temperature.dat
cp GD_Hohlraum/_ne_gdhohlraum_xmic_ne1e20cm3_interp.txt ne.dat
cp GD_Hohlraum/_Zbar_gdhohlraum_xmic_Z_interp.txt zbar.dat
cp GD_Hohlraum/_Q_gdhohlraum_xmic_10ps_IMPACTWcm2.txt flux1.dat
cp GD_Hohlraum/_Q_gdhohlraum_xmic_10ps_separatedsnbWcm2.txt flux2.dat
cp GD_Hohlraum/_Q_gdhohlraum_xmic_10ps_LocalWcm2.txt flux3.dat
cd ..


# Run C7.
mpirun -np $NPROC C7 -p $PROBLEM -m data/segment01.mesh -rs $RS -tf 0.0 -ok $F1ORDER -ot $F0ORDER -no-vis -fa -print -sigma $SIGMA -cl $CL -L $L -xp $XPOINT -minG $MING -s 2 -cfl 1e10 -S0 1.0 -dE 0.00000001 -Em $MAXITER -xn 1 | tee C7E.out

cp results/tmp/C7_1_profiles.* results/fe_analysis/C7_data/
cp results/tmp/C7_1_fe_point.txt results/fe_analysis/C7_data/fe_point_C7.txt
cp results/tmp/C7_1_fe_pointmax.txt results/fe_analysis/C7_data/fe_pointmax_C7.txt
# Store the output file.
#cp C7E.out $DIRanalysis$DIRout"P10_Z"$ZBAR"_"$NAME".output"
# Perform analysis.
cd $DIRanalysis
python C7_analysis.py -N $NPROC -s $SIGMA -cl $CL -fs 20 --labelUseC7 AP1 --labelFluxExt1 Impact --labelFluxExt2 SNBs --pltshow -pn -pZ -pT #--labelFluxExt3 local #-SH #-xp #--Efield #--vlimshow #--AWBSstar #--AWBSoriginal

# Safe figs.

cd $DIRroot
