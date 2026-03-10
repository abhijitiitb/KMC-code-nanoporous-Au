cp $path/all_scripts/001.KMC_MD/header $(pwd)
cp $path/all_scripts/001.KMC_MD/0.AgAu.eam.alloy $(pwd)
cp $path/all_scripts/001.KMC_MD/1.input_minimize_energy 1.input
cp ../step1_KMC/NW.102.xyz $(pwd)

cp NW.102.xyz NW.xyz
sed -i "1,2d" NW.xyz

cat NW.xyz | awk '{print $1, $2, $3, $4}' > tmp && mv tmp NW.xyz

sed -i 's/Pt/2/g' NW.xyz
sed -i 's/Fe/1/g' NW.xyz

count=$(wc -l  < NW.xyz)
sed -i "2s/abc/${count}/g" header

a=1
for i in `seq 1 3`
  do
    a=$(($a+1))
    if [ $a -eq 2 ]; then sed -i "12s/abc/0.0/g" header; sed -i "12s/def/1000.0/g" header; fi
    if [ $a -eq 3 ]; then sed -i "13s/abc/0.0/g" header; sed -i "13s/def/1000.0/g" header; fi
    if [ $a -eq 4 ]; then sed -i "14s/abc/0.0/g" header; sed -i "14s/def/1000.0/g" header; fi
  done 

seq 1 1 $count > sequence
paste -d' '  sequence NW.xyz > tmp && mv tmp NW.xyz
cat header NW.xyz > 2.coordinates
rm sequence header NW.xyz NW.102.xyz 

LAMMPSB_EXE=/apps/scratch/compile/new/v100/lammps-7Aug19/src/lmp_v100
mpirun -np 20 $LAMMPSB_EXE -sf gpu -pk gpu 2 -in 1.input > output                                                 #run step
