cp $path/all_scripts/001.KMC_MD/0.AgAu.eam.alloy $(pwd)
cp $path/all_scripts/001.KMC_MD/1.input_MD 1.input
cp ../step4_MD/dump.lammpstrj $(pwd)
cp ../step4_MD/2.coordinates $(pwd)

sed '4!d' dump.lammpstrj > count
a=`cat count`
tail -$a dump.lammpstrj > NW.xyz
head -29 2.coordinates > header
echo $a

cat header NW.xyz > 2.coordinates
rm dump.lammpstrj count header NW.xyz tmp* max min

LAMMPSB_EXE=/apps/scratch/compile/new/v100/lammps-7Aug19/src/lmp_v100
mpirun -np 20 $LAMMPSB_EXE -sf gpu -pk gpu 2 -in 1.input > output

sed '4!d' dump.lammpstrj > count
a=`cat count`
tail -$a dump.lammpstrj > 18nm_trunc_octa
echo $a
cp 18nm_trunc_octa $path/
mv 18nm_trunc_octa coordinates_final
