cp $path/all_scripts/001.KMC_MD/header $(pwd)
cp $path/all_scripts/001.KMC_MD/0.AgAu.eam.alloy $(pwd)
cp $path/all_scripts/001.KMC_MD/1.input_randomize_water_molecules 1.input
rm coarsening*
cp ../step3_MD/metal-solvent-without-overlap.xyz NW.xyz
sed -i "1,2d" NW.xyz 

sed -i 's/Au/2/g' NW.xyz
sed -i 's/Ag/1/g' NW.xyz
sed -i 's/W/3/g' NW.xyz

count=$(wc -l  < NW.xyz)
sed -i "2s/abc/${count}/g" header
echo $count
cat NW.xyz | awk '{print $1}' | grep "3" | wc -l > tmp
county=`cat tmp`
county=$(($count-$county+1))
echo $county

sed -i "7s/abc/${county}/g" 1.input
sed -i "7s/def/${count}/g" 1.input

a=1
for i in `seq 1 3`
  do
    a=$(($a+1))
    cat NW.xyz | awk -v a=$a '{print $a}' > tmp$a
    sort -k1n tmp$a > tmp && mv tmp tmp$a
    head -1 tmp$a > min && tail -1 tmp$a > max
    b=`cat min` && c=`cat max`
    c=$(echo $c + 3.1 | bc -l)
    if [ $a -eq 2 ]; then sed -i "12s/abc/${b}/g" header; sed -i "12s/def/${c}/g" header; fi
    if [ $a -eq 3 ]; then sed -i "13s/abc/${b}/g" header; sed -i "13s/def/${c}/g" header; fi
    if [ $a -eq 4 ]; then sed -i "14s/abc/${b}/g" header; sed -i "14s/def/${c}/g" header; fi
  done

seq 1 1 $count > sequence
paste -d'  '  sequence NW.xyz > tmp && mv tmp NW.xyz

cat header NW.xyz > 2.coordinates
rm sequence NW.xyz max min header tmp*

LAMMPSB_EXE=/apps/scratch/compile/new/v100/lammps-7Aug19/src/lmp_v100
mpirun -np 20 $LAMMPSB_EXE -sf gpu -pk gpu 2 -in 1.input > output                                                  #run step
rm coarsening.*

for d in `seq 1 5`
  do
    ### check for reached temperature in log.lammps
    loglines=$(wc -l < log.lammps)
    loglines_32=$(($loglines-32))
    sed -n ${loglines_32}p log.lammps | awk '{print $3}' > tmp
    tmp1=`cat tmp`
    printf '%.*f\n' 0 "$tmp1" > tmp2
    tmp3=`cat tmp2`  # reached temperature

    ###  check for set temperature in 1.input
    sed -n 17p 1.input | awk '{print $4}' > tmp4
    tmp5=`cat tmp4`
    printf '%.*f\n' 0 "$tmp5" > tmp6
    oldT=`cat tmp6`
    echo $oldT # set temperature

    if [ $tmp3 -lt 348 ] || [ $tmp3 -gt 358 ]
    then
      delT=$(($tmp3-353))
      newT=$(($oldT-2*$delT))
      echo $newT
      sed -i "17s/${oldT}/${newT}/g" 1.input
      sed -i "22s/${oldT}/${newT}/g" 1.input
      rm tmp*
      LAMMPSB_EXE=/apps/scratch/compile/new/v100/lammps-7Aug19/src/lmp_v100
      mpirun -np 20 $LAMMPSB_EXE -sf gpu -pk gpu 2 -in 1.input > output
      rm coarsening.*

    else
      rm tmp*
      exit
    fi
  done
