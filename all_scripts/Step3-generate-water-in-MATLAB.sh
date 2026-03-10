cp ../step2_MD/dump.lammpstrj $(pwd)
sed '4!d' dump.lammpstrj > count
a=`cat count`
echo $a
tail -$a dump.lammpstrj > 2red.coordinates

rm dump.lammpstrj count

cp $path/all_scripts/MATLAB2/*.m $(pwd)

/usr/local/bin/matlab -r "run ('runn.m');"                             #run step

until [ -f ./metal-solvent-overlap.xyz ]
do
     sleep 50
done

rm metal-solvent-overlap.xyz *.m
