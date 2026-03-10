cp NW.2.xyz NW.xyz
sed -i "1,2d" NW.xyz
sort -k6n NW.xyz > tmp
mv tmp NW.xyz

cat NW.xyz | awk '{print $6}' | sort -k1n | awk '{if ($1 < 3) print $1}'> undercoord
b=$(wc -l < undercoord)

if [ $b -gt 0 ]
then
  sed -i "1,${b}d" NW.xyz
fi
cat NW.xyz | awk '{print $1, $2, $3, $4}' > tmp && mv tmp NW.xyz

sed -i 's/Pt/2/g' NW.xyz
sed -i 's/Fe/1/g' NW.xyz

count=$(wc -l < NW.xyz)
echo $count
seq 1 1 $count > sequence
paste -d' ' sequence NW.xyz > tmp && mv tmp NW.xyz
mv NW.xyz coordinates

cp $path/00.scripts/000.KMC/* $(pwd)
cp kmc.parameters_original kmc.parameters

sed -i "16s/209231/${count}/g" kmc.parameters
sed -i "14s/.TRUE./.FALSE./g" kmc.parameters

$path/01.KMC/amd.par.x
mv latchedslab.xyz NW.102.xyz

rm sequence undercoord kmc.parameters_analyze NW.2.xyz  *over*

