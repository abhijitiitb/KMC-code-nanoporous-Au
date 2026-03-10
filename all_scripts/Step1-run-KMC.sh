cp ../*_eod $(pwd)
eod1=`cat start_eod`
eod2=`cat end_eod`
echo $eod1 $eod2
cp ../18nm_trunc_octa $(pwd)
cp $path/all_scripts/000.KMC/* $(pwd) 
cp kmc.parameters_original kmc.parameters

count=$(wc -l < 18nm_trunc_octa)
echo $count

sed -i "16s/209231/${count}/g" kmc.parameters
number=$(echo $eod1/0.01 | bc)
echo $number
# no to replace atoms
if [ $number -gt 0 ]
then
  sed -i "14s/.TRUE./.FALSE./g" kmc.parameters
fi

#########################################################################################################################
cp 18nm_trunc_octa copy_initial.txt

cat copy_initial.txt | wc -l > total_atoms.txt
total_atoms=$(<total_atoms.txt)
echo $total_atoms

tail -n $total_atoms copy_initial.txt > d.txt
awk '{gsub("1","Fe",$2)}1' d.txt > e.txt
awk '{gsub("2","Pt",$2)}1' e.txt > f.txt

awk '{print $2,$3,$4,$5}' f.txt > i.txt

a=2
z=3
for i in `seq 1 3`
  do
    a=$(($a+1))
    cat d.txt | awk -v a=$a '{print $a}' > tmp$a
    sort -k1n tmp$a > tmp && mv tmp tmp$a
    head -1 tmp$a > min && tail -1 tmp$a > max
    b=`cat min` && c=`cat max`
    if [ $a -eq $z ]; then echo Xmin: $b and Xmax: $c; Xmin=$b Xmax=$c; fi
    if [ $a -eq $(( $z+1 )) ]; then echo Ymin: $b and Ymax: $c; Ymin=$b Ymax=$c; fi
    if [ $a -eq $(( $z+2 )) ]; then echo Zmin: $b and Zmax: $c; Zmin=$b Zmax=$c; fi
  done

echo "$total_atoms $total_atoms $total_atoms $total_atoms 0.00000000E+00 NAtPtms NmPtve PPttentialEnergy" > z.txt
echo "$Xmax $Ymax $Zmax 0.00000000 0.00000000 0.00000000 MD time=0.000E+00 ps" >> z.txt

cat z.txt i.txt > NW.100.xyz

rm  -r *.txt tmp* max min
#########################################################################################################################

cp $path/01.KMC/main_KMC.f90_original $path/01.KMC/main_KMC.f90
sed -i "81s/CALL Initialize_KMC/stop/g" $path/01.KMC/main_KMC.f90
$path/all_scripts/coordtoxyz-to-KMC.sh
$path/01.KMC/amd.par.x > out2_KMC 2>&1
rm  kmc.parameters_analyze *overlap*

mv latchedslab.xyz NW.101.xyz
mv 18nm_trunc_octa coordinates_from_last_MD

cp NW.101.xyz NW.xyz
sed -i "1,2d" NW.xyz
cat NW.xyz | awk '{print $1, $2, $3, $4}' > tmp && mv tmp NW.xyz

sed -i 's/Pt/2/g' NW.xyz
sed -i 's/Fe/1/g' NW.xyz

count1=$(wc -l < NW.xyz)
echo $count1
seq 1 1 $count1 > sequence
paste -d' ' sequence NW.xyz > tmp && mv tmp NW.xyz
mv NW.xyz 18nm_trunc_octa

# calculate the electroactive atoms
cat 18nm_trunc_octa | awk '{print $2}' | awk '{if ($1 == 1) print $1}'> Fe-atoms
b=$(wc -l < Fe-atoms)
echo $b

echo Please enter last eod acheived
eod=$eod1
#read eod
echo $eod
eod=$(echo 1.0-$eod | bc -l)
echo $eod

nEA_atoms=$(echo $b/${eod} | bc) 
echo ${nEA_atoms}
sed -i "17s/156864/${nEA_atoms}/g" kmc.parameters
sed -i "16s/${count}/209231/g" kmc.parameters

count=$(wc -l < 18nm_trunc_octa)
echo $count

sed -i "16s/209231/${count}/g" kmc.parameters

# provide dissolution to be performed
echo Please enter incremental eod to be acheived
eod=$eod2
#read eod           
echo $eod          
                   
cp $path/01.KMC/KMCNanoporous.run.f90_original $path/01.KMC/KMCNanoporous.run.f90
sed -i "213s/0.10/${eod}/g" $path/01.KMC/KMCNanoporous.run.f90
sed -i "244s/0.10/${eod}/g" $path/01.KMC/KMCNanoporous.run.f90
sed -i "245s/0.10/${eod}/g" $path/01.KMC/KMCNanoporous.run.f90
sed -i "255s/0.10/${eod}/g" $path/01.KMC/KMCNanoporous.run.f90

# run KMC
cp $path/01.KMC/main_KMC.f90_original $path/01.KMC/main_KMC.f90
sed -i '78d' $path/01.KMC/main_KMC.f90
sed -i '77d' $path/01.KMC/main_KMC.f90
sed -i '76d' $path/01.KMC/main_KMC.f90
sed -i '75d' $path/01.KMC/main_KMC.f90
$path/all_scripts/coordtoxyz-to-KMC.sh
$path/01.KMC/amd.par.x > out3_KMC 2>&1                                          #run step
rm fort.* *overlap* sequence Fe-atoms
#cp $path/01.KMC/main_KMC.f90_original $path/01.KMC/main_KMC.f90
sed -i "77s/CALL Initialize_KMC/stop/g" $path/01.KMC/main_KMC.f90
sed -i '73d' $path/01.KMC/main_KMC.f90
$path/all_scripts/coordtoxyz-to-KMC.sh

mv NW.1.xyz NW.2.xyz
mv latchedslab.xyz NW.101.xyz

#remove undercordinated atoms and create xyz and 18nm_trunc_octa files###########################################################
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
mv NW.xyz 18nm_trunc_octa

#This code can convert 18nm_trunc_octa file to of type xyz type make necessary changes to code######################################
cp 18nm_trunc_octa copy_initial.txt

cat copy_initial.txt | wc -l > total_atoms.txt
total_atoms=$(<total_atoms.txt)
echo $total_atoms

tail -n $total_atoms copy_initial.txt > d.txt
awk '{gsub("1","Fe",$2)}1' d.txt > e.txt
awk '{gsub("2","Pt",$2)}1' e.txt > f.txt

awk '{print $2,$3,$4,$5}' f.txt > i.txt

a=2
z=3
for i in `seq 1 3`
  do
    a=$(($a+1))
    cat d.txt | awk -v a=$a '{print $a}' > tmp$a
    sort -k1n tmp$a > tmp && mv tmp tmp$a
    head -1 tmp$a > min && tail -1 tmp$a > max
    b=`cat min` && c=`cat max`
    if [ $a -eq $z ]; then echo Xmin: $b and Xmax: $c; Xmin=$b Xmax=$c; fi
    if [ $a -eq $(( $z+1 )) ]; then echo Ymin: $b and Ymax: $c; Ymin=$b Ymax=$c; fi
    if [ $a -eq $(( $z+2 )) ]; then echo Zmin: $b and Zmax: $c; Zmin=$b Zmax=$c; fi
  done

echo "$total_atoms $total_atoms $total_atoms $total_atoms 0.00000000E+00 NAtPtms NmPtve PPttentialEnergy" > z.txt
echo "$Xmax $Ymax $Zmax 0.00000000 0.00000000 0.00000000 MD time=0.000E+00 ps" >> z.txt

cat z.txt i.txt > NW.102.xyz

rm  -r *.txt tmp* max min
#################################################################################################################################
#cp 18nm_trunc_octa $path/18nm_trunc_octa
