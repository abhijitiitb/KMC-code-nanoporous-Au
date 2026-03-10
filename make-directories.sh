echo $(date)           # to print date of execution
echo $(pwd) > a.txt
path=`cat a.txt`
#echo $path
export path            # "path" is a variable that has the folder path from where the simulation was run
rm *.txt

b=.00
e=.20   # increment iEOD

for i in `seq 1 4`
  do
    a=$(echo $b+$e | bc -l)

    mkdir 0$b-0$a
    cd 0$b-0$a
       echo $e > end_eod
       echo $b > start_eod
       cp ../18nm_trunc_octa $(pwd)     
       cp $path/00.scripts/main-KMC-MD.sh $(pwd)
       ./main-KMC-MD.sh
    cd ..
    b=$(echo $a+0.00 | bc -l)
  done

b=0.80 
f=0.10 
mkdir 0.80-0.90
cd 0.80-0.90
   echo $f > end_eod
   echo $b > start_eod
   cp ../18nm_trunc_octa $(pwd)     
   cp $path/00.scripts/main-KMC-MD.sh $(pwd)
   ./main-KMC-MD.sh
cd ..

b=0.90
f=0.05
mkdir 0.90-0.95
cd 0.90-0.95
   echo $f > end_eod
   echo $b > start_eod
   cp ../18nm_trunc_octa $(pwd)     
       cp $path/00.scripts/main-KMC-MD.sh $(pwd)
       ./main-KMC-MD.sh
cd ..

b=0.95
f=0.025
mkdir 0.95-0.975
cd 0.95-0.975
   echo $f > end_eod
   echo $b > start_eod
   cp ../18nm_trunc_octa $(pwd)     
       cp $path/00.scripts/step1-run-KMC.sh $(pwd)
       ./step1-run-KMC.sh
cd ..
