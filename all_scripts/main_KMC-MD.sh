###################################################################
###********************** Step1: perform KMC ****************** ###
###################################################################
mkdir step1_KMC
cd step1_KMC
   cp $path/all_scripts/step1-run-KMC.sh $(pwd)
   ./007.run-KMC.sh > out_step1_KMC.txt   ## run KMC after latching the NP to the desired EOD
cd ..

###################################################################
###************** Step2: perform Energy minimization ***********###
###################################################################

mkdir step2_MD
cd step2_MD
   cp $path/all_scripts/step2-do-energy-minimization-LAMMPS.sh $(pwd)
   ./008.do-energy-minimization-LAMMPS.sh > out_step2_MD.txt   ## Do energy minimization so that NP stabilizes
cd ..

###################################################################
###*************** Step3: create water molecules ***************###
###################################################################

mkdir step3_MD
cd step3_MD
   cp $path/all_scripts/step3-generate-water-in-MATLAB.sh $(pwd)
   ./009.generate-water-in-MATLAB.sh > out_step3_MD.txt
cd ..

###################################################################
##************** Step4: randomize water molecules **************###
###################################################################

mkdir step4_MD
cd step4_MD
   cp $path/all_scripts/step4-randomize-water-LAMMPS.sh $(pwd)
   ./010.randomize-water-LAMMPS.sh > out_step4_MD.txt
cd .. 

###################################################################
##************** Step5: run molecular dynamics *****************###
###################################################################

mkdir step5_MD
cd step5_MD
   cp $path/all_scripts/step5-run-MD.sh $(pwd)
   ./011.run-MD.sh > out_step5_MD.txt
cd ..

###################################################################
##************************** END *******************************###
###################################################################
