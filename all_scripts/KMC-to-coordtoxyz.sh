pwd > $path/01.KMC/path2
cd $path/01.KMC/
cp main_coordtoxyz.f90 main.f90
cp KMCNanoporous.initialize_coordtoxyz.f90 KMCNanoporous.initialize.f90
cp $path/01.KMC/KMCNanoporous.run.f90_original $path/01.KMC/KMCNanoporous.run.f90
touch main.f90
make all

cd `cat path2`
