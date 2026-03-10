pwd > $path/01.KMC/path2
cd $path/01.KMC/
cp main_KMC.f90 main.f90
cp KMCNanoporous.initialize_KMC.f90 KMCNanoporous.initialize.f90
touch main.f90
make all

cd `cat path2`
