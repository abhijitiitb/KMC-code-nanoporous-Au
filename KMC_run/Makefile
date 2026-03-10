#############################################################
#Select compilation option
#COMPILETYPE=gfortdebug
#COMPILETYPE=gfort
#COMPILETYPE=inteldebug
#COMPILETYPE=intel
#COMPILETYPE=pgdebug
#COMPILETYPE=pg
#############################################################
#AMDMode=serial
AMDMode=parallel
OMP=
#OMP= -fopenmp 
#OMP= -openmp 
#MEMMODEL= -mcmodel=small
#MEMMODEL= -mcmodel=medium
#MEMMODEL= -mcmodel=large
MEMMODEL= 
BaseExecName=amd
#############################################################
src=/home/aditya/repository/
#############################################################
ifeq '${AMDMode}' 'serial'
   AMDFILESTmp=main.o
   ExecutableName=$(BaseExecName).x
   COMPILETYPE=gfort
   PARALLEL_MANAGER=  
   PARALLEL_AMD=  
   PARALLEL_KMC= 
endif
ifeq '${AMDMode}' 'parallel'
   AMDFILESTmp=main.o 
   ExecutableName=$(BaseExecName).par.x
   COMPILETYPE=mpidebug
   PARALLEL_MANAGER= mpi_manager.o 
   PARALLEL_MD= mpi_md_domaindecomposition.o #mpi_md_parrep.o  
   PARALLEL_KMC= kmc_processsearch.par.o  
endif
#############################################################
#--- MPI
ifeq '${COMPILETYPE}' 'mpi'
#  F90=/n/srv/abhijit/mpich2-install-gfort/bin/mpif90
#  F90=/n/srv/abhijit/mpich2-install/bin/mpif90
  F90=mpif90 -compile_info
  F90FLAGS= -O0 -fbounds-check  -fcheck-array-temporaries -fimplicit-none 
  LDFLAGS= 
  #cannot combine g77 and gfortran object files
#  F77=/n/srv/abhijit/mpich2-install-gfort/bin/mpif90
  F77=mpif90
  F77FLAGS= -O0
endif
ifeq '${COMPILETYPE}' 'mpidebug'
#  F90=/n/srv/abhijit/mpich2-install-gfort/bin/mpif90
#  F90=/n/srv/abhijit/mpich2-install/bin/mpif90
  F90=mpif90 
#   F90=/opt/mpich2-amd/bin/mpif90
  F90FLAGS=  -O2 
  LDFLAGS= 
  #cannot combine g77 and gfortran object files
#  F77=/n/srv/abhijit/mpich2-install-gfort/bin/mpif90
  F77=mpif90
#  F77=/opt/mpich2-amd/bin/mpif90
  F77FLAGS= -O3
endif
#--- Gfortran
ifeq '${COMPILETYPE}' 'gfortdebug'
  F90=gfortran 
  F90FLAGS= -O0 -fbounds-check  -Wall
  LDFLAGS= -Bstatic 
  #cannot combine g77 and gfortran object files
  F77=gfortran 
  F77FLAGS= -O0 -fbounds-check
endif
ifeq '${COMPILETYPE}' 'gfort'
  F90=gfortran 
  #F90FLAGS= -O2 -fbounds-check 
  F90FLAGS= -O3 
  #LDFLAGS= -Bstatic
  #LDFLAGS= -pg 
  LDFLAGS=  
  F77=gfortran 
  F77FLAGS= -O3 
endif
ifeq '${COMPILETYPE}' 'gfortcheck'
  F90=gfortran 
  F90FLAGS= -O0 -fbounds-check  
  LDFLAGS= -Bstatic 
  F77=gfortran 
  F77FLAGS= -O0 
endif
ifeq '${COMPILETYPE}' 'inteldebug'
  F90=ifort 
  F90FLAGS= -O0 -check -ftrapuv -warn all -CB 
  LDFLAGS= -Bstatic -debug
  F77=ifort 
  F77FLAGS= -O0 -check bounds -check pointers -ftrapuv 
endif 
ifeq '${COMPILETYPE}' 'intel'
  #Intel fast option
  F90=ifort 
  F90FLAGS= -O3 -Bstatic -inline-level=2 -debug all -check bounds -check pointers
  LDFLAGS= -Bstatic 
  F77=ifort 
  F77FLAGS= -O3 -Bstatic
  ##F90FLAGS=-O2 -ipa 
  ##LDFLAGS= -Bstatic 
  ##F90FLAGS= -O2 -Bstatic -tpp7 
endif
ifeq '${COMPILETYPE}' 'pgdebug'
  #for profiler turn PROFILE option on
  #PROFILE=-Mprof=lines
  PROFILE= 
  F90=pgf90 
  F77=pgf77 
  F90FLAGS= -O0 -C -Mchkfpstk -Mchkptr -Mchkstk -Kieee -I../pgapack/include 
  F77FLAGS= -O0 -C -Mchkfpstk -Mchkstk -Kieee -I../pgapack/include 
  LDFLAGS= -Bstatic  
  #LDFLAGS=$(F90FLAGS) $(PROFILE) -Bstatic  
endif
ifeq '${COMPILETYPE}' 'pg'
  #PG optimized
  #PROFILE=-Mprof=func
  PROFILE= 
  F90=pgf90
  F77=pgf77
  #F90FLAGS=-tp=amd64 -fast $(PROFILE) -Minline $(CODEFLAGS) -Mipa
  #F90FLAGS= -tp=p7 -fast $(PROFILE) -Minline $(CODEFLAGS)
  F90FLAGS= -fast $(PROFILE) -Minline $(CODEFLAGS)
  F77FLAGS= -fast $(PROFILE) -Minline $(CODEFLAGS)
  LDFLAGS=$(F90FLAGS) $(PROFILE) -Bstatic 
endif 
#############################################################
COMPILE= -c
COMPILEOMP= $(OMP) -c
EXEC= $(OMP) $(MEMMODEL) -o
#INCLUDES=  -I/usr/lib/ -llapack
#INCLUDES=  -L/opt/Libraries/lapack-3.4.1/ -llapack $(src)BLAS/blas_gfortran.a 
##INCLUDES=  /opt/libraries/LAPACK/lapack_LINUX.a $(src)BLAS/blas_gfortran.a 
INCLUDES= 
#INCLUDES=  -L/opt/intel/lib/intel64
#############################################################
###Files not in use
### trpfunctions.o
%.o: %.mod

##simulatedannealing.o 
AMDCOMMONFILES= amd_typedef.o nrprocedures.o kmc_typdef.o amd_typedef.lvl2.o \
$(PARALLEL_MANAGER)$ randomgen.o amd_db_manipulate.o \
amd_utilities.o io.o check.o \
routines.f OptimPackage.o GA.o  neighborlist.o \
EAMGenerate.o BuckinghamGenerate.o TBSMAGenerate.o StillingerWeberGenerate.o \
CoulombInitialize.o stillingerweber_initialize.o tersoff_initialize.o \
EAMInitialize.o MEAMInitialize.o LJInitialize.o BuckinghamInitialize.o \
Coulomb.o stillingerweber.o tersoff.o EAM.o  \
LJ.o Buckingham.o DFT.o InteratomicPotentials.o patternrecognition.o \
OptimizationAL.o OptimizationNEB.o OptimizationDimer.o \
species.o amd_crystal.o InteratomicPotentialsAssess.o \
CELattice.o\
Dimer.o NEB.o nmodes.o metropolis.o  md.o $(PARALLEL_MD)$ \
kmc_utilities.o kmc_db_manipulate.o kmc_configurations.o \
kmc_relaxlattice.o  kmcneb_processsearch.o kmc_reporterror.o kmc_openfile.o \
kmc_backuprecord.o statistics_spatialaverages.o statistics.o \
machinelearning.o ANN.o Latticekmc2dmodule.o \
LEKMC2ANN.o amd_db_manipulate.lvl2.o tpmd.o \
$(PARALLEL_KMC)$  kmc_generateprocesslistDimer.o kmc_generateprocesslistAMD.o lekmcmd.o \
kmc_generateprocesslistMD.o kmc_generateprocesslistUser.o \
kmc_generateprocesslist.o kmc_initialize.o KMCBinarySearch.o\
kmc_algorithm_updatesystem.domaindecompose.o kmc_algorithm_updatesystem.carve.o kmc_algorithm.o \
alignnanoparticle1.o $(AMDFILESTmp)$
############################################################
BLASFILES=caxpy.o  cher2k.o  csyrk.o   ddot.o    dspmv.o   dtpsv.o   snrm2.o   ssyr2.o   zaxpy.o   zhemv.o   zsyrk.o \
ccopy.o  cher.o    ctbmv.o   dgbmv.o   dspr2.o   dtrmm.o   sasum.o   srot.o    ssyr2k.o  zcopy.o   zher2.o   ztbmv.o \
cdotc.o  cherk.o   ctbsv.o   dgemm.o   dspr.o    dtrmv.o   saxpy.o   srotg.o   ssyr.o    zdotc.o   zher2k.o  ztbsv.o \
cdotu.o  chpmv.o   ctpmv.o   dgemv.o   dswap.o   dtrsm.o   scabs1.o  srotm.o   ssyrk.o   zdotu.o   zher.o    ztpmv.o \
cgbmv.o  chpr2.o   ctpsv.o   dger.o    dsymm.o   dtrsv.o   scasum.o  srotmg.o  stbmv.o   zdrot.o   zherk.o   ztpsv.o \
cgemm.o  chpr.o    ctrmm.o   dnrm2.o   dsymv.o   dzasum.o  scnrm2.o  ssbmv.o   stbsv.o   zdscal.o  zhpmv.o   ztrmm.o \
cgemv.o  crotg.o   ctrmv.o   drot.o    dsyr2.o   dznrm2.o  scopy.o   sscal.o   stpmv.o   zgbmv.o   zhpr2.o   ztrmv.o \
cgerc.o  cscal.o   ctrsm.o   drotg.o   dsyr2k.o  icamax.o  sdot.o    sspmv.o   stpsv.o   zgemm.o   zhpr.o    ztrsm.o \
cgeru.o  csrot.o   ctrsv.o   drotm.o   dsyr.o    idamax.o  sdsdot.o  sspr2.o   strmm.o   zgemv.o   zrotg.o   ztrsv.o \
chbmv.o  csscal.o  dasum.o   drotmg.o  dsyrk.o   isamax.o  sgbmv.o   sspr.o    strmv.o   zgerc.o   zscal.o \
chemm.o  cswap.o   daxpy.o   dsbmv.o   dtbmv.o   izamax.o  sgemm.o   sswap.o   strsm.o   zgeru.o   zswap.o \
chemv.o  csymm.o   dcabs1.o  dscal.o   dtbsv.o   lsame.o   sgemv.o   ssymm.o   strsv.o   zhbmv.o   zsymm.o \
cher2.o  csyr2k.o  dcopy.o   dsdot.o   dtpmv.o   sger.o    ssymv.o   xerbla.o  zhemm.o   zsyr2k.o

############################################################
LAPACKFILES=disnan.o  dlaev2.o    dlansp.o  dlapy2.o  dlarfg.o  dlascl.o \
 dlasr.o   dlassq.o  dorg2l.o  dspev.o   dsteqr.o \
dlae2.o   dlaisnan.o  dlanst.o  dlarf.o   dlartg.o  \
dlaset.o  dlasrt.o  dopgtr.o  dorg2r.o  dsptrd.o  dsterf.o \
dlamch.o  iladlc.o  iladlr.o  lsame.o  xerbla.o
############################################################
#COMMONOBJS = $(AMDFILES) $(AMDCOMMONFILES) $(TADFILES)
COMMONOBJS =  $(AMDCOMMONFILES) 

SRCS := $(wildcard *.f*)
OBJS := $(SRCS:%.f*=%.o)

############################################################
# Compilation section

vpath %.F90 ./:$(src)
vpath %.f90 ./:$(src)
vpath %.f ./:$(src)



all : $(ExecutableName)

#intel : $(ExecutableName)
#	F90=/opt/mpich2-intel/bin/mpif90
#	F77=/opt/mpich2-intel/bin/mpif90
#	INCLUDES=  -L/opt/intel/lib/intel64

#amd : $(ExecutableName)
#	F90=/opt/mpich2-amd/bin/mpif90
#	F77=/opt/mpich2-amd/bin/mpif90
#	INCLUDES=  -L/opt/Libraries/lapack-3.4.1/ -llapack  on garuda
#	INCLUDES=  -L/opt/libraries/LAPACK/ -llapack_LINUX  on serengeti

############################################################


%.o : %.f90
	$(F90) $(COMPILE) $(F90FLAGS) $<

# Compilation section
%.o : %.F90
	$(F90) $(COMPILEOMP) $(F90FLAGS) $<

# Compilation section
%.o : %.f
#	$(F90) $(COMPILE) $(F90FLAGS) $<
	$(F77) $(COMPILE) $(F77FLAGS) $<

$(ExecutableName) : $(COMMONOBJS)
	$(F90) $(LDFLAGS) $(EXEC) $(ExecutableName) $^ $(INCLUDES)

clean :
	rm -f $(src)*.o $(src)*.mod $(src)*.obj *.o *.mod *.obj $(BaseExecName).* 

cleano :
	rm -f *.o *.obj
#############################################################
