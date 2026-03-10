MODULE LatticeKMCModule
!perform KMC - sequential or binary search
!perform NEB for hop and exchange processes
!generate CEKMC models by accessing ANN module
!deposition enabled
   USE ANNModule
   USE OptimizationNEB
   USE VARIABLE_TYPE
   USE NEB_package
   USE Ecuyer_random
   IMPLICIT NONE
   INTEGER,PARAMETER, PRIVATE :: m=10,n=10,natoms=28
   DOUBLE PRECISION, PARAMETER :: F=0.022  ! Deposition rate in monolayers/sec
   DOUBLE PRECISION, PRIVATE :: prefactor=1.e13,time=0.
   DOUBLE PRECISION, DIMENSION(n,m), PRIVATE :: gridn,grids,gride,gridw
   INTEGER, DIMENSION(n,m), PRIVATE :: grid
   DOUBLE PRECISION, DIMENSION(n*m*4+1) :: gridrate
   INTEGER, PRIVATE :: nkmc=10000
   !DOUBLE PRECISION,DIMENSION(1024), PRIVATE :: ActivationEnergy,ETransitionState,EInitialState
   DOUBLE PRECISION, PRIVATE :: ktn,kts,kte,ktw,kt,r,kr,sumrates,coverage
   TYPE(SystemContainer), POINTER, PRIVATE :: AL=>NULL(),AL2=>NULL()
   CHARACTER(len=100), PRIVATE :: filename
   DOUBLE PRECISION, PARAMETER, PRIVATE :: T=300._dp,kB=8.617343e-5_dp
   INTEGER, PRIVATE :: NumberNEBCalls=0
   !INTEGER,DIMENSION(1024) :: IntegerKey
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   
   SUBROUTINE ANNKMC
   !gfortran -o TestANNKMC.x kmc2d.f90 kmc2dmodule.f90
   !Purpose: Find the evolution of a surface containing Ag (Ag/Ag(100) system) using KMC method
   !where the process rates are computed initially usinng NEB method and then using an ANN model
   !This code shows how the ANN model can be constructed on-the-fly
      IMPLICIT NONE
      INTEGER :: ikmc,i,j,dir
      CHARACTER :: crapchar
      REAL :: CPUTime
      
      NumberNEBCalls=0
      
      CALL InitializeKMC2d()
      !CALL PrintLattice()
      CALL SetupAllRates()
      !CALL PrintLattice1()
      
      OPEN(UNIT=540,FILE="ANNKMCSummary")
      WRITE(540,*) "ikmc Time CPUTime NData NDataStored NinteractionsFitted"
      OPEN(UNIT=541,FILE="ANNKMCLattice.m")
      DO ikmc=1,10000
         IF (MOD(ikmc,10000)==0) CALL SetupAllRates()
         CALL SelectKMCProcess(i,j,dir)
         CALL LocalUpdate(i,j,dir)
         IF (MOD(ikmc,100)==0) THEN
            CALL CPU_TIME(CPUTime)
            WRITE(540,*) ikmc,time,CPUTime,ANNList%ANN%ndata,ANNList%ANN%ndata_stored,ninteractionsfitted(ANNList%ANN)
         END IF
         IF (MOD(ikmc,1000)==0) THEN
            WRITE(541,*) "ikmc=[ikmc ",ikmc,";"
            WRITE(541,*) "time=[time ",time,";"
            CALL PrintLattice(541)
         END IF
      END DO
      CLOSE(540)
      CLOSE(541)
      stop
      WRITE(*,*) "_____________________"
      
   
   END SUBROUTINE ANNKMC
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   
   SUBROUTINE SetupAllRates()
      IMPLICIT NONE
      
      INTEGER :: i,j,dir
      INTEGER, DIMENSION(natoms) :: sigma
! kB is Boltzmann constant T is absolute temperature   
      
      DO i=1,n
         DO j=1,m
            gridn(i,j)=0.
            grids(i,j)=0.
            gride(i,j)=0.
            gridw(i,j)=0.
            IF (grid(i,j)==0) CYCLE
            
            IF(grid(iPBCLattice(i-1,n),j)==0) THEN
               dir=1
               CALL CreatesigmaVector(i,j,dir,sigma)
               gridn(i,j)=prefactor*(exp(-GetActivationBarrier(sigma)/(kB*T)))
            END IF
            IF(grid(iPBCLattice(i+1,n),j)==0) THEN 
               dir=2
               CALL CreatesigmaVector(i,j,dir,sigma)
               grids(i,j)=prefactor*(exp(-GetActivationBarrier(sigma)/(kB*T)))
            END IF
            IF(grid(i,iPBCLattice(j+1,m))==0) THEN  
               dir=3
               CALL CreatesigmaVector(i,j,dir,sigma)
               gride(i,j)=prefactor*(exp(-GetActivationBarrier(sigma)/(kB*T)))
            END IF
            IF(grid(i,iPBCLattice(j-1,m))==0) THEN
               dir=4
               CALL CreatesigmaVector(i,j,dir,sigma) 
               gridw(i,j)=prefactor*(exp(-GetActivationBarrier(sigma)/(kB*T))) 
            END IF
         
         END DO
      END DO
      
      !ktn=SUM(SUM(gridn,1),2) !total rate to north
      !kts=SUM(SUM(grids,1),2) !total rate to south
      !kte=SUM(SUM(gride,1),2) !total rate to east
      !ktw=SUM(SUM(gridw,1),2) !total rate to west
      
      WRITE(6,*) "LOOK AT LINE 110 OF LATTICEKMC2DMODULE.F90 --"
      WRITE(6,*) "ktn=SUM(SUM(gridn,1),2) !tota ..."
      WRITE(6,*) " LINES 110-113 NEED TO BE UNCOMMENTED"
      STOP
      
      kt=ktn+kts+kte+ktw
   END SUBROUTINE SetupAllRates
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ninteractionsfitted(ANN)
   !finds how many interactions have been fitted in ANN
      IMPLICIT NONE
      TYPE(NeuralNet), POINTER :: ANN
      INTEGER :: ninteractionsfitted,i,j
      LOGICAL, DIMENSION(:,:), POINTER :: CorrelationNotConsidered
      
      ninteractionsfitted=0
      CorrelationNotConsidered=>ANN%CorrelationNotConsidered
      DO j=1,ANN%natoms
         DO i=1,j
            IF (CorrelationNotConsidered(i,j)) ninteractionsfitted=ninteractionsfitted+1
         END DO
      END DO
   END FUNCTION ninteractionsfitted
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LocalUpdate(i,j,dir)
      IMPLICIT NONE
      INTEGER :: i,j,dir,i1,j1,i2,j2,i3,j3
      
      !Final matrix after the process
      grid(i,j)=0
      !find the affect neighbor
      i1=i
      j1=j
      i2=i
      j2=j
      SELECT CASE (dir)
      CASE(1) !north
         i1=iPBCLattice(i-1,n)
         i2=i-1
      CASE(2) !south
         i1=iPBCLattice(i+1,n)
         i2=i+1
      CASE(3) !east
         j1=iPBCLattice(j+1,m)
         j2=j+1
      CASE(4) !west
         j1=iPBCLattice(j-1,m)
         j2=j-1
      END SELECT
      grid(i1,j1)=1
      
      !affected region of i,j and i1,j1
      IF (natoms==10) THEN
         DO i3=MIN(i-2,i2-2),MAX(i+2,i2+2)
            DO j3=MIN(j-2,j2-2),MAX(j+2,j2+2)
               CALL LocalUpdateSite(iPBCLattice(i3,n),iPBCLattice(j3,m))
            END DO
         END DO
      ELSEIF (natoms==28) THEN
         DO i3=MIN(i-3,i2-3),MAX(i+3,i2+3)
            DO j3=MIN(j-3,j2-3),MAX(j+3,j2+3)
               CALL LocalUpdateSite(iPBCLattice(i3,n),iPBCLattice(j3,m))
            END DO
         END DO
      ELSE
         WRITE(6,*) "This local update has not been implemented"
         STOP
      END IF
      kt=ktn+kts+kte+ktw
   END SUBROUTINE LocalUpdate
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   
   SUBROUTINE LocalUpdateSite(i,j)
      IMPLICIT NONE
      INTEGER :: i,j,dir
      INTEGER, DIMENSION(natoms) :: sigma
      
      ktn=ktn-gridn(i,j)
      kts=kts-grids(i,j)
      kte=kte-gride(i,j)
      ktw=ktw-gridw(i,j)
      
      gridn(i,j)=0.
      grids(i,j)=0.
      gride(i,j)=0.
      gridw(i,j)=0.
            
      IF (grid(i,j)==0) RETURN
      
      IF(grid(iPBCLattice(i-1,n),j)==0) THEN
         dir=1
         CALL CreatesigmaVector(i,j,dir,sigma)
         gridn(i,j)=prefactor*(exp(-GetActivationBarrier(sigma)/(kB*T)))
      END IF
      IF(grid(iPBCLattice(i+1,n),j)==0) THEN 
         dir=2
         CALL CreatesigmaVector(i,j,dir,sigma)
         grids(i,j)=prefactor*(exp(-GetActivationBarrier(sigma)/(kB*T)))
      END IF
      IF(grid(i,iPBCLattice(j+1,m))==0) THEN  
         dir=3
         CALL CreatesigmaVector(i,j,dir,sigma)
         gride(i,j)=prefactor*(exp(-GetActivationBarrier(sigma)/(kB*T)))
      END IF
      IF(grid(i,iPBCLattice(j-1,m))==0) THEN
         dir=4
         CALL CreatesigmaVector(i,j,dir,sigma) 
         gridw(i,j)=prefactor*(exp(-GetActivationBarrier(sigma)/(kB*T))) 
      END IF
      
      ktn=ktn+gridn(i,j)
      kts=kts+grids(i,j)
      kte=kte+gride(i,j)
      ktw=ktw+gridw(i,j)
   END SUBROUTINE LocalUpdateSite
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   
   SUBROUTINE InitializeKMC2d
      IMPLICIT NONE
      INTEGER :: i,j
      
      time=0.
      !CALL ReadActivationBarriers()
      
      OPEN(UNIT=501,FILE="InputCoverage")
      READ(501,*) coverage
      CLOSE(501)
      
      DO i=1,n
         DO j=1,m 
         IF(taus88()<coverage) THEN
            grid(i,j)=1
         ELSE
            grid(i,j)=0
         END IF
         END DO
      END DO
            
   END SUBROUTINE InitializeKMC2d
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    FUNCTION GetActivationBarrier(Sigma)
      IMPLICIT NONE
      TYPE(NeuralNet), POINTER :: ANN
      INTEGER, DIMENSION(natoms) :: Sigma
      REAL(dp), DIMENSION(natoms) :: SigmaReal
      REAL(dp), DIMENSION(natoms+(natoms*(natoms-1))/2) :: SigmaRealCorr
      DOUBLE PRECISION :: GetActivationBarrier,e
      INTEGER :: i,j,pos,n1,errorstatus,ndimensionsx
      CHARACTER :: crapchar
     
      GetActivationBarrier=0. !in eV
      n1=ConvertSigma2Number(sigma)
      SigmaReal=REAL(Sigma,dp)
      SigmaRealCorr(1:natoms)=SigmaReal
      pos=natoms
      DO i=1,natoms-1
         DO j=i+1,natoms
            pos=pos+1
            SigmaRealCorr(pos)=SigmaReal(i)*SigmaReal(j)
         END DO
      END DO
      ANN=>ANNList%ANN
      ndimensionsx=natoms+(natoms*(natoms-1))/2
      CALL ANNPrediction(ANN,SigmaRealCorr,ndimensionsx,e,signature=n1+1,errorstatus=errorstatus,iprint=0)
      !WRITE(6,*) "Errorstatus from ANNPrediction:",errorstatus
      IF (errorstatus==0) THEN
         GetActivationBarrier=e
!         WRITE(6,*) "Prediction vs true value:",GetActivationBarrier,ActivationEnergy(n1+1)
!IF (ABS(GetActivationBarrier-ActivationEnergy(n1+1))>0.03_dp) THEN
!CALL DisplayCorrelations(ANNList%ANN)
!WRITE(6,*) ANN%ndegreesmin,ANN%ndegreesmax
!WRITE(UNIT=6,FMT='("Type any key to continue:")',ADVANCE="NO")
!READ(*,*) crapchar
!END IF
      ELSE
         !GetActivationBarrier=ActivationEnergy(n1+1)
         GetActivationBarrier=NEBLatticeKMCActivationEnergy(sigma)
         CALL StoreDataInANN(ANN,SigmaRealCorr,GetActivationBarrier,ndimensionsx,iprint=0,signature=n1+1)
      END IF
   END FUNCTION GetActivationBarrier
 
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ConvertSigma2Number(sigma)
      IMPLICIT NONE
      INTEGER :: i,sigma(natoms),n1,ConvertSigma2Number
      n1=0
      DO i=1,natoms
         n1=n1+((2**(i-1))*(sigma(natoms+1-i)))
      END DO
      ConvertSigma2Number=n1
   END FUNCTION ConvertSigma2Number
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!    SUBROUTINE ReadActivationBarriers()
!      IMPLICIT NONE
      
!      INTEGER :: i,j1(4),sigma(10),j2,n1
!      DOUBLE PRECISION :: Eact,r1,ETS,EIni
      
!      OPEN(UNIT=315,FILE="AdatomHopEnergies.txt")  
!      DO i=1,1024
!         READ(315,*) j1(1:4),Eact,j2,r1,sigma(1:10),EIni,ETS
!         n1=ConvertSigma2Number(sigma)
!         !IntegerKey(i)=n1
!         ActivationEnergy(n1+1)=Eact
!         ETransitionState(n1+1)=ETS
!         EInitialState(n1+1)=EIni
!      END DO
!      CLOSE(315)
!   END SUBROUTINE ReadActivationBarriers
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreatesigmaVector(i,j,dir,sigma)
      IMPLICIT NONE 
 
      INTEGER, DIMENSION(:) :: sigma
      INTEGER, INTENT(IN) :: i,j,dir
      
      SELECT CASE (natoms)
      CASE(10); CALL CreatesigmaVector10(i,j,dir,sigma)
      CASE(28); CALL CreatesigmaVector28(i,j,dir,sigma)
      CASE DEFAULT
         WRITE(6,*) "Err>> No sigma definition available when natoms is ",natoms
         STOP
      END SELECT
   END SUBROUTINE CreatesigmaVector
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreatesigmaVector10(i,j,dir,sigma)
 
      IMPLICIT NONE 
 
      INTEGER, DIMENSION(:) :: sigma
      INTEGER :: i,j,dir
 
      ! 1 north; 2 south; 3 east; 4 west
      
      IF(dir==1) THEN 
         sigma(1)=grid(i,iPBCLattice(j-1,m))
         sigma(2)=grid(iPBCLattice(i-1,n),iPBCLattice(j-1,m))
         sigma(3)=grid(iPBCLattice(i-2,n),iPBCLattice(j-1,m))
         sigma(4)=grid(iPBCLattice(i-2,n),j)
         sigma(5)=grid(iPBCLattice(i-2,n),iPBCLattice(j+1,m))
         sigma(6)=grid(iPBCLattice(i-1,n),iPBCLattice(j+1,m))
         sigma(7)=grid(i,iPBCLattice(j+1,m))
         sigma(8)=grid(iPBCLattice(i+1,n),iPBCLattice(j+1,m))
         sigma(9)=grid(iPBCLattice(i+1,n),j)
         sigma(10)=grid(iPBCLattice(i+1,n),iPBCLattice(j-1,m)) 
      ELSE IF(dir==2) THEN
         sigma(1)=grid(i,iPBCLattice(j+1,m)) 
         sigma(2)=grid(iPBCLattice(i+1,n),iPBCLattice(j+1,m))
         sigma(3)=grid(iPBCLattice(i+2,n),iPBCLattice(j+1,m))
         sigma(4)=grid(iPBCLattice(i+2,n),j)
         sigma(5)=grid(iPBCLattice(i+2,n),iPBCLattice(j-1,m))
         sigma(6)=grid(iPBCLattice(i+1,n),iPBCLattice(j-1,m))
         sigma(7)=grid(i,iPBCLattice(j-1,m))
         sigma(8)=grid(iPBCLattice(i-1,n),iPBCLattice(j-1,m))
         sigma(9)=grid(iPBCLattice(i-1,n),j)
         sigma(10)=grid(iPBCLattice(i-1,n),iPBCLattice(j+1,m))
      ELSE IF (dir==3) THEN
         sigma(1)=grid(iPBCLattice(i-1,n),j)
         sigma(2)=grid(iPBCLattice(i-1,n),iPBCLattice(j+1,m))
         sigma(3)=grid(iPBCLattice(i-2,n),iPBCLattice(j+2,m))
         sigma(4)=grid(i,iPBCLattice(j+2,m))
         sigma(5)=grid(iPBCLattice(i+1,n),iPBCLattice(j+2,m))
         sigma(6)=grid(iPBCLattice(i+1,n),iPBCLattice(j+1,m))
         sigma(7)=grid(iPBCLattice(i+1,n),j)
         sigma(8)=grid(iPBCLattice(i+1,n),iPBCLattice(j-1,m))
         sigma(9)=grid(i,iPBCLattice(j-1,m))
         sigma(10)=grid(iPBCLattice(i-1,n),iPBCLattice(j-1,m)) 
      ELSE IF (dir==4) THEN
         sigma(1)=grid(iPBCLattice(i+1,n),j)
         sigma(2)=grid(iPBCLattice(i+1,n),iPBCLattice(j-1,m))
         sigma(3)=grid(iPBCLattice(i+1,n),iPBCLattice(j-2,m))
         sigma(4)=grid(i,iPBCLattice(j-2,m))
         sigma(5)=grid(iPBCLattice(i-1,n),iPBCLattice(j-2,m))
         sigma(6)=grid(iPBCLattice(i-1,n),iPBCLattice(j-1,m))
         sigma(7)=grid(iPBCLattice(i-1,n),j)
         sigma(8)=grid(iPBCLattice(i-1,n),iPBCLattice(j+1,m))
         sigma(9)=grid(i,iPBCLattice(j+1,m))
         sigma(10)=grid(iPBCLattice(i+1,n),iPBCLattice(j+1,m))
      END IF      
            
   END SUBROUTINE CreatesigmaVector10
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    
   SUBROUTINE CreatesigmaVector28(i,j,dir,sigma)
 
      IMPLICIT NONE 
 
      INTEGER, DIMENSION(:) :: sigma
      INTEGER :: i,j,dir
 
! 1 north; 2 south; 3 east; 4 west  
      IF(dir==1) THEN 
         sigma(1)=grid(iPBCLattice(i+2,n),iPBCLattice(j-2,m))
         sigma(2)=grid(iPBCLattice(i+1,n),iPBCLattice(j-2,m))
         sigma(3)=grid(iPBCLattice(i,n),iPBCLattice(j-2,m))
         sigma(4)=grid(iPBCLattice(i-1,n),iPBCLattice(j-2,m))
         sigma(5)=grid(iPBCLattice(i-2,n),iPBCLattice(j-2,m))
         sigma(6)=grid(iPBCLattice(i-3,n),iPBCLattice(j-2,m))
         sigma(7)=grid(iPBCLattice(i+2,n),iPBCLattice(j-1,m))
         sigma(8)=grid(iPBCLattice(i+1,n),iPBCLattice(j-1,m))
         sigma(9)=grid(iPBCLattice(i,n),iPBCLattice(j-1,m))
         sigma(10)=grid(iPBCLattice(i-1,n),iPBCLattice(j-1,m))
         sigma(11)=grid(iPBCLattice(i-2,n),iPBCLattice(j-1,m))
         sigma(12)=grid(iPBCLattice(i-3,n),iPBCLattice(j-1,m))
         sigma(13)=grid(iPBCLattice(i+2,n),iPBCLattice(j,m))
         sigma(14)=grid(iPBCLattice(i+1,n),iPBCLattice(j,m))
         sigma(15)=grid(iPBCLattice(i-2,n),iPBCLattice(j,m))
         sigma(16)=grid(iPBCLattice(i-3,n),iPBCLattice(j,m))
         sigma(17)=grid(iPBCLattice(i+2,n),iPBCLattice(j+1,m))
         sigma(18)=grid(iPBCLattice(i+1,n),iPBCLattice(j+1,m))
         sigma(19)=grid(iPBCLattice(i,n),iPBCLattice(j+1,m))
         sigma(20)=grid(iPBCLattice(i-1,n),iPBCLattice(j+1,m))
         sigma(21)=grid(iPBCLattice(i-2,n),iPBCLattice(j+1,m))
         sigma(22)=grid(iPBCLattice(i-3,n),iPBCLattice(j+1,m))
         sigma(23)=grid(iPBCLattice(i+2,n),iPBCLattice(j+2,m))
         sigma(24)=grid(iPBCLattice(i+1,n),iPBCLattice(j+2,m))
         sigma(25)=grid(iPBCLattice(i,n),iPBCLattice(j+2,m))
         sigma(26)=grid(iPBCLattice(i-1,n),iPBCLattice(j+2,m))
         sigma(27)=grid(iPBCLattice(i-2,n),iPBCLattice(j+2,m))
         sigma(28)=grid(iPBCLattice(i-3,n),iPBCLattice(j+2,m)) 
      ELSE IF(dir==2) THEN
         sigma(1)=grid(iPBCLattice(i-2,n),iPBCLattice(j+2,m))
         sigma(2)=grid(iPBCLattice(i-1,n),iPBCLattice(j+2,m))
         sigma(3)=grid(iPBCLattice(i,n),iPBCLattice(j+2,m))
         sigma(4)=grid(iPBCLattice(i+1,n),iPBCLattice(j+2,m))
         sigma(5)=grid(iPBCLattice(i+2,n),iPBCLattice(j+2,m))
         sigma(6)=grid(iPBCLattice(i+3,n),iPBCLattice(j+2,m))
         sigma(7)=grid(iPBCLattice(i-2,n),iPBCLattice(j+1,m))
         sigma(8)=grid(iPBCLattice(i-1,n),iPBCLattice(j+1,m))
         sigma(9)=grid(iPBCLattice(i,n),iPBCLattice(j+1,m))
         sigma(10)=grid(iPBCLattice(i+1,n),iPBCLattice(j+1,m))
         sigma(11)=grid(iPBCLattice(i+2,n),iPBCLattice(j+1,m))
         sigma(12)=grid(iPBCLattice(i+3,n),iPBCLattice(j+1,m))
         sigma(13)=grid(iPBCLattice(i-2,n),iPBCLattice(j,m))
         sigma(14)=grid(iPBCLattice(i-1,n),iPBCLattice(j,m))
         sigma(15)=grid(iPBCLattice(i+2,n),iPBCLattice(j,m))
         sigma(16)=grid(iPBCLattice(i+3,n),iPBCLattice(j,m))
         sigma(17)=grid(iPBCLattice(i-2,n),iPBCLattice(j-1,m))
         sigma(18)=grid(iPBCLattice(i-1,n),iPBCLattice(j-1,m))
         sigma(19)=grid(iPBCLattice(i,n),iPBCLattice(j-1,m))
         sigma(20)=grid(iPBCLattice(i+1,n),iPBCLattice(j-1,m))
         sigma(21)=grid(iPBCLattice(i+2,n),iPBCLattice(j-1,m))
         sigma(22)=grid(iPBCLattice(i+3,n),iPBCLattice(j-1,m))
         sigma(23)=grid(iPBCLattice(i-2,n),iPBCLattice(j-2,m))
         sigma(24)=grid(iPBCLattice(i-1,n),iPBCLattice(j-2,m))
         sigma(25)=grid(iPBCLattice(i,n),iPBCLattice(j-2,m))
         sigma(26)=grid(iPBCLattice(i+1,n),iPBCLattice(j-2,m))
         sigma(27)=grid(iPBCLattice(i+2,n),iPBCLattice(j-2,m))
         sigma(28)=grid(iPBCLattice(i+3,n),iPBCLattice(j-2,m)) 
      ELSE IF (dir==3) THEN
         sigma(1)=grid(iPBCLattice(i-2,n),iPBCLattice(j-2,m))
         sigma(2)=grid(iPBCLattice(i-2,n),iPBCLattice(j-1,m))
         sigma(3)=grid(iPBCLattice(i-2,n),iPBCLattice(j,m))
         sigma(4)=grid(iPBCLattice(i-2,n),iPBCLattice(j+1,m))
         sigma(5)=grid(iPBCLattice(i-2,n),iPBCLattice(j+2,m))
         sigma(6)=grid(iPBCLattice(i-2,n),iPBCLattice(j+3,m))
         sigma(7)=grid(iPBCLattice(i-1,n),iPBCLattice(j-2,m))
         sigma(8)=grid(iPBCLattice(i-1,n),iPBCLattice(j-1,m))
         sigma(9)=grid(iPBCLattice(i-1,n),iPBCLattice(j,m))
         sigma(10)=grid(iPBCLattice(i-1,n),iPBCLattice(j+1,m))
         sigma(11)=grid(iPBCLattice(i-1,n),iPBCLattice(j+2,m))
         sigma(12)=grid(iPBCLattice(i-1,n),iPBCLattice(j+3,m))
         sigma(13)=grid(iPBCLattice(i,n),iPBCLattice(j-2,m))
         sigma(14)=grid(iPBCLattice(i,n),iPBCLattice(j-1,m))
         sigma(15)=grid(iPBCLattice(i,n),iPBCLattice(j+2,m))
         sigma(16)=grid(iPBCLattice(i,n),iPBCLattice(j+3,m))
         sigma(17)=grid(iPBCLattice(i+1,n),iPBCLattice(j-2,m))
         sigma(18)=grid(iPBCLattice(i+1,n),iPBCLattice(j-1,m))
         sigma(19)=grid(iPBCLattice(i+1,n),iPBCLattice(j,m))
         sigma(20)=grid(iPBCLattice(i+1,n),iPBCLattice(j+1,m))
         sigma(21)=grid(iPBCLattice(i+1,n),iPBCLattice(j+2,m))
         sigma(22)=grid(iPBCLattice(i+1,n),iPBCLattice(j+3,m))
         sigma(23)=grid(iPBCLattice(i+2,n),iPBCLattice(j-2,m))
         sigma(24)=grid(iPBCLattice(i+2,n),iPBCLattice(j-1,m))
         sigma(25)=grid(iPBCLattice(i+2,n),iPBCLattice(j,m))
         sigma(26)=grid(iPBCLattice(i+2,n),iPBCLattice(j+1,m))
         sigma(27)=grid(iPBCLattice(i+2,n),iPBCLattice(j+2,m))
         sigma(28)=grid(iPBCLattice(i+2,n),iPBCLattice(j+3,m))
         
      ELSE IF (dir==4) THEN
         sigma(1)=grid(iPBCLattice(i+2,n),iPBCLattice(j+2,m))
         sigma(2)=grid(iPBCLattice(i+2,n),iPBCLattice(j+1,m))
         sigma(3)=grid(iPBCLattice(i+2,n),iPBCLattice(j,m))
         sigma(4)=grid(iPBCLattice(i+2,n),iPBCLattice(j-1,m))
         sigma(5)=grid(iPBCLattice(i+2,n),iPBCLattice(j-2,m))
         sigma(6)=grid(iPBCLattice(i+2,n),iPBCLattice(j-3,m))
         sigma(7)=grid(iPBCLattice(i+1,n),iPBCLattice(j+2,m))
         sigma(8)=grid(iPBCLattice(i+1,n),iPBCLattice(j+1,m))
         sigma(9)=grid(iPBCLattice(i+1,n),iPBCLattice(j,m))
         sigma(10)=grid(iPBCLattice(i+1,n),iPBCLattice(j-1,m))
         sigma(11)=grid(iPBCLattice(i+1,n),iPBCLattice(j-2,m))
         sigma(12)=grid(iPBCLattice(i+1,n),iPBCLattice(j-3,m))
         sigma(13)=grid(iPBCLattice(i,n),iPBCLattice(j+2,m))
         sigma(14)=grid(iPBCLattice(i,n),iPBCLattice(j+1,m))
         sigma(15)=grid(iPBCLattice(i,n),iPBCLattice(j-2,m))
         sigma(16)=grid(iPBCLattice(i,n),iPBCLattice(j-3,m))
         sigma(17)=grid(iPBCLattice(i-1,n),iPBCLattice(j+2,m))
         sigma(18)=grid(iPBCLattice(i-1,n),iPBCLattice(j+1,m))
         sigma(19)=grid(iPBCLattice(i-1,n),iPBCLattice(j,m))
         sigma(20)=grid(iPBCLattice(i-1,n),iPBCLattice(j-1,m))
         sigma(21)=grid(iPBCLattice(i-1,n),iPBCLattice(j-2,m))
         sigma(22)=grid(iPBCLattice(i-1,n),iPBCLattice(j-3,m))
         sigma(23)=grid(iPBCLattice(i-2,n),iPBCLattice(j+2,m))
         sigma(24)=grid(iPBCLattice(i-2,n),iPBCLattice(j+1,m))
         sigma(25)=grid(iPBCLattice(i-2,n),iPBCLattice(j,m))
         sigma(26)=grid(iPBCLattice(i-2,n),iPBCLattice(j-1,m))
         sigma(27)=grid(iPBCLattice(i-2,n),iPBCLattice(j-2,m))
         sigma(28)=grid(iPBCLattice(i-2,n),iPBCLattice(j-3,m))
      END IF       
   END SUBROUTINE CreatesigmaVector28
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   
   FUNCTION iPBCLattice(i,N)
      IMPLICIT NONE
      INTEGER :: i,iPBCLattice,N
      !write(*,*) i
      iPBCLattice=i-N*FLOOR(REAL(i-1)/REAL(N))
      !iPBCLattice
   END FUNCTION iPBCLattice
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintLattice(iunit)
      IMPLICIT NONE
      INTEGER :: i,j,iunit1
      INTEGER, OPTIONAL :: iunit
      
      IF (PRESENT(iunit)) THEN
         iunit1=iunit
      ELSE
         iunit1=6
      END IF
      
      DO i=1,n
         DO j=1,m
            WRITE(UNIT=iunit1,FMT='(i3)',ADVANCE="no") grid(i,j)
         END DO
         WRITE(6,*) ""
      END DO
      
   END SUBROUTINE PrintLattice
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SelectKMCProcess(i,j,dir)
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: i,j,dir
      
      time=time-log(taus88())/kt

      r=taus88()
      kr=r*kt
      sumrates=0._dp
      
      IF (kr<=ktn) THEN
         DO i=1,n
            DO j=1,m
               sumrates=sumrates+gridn(i,j)
               IF (kr<=sumrates) THEN
                  dir=1
                  RETURN
               END IF
            END DO
         END DO
      ELSE
         sumrates=sumrates+ktn
      END IF
      
      IF (kr<=ktn+kts) THEN
         DO i=1,n
            DO j=1,m
               sumrates=sumrates+grids(i,j)
               IF (kr<=sumrates) THEN
                  dir=2
                  RETURN
               END IF
            END DO
         END DO
      ELSE
         sumrates=sumrates+kts
      END IF
      
      IF (kr<=ktn+kts+kte) THEN
         DO i=1,n
            DO j=1,m
               sumrates=sumrates+gride(i,j)
               IF (kr<=sumrates) THEN
                  dir=3
                  RETURN
               END IF
            END DO
         END DO
      ELSE
         sumrates=sumrates+kte
      END IF
      
      DO i=1,n
         DO j=1,m
            sumrates=sumrates+gridw(i,j)
            IF (kr<=sumrates) THEN
               dir=4
               RETURN
            END IF
         END DO
      END DO
      
   END SUBROUTINE SelectKMCProcess
 !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintLattice1()
      IMPLICIT NONE
      INTEGER :: i,j
      
      OPEN(UNIT=10,FILE="PrintLattice1Result.txt")
      WRITE(10,*), "The rates as atom moves to north"
      DO i=1,n
         DO j=1,m
            WRITE(UNIT=10,FMT='(E15.7)',ADVANCE="no") gridn(i,j)
         END DO
         WRITE(10,*) ""   
      END DO
      
      WRITE(10,*), "The rates as atom moves to south"
       DO i=1,n
         DO j=1,m
            WRITE(UNIT=10,FMT='(E15.7)',ADVANCE="no") grids(i,j)
        END DO
        WRITE(10,*) ""
      END DO
      
      WRITE(10,*), "The rates as atom moves to east"
      DO i=1,n
         DO j=1,m
            WRITE(UNIT=10,FMT='(E15.7)',ADVANCE="no") gride(i,j)
         END DO
         WRITE(10,*) ""
      END DO
      
      WRITE(10,*), "The rates as atom moves to west"
      DO i=1,n
         DO j=1,m
            WRITE(UNIT=10,FMT='(E15.7)',ADVANCE="no") gridw(i,j)
         END DO
         WRITE(10,*) ""
      END DO
      CLOSE(10)
      
   END SUBROUTINE PrintLattice1
 !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 FUNCTION NEBLatticeKMCActivationEnergy(sigma)
   IMPLICIT NONE
   INTEGER :: sigma(natoms)
   REAL(dp) :: NEBLatticeKMCActivationEnergy,TSEnergy,springconstant
   INTEGER :: i,ntransitionstates,indx,nadatoms,errorstatus
   
   NEBLatticeKMCActivationEnergy=0._dp
   
   NumberNEBCalls=NumberNEBCalls+1
   
   nadatoms=SUM(sigma)
OPEN(UNIT=539,FILE="initialstate")
WRITE(539,*) nadatoms+501,nadatoms+351,nadatoms+351,nadatoms+351
WRITE(539,*) "    20.45000076    20.45000076    34.95000076     0.00000000     0.00000000     0.00000000  "
IF (sigma(1)==1)  WRITE(UNIT=539,FMT='("Ag           0.97520000         13.28870030         21.41877971   T T T")')
IF (sigma(2)==1)  WRITE(UNIT=539,FMT='("Ag           3.01970000         15.33300000         21.41877971   T T T")')
IF (sigma(3)==1)  WRITE(UNIT=539,FMT='("Ag           5.06420000         17.37770030         21.41877971   T T T")')
IF (sigma(4)==1)  WRITE(UNIT=539,FMT='("Ag           7.10870000         19.42200000         21.41877971   T T T")')
IF (sigma(5)==1)  WRITE(UNIT=539,FMT='("Ag           9.15320000         21.46670000         21.41877971   T T T")')
IF (sigma(6)==1)  WRITE(UNIT=539,FMT='("Ag          11.19800000         23.51120000         21.41877971   T T T")')
IF (sigma(12)==1)  WRITE(UNIT=539,FMT='("Ag          13.24224327         21.46670000         21.41877971   T T T")')
IF (sigma(16)==1)  WRITE(UNIT=539,FMT='("Ag          15.28700000         19.42200000         21.41877971   T T T")')
IF (sigma(22)==1)  WRITE(UNIT=539,FMT='("Ag          17.33120000         17.37770030         21.41877971   T T T")')
IF (sigma(28)==1) WRITE(UNIT=539,FMT='("Ag          19.37570000         15.33300000         21.41877971   T T T")')
IF (sigma(27)==1) WRITE(UNIT=539,FMT='("Ag          17.33120000         13.28870030         21.41877971   T T T")')
IF (sigma(26)==1) WRITE(UNIT=539,FMT='("Ag          15.28700000         11.24400000         21.41877971   T T T")')
IF (sigma(25)==1) WRITE(UNIT=539,FMT='("Ag          13.24224327          9.19970000         21.41877971   T T T")')
IF (sigma(24)==1) WRITE(UNIT=539,FMT='("Ag          11.19800000          7.15220000         21.41877971   T T T")')
IF (sigma(23)==1) WRITE(UNIT=539,FMT='("Ag           9.15320000          5.11070000         21.41877971   T T T")')
IF (sigma(17)==1) WRITE(UNIT=539,FMT='("Ag           7.10870000          7.15220000         21.41877971   T T T")')
IF (sigma(13)==1) WRITE(UNIT=539,FMT='("Ag           5.06420000          9.19970000         21.41877971   T T T")')
IF (sigma(7)==1) WRITE(UNIT=539,FMT='("Ag           3.01970000         11.24400000         21.41877971   T T T")')
IF (sigma(8)==1) WRITE(UNIT=539,FMT='("Ag           5.06420000         13.28870030         21.41877971   T T T")')
IF (sigma(9)==1) WRITE(UNIT=539,FMT='("Ag           7.10870000         15.33300000         21.41877971   T T T")')
IF (sigma(10)==1) WRITE(UNIT=539,FMT='("Ag           9.15320000         17.37770030         21.41877971   T T T")')
IF (sigma(11)==1) WRITE(UNIT=539,FMT='("Ag          11.19800000         19.42200000         21.41877971   T T T")')
IF (sigma(18)==1) WRITE(UNIT=539,FMT='("Ag           9.15324327          9.19970000         21.41877971   T T T")')
IF (sigma(19)==1) WRITE(UNIT=539,FMT='("Ag          11.19800000         11.24400000         21.41877971   T T T")')
IF (sigma(20)==1) WRITE(UNIT=539,FMT='("Ag          13.24200000         13.28870030         21.41877971   T T T")')
IF (sigma(21)==1) WRITE(UNIT=539,FMT='("Ag          15.28700000         15.33300000         21.41877971   T T T")')
IF (sigma(14)==1) WRITE(UNIT=539,FMT='("Ag           7.10870000         11.24400000         21.41877971   T T T")')
IF (sigma(15)==1) WRITE(UNIT=539,FMT='("Ag          13.24200000         17.37770030         21.41877971   T T T")')
WRITE(UNIT=539,FMT='("Ag           9.15324327         13.28870030         21.41877971   T T T")')
CLOSE(539)

filename="chos_adatomhop.xyz"
CALL SYSTEM("cat initialstate AgSlab > "//TRIM(filename))

OPEN(UNIT=539,FILE="finalstate")
WRITE(539,*) nadatoms+501,nadatoms+351,nadatoms+351,nadatoms+351
WRITE(539,*) "    20.45000076    20.45000076    34.95000076     0.00000000     0.00000000     0.00000000  "
IF (sigma(1)==1)  WRITE(UNIT=539,FMT='("Ag           0.97520000         13.28870030         21.41877971   T T T")')
IF (sigma(2)==1)  WRITE(UNIT=539,FMT='("Ag           3.01970000         15.33300000         21.41877971   T T T")')
IF (sigma(3)==1)  WRITE(UNIT=539,FMT='("Ag           5.06420000         17.37770030         21.41877971   T T T")')
IF (sigma(4)==1)  WRITE(UNIT=539,FMT='("Ag           7.10870000         19.42200000         21.41877971   T T T")')
IF (sigma(5)==1)  WRITE(UNIT=539,FMT='("Ag           9.15320000         21.46670000         21.41877971   T T T")')
IF (sigma(6)==1)  WRITE(UNIT=539,FMT='("Ag          11.19800000         23.51120000         21.41877971   T T T")')
IF (sigma(12)==1)  WRITE(UNIT=539,FMT='("Ag          13.24224327         21.46670000         21.41877971   T T T")')
IF (sigma(16)==1)  WRITE(UNIT=539,FMT='("Ag          15.28700000         19.42200000         21.41877971   T T T")')
IF (sigma(22)==1)  WRITE(UNIT=539,FMT='("Ag          17.33120000         17.37770030         21.41877971   T T T")')
IF (sigma(28)==1) WRITE(UNIT=539,FMT='("Ag          19.37570000         15.33300000         21.41877971   T T T")')
IF (sigma(27)==1) WRITE(UNIT=539,FMT='("Ag          17.33120000         13.28870030         21.41877971   T T T")')
IF (sigma(26)==1) WRITE(UNIT=539,FMT='("Ag          15.28700000         11.24400000         21.41877971   T T T")')
IF (sigma(25)==1) WRITE(UNIT=539,FMT='("Ag          13.24224327          9.19970000         21.41877971   T T T")')
IF (sigma(24)==1) WRITE(UNIT=539,FMT='("Ag          11.19800000          7.15220000         21.41877971   T T T")')
IF (sigma(23)==1) WRITE(UNIT=539,FMT='("Ag           9.15320000          5.11070000         21.41877971   T T T")')
IF (sigma(17)==1) WRITE(UNIT=539,FMT='("Ag           7.10870000          7.15220000         21.41877971   T T T")')
IF (sigma(13)==1) WRITE(UNIT=539,FMT='("Ag           5.06420000          9.19970000         21.41877971   T T T")')
IF (sigma(7)==1) WRITE(UNIT=539,FMT='("Ag           3.01970000         11.24400000         21.41877971   T T T")')
IF (sigma(8)==1) WRITE(UNIT=539,FMT='("Ag           5.06420000         13.28870030         21.41877971   T T T")')
IF (sigma(9)==1) WRITE(UNIT=539,FMT='("Ag           7.10870000         15.33300000         21.41877971   T T T")')
IF (sigma(10)==1) WRITE(UNIT=539,FMT='("Ag           9.15320000         17.37770030         21.41877971   T T T")')
IF (sigma(11)==1) WRITE(UNIT=539,FMT='("Ag          11.19800000         19.42200000         21.41877971   T T T")')
IF (sigma(18)==1) WRITE(UNIT=539,FMT='("Ag           9.15324327          9.19970000         21.41877971   T T T")')
IF (sigma(19)==1) WRITE(UNIT=539,FMT='("Ag          11.19800000         11.24400000         21.41877971   T T T")')
IF (sigma(20)==1) WRITE(UNIT=539,FMT='("Ag          13.24200000         13.28870030         21.41877971   T T T")')
IF (sigma(21)==1) WRITE(UNIT=539,FMT='("Ag          15.28700000         15.33300000         21.41877971   T T T")')
IF (sigma(14)==1) WRITE(UNIT=539,FMT='("Ag           7.10870000         11.24400000         21.41877971   T T T")')
IF (sigma(15)==1) WRITE(UNIT=539,FMT='("Ag          13.24200000         17.37770030         21.41877971   T T T")')
WRITE(UNIT=539,FMT='("Ag          11.26165838         15.35124533         21.41877971   T T T")')
CLOSE(539)
CALL SYSTEM("cat "//TRIM(filename)//" finalstate AgSlab > combinedfile")
CALL SYSTEM("mv combinedfile "//TRIM(filename))
CALL SYSTEM("rm initialstate finalstate ")

      CALL ReadCoord(AL,filename)
      CALL Analyze(AL)
      CALL PotentialInitialize(AL)
      CALL WritePotential(AL)
   
      CALL PotentialInitialize(AL%NextNeigh)
      CALL WritePotential(AL%NextNeigh)
      
      OPEN(UNIT=510,FILE="chos_energydiag.xyz")
      springconstant=1._dp
      CALL NEB(state1=AL,state2=AL%NextNeigh,nimages=11,imode=1101,omode=13,interpmode=1, &
         springconstant=springconstant,ITMAX=50,TSEnergy=TSEnergy,ntransitionstates=ntransitionstates, &
         ReuseVL=.FALSE.,TS=AL2,errorstatus=errorstatus)
      
      CALL FLUSH(509)
      CLOSE(510)
      CALL Delete(AL)
      CALL Delete(AL2)
      
      NEBLatticeKMCActivationEnergy=TSEnergy

 END FUNCTION NEBLatticeKMCActivationEnergy
 !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE LatticeKMCModule
