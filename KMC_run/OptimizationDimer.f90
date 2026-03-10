MODULE OptimizationDimer
   !for optimizing different structures:
   !e.g., AL,ChOS
      !omode=1 (SD)
      !omode=2 (CG)
      !omode=3 (LBFGS)
      !omode=4 (dynamical-Euler)
      !omode=5 (dynamical-RK method)
      !omode=6 (QuickMin)
      !omode=7 (FIRE)
   USE VARIABLE_TYPE
   USE NeighborList
   USE PotentialPackage
   USE OptimPackage
   USE utilities
   USE IO
   IMPLICIT NONE
   TYPE(SystemContainer), POINTER, PRIVATE :: AL,AL1,AL2
   TYPE(ChOSContainer), POINTER, PRIVATE :: chos
   INTEGER, PRIVATE :: NAtoms,iprint,TSImage,rmode1,tmode1
   INTEGER, PRIVATE :: rn1,rn2
   LOGICAL, PRIVATE :: CI_NEB1
   REAL(dp), PRIVATE :: SclFac,drmag,TSEnergy,curv
   REAL(dp), DIMENSION(:), POINTER, PRIVATE :: n,t,f,r,r1,r2,dr
   REAL(dp), DIMENSION(:), POINTER, PRIVATE :: AtomCoord,AtomCoord1,AtomCoord2
   REAL(dp), DIMENSION(:), POINTER, PRIVATE :: AtomForce,AtomForce1,AtomForce2
   LOGICAL, DIMENSION(:), POINTER, PRIVATE :: AtomIsMoving,AtomIsMoving1,AtomIsMoving2

   
   INTERFACE Optimize
      MODULE PROCEDURE OptimizeDimer
   END INTERFACE
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE OptimizeDimer(chos,rmode,tmode)
   !performs a dimer search using AL -- note that AL would be
   !obtained either using a MD simulation or by randomly displacing atoms
   
   !rmode1 is the rotation optimization mode
   !rmode1=0 conjugate gradient + Newton minimization
   !rmode1=1 steepest descent search + Newton minimization
   !rmode1=2 LBFGS with Newton minimization
   !rmode1=3 Lagrange multiplier based LBFGS minimization
   
   !tmode1 is the translation mode
   !tmode1=0 is the 
   !tmode1=1 is conjugate gradient based
   !tmode1=2 is LBFGS based
      IMPLICIT NONE
      TYPE(ChOSContainer), POINTER :: chos
      INTEGER :: errorstatus,rmode,tmode,iter
      REAL(dp) :: DimerEnergy
      CHARACTER(len=100) :: filename
      
      rmode1=rmode
      tmode1=tmode
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Bit of housekeeping
      AL1=>chos%AL
      AL2=>AL1%NextNeigh
      NAtoms=AL1%NAtoms
      CALL MakeSize(AL,NAtoms)
      AL%AtomCoord=0.5_dp*(AL1%AtomCoord+AL2%AtomCoord)
      
      AtomCoord1=>AL1%AtomCoord
      AtomCoord2=>AL2%AtomCoord
      AtomCoord=>AL%AtomCoord
      AtomIsMoving1=>AL1%AtomIsMoving
      AtomIsMoving2=>AL2%AtomIsMoving
      !AtomIsMoving=>AL%AtomIsMoving
      AtomForce1=>AL1%AtomForce
      AtomForce2=>AL2%AtomForce
      !AtomForce=>AL%AtomForce
      
      CALL MakeSize(n,3*NAtoms)
      CALL MakeSize(t,3*NAtoms)
      CALL MakeSize(dr,3*NAtoms)
      
      CALL AddVerletList(AL1,NAtoms=NAtoms)
      CALL AddVerletList(AL2,NAtoms=NAtoms)
      n=0._dp
      n=AL1%AtomCoord-AL%AtomCoord !WHERE (AL1%AtomIsMoving) n=AL1%AtomCoord-AL%AtomCoord
      drmag=NORM(n,3*NAtoms)
      
      !OPEN(UNIT=356,FILE="Dimer.xyz")
      DO iter=1,10
         CALL OptimizeDimerRotation(errorstatus,DimerEnergy,curv)
         CALL OptimizeDimerTranslation(errorstatus,DimerEnergy,curv)
      END DO
      !CLOSE(356)
   END SUBROUTINE OptimizeDimer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE OptimizeDimerRotation(errorstatus,DimerEnergy,curv)
   !Performs a dimer rotation
      IMPLICIT NONE
      INTEGER :: iter,errorstatus
      REAL(dp) :: ftol,gtol,xtol,fret,DimerEnergy,curv
      
      SELECT CASE (rmode1) !rotate dimer to min mode
      CASE(1) !steepest descent
         SclFac=0.01_dp !to be used for scaling forces
         ftol=1.0e-9_dp !this should come from potential
         gtol=2.0e-8_dp
         xtol=1.0e-5_dp
         CALL DimerRotationSD(ftol,gtol,xtol,iter,errorstatus,curv)
      CASE(2) !conjugate gradient
         SclFac=0.01_dp !to be used for scaling forces
         ftol=1.0e-9_dp !this should come from potential
         gtol=2.0e-6_dp
         xtol=1.0e-5_dp
         CALL DimerRotationCG(ftol,gtol,xtol,iter,errorstatus)
      CASE DEFAULT
         WRITE(6,*) "This option has not been implemented for dimer rotate"
         STOP
      END SELECT
      DimerEnergy=0.5_dp*(AL1%PotentialEnergy+AL2%PotentialEnergy) + &
         0.25_dp*drmag*DOT_PRODUCT(AL1%AtomForce-AL2%AtomForce,n)
   END SUBROUTINE OptimizeDimerRotation
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DimerRotationSD(ftol,gtol,xtol,iter,errorstatus,curv)
   !employs the steepest-descent method
      IMPLICIT NONE
      INTEGER, PARAMETER :: maxiter=50
      INTEGER :: iter,errorstatus,result
      REAL(dp) :: ftol,gtol,xtol
      REAL(dp), DIMENSION(:), POINTER :: fp,fp1,n1,t1,tprime
      REAL(dp) :: theta,e0,eprev,curv,frotation,theta0
      
      !xtol theta
      !ftol e0
      !gtol frotation
      
      !CALL Figure3()
      !OPEN(UNIT=356,FILE="Tmp.xyz")
      
      DO iter=1,maxiter
         CALL SteepestDescentRotate(theta0,e0,frotation,curv)
      END DO
      
      !CLOSE(356)
   END SUBROUTINE DimerRotationSD
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Figure3()
   !creates Fig.3 of Henkelman's Dimer paper
      IMPLICIT NONE
      INTEGER :: iter
      REAL(dp), DIMENSION(:), POINTER :: fp,fp1,n1,t1,tprime
      REAL(dp) :: theta,e0,eprev,dtheta,theta0
      
      CALL MakeSize(fp,3*NAtoms)
      CALL MakeSize(fp1,3*NAtoms)
      CALL MakeSize(n1,3*NAtoms)
      CALL MakeSize(t1,3*NAtoms)
      CALL MakeSize(tprime,3*NAtoms)
      
      !Obtain n and t for the current orientation
      CALL GetEnergyForceDimerQuick()
      e0=AL1%PotentialEnergy+AL2%PotentialEnergy
      eprev=e0
      CALL GetUnitVectorN(n)
      CALL GetForcePerpendicular(fp,n)
      CALL GetUnitVectorT(t,fp)
      curv=0.5_dp*DOT_PRODUCT(AtomForce2-AtomForce1,n)/drmag
      
      !obtain theta0
      dtheta=0.001_dp
      CALL RotateDimer(dtheta,n,t)
      CALL GetEnergyForceDimerQuick()
      CALL GetUnitVectorN(n1)
      CALL GetForcePerpendicular(fp1,n1)
      CALL GetUnitVectorT(t1,fp1)
      theta0=-0.5_dp*ATAN( &
         (DOT_PRODUCT(fp,t)+DOT_PRODUCT(fp1,t1))* &
         dtheta/(DOT_PRODUCT(fp1,t1)-DOT_PRODUCT(fp,t)))-dtheta*0.5_dp
      eprev=e0
      e0=AL1%PotentialEnergy+AL2%PotentialEnergy
      IF (e0>eprev) theta0=theta0+PI*0.5_dp
      
      !Effect of rotation
      theta=0.46_dp
      !theta=0._dp
      dtheta=0.02_dp*PI
      !dtheta=.0001_dp
      DO iter=1,100
         theta=theta+dtheta
         CALL RotateDimer(theta,n,t)
         CALL GetEnergyForceDimerQuick()
         e0=AL1%PotentialEnergy+AL2%PotentialEnergy
         CALL GetUnitVectorN(n1)
         CALL GetForcePerpendicular(fp1,n1)
         CALL GetUnitVectorT(t1,fp1)
         tprime=t*COS(theta)-n*SIN(theta)
         curv=0.5_dp*DOT_PRODUCT(AtomForce2-AtomForce1,n1)/drmag
         WRITE(6,*) theta,e0,DOT_PRODUCT(AtomForce1-AtomForce2,tprime),&
            (eprev-e0)/drmag/dtheta,NORM(fp1,3*NAtoms),curv
         eprev=e0
      END DO
      stop
   END SUBROUTINE Figure3
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DimerRotationCG(ftol,gtol,xtol,iter,errorstatus)
   !employs the conjugate-gradient method
      IMPLICIT NONE
      INTEGER :: iter,errorstatus
      REAL(dp) :: ftol,gtol,xtol
   END SUBROUTINE DimerRotationCG
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE OptimizeDimerTranslation(errorstatus,DimerEnergy,curv)
   !Performs the dimer translation
      IMPLICIT NONE
      REAL(dp) :: DimerEnergy
      REAL(dp), DIMENSION(:), POINTER :: fr,feff,dx
      REAL(dp) :: curv
      INTEGER :: errorstatus
      
      CALL GetEnergyForceDimerQuick()
      write(*,*) curv
      call flush(6)
      ALLOCATE(fr(3*natoms))
      ALLOCATE(feff(3*natoms))
      ALLOCATE(dx(3*natoms))
      fr=0.5_dp*(AtomForce1+AtomForce2)
      CALL GetUnitVectorN(n)
      feff=-DOT_PRODUCT(fr,n)*n
      IF (curv<0._dp) THEN
         feff=fr+2._dp*feff
      END IF
      dx=0.01_dp*feff/NORM(feff,natoms)
      WHERE(AtomIsMoving1) AtomCoord1=AtomCoord1+dx
      WHERE(AtomIsMoving1) AtomCoord2=AtomCoord2+dx
      CALL WriteXYZ(AL1,fileopen=.FALSE.,fileclose=.FALSE.,iunit=356,iprint=0,NAtoms=AL1%NAtoms)
      CALL GetEnergyForceDimerQuick()
      WRITE(*,*) "Completed rotation"
      call flush(6)
      
      DEALLOCATE(fr,feff,dx)
      
   END SUBROUTINE OptimizeDimerTranslation
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION DimerTranslationEnergy(x,errorstatus)
   !x corresponds to the position of the dimer center
      USE VARIABLE_TYPE
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), INTENT(IN) :: x
      REAL(dp) :: DimerTranslationEnergy
      INTEGER :: errorstatus
      
      DimerTranslationEnergy=0._dp
   END FUNCTION DimerTranslationEnergy
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION DimerTranslationForce(x,errorstatus)
      USE VARIABLE_TYPE
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), INTENT(IN) :: x
      REAL(dp), DIMENSION(SIZE(x)) :: DimerTranslationForce
      INTEGER :: errorstatus
      
      DimerTranslationForce=0._dp
   END FUNCTION DimerTranslationForce
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetEnergyForceDimerQuick()
   !dimer has moved to a new position, find the new energy and force
      IMPLICIT NONE
      INTEGER :: errorstatus
      
      CALL AddVerletListQuick(AL1,NAtoms)
      CALL GetForcesVL(AL1,NAtoms,errorstatus)
      IF (errorstatus/=0) STOP
      CALL AddVerletListQuick(AL2,NAtoms)
      CALL GetForcesVL(AL2,NAtoms,errorstatus)
      IF (errorstatus/=0) STOP
   END SUBROUTINE GetEnergyForceDimerQuick
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RotateDimer(theta,n,t)
   !rotates the dimer by theta given n and t
      IMPLICIT NONE
      REAL(dp) :: theta
      REAL(dp), DIMENSION(:), POINTER :: n,t
   
      WHERE(AtomIsMoving1) AtomCoord1=AtomCoord+drmag*(n*COS(theta)+t*SIN(theta))
      WHERE(AtomIsMoving2) AtomCoord2=2._dp*AtomCoord-AtomCoord1
   END SUBROUTINE RotateDimer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetUnitVectorN(n)
   !finds the n unit vector for a given dimer orientation
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), POINTER :: n
      
      n=0._dp
      WHERE(AtomIsMoving1) n=AtomCoord1-AtomCoord
      n=n/drmag
   END SUBROUTINE GetUnitVectorN
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetForcePerpendicular(fp,n)
   !finds the perpendicular force
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), POINTER :: fp,n
      
      fp=AtomForce1-AtomForce2
      fp=fp-DOT_PRODUCT(fp,n)*n
   END SUBROUTINE GetForcePerpendicular
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetUnitVectorT(t,fp)
   !finds the t vector given fp
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), POINTER :: fp,t
      
      t=0._dp
      WHERE(AtomIsMoving1) t=fp/NORM(fp,3*NAtoms)
   END SUBROUTINE GetUnitVectorT
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SteepestDescentRotate(theta0,e0,frotation,curv)
      IMPLICIT NONE
      INTEGER :: iter
      REAL(dp), DIMENSION(:), POINTER :: fp,fp1,n1,t1
      REAL(dp) :: theta,e0,e1,eprev,dtheta,theta0,frotation,curv
      
      CALL MakeSize(fp,3*NAtoms)
      CALL MakeSize(fp1,3*NAtoms)
      CALL MakeSize(n1,3*NAtoms)
      CALL MakeSize(t1,3*NAtoms)
      
      !Obtain n and t for the current orientation
      write(*,*) "a"; call flush(6)
      CALL GetEnergyForceDimerQuick()
      e0=AL1%PotentialEnergy+AL2%PotentialEnergy
      eprev=e0
      CALL GetUnitVectorN(n)
      CALL GetForcePerpendicular(fp,n)
      CALL GetUnitVectorT(t,fp)
      curv=0.5_dp*DOT_PRODUCT(AtomForce2-AtomForce1,n)/drmag
      write(*,*) "b"; call flush(6)
      !obtain theta0
      dtheta=0.001_dp
      CALL RotateDimer(dtheta,n,t)
      CALL GetEnergyForceDimerQuick()
      CALL GetUnitVectorN(n1)
      CALL GetForcePerpendicular(fp1,n1)
      CALL GetUnitVectorT(t1,fp1)
      theta0=-0.5_dp*ATAN( &
         (DOT_PRODUCT(fp,t)+DOT_PRODUCT(fp1,t1))* &
         dtheta/(DOT_PRODUCT(fp1,t1)-DOT_PRODUCT(fp,t)))-dtheta*0.5_dp
      write(*,*) "c"; call flush(6)
      eprev=e0
      CALL RotateDimer(theta0,n,t)
      CALL GetEnergyForceDimerQuick()
      e0=AL1%PotentialEnergy+AL2%PotentialEnergy
      IF (e0>eprev) THEN
         theta0=theta0+PI*0.5_dp
         CALL RotateDimer(theta0,n,t)
         CALL GetEnergyForceDimerQuick()
         e1=AL1%PotentialEnergy+AL2%PotentialEnergy
         IF (e1>e0) THEN
            theta0=theta0-PI*0.5_dp
            CALL RotateDimer(theta0,n,t)
            CALL GetEnergyForceDimerQuick()
         END IF
      END IF
      write(*,*) "d"; call flush(6)
      CALL GetUnitVectorN(n)
      CALL GetForcePerpendicular(fp,n)
      CALL GetUnitVectorT(t,fp)
      curv=0.5_dp*DOT_PRODUCT(AtomForce2-AtomForce1,n)/drmag
      frotation=DOT_PRODUCT(AtomForce1-AtomForce2,t)
      WRITE(6,*) theta0,e0,frotation,curv
      call flush(6)
      !IF (ABS(theta0)<1.e-3_dp .AND. ABS(frotation)<2.e-3_dp) THEN
      !   SteepestDescentRotate=1
      !END IF
      !CALL WriteXYZ(AL1,fileopen=.FALSE.,fileclose=.FALSE.,iunit=356,iprint=0,NAtoms=AL1%NAtoms)
      DEALLOCATE(fp,fp1,n1,t1)
      write(*,*) "e"; call flush(6)
   END SUBROUTINE SteepestDescentRotate
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE OptimizationDimer
