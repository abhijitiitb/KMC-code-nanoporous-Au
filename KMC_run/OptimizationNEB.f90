MODULE OptimizationNEB
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
   INTEGER, PARAMETER, PRIVATE :: MaxNImages=45
   TYPE(SystemContainer), POINTER, PRIVATE :: AL
   TYPE(ChOSContainer), POINTER, PRIVATE :: chos
   INTEGER, PRIVATE :: NAtoms,NImages,iprint,TSImage,NOptimize,ActiveImage
   LOGICAL, PRIVATE :: CI_NEB1
   REAL(dp), PRIVATE :: SclFac,TSEnergy,SpringConst
   REAL(dp), DIMENSION(3*MaxNEBAtoms), PRIVATE :: t
   REAL(dp), DIMENSION(3*MaxNEBAtoms*(MaxNImages-2)), TARGET, PRIVATE :: MovingAtomCoord
   INTEGER, DIMENSION(3*MaxNEBAtoms), PRIVATE :: AtomCoordIndx
   LOGICAL :: ReuseVL1 !allows reuse of VL and LL belonging to first image for all other images
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !    OPTIMIZATION ROUTINES FOR NEB CALCULATIONS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE OptimizeChOS(chos1,omode1,CI_NEB,ITMAX,errorstatus)
      !given a chain of states containing nimages, the subroutine moves the
      !images to the minimum energy path
      IMPLICIT NONE
      TYPE(ChOSContainer), POINTER :: chos1
      INTEGER, OPTIONAL :: omode1
      INTEGER :: omode,i,j,iter,NEBIteration,ITMAX,errorstatus,image,indx
      REAL(dp) :: ftol,gtol,xtol,fret,InitialEnergy
      LOGICAL :: NotSatisfied,CI_NEB
      REAL(dp), DIMENSION(:), POINTER :: x,AtomCoord
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      REAL(dp), DIMENSION(MaxNImages) :: EnergyChangeImage !maximum 200 images can be stored
      
      errorstatus=0

      chos=>chos1
      IF (.NOT. ASSOCIATED(chos1)) THEN
         WRITE(6,*) "$Err>> ChOS is not associated"
         STOP
      END IF
      NImages=chos%nimages
      IF (NImages>MaxNImages) THEN
         WRITE(6,*) "Err>> Increase the MaxNImages in OptimizationNEB"
         STOP
      END IF
      IF (NImages<3) THEN
         WRITE(6,*) "$Err>> ChOS does contains very few images. NImages",NImages
         STOP
      END IF
      
      CI_NEB1=CI_NEB
      NAtoms=chos%AL%NAtoms
      SpringConst=chos%SpringConst

      IF (PRESENT(omode1)) THEN
         omode=omode1
      ELSE
         omode=3
      END IF
      iter=0
      ftol=1.e-10_dp
      
      WRITE(*,*) "Optim>> Performing NEB calculation"
      
      !Find how many atoms are moving (the dimensionality of the optimization space per image)
      !total degrees of freedom are going to be NOptimize*(NImages-2)
      AL=>chos%AL
      AtomCoord=>AL%AtomCoord
      AtomIsMoving=>AL%AtomIsMoving
      NOptimize=0
      DO indx=1,3*NAtoms
         IF (AtomIsMoving(indx)) THEN
            NOptimize=NOptimize+1
            AtomCoordIndx(NOptimize)=indx
         END IF
      END DO
      !Add verlet list and compute the energies and forces of all remaining images
      CALL AddVerletList(AL,NAtoms=NAtoms)
      CALL GetForcesVL(AL,NAtoms,errorstatus)
      IF (errorstatus/=0) THEN !the atoms are overlapping -- bad chos provided
         errorstatus=3
         RETURN
      END IF
      InitialEnergy=chos%AL%PotentialEnergy

      DO image=2,NImages-1
         AL=>AL%NextNeigh
         CALL AddVerletList(AL,NAtoms=NAtoms)
         CALL GetForcesVL(AL,NAtoms,errorstatus)
         IF (errorstatus/=0) THEN !the atoms are overlapping -- bad chos provided
            errorstatus=3
            RETURN
         END IF
         AtomCoord=>AL%AtomCoord
         j=(image-2)*NOptimize
         DO i=1,NOptimize
            indx=AtomCoordIndx(i)
            MovingAtomCoord(i+j)=AtomCoord(indx)
         END DO
      END DO
      
      EnergyChangeImage=0._dp !note first and last images can be ignored
      
      !Loop over each image and minimize the forces till converged
      NotSatisfied=.TRUE.
      NEBIteration=0
      TSEnergy=0._dp
      TSImage=0
      DO WHILE (NotSatisfied)
         NEBIteration=NEBIteration+1
         AL=>chos%AL
         DO image=2,NImages-1,2 !even numbered images
            !WRITE(6,*) "Image number:",i
            AL=>AL%NextNeigh
            EnergyChangeImage(image)=-AL%PotentialEnergy
            !NAtoms=AL%NAtoms
            x=>MovingAtomCoord((image-2)*NOptimize+1:(image-1)*NOptimize)
WRITE(6,*) "Checking if local version has correct size:",SIZE(x),SIZE(MovingAtomCoord)
STOP
            SclFac=1._dp
            ActiveImage=image
            SELECT CASE(omode)
            CASE(1) !steepest descent
               SclFac=0.01_dp !to be used for scaling forces
               ftol=1.0e-9_dp !this should come from potential
               gtol=2.0e-6_dp
               xtol=1.0e-5_dp
               CALL SteepestDescent(x,ftol,gtol,xtol,iter,fret,ChOSEnergy,ChOSForce, &
                    dx=0.0001_dp,ITMAX=ITMAX,IsPrint1=.TRUE.,errorstatus=errorstatus)
            CASE(2) !conjugate gradient
               SclFac=0.01_dp
               ftol=1.0e-9_dp !this should come from potential
               gtol=2.0e-6_dp
               xtol=1.0e-5_dp
               CALL ConjugateGradient(x,ftol,gtol,xtol,iter,fret,ChOSEnergy,ChOSForce, &
                    dx=0.0001_dp,ITMAX=ITMAX,IsPrint1=.TRUE.,errorstatus=errorstatus)
            CASE(3) !LBFGS
               ftol=1.0e-9_dp !this should come from potential
               gtol=2.0e-6_dp
               xtol=1.0e-5_dp
               CALL Lbfgs(x,ftol,gtol,xtol,iter,fret,ChOSEnergy,ChOSForce,IsPrint1=.TRUE., &
            ITMAX=100,errorstatus=errorstatus)
            CASE(4) !dynamical Euler
               !CALL 
            CASE(5) !Runge Kutta
            CASE(6) !QuickMin
            CASE(7) !FIRE
            CASE(8) !DFP
               ftol=1.0e-9_dp !this should come from potential
               gtol=2.0e-6_dp
               xtol=1.0e-5_dp
               CALL DFP(x,ftol,gtol,xtol,iter,fret,ChOSEnergy,ChOSForce,IsPrint1=.TRUE., &
                    errorstatus=errorstatus)
            CASE DEFAULT
            END SELECT
            IF (errorstatus/=0) STOP !the energy and force calculation had error
            EnergyChangeImage(image)=EnergyChangeImage(image)+AL%PotentialEnergy
            AL=>AL%NextNeigh !to jump to the next odd position
         END DO
         
         AL=>chos%AL%NextNeigh
         DO image=3,NImages-1,2 !odd numbered images
            !WRITE(6,*) "Image number:",i
            AL=>AL%NextNeigh
            EnergyChangeImage(image)=-AL%PotentialEnergy
            !NAtoms=AL%NAtoms
            x=>MovingAtomCoord((image-2)*NOptimize+1:(image-1)*NOptimize)
            SclFac=1._dp
            ActiveImage=image
            SELECT CASE(omode)
            CASE(1) !steepest descent
               SclFac=0.01_dp !to be used for scaling forces
               ftol=1.0e-9_dp !this should come from potential
               gtol=2.0e-6_dp
               xtol=1.0e-5_dp
               CALL SteepestDescent(x,ftol,gtol,xtol,iter,fret,ChOSEnergy,ChOSForce, &
                    dx=0.0001_dp,ITMAX=ITMAX,IsPrint1=.TRUE.,errorstatus=errorstatus)
            CASE(2) !conjugate gradient
               SclFac=0.01_dp
               ftol=1.0e-9_dp !this should come from potential
               gtol=2.0e-6_dp
               xtol=1.0e-5_dp
               CALL ConjugateGradient(x,ftol,gtol,xtol,iter,fret,ChOSEnergy,ChOSForce, &
                    dx=0.0001_dp,ITMAX=ITMAX,IsPrint1=.TRUE.,errorstatus=errorstatus)
            CASE(3) !LBFGS
               ftol=1.0e-9_dp !this should come from potential
               gtol=2.0e-6_dp
               xtol=1.0e-5_dp
               CALL Lbfgs(x,ftol,gtol,xtol,iter,fret,ChOSEnergy,ChOSForce,IsPrint1=.TRUE., &
                    ITMAX=100,errorstatus=errorstatus)
            CASE(4) !dynamical Euler
               !CALL 
            CASE(5) !Runge Kutta
            CASE(6) !QuickMin
            CASE(7) !FIRE
            CASE(8) !DFP
               ftol=1.0e-9_dp !this should come from potential
               gtol=2.0e-6_dp
               xtol=1.0e-5_dp
               CALL DFP(x,ftol,gtol,xtol,iter,fret,ChOSEnergy,ChOSForce,IsPrint1=.TRUE., &
                    errorstatus=errorstatus)
            CASE DEFAULT
            END SELECT
            IF (errorstatus/=0) STOP !the energy and force calculation had error
            EnergyChangeImage(image)=EnergyChangeImage(image)+AL%PotentialEnergy
            AL=>AL%NextNeigh
         END DO
         
         !check if the elastic band is optimized by finding the perpendicular force
         !the AL%AtomForces need not be recalculated as they are already there
         !the tangents      need not be recalculated as they are already there
         
         TSImage=0
         TSEnergy=0._dp
         AL=>chos%AL
         DO image=2,NImages-1
            AL=>AL%NextNeigh
            IF (AL%PotentialEnergy-InitialEnergy>TSEnergy) THEN
               TSEnergy=AL%PotentialEnergy-InitialEnergy
               TSImage=image
            END IF
         END DO
         
         NotSatisfied=MAXVAL(ABS(EnergyChangeImage))>ftol
         WRITE(6,*) "Converged?", .NOT. NotSatisfied
         WRITE(6,*) "Max energy change:",MAXVAL(ABS(EnergyChangeImage))
         !RETURN
      END DO
      
      WRITE(6,*) "NEB has successfully converged"
      WRITE(6,*) "Number of iterations required:",NEBIteration
      
      AL=>chos%AL
      DO i=1,chos%nimages
         WRITE(6,*) AL%PotentialEnergy-chos%AL%PotentialEnergy 
         AL=>AL%NextNeigh
      END DO
      
   END SUBROUTINE OptimizeChOS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ChOSEnergy(x,errorstatus)
      !provides the energy of an image
      USE VARIABLE_TYPE
      IMPLICIT NONE
      INTEGER :: errorstatus,i,j,indx
      REAL(dp), DIMENSION(:), INTENT(IN) :: x
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord
      REAL(dp) :: ChOSEnergy
      
      ChOSEnergy=0._dp
         
      AtomCoord=>AL%AtomCoord
      j=(ActiveImage-2)*NOptimize
      DO i=1,NOptimize
         indx=AtomCoordIndx(i)
         AtomCoord(indx)=x(i)
IF (ActiveImage==3) THEN
WRITE(6,*) "Comparing the coordinates:",x(i),MovingAtomCoord(i+j),MovingAtomCoord(i)
STOP
END IF
      END DO
      CALL AddVerletListQuick(AL,NAtoms) !list generated at each iteration
      CALL GetForcesVL(AL,NAtoms,errorstatus) !get updated energy and force
      IF (errorstatus/=0) STOP
      
      IF (CI_NEB1 .AND. AL%PotentialEnergy>AL%PrevNeigh%PotentialEnergy .AND. &
           AL%PotentialEnergy>AL%NextNeigh%PotentialEnergy) THEN
         !do nothing -- since the TS cannot see the springs
         TSImage=ActiveImage !this will be required in ChOSForce
      ELSE
         t(1:3*NAtoms)=AtomCoord(1:3*NAtoms)-AL%PrevNeigh%AtomCoord(1:3*NAtoms)
         ChOSEnergy=DOT_PRODUCT(t(1:3*NAtoms),t(1:3*NAtoms))
         
         t(1:3*NAtoms)=AL%NextNeigh%AtomCoord(1:3*NAtoms)-AtomCoord(1:3*NAtoms)
         ChOSEnergy=ChOSEnergy+DOT_PRODUCT(t(1:3*NAtoms),t(1:3*NAtoms))
         
         ChOSEnergy=0.5_dp*SpringConst*ChOSEnergy
      END IF
         
      ChOSEnergy=ChOSEnergy+AL%PotentialEnergy
         
   END FUNCTION ChOSEnergy
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ChOSForce(x,errorstatus)
      USE VARIABLE_TYPE
      IMPLICIT NONE
      INTEGER :: errorstatus,image,i,indx
      REAL(dp), DIMENSION(:), INTENT(IN) :: x
      REAL(dp), DIMENSION(SIZE(x)) :: ChOSForce
      REAL(dp), DIMENSION(:), POINTER :: ALChOSForce
      
      errorstatus=0
      image=ActiveImage !if this image is the TS then it has already been declared as TSImage in ChOSEnergy
      
      CALL TangentNEB(AL)
      CALL ImageForce(AL,image)
      ALChOSForce=>AL%ChOSForce
      ChOSForce=0._dp
      DO i=1,NOptimize
         indx=AtomCoordIndx(i)
         ChOSForce(i)=-SclFac*ALChOSForce(indx)
      END DO

   END FUNCTION ChOSForce
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE TangentNEB(AL)
      !finds the tangent at AL assuming the energies have been updated
      !tangent at the first and last image are zero
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(:), POINTER :: ChOSTangent
      REAL(dp) :: Vprev,Vcurr,Vnext,dVmin,dVmax,normt,normc
      INTEGER :: i,r1,r2
      
      ChOSTangent=>AL%ChOSTangent

      !choose one of the following
      !INCLUDE "NEBTangent.f90"
      !INCLUDE "NEBTangent1.f90"
      INCLUDE "NEBTangentImproved.f90"
      
   END SUBROUTINE TangentNEB
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ImageForce(AL,image)
      !finds the force acting on AL
      IMPLICIT NONE
      INTEGER :: i,image,indx
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(:), POINTER :: AtomForce,ChOSTangent,ChOSForce
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomCoordNext,AtomCoordPrev
      REAL(dp) :: fdott,cdott,normt
      
      ChOSTangent=>AL%ChOSTangent
      AtomForce=>AL%AtomForce
      ChOSForce=>AL%ChOSForce
      ChOSForce=0._dp
      
      !perpendicular force
      fdott=DOT_PRODUCT(AtomForce(1:3*NAtoms),ChOSTangent(1:3*NAtoms))
      IF (CI_NEB1 .AND. TSImage==image) THEN
         DO i=1,NOptimize
            indx=AtomCoordIndx(i)
            ChOSForce(indx)=AtomForce(indx)-2._dp*fdott*ChOSTangent(indx)
         END DO
      ELSE
         DO i=1,NOptimize
            indx=AtomCoordIndx(i)
            ChOSForce(indx)=AtomForce(indx)-fdott*ChOSTangent(indx)
         END DO
      END IF

      AtomCoord=>AL%AtomCoord
      AtomCoordPrev=>AL%PrevNeigh%AtomCoord
      AtomCoordNext=>AL%NextNeigh%AtomCoord
      !choose one of the following
      !INCLUDE "NEBSpringForce1.f90"
      INCLUDE "NEBSpringForce2.f90"
      
   END SUBROUTINE ImageForce
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE OptimizeChOSGlobal(chos1,omode1,CI_NEB,ftol,gtol,xtol,ITMAX,errorstatus,ReuseVL)
      !given a chain of states containing nimages, the subroutine moves the
      !images to the minimum energy path
      IMPLICIT NONE
      TYPE(ChOSContainer), POINTER :: chos1
      INTEGER, OPTIONAL :: omode1
      INTEGER :: omode,i,j,iter,ITMAX,errorstatus,indx,image
      LOGICAL :: CI_NEB
      REAL(dp) :: ftol,gtol,xtol,fret
      REAL(dp), DIMENSION(:), POINTER :: x,AtomCoord
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      LOGICAL :: ReuseVL
      
      errorstatus=0
      
      !ReuseVL allows the NEB to use the LL and VL of the first image thus saving the memory requirements
      ReuseVL1=ReuseVL
      
      chos=>chos1
      IF (.NOT. ASSOCIATED(chos1)) THEN
         WRITE(6,*) "$Err>> ChOS is not associated"
         STOP
      END IF
      NImages=chos%nimages
      IF (NImages>MaxNImages) THEN
         WRITE(6,*) "Err>> Increase the MaxNImages in OptimizationNEB"
         STOP
      END IF
      IF (NImages<3) THEN
         WRITE(6,*) "$Err>> ChOS does contains very few images. NImages",NImages
         STOP
      END IF
      
      CI_NEB1=CI_NEB
      NAtoms=chos%AL%NAtoms
      SpringConst=chos%SpringConst

      IF (PRESENT(omode1)) THEN
         omode=omode1
      ELSE
         omode=3
      END IF
      
      WRITE(6,*) "Optim>> Performing global NEB calculation"
      
      !Find how many atoms are moving (the dimensionality of the optimization space per image)
      !total degrees of freedom are going to be NOptimize*(NImages-2)
      AL=>chos%AL !initial state
      AtomCoord=>AL%AtomCoord
      AtomIsMoving=>AL%AtomIsMoving
      NOptimize=0
      DO indx=1,3*NAtoms
         IF (AtomIsMoving(indx)) THEN
            NOptimize=NOptimize+1
            AtomCoordIndx(NOptimize)=indx
         END IF
      END DO

      !Add verlet list and compute the energies and forces
      CALL AddVerletList(AL,NAtoms=NAtoms)
      CALL GetForcesVL(AL,NAtoms,errorstatus)
      IF (errorstatus/=0) THEN !the atoms are overlapping -- bad chos provided
         errorstatus=3
         RETURN
      END IF
      
      DO image=2,NImages-1
         WRITE(*,*) "Image #",image
         AL=>AL%NextNeigh
         
         IF (ReuseVL1) THEN
            AL%VL=>chos%AL%VL
            AL%LL=>chos%AL%LL
         END IF
         CALL AddVerletList(AL,NAtoms=NAtoms)
         CALL GetForcesVL(AL,NAtoms,errorstatus)
         IF (errorstatus/=0) THEN !the atoms are overlapping -- bad chos provided
            errorstatus=3
            RETURN
         END IF
         AtomCoord=>AL%AtomCoord
         j=(image-2)*NOptimize
         DO i=1,NOptimize
            indx=AtomCoordIndx(i)
            MovingAtomCoord(i+j)=AtomCoord(indx)
         END DO
      END DO
      
      AL=>AL%NextNeigh !final state
      
      IF (ReuseVL1) THEN
         AL%VL=>chos%AL%VL
         AL%LL=>chos%AL%LL
      END IF
      
      CALL AddVerletList(AL,NAtoms=NAtoms)
      CALL GetForcesVL(AL,NAtoms,errorstatus)
      IF (errorstatus/=0) THEN !the atoms are overlapping -- bad chos provided
         errorstatus=3
         RETURN
      END IF
      
      SclFac=1._dp
      x=>MovingAtomCoord(1:(NImages-2)*NOptimize)

      SELECT CASE(omode)
      CASE(1) !steepest descent
         SclFac=0.01_dp !to be used for scaling forces
         !ftol=1.0e-9_dp !this should come from potential
         !gtol=2.0e-6_dp
         !xtol=1.0e-5_dp
         CALL SteepestDescent(x,ftol,gtol,xtol,iter,fret,ChOSEnergyGlobal,ChOSForceGlobal, &
            dx=0.0001_dp,ITMAX=ITMAX,IsPrint1=.TRUE.,errorstatus=errorstatus)
      CASE(2) !conjugate gradient
         SclFac=0.01_dp
         !ftol=1.0e-9_dp !this should come from potential
         !gtol=2.0e-6_dp
         !xtol=1.0e-5_dp
         CALL ConjugateGradient(x,ftol,gtol,xtol,iter,fret,ChOSEnergyGlobal,ChOSForceGlobal, &
            dx=0.0001_dp,ITMAX=ITMAX,IsPrint1=.TRUE.,errorstatus=errorstatus)
      CASE(3) !LBFGS
         !ftol=1.0e-9_dp !this should come from potential
         !gtol=2.0e-3_dp
         !xtol=1.0e-5_dp
         !ITMAX=100
         CALL Lbfgs(x,ftol,gtol,xtol,iter,fret,ChOSEnergyGlobal,ChOSForceGlobal, &
              IsPrint1=.TRUE.,ITMAX=ITMAX,errorstatus=errorstatus,memcorrections=4)
      CASE(4) !dynamical Euler
         !CALL 
      CASE(5) !Runge Kutta
      CASE(6) !QuickMin
      CASE(7) !FIRE
      CASE(8) !DFP
         !ftol=1.0e-9_dp !this should come from potential
         !gtol=2.0e-6_dp
         !xtol=1.0e-5_dp
         CALL DFP(x,ftol,gtol,xtol,iter,fret,ChOSEnergyGlobal,ChOSForceGlobal, &
            IsPrint1=.TRUE.,errorstatus=errorstatus)
      CASE DEFAULT
      END SELECT
      
      IF (errorstatus/=0) RETURN
      
      WRITE(6,*) "Global NEB has successfully converged"
      
      !print the MEP energies
      AL=>chos%AL
      DO i=1,NImages
         WRITE(6,*) AL%PotentialEnergy-chos%AL%PotentialEnergy 
         AL=>AL%NextNeigh
      END DO
      
   END SUBROUTINE OptimizeChOSGlobal
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ChOSEnergyGlobal(x,errorstatus)
      !provides the sum of energy of all images
      USE VARIABLE_TYPE
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), INTENT(IN) :: x
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomCoordNext,AtomCoordPrev
      REAL(dp) :: ChOSEnergyGlobal,InitialEnergy
      INTEGER :: i,j,errorstatus,image,indx
      
      !add the correct positions and label the highest energy point as the TS
      TSImage=0
      TSEnergy=0._dp
      AL=>chos%AL%NextNeigh
      InitialEnergy=chos%AL%PotentialEnergy
!IntelligentPrint=.TRUE.
      DO image=2,NImages-1
         AtomCoord=>AL%AtomCoord
         j=(image-2)*NOptimize
         DO i=1,NOptimize
            indx=AtomCoordIndx(i)
            AtomCoord(indx)=x(i+j)
         END DO
         
         IF (ReuseVL1) THEN
            AL%VL=>chos%AL%VL
            AL%LL=>chos%AL%LL
         END IF
 
         CALL AddVerletList(AL,NAtoms=AL%NAtoms) !list generated at each iteration
         CALL GetForcesVL(AL,NAtoms,errorstatus) !get updated energy and force
         
         IF (errorstatus/=0) RETURN
         !Currently a very simplistic approach for determining the TS is used
         !The image with the highest energy is termed as the TS
         !However, if there are two or more images that are close to each other - the defn of 
         !TS is not clear. Such cases need to be thought out.

         IF (AL%PotentialEnergy-InitialEnergy>TSEnergy) THEN !highest energy image is the TS
            TSEnergy=AL%PotentialEnergy-InitialEnergy
            TSImage=image
         END IF
         AL=>AL%NextNeigh
      END DO
      
      !t=PBCDistance(AL,x,AL%AtomCoord)
      !write(*,*) "max displacement:",maxval(abs(t))
      !WHERE (AL%AtomIsMoving .AND. ABS(t)<0.20_dp) AL%AtomCoord=x
      !WHERE (AL%AtomIsMoving .AND. ABS(t)>=0.20_dp) AL%AtomCoord=AL%AtomCoord+0.20_dp*t/ABS(t)
      
      ChOSEnergyGlobal=0._dp
      AL=>chos%AL%NextNeigh
      DO image=2,NImages-1
         !WRITE(6,*) "Im",image,AL%PotentialEnergy
         ChOSEnergyGlobal=ChOSEnergyGlobal+AL%PotentialEnergy
         AtomCoord=>AL%AtomCoord
         AtomCoordPrev=>AL%PrevNeigh%AtomCoord
         AtomCoordNext=>AL%NextNeigh%AtomCoord
         IF (CI_NEB1 .AND. TSImage==image) THEN
            !do nothing -- since the TS cannot see the springs
         ELSE
            t(1:3*NAtoms)=AtomCoord(1:3*NAtoms)-AtomCoordPrev(1:3*NAtoms) !note that t is not a unit vector
            ChOSEnergyGlobal=ChOSEnergyGlobal+0.5_dp*SpringConst*DOT_PRODUCT(t(1:3*NAtoms),t(1:3*NAtoms))
         
            t(1:3*NAtoms)=AtomCoordNext(1:3*NAtoms)-AtomCoord(1:3*NAtoms)
            ChOSEnergyGlobal=ChOSEnergyGlobal+0.5_dp*SpringConst*DOT_PRODUCT(t(1:3*NAtoms),t(1:3*NAtoms))
         END IF
         AL=>AL%NextNeigh
      END DO
      
   END FUNCTION ChOSEnergyGlobal
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ChOSForceGlobal(x,errorstatus)
      USE VARIABLE_TYPE
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), INTENT(IN) :: x
      REAL(dp), DIMENSION(SIZE(x)) :: ChOSForceGlobal
      REAL(dp), DIMENSION(:), POINTER :: ChOSTangent,ALChOSForce
      INTEGER :: image,errorstatus,i,j,indx
      
      errorstatus=0
      
      ChOSTangent=>AL%ChOSTangent
      ChOSForceGlobal=0._dp

      AL=>chos%AL%NextNeigh
      DO image=2,NImages-1
         CALL TangentNEB(AL)
         CALL ImageForce(AL,image)
         ALChOSForce=>AL%ChOSForce
         j=(image-2)*NOptimize
         DO i=1,NOptimize
            indx=AtomCoordIndx(i)
            ChOSForceGlobal(i+j)=-SclFac*ALChOSForce(indx) !to get the actual energy gradient
         END DO
         AL=>AL%NextNeigh
      END DO

   END FUNCTION ChOSForceGlobal
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE OptimizationNEB
