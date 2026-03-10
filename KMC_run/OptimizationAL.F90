MODULE OptimizationAL
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
   
   TYPE(SystemContainer), POINTER, PRIVATE :: AL1,ALa,ALb

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !TYPE(SystemContainer), POINTER, PRIVATE :: ALCopy
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   INTEGER, PRIVATE :: iprint
   REAL(dp), PRIVATE :: SclFac,MaxDisplacement1
   REAL(dp), DIMENSION(3*MaxNAtoms), TARGET, PRIVATE :: MovingAtomCoord,AtomCoordCopy
   REAL(dp), DIMENSION(3*MaxNAtoms), PRIVATE :: t
   INTEGER, DIMENSION(3*MaxNAtoms), PRIVATE :: AtomCoordLocation
   INTEGER, PRIVATE :: NOptimize,NAtoms1
   REAL(dp), DIMENSION(:), POINTER, PRIVATE :: AtomCoord,AtomForce
   LOGICAL, DIMENSION(:), POINTER, PRIVATE :: AtomIsMoving
   LOGICAL, PRIVATE :: PrintXYZOptimized,ByPassStop
   INTEGER, DIMENSION(:), POINTER, PRIVATE :: NearestMatch
   REAL(dp), DIMENSION(:), POINTER, PRIVATE :: NearestMatchDistance
   LOGICAL, DIMENSION(:), POINTER :: collision
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE OptimizeAL(AL,omode1,NAtoms,ftol,gtol,xtol,ITMAX,MaxDisplacement,FileXYZOptimized,errorstatus)
      !minimize energy of AL by moving in the direction of the forces
      !FileXYZOptimized is used to print the atom positions during optimization
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER, OPTIONAL :: omode1
      INTEGER :: omode,iter,errorstatus1,ITMAX
      INTEGER, OPTIONAL :: errorstatus
      INTEGER :: NAtoms,loc
      REAL(dp) :: ftol,gtol,xtol,fret,MaxDisplacement
      REAL(dp), DIMENSION(:), POINTER :: x
      CHARACTER(len=100), OPTIONAL :: FileXYZOptimized
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Bit of housekeeping
      AL1=>AL
      AtomCoordCopy(1:3*NAtoms)=AL%AtomCoord(1:3*NAtoms)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Created temprorailiy - this can be deleted
      !NULLIFY(ALCopy)
      !CALL MakeSize(ALCopy,AL%NAtoms)
      !CALL Copy(AL,ALCopy)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      NAtoms1=NAtoms
      MaxDisplacement1=MaxDisplacement
      CALL AddVerletList(AL1,NAtoms=NAtoms1,ListType=AL%VL%ListType)
      
      IF (PRESENT(omode1)) THEN
         omode=omode1
      ELSE
         omode=3
      END IF
      
      !iprint=10
      
      ByPassStop=.FALSE.
      IF (PRESENT(errorstatus)) BypassStop=.TRUE.
      
      AtomCoord=>AL1%AtomCoord
      AtomForce=>AL1%AtomForce
      AtomIsMoving=>AL1%AtomIsMoving

      PrintXYZOptimized=.FALSE.
      IF (PRESENT(FileXYZOptimized)) THEN
         OPEN(UNIT=363,FILE=FileXYZOptimized)
         PrintXYZOptimized=.TRUE.
         CALL WriteXYZ(AL1,fileopen=.FALSE.,fileclose=.FALSE.,iunit=363)
      END IF
      
      NOptimize=0 !this is the dimensionality of optimization space
      DO loc=1,3*NAtoms !all coordinate dimensions
         IF (AtomIsMoving(loc)) THEN
            NOptimize=NOptimize+1
            MovingAtomCoord(NOptimize)=AtomCoord(loc) !this is the array to be optimized
            AtomCoordLocation(NOptimize)=loc !this is not atom index
         END IF
      END DO
      t(1:NOptimize)=MovingAtomCoord(1:NOptimize) !t stores the old state
      x=>MovingAtomCoord(1:NOptimize)
      
      iter=0
      SclFac=1._dp
!From Voter:
!ftol=1.e-9 is good for eam - keep these tight to prevent balancing on saddle
!gtol=1.e-5 is good for eam
!xtol=1.e-5 is good for eam
      SELECT CASE(omode)
      CASE(1) !steepest descent
         SclFac=0.01_dp !to be used for scaling forces
         !ftol=1.0e-9_dp !this should come from potential
         !gtol=2.0e-6_dp
         !xtol=1.0e-5_dp
         CALL SteepestDescent(x,ftol,gtol,xtol,iter,fret,ALEnergy,ALForce, &
            dx=0.0001_dp,ITMAX=ITMAX,IsPrint1=.FALSE.,errorstatus=errorstatus1)
      CASE(2) !conjugate gradient
         SclFac=0.01_dp
         !ftol=1.0e-9_dp !this should come from potential
         !gtol=2.0e-6_dp
         !xtol=1.0e-5_dp
         CALL ConjugateGradient(x,ftol,gtol,xtol,iter,fret,ALEnergy,ALForce, &
            dx=0.0001_dp,ITMAX=ITMAX,IsPrint1=.FALSE.,errorstatus=errorstatus1)
      CASE(3) !LBFGS
         !SclFac=0.01_dp
         !ftol=1.0e-9_dp !this should come from potential
         !gtol=2.0e-4_dp
         !xtol=1.0e-5_dp
         CALL Lbfgs(x,ftol,gtol,xtol,iter,fret,ALEnergy,ALForce,IsPrint1=.FALSE., &
            ITMAX=ITMAX,errorstatus=errorstatus1)
      CASE(4) !dynamical Euler
         !CALL 
      CASE(5) !Runge Kutta
      CASE(6) !QuickMin
      CASE(7) !FIRE
      CASE(8) !DFP
         SclFac=0.01_dp
         !ftol=1.0e-9_dp !this should come from potential
         !gtol=2.0e-6_dp
         !xtol=1.0e-5_dp
         CALL DFP(x,ftol,gtol,xtol,iter,fret,ALEnergy,ALForce,IsPrint1=.TRUE., &
            errorstatus=errorstatus1)
      !CASE(9) !mix of steepest-descent and conjugate gradient
      !   CALL SteepestDescent(x,ftol,gtol,xtol,iter,fret,ALEnergy,ALForce, &
      !      dx=0.0001_dp,ITMAX=4,IsPrint1=.TRUE.,errorstatus=errorstatus1)
      !   CALL ConjugateGradient(x,ftol,gtol,xtol,iter,fret,ALEnergy,ALForce, &
      !      dx=0.0001_dp,IsPrint1=.TRUE.,errorstatus=errorstatus1)
      !CASE(10) !mix of steepest-descent and L-BFGS
      !   CALL SteepestDescent(x,ftol,gtol,xtol,iter,fret,ALEnergy,ALForce, &
      !      dx=0.0001_dp,ITMAX=4,IsPrint1=.TRUE.,errorstatus=errorstatus1)
      !   CALL Lbfgs(x,ftol,gtol,xtol,iter,fret,ALEnergy,ALForce,IsPrint1=.TRUE., &
      !      ITMAX=ITMAX,errorstatus=errorstatus1)
      CASE DEFAULT
      END SELECT
      
      IF (PRESENT(errorstatus)) errorstatus=errorstatus1
      IF (PrintXYZOptimized) CLOSE(363)
      
      
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !CALL Delete(ALCopy)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE OptimizeAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ALEnergy(x,errorstatus)
      USE VARIABLE_TYPE
      IMPLICIT NONE
      INTEGER :: errorstatus,indx,loc
      REAL(dp), DIMENSION(:), INTENT(IN) :: x
      REAL(dp) :: ALEnergy,dx1,dx2,coordx,coordt,maxdx
      CHARACTER(len=100) :: filename
      
      !Originally was working well with 0.05
      !using 0.03 slows down the code -- however, it prevents code from crashing for the example that was tried
      !using L-BFGS. It appeared that L-BFGS, when it saw significant reduction in energy using large displacements,
      !tried to bring atoms closer than 1 Ang which made to code stop in CALL GetVerletForce
      
      DO indx=1,NOptimize
         coordx=x(indx) !new coordinate - non-PBC coordinate
         coordt=t(indx) !old coordinate
         dx1=coordx-coordt !displacement occurred
         dx2=MIN(ABS(dx1),MaxDisplacement1) !this prevents the atom to be displaced more than MaxDisplacement1
         IF (dx2>0._dp) THEN
            !dx1=dx2*dx1/ABS(dx1) !dx2 is a positive number
            coordt=coordt+dx1
            t(indx)=coordt
            MovingAtomCoord(indx)=coordt
            loc=AtomCoordLocation(indx)
            AtomCoord(loc)=coordt
         END IF
      END DO
      !IF (iprint==10) OPEN(UNIT=312,FILE="optimization.xyz")
      !IF (iprint==10) CALL WriteXYZ(AL,fileopen=.FALSE.,fileclose=.FALSE.,iunit=312)
      !IF (iprint==10) CLOSE(312)
      
      CALL AddVerletListQuick(AL1,NAtoms1)
      CALL GetForcesVL(AL1,NAtoms1,errorstatus)
      
      IF (PrintXYZOptimized) &
         CALL WriteXYZ(AL1,fileopen=.FALSE.,fileclose=.FALSE.,iunit=363)
      
      IF (errorstatus/=0) THEN
         WRITE(6,*) "Err>> Printing the atom configuration to OptimizationError.xyz"
         AL1%AtomCoord(1:3*AL1%NAtoms)=AtomCoordCopy(1:3*AL1%NAtoms)
         filename="OptimizationError.xyz"
         CALL WriteXYZ(AL1,filename=filename)
         IF (PrintXYZOptimized) CLOSE(363)
         IF (.NOT. ByPassStop) STOP
      END IF
      ALEnergy=AL1%PotentialEnergy
      !WRITE(*,*) ALEnergy
   END FUNCTION ALEnergy
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ALForce(x,errorstatus)
      USE VARIABLE_TYPE
      IMPLICIT NONE
      INTEGER :: errorstatus,indx,loc,rn1,rn2
      REAL(dp), DIMENSION(:), INTENT(IN) :: x
      REAL(dp), DIMENSION(SIZE(x)) :: ALForce
      
      DO indx=1,NOptimize
         loc=AtomCoordLocation(indx)
         ALForce(indx)=-SclFac*AtomForce(loc) !to get actual energy gradient
      END DO
   END FUNCTION ALForce
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !SUBROUTINES TO CONVERT AN OFF-LATTICE STRUCTURE TO A LATTICE STRUCTURE
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Translate2Align(ALa1,ALb1,x,v,matcha2b,matchb2a,matchseparation)
   !move AL1 by an amount t so that it matches best with AL
   !The main idea is that while finding the lattice sites
   ! we wish to match where atoms are present with respect
   ! to some reference lattice structure. First we need to have
   ! a god match between the two. To do so we translate the
   ! reference to get a good match
   !v is the repeat unit of the lattice
   !x is the final amount by which ALb1 should be translated
   !matcha2b is the index in b that matches with a
   !matchb2a is the index in a that matches with b
   !matchseparation is the distance between a and b with indices given in matcha2b
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: ALa1,ALb1
      REAL(dp) :: dx,sm1,fret,ftol,gtol,xtol,x(3),v(3),xmin(3)
      REAL(dp), DIMENSION(:), POINTER :: freq1,freq2
      INTEGER :: iter,errorstatus1,ITMAX,dimn
      INTEGER, DIMENSION(:), POINTER, OPTIONAL :: matcha2b,matchb2a
      REAL(dp), DIMENSION(:), POINTER, OPTIONAL :: matchseparation
      LOGICAL :: IntelliCheck
      

      ALa=>ALa1
      ALb=>ALb1 !this should be reference state
      ALLOCATE(NearestMatch(ALa%NAtoms),NearestMatchDistance(ALa%NAtoms))
      ALLOCATE(collision(ALa%NAtoms))

      x=0._dp
      !optimize along x,y  and z
      DO dimn=1,3
         dx=0.4; CALL LocateMin(dimn=dimn,xrange=(/v(dimn),0._dp/),x=x,dx=dx,xmin=x)
         dx=0.2; CALL LocateMin(dimn=dimn,xrange=(/0.4_dp,-0.4_dp/),x=x,dx=dx,xmin=x)
         dx=0.05; CALL LocateMin(dimn=dimn,xrange=(/0.2_dp,-0.2_dp/),x=x,dx=dx,xmin=x)
         dx=0.01; CALL LocateMin(dimn=dimn,xrange=(/0.05_dp,-0.05_dp/),x=x,dx=dx,xmin=x)
      END DO
      IntelliCheck=IntelligentPrint
      IntelligentPrint=.TRUE.
      sm1=MapDistance(x,errorstatus1)
      IF (.NOT. IntelliCheck) IntelligentPrint=.FALSE.
      WRITE(*,*) "Translate by ",x

      IF (PRESENT(matcha2b)) THEN
         matcha2b=NearestMatch
      END IF
      IF (PRESENT(matchb2a)) THEN
         matchb2a=0
         DO iter=1,ALa1%NAtoms
            IF (NearestMatch(iter)>0) matchb2a(NearestMatch(iter))=iter
         END DO
      END IF
      DEALLOCATE(NearestMatch)
 
      IF (PRESENT(matchseparation)) THEN
         matchseparation=NearestMatchDistance
      END IF
      DEALLOCATE(NearestMatchDistance)
      
      DEALLOCATE(collision)
      RETURN

      x=0._dp
      ftol=1.e-3
      gtol=1.e-4
      xtol=5.e-2
      ITMAX=55555
      !CALL SteepestDescent(x,ftol,gtol,xtol,iter,fret,MapDistance,MapDistanceGrad, &
      !    dx=0.0001_dp,ITMAX=ITMAX,IsPrint1=.TRUE.,errorstatus=errorstatus1)
      !CALL ConjugateGradient(x,ftol,gtol,xtol,iter,fret,MapDistance,MapDistanceGrad, &
      !     dx=0.0001_dp,ITMAX=ITMAX,IsPrint1=.FALSE.,errorstatus=errorstatus1)
      !CALL Lbfgs(x,ftol,gtol,xtol,iter,fret,MapDistance,MapDistanceGrad,IsPrint1=.FALSE., &
      !    ITMAX=ITMAX,errorstatus=errorstatus1)
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE LocateMin(dimn,xrange,x,dx,xmin)
         IMPLICIT NONE
         INTEGER :: dimn !this is the dimension where we optimize
         REAL(dp) :: x(3),t(3) !this is the input value of x or translation vector
         REAL(dp) :: dx !this is by how much we need to increment x
         REAL(dp) :: xrange(2) !this is maximum and minimum increment to x
         REAL(dp) :: xmin(3) !this is the value of minimum found
         REAL(dp) :: val(1000),xloc(1000) !hold the value and location along x(dimn)
         INTEGER :: errorstatus,i,npts
         REAL(dp) :: minval

         t=x
         t(dimn)=t(dimn)+xrange(2)
         npts=CEILING((xrange(1)-xrange(2))/dx)+1
         minval=100000000._dp
         DO i=1,npts
            xloc(i)=t(dimn)
            val(i)=MapDistance(t,errorstatus)
            IF (val(i)<minval) THEN
               xmin=t
               minval=val(i)
            END IF
            !WRITE(*,*) i,t,val(i)
            t(dimn)=t(dimn)+dx
         END DO
         WRITE(*,*) "Minimum found at ",xmin," minimum value ",minval
      END SUBROUTINE LocateMin
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE Translate2Align
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION MapDistance(x,errorstatus)
   !THIS IS A CUSTOMIZED FUNCTION GO THROUGH THE FUNCTION CAREFULLY
     USE VARIABLE_TYPE
     IMPLICIT NONE
     INTEGER :: errorstatus
     REAL(dp), DIMENSION(:), INTENT(IN) :: x
     REAL(dp) :: MapDistance
     REAL(dp) :: d,cutoff,chg1,chg2
     INTEGER :: iatom1,iatom2,count,N0_5,N5_1,N1_5,N5
     INTEGER, PARAMETER :: chunksize=20

     errorstatus=0
     DO iatom2=1,ALb%NAtoms
        ALb%AtomCoord(3*iatom2-2:3*iatom2)=ALb%AtomCoord(3*iatom2-2:3*iatom2)+x
     END DO

     collision=.FALSE.
!$OMP PARALLEL
!$OMP DO SCHEDULE(DYNAMIC,chunksize)
     DO iatom1=1,ALa%NAtoms
        CALL MapDistanceWorker(iatom1,NearestMatch,NearestMatchDistance)
     END DO
!$OMP END DO NOWAIT
!$OMP BARRIER
     !check for collisions
!$OMP DO SCHEDULE(DYNAMIC,chunksize)
     DO iatom1=1,ALa%NAtoms
        DO iatom2=1,ALa%NAtoms
           collision(iatom1)=NearestMatch(iatom1)==NearestMatch(iatom2) .AND. &
              iatom1/=iatom2
           IF (collision(iatom1)) EXIT
        END DO
     END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
     !Check for collisions
     MapDistance=0._dp
     IF (ANY(collision)) THEN
        WRITE(6,*) "Collision occurred ... penalizing map"
        MapDistance=1000._dp
     END IF

     DO iatom1=1,ALa%Natoms
        !chg1=ALa%AtomCharge(iatom1)
        !IF (NearestMatchDistance(iatom1)<1000._dp .AND. chg1>0._dp)  &
          MapDistance=MapDistance+NearestMatchDistance(iatom1)
     END DO
     IF (IntelligentPrint) THEN
         N0_5=0 !number of points between 0 and 0.5
         N5_1=0  !between 0.5 and 1
         N1_5=0  !between 1 and 1.5
         N5=0 !beyond
         write(*,*) "Match found"
         DO iatom1=1,ALa%NAtoms
            !iatom2=NearestMatch(iatom1)
            !IF (iatom2==0) THEN
            !    WRITE(UNIT=6,FMT='(2I5)') iatom1,iatom2
            !ELSE
            !    WRITE(UNIT=6,FMT='("Index: ",2I7,3ES15.5,"    <-->   ",3ES15.5)') &
            !       iatom1,iatom2,ALa%AtomCoord(3*iatom1-2:3*iatom1), &
            !       ALb%AtomCoord(3*iatom2-2:3*iatom2)
            !END IF
            d=NearestMatchDistance(iatom1)
            IF (d<0.5_dp) THEN
               N0_5=N0_5+1
            ELSEIF (d<1._dp) THEN
               N5_1=N5_1+1
            ELSEIF (d<1.5_dp) THEN
               N1_5=N1_5+1
            ELSE
               N5=N5+1
            END IF
         END DO
         WRITE(*,*) "Number of atoms with match found in range [<0.5]:",N0_5
         WRITE(*,*) "Number of atoms with match found in range [0.5-1]:",N5_1
         WRITE(*,*) "Number of atoms with match found in range [1-1.5]:",N1_5
         WRITE(*,*) "Number of atoms with match found in range [>1.5]:",N5
     END IF
     !write(*,*) MapDistance,x,MAXVAL(NearestMatchDistance)
     !Set the coordinates back to original value
     DO iatom2=1,ALb%NAtoms
        ALb%AtomCoord(3*iatom2-2:3*iatom2)=ALb%AtomCoord(3*iatom2-2:3*iatom2)-x
     END DO
   END FUNCTION MapDistance
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MapDistanceWorker(iatom1,NearestMatch,NearestMatchDistance)
   !THIS IS A CUSTOMIZED FUNCTION GO THROUGH THE FUNCTION CAREFULLY
     USE VARIABLE_TYPE
     IMPLICIT NONE
     INTEGER :: errorstatus
     REAL(dp) :: d,coord1(3),coord2(3),chg1,chg2
     REAL(dp), DIMENSION(:) :: NearestMatchDistance
     INTEGER :: iatom1,iatom2,count
     INTEGER, DIMENSION(:) :: NearestMatch

        coord1=ALa%AtomCoord(3*iatom1-2:3*iatom1)
        chg1=ALa%AtomCharge(iatom1)
        NearestMatchDistance(iatom1)=1000._dp
        DO iatom2=1,ALb%NAtoms
           coord2=ALb%AtomCoord(3*iatom2-2:3*iatom2)
           coord2=PBCdistance(ALa,coord1,coord2)
           coord2=coord2*coord2
           d=SQRT(coord2(1)+coord2(2)+coord2(3))
           chg2=ALb%AtomCharge(iatom2)
           IF (d<=NearestMatchDistance(iatom1) .AND. chg1*chg2>=0._dp) THEN !chg1 and chg2 should be same sign
              NearestMatchDistance(iatom1)=d
              NearestMatch(iatom1)=iatom2
           END IF
        END DO
   END SUBROUTINE MapDistanceWorker
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION MapDistanceGrad(x,errorstatus)
      USE VARIABLE_TYPE
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(dp), DIMENSION(:), INTENT(IN) :: x
      REAL(dp), DIMENSION(SIZE(x)) :: MapDistanceGrad,tx
      REAL(dp) :: fp,fn

      MapDistanceGrad=0._dp
      tx=x
      tx(1)=x(1)+0.01_dp
      fp=MapDistance(tx,errorstatus)
      tx(1)=x(1)-0.01_dp
      fn=MapDistance(tx,errorstatus)
      MapDistanceGrad(1)=(fp-fn)/0.02_dp

      tx=x
      tx(2)=x(2)+0.01_dp
      fp=MapDistance(tx,errorstatus)
      tx(2)=x(2)-0.01_dp
      fn=MapDistance(tx,errorstatus)
      MapDistanceGrad(2)=(fp-fn)/0.02_dp
      
      tx=x
      t(3)=x(3)+0.01_dp
      fp=MapDistance(tx,errorstatus)
      tx(3)=x(3)-0.01_dp
      fn=MapDistance(tx,errorstatus)
      MapDistanceGrad(3)=(fp-fn)/0.02_dp
   END FUNCTION MapDistanceGrad
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE OptimizationAL
