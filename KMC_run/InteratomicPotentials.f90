MODULE PotentialPackage
   !Master module for performing energy based calculations

   USE VARIABLE_TYPE
   USE IO
   USE db_manipulate
   USE NeighborList
   USE EAMInitialize
   USE EAMSubroutines
   USE LJInitialize 
   USE LJSubroutines
   USE CoulombInitialize
   USE CoulombSubroutines
   USE BuckinghamInitialize
   USE BuckinghamSubroutines
   USE StillingerWeberInitialize
   USE StillingerWeberSubroutines
   USE TersoffInitialize
   USE TersoffSubroutines
   USE DFTSubroutines
   IMPLICIT NONE
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PotentialInitialize(AL,IsPrint)
      !initialize potentials - this has to be called for each new AL
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(InteracPotential), POINTER :: Potential
      !do not deallocate AL%Potential unless sure ... 
      !   ...other ALs might be using this
      LOGICAL, OPTIONAL :: IsPrint
      LOGICAL :: IsPrint1
      
      IF (PRESENT(IsPrint)) THEN
         IsPrint1=IsPrint
      ELSE
         IsPrint1=.TRUE.
      END IF
      
      IF (ASSOCIATED(AL%Potential)) THEN
         IF (.NOT. AL%Potential%UnInitialized) RETURN
      END IF
      
      IF (IsPrint1) WRITE(6,*) "Ini>> Initializing all interatomic potentials ..."
      ALLOCATE(AL%Potential)
      Potential=>AL%Potential

      CALL ReadPotentialInfo(AL)
         
      IF (ASSOCIATED(Potential%EAM)) CALL EAMSetUp(AL,IsPrint1)
      IF (ASSOCIATED(Potential%LJ)) CALL LJSetUp(AL)
      IF (ASSOCIATED(Potential%Buckingham)) CALL BuckinghamSetUp(AL)
      IF (ASSOCIATED(Potential%SW)) CALL SWSetUp(AL)
      IF (ASSOCIATED(Potential%Tersoff)) CALL TersoffSetUp(AL)
      !CALL BOSetUp(AL)
      
      !Coulomb potential requires neighborlist, while neighborlist needs
      !all interactions, so we send Coulomb to the end of potential initialize
      IF (ASSOCIATED(Potential%Coulomb)) THEN
         IF (ASSOCIATED(Potential%Coulomb%Ewald)) CALL EwaldSetUp(AL)
         IF (ASSOCIATED(Potential%Coulomb%Wolf)) CALL WolfSetUp(AL)
         WRITE(6,*) "Coulomb potential is now available"
      END IF

      !IF (ASSOCIATED(Potential%DFT)) CALL 
      
      Potential%UnInitialized=.FALSE.

   END SUBROUTINE PotentialInitialize
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetEnergyVL(AL,NAtoms)
      !VL should already be filled
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(InteracPotential), POINTER :: Potential
      INTEGER :: NAtoms
      
      CALL PotentialInitialize(AL)

      Potential=>AL%Potential
      AL%PotentialEnergy=0._dp !initialize
      
      IF (ASSOCIATED(Potential%EAM)) CALL EAMEnTaylor(AL,NAtoms)
      IF (ASSOCIATED(Potential%LJ)) CALL LJEnTaylor(AL,NAtoms)
      IF (ASSOCIATED(Potential%Coulomb)) THEN
         IF (ASSOCIATED(Potential%Coulomb%Ewald)) CALL CoulombEnEwald(AL,NAtoms)
         IF (ASSOCIATED(Potential%Coulomb%Wolf)) CALL CoulombEnWolf(AL,NAtoms)
      END IF
      IF (ASSOCIATED(Potential%Buckingham)) CALL BuckinghamEnTaylor(AL,NAtoms)
      IF (ASSOCIATED(Potential%SW)) CALL SWEnTaylor(AL,NAtoms)
      IF (ASSOCIATED(Potential%Tersoff)) CALL TersoffEnTaylor(AL,NAtoms)

   END SUBROUTINE GetEnergyVL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetForcesVL(AL,NAtoms,errorstatus)
      !VL should already be filled
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(InteracPotential), POINTER :: Potential
      REAL(dp), DIMENSION(:), POINTER :: AtomForce
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      INTEGER :: errorstatus,errorstatus1,NAtoms
      
      CALL PotentialInitialize(AL)

      Potential=>AL%Potential
      AL%PotentialEnergy=0._dp !initialize
      AL%AtomForce=0._dp
      
      !get all forces and energies
      errorstatus=0
      errorstatus1=0
      IF (ASSOCIATED(Potential%EAM)) CALL EAMForceLinearInterp(AL,NAtoms,errorstatus1)
      errorstatus=MAX(errorstatus,errorstatus1)
      IF (ASSOCIATED(Potential%LJ)) CALL LJForceLinearInterp(AL,NAtoms,errorstatus1)
      !IF (ASSOCIATED(Potential%Buckingham)) CALL BuckinghamForceLinearInterp(AL,NAtoms,errorstatus1)
      !errorstatus=MAX(errorstatus,errorstatus1)
      IF (ASSOCIATED(Potential%Coulomb)) THEN
         IF (ASSOCIATED(Potential%Coulomb%Ewald)) CALL CoulombForceEwald(AL,NAtoms,errorstatus1)
         errorstatus=MAX(errorstatus,errorstatus1)
      END IF
      IF (ASSOCIATED(Potential%SW)) CALL SWForceLinearInterp(AL,NAtoms,errorstatus1)
      IF (ASSOCIATED(Potential%Tersoff)) CALL TersoffForceLinearInterp(AL,NAtoms,errorstatus1)
      
      IF (ASSOCIATED(Potential%DFT)) CALL DFTForces(AL)
      
      AtomIsMoving=>AL%AtomIsMoving
      AtomForce=>AL%AtomForce
      WHERE (.NOT. AtomIsMoving(1:3*NAtoms)) AtomForce(1:3*NAtoms)=0._dp !makes the force on non-moving dimensions as zero

   END SUBROUTINE GetForcesVL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetHessianVL(AL,NAtoms,curvature,ndimn,AtomMap,nmove)
   !finds the scaled Hessian matrix for use with normal mode analysis
   !curvature contains the scaled Hessian matrix terms
   !ndim will tell the number of elements in curvature
   !AtomMap contains the atom mapping -- see below
   !nmove is number of moving atoms
      IMPLICIT NONE
      INTEGER, PARAMETER :: maxnmove=1000 !size of curvature is 3maxnmove*(3maxnmove+1)
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION( (3*maxnmove*(3*maxnmove+1))/2 ) :: curvature
      REAL(dp) :: hess,ep,e0,en
      REAL(dp), DIMENSION(:), ALLOCATABLE :: fp,fn
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomForce,AtomMass
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      REAL(dp) :: displacement,coord,mass_iatom,mass_jatom,effectivemass
      INTEGER, DIMENSION(maxnmove) :: AtomMap
      INTEGER :: NAtoms,ndimn,NMove,imove,jmove,iatom,jatom
      INTEGER :: dir,dir1,idx,jdx,idy,jdy,ijy,errorstatus
      
      CALL PotentialInitialize(AL)
      CALL AddVerletListQuick(AL,NAtoms) !initialize the VL
      
      !Collect the moving atoms and create AtomMap
      Nmove=0
      idx=-2
      AtomCoord=>AL%AtomCoord
      AtomForce=>AL%AtomForce
      AtomIsMoving=>AL%AtomIsMoving
      AtomMass=>AL%AtomMass
      DO iatom=1,NAtoms
         idx=idx+3
         IF (ANY(AtomIsMoving(idx:idx+2))) THEN
            Nmove=Nmove+1
            AtomMap(Nmove)=iatom !this indexing is very important later
         END IF
      END DO
      
      !Check the number of dimensions present
      IF (maxnmove<Nmove) THEN
         WRITE(6,*) "Err>> Curvature array is small"
         STOP
      END IF
      ndimn=3*Nmove
      ndimn=(ndimn*(ndimn+1))/2
      
      !initialize and allocate
      curvature(1:ndimn)=0._dp
      ALLOCATE(fp(3*NAtoms))
      ALLOCATE(fn(3*NAtoms))

      !find the curvature
      displacement=0.0001_dp !this depends on the potential type, but we
        !are going to assume that the potential is able to resolve this 
        !difference
      
      DO imove=1,Nmove
         iatom=AtomMap(imove) !the ith atom
         DO dir=1,3
            idx=3*(iatom-1)+dir
            mass_iatom=AtomMass(idx)
            coord=AtomCoord(idx)
            idy=3*(imove-1)+dir
            
            !step 1
            AtomCoord(idx)=coord+displacement !positive displacement
            !CALL UpdateOldVL(AL,NAtoms,iatom,errorstatus) !update the Verlet neighborlist
            !IF (errorstatus/=0) 
            CALL AddVerletListQuick(AL,NAtoms)
            CALL GetForcesVL(AL,NAtoms,errorstatus)
            fp(1:3*NAtoms)=AtomForce(1:3*NAtoms)
            ep=AL%PotentialEnergy
            
            !step 2
            AtomCoord(idx)=coord-displacement !negative displacement
            !CALL UpdateOldVL(AL,NAtoms,iatom,errorstatus) !update the neighborlist
            !IF (errorstatus/=0) 
            CALL AddVerletListQuick(AL,NAtoms)
            CALL GetForcesVL(AL,NAtoms,errorstatus)
            fn(1:3*NAtoms)=AtomForce(1:3*NAtoms)
            en=AL%PotentialEnergy
            
            !step 3
            AtomCoord(idx)=coord !move to original location
            !CALL UpdateOldVL(AL,NAtoms,iatom,errorstatus) !update the neighborlist
            !IF (errorstatus/=0) 
            CALL AddVerletListQuick(AL,NAtoms)
            CALL GetForcesVL(AL,NAtoms,errorstatus)
            e0=AL%PotentialEnergy
           
            DO jmove=1,NMove
               jatom=AtomMap(jmove)
               DO dir1=1,3
                  jdx=3*(jatom-1)+dir1
                  mass_jatom=AtomMass(jdx)
                  jdy=3*(jmove-1)+dir1
                  ijy=(MAX(idy,jdy)*( MAX(idy,jdy) - 1 ))/2+ MIN(idy,jdy) !location in curvature matrix (upper diag matrix from NETLIB)
                  !ijx=(MAX(idx,jdx)*( MAX(idx,jdx) - 1 ))/2+ MIN(idx,jdx) !location in curvature matrix (upper diag matrix from NETLIB)
                  
                  effectivemass=SQRT(mass_iatom*mass_jatom)
                  
                  hess=(fn(jdx)-fp(jdx))/effectivemass !*0.5_dp/displacement --- this will be multiplied later
                  !IF (hess<0._dp) THEN
                  !   WRITE(6,*) "hess is negative",hess,iatom,jatom,dir,dir1
                     !WRITE(6,*) ep,e0,en
                     !stop
                  !ELSEIF (hess>0._dp) THEN
                  !   WRITE(6,*) "hess is positive",fn(jdx),fp(jdx)
                  !END IF
                  
                  
                  IF (idx==jdx) THEN
                     curvature(ijy)=hess
                  ELSE
                     curvature(ijy)=curvature(ijy)+0.5_dp*hess !two contributions
                  END IF
                  
                  !if(iatom==1) then
                  !   write(*,*) -fp(jdx)/27.21_dp,-fn(jdx)/27.21_dp
                  !end if
               END DO
            END DO
         END DO
      END DO
      
      !finalize and deallocate
      curvature(1:ndimn)=curvature(1:ndimn)*(0.5_dp/displacement)
      DEALLOCATE(fp)
      DEALLOCATE(fn)
      
   END SUBROUTINE GetHessianVL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ChkForces(AL)
   
      !either VL or LL based - currently VL based
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: i,j,NAtoms,errorstatus
      REAL(dp) :: dr,df,en,ep,t1,t2,maxerr
      LOGICAL :: IsError
      
      WRITE(*,*) "Check>> Ensuring forces are correctly calculated"
      NAtoms=AL%NAtoms
      
      CALL AddVerletList(AL,NAtoms=NAtoms)
      CALL GetForcesVL(AL,NAtoms,errorstatus)
      IF (errorstatus/=0) STOP

      dr=0.001_dp
      
      CALL CPU_TIME(t1)
      CALL AddVerletList(AL,.TRUE.,NAtoms=NAtoms)
      DO i=1,100
         CALL GetForcesVL(AL,NAtoms,errorstatus) !get the forces
      END DO
      CALL CPU_TIME(t2)
      WRITE(UNIT=*,FMT='(" ...CPU time/atom/force call=",ES12.3)') (t2-t1)/100./REAL(AL%NAtoms)
      
      IsError=.FALSE.
      
      maxerr=0._dp
      DO i=1,NAtoms
         DO j=1,3 !dimension
            IF (.NOT. AL%AtomIsMoving(3*i+j-3)) CYCLE
            
            AL%AtomCoord(3*i+j-3)=AL%AtomCoord(3*i+j-3)+dr
            CALL AddVerletList(AL,.FALSE.,NAtoms=NAtoms) !dont reinitialize
            CALL GetEnergyVL(AL,NAtoms) !only energy
            ep=AL%PotentialEnergy
            
            AL%AtomCoord(3*i+j-3)=AL%AtomCoord(3*i+j-3)-2._dp*dr
            CALL AddVerletList(AL,.FALSE.,NAtoms=NAtoms) !dont reinitialize
            CALL GetEnergyVL(AL,NAtoms) !only energy
            en=AL%PotentialEnergy
            
            df=ABS((en-ep)/2._dp/dr-AL%AtomForce(3*i+j-3))
            maxerr=MAX(maxerr,ABS(100*df/((en-ep)/2._dp/dr)))
            AL%AtomCoord(3*i+j-3)=AL%AtomCoord(3*i+j-3)+dr
            
            IF ((en-ep)/2._dp/dr /=0._dp) THEN
               IF (ABS(100._dp*df/((en-ep)/2._dp/dr))>1._dp) THEN !more than 0.1% error
               !IF (df>1.e-6_dp) THEN !more than 0.1% error
                  IsError=.TRUE.
                  WRITE(*,*) "$Err>> Forces computed incorrectly using Verlet list"
                  WRITE(*,*) " >> Offending atom:",i
                  WRITE(*,*) " >> Offending dimension:",j
                  WRITE(*,*) " >> % error:",100*df/((en-ep)/2._dp/dr)
                  WRITE(*,*) " >> Approximate force:",(en-ep)/2._dp/dr
                  WRITE(*,*) " >> Calculated force:",AL%AtomForce(3*i+j-3)
               END IF
            END IF
         END DO
      END DO
      WRITE(UNIT=6,FMT='("... maximum error in computing forces:",ES10.3," %")') maxerr
      
      IF (IsError) THEN
         WRITE(*,*) "...forces are incorrect."
      ELSE
         WRITE(*,*) "...forces are correct."
      END IF
      
   END SUBROUTINE ChkForces
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

END MODULE PotentialPackage
