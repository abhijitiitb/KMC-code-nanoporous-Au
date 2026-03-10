MODULE db_manipulate
   USE VARIABLE_TYPE
   USE nrprocedures
   USE randomgen
   IMPLICIT NONE
   
   INTERFACE MakeSize
      MODULE PROCEDURE  MakeSize1Dsp, &
         MakeSize1Ddp,MakeSize1Dint,MakeSize1Dlog,MakeSize1Dchar,&
         MakeSize2Ddp,MakeSize2Dint,MakeSize2Dlog,MakeSize2Dchar,&
         MakeSizeAL,MakeSizeLL,MakeSizeVL
   END INTERFACE MakeSize
   INTERFACE Copy
      MODULE PROCEDURE CopyAL,CopyALFull,CopyPotential !,CopyALSection
   END INTERFACE Copy
   INTERFACE Delete
      MODULE PROCEDURE DeleteAL,DeleteVL,DeleteLL,DeleteChos,DeleteMD,DeletePotential
   END INTERFACE Delete
   INTERFACE Merge
      MODULE PROCEDURE MergeAL
   END INTERFACE Merge
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddALToMD(AL,MDSiml,dt,integrator,coeff, &
      temperature,TemperatureRamp,warnall,IsPrint)
      !add AL to MD simulation
      !  used when MDSiml is empty
      !  and AL contains all the atomic species, position, 
      !  velocity information   
      
      !Last checked: July 03, 2009
      !Last modified: 
      !Last modified by: AC
      
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(MDContainer), POINTER :: MDSiml
      
      REAL(dp), OPTIONAL :: dt,coeff,temperature,TemperatureRamp
      REAL(dp) :: dt1,coeff1,temperature1,TemperatureRamp1
      CHARACTER(len=10), OPTIONAL :: integrator
      CHARACTER(len=10) :: integrator1
      LOGICAL, OPTIONAL :: warnall,IsPrint
      LOGICAL :: warnall1,IsPrint1
      
      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(*,*) "$Err>> AtomList is not associated in AddALToMD"
         STOP
      END IF
      
      IF (PRESENT(IsPrint)) THEN
         IsPrint1=IsPrint
      ELSE
         IsPrint1=.TRUE.
      END IF
      
      IF (PRESENT(warnall)) THEN
         warnall1=warnall
      ELSE
         warnall1=.FALSE.
      END IF
      
      IF (PRESENT(integrator)) THEN
         integrator1=integrator
      ELSE
         integrator1="verlet"
      END IF
      
      IF (PRESENT(dt)) THEN
         dt1=dt
      ELSE
         SELECT CASE (integrator1(1:4))
         CASE("verl"); dt1=1.e-15_dp
         CASE("velo"); dt1=1.e-15_dp
         CASE("leap"); dt1=1.e-15_dp
         CASE("lang"); dt1=4.e-15_dp
         END SELECT
      END IF
      
      IF (PRESENT(coeff)) THEN
         coeff1=coeff
      ELSE
         IF (integrator1(1:4)=="lang") THEN
            coeff1=5.e12_dp
         ELSE
            coeff1=0._dp
         END IF
      END IF
      
      IF (PRESENT(temperature)) THEN !mainly important for Langevin dynamics
         temperature1=temperature
      ELSE
         IF (integrator1(1:4)=="lang") THEN
            WRITE(*,*) "$Err>> Provide temperature in AddALToMD for Langevin dynamics"
         ELSE
            temperature1=0._dp
         END IF
      END IF
      
      IF (PRESENT(TemperatureRamp)) THEN
         TemperatureRamp1=TemperatureRamp
      ELSE
         TemperatureRamp1=0._dp
      END IF
      
      IF (ASSOCIATED(MDSiml)) THEN
         IF (warnall1) THEN
            WRITE(*,*) "Warning>> Replacing AtomList in MDSiml in AddALToMD"
         END IF
         NULLIFY(MDSiml%AL)
      ELSE
         ALLOCATE(MDSiml)
      END IF

      MDSiml%AL=>AL
      MDSiml%AL%Time=0._dp
      MDSiml%dt=dt1
      MDSiml%Temperature=temperature1
      MDSiml%LangevinCoeff=coeff1
      MDSiml%AL=>AL
      MDSiml%Integrator=integrator1
      MDSiml%TemperatureRamp=TemperatureRamp1
      
      
      !MDSiml%History=>MDSiml%AL
      !MDSiml%HistorySizeMax=100
      !MDSiml%HistoryFreq=30
      !MDSiml%HistoryOn=.FALSE.

      IF (.NOT. MDSiml%HistoryOn) THEN
         MDSiml%HistorySizeMax=0
         !MDSiml%HistoryFreq=0
      END IF
      MDSiml%HistorySize=0
      
      IF (IsPrint1) THEN
         WRITE(UNIT=*,FMT='("DB>> Adding AL to MD container")') 
         WRITE(UNIT=*,FMT='(" ... MD time:",ES10.3)') MDSiml%AL%Time
         WRITE(UNIT=*,FMT='(" ... MD time increment:",ES12.3)') MDSiml%dt
         WRITE(*,*) "... MD integrator: ",MDSiml%Integrator
         IF (MDSiml%Integrator(1:4)=="lang") THEN
            WRITE(UNIT=*,FMT='(" ... MD temperature:",F10.3)') MDSiml%Temperature
            WRITE(UNIT=*,FMT='(" ... Langevin coeff:",ES14.3," /s")') MDSiml%LangevinCoeff
         END IF
      END IF
   END SUBROUTINE AddALToMD
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddALToMDHistory(MDSiml,iprint)
      IMPLICIT NONE
      !to add a snapshot to the history
      !AL and MDHistory are not linked 
      !AL contains the current state of the system
      !MDHistory contains the list of past ALs
      ! with the oldest at the beginning of the list 
      ! and the latest at the end of the list
      TYPE(MDContainer), POINTER :: MDSiml
      TYPE(SystemContainer), POINTER :: AL,AL1
      INTEGER, OPTIONAL :: iprint
      INTEGER :: iprint1,NAtoms
      
      IF (PRESENT(iprint)) THEN
         iprint1=iprint
      ELSE
         iprint1=0
      END IF
      
      IF (iprint1>0) WRITE(UNIT=6,FMT='("Storing MD snapshot in history @ time ",ES12.3, &
        " ps ...")',ADVANCE="NO") MDSiml%AL%Time*1.e12_dp
      
      IF (ASSOCIATED(MDSiml%History) .AND. MDSiml%HistorySize==0) THEN
         WRITE(*,*) "History is associated"
         STOP
      END IF
      
      IF (.NOT. ASSOCIATED(MDSiml%AL)) THEN
         WRITE(6,*) "$Err>> Atomlist to be added to MD history is not associated"
         STOP
      END IF
      
      !Add a new image
      NULLIFY(AL1)
      NAtoms=MDSiml%AL%NAtoms
      CALL MakeSize(AL1,NAtoms)
      CALL Copy(MDSiml%AL,AL1) !copy current image to AL
      
      !Check if we need to delete an image from MDHistory
      NULLIFY(AL)
      IF (MDSiml%HistorySize==MDSiml%HistorySizeMax) THEN
         AL=>MDSiml%History
         IF (ASSOCIATED(AL%NextNeigh)) THEN
            MDSiml%History=>AL%NextNeigh
            NULLIFY(MDSiml%History%PrevNeigh)
            NULLIFY(AL%NextNeigh)
            CALL Delete(AL)
            MDSiml%HistorySize=MDSiml%HistorySize-1
         ELSE
            IF (MDSiml%HistorySizeMax<=1) THEN
               WRITE(UNIT=6,FMT='("$Err>> Maximum number of snapshots stored" &
                  " in MDHistory is too small:",I3)') MDSiml%HistorySizeMax
               STOP
            END IF
            WRITE(6,*) "$Err>> Unable to find the history of MD trajectory"
            STOP
         END IF
      END IF
      
      !Go to the end of the history
      AL=>MDSiml%History
      IF (ASSOCIATED(AL)) THEN
         DO WHILE (ASSOCIATED(AL%NextNeigh))
            AL=>AL%NextNeigh
         END DO
      END IF
      
      !Add the latest image
      IF (ASSOCIATED(AL)) THEN
         AL%NextNeigh=>AL1
         AL1%PrevNeigh=>AL
      ELSE
         MDSiml%History=>AL1
         NULLIFY(AL1%PrevNeigh)
      END IF
      
      NULLIFY(AL1%NextNeigh)
      NULLIFY(AL1)
      
      !ALLOCATE(MDSiml%History%PrevNeigh) !old version
      !MDSiml%AL%PrevNeigh%NextNeigh=>MDSiml%AL
      !CALL Copy(MDSiml%AL,MDSiml%AL%PrevNeigh)
      !MDSiml%AL=>MDSiml%AL%NextNeigh
      
      IF (MDSiml%HistorySize==MDSiml%HistorySizeMax+1) THEN
         AL=>MDSiml%History
         MDSiml%History=>MDSiml%History%NextNeigh
         NULLIFY(AL%NextNeigh)
         NULLIFY(MDSiml%History%PrevNeigh)
         CALL Delete(AL)
      ELSE
         MDSiml%HistorySize=MDSiml%HistorySize+1
      END IF
      
      IF (iprint1>0) WRITE(UNIT=6,FMT='("[DONE]")')
      
   END SUBROUTINE AddALToMDHistory
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddChOSNEB(s1,s2,chos,interpmode,nimages) 
      !total images = n+1
      IMPLICIT NONE
      !see imode of NEB for description
      TYPE(SystemContainer), POINTER :: s1,s2 !start state
      TYPE(SystemContainer), POINTER :: ALPrev,ALCurr
      TYPE(ChOSContainer), POINTER :: chos,c
      REAL(dp), DIMENSION(3*s1%NAtoms) :: displacement
      REAL(dp), DIMENSION(s1%NAtoms) :: spacing
      INTEGER :: NAtoms,i,interpmode,iatom
      INTEGER :: nimages
      REAL(dp), DIMENSION(:), POINTER :: dr,r
      
      IF (ASSOCIATED(chos)) THEN
         WRITE(*,*) "Err>> chos should have been empty"
         STOP
      END IF
      
      IF (.NOT. ASSOCIATED(s1)) THEN
         WRITE(6,*) "$Err>> NEB initial state is not initialized"
         STOP
      END IF
      
      IF (.NOT. ASSOCIATED(s2)) THEN
         WRITE(6,*) "$Err>> NEB final state is not initialized"
         STOP
      END IF
      
      IF (s1%NAtoms/=s2%NAtoms) THEN
         WRITE(6,*) "$Err>> Number of atoms in initial and final states in NEB do not match"
         STOP
      END IF
      
      NAtoms=s1%NAtoms
      
      CALL MakeSize(dr,3*NAtoms)
      CALL MakeSize(r,3*NAtoms)
      
      dr=0._dp !diplacement from start location
      CALL AddChosImage(s1,chos,dr,AddChOSvariablesInAL=.TRUE.) !first image
      c=>chos
      
      IF (interpmode==1) THEN !linear interpolation
         IF (nimages<3) THEN
            WRITE(*,*) "$Err>> Too few images for NEB"
            STOP
         END IF
         !dr=(s2%AtomCoord(1:3*NAtoms)-s1%AtomCoord(1:3*NAtoms))/REAL(nimages-1)
         dr=AMDDbPBCdistance2(s2,s1)/REAL(nimages-1)
         r=0._dp
         DO i=2,nimages
            r=r+dr
            CALL AddChosImage(s1,c,r,AddChOSvariablesInAL=.TRUE.)
         END DO
         chos%nimages=nimages
      ELSEIF (interpmode==2) THEN 
         !growing string method- uses 4 images to begin with
         dr=0.01_dp*(s2%AtomCoord(1:3*NAtoms)-s1%AtomCoord(1:3*NAtoms))  !as if 100 images were present
         CALL AddChosImage(s1,c,dr,AddChOSvariablesInAL=.TRUE.)
         dr=99._dp*dr
         CALL AddChosImage(s1,c,dr,AddChOSvariablesInAL=.TRUE.)
         dr=(100._dp/99._dp)*dr
         CALL AddChosImage(s1,c,dr,AddChOSvariablesInAL=.TRUE.)
         chos%nimages=4
      ELSEIF (interpmode==3) THEN !mixed linear-nonlinear interpolation, atoms that overlap are shifted to a new position
         IF (nimages<3) THEN
            WRITE(*,*) "$Err>> Too few images for NEB"
            STOP
         END IF
         !dr=(s2%AtomCoord(1:3*NAtoms)-s1%AtomCoord(1:3*NAtoms))/REAL(nimages-1)
         dr=AMDDbPBCdistance2(s2,s1)/REAL(nimages-1)
         r=0._dp
         ALPrev=>chos%AL
         DO i=2,nimages
            r=r+dr
            !Step 1. Begin with linear interpolation
            CALL AddChosImage(s1,c,r,AddChOSvariablesInAL=.TRUE.,ALPrev=ALPrev,ALCurr=ALCurr) !ALCurr is the AL just added
            !Step 2. Get the movement from the initial position/state
            displacement=AMDDbPBCdistance2(ALCurr,s1) !displacement with respect to the previous 
            DO iatom=1,NAtoms
               spacing(iatom)=SQRT(displacement(3*iatom-2)*displacement(3*iatom-2)+ &
                  displacement(3*iatom-1)*displacement(3*iatom-1)+ &
                  displacement(3*iatom)*displacement(3*iatom) )
            END DO
            !Step 3. Correct the positions
            CALL NudgeImage(ALCurr,ALPrev,spacing,d=1._dp) !nudges ALCurr given ALPrev such that the spacing betwen atoms is more than d
         END DO
         chos%nimages=nimages
      END IF

      DEALLOCATE(dr)
      DEALLOCATE(r)
      
   END SUBROUTINE AddChOSNEB
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION AMDDbPBCdistance2(AL1,AL2)
   !returns PBC(AL1Coord-AL2Coord)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL1,AL2 
      REAL(dp), DIMENSION(3*AL1%NAtoms) :: AMDDbPBCdistance2
      REAL(dp), DIMENSION(3*MaxNAtoms) :: BoxSize
      INTEGER :: NAtoms
      
      IF (AL1%NAtoms/=AL2%NAtoms) THEN
         WRITE(6,*) "$Err>> Size of Atom lists do not match to get distance"
         STOP
      END IF
      IF (ANY(AL1%BoxSize/=AL2%BoxSize)) THEN
         WRITE(6,*) "$Err>> Size of Atom list boxes do not match to get distance"
         STOP
      END IF
      
      NAtoms=AL1%NAtoms
      BoxSize(1:3*NAtoms)=PACK(SPREAD(AL1%BoxSize(1:3)-AL1%BoxSize(4:6),2,NAtoms),.TRUE.)
      AMDDbPBCdistance2(1:3*NAtoms)=AL1%AtomCoord(1:3*NAtoms)-AL2%AtomCoord(1:3*NAtoms)
      AMDDbPBCdistance2(1:3*NAtoms)=AMDDbPBCdistance2(1:3*NAtoms)-BoxSize(1:3*NAtoms)* &
         NINT(AMDDbPBCdistance2(1:3*NAtoms)/BoxSize(1:3*NAtoms))
   END FUNCTION AMDDbPBCdistance2
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddChosImage(AL,chos,dr,AddChOSvariablesInAL,ALPrev,ALCurr)
      !add image that looks like AL but displaced by dr
      !image is always added to the end of linked list
      TYPE(ChOSContainer), POINTER :: chos
      TYPE(SystemContainer), POINTER :: AL,AL1,AL2
      TYPE(SystemContainer), POINTER, OPTIONAL :: ALPrev,ALCurr
      REAL(dp), DIMENSION(:), POINTER :: dr
      LOGICAL, OPTIONAL :: AddChOSvariablesInAL
      LOGICAL :: AddChOSvariablesInAL1
      INTEGER :: NAtoms
      
      IF (.NOT. ASSOCIATED(chos)) THEN
         ALLOCATE(chos)
         chos%nimages=0
      END IF
      
      IF (PRESENT(AddChOSvariablesInAL)) THEN
         AddChOSvariablesInAL1=AddChOSvariablesInAL
      ELSE
         AddChOSvariablesInAL1=.FALSE.
      END IF
      
      NAtoms=AL%NAtoms
      NULLIFY(AL1)
      CALL MakeSize(AL1,NAtoms,AddChOSvariablesInAL=AddChOSvariablesInAL1)
      CALL Copy(AL,AL1)
      IF (SIZE(AL1%AtomCoord)==SIZE(dr)) THEN
         WHERE(AL1%AtomIsMoving(1:3*NAtoms)) &
            AL1%AtomCoord(1:3*NAtoms)=AL1%AtomCoord(1:3*NAtoms)+dr(1:3*NAtoms)
      ELSE
         WRITE(6,*) "$Err>> Size of arrays do not match in AddChOSImage"
         STOP
      END IF
      NULLIFY(AL1%PrevNeigh)
      NULLIFY(AL1%NextNeigh)
      
      IF (ASSOCIATED(chos%AL)) THEN
         AL2=>chos%AL
         !reach the end of the list
         DO WHILE (ASSOCIATED(AL2%NextNeigh))
            AL2=>AL2%NextNeigh
         END DO
         AL2%NextNeigh=>AL1
         AL1%PrevNeigh=>AL2
         IF (PRESENT(ALPrev)) ALPrev=>AL2
         IF (PRESENT(ALCurr)) ALCurr=>AL1
         NULLIFY(AL2)
         NULLIFY(AL1)
      ELSE
         chos%AL=>AL1
         NULLIFY(AL1)
      END IF
      
      chos%nimages=chos%nimages+1
      WRITE(UNIT=6,FMT='(" ... number of images:",I3)') chos%nimages
   END SUBROUTINE AddChosImage
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION AMDDbPBCSpacing(AL,atom1,atom2)
   !identical to the function AMDDbPBCSpacing defined in amd_utilities.f90
   !difference vector from atom2
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: atom1,atom2
      REAL(dp), DIMENSION(3) :: r1,r2,vec1,vec2,vec3,AMDDbPBCSpacing,BoxSize
      TYPE(SystemContainer), POINTER :: AL
      
      IF (atom1==atom2) THEN
         AMDDbPBCSpacing=0._dp
         RETURN
      END IF
      
      r1=AL%AtomCoord(3*atom1-2:3*atom1)
      r2=AL%AtomCoord(3*atom2-2:3*atom2)
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      
      AMDDbPBCSpacing=r1-r2 !vector from atom2
      AMDDbPBCSpacing=AMDDbPBCSpacing-BoxSize*NINT(AMDDbPBCSpacing/BoxSize)
   END FUNCTION AMDDbPBCSpacing
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE NudgeImage(ALCurr,ALPrev,spacing,d)
   !nudges ALCurr given ALPrev such that the spacing betwen atoms in ALCurr is more than d
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: ALCurr,ALPrev
      REAl(dp), DIMENSION(ALCurr%NAtoms) :: spacing !contains the displacement from the original state
      REAL(dp) :: d,dcurr,dprev,sgn,vec_curr(3),vec_prev(3)
      INTEGER :: iatom,jatom,NAtoms,iatom1,jatom1
      
      !Step 1. Find which the movement of the atoms from previous image
      NAtoms=ALCurr%NAtoms
      
      !Step 2. Find whether two atoms are within distance d -- the coordinates are x1 and x2
      DO iatom=1,NAtoms
         IF (spacing(iatom)>0.3_dp) THEN !iatom has moved a lot -- hopefully this is just a handful of atoms
            DO jatom=1,NAtoms !compare with all atoms and check if it is overlapping
               IF (iatom==jatom) CYCLE
               vec_curr=AMDDbPBCSpacing(ALCurr,iatom,jatom) !pbc vector from jatom to iatom
               dcurr=SQRT(vec_curr(1)*vec_curr(1)+vec_curr(2)*vec_curr(2)+vec_curr(3)*vec_curr(3)) !current spacing
               IF (dcurr<d) THEN !we need to move the atoms apart
                  !iatom1 is the atom which has moved the most
                  iatom1=iatom
                  jatom1=jatom
                  sgn=1._dp
                  IF (spacing(iatom)<spacing(jatom)) THEN
                     iatom1=jatom
                     jatom1=iatom
                     sgn=-1._dp !go in the opposite direction
                  END IF
                  vec_prev=AMDDbPBCSpacing(ALPrev,iatom1,jatom1) !vector from jatom to iatom
                  dprev=SQRT(vec_prev(1)*vec_prev(1)+vec_prev(2)*vec_prev(2)+vec_prev(3)*vec_prev(3)) !current spacing
                  !displace the atom which has moved the most
                  ALCurr%AtomCoord(3*iatom1-2:3*iatom1)=ALCurr%AtomCoord(3*jatom1-2:3*jatom1)+sgn*vec_curr*dprev/dcurr
                  !hopefully this will take care of the overlapping position and will not be undone
               END IF
            END DO
         END IF
      END DO
      
      !Step 3. Find the distance between these two atoms in previous image dp and in the current image dc
      
      !Step 4. Get new separation dp*(x1-x2)/dc
      
      !Step 5. Assign the atom which has moved the most displacement
   END SUBROUTINE NudgeImage
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CountImages(AL)
   !counts the number of images that are present in the chos formed by AL
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1
      INTEGER :: CountImages
      
      CountImages=0
      AL1=>AL
      DO WHILE (ASSOCIATED(AL1))
         CountImages=CountImages+1
         AL1=>AL1%NextNeigh
      END DO
   END FUNCTION CountImages
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddChOSDimer(s1,chos)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: s1 !start state
      !TYPE(SystemContainer), POINTER, OPTIONAL :: s2 !end state
      TYPE(ChOSContainer), POINTER :: chos
      INTEGER :: imode,NAtoms,i
      REAL(dp), DIMENSION(:), POINTER :: dr
      REAL(sp), DIMENSION(:), POINTER :: dr_rand
      
      IF (ASSOCIATED(chos)) THEN
         WRITE(*,*) "Err>> chos should have been empty"
         STOP
      END IF
      WRITE(6,*) "Adding dimer ..."
      
      NAtoms=s1%NAtoms
      CALL MakeSize(dr,3*NAtoms)
      CALL MakeSize(dr_rand,3*NAtoms)
      
      CALL ran1_v(dr_rand)
      dr=(REAL(dr_rand,dp)-0.5_dp)*0.01_dp !maximum distance for atom from original position is 0.005 Ang
      CALL AddChosImage(s1,chos,dr)
      dr=-1._dp*dr
      CALL AddChosImage(s1,chos,dr)
      
      DEALLOCATE(dr)
      DEALLOCATE(dr_rand)
   END SUBROUTINE AddChOSDimer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteChos(chos,entire1)
      IMPLICIT NONE
      TYPE(ChOSContainer), POINTER :: chos
      TYPE(SystemContainer), POINTER :: c
      LOGICAL, OPTIONAL :: entire1
      LOGICAL :: entire
      
      IF (.NOT. ASSOCIATED(chos)) RETURN
      
      IF (PRESENT(entire1)) THEN
         entire=entire1
      ELSE
         entire=.FALSE. !remove only current image
      END IF
      
      IF (entire) THEN
         CALL Delete(chos%AL)
         DEALLOCATE(chos)
      ELSE
         c=>chos%AL
         IF (ASSOCIATED(c%PrevNeigh)) THEN !this might not work properly -- are the pointers for AL%PrevNeigh and AL%NextNeigh working properly in DeleteAL(AL,.FALSE.)??
            IF (ASSOCIATED(c%NextNeigh)) THEN
               c%PrevNeigh%NextNeigh=>c%NextNeigh
               c%NextNeigh%PrevNeigh=>c%PrevNeigh
            ELSE
               NULLIFY(c%PrevNeigh%NextNeigh)
            END IF
         ELSE
            IF (ASSOCIATED(c%NextNeigh)) THEN
               NULLIFY(c%NextNeigh%PrevNeigh)
            END IF
         END IF
         CALL Delete(c,.FALSE.)
      END IF
   END SUBROUTINE DeleteChos
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteAL(AL,entire)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1
      LOGICAL, OPTIONAL :: entire
      LOGICAL :: entire1
      
      IF (.NOT. ASSOCIATED(AL)) RETURN
      !IF (ASSOCIATED(AL%SpeciesDirectory)) DEALLOCATE(AL%SpeciesDirectory)
      
      IF (PRESENT(entire)) THEN
         entire1=entire
      ELSE
         entire1=.TRUE.
      END IF
      
      AL1=>AL
      IF (entire1) THEN
         DO WHILE (ASSOCIATED(AL1))
            AL=>AL%NextNeigh
            NULLIFY(AL1%PrevNeigh)
            NULLIFY(AL1%NextNeigh)
            CALL DeleteALImage(AL1)
            AL1=>AL
         END DO
      ELSE
         AL=>AL%NextNeigh
         CALL DeleteALImage(AL1)
      END IF
   END SUBROUTINE DeleteAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteALImage(AL)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      
      IF (ASSOCIATED(AL%AtomCharge)) DEALLOCATE(AL%AtomCharge)
      IF (ASSOCIATED(AL%AtomSpecies)) DEALLOCATE(AL%AtomSpecies)
      IF (ASSOCIATED(AL%AtomCoord)) DEALLOCATE(AL%AtomCoord)
      IF (ASSOCIATED(AL%AtomVelocity)) DEALLOCATE(AL%AtomVelocity)
      IF (ASSOCIATED(AL%AtomDriftVelocity)) DEALLOCATE(AL%AtomDriftVelocity)
      IF (ASSOCIATED(AL%AtomMass)) DEALLOCATE(AL%AtomMass)
      IF (ASSOCIATED(AL%AtomForce)) DEALLOCATE(AL%AtomForce)
      IF (ASSOCIATED(AL%AtomShiftCoord)) DEALLOCATE(AL%AtomShiftCoord)
      IF (ASSOCIATED(AL%AtomIsMoving)) DEALLOCATE(AL%AtomIsMoving)
      IF (ASSOCIATED(AL%ChOSTangent)) DEALLOCATE(AL%ChOSTangent)
      IF (ASSOCIATED(AL%ChOSForce)) DEALLOCATE(AL%ChOSForce)
      
      !IF (ASSOCIATED(AL%AtomEnergyHessian)) DEALLOCATE(AL%AtomEnergyHessian)
      
      IF (ASSOCIATED(AL%VL)) CALL Delete(AL%VL)
      IF (ASSOCIATED(AL%LL)) CALL Delete(AL%LL)
      IF (ASSOCIATED(AL%Potential)) CALL Delete(AL%Potential)
      DEALLOCATE(AL)
         
   END SUBROUTINE DeleteALImage
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteVL(VL)
      IMPLICIT NONE
      TYPE(VerletListContainer), POINTER :: VL
      
      IF (.NOT. ASSOCIATED(VL)) RETURN
      
      IF (ASSOCIATED(VL%ListRange)) DEALLOCATE(VL%ListRange)
      IF (ASSOCIATED(VL%ListDomainAtom)) DEALLOCATE(VL%ListDomainAtom)
      IF (ASSOCIATED(VL%List)) DEALLOCATE(VL%List)
      IF (ASSOCIATED(VL%drmag)) DEALLOCATE(VL%drmag)
      IF (ASSOCIATED(VL%OrigCoord)) DEALLOCATE(VL%OrigCoord) !earlier OrigCoord was of fixed size, pre-allocated
      IF (ASSOCIATED(VL%dr)) DEALLOCATE(VL%dr)
      
      DEALLOCATE(VL)
   END SUBROUTINE DeleteVL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteLL(LL)
      IMPLICIT NONE
      TYPE(LinkedListContainer), POINTER :: LL
      
      IF (.NOT. ASSOCIATED(LL)) RETURN
      IF (ASSOCIATED(LL%ListRange)) DEALLOCATE(LL%ListRange)
      IF (ASSOCIATED(LL%List)) DEALLOCATE(LL%List)
      IF (ASSOCIATED(LL%AtomCellIndx)) DEALLOCATE(LL%AtomCellIndx)
      DEALLOCATE(LL)
   END SUBROUTINE DeleteLL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteMD(MDSiml,AL)
      !removes the MD container & return the Atom list
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      TYPE(SystemContainer), POINTER, OPTIONAL :: AL
      TYPE(SystemContainer), POINTER :: AL1
      INTEGER :: i

      IF (.NOT. ASSOCIATED(MDSiml)) RETURN
      
      IF (PRESENT(AL)) THEN
         AL=>MDSiml%AL
         NULLIFY(MDSiml%AL)
      ELSE
         CALL Delete(MDSiml%AL)
      END IF

      CALL DeleteMDHistory(MDSiml) !delete the history
      
      DEALLOCATE(MDSiml)
   END SUBROUTINE DeleteMD
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteMDHistory(MDSiml)
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: i
      
      NULLIFY(AL)
      DO i=1,MDSiml%HistorySize
         AL=>MDSiml%History
         IF (ASSOCIATED(MDSiml%History%NextNeigh)) THEN
            MDSiml%History=>MDSiml%History%NextNeigh
            NULLIFY(MDSiml%History%PrevNeigh)
         ELSE
            NULLIFY(MDSiml%History)
         END IF
         NULLIFY(AL%NextNeigh)
         CALL Delete(AL)
      END DO
      MDSiml%HistorySize=0
   END SUBROUTINE DeleteMDHistory
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeletePotential(Potential)
      IMPLICIT NONE
      TYPE(InteracPotential), POINTER :: Potential
      
      IF (.NOT. ASSOCIATED(Potential)) RETURN
      
      CALL DeleteLJ(Potential%LJ)
      CALL DeleteEAM(Potential%EAM)
      CALL DeleteCoulomb(Potential%Coulomb)
      CALL DeleteSW(Potential%SW)
      DEALLOCATE(Potential)
   END SUBROUTINE DeletePotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteLJ(LJ)
      IMPLICIT NONE
      TYPE(LJPotential), POINTER :: LJ
      INTEGER :: nx,ny,ix,iy
      
      IF (.NOT. ASSOCIATED(LJ)) RETURN
      DEALLOCATE(LJ)
   END SUBROUTINE DeleteLJ
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteEAM(EAM)
      IMPLICIT NONE
      TYPE(EAMPotential), POINTER :: EAM
      INTEGER :: nx,ny,ix,iy
      
      IF (.NOT. ASSOCIATED(EAM)) RETURN
      DEALLOCATE(EAM)
   END SUBROUTINE DeleteEAM
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteCoulomb(Coulomb)
      IMPLICIT NONE
      TYPE(CoulombPotential), POINTER :: Coulomb
      
      IF (.NOT. ASSOCIATED(Coulomb)) RETURN
      
      CALL DeleteCoulombEwald(Coulomb%Ewald)
      CALL DeleteCoulombWolf(Coulomb%Wolf)
      CALL DeleteCoulombPPPM(Coulomb%PPPM)
      CALL DeleteCoulombFMM(Coulomb%FMM)
   END SUBROUTINE DeleteCoulomb
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteCoulombEwald(Ewald)
      IMPLICIT NONE
      TYPE(EwaldParameters), POINTER :: Ewald
      
      IF (.NOT. ASSOCIATED(Ewald)) RETURN
      
      !IF (ASSOCIATED(Ewald%EIKX)) DEALLOCATE(Ewald%EIKX)
      !IF (ASSOCIATED(Ewald%EIKY)) DEALLOCATE(Ewald%EIKY)
      !IF (ASSOCIATED(Ewald%EIKZ)) DEALLOCATE(Ewald%EIKZ)
      IF (ASSOCIATED(Ewald%KVEC)) DEALLOCATE(Ewald%KVEC)
      !IF (ASSOCIATED(Ewald%RVEC)) DEALLOCATE(Ewald%RVEC)
      IF (ASSOCIATED(Ewald%KVector)) DEALLOCATE(Ewald%KVector)
      
      DEALLOCATE(Ewald)
   END SUBROUTINE DeleteCoulombEwald
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteCoulombWolf(Wolf)
      IMPLICIT NONE
      TYPE(WolfParameters), POINTER :: Wolf
      
      IF (.NOT. ASSOCIATED(Wolf)) RETURN
   END SUBROUTINE DeleteCoulombWolf
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteCoulombPPPM(PPPM)
      IMPLICIT NONE
      TYPE(PPPMParameters), POINTER :: PPPM
      
      IF (.NOT. ASSOCIATED(PPPM)) RETURN
   END SUBROUTINE DeleteCoulombPPPM
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteCoulombFMM(FMM)
      IMPLICIT NONE
      TYPE(FMMParameters), POINTER :: FMM
      
      IF (.NOT. ASSOCIATED(FMM)) RETURN
   END SUBROUTINE DeleteCoulombFMM
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteSW(SW)
      IMPLICIT NONE
      TYPE(StillingerWeberPotential), POINTER :: SW
      
      IF (.NOT. ASSOCIATED(SW)) RETURN
      DEALLOCATE(SW)
   END SUBROUTINE DeleteSW
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteListOfScalars(arr)
      IMPLICIT NONE
      TYPE(ListOfScalars), POINTER :: arr,t
      
      DO WHILE (ASSOCIATED(arr))
         t=>arr
         arr=>arr%NextNeigh
         IF (ASSOCIATED(arr)) NULLIFY(arr%PrevNeigh)
         NULLIFY(t%NextNeigh)
         DEALLOCATE(t)
      END DO
   END SUBROUTINE DeleteListOfScalars
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MergeAL(AL0,AL1,AL2)
      !stronger version of CopyAL, also similar to Add in amd_crystal.f90
      !merge AL0 and AL1 into AL2
      !performs -
      !  merging
      !  box resizing
      !  species renumbering
      !  overlap check

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL1,AL2,AL0
      INTEGER :: i,AtomicNo,indx,p1,p2,r1,r2
      INTEGER :: NAtoms !,NSpecies
      LOGICAL :: RenumberingRequired
      
      !Atom count merge
      NAtoms=AL0%NAtoms+AL1%NAtoms
      CALL MakeSize(AL2,NAtoms)

      !Box resizing
      AL2%BoxSize(1:3)=MAX(AL0%BoxSize(1:3),AL1%BoxSize(1:3))
      AL2%BoxSize(4:6)=MIN(AL0%BoxSize(4:6),AL1%BoxSize(4:6))

      CALL MergeALSlow()

   CONTAINS
   !oooooooooooooooooooooooooooooo
   SUBROUTINE MergeALSlow()

      !Species renumbering-add atoms one by one
      !--------------------
      !incase AL0 has fewer species than AL1 then we need to add them
      !DO i=1,AL0%NSpecies
      !   !AtomicNo=AL0%SpeciesDirectory(i)
      !   AtomicNo=SpeciesDirectory_global(i)
      !   indx= GetSpecies(AL2,AtomicNo,ForceAdd1=.TRUE.) !forcibly adding the species into SpeciesDirectory
      !END DO
      !DO i=1,AL1%NSpecies
      !   !AtomicNo=AL1%SpeciesDirectory(i)
      !   AtomicNo=SpeciesDirectory(i)
      !   indx= GetSpecies(AL2,AtomicNo,ForceAdd1=.TRUE.)
      !END DO
      !NSpeciesType=AL2%NSpeciesType

      !merge
      !IF (AL1%BoundaryCondX/=AL%BoundaryCondX .OR.&
      !     AL1%BoundaryCondY/=AL%BoundaryCondY .OR. &
      !     AL1%BoundaryCondZ/=AL%BoundaryCondZ ) THEN
      !   WRITE(*,*) "$Err>> Boundary conditions do not match"
      !   STOP
      !END IF

      r1=1
      r2=AL0%NAtoms
      p1=1
      p2=3*AL0%NAtoms
      AL2%AtomCharge(r1:r2)=AL0%AtomCharge
      AL2%AtomSpecies(r1:r2)=AL0%AtomSpecies

      AL2%AtomCoord(p1:p2)=AL0%AtomCoord
      AL2%AtomVelocity(p1:p2)=AL0%AtomVelocity
      AL2%AtomMass(p1:p2)=AL0%AtomMass
      AL2%AtomForce(p1:p2)=AL0%AtomForce
      AL2%AtomShiftCoord(p1:p2)=AL0%AtomShiftCoord
      AL2%AtomIsMoving(p1:p2)=AL0%AtomIsMoving

      !now add AL1
      r1=r2+1
      r2=r2+AL1%NAtoms
      p1=p2+1
      p2=p2+3*AL1%NAtoms
      AL2%AtomCharge(r1:r2)=AL1%AtomCharge
      AL2%AtomSpecies(r1:r2)=AL1%AtomSpecies
      !We need to update the atom species ... (older version)
      !in new version atom species is fixed from beginning.
      !DO i=1,AL1%NAtoms
      !   indx=AL1%AtomSpecies(i)
      !   !AtomicNo=AL1%SpeciesDirectory(indx)
      !   AtomicNo=SpeciesDirectory(indx)
      !   !now check which element of AL2%SpeciesDirectory is AtomicNo
      !   indx=GetSpecies(AL2,AtomicNo)
      !   AL2%AtomSpecies(AL0%NAtoms+i)=indx
      !END DO
      AL2%AtomCoord(p1:p2)=AL1%AtomCoord
      AL2%AtomVelocity(p1:p2)=AL1%AtomVelocity
      AL2%AtomDriftVelocity(p1:p2)=AL1%AtomDriftVelocity
      AL2%AtomMass(p1:p2)=AL1%AtomMass
      AL2%AtomForce(p1:p2)=AL1%AtomForce
      AL2%AtomShiftCoord(p1:p2)=AL1%AtomShiftCoord
      AL2%AtomIsMoving(p1:p2)=AL1%AtomIsMoving

   END SUBROUTINE MergeALSlow
   !oooooooooooooooooooooooooooooo
   END SUBROUTINE MergeAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CopyAL(AL,AL1,icopy) !copy AL into AL1
      !icopy=0: copy all
      !icopy=2: copy coordinates
      !icopy=3: copy velocity and force
      !icopy=5: copy mass, species and charge
      !icopy=7: copy neighborlist -- not in use
      !icopy=11: copy shift coord
      !icopy=13: copy is moving
      !icopy=17: copy potential
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1
      INTEGER, OPTIONAL :: icopy
      INTEGER :: icopy1
      
      icopy1=0
      IF (PRESENT(icopy)) icopy1=icopy
      
      CALL CopyALFull(AL,AL1,1,AL%NAtoms,icopy1)
   END SUBROUTINE CopyAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CopyALFull(AL,AL1,r1,r2,icopy) 
      !copy AL into AL1 (empty list)
      !should not be used for merging two lists 
      !icopy=0: copy all
      !icopy=2: copy coordinates
      !icopy=3: copy velocity, drift velocity and force
      !icopy=5: copy mass, species and charge
      !icopy=7: copy neighborlist -- not in use
      !icopy=11: copy shift coord
      !icopy=13: copy is moving
      !icopy=17: copy potential
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1
      INTEGER :: i,r1,r2 !range in AL copied
      INTEGER :: NAtoms,p1,p2,icopy1
      INTEGER, OPTIONAL :: icopy

      icopy1=0
      IF (PRESENT(icopy)) icopy1=icopy
      
      !Rule 1: AL1 should be already allocated
      IF (.NOT. ASSOCIATED(AL1)) THEN
         WRITE(*,*) "$Err>> CopyAL: Atom list is not associated"
         STOP
      END IF
      
      !Rule 2: AL should have atoms in it
      NAtoms=AL%NAtoms
      IF (NAtoms==0) THEN
         WRITE(*,*) "Err>> CopyAL: attempt to copy empty AtomList"
         RETURN
      END IF
      
      !Rule 3: AL1 should exactly fit the range of atoms specified
      IF (AL1%NAtoms<r2-r1+1) THEN
         WRITE(*,*) "Err>> CopyAL: range of AL1 copied exceeded"
         STOP
      END IF
      
      !range of direct copy (xyz of one atom are together in array)
      p1=3*r1-2
      p2=3*r2
      
      AL1%Time=AL%Time
      
      IF (MOD(icopy1,5)==0) AL1%AtomCharge(r1:r2)=AL%AtomCharge(1:NAtoms)
      IF (MOD(icopy1,5)==0) AL1%AtomSpecies(r1:r2)=AL%AtomSpecies(1:NAtoms)

      IF (MOD(icopy1,2)==0) AL1%AtomCoord(p1:p2)=AL%AtomCoord(1:3*NAtoms)
      IF (MOD(icopy1,3)==0) AL1%AtomVelocity(p1:p2)=AL%AtomVelocity(1:3*NAtoms)
      IF (MOD(icopy1,3)==0) AL1%AtomDriftVelocity(p1:p2)=AL%AtomDriftVelocity(1:3*NAtoms)
      IF (MOD(icopy1,5)==0) AL1%AtomMass(p1:p2)=AL%AtomMass(1:3*NAtoms)
      IF (MOD(icopy1,3)==0) AL1%AtomForce(p1:p2)=AL%AtomForce(1:3*NAtoms)
      IF (MOD(icopy1,11)==0) AL1%AtomShiftCoord(p1:p2)=AL%AtomShiftCoord(1:3*NAtoms)
      
      IF (MOD(icopy1,13)==0) AL1%AtomIsMoving(p1:p2)=AL%AtomIsMoving(1:3*NAtoms)

      AL1%BoundaryCondX=AL%BoundaryCondX
      AL1%BoundaryCondY=AL%BoundaryCondY
      AL1%BoundaryCondZ=AL%BoundaryCondZ
      
      AL1%PotentialEnergy=AL%PotentialEnergy
      AL1%KineticEnergy=AL%KineticEnergy
      
      AL1%BoxSize(1:3)=AL%BoxSize(1:3)
      AL1%BoxSize(4:6)=AL%BoxSize(4:6)

      IF (MOD(icopy1,13)==0) THEN
         AL1%NMove=0
         p1=3*r1-2
         p2=3*r1 !starting range
         DO i=r1,r2
            WHERE (AL1%AtomIsMoving(p1:p2)) AL1%NMove=AL1%NMove+1
            p1=p1+3
            p2=p2+3
         END DO
      ELSE
         AL1%NMove=AL%NMove
      END IF
      
      !species terms & potentials
      !for now its assumed that single species array for all ALs
      !potentials can also be directly copied in this way
      
      !IF (AL1%NSpeciesType==0) THEN
      !   AL1%NSpeciesType=AL%NSpeciesType
         !CALL MakeSize(AL1%SpeciesDirectory,AL1%NSpeciesType)
         !AL1%SpeciesDirectory=AL%SpeciesDirectory
         IF (MOD(icopy1,17)==0) CALL Copy(AL%Potential,AL1%Potential)
         !AL1%Potential=>AL%Potential
      !ELSE
      !   WRITE(*,*) "$Err>> CopyAL: AL1 seems to already contain atoms"
      !   STOP
      !END IF
      
   END SUBROUTINE CopyALFull
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CopyPotential(Potential,Potential1)
      IMPLICIT NONE
      TYPE(InteracPotential), POINTER :: Potential,Potential1
      TYPE(LJPotential), POINTER :: LJ,LJ1
      TYPE(EAMPotential), POINTER :: EAM,EAM1
      TYPE(CoulombPotential), POINTER :: Coulomb
      TYPE(EwaldParameters), POINTER :: Ewald,Ewald1
      TYPE(StillingerWeberPotential), POINTER :: SW,SW1
      INTEGER :: nx,ny
      
      IF (ASSOCIATED(Potential)) THEN
         
         ALLOCATE(Potential1)
         
         Potential1%UnInitialized=Potential%UnInitialized
         !Potential1%NSpecies=Potential%NSpecies
         
         Potential1%VLBuffer=Potential%VLBuffer
         Potential1%MaxPotentialCutoff=Potential%MaxPotentialCutoff
         
         IF (ASSOCIATED(Potential%LJ)) THEN
         
            ALLOCATE(Potential1%LJ)
            
            LJ=>Potential%LJ
            LJ1=>Potential1%LJ
            
            LJ1%TableLocationPP=LJ%TableLocationPP
            LJ1%SizePP=LJ%SizePP
            LJ1%RangePP=LJ%RangePP
            LJ1%IncrPPinv=LJ%IncrPPinv
            LJ1%PP=LJ%PP
               
         END IF
         
         IF (ASSOCIATED(Potential%EAM)) THEN
         
            ALLOCATE(Potential1%EAM)
            
            EAM=>Potential%EAM
            EAM1=>Potential1%EAM
            
            EAM1%TableLocationPP=EAM%TableLocationPP
            EAM1%TableLocationDT=EAM%TableLocationDT
            EAM1%TableLocationEE=EAM%TableLocationEE
            EAM1%SizePP=EAM%SizePP
            EAM1%SizeDT=EAM%SizeDT
            EAM1%SizeEE=EAM%SizeEE
            EAM1%RangePP=EAM%RangePP
            EAM1%RangeDT=EAM%RangeDT
            EAM1%RangeEE=EAM%RangeEE
            EAM1%IncrPPinv=EAM%IncrPPinv
            EAM1%IncrDTinv=EAM%IncrDTinv
            EAM1%IncrEEinv=EAM%IncrEEinv
            EAM1%PP=EAM%PP
            EAM1%DT=EAM%DT
            EAM1%EE=EAM%EE
         
         END IF
         
         IF (ASSOCIATED(Potential%Coulomb)) THEN
         
            Coulomb=>Potential%Coulomb
            ALLOCATE(Potential1%Coulomb)
            
            IF (ASSOCIATED(Coulomb%Ewald)) THEN
            
               ALLOCATE(Potential1%Coulomb%Ewald)
               Ewald=>Potential%Coulomb%Ewald
               Ewald1=>Potential1%Coulomb%Ewald
               
               Ewald1%NAtoms=Ewald%NAtoms
               Ewald1%ALPHA=Ewald%ALPHA
               Ewald1%ReCutOff=Ewald%ReCutOff
               Ewald1%TimeRatio=Ewald%TimeRatio
               Ewald1%Precision=Ewald%Precision
               Ewald1%ImCutOff=Ewald%ImCutOff
               !Ewald1%TOTr=Ewald%TOTr
               Ewald1%TOTk=Ewald%TOTk
            
               !IF (ASSOCIATED(Ewald%EIKX)) THEN
               !   nx=SIZE(Ewald%EIKX,1)
               !   ny=SIZE(Ewald%EIKX,2)
               !   ALLOCATE(Ewald1%EIKX(nx,ny))
               !   Ewald1%EIKX=Ewald%EIKX
               !END IF
               
               !IF (ASSOCIATED(Ewald%EIKY)) THEN
               !   nx=SIZE(Ewald%EIKY,1)
               !   ny=SIZE(Ewald%EIKY,2)
               !   ALLOCATE(Ewald1%EIKY(nx,ny))
               !   Ewald1%EIKY=Ewald%EIKY
               !END IF
               
               !IF (ASSOCIATED(Ewald%EIKZ)) THEN
               !   nx=SIZE(Ewald%EIKZ,1)
               !   ny=SIZE(Ewald%EIKZ,2)
               !   ALLOCATE(Ewald1%EIKZ(nx,ny))
               !   Ewald1%EIKZ=Ewald%EIKZ
               !END IF
               
               IF (ASSOCIATED(Ewald%KVEC)) THEN
                  nx=SIZE(Ewald%KVEC)
                  ALLOCATE(Ewald1%KVEC(nx))
                  Ewald1%KVEC=Ewald%KVEC
               END IF
               
               !IF (ASSOCIATED(Ewald%RVEC)) THEN
               !   nx=SIZE(Ewald%RVEC)
               !   ALLOCATE(Ewald1%RVEC(nx))
               !   Ewald1%RVEC=Ewald%RVEC
               !END IF
               
               IF (ASSOCIATED(Ewald%KVector)) THEN
                  nx=SIZE(Ewald%KVector)
                  ALLOCATE(Ewald1%KVector(nx))
                  Ewald1%KVector=Ewald%KVector
               END IF
            END IF
            
         END IF
         
         IF (ASSOCIATED(Potential%SW)) THEN
         
            ALLOCATE(Potential1%SW)
            
            SW=>Potential%SW
            SW1=>Potential1%SW
            
            SW1%TableLocationPP=SW%TableLocationPP
            SW1%SizePP=SW%SizePP
            SW1%SizeF3=SW%SizeF3
            SW1%RangePP=SW%RangePP
            SW1%RangeF3=SW%RangeF3
            SW1%IncrPPinv=SW%IncrPPinv
            SW1%IncrF3inv=SW%IncrF3inv
            SW1%PP=SW%PP
            SW1%F3=SW%F3
            SW1%epssqrt=SW%epssqrt
            SW1%lam=SW%lam
            SW1%ConvFacPP=SW%ConvFacPP
            SW1%EnabledPP=SW%EnabledPP
            
            SW1%A=SW%A
            SW1%B=SW%B
            SW1%g=SW%g
            SW1%q=SW%q
            SW1%epsln=SW%epsln
            SW1%sigma=SW%sigma
            SW1%ac=SW%ac
            SW1%lmd=SW%lmd
            SW1%gamma=SW%gamma
            SW1%cutoff=SW%cutoff
         END IF
         
         !WRITE(6,*) "Following potentials are associated:"
         !WRITE(6,*) "LJ:",ASSOCIATED(Potential1%LJ)
         !WRITE(6,*) "EAM:",ASSOCIATED(Potential1%EAM)
         !WRITE(6,*) "Coulomb:",ASSOCIATED(Potential1%Coulomb)
         !WRITE(6,*) "Buckingham:",ASSOCIATED(Potential1%Buckingham)
         !WRITE(6,*) "SW:",ASSOCIATED(Potential1%SW)
         !WRITE(6,*) "Tersoff:",ASSOCIATED(Potential1%Tersoff)
         !WRITE(6,*) "FENE:",ASSOCIATED(Potential1%FENE)
         !WRITE(6,*) "AIREBO:",ASSOCIATED(Potential1%AIREBO)
         !WRITE(6,*) "MEAM:",ASSOCIATED(Potential1%MEAM)
         !WRITE(6,*) "TB:",ASSOCIATED(Potential1%TB)
      ELSE
         RETURN
      END IF
   END SUBROUTINE CopyPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!    SUBROUTINE CopyALSection(AL,AL1,r1,r2,q1,q2) !copy AL into AL1
!       IMPLICIT NONE
!       !only copy selection of AL
!       !unlike CopyALFull AL1 can be of different size
!       TYPE(SystemContainer), POINTER :: AL,AL1
!       INTEGER :: r1,r2 !range in AL copied to AL1
!       INTEGER :: q1,q2 !range in AL1 into which AL is copied
!       
!       
!    END SUBROUTINE CopyALSection
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ResizeToLargerAL(AL,NAtoms)
   !resize an atom list to contain more atoms
   !NAtoms is the new size
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: NAtoms
      REAL(dp), DIMENSION(:), POINTER :: arrdp
      REAL(sp), DIMENSION(:), POINTER :: arrsp
      INTEGER, DIMENSION(:), POINTER :: arrint
      LOGICAL, DIMENSION(:), POINTER :: arrlog
      REAL(dp), DIMENSION(:), POINTER :: arrdp_ptr
      REAL(sp), DIMENSION(:), POINTER :: arrsp_ptr
      INTEGER, DIMENSION(:), POINTER :: arrint_ptr
      LOGICAL, DIMENSION(:), POINTER :: arrlog_ptr
      
      !check if more than one image present
      IF (ASSOCIATED(AL%PrevNeigh) .OR. ASSOCIATED(AL%NextNeigh)) &
         WRITE(6,*) "Warning>> AL found to have neighbors in ResizeToLargerAL"
         
      !check if LL and VL present
      IF (ASSOCIATED(AL%LL)) WRITE(6,*) "Warning>> AL has LL in ResizeToLargerAL"
      IF (ASSOCIATED(AL%LL)) WRITE(6,*) "Warning>> AL has VL in ResizeToLargerAL"
      
      !create new arrays
      IF (ASSOCIATED(AL%AtomCharge)) CALL resizedp(AL%AtomCharge,NAtoms,AL%NAtoms)
      IF (ASSOCIATED(AL%AtomSpecies)) CALL resizeint(AL%AtomSpecies,NAtoms,AL%NAtoms)
      IF (ASSOCIATED(AL%AtomCoord)) CALL resizedp(AL%AtomCoord,3*NAtoms,3*AL%NAtoms)
      IF (ASSOCIATED(AL%AtomVelocity)) CALL resizedp(AL%AtomVelocity,3*NAtoms,3*AL%NAtoms)
      IF (ASSOCIATED(AL%AtomDriftVelocity)) CALL resizedp(AL%AtomDriftVelocity,3*NAtoms,3*AL%NAtoms)
      IF (ASSOCIATED(AL%AtomMass)) CALL resizedp(AL%AtomMass,3*NAtoms,3*AL%NAtoms)
      IF (ASSOCIATED(AL%AtomForce)) CALL resizedp(AL%AtomForce,3*NAtoms,3*AL%NAtoms)
      IF (ASSOCIATED(AL%AtomShiftCoord)) CALL resizedp(AL%AtomShiftCoord,3*NAtoms,3*AL%NAtoms)
      IF (ASSOCIATED(AL%AtomIsMoving)) CALL resizelog(AL%AtomIsMoving,3*NAtoms,3*AL%NAtoms)
      IF (ASSOCIATED(AL%ChOSTangent)) CALL resizedp(AL%ChOSTangent,3*NAtoms,3*AL%NAtoms)
      IF (ASSOCIATED(AL%ChOSForce)) CALL resizedp(AL%ChOSForce,3*NAtoms,3*AL%NAtoms)
      
      AL%NAtoms=NAtoms
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE resizedp(ptr,rn,ro)
      !ptr should be AL%element pointer, rn and ro are the new and old size
         IMPLICIT NONE
         REAL(dp), DIMENSION(:), POINTER :: ptr
         INTEGER :: rn,ro,r
         
         r=MIN(rn,ro)
         ALLOCATE(arrdp(rn))
         arrdp_ptr=>ptr
         arrdp(1:r)=arrdp_ptr(1:r)
         ptr=>arrdp
         DEALLOCATE(arrdp_ptr)
      END SUBROUTINE resizedp
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE resizeint(ptr,rn,ro)
      !ptr should be AL%element pointer, rn and ro are the new and old size
         IMPLICIT NONE
         INTEGER, DIMENSION(:), POINTER :: ptr
         INTEGER :: rn,ro,r
         
         r=MIN(rn,ro)
         ALLOCATE(arrint(rn))
         arrint_ptr=>ptr
         arrint(1:r)=arrint_ptr(1:r)
         ptr=>arrint
         DEALLOCATE(arrint_ptr)
      END SUBROUTINE resizeint
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE resizelog(ptr,rn,ro)
      !ptr should be AL%element pointer, rn and ro are the new and old size
         IMPLICIT NONE
         LOGICAL, DIMENSION(:), POINTER :: ptr
         INTEGER :: rn,ro,r
         
         r=MIN(rn,ro)
         ALLOCATE(arrlog(rn))
         arrlog_ptr=>ptr
         arrlog(1:r)=arrlog_ptr(1:r)
         ptr=>arrlog
         DEALLOCATE(arrlog_ptr)
      END SUBROUTINE resizelog
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE ResizeToLargerAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ResizeToLargerAL1(AL,NAtoms)
   !resize an atom list to contain more atoms -- older version
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1
      INTEGER :: NAtoms
      
      NULLIFY(AL1)
      CALL MakeSize(AL1,NAtoms)
      CALL Copy(AL,AL1,1,AL%NAtoms)
      CALL Delete(AL)
      AL=>AL1
      NULLIFY(AL1)
   END SUBROUTINE ResizeToLargerAL1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MakeSizeAL(AL,n,AddChOSvariablesInAL,imake)
      !create/modify AtomList where n is number of atoms
      !imake=0: make all
      !imake=2: make coordinates
      !imake=3: make velocity, driftvelocity and force
      !imake=5: make mass, species and charge
      !imake=7: make neighborlist -- not in use
      !imake=11: make shift coord
      !imake=13: make is moving
      !imake=17: make potential
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: n,imake1
      INTEGER, OPTIONAL :: imake
      LOGICAL :: AddChOSvariablesInAL1
      LOGICAL, OPTIONAL :: AddChOSvariablesInAL
      
      IF (.NOT. ASSOCIATED(AL)) ALLOCATE(AL)
      
      IF (PRESENT(AddChOSvariablesInAL)) THEN
         AddChOSvariablesInAL1=AddChOSvariablesInAL
      ELSE
         AddChOSvariablesInAL1=.FALSE.
      END IF
      
      imake1=0
      IF (PRESENT(imake)) imake1=imake
      
      AL%PotentialEnergy=0._dp
      
      AL%KineticEnergy=0._dp
      AL%NAtoms=n
      AL%NMove=n !note NMove is 3D
         !all atoms in all directions are moving
      !AL%NSpeciesType=0
      
      AL%BoundaryCondX="periodic"
      AL%BoundaryCondY="periodic"
      AL%BoundaryCondZ="periodic"
      
      IF (MOD(imake1,5)==0) CALL MakeSize(AL%AtomCharge,n)
      IF (MOD(imake1,5)==0) CALL MakeSize(AL%AtomSpecies,n)
      
      IF (MOD(imake1,2)==0) CALL MakeSize(AL%AtomCoord,3*n)
      IF (MOD(imake1,5)==0) CALL MakeSize(AL%AtomMass,3*n)
      IF (MOD(imake1,3)==0) CALL MakeSize(AL%AtomForce,3*n)
      IF (MOD(imake1,13)==0) CALL MakeSize(AL%AtomIsMoving,3*n)
      
      IF (MOD(imake1,3)==0) CALL MakeSize(AL%AtomVelocity,3*n)
      IF (MOD(imake1,3)==0) CALL MakeSize(AL%AtomDriftVelocity,3*n)
      IF (MOD(imake1,3)==0) AL%AtomDriftVelocity=0._dp
      IF (MOD(imake1,11)==0) CALL MakeSize(AL%AtomShiftCoord,3*n)
         
      IF (AddChOSvariablesInAL1) THEN
         CALL MakeSize(AL%ChOSTangent,3*n)
         CALL MakeSize(AL%ChOSForce,3*n)
      !ELSE
      END IF
   END SUBROUTINE MakeSizeAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MakeSize1Dint(arr1,n)
      IMPLICIT NONE
      INTEGER, DIMENSION(:), POINTER :: arr1
      INTEGER :: n
      LOGICAL :: IsAllocated
      
      IsAllocated=ASSOCIATED(arr1)
      
      IF (IsAllocated) THEN
         IF (SIZE(arr1)/=n) THEN
            DEALLOCATE(arr1)
            IsAllocated=.FALSE.
         END IF
      END IF
      IF (.NOT. IsAllocated) THEN
         ALLOCATE(arr1(n))
      END IF
      
      arr1=0
   END SUBROUTINE MakeSize1Dint
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MakeSize2Dint(arr1,n)
      IMPLICIT NONE
      INTEGER, DIMENSION(:,:), POINTER :: arr1
      INTEGER :: n
      LOGICAL :: IsAllocated
      
      IsAllocated=ASSOCIATED(arr1)
      IF (IsAllocated) THEN
         IF (SIZE(arr1)/=n*n) THEN
            DEALLOCATE(arr1)
            IsAllocated=.FALSE.
         END IF
      END IF
      IF (.NOT. IsAllocated) THEN
         ALLOCATE(arr1(n,n))
      END IF
      
      arr1=0
   END SUBROUTINE MakeSize2Dint
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MakeSize1Dsp(arr1,n)
      IMPLICIT NONE
      REAL(sp), DIMENSION(:), POINTER :: arr1
      INTEGER :: n
      LOGICAL :: IsAllocated
      
      IsAllocated=ASSOCIATED(arr1)
      IF (IsAllocated) THEN
         IF (SIZE(arr1)/=n) THEN
            DEALLOCATE(arr1)
            IsAllocated=.FALSE.
         END IF
      END IF
      IF (.NOT. IsAllocated) THEN
         ALLOCATE(arr1(n))
      END IF
      
      arr1=0._dp
   END SUBROUTINE MakeSize1Dsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MakeSize1Ddp(arr1,n)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), POINTER :: arr1
      INTEGER :: n
      LOGICAL :: IsAllocated
      
      IsAllocated=ASSOCIATED(arr1)
      IF (IsAllocated) THEN
         IF (SIZE(arr1)/=n) THEN
            DEALLOCATE(arr1)
            IsAllocated=.FALSE.
         END IF
      END IF
      IF (.NOT. IsAllocated) THEN
         ALLOCATE(arr1(n))
      END IF
      
      arr1=0._dp
   END SUBROUTINE MakeSize1Ddp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MakeSize2Ddp(arr1,n,n1)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:,:), POINTER :: arr1
      INTEGER, INTENT(IN) :: n
      INTEGER, OPTIONAL :: n1
      INTEGER :: n2
      LOGICAL :: IsAllocated
      
      IF (PRESENT(n1)) THEN
         n2=n1
      ELSE
         n2=n
      END IF
      
      IsAllocated=ASSOCIATED(arr1)
      IF (IsAllocated) THEN
         IF (SIZE(arr1)/=n*n2) THEN
            DEALLOCATE(arr1)
            IsAllocated=.FALSE.
         END IF
      END IF
      IF (.NOT. IsAllocated) THEN
         ALLOCATE(arr1(n,n2))
      END IF
      
      arr1=0._dp
   END SUBROUTINE MakeSize2Ddp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MakeSize1Dlog(arr1,n)
      IMPLICIT NONE
      LOGICAL, DIMENSION(:), POINTER :: arr1
      INTEGER :: n
      LOGICAL :: IsAllocated
      
      IsAllocated=ASSOCIATED(arr1)
      IF (IsAllocated) THEN
         IF (SIZE(arr1)/=n) THEN
            DEALLOCATE(arr1)
            IsAllocated=.FALSE.
         END IF
      END IF
      IF (.NOT. IsAllocated) THEN
         ALLOCATE(arr1(n))
      END IF
      
      arr1=.TRUE.
   END SUBROUTINE MakeSize1Dlog
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MakeSize2Dlog(arr1,n)
      IMPLICIT NONE
      LOGICAL, DIMENSION(:,:), POINTER :: arr1
      INTEGER :: n
      LOGICAL :: IsAllocated
      
      IsAllocated=ASSOCIATED(arr1)
      IF (IsAllocated) THEN
         IF (SIZE(arr1)/=n*n) THEN
            DEALLOCATE(arr1)
            IsAllocated=.FALSE.
         END IF
      END IF
      IF (.NOT. IsAllocated) THEN
         ALLOCATE(arr1(n,n))
      END IF
      
      arr1=.TRUE.
   END SUBROUTINE MakeSize2Dlog
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MakeSize1Dchar(arr1,n)
      IMPLICIT NONE
      CHARACTER(*), DIMENSION(:), POINTER :: arr1
      INTEGER :: n
      LOGICAL :: IsAllocated
      
      IsAllocated=ASSOCIATED(arr1)
      IF (IsAllocated) THEN
         IF (SIZE(arr1)/=n) THEN
            DEALLOCATE(arr1)
            IsAllocated=.FALSE.
         END IF
      END IF
      IF (.NOT. IsAllocated) THEN
         ALLOCATE(arr1(n))
      END IF
      
   END SUBROUTINE MakeSize1Dchar
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MakeSize2Dchar(arr1,n)
      IMPLICIT NONE
      CHARACTER(*), DIMENSION(:,:), POINTER :: arr1
      INTEGER :: n
      LOGICAL :: IsAllocated
      
      IsAllocated=ASSOCIATED(arr1)
      IF (IsAllocated) THEN
         IF (SIZE(arr1)/=n*n) THEN
            DEALLOCATE(arr1)
            IsAllocated=.FALSE.
         END IF
      END IF
      IF (.NOT. IsAllocated) THEN
         ALLOCATE(arr1(n,n))
      END IF
      
   END SUBROUTINE MakeSize2Dchar
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   SUBROUTINE MakeSizeLL(LL,n)

      IMPLICIT NONE
      TYPE(LinkedListContainer), POINTER :: LL
      INTEGER, INTENT(IN) :: n

      IF (.NOT. ASSOCIATED(LL)) ALLOCATE(LL)

      LL%NAtoms=n

      CALL MakeSize(LL%AtomCellIndx,n)
      
   END SUBROUTINE MakeSizeLL

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   SUBROUTINE MakeSizeVL(VL,n,nmax)

      IMPLICIT NONE
      TYPE(VerletListContainer), POINTER :: VL
      INTEGER, INTENT(IN) :: n,nmax
      
      IF (.NOT. ASSOCIATED(VL)) ALLOCATE(VL)

      VL%NAtoms=n
      
      IF (n>nmax) RETURN
      
      IF (VL%MaxAtomPerAtom==0) VL%MaxAtomPerAtom=80 !before 80

      CALL MakeSize(VL%ListRange,n)
      CALL MakeSize(VL%OrigCoord,3*n) !earlier we had OrigCoord of fixed size, pre-allocated
      
      VL%ListRange=0
      !VL%ListDomainAtom=.TRUE.
      VL%OrigCoord=0._dp

   END SUBROUTINE MakeSizeVL

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SortXYZ(AL,dimn,sgn)
      !sort only along 1 dimension of XYZ
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: dimn,NAtoms,sgn
      
      WRITE(UNIT=*,FMT='(" ...sorting spatial dimension ",I1," in")',ADVANCE="NO") dimn
      
      NAtoms=AL%NAtoms
      SELECT CASE(dimn)
      CASE(1); CALL SortAL(AL,AL%AtomCoord(1:3*NAtoms-2:3),sgn)
      CASE(2); CALL SortAL(AL,AL%AtomCoord(2:3*NAtoms-1:3),sgn)
      CASE(3); CALL SortAL(AL,AL%AtomCoord(3:3*NAtoms  :3),sgn)
      CASE DEFAULT
         WRITE(*,*) "$Err>> Unknown dimension has been passed to sort"
         STOP
      END SELECT
   END SUBROUTINE SortXYZ
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SortPXYZ()
   END SUBROUTINE SortPXYZ
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SortSpeciesType()
   END SUBROUTINE SortSpeciesType
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SortForce()
   END SUBROUTINE SortForce
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SortAL(AL,arr_real,sgn) 
      !sort AL
      !arr_real is an array in AL
      !sgn=1 ascending, sgn=-1 descending
      !currently all atoms are sorted
      IMPLICIT NONE
         !for the atom # based on the values obtained from sortobj
         !sortobj is a user fn that returns obj value for AL
         !e.g. sort along z ascending returns AL%AtomCoord(3:3*NAtoms,3)
         ! while sort along z descending returns -1*AL%AtomCoord(3:3*NAtoms,3)
         !sort based on Atomic# can be written in a similar fashion
         !NOTE : IF 
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: sgn,NAtoms,i,j
      REAL(dp), DIMENSION(:) :: arr_real
      INTEGER, DIMENSION(:), POINTER :: indx
      REAL(dp), DIMENSION(:), POINTER :: arr_tmp_r
      INTEGER, DIMENSION(:), POINTER :: arr_tmp_i
      LOGICAL, DIMENSION(:), POINTER :: arr_tmp_l

      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(*,*) "$Err>> AL not associated, cannot be sorted"
         STOP
      END IF
      
      NAtoms=AL%NAtoms
      IF (NAtoms==0) THEN
         WRITE(*,*) "$Err>> AL should contain atoms before sorting"
         STOP
      END IF
      IF (SIZE(arr_real)/=NAtoms) THEN
         WRITE(*,*) "$Err>> Array to be sorted should be of same size as number of atoms"
         STOP
      END IF
      
      IF (ABS(sgn)/=1) THEN
         WRITE(*,*) "$Err>> Value of sgn passed to sort is unclear"
         STOP
      END IF
      
      CALL MakeSize(arr_tmp_r,SIZE(arr_real))
      CALL MakeSize(indx,AL%NAtoms)
      arr_tmp_r=arr_real
      DO i=1,NAtoms
         indx(i)=i
      END DO
      
      IF (sgn>1) THEN
         WRITE(UNIT=*,FMT='("ascending order using")',ADVANCE="NO")
      ELSE
         WRITE(UNIT=*,FMT='("descending order using")',ADVANCE="NO")
      END IF
      IF (NAtoms<20) THEN
         WRITE(UNIT=*,FMT='(" Regular sort")') 
         CALL sort_pick_abhijit(arr_tmp_r,indx)
      ELSEIF (NAtoms<1000) THEN
         WRITE(UNIT=*,FMT='(" Heap sort")') 
         CALL sort_heap_abhijit(arr_tmp_r,indx)
      ELSE
         WRITE(UNIT=*,FMT='(" Quick sort")') 
         CALL sort_quick_abhijit(arr_tmp_r,indx)
      END IF
      
      !check
      !DO i=1,NAtoms
      !   WRITE(UNIT=*,FMT='(i3,f8.3,"       ",i3,f8.3)') i,arr_real(i),indx(i),arr_tmp_r(i)
      !END DO
      DO i=1,NAtoms
         j=indx(i)
         IF (arr_real(j)/=arr_tmp_r(i)) THEN
            WRITE(*,*) ""
            WRITE(*,*) "$Err>> sorted index is incorrect"
            STOP
         END IF
      END DO
      
      !rearrange using the updated indx & store 
      !the result in arr if real, in arr1 if integer, in arr2 if logical
      CALL MakeSize(arr_tmp_r,NAtoms)
      CALL MakeSize(arr_tmp_i,NAtoms)
      
      DO i=1,NAtoms
         j=indx(i)
         IF (sgn>1) THEN
            arr_tmp_r(i)=AL%AtomCharge(j)
            arr_tmp_i(i)=AL%AtomSpecies(j)
         ELSE
            arr_tmp_r(NAtoms+1-i)=AL%AtomCharge(j)
            arr_tmp_i(NAtoms+1-i)=AL%AtomSpecies(j)
         END IF
      END DO
      AL%AtomCharge=arr_tmp_r
      AL%AtomSpecies=arr_tmp_i
      
      CALL MakeSize(arr_tmp_r,3*NAtoms)
      CALL MakeSize(arr_tmp_l,3*NAtoms)
      !--------------
      DO i=1,NAtoms
         j=indx(i)
         IF (sgn>1) THEN
            arr_tmp_r(3*i-2:3*i)=AL%AtomCoord(3*j-2:3*j)
         ELSE
            arr_tmp_r(3*(NAtoms+1-i)-2:3*(NAtoms+1-i))=AL%AtomCoord(3*j-2:3*j)
         END IF
      END DO
      AL%AtomCoord=arr_tmp_r
      !--------------
      DO i=1,NAtoms
         j=indx(i)
         IF (sgn>1) THEN
            arr_tmp_r(3*i-2:3*i)=AL%AtomVelocity(3*j-2:3*j)
         ELSE
            arr_tmp_r(3*(NAtoms+1-i)-2:3*(NAtoms+1-i))=AL%AtomVelocity(3*j-2:3*j)
         END IF
      END DO
      AL%AtomVelocity=arr_tmp_r
      !--------------
      DO i=1,NAtoms
         j=indx(i)
         IF (sgn>1) THEN
            arr_tmp_r(3*i-2:3*i)=AL%AtomMass(3*j-2:3*j)
         ELSE
            arr_tmp_r(3*(NAtoms+1-i)-2:3*(NAtoms+1-i))=AL%AtomMass(3*j-2:3*j)
         END IF
      END DO
      AL%AtomMass=arr_tmp_r
      !--------------
      DO i=1,NAtoms
         j=indx(i)
         IF (sgn>1) THEN
            arr_tmp_r(3*i-2:3*i)=AL%AtomForce(3*j-2:3*j)
         ELSE
            arr_tmp_r(3*(NAtoms+1-i)-2:3*(NAtoms+1-i))=AL%AtomForce(3*j-2:3*j)
         END IF
      END DO
      AL%AtomForce=arr_tmp_r
      !--------------
      DO i=1,NAtoms
         j=indx(i)
         IF (sgn>1) THEN
            arr_tmp_r(3*i-2:3*i)=AL%AtomShiftCoord(3*j-2:3*j)
         ELSE
            arr_tmp_r(3*(NAtoms+1-i)-2:3*(NAtoms+1-i))=AL%AtomShiftCoord(3*j-2:3*j)
         END IF
      END DO
      AL%AtomShiftCoord=arr_tmp_r
      !--------------
      DO i=1,NAtoms
         j=indx(i)
         IF (sgn>1) THEN
            arr_tmp_l(3*i-2:3*i)=AL%AtomIsMoving(3*j-2:3*j)
         ELSE
            arr_tmp_l(3*(NAtoms+1-i)-2:3*(NAtoms+1-i))=AL%AtomIsMoving(3*j-2:3*j)
         END IF
      END DO
      AL%AtomIsMoving=arr_tmp_l
      !--------------
      DEALLOCATE(arr_tmp_r,arr_tmp_i,arr_tmp_l,indx)
   END SUBROUTINE SortAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !FUNCTION GetSpecies(AL,atomicnumber) !,ForceAdd1)
   FUNCTION GetSpecies(atomicnumber) !,ForceAdd1)
      !if atomicnumber is present in AL%SpeciesDirectory
      !GetSpecies returns its position
      !else a new spot is created for atomicnumber
      IMPLICIT NONE
      !TYPE(SystemContainer), POINTER :: AL
      INTEGER :: atomicnumber,GetSpecies,i
      !INTEGER, DIMENSION(:), POINTER :: an
      !LOGICAL, OPTIONAL :: ForceAdd1
      !LOGICAL :: ForceAdd
      
      !IF (PRESENT(ForceAdd1)) THEN
      !   ForceAdd=ForceAdd1
      !   WRITE(*,*) "Force add is not allowed"
      !   STOP
      !ELSE
      !   ForceAdd=.FALSE.
      !END IF
      !IF (.NOT. ForceAdd .AND. .NOT. ASSOCIATED(SpeciesDirectory)) THEN
      IF (.NOT. ASSOCIATED(SpeciesDirectory_global)) THEN
         WRITE(*,*) "$Err>>Species array for Atom List not initialized"
         STOP
      END IF
      
      !DO i=1,AL%NSpeciesType
      DO i=1,NSpecies_global
         IF (SpeciesDirectory_global(i)==atomicnumber) THEN
            GetSpecies=i
            RETURN
         END IF
      END DO
      
      !older version -- where new species could be added
      !CALL MakeSize(an,AL%NSpeciesType)
      !an=AL%SpeciesDirectory
      !AL%NSpeciesType=AL%NSpeciesType+1
      !CALL MakeSize(AL%SpeciesDirectory,AL%NSpeciesType)
      !AL%SpeciesDirectory(1:AL%NSpeciesType-1) = an
      !AL%SpeciesDirectory(AL%NSpeciesType)=atomicnumber
      
      !DEALLOCATE(an)
   END FUNCTION GetSpecies
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CompareALToAL(AL1,AL2,xtol,dx1,match_algorithm)
   
      !returns TRUE if AL1 is identical to AL2
      !if AL1 and AL2 are minimized then even better chances of finding match
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL1,AL2
      REAL(dp), OPTIONAL :: xtol,dx1
      CHARACTER(len=10), OPTIONAL :: match_algorithm
      REAL(dp) :: xtol1
      CHARACTER(len=10) :: match_algorithm1
      LOGICAL :: CompareALToAL
      
      IF (PRESENT(xtol)) THEN
         xtol1=xtol
      ELSE
         xtol1=0.5_dp !Angstroms
      END IF
      
      IF (PRESENT(match_algorithm)) THEN
         match_algorithm1=match_algorithm
      ELSE
         match_algorithm1="quick"
      END IF
      
      SELECT CASE (match_algorithm1(1:5))
      CASE('quick')
         IF (present(dx1)) THEN
            CompareALToAL= CompareALToALQuick(AL1,AL2,xtol1,dx1)
         ELSE
            CompareALToAL= CompareALToALQuick(AL1,AL2,xtol1)
         END IF
      CASE('rigor')
         CompareALToAL= CompareALToALRigourous(AL1,AL2,xtol1)
      CASE DEFAULT
         WRITE(*,*) "$Err>> Unknown option provided to CompareALToAL"
         STOP
      END SELECT
   END FUNCTION CompareALToAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CompareALToALRigourous(AL1,AL2,xtol)
   !AL1 is the reference state and may comprise of more than one image
   !AL2 is the state which has be compared to the reference state and may comprise of more than one image
   !rigourous comparison of atoms in AL1 and AL2
      IMPLICIT NONE
      REAL(dp) :: xtol
      TYPE(SystemContainer), POINTER :: AL1,AL2
      LOGICAL :: CompareALToALRigourous
      
      WRITE(*,*) "Not implemented yet"
      STOP
      CompareALToALRigourous=.FALSE.
   END FUNCTION CompareALToALRigourous
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CompareALToALQuick(AL1,AL2,xtol,dx1)
   !AL1 is the reference state and may comprise of more than one image
   !AL2 is the state which has be compared to the reference state and may comprise of more than one image
   !Direct comparison of AL atom positions
   !useful for MD simulations
      IMPLICIT NONE
      REAL(dp) :: dx,xtol
      REAL(dp), OPTIONAL :: dx1
      INTEGER :: NAtoms
      TYPE(SystemContainer), POINTER :: AL1,AL2
      LOGICAL :: CompareALToALQuick
      
      NAtoms=AL1%NAtoms
      CompareALToALQuick= NAtoms==AL2%NAtoms
      IF (.NOT. CompareALToALQuick) RETURN
      
      CompareALToALQuick= ALL(AL1%AtomSpecies==AL2%AtomSpecies)
      
      !dx=MAXVAL(ABS(AL1%AtomCoord-AL2%AtomCoord))
      
      dx=PBCdist()
      CompareALToALQuick=dx<xtol
      !WRITE(6,*) TRIM(TxtHeader)//"CompareAL>> Difference found=",dx,xtol,CompareALToALQuick
      
      IF (PRESENT(dx1)) dx1=dx
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      FUNCTION PBCdist()
   
         IMPLICIT NONE
         REAL(dp), DIMENSION(3*MaxNAtoms) :: BoxSize,vec1
         REAL(dp) :: PBCdist
         !INTEGER :: i
         
         IF (MaxNAtoms<NAtoms) THEN
            WRITE(6,*) "Err>> Cannot compare the two ALs. Increase MaxNAtomsAL"
            STOP
         END IF
         BoxSize(1:3*NAtoms)=PACK(SPREAD(AL1%BoxSize(1:3)-AL1%BoxSize(4:6),2,NAtoms),.TRUE.)
      
         vec1(1:3*NAtoms)=ABS(AL1%AtomCoord(1:3*NAtoms)-AL2%AtomCoord(1:3*NAtoms))
         vec1(1:3*NAtoms)=MIN(BoxSize(1:3*NAtoms)-vec1(1:3*NAtoms),vec1(1:3*NAtoms))
         
         !DO i=1,NAtoms
         !   IF (ANY(vec1(3*i-2:3*i)>xtol)) THEN
         !      WRITE(6,*) "Atom index",i
         !      WRITE(6,*) "AL1Coord:",AL1%AtomCoord(3*i-2:3*i)
         !      WRITE(6,*) "AL2Coord:",AL2%AtomCoord(3*i-2:3*i)
         !   END IF
         !END DO
         PBCdist=MAXVAL(vec1(1:3*NAtoms))
      
      END FUNCTION PBCdist
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   END FUNCTION CompareALToALQuick
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE db_manipulate
