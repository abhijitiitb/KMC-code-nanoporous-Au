MODULE KMCRelax
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur

   USE KMC_VARIABLE_TYPE
   USE VARIABLE_TYPE
   USE db_manipulate
   USE KMCUtilities, ONLY : GetPBCCoord,GetPBCSpacing,PrintAtoms
   USE KMCDbManipulate, ONLY : RefreshCell
   USE OptimizationAL
   USE NEB_package
   
   IMPLICIT NONE
   REAL, PRIVATE :: a=0.,UnitCell(3)=0.
   REAL, ALLOCATABLE, DIMENSION(:,:), PRIVATE :: UnitCellAtomCoord(:,:)
   INTEGER, PRIVATE :: SizeUnitCellAtomCoord=0
   TYPE(SystemContainer), POINTER, PRIVATE :: AL
   LOGICAL, PRIVATE :: Initialized=.FALSE.
   SAVE

   CONTAINS

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RelaxFullKMCSystem()
   !converts the KMC system into an AMD atom list and then does the optimization
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      INTEGER :: idx,ix,iy,iz,CellIndex(3),NAtoms,NMove(3)
      INTEGER :: i,NoImages(3),MinDistanceIndx,idy
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      REAL :: NewCoord(3),MinDistance,Distance,TmpRealDP=0.,vec(3)
      REAL :: TmpReal=0.,coord(3),MinBox(3),MaxBox(3)
      CHARACTER(len=100) :: TmpName
      CHARACTER :: CrapChar
      
      
      IF (KMCRelaxOption==0) RETURN

      IF (KMCRelaxOption==2) THEN !snap atoms to lattice
         atom=>KMC_AL
         DO idx=1,NMoveKMC
            NewCoord=atom%coord
            NoImages=INT(NewCoord/UnitCell)
            NewCoord=MOD(NewCoord,UnitCell)
            WHERE(NewCoord<0.) NoImages=NoImages-1
            WHERE(NewCoord<0.) NewCoord=NewCoord+UnitCell
            MinDistance=a
            MinDistanceIndx=0
            DO idy=1,SizeUnitCellAtomCoord
               vec=UnitCellAtomCoord(idy,:)-NewCoord
               Distance=SQRT(DOT_PRODUCT(vec,vec))
               IF (Distance<MinDistance) THEN
                  MinDistance=Distance
                  MinDistanceIndx=idy
               END IF
            END DO
            NewCoord=UnitCellAtomCoord(MinDistanceIndx,:)
            NewCoord=NewCoord+REAL(NoImages)*UnitCell
            !WRITE(UNIT=UnitScrap,FMT='(3f10.3)') atom%Coord
            atom%coord=NewCoord
            IF (MinDistance>KMCMaxLatticeSpacing) THEN
               WRITE(*,*) 'Ini>> (WARNING) Atom shifted by ',MinDistance
               WRITE(UnitScrap,*) 'Atom shifted by ',MinDistance
            END IF
            atom=>atom%NextNeigh
         END DO
      ELSEIF (KMCRelaxOption==1) THEN !energy minimization of affected cells
         IF (.NOT. Initialized) THEN
            NULLIFY(AL)
            CALL MakeSize(AL,MaxNAtoms)
            CALL PotentialInitialize(AL)
            Initialized=.TRUE.
         END IF
         NAtoms=0
         NMove=0
         MaxBox=KMCBoxSize
         MinBox=0.
         atom=>KMC_AL
         AtomCoord=>AL%AtomCoord
         AtomIsMoving=>AL%AtomIsMoving
         AtomSpecies=>AL%AtomSpecies
         DO WHILE (ASSOCIATED(atom))
            NAtoms=NAtoms+1
            coord=atom%Coord
            AtomCoord(3*NAtoms-2:3*NAtoms)=coord
            AtomSpecies(NAtoms)=atom%Species
            IF (atom%IsMoving) THEN
               AtomIsMoving(3*NAtoms-2:3*NAtoms)=.TRUE.
               NMove=NMove+1
            ELSE
               AtomIsMoving(3*NAtoms-2:3*NAtoms)=.FALSE.
            END IF
            atom=>atom%NextNeigh
         END DO
         AL%NAtoms=NAtoms
         AL%NMove=NMove
         AL%BoxSize(1:3)=REAL(MaxBox,dp)
         AL%BoxSize(4:6)=REAL(MinBox,dp)
         CALL AddVerletList(AL,NAtoms=NAtoms)
         WRITE(*,*) "Made global optimization go slow ...",KMCStep,KMCBlock
         !TmpName="KMCGlobalRelax1.xyz"
         !CALL WriteXYZ(AL,TmpName)
         CALL OptimizeAL(AL,4,NAtoms=NAtoms,ftol=1.e-4_dp,gtol=1.e-3_dp,xtol=1.e-2_dp, &
            ITMAX=50,MaxDisplacement=0.05_dp)
         !CALL OptimizeAL(AL,3,NAtoms=NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
         !   ITMAX=500,MaxDisplacement=0.05_dp)
         CALL PBC(AL)
         
         !now read the optimized positions
         NAtoms=0
         atom=>KMC_AL
         DO WHILE (ASSOCIATED(atom))
            NAtoms=NAtoms+1
            atom%Coord=AtomCoord(3*NAtoms-2:3*NAtoms)
            CALL RefreshCell(atom)
            atom=>atom%NextNeigh
         END DO
         !NULLIFY(AtomCoord)
         !NULLIFY(AtomIsMoving)
         
         !CALL Delete(AL) !not needed
         !Initialized=.FALSE.
         
      END IF
   END SUBROUTINE RelaxFullKMCSystem
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RelaxKMCSystem(CenterCoord)
   !Performs local optimization/lattice fit for a subsystem that is a part of the large system
   !The subsystem consists of cells that listed in AffectedCells
   !CenterCoord is the coordinate around which the sub-system is carved
   !CenterCoord is used to move atoms relative to the center
   ! at the end of the optimization the atoms are moved back by the distance CenterCoord
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCAtomList), POINTER :: KMC_AL1,KMC_AL2
      INTEGER :: idx,ix,iy,iz,CellIndex(3),NAtoms,NMove(3),NAtoms1,NAtoms2,NAtoms3
      INTEGER :: i,NoImages(3),MinDistanceIndx,idy
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      REAL :: NewCoord(3),MinDistance,Distance,TmpRealDP=0.,vec(3)
      REAL :: TmpReal=0.,coord(3),MinBox(3),MaxBox(3),CenterCoord(3)
      CHARACTER(len=100) :: TmpName
      CHARACTER :: CrapChar
      LOGICAL :: CellIsMoving

      IF (KMCRelaxOption==0) RETURN

      IF (KMCRelaxOption==2) THEN !snap atoms to lattice
         atom=>KMC_AL
         DO idx=1,NMoveKMC
            NewCoord=atom%coord
            NoImages=INT(NewCoord/UnitCell)
            NewCoord=MOD(NewCoord,UnitCell)
            WHERE(NewCoord<0.) NoImages=NoImages-1
            WHERE(NewCoord<0.) NewCoord=NewCoord+UnitCell
            MinDistance=a
            MinDistanceIndx=0
            DO idy=1,SizeUnitCellAtomCoord
               vec=UnitCellAtomCoord(idy,:)-NewCoord
               Distance=SQRT(DOT_PRODUCT(vec,vec))
               IF (Distance<MinDistance) THEN
                  MinDistance=Distance
                  MinDistanceIndx=idy
               END IF
            END DO
            NewCoord=UnitCellAtomCoord(MinDistanceIndx,:)
            NewCoord=NewCoord+REAL(NoImages)*UnitCell
            !WRITE(UNIT=UnitScrap,FMT='(3f10.3)') atom%Coord
            atom%coord=NewCoord
            IF (MinDistance>KMCMaxLatticeSpacing) THEN
               WRITE(*,*) 'Ini>> (WARNING) Atom shifted by ',MinDistance
               WRITE(UnitScrap,*) 'Atom shifted by ',MinDistance
            END IF
            atom=>atom%NextNeigh
         END DO
      ELSEIF (KMCRelaxOption==1) THEN !energy minimization of affected cells
         !CALL MakeSize(AL,MaxNAtoms)
         IF (.NOT. Initialized) THEN
            NULLIFY(AL)
            CALL MakeSize(AL,MaxSlaveNAtoms)
            CALL PotentialInitialize(AL)
            Initialized=.TRUE.
         END IF
         NAtoms=0
         NMove=0
         MaxBox=0.
         MinBox=10000.
         AtomCoord=>AL%AtomCoord
         AtomIsMoving=>AL%AtomIsMoving
         AtomSpecies=>AL%AtomSpecies
         DO i=1,NAffectedCells
            CellIndex=AffectedCells(3*i-2:3*i)
            ix=CellIndex(1)
            iy=CellIndex(2)
            iz=CellIndex(3)
            CellIsMoving=IsCellAffected(ix,iy,iz)==1 .OR. IsCellAffected(ix,iy,iz)==2
            KMC_AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
            DO WHILE (ASSOCIATED(KMC_AL1))
               atom=>KMC_AL1%Atom
               NAtoms=NAtoms+1
               
               coord=GetPBCSpacing(atom%Coord,CenterCoord) !+KMCBoxSize*0.5 !so that CenterCoord moves to KMCBoxSize/2.
               AtomCoord(3*NAtoms-2:3*NAtoms)=coord
               AtomSpecies(NAtoms)=atom%Species
               MinBox=MIN(MinBox,coord)
               MaxBox=MAX(MaxBox,coord)
               IF (atom%IsMoving .AND. CellIsMoving) THEN
                  AtomIsMoving(3*NAtoms-2:3*NAtoms)=.TRUE.
                  NMove=NMove+1
               ELSE
                  AtomIsMoving(3*NAtoms-2:3*NAtoms)=.FALSE.
               END IF
               KMC_AL1=>KMC_AL1%NextNeigh
            END DO
         END DO
         
         AL%NAtoms=NAtoms
         AL%NMove=NMove
         !AL%BoxSize(1:3)=REAL(MaxBox+MaxBox-MinBox,dp)
         !AL%BoxSize(4:6)=REAL(MinBox-MaxBox+MinBox,dp)
         AL%BoxSize(1:3)=REAL(MaxBox+10.,dp)
         AL%BoxSize(4:6)=REAL(MinBox-10.,dp)
         
         TmpName="KMCLocalRelax1.xyz"
         CALL WriteXYZ(AL,TmpName)
         !stop
         !CALL PotentialInitialize(AL)
         CALL AddVerletList(AL,NAtoms=NAtoms)
         WRITE(*,*) "MADE OPTIMIZATION GO SLOW ...",KMCStep,KMCBlock
         CALL OptimizeAL(AL,4,NAtoms=NAtoms,ftol=1.e-4_dp,gtol=1.e-3_dp,xtol=1.e-2_dp, &
            ITMAX=30,MaxDisplacement=0.2_dp)
         !CALL OptimizeAL(AL,3,NAtoms=NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
         !   ITMAX=500,MaxDisplacement=0.2_dp)
         !CALL SYSTEM("cat KMCRel.xyz KMCLocalRelax1.xyz > tmp00.xyz")
         !CALL SYSTEM("mv tmp00.xyz KMCRel.xyz")
         TmpName="KMCLocalRelax2.xyz"
         CALL WriteXYZ(AL,TmpName)
         !READ(*,*) CrapChar
         !stop
         
         !now read the optimized positions
         NAtoms=0
         !NAtoms1=0
         !NAtoms2=0
         !NAtoms3=0
         DO i=1,NAffectedCells
            CellIndex=AffectedCells(3*i-2:3*i)
            ix=CellIndex(1)
            iy=CellIndex(2)
            iz=CellIndex(3)
            CellIsMoving=IsCellAffected(ix,iy,iz)==1 !values 0, 2 and 3 are not copied
            KMC_AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
            DO WHILE (ASSOCIATED(KMC_AL1))
               atom=>KMC_AL1%Atom
               NAtoms=NAtoms+1
               !IF (IsCellAffected(ix,iy,iz)==1) NAtoms1=NAtoms1+1
               !IF (IsCellAffected(ix,iy,iz)==2) NAtoms2=NAtoms2+1
               !IF (IsCellAffected(ix,iy,iz)==3) NAtoms3=NAtoms3+1
               IF (CellIsMoving) THEN
                  coord=AtomCoord(3*NAtoms-2:3*NAtoms)+CenterCoord !-KMCBoxSize*0.5
                  coord=GetPBCCoord(coord)
                  IF (ANY(ABS(GetPBCSpacing(coord,atom%Coord))>1.5)) THEN
                     WRITE(6,*) "Atom has moved too much",coord,atom%Coord
                     STOP
                  END IF
                  atom%Coord=coord
               END IF
               KMC_AL1=>KMC_AL1%NextNeigh
            END DO
         END DO
         
         DO i=1,NAffectedCells
            CellIndex=AffectedCells(3*i-2:3*i)
            ix=CellIndex(1)
            iy=CellIndex(2)
            iz=CellIndex(3)
            !CellIsMoving=IsCellAffected(ix,iy,iz)==1
            !IF (CellIsMoving) THEN
               KMC_AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
               DO WHILE (ASSOCIATED(KMC_AL1))
                  atom=>KMC_AL1%Atom
                  KMC_AL2=>KMC_AL1%NextNeigh !it is possible that KMC_AL1 points to some other cell after RefreshCell is done when the atom "atom" has moved to a different cell
                  CALL RefreshCell(atom)
                  KMC_AL1=>KMC_AL2
               END DO
            !END IF
         END DO
         
         !NULLIFY(AtomCoord)
         !NULLIFY(AtomIsMoving)
         !CALL Delete(AL) !Not needed
         
         
         !IF (NAtomKMC>6000) THEN
         !   WRITE(*,*) "...bypassing full system relax due to its size"
         !   RETURN
         !END IF
         !CALL OpenFile(UnitRelaxCoord)
         !WRITE(UnitRelaxCoord,FMT='(i5)') NAtomKMC
         !WRITE(UnitRelaxCoord,FMT='(3f20.8)') (KMCCellSize(idx),idx=1,3)
         !atom=>Site
         !DO WHILE (ASSOCIATED(atom))
         !   WRITE(UNIT=UnitRelaxCoord,FMT='(3f20.8,i5)') &
         !      (atom%coord(idx),idx=1,3),atom%iSpecies
         !   atom=>atom%NextSite
         !END DO
         !CLOSE(UnitRelaxCoord)
         
         !!Print .input file
         !CALL OpenFile(UnitRelaxInput)
         !TmpName='Rlx'
         !CALL GenerateTADInput(UnitRelaxInput,TmpName,-1,TempLo,TempLo,1._dp, &
         !     0.01_dp,NAtomKMC,NMoveKMC,-100,.FALSE.,TmpReal,TmpRealDP,TmpRealDP, &
         !     TmpRealDP,TmpRealDP,TmpRealDP,.FALSE.,0_i2b)
         !CLOSE(UnitRelaxInput)
         !CALL TADDriver('r') !relax type run
         !OPEN(UNIT=110,FILE=FileRelaxedPositions)
         !READ(110,*) CrapChar
         !atom=>Site
         !DO i=1,NMoveKMC
         !   READ(110,*) (atom%coord(idx), idx=1,3)
         !   atom=>atom%NextNeigh
         !END DO
         !CLOSE(110)
      END IF
   END SUBROUTINE RelaxKMCSystem
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE InitializeUnitCell
      IMPLICIT NONE

      SELECT CASE (UnitCellType)
      CASE(1); CALL LoadFCC100()
      CASE(2); CALL LoadBCC()
      CASE(3); CALL LoadSC()
      CASE(4); CALL LoadGraphite()
      CASE(5); CALL LoadHCP()
      CASE(6); CALL LoadDiamond()
      END SELECT
   END SUBROUTINE InitializeUnitCell
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LoadFCC100()
      IMPLICIT NONE

      UnitCell=LatticeConst
      a=LatticeConst(1)
      KMCMaxLatticeSpacing=a/SQRT(2.)
      !BOCutoffForCfg=8.3 !Coordination - bulk#
      !BOCutoffForCfg=6.45 !Coordination - surface#
      !BOCutoffForCfg=4.90
      SizeUnitCellAtomCoord=14
      ALLOCATE(UnitCellAtomCoord(SizeUnitCellAtomCoord,3)); NAllottedOthers=NAllottedOthers+&
        SizeUnitCellAtomCoord*3

      UnitCellAtomCoord(1,1)=0.; UnitCellAtomCoord(1,2)=0.; UnitCellAtomCoord(1,3)=0.;
      UnitCellAtomCoord(2,1)=a; UnitCellAtomCoord(2,2)=0.; UnitCellAtomCoord(2,3)=0.; 
      UnitCellAtomCoord(3,1)=0.; UnitCellAtomCoord(3,2)=a; UnitCellAtomCoord(3,3)=0.; 
      UnitCellAtomCoord(4,1)=a; UnitCellAtomCoord(4,2)=a; UnitCellAtomCoord(4,3)=0.;
      UnitCellAtomCoord(5,1)=a/2.; UnitCellAtomCoord(5,2)=a/2.; UnitCellAtomCoord(5,3)=0.;
      
      UnitCellAtomCoord(6,1)=0.; UnitCellAtomCoord(6,2)=0.; UnitCellAtomCoord(6,3)=a; 
      UnitCellAtomCoord(7,1)=a; UnitCellAtomCoord(7,2)=0.; UnitCellAtomCoord(7,3)=a; 
      UnitCellAtomCoord(8,1)=0.; UnitCellAtomCoord(8,2)=a; UnitCellAtomCoord(8,3)=a; 
      UnitCellAtomCoord(9,1)=a; UnitCellAtomCoord(9,2)=a; UnitCellAtomCoord(9,3)=a; 
      UnitCellAtomCoord(10,1)=a/2.; UnitCellAtomCoord(10,2)=a/2.; UnitCellAtomCoord(10,3)=a;

      UnitCellAtomCoord(11,1)=a/2.; UnitCellAtomCoord(11,2)=0.; UnitCellAtomCoord(11,3)=a/2.; 
      UnitCellAtomCoord(12,1)=0; UnitCellAtomCoord(12,2)=a/2.; UnitCellAtomCoord(12,3)=a/2.; 
      UnitCellAtomCoord(13,1)=a/2.; UnitCellAtomCoord(13,2)=a; UnitCellAtomCoord(13,3)=a/2.; 
      UnitCellAtomCoord(14,1)=a; UnitCellAtomCoord(14,2)=a/2.; UnitCellAtomCoord(14,3)=a/2.; 
   END SUBROUTINE LoadFCC100
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LoadBCC()
      IMPLICIT NONE
      INTEGER :: i,j

      UnitCell=LatticeConst
      a=LatticeConst(1)
      KMCMaxLatticeSpacing=SQRT(3.)*a/2. !required to compute BO
      !BOCutoffForCfg=7.
      SizeUnitCellAtomCoord=39
      ALLOCATE(UnitCellAtomCoord(SizeUnitCellAtomCoord,3)); NAllottedOthers=NAllottedOthers+&
        SizeUnitCellAtomCoord*3


      UnitCellAtomCoord(1,1)=0.; UnitCellAtomCoord(1,2)=0.; UnitCellAtomCoord(1,3)=0.
      UnitCellAtomCoord(2,1)=a; UnitCellAtomCoord(2,2)=0.; UnitCellAtomCoord(2,3)=0.
      UnitCellAtomCoord(3,1)=0.; UnitCellAtomCoord(3,2)=a; UnitCellAtomCoord(3,3)=0.
      UnitCellAtomCoord(4,1)=a; UnitCellAtomCoord(4,2)=a; UnitCellAtomCoord(4,3)=0.

      UnitCellAtomCoord(5,1)=0.; UnitCellAtomCoord(5,2)=0.; UnitCellAtomCoord(5,3)=a 
      UnitCellAtomCoord(6,1)=a; UnitCellAtomCoord(6,2)=0.; UnitCellAtomCoord(6,3)=a 
      UnitCellAtomCoord(7,1)=0.; UnitCellAtomCoord(7,2)=a; UnitCellAtomCoord(7,3)=a
      UnitCellAtomCoord(8,1)=a; UnitCellAtomCoord(8,2)=a; UnitCellAtomCoord(8,3)=a

      UnitCellAtomCoord(9,1)=a/2.; UnitCellAtomCoord(9,2)=a/2.; UnitCellAtomCoord(9,3)=a/2. 
      
      i=9
      !FCC points (10-15)
      i=i+1; UnitCellAtomCoord(i,1)=a/2.; UnitCellAtomCoord(i,2)=a/2.; UnitCellAtomCoord(i,3)=0.; !10
      i=i+1; UnitCellAtomCoord(i,1)=a/2.; UnitCellAtomCoord(i,2)=0.; UnitCellAtomCoord(i,3)=a/2.;
      i=i+1; UnitCellAtomCoord(i,1)=0.; UnitCellAtomCoord(i,2)=a/2.; UnitCellAtomCoord(i,3)=a/2.; !12
      i=i+1; UnitCellAtomCoord(i,1)=a/2.; UnitCellAtomCoord(i,2)=a/2.; UnitCellAtomCoord(i,3)=a;
      i=i+1; UnitCellAtomCoord(i,1)=a/2.; UnitCellAtomCoord(i,2)=a; UnitCellAtomCoord(i,3)=a/2.; !14
      i=i+1; UnitCellAtomCoord(i,1)=a; UnitCellAtomCoord(i,2)=a/2.; UnitCellAtomCoord(i,3)=a/2.;
      
      !Midway between BCC corner and center sites (16-23)
      DO j=1,8
         i=i+1
         UnitCellAtomCoord(i,:)=0.333*UnitCellAtomCoord(j,:)+0.667*UnitCellAtomCoord(9,:)
      END DO
      DO j=1,8 !(24-31)
         i=i+1
         UnitCellAtomCoord(i,:)=0.667*UnitCellAtomCoord(j,:)+0.333*UnitCellAtomCoord(9,:)
      END DO
      
      !Sites between two BCC corner, BCC center and FCC edge center sites (32-39)
      i=i+1; UnitCellAtomCoord(i,:)=0.25*(UnitCellAtomCoord(9,:)+UnitCellAtomCoord(1,:)+ &
        UnitCellAtomCoord(2,:)+UnitCellAtomCoord(11,:))
      i=i+1; UnitCellAtomCoord(i,:)=0.25*(UnitCellAtomCoord(9,:)+UnitCellAtomCoord(1,:)+ &
        UnitCellAtomCoord(3,:)+UnitCellAtomCoord(12,:))
      i=i+1; UnitCellAtomCoord(i,:)=0.25*(UnitCellAtomCoord(9,:)+UnitCellAtomCoord(4,:)+ &
        UnitCellAtomCoord(2,:)+UnitCellAtomCoord(15,:))
      i=i+1; UnitCellAtomCoord(i,:)=0.25*(UnitCellAtomCoord(9,:)+UnitCellAtomCoord(4,:)+ &
        UnitCellAtomCoord(3,:)+UnitCellAtomCoord(14,:))
        
      i=i+1; UnitCellAtomCoord(i,:)=0.25*(UnitCellAtomCoord(9,:)+UnitCellAtomCoord(5,:)+ &
        UnitCellAtomCoord(6,:)+UnitCellAtomCoord(11,:))
      i=i+1; UnitCellAtomCoord(i,:)=0.25*(UnitCellAtomCoord(9,:)+UnitCellAtomCoord(5,:)+ &
        UnitCellAtomCoord(7,:)+UnitCellAtomCoord(12,:))
      i=i+1; UnitCellAtomCoord(i,:)=0.25*(UnitCellAtomCoord(9,:)+UnitCellAtomCoord(8,:)+ &
        UnitCellAtomCoord(6,:)+UnitCellAtomCoord(15,:))
      i=i+1; UnitCellAtomCoord(i,:)=0.25*(UnitCellAtomCoord(9,:)+UnitCellAtomCoord(8,:)+ &
        UnitCellAtomCoord(7,:)+UnitCellAtomCoord(14,:))
        
      !DO i=1,39
      !   WRITE(UNIT=UnitScrap,FMT='(3f20.8,i5)') UnitCellAtomCoord(i,:),1
      !END DO

   END SUBROUTINE LoadBCC
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LoadSC()
      IMPLICIT NONE
   END SUBROUTINE LoadSC
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LoadGraphite()
      IMPLICIT NONE
   END SUBROUTINE LoadGraphite
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LoadHCP()
      IMPLICIT NONE
   END SUBROUTINE LoadHCP
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LoadDiamond()
      IMPLICIT NONE
   END SUBROUTINE LoadDiamond
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE KMCRelax

