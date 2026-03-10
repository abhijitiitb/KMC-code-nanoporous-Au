MODULE KMCNEBProcessSearch
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur

   USE KMC_VARIABLE_TYPE
   USE VARIABLE_TYPE
   USE db_manipulate
   USE utilities, ONLY : INT2CHAR,NINT2CHAR
   USE KMCUtilities, ONLY : GetPBCCoord,GetPBCSpacing,PrintAtoms
   USE KMCDbManipulate, ONLY : RefreshCell
   USE OptimizationAL
   USE NEB_package
   
   IMPLICIT NONE
   REAL, DIMENSION(3*MaxNAtoms) :: PrcInitialCoord=0.,PrcSaddleCoord=0.,PrcFinalCoord=0.
   LOGICAL, DIMENSION(MaxNAtoms) :: ProcessAtom
   TYPE(SystemContainer), POINTER :: NEB_AL1,NEB_AL2
   LOGICAL :: NEBInitialized=.FALSE.
   SAVE

   CONTAINS

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddKMCProcessInitialState(CenterCoord,ifile)
   !similar to RelaxKMCSystem 
   !sets up the Atom List (initial state)
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCAtomList), POINTER :: KMC_AL1,KMC_AL2
      INTEGER :: idx,ix,iy,iz,CellIndex(3),NAtoms,NMove(3)
      INTEGER :: i,NoImages(3),MinDistanceIndx,idy
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      REAL :: NewCoord(3),MinDistance,Distance,TmpRealDP=0.,vec(3)
      REAL :: TmpReal=0.,coord(3),MinBox(3),MaxBox(3),CenterCoord(3)
      CHARACTER(len=100) :: TmpName
      CHARACTER :: CrapChar
      LOGICAL :: CellIsMoving
      INTEGER, OPTIONAL :: ifile
      CHARACTER(len=NINT2CHAR) :: cfile

      IF (.NOT. NEBInitialized) THEN !initialize NEB_AL
         CALL MakeSize(NEB_AL1,MaxNAtoms)
         CALL MakeSize(NEB_AL2,MaxNAtoms)
         CALL PotentialInitialize(NEB_AL1)
         CALL PotentialInitialize(NEB_AL2)
         NEBInitialized=.TRUE.
      END IF

      NAtoms=0
      NMove=0
      MaxBox=0.
      MinBox=10000.
      AtomCoord=>NEB_AL1%AtomCoord
      AtomIsMoving=>NEB_AL1%AtomIsMoving
      AtomSpecies=>NEB_AL1%AtomSpecies
      
      WRITE(6,*) "Number of affected cells:",NAffectedCells
      
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
            
            coord=GetPBCSpacing(atom%Coord,CenterCoord)
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
         
      NEB_AL1%NAtoms=NAtoms
      NEB_AL1%NMove=NMove
      NEB_AL1%BoxSize(1:3)=REAL(MaxBox+10.,dp)
      NEB_AL1%BoxSize(4:6)=REAL(MinBox-10.,dp)
      TmpName="KMCNEBInitial.xyz"
      IF (PRESENT(ifile)) THEN
         cfile=INT2CHAR(ifile)
         TmpName="KMCNEBInitial."//cfile(1:LEN_TRIM(cfile))//".xyz"
      END IF
      CALL WriteXYZ(NEB_AL1,TmpName)
      !WRITE(*,*) "3.Size atomcoord:",SIZE(NEB_AL1%AtomCoord)
      !WRITE(*,*) "4.Size atomcoord:",SIZE(NEB_AL2%AtomCoord)
   END SUBROUTINE AddKMCProcessInitialState
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddKMCProcessFinalState(CenterCoord,ifile)
   !similar to RelaxKMCSystem 
   !sets up the Atom List (final state)
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCAtomList), POINTER :: KMC_AL1,KMC_AL2
      INTEGER :: idx,ix,iy,iz,CellIndex(3),NAtoms,NMove(3)
      INTEGER :: i,NoImages(3),MinDistanceIndx,idy
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      REAL :: NewCoord(3),MinDistance,Distance,TmpRealDP=0.,vec(3)
      REAL :: TmpReal=0.,coord(3),MinBox(3),MaxBox(3),CenterCoord(3)
      CHARACTER(len=100) :: TmpName
      CHARACTER :: CrapChar
      LOGICAL :: CellIsMoving
      INTEGER, OPTIONAL :: ifile
      CHARACTER(len=NINT2CHAR) :: cfile

      NAtoms=0
      NMove=0
      MaxBox=0.
      MinBox=10000.
      AtomCoord=>NEB_AL2%AtomCoord
      AtomIsMoving=>NEB_AL2%AtomIsMoving
      AtomSpecies=>NEB_AL2%AtomSpecies
      !WRITE(*,*) "5.Size atomcoord:",SIZE(NEB_AL1%AtomCoord)
      !WRITE(*,*) "6.Size atomcoord:",SIZE(NEB_AL2%AtomCoord)
      DO i=1,NAffectedCells
         CellIndex=AffectedCells(3*i-2:3*i)
         ix=CellIndex(1)
         iy=CellIndex(2)
         iz=CellIndex(3)
         CellIsMoving=IsCellAffected(ix,iy,iz)==1 .OR. IsCellAffected(ix,iy,iz)==2
         KMC_AL1=>Cell(ix,iy,iz)%AL
         DO WHILE (ASSOCIATED(KMC_AL1))
            atom=>KMC_AL1%Atom
            NAtoms=NAtoms+1
            
            coord=GetPBCSpacing(atom%Coord,CenterCoord)
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
         
      NEB_AL2%NAtoms=NAtoms
      NEB_AL2%NMove=NMove
      NEB_AL2%BoxSize(1:3)=REAL(MaxBox+10.,dp)
      NEB_AL2%BoxSize(4:6)=REAL(MinBox-10.,dp)
      TmpName="KMCNEBFinal.xyz"
      IF (PRESENT(ifile)) THEN
         cfile=INT2CHAR(ifile)
         TmpName="KMCNEBFinal."//cfile(1:LEN_TRIM(cfile))//".xyz"
      END IF
      CALL WriteXYZ(NEB_AL2,TmpName)
      !STOP   
      !WRITE(*,*) "7.Size atomcoord:",SIZE(NEB_AL1%AtomCoord)
      !WRITE(*,*) "8.Size atomcoord:",SIZE(NEB_AL2%AtomCoord)
   END SUBROUTINE AddKMCProcessFinalState
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ProcessNEB(MinDisplacement,ActivationEnergy,MaxDisplacement,ifile,ReuseVL)
   !uses NEB calculation to find the activation barrier
   !finally provides the initial, saddle and final coordinates
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL2,AL3
      REAL :: MinDisplacement,MinDisplacement1 !MinDisplacement is the minimum displacement allowed for any process atom
      REAL :: MaxDisplacement !MaxDisplacement is the maximum displacement for any atom
      INTEGER :: ntransitionstates,nimages,NAtoms,image,iatom,nprcatom
      INTEGER :: errorstatus
      REAL :: ActivationEnergy
      REAL(dp) :: springconstant,TSEnergy,distancev(3),distance,maxdistance(MaxNAtoms)
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomCoord1
      TYPE(ChOSContainer), POINTER :: chos
      INTEGER, OPTIONAL :: ifile
      CHARACTER(len=NINT2CHAR) :: cfile
      CHARACTER(len=100) :: TmpName
      LOGICAL :: PrintResults,ReuseVL

      IF (.NOT. NEBInitialized) THEN !initialize NEB_AL
         CALL MakeSize(NEB_AL1,MaxNAtoms)
         CALL MakeSize(NEB_AL2,MaxNAtoms)
         CALL PotentialInitialize(NEB_AL1)
         CALL PotentialInitialize(NEB_AL2)
         NEBInitialized=.TRUE.
      END IF
      
      MinDisplacement1=MinDisplacement
      PrintResults=.FALSE.
      IF (PRESENT(ifile)) THEN
         cfile=INT2CHAR(ifile)
         TmpName="NEBDetails."//cfile(1:LEN_TRIM(cfile))
         WRITE(6,*) "TmpName:",TmpName
         OPEN(UNIT=353,FILE=TRIM(TmpName))
         READ(353,*) MinDisplacement1
         CLOSE(353)
         !create NEB_AL1
         TmpName="KMCNEBInitial."//cfile(1:LEN_TRIM(cfile))//".xyz"
         CALL ReadCoord(NEB_AL1,TmpName)
      
         !create NEB_AL2
         TmpName="KMCNEBFinal."//cfile(1:LEN_TRIM(cfile))//".xyz"
         CALL ReadCoord(NEB_AL2,TmpName)
         PrintResults=.TRUE.
      END IF
      
      springconstant=1._dp
      NULLIFY(AL2)
      NULLIFY(chos)
      nimages=9
      NAtoms=NEB_AL1%NAtoms
      
      TmpName="chos2.xyz"
      IF (PrintResults) TmpName="chos2."//cfile(1:LEN_TRIM(cfile))//".xyz"
      CALL NEB(NEB_AL1,NEB_AL2,nimages=nimages,imode=1101,omode=13,interpmode=1, &
         springconstant=springconstant,ITMAX=100,TSEnergy=TSEnergy,ntransitionstates=ntransitionstates, &
         TS=AL2,NEBChos=chos,ReuseVL=ReuseVL,ChOSFile=TmpName,errorstatus=errorstatus)

      IF (ntransitionstates/=1 .OR. TSEnergy>6._dp) THEN
         WRITE(6,*) "Number of transition states found:",ntransitionstates
         WRITE(6,*) "TSEnergy found:",TSEnergy
         WRITE(6,*) "Increasing the number of images to get a better estimate"
         CALL Delete(AL2)
         CALL DeleteChOS(chos,.TRUE.)
         NULLIFY(AL2)
         NULLIFY(chos)
         nimages=15
         CALL NEB(NEB_AL1,NEB_AL2,nimages=nimages,imode=1101,omode=13,interpmode=1, &
         springconstant=springconstant,ITMAX=100,TSEnergy=TSEnergy,ntransitionstates=ntransitionstates, &
         TS=AL2,NEBChos=chos,ReuseVL=ReuseVL,ChOSFile=TmpName,errorstatus=errorstatus)
         WRITE(UnitScrap,*) "After repeat the energy found to be ..."
         WRITE(UnitScrap,*) "Number of transition states found:",ntransitionstates
         WRITE(UnitScrap,*) "TSEnergy:",TSEnergy
      END IF
      
!  WRITE(6,*) "Completed step 1 of ifile",ifile ; CALL FLUSH(6)
 
      IF (ntransitionstates==0 .OR. TSEnergy<0.05 .OR. .NOT. ASSOCIATED(AL2)) THEN
!such a case can occur when there is no transition state present and the system can fall from one
!state to another. for e.g., if a surface adatom come close to an edge then it has to stick to it
!the process is not first order kinetic in such a case. Such cases can be tackled by performing a
!molecular dynamics before each KMC step. The MD fluctuations will cause the system to fall into 
!the correct state
         WRITE(UnitScrap,*) "Ignoring the process due to ntransitionstates, TSEnergy<0.05", &
            ntransitionstates,TSEnergy
         IF (ASSOCIATED(AL2)) CALL Delete(AL2)
         IF (ASSOCIATED(chos)) CALL Delete(chos)
         ActivationEnergy=0.
         RETURN
      END IF

      IF (ntransitionstates>1 .OR. TSEnergy>6.0_dp) THEN
         WRITE(6,*) "Err>> Even after repeating NEB something is fishy"
         STOP
      END IF
      
!   WRITE(6,*) "Completed step 2 of ifile",ifile ; CALL FLUSH(6)    
   
      !now store the relative coordinates
      AL3=>chos%AL
      AtomCoord=>NEB_AL1%AtomCoord
      maxdistance=0._dp
      DO image=2,nimages
         AL3=>AL3%NextNeigh
         AtomCoord1=>AL3%AtomCoord
         DO iatom=1,NAtoms
            distancev=PBCdistance3(NEB_AL1,AtomCoord1(3*iatom-2:3*iatom),AtomCoord(3*iatom-2:3*iatom))!PBCdistance3 required updation
            distance=SQRT(distancev(1)*distancev(1)+distancev(2)*distancev(2)+distancev(3)*distancev(3))
            maxdistance(iatom)=MAX(maxdistance(iatom),distance)
         END DO
      END DO
!    WRITE(6,*) "Completed step 3 of ifile",ifile ; CALL FLUSH(6)   
      ActivationEnergy=REAL(TSEnergy)
      
      ProcessAtom(1:NAtoms)=.FALSE.
      MaxDisplacement=0.
      nprcatom=0
      DO iatom=1,NAtoms
         distance=maxdistance(iatom)
         MaxDisplacement=MAX(MaxDisplacement,distance)
         IF (distance>REAL(MinDisplacement1,dp)) THEN !it is a moving atom
            nprcatom=nprcatom+1
            ProcessAtom(iatom)=.TRUE.
            PrcInitialCoord(3*iatom-2:3*iatom)=NEB_AL1%AtomCoord(3*iatom-2:3*iatom) !relative coord
            PrcSaddleCoord(3*iatom-2:3*iatom)=AL2%AtomCoord(3*iatom-2:3*iatom)
            PrcFinalCoord(3*iatom-2:3*iatom)=NEB_AL2%AtomCoord(3*iatom-2:3*iatom)
         END IF
      END DO
!    WRITE(6,*) "Completed step 4 of ifile",ifile ; CALL FLUSH(6)   
      IF (PrintResults) THEN
         TmpName="NEBResults."//cfile(1:LEN_TRIM(cfile))
         OPEN(UNIT=353,FILE=TmpName(1:LEN_TRIM(TmpName)))
         WRITE(353,*) NAtoms,nprcatom,ActivationEnergy
         DO iatom=1,NAtoms
            IF (ProcessAtom(iatom)) THEN
               WRITE(353,*) iatom,NEB_AL1%AtomCoord(3*iatom-2:3*iatom), & 
                  AL2%AtomCoord(3*iatom-2:3*iatom),NEB_AL2%AtomCoord(3*iatom-2:3*iatom)
            END IF
         END DO
         !WRITE(353,*) MaxDisplacement
         CLOSE(353)
      END IF
!    WRITE(6,*) "Completed step 5 of ifile",ifile ; CALL FLUSH(6)         
      TmpName="choslocal.xyz"
      IF (PrintResults) TmpName="choslocal."//cfile(1:LEN_TRIM(cfile))//".xyz"
      
      OPEN(UNIT=353,FILE=TmpName(1:LEN_TRIM(TmpName)))
      WRITE(353,*) nprcatom
      WRITE(353,*) " "
      DO iatom=1,NAtoms
         IF (ProcessAtom(iatom)) THEN
            WRITE(353,*) "Ag ",NEB_AL1%AtomCoord(3*iatom-2:3*iatom),maxdistance(iatom)
         END IF
      END DO
      WRITE(353,*) nprcatom
      WRITE(353,*) " "
      DO iatom=1,NAtoms
         IF (ProcessAtom(iatom)) THEN
            WRITE(353,*) "Ag ",AL2%AtomCoord(3*iatom-2:3*iatom)
         END IF
      END DO
      WRITE(353,*) nprcatom
      WRITE(353,*) " "
      DO iatom=1,NAtoms
         IF (ProcessAtom(iatom)) THEN
            WRITE(353,*) "Ag ",NEB_AL2%AtomCoord(3*iatom-2:3*iatom)
         END IF
      END DO
      CLOSE(353)
      
!    WRITE(6,*) "Completed step 6 of ifile",ifile ; CALL FLUSH(6)   
    
      CALL Delete(AL2)
      CALL DeleteChOS(chos,.TRUE.)
      !WRITE(*,*) "11.Size atomcoord:",SIZE(NEB_AL1%AtomCoord)
      !WRITE(*,*) "12.Size atomcoord:",SIZE(NEB_AL2%AtomCoord)
     WRITE(6,*) "Completed step 7 of ifile",ifile ; CALL FLUSH(6)
   END SUBROUTINE ProcessNEB
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE KMCNEBProcessSearch

