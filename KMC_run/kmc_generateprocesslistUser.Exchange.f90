!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateExchangeProcess(atom1,displacement1,atom2,displacement2,istage1)
   !procedure for diffusion on a (100) surface of a FCC
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom1,atom2
      LOGICAL :: ProcessIsPossible
      REAL, DIMENSION(3) :: displacement1,displacement2
      LOGICAL :: NewProcessFound
      TYPE(KMCProcess), POINTER :: prc,MatchingProcess
      INTEGER :: i,istage1
      
      IF (istage1==2) THEN
         WRITE(6,*) "Err>> Stage 2 attempted in CreateHopProcess"
         STOP
      END IF
      
      ProcessIsPossible= CheckDisplacement(atom2,displacement2,tol=1.) !check if this displacement is possible at all
      !since the exchange displacement is not known to us a priori we shall give a larger search tolerance
      IF (ProcessIsPossible) ProcessIsPossible=ExcProcessAbsent(atom1,displacement1,atom2,displacement2)
      
      IF (ProcessIsPossible) THEN
         !Forward process
         WRITE(6,*) "___________________________________________"
         CALL ResetAffectedCellInfo()
         CALL MarkAffectedCells(atom1,displacement1,NumberNeighborCellsInEnv)
         CALL MarkAffectedCells(atom2,displacement2,NumberNeighborCellsInEnv)
         
         NULLIFY(prc)
         IF (istage1==0 .OR. istage1==3) THEN
            ALLOCATE(prc); NAllottedKMCProcess=NAllottedKMCProcess+1
         END IF
         
         !IF (atom1%Cfg%ConfigurationIndex==1) THEN
         !   WRITE(6,*) "Exchange process added with atom2 coord:",atom2%Coord
         !   WRITE(6,*) "Displacement1:",displacement1
         !   WRITE(6,*) "Displacement2:",displacement2
         !END IF
         
         IF (EnabledNEBWithUserKMCProcess) THEN
            CALL NEBExchangeProcess(prc,atom1,displacement1,atom2,displacement2,istage1)
            IF ((istage1==0 .OR. istage1==3)) THEN
               IF (prc%FwdBarrier==0.) THEN !delete such a process
                  CALL DeleteProcess(prc)
                  RETURN
               END IF
            END IF
         ELSE
            CALL ExchangeProcess(prc,atom1,displacement1,atom2,displacement2)
         END IF
         
         IF (istage1==0 .OR. istage1==3) THEN
            NewProcessFound= .NOT. CompareProcess2PL(prc,MatchingProcess,MinDisplacement=0.1)
         
            IF (NewProcessFound) THEN
               CALL AddProcessToPL(prc)
               WRITE(6,*) "Exchange rate:",prc%Rate
               CALL AssignProcessFull(prc)
            ELSE
               IF (prc%NumberProcessAtoms>=MatchingProcess%NumberProcessAtoms) THEN
                  CALL DeleteProcess(prc)
               ELSE!store the process with the fewer number of process atoms
               !the code assumes that both processes are the same but some of the atoms did not show up because 
               !they were moving very less. it will hold on to process atoms that are known to move a lot while 
               !delete those who move less. So the process with fewer process atoms is stored
               !REPLACE MatchingProcess WITH prc
                  CALL ReplaceProcess(MatchingProcess,prc) !MatchingProcess is now floating
                  CALL DeleteProcess(prc)
               END IF
            END IF
         END IF
         WRITE(6,*) "Done."
         
      !ELSE
      !   IF (atom1%Cfg%ConfigurationIndex==1) THEN
      !      WRITE(6,*) "Exchange process NOT added with atom2 coord:",atom2%Coord
      !      WRITE(6,*) "Displacement1:",displacement1
      !      WRITE(6,*) "Displacement2:",displacement2
      !   END IF
      END IF
   
   END SUBROUTINE CreateExchangeProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ExcProcessAbsent(atom1,displacement1,atom2,displacement2)
   !determines whether the process is already present with the atom by checking
   !the KMCProcessSubscriptionAction - the following details are checked
   !the final displacement of two atoms is significant and is given approx by displacement1
   !and displacement2, that all other atoms have small final displacements (less than 1 Ang)
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom1,atom2
      REAL, DIMENSION(3) :: displacement1,displacement2,coord,coord1
      TYPE(KMCProcessSubscriptionInfo), POINTER :: prclist
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: prcatom
      LOGICAL :: ExcProcessAbsent,Found
      INTEGER :: cfgno1,cfgno2,NPrcAtomsMoving,NPrcAtomsMatch,CurrentKMCStep
      REAL :: rmag2
      
      Found=.FALSE. !so far a process matching the exchange process is not found
      
      CurrentKMCStep=KMCStep+KMCBlock*NKMCStepsPerBlock
      
      cfgno1=atom1%Cfg%ConfigurationIndex
      cfgno2=atom2%Cfg%ConfigurationIndex
      prclist=>atom1%PrcSubscriberInfo
      DO WHILE (ASSOCIATED(prclist))
         prc=>prclist%Process
         
         IF (prc%FirstAdded<CurrentKMCStep .AND. prclist%IsActive) THEN 
            NPrcAtomsMoving=0 !# atoms found moving
            NPrcAtomsMatch=0 !#atoms that are moving and their displacement match
            prcatom=>prc%ReactantAL
            DO WHILE (ASSOCIATED(prcatom))
               coord=prcatom%RelFinalPos-prcatom%RelInitialPos !displacement
               rmag2=coord(1)*coord(1)+coord(2)*coord(2)+coord(3)*coord(3) !magnitude displacement
               IF (rmag2>1.) THEN !atom is moving in this process
                  NPrcAtomsMoving=NPrcAtomsMoving+1
                  coord1=coord-displacement1
                  rmag2=coord1(1)*coord1(1)+coord1(2)*coord1(2)+coord1(3)*coord1(3) !deviation from displacement
                  IF (rmag2<0.25 .AND. cfgno1==prcatom%CfgIndex) & !displacement is matching
                     NPrcAtomsMatch=NPrcAtomsMatch+1
                  coord1=coord-displacement2
                  rmag2=coord1(1)*coord1(1)+coord1(2)*coord1(2)+coord1(3)*coord1(3) !deviation from displacement
                  IF (rmag2<0.25 .AND. cfgno2==prcatom%CfgIndex) & !displacement is matching
                     NPrcAtomsMatch=NPrcAtomsMatch+1
               END IF
               prcatom=>prcatom%NextNeigh
            END DO
            Found=NPrcAtomsMoving==2 .AND. NPrcAtomsMatch==2
            
            IF (Found) EXIT !all atoms matched with the process
         END IF
         prclist=>prclist%NextNeigh
      END DO
      ExcProcessAbsent=.NOT. Found
!      WRITE(UNIT=6,FMT='("Exc process absent: ",L1," Press key to continue ")') &
!         ExcProcessAbsent
   END FUNCTION ExcProcessAbsent
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ExchangeProcess(prc,atom1,displacement1,atom2,displacement2)
   !Adds the atoms positions in the initial state
      IMPLICIT NONE
      REAL(dp), PARAMETER :: J0=0.243,Jb=0.2,Js=-0.05
      TYPE(KMCAtom), POINTER :: atom1,atom2 !these are the inputs
      TYPE(KMCAtom), POINTER :: atom3,atom4
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(KMCProcess), POINTER :: prc
      INTEGER :: i,cfgindex,CellIndex(3),AtomIndex1,AtomIndex2
      INTEGER :: BOInitial,BOSaddle
      REAL, DIMENSION(3) :: displacement1,displacement2,initialcoord,saddlecoord,finalcoord,atomcoord1,atomcoord2
      REAL, DIMENSION(3) :: coord,unitcoord,unitcoordnormal,v1,v2,d1,d2
      REAL :: a,tol,Eact,s1,s2,s3
      
      AtomIndex1=atom1%Index
      AtomIndex2=atom2%Index
      atomcoord1=atom1%Coord
      prc%CfgType=>atom1%Cfg
      
      !STEP 1:
      BOInitial=0
      a=LatticeConst(1)
      tol=RadialTolerance
      
      coord=GetPBCCoord(atomcoord1+a*(/.5,.5,0./))
      atom4=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom4)) BOInitial=BOInitial+1
      
      coord=GetPBCCoord(atomcoord1+a*(/.5,-.5,0./))
      atom4=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom4)) BOInitial=BOInitial+1
   
      coord=GetPBCCoord(atomcoord1+a*(/-.5,.5,0./))
      atom4=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom4)) BOInitial=BOInitial+1
      
      coord=GetPBCCoord(atomcoord1+a*(/-.5,-.5,0./))
      atom4=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom4)) BOInitial=BOInitial+1
      
      IF (BOInitial==4) THEN
         WRITE(6,*) "BOInitial is four adatoms are present -- exchange cannot take place -- bug in the code"
         STOP
      END IF
      
      !STEP 2:
      !atom1%coord=GetPBCCoord(atom1%coord+0.5*displacement1) !displace to saddle point
      !atom2%coord=GetPBCCoord(atom2%coord+0.5*displacement2)
      !CALL RefreshCell(atom)
      
      !STEP 3:
      BOSaddle=0
      saddlecoord=atom1%coord+(/displacement1(1),displacement1(2),0./)
      unitcoord=2.*(/displacement1(1),displacement1(2),0./)/a
      unitcoordnormal(1)= unitcoord(2)
      unitcoordnormal(2)=-unitcoord(1)
      unitcoordnormal(3)=0.
      
      coord=GetPBCCoord(saddlecoord+a*0.5*unitcoordnormal)
      atom3=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom3)) BOSaddle=BOSaddle+1
      
      coord=GetPBCCoord(saddlecoord-a*0.5*unitcoordnormal)
      atom3=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom3)) BOSaddle=BOSaddle+1
      
      !STEP 4:
      !atom%coord=GetPBCCoord(atom%coord-0.5*displacement) !back to initial spot
      !CALL RefreshCell(atom)
      
      !STEP 5:
      Eact=J0+Jb*REAL(BOInitial,dp)+Js*REAL(BOSaddle,dp)
      WRITE(6,*) "BO:",BOInitial,BOSaddle
      prc%Frequency=1.e13
      prc%FwdBarrier=Eact
      prc%Rate=prc%Frequency*EXP(-Eact/kBT)
      
      !STEP 6:
      v1=(/0.,0.,0./)
      v2=v1+displacement1+displacement2
      s3=a/SQRT(2.)+RadialTolerance
      
      DO i=1,NAffectedCells
         CellIndex=AffectedCells(3*i-2:3*i)
         AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
         DO WHILE (ASSOCIATED(AL1))
            initialcoord=GetPBCSpacing(AL1%Atom%Coord,atomcoord1) !relative position
            d1=initialcoord-v1
            d2=initialcoord-v2
            s1=SQRT(d1(1)*d1(1)+d1(2)*d1(2)+d1(3)*d1(3))
            s2=SQRT(d2(1)*d2(1)+d2(2)*d2(2)+d2(3)*d2(3))
            IF (s1<=s3 .OR. s2<=s3) THEN
               IF (AtomIndex1==AL1%Atom%Index) THEN
                  finalcoord=initialcoord+displacement1
               ELSEIF (AtomIndex2==AL1%Atom%Index) THEN
                  finalcoord=initialcoord+displacement2
               ELSE
                  finalcoord=initialcoord
               END IF
               saddlecoord=0.5*(initialcoord+finalcoord)
               cfgindex=AL1%atom%Cfg%ConfigurationIndex
               CALL AddAtomToProcess(prc,initialcoord,finalcoord,saddlecoord,cfgindex,AL1%Atom%Species)
            END IF
            AL1=>AL1%NextNeigh
         END DO
      END DO
   END SUBROUTINE ExchangeProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE NEBExchangeProcess(prc,atom1,displacement1,atom2,displacement2,istage1)
   !This subroutine will generate a process that contains participant atoms which relax by a distance of MaxRelaxDist or more
   !atom (input) is the main atom which is going to hop
   !prc (output) is the process which shall be filled in
   !displacement(input) is the displacement of the atom
      IMPLICIT NONE
      INTEGER :: NAffectedCells1,CellIndex(3),i,ix,iy,iz,NAtoms,cfgindex,NProcessAtoms
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCAtom), POINTER :: atom1,atom2,atom3
      TYPE(KMCAtomList), POINTER :: AL1
      REAL :: coord(3),displacement1(3),displacement2(3),initialcoord(3),saddlecoord(3),finalcoord(3),s1,s2
      REAL :: Eact,tol,MaxDisplacement
      LOGICAL :: CellIsMoving
      INTEGER :: istage1,iatom,nprcatom
      CHARACTER(len=NINT2CHAR) :: cfile
      CHARACTER(len=30) :: TmpName
      
      CHARACTER(len=2) :: CrapChar2
      TYPE(KMCAtom), POINTER :: atom10
      TYPE(KMCAtomList), POINTER :: KMC_AL1
      REAL :: relcoordstg1(3),relcoordstg3(3)
      INTEGER :: NAtomsStg1,NAtomsStg3
      
      WRITE(6,*) "NEBExcProcess istage and ifile:",istage1,ifile+1
      
      IF (istage1==0 .OR. istage1==3) prc%CfgType=>atom1%Cfg
      tol=RadialTolerance
      
      !the cells affected by the main atom and its resulting displacement are already known
      !add few layers
      NAffectedCells1=NAffectedCells
      DO i=1,NAffectedCells1
         CellIndex=AffectedCells(3*i-2:3*i)
         CALL ExtendAffectedCellLayer(CellIndex,value=1,range=1) !add a layer of moving atoms
      END DO
      NAffectedCells1=NAffectedCells
      DO i=1,NAffectedCells1
         CellIndex=AffectedCells(3*i-2:3*i)
         CALL ExtendAffectedCellLayer(CellIndex,value=2,range=1) !add a layer of pseudo-moving atoms
      END DO
      NAffectedCells1=NAffectedCells
      DO i=1,NAffectedCells1
         CellIndex=AffectedCells(3*i-2:3*i)
         CALL ExtendAffectedCellLayer(CellIndex,value=3,range=1) !add a layer of non-moving atoms
      END DO
      NAffectedCells1=NAffectedCells
      DO i=1,NAffectedCells1
         CellIndex=AffectedCells(3*i-2:3*i)
         CALL ExtendAffectedCellLayer(CellIndex,value=3,range=1) !add a layer of non-moving atoms
      END DO
      
      ifile=ifile+1
         
      coord=atom1%Coord
      WRITE(6,*) "Center coord:",coord
      
      IF (istage1==0) THEN
         CALL AddKMCProcessInitialState(coord)
         atom1%Coord=atom1%Coord+displacement1
         atom2%Coord=atom2%Coord+displacement2
         CALL AddKMCProcessFinalState(coord)
         atom1%Coord=atom1%Coord-displacement1
         atom2%Coord=atom2%Coord-displacement2
      END IF
      IF (istage1==1) THEN
         CALL AddKMCProcessInitialState(coord,ifile)
         atom1%Coord=atom1%Coord+displacement1
         atom2%Coord=atom2%Coord+displacement2
         CALL AddKMCProcessFinalState(coord,ifile)
         atom1%Coord=atom1%Coord-displacement1
         atom2%Coord=atom2%Coord-displacement2
         cfile=INT2CHAR(ifile)
         TmpName="NEBDetails."//cfile(1:LEN_TRIM(cfile))
         OPEN(UNIT=353,FILE=TmpName(1:LEN_TRIM(TmpName)))
         WRITE(353,*) 0.15 !MinDisplacement
         CLOSE(353)
      END IF
      
      
      !Do NEB calculation
      IF (istage1==0) THEN
         CALL ProcessNEB(MinDisplacement=0.15,ActivationEnergy=Eact, &
            MaxDisplacement=MaxDisplacement,ReuseVL=.FALSE.)
      ELSEIF (istage1==3) THEN !read the NEB results for ifile (ifile should be global in module decl.)
         cfile=INT2CHAR(ifile)
         TmpName="NEBResults."//TRIM(cfile)
         WRITE(6,*) "Reading NEB results from file:",TRIM(TmpName)
         OPEN(UNIT=353,FILE=TRIM(TmpName))
         READ(353,*) NAtoms,nprcatom,Eact
         WRITE(6,*) "Total number of atoms:",NAtoms
         WRITE(6,*) "Number of process atoms:",nprcatom
         WRITE(6,*) "Recovered activation barrier is:",Eact
         
         !*************************************************************
         !check whether the coordinates are correct
         NAtomsStg1=0
         NAtomsStg3=0
         OPEN(UNIT=391,FILE="KMCNEBInitial."//cfile(1:LEN_TRIM(cfile))//".xyz")
         READ(391,*) NAtomsStg1
         READ(391,*) CrapChar2
         
         DO i=1,NAffectedCells
            CellIndex=AffectedCells(3*i-2:3*i)
            KMC_AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
            DO WHILE (ASSOCIATED(KMC_AL1))
               atom10=>KMC_AL1%Atom
               NAtomsStg3=NAtomsStg3+1
               relcoordstg3=GetPBCSpacing(atom10%Coord,coord)
               IF (NAtomsStg3>NAtomsStg1) THEN
                  WRITE(6,*) "Number of atoms in affected cells much more than that in"// &
                    "KMCNEBInitial."//cfile(1:LEN_TRIM(cfile))//".xyz"
                  STOP
               END IF
               READ(391,*) CrapChar2,relcoordstg1
               IF (ANY(ABS(relcoordstg1-relcoordstg3)>0.0001)) THEN
                  WRITE(UnitScrap,*) "Mismatch in coordinates",ifile," @ position ",NAtomsStg3
                  WRITE(UnitScrap,*) "Expected rel coord:",relcoordstg1
                  WRITE(UnitScrap,*) "Found rel coord:",relcoordstg3
                  WRITE(UnitScrap,*) "Main atom coord:",coord
                  STOP
               END IF
               KMC_AL1=>KMC_AL1%NextNeigh
            END DO
         END DO
         IF (NAtomsStg1/=NAtomsStg3) THEN
            WRITE(UnitScrap,*) "Number of atoms in stg 1 and 3 are different",ifile,NAtomsStg1,NAtomsStg3
                  CALL FLUSH(UnitScrap)
            STOP
         END IF
         !*************************************************************
         
         ProcessAtom(1:NAtoms)=.FALSE.
         DO i=1,nprcatom
            READ(353,*) iatom,initialcoord,saddlecoord,finalcoord
            ProcessAtom(iatom)=.TRUE.
            PrcInitialCoord(3*iatom-2:3*iatom)=initialcoord
            PrcSaddleCoord(3*iatom-2:3*iatom)=saddlecoord
            PrcFinalCoord(3*iatom-2:3*iatom)=finalcoord
         END DO
         !READ(353,*) MaxDisplacement
         CLOSE(353)
         
         !delete files
         CALL SYSTEM("rm NEBResults."//TRIM(cfile)// & 
            " NEBDetails."//TRIM(cfile)// &
            " choslocal."//TRIM(cfile)//".xyz"// &
            " KMCNEBInitial."//TRIM(cfile)//".xyz"// &
            " KMCNEBFinal."//TRIM(cfile)//".xyz"// &
            " chos2."//TRIM(cfile)//".xyz")
      END IF
      
      IF (istage1==0 .OR. istage1==3) THEN
         prc%Frequency=1.e13
         prc%FwdBarrier=Eact
         prc%Rate=prc%Frequency*EXP(-Eact/kBT)
         !prc%MaxPrcAtomDisplacement=MaxDisplacement
      
         !now obtain the process atoms
         NAtoms=0
         NProcessAtoms=0
         DO i=1,NAffectedCells
            CellIndex=AffectedCells(3*i-2:3*i)
            ix=CellIndex(1)
            iy=CellIndex(2)
            iz=CellIndex(3)
            AL1=>Cell(ix,iy,iz)%AL
            DO WHILE (ASSOCIATED(AL1))
               atom3=>AL1%Atom
               NAtoms=NAtoms+1
               IF (ProcessAtom(NAtoms)) THEN
                  NProcessAtoms=NProcessAtoms+1
                  initialcoord=PrcInitialCoord(3*NAtoms-2:3*NAtoms)
                  saddlecoord=PrcSaddleCoord(3*NAtoms-2:3*NAtoms)
                  finalcoord=PrcFinalCoord(3*NAtoms-2:3*NAtoms)
                  cfgindex=atom3%Cfg%ConfigurationIndex
                  CALL AddAtomToProcess(prc,initialcoord,finalcoord,saddlecoord,cfgindex,atom3%Species)
               END IF
               AL1=>AL1%NextNeigh
            END DO
         END DO
         WRITE(6,*) "Created temporary process with number of atoms:",NProcessAtoms
      END IF
   END SUBROUTINE NEBExchangeProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
