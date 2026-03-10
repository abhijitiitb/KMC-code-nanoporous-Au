!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
!Potential issues:
! a) IMPROPER OPTIMIZATION- 
! Due to improper optimization in the localized optimization code (this improper optimization is partly
! because of the way the outer non-moving shell is chosen. Since this shell will keep varying from iteration 
! to iteration, there can be a signficant deviation from the global optimized structure) 
! Solution: For the time being in NEBHopProcess atom3%Coord is stored and not the optimized position from the
! NEB initial state (see subroutine)
! b) DUPLICATE PROCESS - A cfg has access to two processes which are both the same. Two processes are different
! if the participating atoms are slightly off from each other. The subroutine CompareProcess2Process
! would correctly find that the two processes are different. However, the local configurations of the
! current state are some where near the average of the coordinates of the two processes' atoms. Hence,
! the cfg incorrectly would get both the processes
! Solution; CompareProcess2Process now uses 2*RadialTolerance based on propagation of errors
! c) MISTAKEN IDENTITY - In MD, two processes were found at different times and these processes where the same. 
! However, upon doing NEB it was found that they have slightly different participating atoms. This resulted in the
! same process being reported twice. Thus the less displaced atoms should get lower weightage
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateHopProcess(atom,displacement,istage1)
   !Last checked: Jan 12, 2011
   !procedure for diffusion on a (100) surface of a FCC
   !single atom
   !istage1 can 0,1 or 3 but not 2
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      LOGICAL :: ProcessIsPossible
      REAL :: displacement(3)
      LOGICAL :: NewProcessFound
      TYPE(KMCProcess), POINTER :: prc,MatchingProcess
      INTEGER :: i,istage1
      
      IF (istage1==2) THEN
         WRITE(6,*) "Err>> Stage 2 attempted in CreateHopProcess"
         STOP
      END IF
      
      ProcessIsPossible=CheckDisplacement(atom,displacement,tol=1.0) !check if this displacement is possible at all
      IF (ProcessIsPossible) ProcessIsPossible=HopProcessAbsent(atom,displacement)

      IF (ProcessIsPossible) THEN
         WRITE(6,*) "___________________________________________"
         CALL ResetAffectedCellInfo()
         CALL MarkAffectedCells(atom,displacement,NumberNeighborCellsInEnv)

         NULLIFY(prc)
         IF (istage1==0 .OR. istage1==3) THEN
            ALLOCATE(prc); NAllottedKMCProcess=NAllottedKMCProcess+1
         END IF
         
         IF (EnabledNEBWithUserKMCProcess) THEN
            CALL NEBHopProcess(atom,prc,displacement,istage1)
            IF ((istage1==0 .OR. istage1==3)) THEN
               IF (prc%FwdBarrier==0.) THEN !delete such a process
                  CALL DeleteProcess(prc)
                  RETURN
               END IF
            END IF
         ELSE
            CALL HopProcess(atom,prc,displacement)
         END IF
         
         IF (istage1==0 .OR. istage1==3) THEN
            NewProcessFound= .NOT. CompareProcess2PL(prc,MatchingProcess,MinDisplacement=0.1)
            IF (NewProcessFound) THEN

               CALL AddProcessToPL(prc)
               WRITE(6,*) "Hop rate:",prc%Rate
               CALL AssignProcessFull(prc)
               WRITE(6,*) "Process is added @ KMC step:",KMCStep+KMCBlock*NKMCStepsPerBlock
            ELSE
               IF (prc%NumberProcessAtoms>=MatchingProcess%NumberProcessAtoms) THEN
                  WRITE(6,*) "Process is already present ... code was unable to assign the existing process"
                  WRITE(UnitScrap,*) "Process is already present ... yet it was not allotted earlier"
                  WRITE(UnitScrap,*) "Matching process index:",MatchingProcess%Index
                  WRITE(UnitScrap,*) "Matching process added first to db @ iteration:",MatchingProcess%FirstAdded
! IF (MatchingProcess%FirstAdded<KMCStep+NKMCStepsPerBlock*KMCBlock) THEN
!    WRITE(UnitScrap,*) "Culprit atom is:",atom%Index
!    CALL PrintAtomProcessAssignments()
! END IF
                  CALL PrintProcess(prc)
!CALL CheckKMC()
!WRITE(6,*) "No issues with KMC database found"
                  CALL DeleteProcess(prc)
               ELSE !store the process with the fewer number of process atoms
               !the code assumes that both processes are the same but some of the atoms did not show up because 
               !they were moving very less. it will hold on to process atoms that are known to move a lot while 
               !delete those who move less. So the process with fewer process atoms is stored
               !REPLACE MatchingProcess WITH prc
          WRITE(6,*) "About to replace process"
                  CALL ReplaceProcess(MatchingProcess,prc)
                  CALL DeleteProcess(prc) !prc now contains the information of MatchingProcess
               END IF
            END IF
         END IF
         WRITE(6,*) "Done."
      END IF
   END SUBROUTINE CreateHopProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION HopProcessAbsent(atom,displacement)
   !determines whether the process is already present with the atom by checking
   !the KMCProcessSubscriptionAction - the following details are checked
   !the final displacement of atom atom is significant and is given approx by displacement
   !all other atoms have small final displacements (less than 1 Ang)
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      REAL, DIMENSION(3) :: displacement,coord
      TYPE(KMCProcessSubscriptionInfo), POINTER :: prclist
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: prcatom
      LOGICAL :: HopProcessAbsent,Found
      INTEGER :: cfgno,NPrcAtomsMoving,NPrcAtomsMatch,CurrentKMCStep
      REAL :: rmag2

      Found=.FALSE. !so far a process matching the hop process is not found
      
      CurrentKMCStep=KMCStep+KMCBlock*NKMCStepsPerBlock
!WRITE(UnitScrap,*) "Checking if proposed process already exists @ KMC step:",KMCStep+KMCBlock*NKMCStepsPerBlock
      
      cfgno=atom%Cfg%ConfigurationIndex
      prclist=>atom%PrcSubscriberInfo
!WRITE(UnitScrap,*) "Investigating hop for atom",atom%Index," displacement ",displacement
      DO WHILE (ASSOCIATED(prclist))
         prc=>prclist%Process
    
         IF (prc%FirstAdded<CurrentKMCStep .AND. prclist%IsActive) THEN 
!the process is old and we should check the displacement
!if this process is new it will be added only at the current KMC step
!if CurrentKMCStep==0 then the process was read from previous files
            NPrcAtomsMoving=0 !# atoms found moving
            NPrcAtomsMatch=0 !#atoms that are moving and their displacement match
            prcatom=>prc%ReactantAL
            DO WHILE (ASSOCIATED(prcatom)) !find such a process in subscribed processes
               coord=prcatom%RelFinalPos-prcatom%RelInitialPos !displacement
               rmag2=coord(1)*coord(1)+coord(2)*coord(2)+coord(3)*coord(3) !magnitude displacement
               IF (rmag2>1.) THEN !atom is moving in this process
                  NPrcAtomsMoving=NPrcAtomsMoving+1
                  coord=coord-displacement
                  rmag2=coord(1)*coord(1)+coord(2)*coord(2)+coord(3)*coord(3) !deviation from displacement
                  IF (rmag2<Rad2Correction .AND. cfgno==prcatom%CfgIndex) & !displacement is matching
                     NPrcAtomsMatch=NPrcAtomsMatch+1
                  END IF
               prcatom=>prcatom%NextNeigh
            END DO
            Found=NPrcAtomsMoving==1 .AND. NPrcAtomsMatch==1
            
!WRITE(UnitScrap,*) "HopProcessAbsent: Comparing to process # ",prc%Index," found ",Found,prc%FirstAdded
            
            IF (Found) EXIT !all atoms matched with the process
         END IF
         prclist=>prclist%NextNeigh
      END DO
      HopProcessAbsent=.NOT. Found
!WRITE(UNIT=UnitScrap,FMT='("Hop process absent: ",L1," Press key to continue ")') &
 !  HopProcessAbsent
   END FUNCTION HopProcessAbsent
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE HopProcess(atom,prc,displacement)
   !Adds the atoms positions in the initial state
      IMPLICIT NONE
      REAL(dp), PARAMETER :: J0=0.243_dp,Jb=0.2_dp,Js=-0.05_dp
      TYPE(KMCAtom), POINTER :: atom,atom1
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(KMCProcess), POINTER :: prc
      INTEGER :: i,cfgindex,CellIndex(3),AtomIndex
      INTEGER :: BOInitial,BOSaddle
      REAL, DIMENSION(3) :: displacement,initialcoord,saddlecoord,finalcoord,atomcoord
      REAL, DIMENSION(3) :: coord,unitcoord,unitcoordnormal,v1,v2,d1,d2
      REAL :: a,tol,Eact,s1,s2,s3
      
      AtomIndex=atom%Index
      atomcoord=atom%Coord
      prc%CfgType=>atom%Cfg
      
      !STEP 1:
      BOInitial=0
      a=LatticeConst(1)
      tol=RadialTolerance
      
      coord=GetPBCCoord(atomcoord+0.5*a*(/1.,1.,0./))
      atom1=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom1)) BOInitial=BOInitial+1
      
      coord=GetPBCCoord(atomcoord+0.5*a*(/1.,-1.,0./))
      atom1=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom1)) BOInitial=BOInitial+1
      
      coord=GetPBCCoord(atomcoord+0.5*a*(/-1.,1.,0./))
      atom1=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom1)) BOInitial=BOInitial+1
      
      coord=GetPBCCoord(atomcoord+0.5*a*(/-1.,-1.,0./))
      atom1=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom1)) BOInitial=BOInitial+1
      
      IF (BOInitial==4) THEN
         WRITE(6,*) "BOInitial is four -- hop cannot occur -- bug in the code"
         STOP
      END IF
      
      !STEP 2:
      atom%coord=GetPBCCoord(atom%coord+0.5*displacement) !displace to saddle point
      CALL RefreshCell(atom)
      
      !STEP 3:
      BOSaddle=0
      saddlecoord=atom%coord
      unitcoord=SQRT(2.)*displacement/a
      unitcoordnormal(1)= unitcoord(2)
      unitcoordnormal(2)=-unitcoord(1)
      unitcoordnormal(3)=0.
      
      coord=GetPBCCoord(saddlecoord+a*(0.5*unitcoord+unitcoordnormal)/SQRT(2.))
      atom1=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom1)) BOSaddle=BOSaddle+1
      
      coord=GetPBCCoord(saddlecoord+a*(0.5*unitcoord-unitcoordnormal)/SQRT(2.))
      atom1=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom1)) BOSaddle=BOSaddle+1
      
      coord=GetPBCCoord(saddlecoord+a*(-0.5*unitcoord+unitcoordnormal)/SQRT(2.))
      atom1=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom1)) BOSaddle=BOSaddle+1
      
      coord=GetPBCCoord(saddlecoord+a*(-0.5*unitcoord-unitcoordnormal)/SQRT(2.))
      atom1=>SearchAtom(coord,tol)
      IF (ASSOCIATED(atom1)) BOSaddle=BOSaddle+1
      
      !STEP 4:
      atom%coord=GetPBCCoord(atom%coord-0.5*displacement) !back to initial spot
      CALL RefreshCell(atom)
      
      !STEP 5:
      Eact=J0+Jb*REAL(BOInitial,dp)+Js*REAL(BOSaddle,dp)
      WRITE(6,*) "BO:",BOInitial,BOSaddle
      prc%Frequency=1.e13_dp
      prc%FwdBarrier=Eact
      prc%Rate=prc%Frequency*EXP(-Eact/kBT)
      
      !STEP 6:
      v1=(/0.,0.,0./)
      v2=v1+displacement
      s3=a/SQRT(2.)+RadialTolerance
      
      DO i=1,NAffectedCells
         CellIndex=AffectedCells(3*i-2:3*i)
         AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
         DO WHILE (ASSOCIATED(AL1))
            initialcoord=GetPBCSpacing(AL1%Atom%Coord,atomcoord) !relative position
            d1=initialcoord-v1
            d2=initialcoord-v2
            s1=SQRT(d1(1)*d1(1)+d1(2)*d1(2)+d1(3)*d1(3))
            s2=SQRT(d2(1)*d2(1)+d2(2)*d2(2)+d2(3)*d2(3))
            IF (s1<=s3 .OR. s2<=s3) THEN
               IF (AtomIndex==AL1%Atom%Index) THEN
                  finalcoord=initialcoord+displacement
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
   END SUBROUTINE HopProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE NEBHopProcess(atom,prc,displacement,istage1)
   !This subroutine will generate a process that contains participant atoms which relax by a distance of MaxRelaxDist or more
   !atom (input) is the main atom which is going to hop
   !prc (output) is the process which shall be filled in
   !displacement(input) is the displacement of the atom
      IMPLICIT NONE
      INTEGER :: NAffectedCells1,CellIndex(3),i,ix,iy,iz,NAtoms,cfgindex,NProcessAtoms
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCAtom), POINTER :: atom,atom3,atom1
      TYPE(KMCAtomList), POINTER :: AL1
      REAL :: coord(3),displacement(3),initialcoord(3),saddlecoord(3),finalcoord(3),s1,s2,coord1(3)
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
      
!WRITE(UnitScrap,*) "Hop process for atom index:",atom%Index,atom%Coord,istage1,ifile+1
      WRITE(6,*) "NEBHopProcess istage and ifile:",istage1,ifile+1
! IF (ifile>169) CALL PrintAtomProcessAssignments()
      IF (istage1==0 .OR. istage1==3) prc%CfgType=>atom%Cfg
      tol=RadialTolerance
      
      !the cells affected by the main atom and its resulting displacement are already known
      !add few layers
      NAffectedCells1=NAffectedCells
      DO i=1,NAffectedCells1
         CellIndex=AffectedCells(3*i-2:3*i)
         CALL ExtendAffectedCellLayer(CellIndex,value=1,range=1) !add a layer of psuedo-moving atoms
      END DO
      NAffectedCells1=NAffectedCells
      DO i=1,NAffectedCells1
         CellIndex=AffectedCells(3*i-2:3*i)
         CALL ExtendAffectedCellLayer(CellIndex,value=2,range=1) !add a layer of non-moving atoms
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
         
      coord=atom%Coord
      WRITE(6,*) "Center coord:",coord
      
      IF (istage1==0) THEN
         CALL AddKMCProcessInitialState(coord)
         atom%Coord=atom%Coord+displacement
         CALL AddKMCProcessFinalState(coord)
         atom%Coord=atom%Coord-displacement
      END IF
      IF (istage1==1) THEN
         CALL AddKMCProcessInitialState(coord,ifile)
         atom%Coord=atom%Coord+displacement
         CALL AddKMCProcessFinalState(coord,ifile)
         atom%Coord=atom%Coord-displacement
         cfile=INT2CHAR(ifile)
         TmpName="NEBDetails."//cfile(1:LEN_TRIM(cfile))
         OPEN(UNIT=353,FILE=TmpName(1:LEN_TRIM(TmpName)))
         WRITE(353,*) 0.1 !MinDisplacement
         CLOSE(353)
      END IF
      
      !Do NEB calculation
      IF (istage1==0) THEN
         CALL ProcessNEB(MinDisplacement=0.1,ActivationEnergy=Eact,MaxDisplacement=MaxDisplacement, &
            ReuseVL=.FALSE.)
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
                  WRITE(UnitScrap,*) "For reference check ","KMCNEBInitial."//cfile(1:LEN_TRIM(cfile))//"-1.xyz"
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
            CellIsMoving=IsCellAffected(ix,iy,iz)==1
            AL1=>Cell(ix,iy,iz)%AL
            DO WHILE (ASSOCIATED(AL1))
               atom3=>AL1%Atom
               NAtoms=NAtoms+1
               IF (ProcessAtom(NAtoms)) THEN
                  NProcessAtoms=NProcessAtoms+1
      !WARNING: Using initialcoord from the PrcInitialCoord could be potentially dangerous.
      !The NEB calculation optimizes the initial state. During this process the atoms could be displaced
      !substantially resulting in a process initial coord that is quite different from what
      !was originally expected!!
      !initialcoord=PrcInitialCoord(3*NAtoms-2:3*NAtoms) !if during the optimization at the starting of the NEB calculation, 
      !the atoms move too much then initial coordinates will appear quite different from actual local environments
                  !coord1=coord+initialcoord
                  !coord1=GetPBCCoord(coord1)
                  !atom1=>SearchAtom(coord1,tol)
                  !IF (.NOT. ASSOCIATED(atom1)) THEN
                  !   WRITE(6,*) "Atom has relaxed too much during NEB:"
                  !   WRITE(6,*) "Optimized atom coord:",coord1
                  !   WRITE(6,*) "Expected atom coord:",atom3%Coord
                  !   STOP
                  !END IF
                  !IF (atom1%Cfg%ConfigurationIndex/=atom3%Cfg%ConfigurationIndex) THEN
                  !   WRITE(6,*) "Atom has moved too much during NEB and the "
                  !   WRITE(6,*) "atom which is expected to be present has a different cfg index"
                  !   WRITE(6,*) "Found cfg index of atom:",atom1%Cfg%ConfigurationIndex
                  !   WRITE(6,*) "Expected cfg index of atom:",atom3%Cfg%ConfigurationIndex
                  !   STOP
                  !END IF
                  initialcoord=GetPBCSpacing(atom3%Coord,coord)
                  saddlecoord=PrcSaddleCoord(3*NAtoms-2:3*NAtoms)
                  finalcoord=PrcFinalCoord(3*NAtoms-2:3*NAtoms)
                  cfgindex=atom3%Cfg%ConfigurationIndex
                  CALL AddAtomToProcess(prc,initialcoord,finalcoord,saddlecoord,cfgindex,atom3%Species)
               END IF
               AL1=>AL1%NextNeigh
            END DO
         END DO
         WRITE(6,*) "Created temporary process with number of atoms:",NProcessAtoms
      
         IF (NProcessAtoms>40) THEN
            WRITE(6,*) "Number of process atoms is too large"
            STOP
         END IF
      END IF
   END SUBROUTINE NEBHopProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
