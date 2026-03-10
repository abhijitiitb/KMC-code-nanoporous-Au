MODULE KMCUpdateSystemCarve
!contains update subroutines for LEKMC that carve out an affected region by adding cells one-by-one
!The main step is the Updating of system. Following steps are performed
!Step 1. Find cells where the selected process atom are/will be located. This group of cells is C1
!Step 2. Find cells neighboring C1 within environment cutoff where config can change. The union of step 1 and step 2 gives C2, where configs can change. (see SBRTN FindAffectedCells)
!Step 3. Due to change in config of C2 relaxation is performed by adding buffer layer of type 2 (pseudo-moving) and 3 (non-moving). C2 plus buffer layer cells are C3. Later in step 5 the energy relaxation is performed and the environments in C3 are updated. Note that buffer atoms do not move yet their environments can change.
!Step 4. Find processes that involve atoms in C2. The cells involving this process are added to C2 and together they are called C2, C3 and the new cells are called C4. C4 is the region where processes could have been altered. This information needs to be saved some where else since we will use it later.
!Step 5. Coming back to C3 (we ignore the new cells forming C4), energy relaxation is performed and the cfg in C3 are updated
!Step 6. Coming back to C4, processes are removed from C4. And new process assignments are made.
!Step 7. Add a layer to C4 to allow discovery of new previously unseen processes
   USE KMC_VARIABLE_TYPE
   USE KMCDbManipulate
   USE KMCRecords
   USE GenerateProcess
   USE Ecuyer_random
   USE KMCRelax
   IMPLICIT NONE
   INTEGER :: AffectedCells1(3*SizeAffectedCells)=0
   CONTAINS
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE UpdateSystem(prc,atom)
   !works with list of affected cells and updates the following information
   ! 0) Find which cells are affected
   ! 1) Processes associated with atoms lose their activity by one
   ! 2) Atom configurations are recomputed
   ! 3) Existing processes are checked to increase their activity
   ! 4) Check for new processes
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCProcess), POINTER :: prc ,prc1
      REAL :: coord(3)
      INTEGER :: NAffectedCellsC3,NAffectedCellsC4,i,j,CellIndex1(3),ix,iy,iz
      INTEGER :: NewC4AffectedCells(3*SizeAffectedCells),NewC4AffectedCellsValue(SizeAffectedCells)

      coord=atom%Coord
WRITE(6,*) "Process selected index:",prc%Index
      CALL FindAffectedCells(prc,atom) !just find the cell list
      
      CALL AddCellsForEnergyMinimization()
      NAffectedCellsC3=NAffectedCells !cells with affected environments - also this will be used in energy min
      
      CALL ExtendAffectedCellsBasedOnProcess()
      NAffectedCellsC4=NAffectedCells !cells with affected environments and processes

      NAffectedCells=NAffectedCellsC3 !focus on only the first NAffectedCellsC3 cells
      CALL RemoveProcessesInAffected() !remove association between atoms and cfg without modifying affected cells

      AffectedCells1(1:3*NAffectedCellsC4)=AffectedCells(1:3*NAffectedCellsC4) !this is used to search process
      !Find which cells lie between NAffectedCellsC3+1 and NAffectedCellsC4
      j=0
      DO i=NAffectedCellsC3+1,NAffectedCellsC4
         ix=AffectedCells(3*i-2)
         iy=AffectedCells(3*i-1)
         iz=AffectedCells(3*i)
         j=j+1
         NewC4AffectedCells(3*j-2)=ix
         NewC4AffectedCells(3*j-1)=iy
         NewC4AffectedCells(3*j)=iz
         NewC4AffectedCellsValue(j)=IsCellAffected(ix,iy,iz)
         IsCellAffected(ix,iy,iz)=0 !reset value so that relaxation step (step 5) can be performed
      END DO
      CALL UpdateSystemCfg(prc,atom) !update atoms positions, do energy minimization, and find new cfg in affected cells

      !Reset the IsCellAffected values of cells added during relaxation (these correspond to non-moving atoms)
      DO i=NAffectedCellsC3+1,NAffectedCells
         ix=AffectedCells(3*i-2)
         iy=AffectedCells(3*i-1)
         iz=AffectedCells(3*i)
         IsCellAffected(ix,iy,iz)=0
      END DO

      !now recover back the information we had regd C4
      NAffectedCells=NAffectedCellsC4
      AffectedCells(1:3*NAffectedCellsC4)=AffectedCells1(1:3*NAffectedCellsC4)
      j=0
      DO i=NAffectedCellsC3+1,NAffectedCellsC4
         j=j+1
         ix=NewC4AffectedCells(3*j-2)
         iy=NewC4AffectedCells(3*j-1)
         iz=NewC4AffectedCells(3*j)
         AffectedCells(3*i-2)=ix
         AffectedCells(3*i-1)=iy
         AffectedCells(3*i)=iz
         IsCellAffected(ix,iy,iz)=NewC4AffectedCellsValue(j)
      END DO
      
      !configuration indices might have changed in NAffectedCells so delete the 
      !processes involving atoms in these cells, 
      !but first we will find which are the cells affected by the missing processes
      
      !now we remove the processes from the affected cells where cfg could be different
      DO i=NAffectedCellsC3+1,NAffectedCellsC4
         CellIndex1=AffectedCells(3*i-2:3*i)
         CALL ExtendAffectedCellLayer(CellIndex1,value=1,range=1) !this is needed to capture any new process that can appear since our local configurations have changed
      END DO
      
      !finally we assign and search processes
      CALL AssignProcessInAllAffectedCells()
      
      !prc1=>KMC_PL
      !DO WHILE (ASSOCIATED(prc1))
      !   CALL AssignProcessInAffectedCells(prc1) !look at the existing list and assign processes
      !   prc1=>prc1%NextNeigh
      !END DO
      CALL CheckForNewProcesses() !find if new processes can occur in affected region -- and assign them
      !write(*,*) "d"
      !CALL CheckSubscriptionInfoAllAtoms()
      !CALL CheckSubscriberAllProcesses()
WRITE(6,*) "Successfully updated system after selecting process"
!CALL PrintAtomProcessAssignments()

   END SUBROUTINE UpdateSystem
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE FindAffectedCells(prc,atom)
   !adds which cells are affected by a process (see steps 1 and 2 on line 1)
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: prcatom
      REAL, DIMENSION(3) :: coord,coord1
      INTEGER :: CellIndex(3),NAffectedCells1,cutoff,i
      
      !WRITE(*,*) "Atom index(A)",atom%Index,atom%Coord
      !WRITE(*,*) "Process:",prc%Index
      coord=atom%Coord
      prcatom=>prc%ReactantAL
      CALL ResetAffectedCellInfo() !NAffectedCells=0
      
      !Get the cells of the selected process atoms
      DO WHILE (ASSOCIATED(prcatom))
         coord1=GetPBCCoord(coord+prcatom%RelInitialpos)
         CellIndex=GetCellLocation(coord1)
         !write(*,*) coord1,CellIndex
         CALL AddCellToAffectedList(CellIndex,value=1)
         
         coord1=GetPBCCoord(coord+prcatom%RelSaddlepos)
         CellIndex=GetCellLocation(coord1)
         !write(*,*) coord1,CellIndex
         CALL AddCellToAffectedList(CellIndex,value=1)
         
         coord1=GetPBCCoord(coord+prcatom%RelFinalpos)
         CellIndex=GetCellLocation(coord1)
         !write(*,*) coord1,CellIndex
         CALL AddCellToAffectedList(CellIndex,value=1)
         prcatom=>prcatom%NextNeigh
      END DO
      
      !now find cells within one environment cutoff
      NAffectedCells1=NAffectedCells
      cutoff=MAXVAL(CEILING(EnvCutOffRadius/KMCCellSize))
      cutoff=MAX(cutoff,1)
      DO i=1,NAffectedCells1
         CellIndex=AffectedCells(3*i-2:3*i)
         CALL ExtendAffectedCellLayer(CellIndex,value=1,range=cutoff)
      END DO
      
      !OPEN(unit=493,FILE="affectedlist")
      !CALL PrintAffectedCells(493)
   END SUBROUTINE FindAffectedCells
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ExtendAffectedCellsBasedOnProcess()
   !configurations in affected cells will probably change. 
   !so check which processes are in turn affected
   !See step 4 on line 1. We will construct the cells C4
      IMPLICIT NONE
      INTEGER :: NAffectedCells1,i,j,CellIndex(3),CellIndex1(3), CellIndex2(3)
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(KMCAtom), POINTER :: activeatom
      TYPE(KMCProcessSubscriptionInfo), POINTER :: ptr
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: prcatom
      REAL :: coord(3),coord1(3)
      
      NAffectedCells1=NAffectedCells !so far these are cells with affected environments
      DO i=1,NAffectedCells1
         CellIndex=AffectedCells(3*i-2:3*i)
         
         !CALL ExtendAffectedCellLayer(CellIndex,value=1)
!add a tolerance of 1 cell -- this can be important when a coordinate lies on a cell edge 
!and its cell is unclear
         
         AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
         DO WHILE (ASSOCIATED(AL1))
            ptr=>AL1%atom%PrcSubscriberInfo !list of processes this atom is participating in
            DO WHILE (ASSOCIATED(ptr))
               activeatom=>ptr%ActiveAtom !get the head atom
               prc=>ptr%Process
               coord=activeatom%Coord
               prcatom=>prc%ReactantAL
               DO WHILE (ASSOCIATED(prcatom))
                  coord1=GetPBCCoord(coord+prcatom%RelInitialpos)
                  CellIndex1=GetCellLocation(coord1)
                  CALL AddCellToAffectedList(CellIndex1,value=2) !add these as moving buffer atoms
                  !CALL ExtendAffectedCellLayer(CellIndex1,value=1,range=1) !this is needed to capture any new process that can appear since our local configurations have changed
                  prcatom=>prcatom%NextNeigh
               END DO
               ptr=>ptr%NextNeigh
            END DO
            AL1=>AL1%NextNeigh
         END DO
      END DO
      
      !CALL PrintAffectedCells(493)
      !CLOSE(493)
      
   END SUBROUTINE ExtendAffectedCellsBasedOnProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RemoveProcessesInAffected()
      IMPLICIT NONE
      TYPE(KMCAtomList), POINTER :: AL1
      INTEGER :: i,CellIndex(3),NAffectedCells1
      
      NAffectedCells1=NAffectedCells
      DO i=1,NAffectedCells1
         CellIndex=AffectedCells(3*i-2:3*i)
         AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
         DO WHILE (ASSOCIATED(AL1))
            CALL RemoveAllProcessesFromAtom(AL1%Atom)
            AL1=>AL1%NextNeigh
         END DO
      END DO
   END SUBROUTINE RemoveProcessesInAffected
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddCellsForEnergyMinimization()
   !The environments of all these cells could be affected
      IMPLICIT NONE
      INTEGER :: NAffectedCells1,CellIndex(3),cutoff,i
      
      !For hop we have 2 and 3
      ! for exchange we have 1, 2, 3
      NAffectedCells1=NAffectedCells
      cutoff=MAXVAL(CEILING(EnvCutOffRadius/KMCCellSize))
      cutoff=MAX(cutoff,1)
      DO i=1,NAffectedCells1
         CellIndex=AffectedCells(3*i-2:3*i)
         CALL ExtendAffectedCellLayer(CellIndex,value=2,range=2*cutoff) !add a layer of psuedo-moving atoms
      END DO
   END SUBROUTINE AddCellsForEnergyMinimization
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE UpdateSystemCfg(prc,atom)
   !prc and atom are process and atom, respectively, that is selected
   !See step 3 on line 1.
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom,atom1
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: prcatom
      REAL, DIMENSION(3) :: coord,coord1
      REAL :: tol
      INTEGER :: i,CellIndex(3),NAffectedCells1,oldcfg,NAffectedCellsC3
      INTEGER :: ix,iy,iz,cellin(3),cutoff
 
      !move atoms to new position and update the configurations
      coord=atom%Coord
      tol=RadialTolerance
     
      IF (prc%CfgType%ConfigurationIndex/=atom%Cfg%ConfigurationIndex) THEN
         WRITE(6,*) "Process cfg does not match atom cfg"
         WRITE(6,*) "Prc Cfg index:",prc%CfgType%ConfigurationIndex
         WRITE(6,*) "Atom Cfg index:",atom%Cfg%ConfigurationIndex
         STOP
      END IF
      
      prcatom=>prc%ReactantAL
      DO WHILE (ASSOCIATED(prcatom))
         coord1=GetPBCCoord(coord+prcatom%RelInitialpos)
         atom1=>SearchAtom(coord1,tol)
         
         IF (.NOT. ASSOCIATED(atom1)) THEN
            WRITE(6,*) "Atom1 in UpdateSystemCfg is not associated"
            STOP
         END IF
         atom1%Coord=GetPBCCoord(coord+prcatom%RelFinalpos)
         CALL RefreshCell(atom1)
         prcatom=>prcatom%NextNeigh
      END DO
      
      !CALL CheckOverlap() !--added to make sure there is no overlap between atoms
      
      !NAffectedCells1=NAffectedCells
      !DO i=1,NAffectedCells1
      !   CellIndex=AffectedCells(3*i-2:3*i)
      !   CALL ExtendAffectedCellLayer(CellIndex,value=2) !add a layer of non-moving atoms
      !END DO
      
      NAffectedCellsC3=NAffectedCells !we will reset NAffectedCells to this value
      !note that the pseudo layer atoms are not moving but their cfg could be modified
 
      NAffectedCells1=NAffectedCells
      DO i=1,NAffectedCells1
         CellIndex=AffectedCells(3*i-2:3*i)
         CALL ExtendAffectedCellLayer(CellIndex,value=3,range=2) !add a layer of non-moving atoms
      END DO
  
      !Check atoms in correct cells
      !DO ix=1,NCells(1)
      !   DO iy=1,NCells(2)
      !      DO iz=1,NCells(3)
      !         AL1=>Cell(ix,iy,iz)%AL
      !         DO WHILE (ASSOCIATED(AL1))
      !            atom1=>AL1%Atom
      !            cellin=GetCellLocation(atom1%Coord)
      !            IF (ANY((/ix,iy,iz/)/=cellin)) THEN
      !               WRITE(*,*) "Cells do not match"
      !               STOP
      !            END IF
      !            AL1=>AL1%NextNeigh
      !         END DO
      !      END DO
      !   END DO
      !END DO
      !atom1=>KMC_AL
      !DO WHILE (ASSOCIATED(atom1))
      !   cellin=GetCellLocation(atom1%Coord)
      !   IF (ANY(cellin/=atom1%CellIndex)) THEN
      !      WRITE(6,*) "Atom in wrong cell"
      !      STOP
      !   END IF
      !   atom1=>atom1%NextNeigh
      !END DO
      !STOP
      
      CALL RelaxKMCSystem(atom%Coord)
      !CALL RelaxFullKMCSystem()
      
      NAffectedCells=NAffectedCellsC3 !THE NON-MOVING BUFFER ATOMS WILL NOW BE OVERWRITTEN. reset the value since we dont needs cells belonging to non-moving buffer -- they will not be able to see the change in environments. Also NAffectedCellsC3 will be used to find processes that are affected.

      DO i=1,NAffectedCells !in all cells where the relaxation was performed
         CellIndex=AffectedCells(3*i-2:3*i)
         !CALL ConfirmCellIsAffected(CellIndex)  !check added ------<<<<<<<<<------------
         AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
         DO WHILE (ASSOCIATED(AL1))
            atom1=>AL1%Atom
         !IF (ASSOCIATED(atom1%PrcSubscriberInfo)) THEN
         !   WRITE(*,*) "atom is still associated with process---"
         !   STOP
         !END IF
         
            !oldcfg=atom1%Cfg%ConfigurationIndex
            IF (atom1%IsMoving) THEN
!               IF (atom1%Index==1 .OR. atom1%Index==2)CALL PrintEnv(AL1%Atom)
!oldcfg=atom1%Cfg%ConfigurationIndex
               CALL DeleteLocalCfg(atom1)
               CALL AddLocalCfg(atom1)
!               IF (atom1%Index==1 .OR. atom1%Index==2)CALL PrintEnv(AL1%Atom)
            END IF
            CALL AssignConfig(atom1)
!WRITE(*,*) "Modifying cfg of atom ",atom1%Index," @ KMC step",KMCStep+KMCBlock*NKMCStepsPerBlock,oldcfg, &
!atom1%Cfg%ConfigurationIndex
            
!      IF (NumberCfgTypes>100000.AND. NumberProcesses>100000) THEN
!         CALL GenSystemRecord()
!         WRITE(6,*) "NumberCfgTypes or NumberProcesses greater than allowed number"
!         STOP
!      END IF
            !IF (oldcfg/=atom1%Cfg%ConfigurationIndex) &
               !write(*,*) atom1%Index,oldcfg,atom1%Cfg%ConfigurationIndex,atom1%Coord
            !IF (atom1%IsMoving)write(*,*) atom1%Index,atom1%Coord,oldcfg,atom1%Cfg%ConfigurationIndex
            !IF (atom1%Index==427)CALL PrintEnv(AL1%Atom)
            AL1=>AL1%NextNeigh
         END DO
      END DO
      
      
      !atom1=>KMC_AL
      !DO WHILE (ASSOCIATED(atom1))
      !   IF (atom1%IsMoving) THEN
      !      IF (.NOT. CompareAtomEnvToCfg(atom1,atom1%Cfg)) THEN
      !         WRITE(6,*) "Atom config is incorrect"
      !         STOP
      !      END IF
      !   END IF
      !   atom1=>atom1%NextNeigh
      !END DO
      !DO ix=1,NCells(1)
      !   DO iy=1,NCells(2)
      !      DO iz=1,NCells(3)
      !         AL1=>Cell(ix,iy,iz)%AL
      !         DO WHILE (ASSOCIATED(AL1))
      !            atom1=>AL1%Atom
      !            cellin=GetCellLocation(atom1%Coord)
      !            IF (ANY((/ix,iy,iz/)/=cellin)) THEN
      !               WRITE(*,*) "Cells do not match"
      !               STOP
      !            END IF
      !            AL1=>AL1%NextNeigh
      !         END DO
      !      END DO
      !   END DO
      !END DO
      !WRITE(6,*) "Everythings fine"
      !STOP
   END SUBROUTINE UpdateSystemCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckForNewProcesses()
      IMPLICIT NONE
      INTEGER :: i,CellIndex(3),NAffectedCells1,istage1 !,ifile1,nfiles1
      REAL :: ActivationEnergy,MaxDisplacement
      
      !copy this information because NAffectedCells and AffectedCells will be changed
      NAffectedCells1=NAffectedCells
      AffectedCells1(1:3*NAffectedCells)=AffectedCells(1:3*NAffectedCells)
      
      nfiles=0
      !stage1 is prepare the 
      DO istage1=1,3
         ifile=0
         
!         WRITE(UnitScrap,*) "Number of files:",nfiles, " @ begin of stage ",istage1
!         CALL FLUSH(UnitScrap)
         
         IF (istage1==2 .AND. EnabledNEBWithUserKMCProcess) THEN
            !CALL PerformNEB() !for serial
!            WRITE(*,*) "Number of files:",nfiles, " @stage ",istage1
            CALL PerformNEBParallel(nfiles) !for parallel
!   WRITE(UnitScrap,*) "Completed all NEB calculations"
!   CALL FLUSH(UnitScrap)
         ELSE
            nfiles=0
            DO i=1,NAffectedCells1
               istage=istage1 !need to reset to istage1 each iteration
               CellIndex=AffectedCells1(3*i-2:3*i)
               CALL GenerateProcessListCell(CellIndex)
            END DO
            nfiles=ifile
!            WRITE(*,*) "Number of files:",nfiles, " @ end of stage ",istage1
         END IF
      END DO
   END SUBROUTINE CheckForNewProcesses
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE KMCUpdateSystemCarve
