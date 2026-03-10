MODULE KMCUpdateSystemDomainDecomposition
!contains update subroutines for LEKMC that create an affected region using domains 
!this version is newer than kmc_algorithm_updatesystem.carve.f90

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
      TYPE(KMCProcess), POINTER :: prc 
      REAL :: coord(3)
      INTEGER :: i,j,CellIndex1(3),ix,iy,iz,indxmax,indxmin,NAffectedCellsC4,NAffectedCellsC3,range
      INTEGER :: ixmax,ixmin,iymax,iymin,izmax,izmin
      INTEGER :: NewC4AffectedCells(3*SizeAffectedCells),NewC4AffectedCellsValue(SizeAffectedCells)
      
      coord=atom%Coord
      WRITE(6,*) "Process selected index:",prc%Index
      CALL FindDomain(atom,domainsize=1) !keeping the process atom as center, curve out the domain and the buffers C3 and C4
      CALL FindBufferC3(atom,domainsize=1,C3thickness=2)
      NAffectedCellsC3=NAffectedCells
      CALL RemoveProcessesInAffected()
      !AffectedCells1(1:3*NAffectedCellsC4)=AffectedCells(1:3*NAffectedCellsC4) !this is used to search process
      !Find which cells lie between NAffectedCellsC3+1 and NAffectedCellsC4
      CALL UpdateSystemCfg(prc,atom) 
      
      !DO i=NAffectedCellsC3+1,NAffectedCellsC4
         !CellIndex1=AffectedCells(3*i-2:3*i)
         !CALL ExtendAffectedCellLayer(CellIndex1,value=1,range=1) !this is needed to capture any new process that can appear since our local configurations have changed
      !END DO
   
      CALL AssignProcessInAllAffectedCells()   
      CALL CheckForNewProcesses() !find if new processes can occur in affected region -- and assign them
      WRITE(6,*) "Successfully updated system after selecting process"
      !CALL PrintAtomProcessAssignments()

   END SUBROUTINE UpdateSystem
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE FindDomain(atom,domainsize)
   !takes main atom as center atom and carves out a domain based on the range
   !range is being fixed in this subroutine -- this is an internal parameter
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      INTEGER :: i,j,ix,iy,iz,CellIndex(3),CellIndex1(3),domainsize
      REAL ::coord(3)
      
      CALL ResetAffectedCellInfo()
      coord=atom%Coord
      WRITE(*,*)"The affected atom coordinate is", coord
      CellIndex=GetCellLocation(coord)
      write(*,*)CellIndex
      
      !AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL !main center atom
      domainsize=1 !domainsize can be altered hence domain size will change
      DO ix=-domainsize,domainsize
         DO iy=-domainsize,domainsize
            DO iz=-domainsize,domainsize
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               !CALL PrintAffectedCellsDD(CellIndex1)
               CALL AddCellToAffectedList(CellIndex1,value=1)
            ENDDO
         ENDDO
      ENDDO
      !OPEN(unit=491,FILE="affectedlistdomain") 
      !CALL PrintAffectedCells(491) !list of all the domain cells 
   END SUBROUTINE FindDomain
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx
   SUBROUTINE FindBufferC3(atom,domainsize,C3thickness)
   !if domain 1 is done it carves out layers based on the range 
   !the value of range can be different from the one used in FindDomain
      IMPLICIT NONE
      INTEGER :: range,ix,iy,iz,CellIndex1(3),CellIndex(3),NAffectedCellsC3,domainsize,C3thickness
      REAL :: coord(3)
      TYPE(KMCAtom), POINTER :: atom
      
      coord=atom%Coord
      CellIndex=GetCellLocation(coord)
      
      range=domainsize+C3thickness !range of C3
      !Creating Buffer C3
      DO ix=-range,-domainsize-1
         DO iy=-range,range
            DO iz=-range,range
                  CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
                  !CALL PrintAffectedCellsBuffer(CellIndex1)
                  CALL AddCellToAffectedList(CellIndex1,value=2)
            ENDDO
         ENDDO
      ENDDO
      DO ix=domainsize+1,range
         DO iy=-range,range
            DO iz=-range,range
                  CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
                  !write(*,*)CellIndex1
                  !CALL PrintAffectedCellsBuffer(CellIndex1)
                  CALL AddCellToAffectedList(CellIndex1,value=2)
            ENDDO
         ENDDO
      ENDDO
      DO ix=-domainsize,domainsize      !Domain Size=2; Changes according to domain size line 158
         DO iy=-range,-domainsize-1
            DO iz=-range,range
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               !CALL PrintAffectedCellsBuffer(CellIndex1)
               CALL AddCellToAffectedList(CellIndex1,value=2)
            ENDDO
         ENDDO
      ENDDO
      DO ix=-domainsize,domainsize
         DO iy=domainsize+1,range
            DO iz=-range,range
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               !CALL PrintAffectedCellsBuffer(CellIndex1)
               CALL AddCellToAffectedList(CellIndex1,value=2)
            ENDDO
         ENDDO
      ENDDO
      DO ix=-domainsize,domainsize
         DO iy=-domainsize,domainsize
            DO iz=-range,-domainsize-1
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               !CALL PrintAffectedCellsBuffer(CellIndex1)
               CALL AddCellToAffectedList(CellIndex1,value=2)
            ENDDO
         ENDDO
      ENDDO
      DO ix=-domainsize,domainsize
         DO iy=-domainsize,domainsize
            DO iz=domainsize+1,range
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               !CALL PrintAffectedCellsBuffer(CellIndex1)
               CALL AddCellToAffectedList(CellIndex1,value=2)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(35)
      NAffectedCellsC3=NAffectedCells
      WRITE(*,*)"Affected Cells in domain C3 is  : ", NAffectedCellsC3
      !OPEN(unit=493,FILE="affectedlistbufferC3")
      !CALL PrintAffectedCells(493)
      !CLOSE(493)
      !CALL UpdateSystemCfg(prc,atom)
    ENDSUBROUTINE FindBufferC3
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
   SUBROUTINE FindBufferC4(atom,domainsize,C3thickness,C4thickness)
   !adds a region called C4 which comprises of non-moving buffer atom once C3 region has been created
   !C4 is required only for optimization purposes
   !value of range can be set independently of what has been used for C1 and C3 regions
      IMPLICIT NONE
      INTEGER :: range,ix,iy,iz,CellIndex(3),CellIndex1(3),NAffectedCellsC4,domainsize,C3thickness,C4thickness
      REAL :: coord(3)
      TYPE(KMCAtom), POINTER :: atom
      coord=atom%Coord
      CellIndex=GetCellLocation(coord)
      
      !Creating Buffer C4
      range=domainsize+C3thickness+C4thickness  
      
      DO ix=-range,-domainsize-C3thickness-1
         DO iy=-range,range
            DO iz=-range,range
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               !CALL PrintAffectedCellsBuffer1(CellIndex1)
               CALL AddCellToAffectedList(CellIndex1,value=3)
            ENDDO
         ENDDO
      ENDDO
      DO ix=domainsize+C3thickness+1,range
         DO iy=-range,range
            DO iz=-range,range
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               !CALL PrintAffectedCellsBuffer1(CellIndex1)
               CALL AddCellToAffectedList(CellIndex1,value=3)
            ENDDO
         ENDDO
      ENDDO
      
      DO ix=-domainsize-C3thickness,domainsize+C3thickness      !Domain Size=2; Changes according to domain size line 158
         DO iy=-range,-domainsize-C3thickness-1
            DO iz=-range,range
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               !CALL PrintAffectedCellsBuffer1(CellIndex1)
               CALL AddCellToAffectedList(CellIndex1,value=3)
            ENDDO
         ENDDO
      ENDDO
      DO ix=-domainsize-C3thickness,domainsize+C3thickness
         DO iy=domainsize+C3thickness+1,range
            DO iz=-range,range
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               !CALL PrintAffectedCellsBuffer1(CellIndex1)
               CALL AddCellToAffectedList(CellIndex1,value=3)
            ENDDO
         ENDDO
      ENDDO

      DO ix=-domainsize-C3thickness,domainsize+C3thickness
         DO iy=-domainsize-C3thickness,domainsize+C3thickness
            DO iz=-range,-domainsize-C3thickness-1
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               !CALL PrintAffectedCellsBuffer1(CellIndex1)
               CALL AddCellToAffectedList(CellIndex1,value=3)
            ENDDO
         ENDDO
      ENDDO
      DO ix=-domainsize-C3thickness,domainsize+C3thickness
         DO iy=-domainsize-C3thickness,domainsize+C3thickness
            DO iz=domainsize+C3thickness+1,range
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               !CALL PrintAffectedCellsBuffer1(CellIndex1)
               CALL AddCellToAffectedList(CellIndex1,value=3)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(45)
      NAffectedCellsC4=NAffectedCells
      WRITE(*,*)"Affected Cells in domain C4 is  : ", NAffectedCellsC4
      !OPEN(unit=494,FILE="affectedlistbufferC4")
      !CALL PrintAffectedCells(494)      !Working.. Last checked 6th Feb '12'
      !CLOSE(494)
   ENDSUBROUTINE FindBufferC4
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintAffectedCellsDD(CellIndex)
      IMPLICIT NONE
      !Prints the cell info of the domain C1
      INTEGER :: CellIndex(3)
      TYPE(KMCAtomList), POINTER :: AL1
      AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
      open(unit=25,file="domain1.xyz")
      DO WHILE(ASSOCIATED(AL1))
      write(25,*)"Ag ", AL1%Atom%Coord ,"T",   "T",   "T"
      !write(*,*)AL1%Atom%Index
      AL1=>AL1%NextNeigh
      ENDDO
   ENDSUBROUTINE PrintAffectedCellsDD
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintAffectedCellsBuffer(CellIndex) !for C3
      IMPLICIT NONE
      INTEGER :: CellIndex(3)
      TYPE(KMCAtomList), POINTER :: AL1
      open(unit=35,file="bufferc31.xyz")
      AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
      DO WHILE(ASSOCIATED(AL1))
         write(35,*)"Ag  ",AL1%Atom%Coord,"T","T","T"
         AL1=>AL1%NextNeigh
      ENDDO
   ENDSUBROUTINE PrintAffectedCellsBuffer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintAffectedCellsBuffer1(CellIndex) !Creating C4;can be modified later
      IMPLICIT NONE
      INTEGER :: i,j,k,CellIndex(3)
      TYPE(KMCAtomList), POINTER :: AL1
      open(unit=45,file="bufferc41.xyz")
      AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
      DO WHILE(ASSOCIATED(AL1))
         write(45,*)"Ag  ",AL1%Atom%Coord,"F","F","F"
         AL1=>AL1%NextNeigh
      ENDDO
   ENDSUBROUTINE PrintAffectedCellsBuffer1
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
   SUBROUTINE UpdateSystemCfg(prc,atom)
   !prc and atom are process and atom, respectively, that is selected
   !See step 3 on line 1.
   !move atoms to new position and update the configurations......will be same for Domain decomposition
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom,atom1
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: prcatom
      REAL, DIMENSION(3) :: coord,coord1
      INTEGER :: range,ix,iy,iz,i
      REAL :: tol
      INTEGER :: CellIndex(3),NAffectedCells1,CellIndex1(3),NAffectedCellsC3
      NAffectedCellsC3=NAffectedCells
      
      coord=atom%Coord
      
      tol=RadialTolerance
      !write(*,*)prc%CfgType%ConfigurationIndex,atom%Cfg%ConfigurationIndex
      IF (prc%CfgType%ConfigurationIndex/=atom%Cfg%ConfigurationIndex) THEN !Test
         WRITE(6,*) "Process cfg does not match atom cfg"
         WRITE(6,*) "Prc Cfg index:",prc%CfgType%ConfigurationIndex
         WRITE(6,*) "Atom Cfg index:",atom%Cfg%ConfigurationIndex
         STOP
      END IF
      !write(*,*)"The Participating Atoms are.."
      prcatom=>prc%ReactantAL !participating atoms
      DO WHILE (ASSOCIATED(prcatom))
         coord1=GetPBCCoord(coord+prcatom%RelInitialpos)
         !write(*,*)coord1
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
      CALL FindBufferC4(atom,domainsize=1,C3thickness=2,C4thickness=1)
         
      CALL RelaxKMCSystem(coord) 
      !CALL RelaxFullKMCSystem()
      NAffectedCells1=NAffectedCellsC3
      DO i=1,NAffectedCells1 !cells C1+C3
         CellIndex=AffectedCells(3*i-2:3*i)
         !CALL PrintAffectedCellsBuffer(CellIndex)
         !CALL ConfirmCellIsAffected(CellIndex)  !check added ------<<<<<<<<<------------
         AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
         DO WHILE (ASSOCIATED(AL1))
            atom1=>AL1%Atom
   
            IF (atom1%IsMoving) THEN             
               CALL DeleteLocalCfg(atom1)
               CALL AddLocalCfg(atom1)
!               IF (atom1%Index==1 .OR. atom1%Index==2)CALL PrintEnv(AL1%Atom)
            END IF
            CALL AssignConfig(atom1)
            AL1=>AL1%NextNeigh
         END DO
      END DO
      write(*,*)"System Updated......" !Running fine till here....last Checked 30 Jan
   
   END SUBROUTINE UpdateSystemCfg
   !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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
            !CALL PerformNEB()
            !WRITE(*,*) "Number of files:",nfiles, " @stage ",istage1
            CALL PerformNEBParallel(nfiles)
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

END MODULE KMCUpdateSystemDomainDecomposition
