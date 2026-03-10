MODULE UserProcess
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
!Provided by user find the processes possible in the current system
   USE KMC_VARIABLE_TYPE
   USE KMCDbManipulate
   USE KMCConfigurations
   USE KMCRelax
   
   USE KMCProcessSearchParallel !for parallel NEB
   
   USE KMCNEBProcessSearch, ONLY : AddKMCProcessInitialState,AddKMCProcessFinalState,ProcessNEB, &
      PrcInitialCoord,PrcSaddleCoord,PrcFinalCoord,ProcessAtom !for serial NEB
   IMPLICIT NONE
   INTEGER :: ifile=0,istage=0,nfiles=0
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PerformNEB()
      IMPLICIT NONE
      INTEGER :: ifile1
      REAL :: ActivationEnergy,MaxDisplacement
      
      DO ifile1=1,nfiles
         CALL ProcessNEB(MinDisplacement=0.,ActivationEnergy=ActivationEnergy, &
            MaxDisplacement=MaxDisplacement,ifile=ifile1,ReuseVL=.FALSE.)
      END DO
   END SUBROUTINE PerformNEB
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateUserProcessListFull()
   !used for create process list of full system
   !Last checked: Jan 09, 2011
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      REAL :: ActivationEnergy,MaxDisplacement
      !TYPE(KMCProcessList), POINTER :: ptr
      
      IF (EnabledNEBWithUserKMCProcess) THEN !do in three stages
         DO istage=1,3
            ifile=0
            WRITE(*,*) "istage",istage
            IF (istage==2) THEN !perform NEB calculations
               !DO ifile=1,nfiles
               !   CALL ProcessNEB(MinDisplacement=0.,ActivationEnergy=ActivationEnergy, &
               !      MaxDisplacement=MaxDisplacement,ifile=ifile,ReuseVL=.FALSE.)
               !END DO
               !CALL PerformNEB() !for serial
               CALL PerformNEBParallel(nfiles) !for parallel
            ELSE
               atom=>KMC_AL
               DO WHILE (ASSOCIATED(atom))
                  CALL ResetAffectedCellInfo() !NAffectedCells=0
                  CALL GenerateUserHopProcessList(atom)
                  !CALL ResetAffectedCellInfo() !NAffectedCells=0
                  !CALL GenerateUserExchangeProcessList(atom)
                  atom=>atom%NextNeigh
               END DO
               nfiles=ifile
               WRITE(6,*) "Number of files:",nfiles, " @ end of stage ",istage
               !WRITE(UnitScrap,*) "Printing the NewProcessList @ end of stage ",istage
               !ptr=>NewProcessList
               !DO WHILE (ASSOCIATED(ptr))
               !   WRITE(Unit=UnitScrap,FMT='(I5,",")',ADVANCE="NO") ptr%Process%Index
               !   ptr=>ptr%NextNeigh
               !END DO
            END IF
         END DO
      ELSE !do in one stage
         istage=0
         DO WHILE (ASSOCIATED(atom))
            CALL ResetAffectedCellInfo() !NAffectedCells=0
            CALL GenerateUserHopProcessList(atom)
            !CALL ResetAffectedCellInfo() !NAffectedCells=0
            !CALL GenerateUserExchangeProcessList(atom)
            atom=>atom%NextNeigh
         END DO
      END IF
   END SUBROUTINE GenerateUserProcessListFull
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateUserProcessListCell(CellIndex)
   !used for create process list for atoms in a cell
      IMPLICIT NONE
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(KMCAtom), POINTER :: atom
      INTEGER :: CellIndex(3),NAffectedAtomCoord,i !,ifile1,nfiles1
      REAL :: AffectedAtomCoord(3000),coord(3),tol
      REAL :: ActivationEnergy,MaxDisplacement
      
      !ifile=0 set in kmc_algorithm
      !istage is set in kmc_algorithm
      !ifile=ifile1
      !nfiles=nfiles1
      
      IF (.NOT. EnabledNEBWithUserKMCProcess) THEN
         IF (istage>1) RETURN !otherwise completed in 1 stage
         istage=0 !when istage is 1
      END IF
      
      tol=RadialTolerance
      NAffectedAtomCoord=0
      AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL      
      DO WHILE (ASSOCIATED(AL1))
         NAffectedAtomCoord=NAffectedAtomCoord+1
         AffectedAtomCoord(3*NAffectedAtomCoord-2:3*NAffectedAtomCoord)=AL1%Atom%Coord
         AL1=>AL1%NextNeigh
      END DO
      
      IF (istage==2) THEN !perform NEB calculations
         !WRITE(*,*) "Number of files:",nfiles
         !DO ifile1=1,nfiles
         !   CALL ProcessNEB(MinDisplacement=0.,ActivationEnergy=ActivationEnergy, &
         !      MaxDisplacement=MaxDisplacement,ifile=ifile1,ReuseVL=.TRUE.)
         !END DO
         WRITE(6,*) "Err>> Stage2 attempted in GenerateUserProcessListCell"
         STOP
      ELSE
         DO i=1,NAffectedAtomCoord
            coord=AffectedAtomCoord(3*i-2:3*i)
            atom=>SearchAtom(coord,tol)
            CALL ResetAffectedCellInfo() !NAffectedCells=0
            CALL GenerateUserHopProcessList(atom)
            !CALL ResetAffectedCellInfo() !NAffectedCells=0
            !CALL GenerateUserExchangeProcessList(atom)
         END DO
         nfiles=ifile !keeps getting updated as each cell is accounted for during stage=0, 1 and 3
         !ifile1=ifile
         !nfiles1=nfiles
      END IF
   END SUBROUTINE GenerateUserProcessListCell
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateUserHopProcessList(atom)
   !Creates only hop process on (100) surface
   !Last checked: Jan 09, 2011
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      REAL :: displacement(3),a
      
      NAffectedCells=0
      a=LatticeConst(1)
      displacement(1)=0.5*a
      displacement(2)=displacement(1)
      displacement(3)=0.
      IF (atom%IsMoving) THEN
         CALL CreateHopProcess(atom,displacement,istage)
         displacement(2)=-displacement(2)
         CALL CreateHopProcess(atom,displacement,istage)
         displacement(1)=-displacement(1)
         CALL CreateHopProcess(atom,displacement,istage)
         displacement(2)=-displacement(2)
         CALL CreateHopProcess(atom,displacement,istage)
      END IF
   END SUBROUTINE GenerateUserHopProcessList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateUserExchangeProcessList(atom1)
   !Creates only exchange process on (100) surface
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom1,atom2
      REAL :: displacement1(3),displacement2(3),a,coord2(3),tol
      
      a=LatticeConst(1)
      tol=RadialTolerance
      
      NAffectedCells=0
      coord2=atom1%Coord+a*(/0.5,0.,-0.44/)
      coord2=GetPBCCoord(coord2)
      atom2=>SearchAtom(coord2,tol)
      IF (ASSOCIATED(atom2)) THEN
         displacement1=a*(/0.5,0.,-0.44/)
         displacement2=a*(/0.5,0.,0.44/)
         IF (atom1%IsMoving .AND. atom2%IsMoving) THEN
            CALL CreateExchangeProcess(atom1,displacement1,atom2,displacement2,istage)
         END IF
      ELSE
         IF (atom1%cfg%ConfigurationIndex==1) THEN
            WRITE(6,*) "Unable to find the atom coordinate:", Coord2
            CALL PrintAtoms(PrintAtomsFormat,.FALSE.)
         END IF
      END IF
      
      NAffectedCells=0
      coord2=atom1%Coord+a*(/0.,0.5,-0.44/)
      coord2=GetPBCCoord(coord2)
      atom2=>SearchAtom(coord2,tol)
      IF (ASSOCIATED(atom2)) THEN
         displacement1=a*(/0.,0.5,-0.44/)
         displacement2=a*(/0.,0.5,0.44/)
         IF (atom1%IsMoving .AND. atom2%IsMoving) THEN
            CALL CreateExchangeProcess(atom1,displacement1,atom2,displacement2,istage)
         END IF
      ELSE
         IF (atom1%cfg%ConfigurationIndex==1) THEN
            WRITE(6,*) "Unable to find the atom coordinate:", Coord2
            CALL PrintAtoms(PrintAtomsFormat,.FALSE.)
         END IF
      END IF
      
      NAffectedCells=0
      coord2=atom1%Coord+a*(/-0.5,0.,-0.44/)
      coord2=GetPBCCoord(coord2)
      atom2=>SearchAtom(coord2,tol)
      IF (ASSOCIATED(atom2)) THEN
         displacement1=a*(/-0.5,0.,-0.44/)
         displacement2=a*(/-0.5,0.,0.44/)
         IF (atom1%IsMoving .AND. atom2%IsMoving) THEN
            CALL CreateExchangeProcess(atom1,displacement1,atom2,displacement2,istage)
         END IF
      ELSE
         IF (atom1%cfg%ConfigurationIndex==1) THEN
            WRITE(6,*) "Unable to find the atom coordinate:", Coord2
            CALL PrintAtoms(PrintAtomsFormat,.FALSE.)
         END IF
      END IF
      
      NAffectedCells=0
      coord2=atom1%Coord+a*(/0.,-.5,-0.44/)
      coord2=GetPBCCoord(coord2)
      atom2=>SearchAtom(coord2,tol)
      IF (ASSOCIATED(atom2)) THEN
         displacement1=a*(/0.,-0.5,-0.44/)
         displacement2=a*(/0.,-0.5,0.44/)
         IF (atom1%IsMoving .AND. atom2%IsMoving) THEN
            CALL CreateExchangeProcess(atom1,displacement1,atom2,displacement2,istage)
         END IF
      ELSE
         IF (atom1%cfg%ConfigurationIndex==1) THEN
            WRITE(6,*) "Unable to find the atom coordinate:", Coord2
            CALL PrintAtoms(PrintAtomsFormat,.FALSE.)
         END IF
      END IF
      
   END SUBROUTINE GenerateUserExchangeProcessList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DisplaceAtoms(atom,displacement)
   !displaces atom to a new position
   !Last checked: Jan 12, 2011
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      REAL :: displacement(3)
      
      atom%Coord=GetPBCCoord(atom%Coord+displacement)
      CALL RefreshCell(atom)
   END SUBROUTINE DisplaceAtoms
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE UpdateCfg()
   !Updates the configuration of atoms in affected cells
   !Last checked: Jan 12, 2011
      IMPLICIT NONE
      INTEGER :: i,CellIndex(3)
      TYPE(KMCAtomList), POINTER :: AL1
      
      DO i=1,NAffectedCells
         CellIndex=AffectedCells(3*i-2:3*i)
         AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
         DO WHILE (ASSOCIATED(AL1))
            CALL AssignConfig(AL1%Atom)
            AL1=>AL1%NextNeigh
         END DO
      END DO
   END SUBROUTINE UpdateCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   INCLUDE "kmc_generateprocesslistUser.Hop.f90"
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   INCLUDE "kmc_generateprocesslistUser.Exchange.f90"  
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
END MODULE UserProcess
