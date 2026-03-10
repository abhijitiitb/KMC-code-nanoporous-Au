MODULE checkdb
   USE VARIABLE_TYPE
   USE utilities
   IMPLICIT NONE
   
   INTERFACE Check
      MODULE PROCEDURE ChkNonZeroMag
   END INTERFACE Check
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ChkNonZeroMag(r)
   
      IMPLICIT NONE
      REAL(dp) :: r
      LOGICAL :: ChkNonZeroMag
      
      ChkNonZeroMag= (ABS(r)>1.e-3)
      
   END FUNCTION ChkNonZeroMag
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CheckLL(AL)
   
      !Checks if the atom supposed to be present in a cell
      !is present in the linked list for the cells
      
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(LinkedListContainer), POINTER :: LL
      LOGICAL :: CheckLL,Found
      REAL(dp) :: CellSize(3)
      INTEGER :: i,r1,r2,Cell(3),indx,NCell(3)
      INTEGER, DIMENSION(:), POINTER :: Counter
      
      !WRITE(*,*) "Check>> Checking LL for correctness"
      LL=>AL%LL
      
      CheckLL=.TRUE.
      
      CALL PBC(AL)
      
      CellSize=LL%CellSize
      NCell=LL%NumberCell3D

      CALL MakeSize(Counter,LL%NumberCell)
      
      DO i=1,AL%NAtoms
         Cell=CEILING((AL%AtomCoord(3*i-2:3*i)-AL%BoxSize(4:6))/CellSize)
         WHERE (Cell==0) Cell=NCell
         
         indx=((Cell(1)-1)*NCell(2)+Cell(2)-1)*NCell(3)+Cell(3)
         
         r1=LL%MaxAtomPerCell*(indx-1)+1
         r2=r1-1+LL%ListRange(indx)
         
         Found=ANY(LL%List(r1:r2)==i)
         IF (.NOT. Found) THEN
            WRITE(UNIT=6,FMT='(" >> Unable to find atom:",i5)') i
            WRITE(UNIT=6,FMT='(" >> offending cell:",i5)') indx
            WRITE(UNIT=6,FMT='(" >> supposed cell:",i5)') LL%AtomCellIndx(i)
            WRITE(UNIT=6,FMT='(" >> atom position:",3ES15.8)') AL%AtomCoord(3*i-2:3*i)
            WRITE(UNIT=6,FMT='(" >> cell position found:",3I5)') Cell
            CheckLL=.FALSE.
            CALL PrintList(indx)
            STOP
         END IF
         Counter(indx)=Counter(indx)+1
         
      END DO
      
      IF (ANY(Counter/=LL%ListRange)) THEN
         WRITE(*,*) "$ChkErr>> Linked list count is messed up"
         STOP
      END IF
      
      DEALLOCATE(Counter)
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE PrintList(i)
         IMPLICIT NONE
         INTEGER :: i,j
         
         WRITE(UNIT=*,FMT='("Printing list for cell:",i5)') i
         DO j=r1,r2
            WRITE(UNIT=*,FMT='(i5)') LL%List(j)
         END DO
      END SUBROUTINE PrintList
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      
   END FUNCTION CheckLL
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   FUNCTION CheckVL(AL)
      !can check for both half and full list -- generates the list without a half list
   
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(LinkedListContainer), POINTER :: LL
      INTEGER :: i,j,k,n,r0,r1,r2,NAtoms,nx,ny,nz
      REAL(dp) :: CutOff,drmag,dr(3),BoxSize(3)
      INTEGER :: ri0,ri1,ri2,rj0,rj1,rj2
      LOGICAL :: CheckVL,Found,halflist,LargePotentialCutOff(3)
      INTEGER, DIMENSION(:), POINTER :: counter
      
      !FOR SOME REASON THIS FUNCTION IS GIVING PROBLEMS WITH MPI
      LL=>AL%LL
      VL=>AL%VL
      
      NAtoms=AL%NAtoms
      
      IF (VL%ListType==0) THEN
         WRITE(6,*) "$Err>> VL type (half/full) is not clear"
         STOP
      ELSEIF (VL%ListType==1) THEN
         halflist=.FALSE.
         WRITE(6,*) "VL List type: FULL"
      ELSEIF (VL%ListType==2) THEN
         halflist=.TRUE.
         WRITE(6,*) "VL List type: HALF"
      ELSE
         WRITE(6,*) "$Err>> Incorrect VL type prvoided"
         STOP
      END IF
      ALLOCATE(counter(NAtoms)) !this will store the number of atoms in the list as it is being prepared
      
      CheckVL=.TRUE.
      
      CALL PBC(AL)
      
      CutOff=VL%CutOff+VL%Buffer !note we add the buffer, though an atom from outside the
        !cutoff could have entered the buffer region while we had not updated the linked list
      counter=0
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      
      !check if the atoms are correctly present
      CutOff=VL%CutOff+VL%Buffer
      LargePotentialCutOff(1)=CutOff*2._dp>BoxSize(1)
      LargePotentialCutOff(2)=CutOff*2._dp>BoxSize(2)
      LargePotentialCutOff(3)=CutOff*2._dp>BoxSize(3)
      
      nx=FLOOR(CutOff/BoxSize(1))
      ny=FLOOR(CutOff/BoxSize(2))
      nz=FLOOR(CutOff/BoxSize(3))
      
      DO i=1,NAtoms-1
         
         INCLUDE "checkvl.f90"
         !WRITE(UNIT=6,FMT='("List is correct for atom # ",I5)') i
         
      END DO
      
      !check if the count is correct
      DO i=1,NAtoms
         IF (counter(i)/=VL%ListRange(i)) THEN
            WRITE(6,*) "$ChkErr>> Incorrect # neighbors in Verlet list"
            WRITE(6,*) " >> found:",counter(i)
            WRITE(6,*) " >> supposed to be:",VL%ListRange(i)
            CALL PrintList(i)
            CheckVL=.FALSE.
         END IF
      END DO
      
      IF (.NOT. CheckVL) STOP
      DEALLOCATE(counter)
      
      CONTAINS
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      
      SUBROUTINE PrintList(iatom)
         IMPLICIT NONE
         INTEGER :: j,iatom,jatom
         
         !RETURN
         
         WRITE(UNIT=6,FMT='("Printing Verlet list for atom:",i5)') iatom
         WRITE(UNIT=6,FMT='(" >> Coord  atom:",3ES15.8)') &
           AL%AtomCoord(3*iatom-2:3*iatom)
         WRITE(UNIT=6,FMT='(" >> Atom belongs to cell :",i5)') LL%AtomCellIndx(iatom)
         DO j=ri1,ri2
            jatom=VL%List(j)
            dr=AL%AtomCoord(3*jatom-2:3*jatom)-AL%AtomCoord(3*iatom-2:3*iatom)
            dr=dr-BoxSize*NINT(dr/BoxSize)
            drmag=SQRT(DOT_PRODUCT(dr,dr))
            !!!IF (ABS(drmag-VL%drmag(j))>0.1_dp) 
            
            WRITE(6,*) "Ag ",VL%dr(3*j-2:3*j)
            !WRITE(UNIT=*,FMT='(i5," ",5ES10.3,"    ",L)') VL%List(j), &
            !   VL%drmag(j),drmag,AL%AtomCoord(3*jatom-2:3*jatom)
         END DO
      
      END SUBROUTINE PrintList
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      
   END FUNCTION CheckVL
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
END MODULE checkdb
