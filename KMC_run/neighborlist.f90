MODULE NeighborList

   !returns the verlet list for the atoms provided
   !uses modified version of A&T
   
   !To do:
   !Ability to rebuild the lists if memory issues occur

   USE VARIABLE_TYPE
   USE db_manipulate
   USE checkdb
   USE utilities, ONLY : PBC,PBCdistance

   IMPLICIT NONE
   REAL(dp), PRIVATE, SAVE :: CutOff,CutOff2,Buffer
   LOGICAL, PRIVATE :: CreateHalfList,ByPassVLCreation
   REAL(dp), DIMENSION(:), POINTER, PRIVATE :: OrigCoord,AtomCoord
   REAL(dp), DIMENSION(3), PRIVATE :: BoxSize=0._dp
   INTEGER, PARAMETER, PRIVATE :: MaxAtomsVL=50000000
   
   INTERFACE UpdateOldVL
      MODULE PROCEDURE  UpdateOldVLFull,UpdateOldVLLocal
   END INTERFACE UpdateOldVL
   
   CONTAINS

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckAtomOverlap(AL,distance,NAtoms) !for first NAtoms
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      REAL(dp) :: distance !distance between two atoms
      REAL(dp) :: overlap,mindistance
      INTEGER :: iatom,LR,MaxAtomPerAtom,NAtoms
      INTEGER, DIMENSION(:), POINTER :: VLListRange
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag
      
      CALL AddVerletList(AL,NAtoms=AL%NAtoms)
      VL=>AL%VL
      VLdrmag=>VL%drmag
      VLListRange=>VL%ListRange
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      overlap=1000._dp
      DO iatom=1,NAtoms
         LR=VLListRange(iatom)
         IF (LR>0) THEN
            mindistance=MINVAL(VLdrmag(MaxAtomPerAtom*(iatom-1)+1:MaxAtomPerAtom*(iatom-1)+LR))
            overlap=MIN(overlap,mindistance)
         END IF
      END DO
      WRITE(6,*) "Overlap of atoms is ",overlap
      IF (overlap<distance) THEN
         WRITE(*,*) "Atom overlap is ",overlap
         STOP
      END IF
   END SUBROUTINE CheckAtomOverlap
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MinimumNearestNeighborDistanceVL(AL,NAtoms)
   !finds the nearest neighbor from VL
   !subroutine is not completely tested
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag
      INTEGER, DIMENSION(:), POINTER :: VLListRange
      INTEGER :: NAtoms,iatom,MaxAtomPerAtom,i
      REAL(dp) :: drmagmin,mindistance
      
      VL=>AL%VL
      VLdrmag=>VL%drmag
      VLListRange=>VL%ListRange
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      drmagmin=1000._dp
      DO iatom=1,NAtoms
         mindistance=MINVAL(VLdrmag( (iatom-1)*MaxAtomPerAtom+1: (iatom-1)*MaxAtomPerAtom+ &
            VLListRange(iatom)))
         drmagmin=MIN(mindistance,drmagmin)
         IF (mindistance<0.1_dp) THEN
            WRITE(6,*) "iatom distance:",mindistance,iatom
            DO i=1,VLListRange(iatom)
               WRITE(6,*) "VLdrmag:",VLdrmag((iatom-1)*MaxAtomPerAtom+i),VL%List((iatom-1)*MaxAtomPerAtom+i)
            END DO
         END IF
      END DO
      WRITE(6,*) ">> Smallest separation between atoms is:",drmagmin
   END SUBROUTINE MinimumNearestNeighborDistanceVL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddVLChOS(chos)
      IMPLICIT NONE
      TYPE(ChOSContainer), POINTER :: chos
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: nimages,NAtoms
      
      nimages=0
      AL=>chos%AL
      
      DO WHILE (ASSOCIATED(AL))
         nimages=nimages+1
         NAtoms=AL%NAtoms
         CALL AddVerletList(AL,IsUninitialized=.TRUE.,NAtoms=NAtoms)
         AL=>AL%NextNeigh
      END DO
      
      IF (nimages/=chos%nimages) THEN
         WRITE(6,*) "$Err>> Number of images in ChOS is incorrect"
         STOP
      END IF
   END SUBROUTINE AddVLChOS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddVerletList(AL,IsUninitialized,ListType,NAtoms)

      !if Verlet and linked list already exist but require
      !minor updating then pass IsUninitialized1=.FALSE.
      !CheckVL1 - check if any atom outside buffer
      
      !Both half and full list can be created by selecting setting the option
      !AL%VL%ListType to 2

      IMPLICIT NONE
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(LinkedListContainer), POINTER :: LL
      TYPE(SystemContainer), POINTER :: AL
      TYPE(InteracPotential), POINTER :: Potential
      INTEGER :: iunit,i,NAtoms,rn1,rn2
      INTEGER, OPTIONAL :: ListType
      LOGICAL, OPTIONAL :: IsUninitialized
      LOGICAL :: IsNewVL,IsUninitialized1,UpdateVL
      
      iunit=6
      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(*,*) "$Err>> AtomList is not associated in AddVerletList"
         STOP
      END IF
      
      
      IF (PRESENT(IsUninitialized)) THEN
         IsUninitialized1=IsUninitialized
      ELSE
         IsUninitialized1=.TRUE.
      END IF
      
      ByPassVLCreation=.FALSE.
      IF (NAtoms>MaxAtomsVL) THEN !if too many atoms then VL is memory-expensive so use LL
         WRITE(UNIT=6,FMT='("Number of atoms is greater than ",i8," ... switching over to LL mode")') MaxAtomsVL
         ByPassVLCreation=.TRUE.
      END IF
      
      IF (IsUninitialized1) THEN
      
         Potential=>AL%Potential

         IF (.NOT. ASSOCIATED(Potential)) THEN
            WRITE(6,*) "$Err>> AtomList is not associated in AddVerletList"
            STOP
         END IF
         
!         CutOff=MAXVAL(Potential%MaxPotentialCutoff)
         Cutoff=5.50_dp !For FCC nn    !earlier value 5.50_dp
         
!         WRITE(6,*)CutOff
         Buffer=Potential%VLBuffer
!         WRITE(6,*)Buffer

         IF (CutOff<=0._dp .OR. Buffer<=0._dp) THEN
            WRITE(UNIT=iunit, &
              FMT='(" $Err>> Unreasonable cutoff or buffer passed to AddVerletList",2ES12.3)') &
              CutOff,Buffer
            STOP
         END IF

         !WRITE(6,*) ">> Setting up new Verlet & Linked list"
         !WRITE(UNIT=6,FMT='(" ... cutoff is ",F8.3)') CutOff
         !WRITE(UNIT=6,FMT='(" ... Verlet list buffer is ",F8.3)') Buffer

         CALL LLInitialize(AL,NAtoms)
         LL=>AL%LL

         CALL VLInitialize(AL,NAtoms)
         VL=>AL%VL
         
         VL%CutOff=CutOff
         VL%Buffer=Buffer

      ELSE

         LL=>AL%LL
         VL=>AL%VL

         Buffer=VL%Buffer
         CutOff=VL%CutOff

      END IF
      
      !NAtoms=AL%NAtoms
      rn1=1
      rn2=3*NAtoms
      CutOff=CutOff+Buffer !now cutoff includes the buffer
      CutOff2=CutOff*CutOff
      AtomCoord=>AL%AtomCoord
      IF (.NOT. ByPassVLCreation) OrigCoord=>VL%OrigCoord
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      
      UpdateVL=.FALSE.
      
      !AL%VL%ListType=1
      !IF (AL%VL%ListType==0) AL%VL%ListType=2 !half list - default
      IF (PRESENT(ListType)) THEN
         IF (.NOT.(ListType==1 .OR. ListType==2)) THEN
            WRITE(UNIT=6,FMT='("$Err>> Verlet List type does not make sense ",I3)') ListType
            STOP
         END IF
         IF (VL%ListType>0 .AND. VL%ListType/=ListType) THEN
            WRITE(UNIT=6,FMT='("$Warning>> Changing Verlet list from type",I2," to type",I2)') &
              VL%ListType,ListType
         END IF
         VL%ListType=ListType
      ELSE
         IF (VL%ListType==0) VL%ListType=VLListTypeDefault !half list - default
      END IF
      
      CreateHalfList= VL%ListType==2 !or else it must have been type full list
      
      IF (IsUninitialized1) THEN
         UpdateVL=.TRUE.
      ELSE !if VL is already initialized then useful to get this
         UpdateVL=IsVLMismatched(AL,NAtoms=NAtoms) !automatic check to see if atoms beyond buffer
      END IF
      
      IF (UpdateVL) THEN !rebuild entire LL and VL
         !rather than fill the entire list again
         !one should identify which atoms moved the most
         CALL LLFill(AL,LL,NAtoms)
         
         IF (.NOT. ByPassVLCreation) THEN
            CALL VLFill(AL,LL,VL)
         
            !VL%OrigCoord=AL%AtomCoord
            OrigCoord(rn1:rn2)=AtomCoord(rn1:rn2)
         END IF
         
      ELSE !only find the updated distances without updating LL
         IF (.NOT. ByPassVLCreation) CALL UpdateOldVL(AL,VL,NAtoms)
      END IF
      
      AL%VL=>VL

   END SUBROUTINE AddVerletList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddVerletListQuick(AL,NAtoms)

      !faster version of AddVerletList -- assumes VL was created before
      !and makes fewer checks

      IMPLICIT NONE
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(LinkedListContainer), POINTER :: LL
      TYPE(SystemContainer), POINTER :: AL
      TYPE(InteracPotential), POINTER :: Potential
      INTEGER :: iunit,i,CheckVLFreq1,NAtoms,rn1,rn2
      LOGICAL :: IsNewVL,IsUninitialized,UpdateVL,CheckVL
      
      iunit=6 
      
      ByPassVLCreation=.FALSE.
      IF (NAtoms>MaxAtomsVL) THEN !if too many atoms then VL is memory-expensive so use LL
         ByPassVLCreation=.TRUE.
      END IF
      
      LL=>AL%LL
      VL=>AL%VL

      Buffer=VL%Buffer
      CutOff=VL%CutOff
      
      CutOff=CutOff+Buffer !now cutoff includes the buffer
      CutOff2=CutOff*CutOff
      AtomCoord=>AL%AtomCoord
      IF (.NOT. ByPassVLCreation) OrigCoord=>VL%OrigCoord
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      
      UpdateVL=.FALSE.
      
      !AL%VL%ListType=1
      !IF (AL%VL%ListType==0) AL%VL%ListType=2 !half list - default
      IF (AL%VL%ListType==0) THEN !this will happen if VL is uninitialized
         WRITE(6,*) "$Err>> Verlet list type is zero in AddVerletListQuick"
         STOP
      END IF
      
      CreateHalfList= VL%ListType==2 !or else it must have been type full list
      
      IF (.NOT. ByPassVLCreation) UpdateVL=IsVLMismatchedQuick(AL,NAtoms)
      
      IF (UpdateVL) THEN !rebuild entire LL and VL
         CALL LLFill(AL,LL,NAtoms)
         IF (.NOT. ByPassVLCreation) THEN
            CALL VLFill(AL,LL,VL)
            rn1=1
            rn2=3*NAtoms
            OrigCoord(rn1:rn2)=AtomCoord(rn1:rn2)
         END IF
      ELSE !only find the updated distances without updating LL
         IF (.NOT. ByPassVLCreation) CALL UpdateOldVL(AL,VL,NAtoms)
      END IF

   END SUBROUTINE AddVerletListQuick
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION IsVLMismatched(AL,mask,NAtoms) !AtomCoord and OrigCoord values are set in AddVerletList(Quick)

      !use mask to identify which atoms moved

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      REAL(dp) :: dr(3),dr2mag,dr2mag_max,Buffer2
      INTEGER :: r1,r2,i,NAtoms
      LOGICAL :: IsVLMismatched
      LOGICAL, DIMENSION(:), POINTER, OPTIONAL :: mask

      VL=>AL%VL

      Buffer2=Buffer*Buffer/4._dp !added this factor of 0.25 to make it more sensitive -- Aug 25
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)

      r1=-2
      r2=0
      
      dr2mag_max=0._dp

      DO i=1,NAtoms

         r1=r1+3
         r2=r2+3

         !dr=ABS(PBCdistance(AL,AL%AtomCoord(r1:r2),VL%OrigCoord(r1:r2)))
         dr=AtomCoord(r1:r2)-OrigCoord(r1:r2)
         dr=dr-BoxSize*NINT(dr/BoxSize)

         dr2mag=DOT_PRODUCT(dr,dr)
         dr2mag_max=MAX(dr2mag,dr2mag_max)

         IF (PRESENT(mask)) THEN
            mask(i)=dr2mag>Buffer2 !store which atom went outside buffer
            IsVLMismatched=IsVLMismatched .OR. mask(i)
         ELSE
            IsVLMismatched=dr2mag>Buffer2
            IF (IsVLMismatched) RETURN
         END IF

      END DO

   END FUNCTION IsVLMismatched
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION IsVLMismatchedQuick(AL,NAtoms)

      !use mask to identify which atoms moved

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      REAL(dp) :: dr(3),dr2mag,Buffer2
      INTEGER :: r1,r2,i,NAtoms
      LOGICAL :: IsVLMismatchedQuick

      VL=>AL%VL

      Buffer2=Buffer*Buffer/4._dp !added this factor of 0.25 to make it more sensitive -- Aug 25
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)

      r1=-2
      r2=0

      DO i=1,NAtoms

         r1=r1+3
         r2=r2+3

         !dr=ABS(PBCdistance(AL,AL%AtomCoord(r1:r2),VL%OrigCoord(r1:r2)))
         dr=AtomCoord(r1:r2)-OrigCoord(r1:r2)
         dr=dr-BoxSize*NINT(dr/BoxSize)

         dr2mag=DOT_PRODUCT(dr,dr)

         IsVLMismatchedQuick=dr2mag>Buffer2
         IF (IsVLMismatchedQuick) RETURN

      END DO
   END FUNCTION IsVLMismatchedQuick
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LLInitialize(AL,NAtoms)

      !Sets up the LL for an AL

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(LinkedListContainer), POINTER :: LL
      REAL(dp), DIMENSION(3) :: BoxSize
      INTEGER :: NAtoms
      
      !WRITE(*,*) ">> Initializing linked list"

      CALL MakeSize(AL%LL,NAtoms) !set size of LL%AtomCellIndx to NAtoms
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)

      LL=>AL%LL

      LL%CellSize=CutOff/2._dp !use interaction potential cutoff as reference
      LL%NumberCell3D=FLOOR((AL%BoxSize(1:3)-AL%BoxSize(4:6))/ &
           LL%CellSize)

      !assuming fcc of lattice constant a=3 Ang 
      !  (2 atoms per unitcell)
      !WRITE(UNIT=6,FMT='(" ...Number of cells for LL:",3I3)') LL%NumberCell3D
      LL%NumberCell=PRODUCT(LL%NumberCell3D)
      LL%CellSize=(AL%BoxSize(1:3)-AL%BoxSize(4:6))/ &
           REAL(LL%NumberCell3D,dp)
           
      !LL%MaxAtomPerCell=CEILING(2.5_dp*PRODUCT(LL%CellSize)/27._dp)*3
      LL%MaxAtomPerCell=CEILING(4._dp*PRODUCT(LL%CellSize)/8._dp) !assuming a FCC unit cell of lattice const. 2.0 Ang
      !LL%MaxAtomPerCell=100

      LL%MaxAtom=LL%NumberCell*LL%MaxAtomPerCell
      
      !conclude
      CALL MakeSize(LL%ListRange,LL%NumberCell)
      CALL MakeSize(LL%List,LL%MaxAtom)
      
      LL%ListRange=0
      LL%List=0
      
      LL%UnInitialized=.FALSE.

   END SUBROUTINE LLInitialize
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LLInitializeCustomized(AL,LL,CellSize,NAtoms)

      !Sets up the LL for an AL (used by parallel MD domain decomposition)
      !NAtoms is the target number of atoms that are possible in AL

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(LinkedListContainer), POINTER :: LL
      REAL(dp), DIMENSION(3) :: BoxSize,CellSize
      INTEGER :: NAtoms
      
      IF (ASSOCIATED(LL)) THEN
         WRITE(6,*) "Err>> LL is already associated in LLInitializedCustomized"
         STOP
      END IF

      CALL MakeSize(AL%LL,NAtoms) !set size of LL%AtomCellIndx to NAtoms
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      
      LL=>AL%LL

      LL%CellSize=CellSize !use interaction potential cutoff as reference
      LL%NumberCell3D=FLOOR(BoxSize/CellSize)
      LL%NumberCell=PRODUCT(LL%NumberCell3D)
      IF (ANY(ABS(CellSize*REAL(LL%NumberCell3D)-BoxSize)>0.0001_dp)) THEN
         WRITE(6,*) "Err>> Box size is not an integer multiple of cell size in LLInitializedCustomized"
         STOP
      END IF
           
      LL%MaxAtomPerCell=CEILING(4._dp*PRODUCT(LL%CellSize)/8._dp) !see LLInitialize
      LL%MaxAtom=LL%NumberCell*LL%MaxAtomPerCell
      
      !conclude
      CALL MakeSize(LL%ListRange,LL%NumberCell)
      CALL MakeSize(LL%List,LL%MaxAtom)
      
      LL%ListRange=0
      LL%List=0
      
      LL%UnInitialized=.FALSE.

   END SUBROUTINE LLInitializeCustomized
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE VLInitialize(AL,NAtoms)

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      INTEGER :: NAtoms
      
      !WRITE(*,*) ">> Initializing Verlet list for # atoms",NAtoms
      
      CALL MakeSize(AL%VL,NAtoms,MaxAtomsVL) !if NAtoms>MaxAtomVL creation of arrays will be bypassed
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)

      VL=>AL%VL

      VL%MaxAtom=VL%NAtoms*VL%MaxAtomPerAtom
      !WRITE(UNIT=*,FMT='(" ...Maximum atoms per atom:",I5)') VL%MaxAtomPerAtom
         !where MaxAtomPerAtom is max atoms in list per atom
      IF (VL%MaxAtom>1000000000) THEN
         WRITE(6,*) "$Err>> Verlet list will be taking up too much space"
         WRITE(6,*) "Number of atoms:",VL%NAtoms,NAtoms
         WRITE(6,*) "Number of atoms per atom:",VL%MaxAtomPerAtom
         STOP
      ELSE
         WRITE(6,*) "Verlet list will be taking below amount of space"
         WRITE(6,*) "Number of atoms:",VL%NAtoms,NAtoms
         WRITE(6,*) "Number of atoms per atom:",VL%MaxAtomPerAtom
      END IF

      VL%UnInitialized=.FALSE.
      
      IF (ByPassVLCreation) RETURN
      
      CALL MakeSize(VL%List,VL%MaxAtom)
      VL%List=0
      CALL MakeSize(VL%ListDomainAtom,VL%MaxAtom)
      VL%ListDomainAtom=.TRUE.
      CALL MakeSize(VL%drmag,VL%MaxAtom)
      VL%drmag=0._dp
      CALL MakeSize(VL%dr,3*VL%MaxAtom)
      VL%dr=0._dp

   END SUBROUTINE VLInitialize
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LLFill(AL,LL,NAtoms)

   !assembles atoms based on their cell location

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(LinkedListContainer), POINTER :: LL
      INTEGER :: NAtoms,i,LLFactor(3),MaxAtomPerCell,NumberCell3D(3)
      INTEGER, DIMENSION(:), POINTER :: AtomCellIndx,ListRange,List
      REAL(dp) :: xc(3),CellSize(3)
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord
      
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      
      !NAtoms=AL%NAtoms
      CALL PBC(AL) !move atoms to PBC
      xc=AL%BoxSize(4:6) !lower corner
      
      !get the cell index
      LLFactor(1)=LL%NumberCell3D(2)*LL%NumberCell3D(3) !x
      LLFactor(2)=LL%NumberCell3D(3) !y
      LLFactor(3)=1 !z

      AtomCellIndx=>LL%AtomCellIndx
      AtomCoord=>AL%AtomCoord
      MaxAtomPerCell=LL%MaxAtomPerCell
      List=>LL%List
      ListRange=>LL%ListRange
      CellSize=LL%CellSize
      NumberCell3D=LL%NumberCell3D
      
      !For each atom add its cell # to AtomCellIndx
      List=0
      ListRange=0 !no atom filled in so far ... 
      
      DO i=1,NAtoms
         CALL LLFillDriver(i)
      END DO
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE LLFillDriver(i)
         IMPLICIT NONE
         INTEGER :: i,r1,r2,cellindx,pos
         INTEGER :: icell(3)
      
         r1=3*i-2
         r2=r1+2
         icell=CEILING((AtomCoord(r1:r2)-xc)/CellSize)
         !WHERE (icell==0) icell=1 !NumberCell3D !for periodic BC
         icell=MAX(icell,1)

         !Convert to scalar number 
         icell(1:2)=icell(1:2)-1
         cellindx=SUM (icell*LLFactor)
         
         !LL%AtomCellIndx(i)=cellindx
         AtomCellIndx(i)=cellindx
         !IF (cellindx<0) THEN
         !   WRITE(6,*) "$Err>> LL fill error ... cell index is negative",icell,AtomCoord(r1:r2)
         !END IF
      
         !Next invert the info in AtomCellIndx and store based on cells
         !LLIndx gives number of atoms inside a cell
         !for ith cell these atoms are recorded in 
         !LL(MaxAtomPerCell*(i-1)+1:MaxAtomPerCell*(i-1)+LLIndx(i))
      
         !add atom to LL%List
         !pos=LL%ListRange(cellindx)+1
         pos=ListRange(cellindx)+1

         IF (pos>MaxAtomPerCell) THEN
            WRITE(6,*) "$Err>> Max # atoms per cell is smaller than needed",pos
            WRITE(UNIT=6,FMT='(" ...# atoms in cell so far:",I8)') pos
            WRITE(UNIT=6,FMT='(" ...maximum #atoms allowed:",I8)') MaxAtomPerCell
            STOP
         END IF

         !LL%List((cellindx-1)*LL%MaxAtomPerCell+pos)=i
         !LL%ListRange(cellindx)=pos !update index
         List((cellindx-1)*MaxAtomPerCell+pos)=i
         ListRange(cellindx)=pos !update index
      
      END SUBROUTINE LLFillDriver
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   END SUBROUTINE LLFill
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE VLFill(AL,LL,VL) 

      !scan cells and add atoms - assumes periodic BC

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(LinkedListContainer), POINTER :: LL
      TYPE(VerletListContainer), POINTER :: VL
      INTEGER, DIMENSION(:), POINTER :: LLListRange,LLList,VLListRange,VLList
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom
      REAL(dp), DIMENSION(:), POINTER :: VLdr,VLdrmag
      INTEGER :: maincellindx,j1,j2,j3,neighcellindx,MaxAtomPerAtom,MaxAtomPerCell
      INTEGER :: nx,ny,nz
      INTEGER, DIMENSION(3) :: c1,c2,NCell,BoxPosition,BoxPosition1,AllowPBC
      LOGICAL :: CellOutsideBox,CellInsideBox,LargePotentialCutOff(3),DoubleCellImage(3)
      
      IF (ByPassVLCreation) RETURN
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)

      VLList=>VL%List
      VLListRange=>VL%ListRange
      VLList=0
      VLListRange=0 !no atom filled so far
      LLList=>LL%List
      LLListRange=>LL%ListRange
      VLdr=>VL%dr
      VLdrmag=>VL%drmag
      VLListDomainAtom=>VL%ListDomainAtom
      VLListDomainAtom=.TRUE.
      
      c1=(/1,1,0/) !starting cell #

      NCell=LL%NumberCell3D
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      MaxAtomPerCell=LL%MaxAtomPerCell
      
      CutOff=VL%CutOff+VL%Buffer
      CutOff2=CutOff*CutOff
      
      !CutOff= 20. !artificial
      !CutOff2=CutOff*CutOff
      
      CreateHalfList=VL%ListType==2
      
      nx=CEILING(CutOff/LL%CellSize(1))
      ny=CEILING(CutOff/LL%CellSize(2))
      nz=CEILING(CutOff/LL%CellSize(3))
      
      LargePotentialCutOff(1)=CutOff*2._dp>BoxSize(1)
      LargePotentialCutOff(2)=CutOff*2._dp>BoxSize(2)
      LargePotentialCutOff(3)=CutOff*2._dp>BoxSize(3)
      
      
     !write(*,*) "LargePoten",LargePotentialCutOff
     !       DoubleCellImage(1)=.NOT. LargePotentialCutOff(1) .AND. 2*nx+1>NCell(1) !is it possible that a neighcell can be observed twice
     !       DoubleCellImage(2)=.NOT. LargePotentialCutOff(2) .AND. 2*ny+1>NCell(2) !if so we cannot allow that since the neighcell atoms
     !       DoubleCellImage(3)=.NOT. LargePotentialCutOff(3) .AND. 2*nz+1>NCell(3) !shall be added twice.
      !write(*,*) "Double",DoubleCellImage
      !write(*,*) "Cellsize",LL%CellSize
      !write(*,*) "NCell",NCell
      !write(*,*) "Cutoff",CutOff
      !write(*,*) "nx,ny,",nx,ny,nz
      !A PROBLEM -- if the Cutoff is slighly less than half the box size then
      ! LargePotentialCutOff is FALSE, but nx can make the same atom show up twice !!
      ! For ex. nx=10, CutOff=4.9 and BoxSize=10, and atom can show up in a list only once
      ! However, in the present version and atom can INCORRECTLY show up more than once!!
      
      IF (CreateHalfList) THEN
      
         !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         IF (ANY(LargePotentialCutOff)) THEN !at least one is long ranged
         !HALF LIST, AT LEAST ONE LARGE POTENTIAL
         
            DoubleCellImage(1)=.NOT. LargePotentialCutOff(1) .AND. 2*nx+1>NCell(1) !is it possible that a neighcell can be observed twice
            DoubleCellImage(2)=.NOT. LargePotentialCutOff(2) .AND. 2*ny+1>NCell(2) !if so we cannot allow that since the neighcell atoms
            DoubleCellImage(3)=.NOT. LargePotentialCutOff(3) .AND. 2*nz+1>NCell(3) !shall be added twice.
         
            IF (ANY(DoubleCellImage)) THEN
            !HALF LIST, LONG INTERACTION ALONG SOME DIMN, SHORT ALONG OTHERS WITH AT LEAST ONE HAVING DOUBLE IMAGE
!write(*,*) "half-longinterac-doubleimage"
               DO maincellindx=1,LL%NumberCell
      
                  c1(3)=c1(3)+1
                  c1(2)=c1(2)+(c1(3)-1)/NCell(3)
                  c1(1)=c1(1)+(c1(2)-1)/NCell(2)
                  c1=MOD(c1-1,NCell)+1
                  
                  IF (LLListRange(maincellindx)==0) CYCLE !no atoms present
                  
                  DO j3=-nz,nz
                     DO j2=-ny,ny
                        DO j1=-nx,nx
                           
                           DoubleCellImage(1)=.NOT. LargePotentialCutOff(1) .AND. 2*ABS(j1)+1>NCell(1) !can this cell be observed twice
                           DoubleCellImage(2)=.NOT. LargePotentialCutOff(2) .AND. 2*ABS(j2)+1>NCell(2)
                           DoubleCellImage(3)=.NOT. LargePotentialCutOff(3) .AND. 2*ABS(j3)+1>NCell(3)
                           
                           IF (ANY(DoubleCellImage)) THEN

                              c2=c1+(/j1,j2,j3/) !present cell + displacement
                              BoxPosition=0
                              BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                              WHERE (DoubleCellImage) BoxPosition=BoxPosition1

                              IF (ALL(BoxPosition==0)) THEN !cell inside box and no double possible
                                 BoxPosition=0
                                 AllowPBC=1
                                 WHERE (LargePotentialCutOff)
                                    BoxPosition=BoxPosition1
                                    AllowPBC=0
                                 END WHERE
                              
                                 !wrap the box using periodic BC
                                 c2=c2-NCell*BoxPosition1 !FLOOR(REAL(c2-1)/REAL(NCell)) 
                                 
                                 neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
                                 
                                 CellOutsideBox=ANY(BoxPosition/=0)
                                 IF (CellOutsideBox) THEN
                                    CALL VLAddHalfCellOutBoxLongCutOff()
                                 ELSE
                                    IF (maincellindx==neighcellindx) THEN
                                       CALL VLAddHalfAACellInBoxLongCutOff()
                                    ELSEIF (maincellindx<neighcellindx) THEN !prevent double
                                       CALL VLAddHalfABCellInBoxLongCutOff()
                                    END IF
                                 END IF
                              END IF

                           ELSE !no double is possible -- only pure long or pure short

                              c2=c1+(/j1,j2,j3/) !present cell + displacement
                              BoxPosition=0
                              BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                              AllowPBC=1
                              WHERE (LargePotentialCutOff)
                                 BoxPosition=BoxPosition1
                                 AllowPBC=0
                              END WHERE
                              
                              !wrap the box using periodic BC
                              c2=c2-NCell*BoxPosition1 !FLOOR(REAL(c2-1)/REAL(NCell)) 
                              
                              neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
                              
                              CellOutsideBox=ANY(BoxPosition/=0)
                              IF (CellOutsideBox) THEN
                                 CALL VLAddHalfCellOutBoxLongCutOff()
                              ELSE
                                 IF (maincellindx==neighcellindx) THEN
                                    CALL VLAddHalfAACellInBoxLongCutOff()
                                 ELSEIF (maincellindx<neighcellindx) THEN !prevent double
                                    CALL VLAddHalfABCellInBoxLongCutOff()
                                 END IF
                              END IF

                           END IF

                        END DO
                     END DO
                  END DO
                  
               END DO

            ELSE !cell cannot be seen double, either it is long potential along the direction or it is short with cell observed only once
!write(*,*) "large potential no double -- half list"
               DO maincellindx=1,LL%NumberCell
      
                  c1(3)=c1(3)+1
                  c1(2)=c1(2)+(c1(3)-1)/NCell(3)
                  c1(1)=c1(1)+(c1(2)-1)/NCell(2)
                  c1=MOD(c1-1,NCell)+1
                  
                  IF (LLListRange(maincellindx)==0) CYCLE !no atoms present
                  
                  DO j3=-nz,nz
                     DO j2=-ny,ny
                        DO j1=-nx,nx
                           
                           c2=c1+(/j1,j2,j3/) !present cell + displacement
                           BoxPosition=0
                           BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                           AllowPBC=1
                           WHERE (LargePotentialCutOff)
                              BoxPosition=BoxPosition1
                              AllowPBC=0
                           END WHERE
                           
                           !wrap the box using periodic BC
                           c2=c2-NCell*BoxPosition1 !FLOOR(REAL(c2-1)/REAL(NCell)) 
                           
                           neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
      
!IF ((maincellindx==1 .OR. neighcellindx==141) .AND. (maincellindx==141 .OR. neighcellindx==1)) THEN
!write(*,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"                     
!write(*,*) "maincellindx:",maincellindx
!write(*,*) "neighcellindx:",neighcellindx
!write(*,*)"j1,j2,:", j1,j2,j3
!write(*,*) "BoxPOsition:",BoxPosition
!write(*,*) "AllowPBC:",AllowPBC
!END IF
                           CellOutsideBox=ANY(BoxPosition/=0)
                           IF (CellOutsideBox) THEN
                              CALL VLAddHalfCellOutBoxLongCutOff()
                           ELSE
                              IF (maincellindx==neighcellindx) THEN
                                 CALL VLAddHalfAACellInBoxLongCutOff()
                              ELSEIF (maincellindx<neighcellindx) THEN !prevent double
                                 CALL VLAddHalfABCellInBoxLongCutOff()
                              END IF
                           END IF
                           
                        END DO
                     END DO
                  END DO
                  
               END DO
            END IF
         !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         
         ELSE !completely short ranged interaction
         !HALF LIST, ALL SHORT RANGED POTENTIAL
         
            DoubleCellImage(1)=2*nx+1>NCell(1) !is it possible that a neighcell can be observed twice
            DoubleCellImage(2)=2*ny+1>NCell(2) !if so we cannot allow that since the neighcell atoms
            DoubleCellImage(3)=2*nz+1>NCell(3) !shall be added twice.
            
            IF (ANY(DoubleCellImage)) THEN !it is possible that cell is observed twice in the do loop
            !HALF LIST, SHORT POTENTIAL -- AT LEAST ONE CAN BE SEEN DOUBLE
!write(*,*) "half-double-short"            
               DO maincellindx=1,LL%NumberCell
               
                  c1(3)=c1(3)+1
                  c1(2)=c1(2)+(c1(3)-1)/NCell(3)
                  c1(1)=c1(1)+(c1(2)-1)/NCell(2)
                  c1=MOD(c1-1,NCell)+1
         
                  IF (LLListRange(maincellindx)==0) CYCLE !no atoms present
                  
                  DO j3=-nz,nz
                     DO j2=-ny,ny
                        DO j1=-nx,nx
                        
                           DoubleCellImage(1)=2*ABS(j1)+1>NCell(1) !can this cell be observed twice
                           DoubleCellImage(2)=2*ABS(j2)+1>NCell(2)
                           DoubleCellImage(3)=2*ABS(j3)+1>NCell(3)
                           
                           IF (ANY(DoubleCellImage)) THEN
                           !since potential is short each atom will appear only once in the VL, 
                           ! and the neighcell will appear twice in the do loop
                           ! when the neighcell is called once it will be inside the box
                           ! another time it is outside the box, in both cases only a subset
                           ! is added (intersection of the two subsets should be a null space)
                           !so we can allow the cell to be considered only once when it is inside the
                           ! box
                        
                              c2=c1+(/j1,j2,j3/) !present cell + displacement
                              !BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                  
                              BoxPosition=0
                              BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                              WHERE (DoubleCellImage) BoxPosition=BoxPosition1
                              
                              !wrap the box using periodic BC
                              c2=c2-NCell*BoxPosition1
                              neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
                              
                              CellInsideBox=ALL(BoxPosition==0)
                              
                              IF (CellInsideBox) THEN
                                 IF (maincellindx==neighcellindx) THEN
                                    CALL VLAddHalfAAShortCutOff()
                                 ELSEIF (maincellindx<neighcellindx) THEN
                                    CALL VLAddHalfABShortCutOff()
                                 END IF
                              END IF
                              
                           ELSE !add this cell using min image convention
                              
                              c2=c1+(/j1,j2,j3/) !present cell + displacement
                              
                              !wrap the box using periodic BC
                              c2=c2-NCell*FLOOR(REAL(c2-1)/REAL(NCell)) 
                              neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
                              IF (maincellindx==neighcellindx) THEN
                                 CALL VLAddHalfAAShortCutOff()
                              ELSEIF (maincellindx<neighcellindx) THEN
                                 CALL VLAddHalfABShortCutOff()
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
         
               END DO
            
            
            ELSE !a neighcell will not be encountered twice
            !HALF LIST, SHORT POTENTIAL -- DOUBLE IS NOT POSSIBLE
            
               DO maincellindx=1,LL%NumberCell
      
                  c1(3)=c1(3)+1
                  c1(2)=c1(2)+(c1(3)-1)/NCell(3)
                  c1(1)=c1(1)+(c1(2)-1)/NCell(2)
                  c1=MOD(c1-1,NCell)+1
                  
                  IF (LLListRange(maincellindx)==0) CYCLE !no atoms present
                  
                  DO j3=-nz,nz
                     DO j2=-ny,ny
                        DO j1=-nx,nx
                        
                           c2=c1+(/j1,j2,j3/) !present cell + displacement
                           !BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                  
                           !wrap the box using periodic BC
                           c2=c2-NCell*FLOOR(REAL(c2-1)/REAL(NCell)) 
                  
                           neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
                  
                           IF (maincellindx==neighcellindx) THEN
                              CALL VLAddHalfAAShortCutOff()
                           ELSEIF (maincellindx<neighcellindx) THEN
                              CALL VLAddHalfABShortCutOff()
                           END IF
                        END DO
                     END DO
                  END DO
         
               END DO
            
            END IF
            
         END IF
         
         !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         
      ELSE !full list
      
         !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         IF (ANY(LargePotentialCutOff)) THEN !at least one is long ranged     
         
            DoubleCellImage(1)=.NOT. LargePotentialCutOff(1) .AND. 2*nx+1>NCell(1) !is it possible that a neighcell can be observed twice
            DoubleCellImage(2)=.NOT. LargePotentialCutOff(2) .AND. 2*ny+1>NCell(2) !if so we cannot allow that since the neighcell atoms
            DoubleCellImage(3)=.NOT. LargePotentialCutOff(3) .AND. 2*nz+1>NCell(3) !shall be added twice.
            
            IF (ANY(DoubleCellImage)) THEN
!write(*,*) "large potential double will occur -- full list"
               DO maincellindx=1,LL%NumberCell
      
                  c1(3)=c1(3)+1
                  c1(2)=c1(2)+(c1(3)-1)/NCell(3)
                  c1(1)=c1(1)+(c1(2)-1)/NCell(2)
                  c1=MOD(c1-1,NCell)+1
                  
                  IF (LLListRange(maincellindx)==0) CYCLE !no atoms present
               
                  DO j3=-nz,nz
                     DO j2=-ny,ny
                        DO j1=-nx,nx
                           
                           DoubleCellImage(1)=.NOT. LargePotentialCutOff(1) .AND. 2*ABS(j1)+1>NCell(1) !can this cell be observed twice
                           DoubleCellImage(2)=.NOT. LargePotentialCutOff(2) .AND. 2*ABS(j2)+1>NCell(2)
                           DoubleCellImage(3)=.NOT. LargePotentialCutOff(3) .AND. 2*ABS(j3)+1>NCell(3)

                           IF (ANY(DoubleCellImage)) THEN !this cell will show up double

!since potential is short where DoubleCellImage is TRUE
! each atom will appear only once in the VL for such cells, 
! and the neighcell will appear twice in the do loop.
! When the neighcell is called once it will be inside the box
! another time it is outside the box, in both cases only a subset
! is added (intersection of the two subsets should be a null space).
! So we can allow the cell to be considered only once when it is inside the
! box
                              c2=c1+(/j1,j2,j3/) !present cell + displacement
                              !BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                  
                              BoxPosition=0
                              BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                              WHERE (DoubleCellImage) BoxPosition=BoxPosition1

                              IF (ANY(BoxPosition/=0)) CYCLE !the neighcell has a double image which can be ignored

                              BoxPosition=0
                              AllowPBC=1
                              WHERE (LargePotentialCutOff) 
                                 BoxPosition=BoxPosition1
                                 AllowPBC=0
                              END WHERE
                              
                              !wrap the box using periodic BC
                              c2=c2-NCell*BoxPosition1
                              neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
                              
                              CellOutsideBox=ANY(BoxPosition/=0)
                              IF (CellOutsideBox) THEN
                                 CALL VLAddFullCellOutBoxLongCutOff()
                              ELSE
                                 IF (maincellindx==neighcellindx) THEN
                                    CALL VLAddFullAACellInBoxLongCutOff()
                                 ELSEIF (maincellindx<neighcellindx) THEN
                                    CALL VLAddFullABCellInBoxLongCutOff()
                                 END IF
                              END IF

                           ELSE

                              c2=c1+(/j1,j2,j3/) !present cell + displacement
                              BoxPosition=0
                              BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                              AllowPBC=1
                              WHERE (LargePotentialCutOff) 
                                 BoxPosition=BoxPosition1
                                 AllowPBC=0
                              END WHERE
                              
                              !wrap the box using periodic BC
                              c2=c2-NCell*BoxPosition1 !FLOOR(REAL(c2-1)/REAL(NCell)) 
                              
                              neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
                              
                              CellOutsideBox=ANY(BoxPosition/=0)
                              IF (CellOutsideBox) THEN
                                 CALL VLAddFullCellOutBoxLongCutOff()
                              ELSE
                                 IF (maincellindx==neighcellindx) THEN
                                    CALL VLAddFullAACellInBoxLongCutOff()
                                 ELSEIF (maincellindx<neighcellindx) THEN
                                    CALL VLAddFullABCellInBoxLongCutOff()
                                 END IF
                              END IF

                           END IF
                        END DO
                     END DO
                  END DO
         
               END DO

            ELSE

               DO maincellindx=1,LL%NumberCell
      
                  c1(3)=c1(3)+1
                  c1(2)=c1(2)+(c1(3)-1)/NCell(3)
                  c1(1)=c1(1)+(c1(2)-1)/NCell(2)
                  c1=MOD(c1-1,NCell)+1
                  
                  IF (LLListRange(maincellindx)==0) CYCLE !no atoms present
               
                  DO j3=-nz,nz
                     DO j2=-ny,ny
                        DO j1=-nx,nx
                           
                           c2=c1+(/j1,j2,j3/) !present cell + displacement
                           BoxPosition=0
                           BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                           AllowPBC=1
                           WHERE (LargePotentialCutOff) 
                              BoxPosition=BoxPosition1
                              AllowPBC=0
                           END WHERE
                           
                           !wrap the box using periodic BC
                           c2=c2-NCell*BoxPosition1 !FLOOR(REAL(c2-1)/REAL(NCell)) 
                           
                           neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
                           
                           CellOutsideBox=ANY(BoxPosition/=0)
                           IF (CellOutsideBox) THEN
                              CALL VLAddFullCellOutBoxLongCutOff()
                           ELSE
                              IF (maincellindx==neighcellindx) THEN
                                 CALL VLAddFullAACellInBoxLongCutOff()
                              ELSEIF (maincellindx<neighcellindx) THEN
                                 CALL VLAddFullABCellInBoxLongCutOff()
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
         
               END DO

            END IF
         !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         
         ELSE !completely short ranged interaction
         
            DoubleCellImage(1)=2*nx+1>NCell(1) !is it possible that a neighcell can be observed twice
            DoubleCellImage(2)=2*ny+1>NCell(2) !if so we cannot allow that since the neighcell atoms
            DoubleCellImage(3)=2*nz+1>NCell(3) !shall be added twice.
         
            IF (ANY(DoubleCellImage)) THEN
            
               DO maincellindx=1,LL%NumberCell
      
                  c1(3)=c1(3)+1
                  c1(2)=c1(2)+(c1(3)-1)/NCell(3)
                  c1(1)=c1(1)+(c1(2)-1)/NCell(2)
                  c1=MOD(c1-1,NCell)+1
         
                  IF (LLListRange(maincellindx)==0) CYCLE !no atoms present
         
                  DO j3=-nz,nz
                     DO j2=-ny,ny
                        DO j1=-nx,nx
                        
                           DoubleCellImage(1)=2*ABS(j1)+1>NCell(1) !can this cell be observed twice
                           DoubleCellImage(2)=2*ABS(j2)+1>NCell(2)
                           DoubleCellImage(3)=2*ABS(j3)+1>NCell(3)
                           
                           IF (ANY(DoubleCellImage)) THEN
                           !since potential is short each atom will appear only once in the VL, 
                           ! and the neighcell will appear twice in the do loop.
                           ! When the neighcell is called once it will be inside the box
                           ! another time it is outside the box, in both cases only a subset
                           ! is added (intersection of the two subsets should be a null space)
                           ! So we can allow the cell to be considered only once when it is inside the
                           ! box
                        
                              c2=c1+(/j1,j2,j3/) !present cell + displacement
                              !BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                  
                              BoxPosition=0
                              BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                              WHERE (DoubleCellImage) BoxPosition=BoxPosition1
                              
                              !wrap the box using periodic BC
                              c2=c2-NCell*BoxPosition1
                              neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
                              
                              CellInsideBox=ANY(BoxPosition/=0)
                              
                              IF (CellInsideBox) THEN
                                 IF (maincellindx==neighcellindx) THEN
                                    CALL VLAddFullAAShortCutOff()
                                 ELSEIF (maincellindx<neighcellindx) THEN
                                    CALL VLAddFullABShortCutOff()
                                 END IF
                              END IF
                              
                           ELSE !add this cell using min image convention
                           
                              c2=c1+(/j1,j2,j3/) !present cell + displacement
                              !BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                  
                              !wrap the box using periodic BC
                              c2=c2-NCell*FLOOR(REAL(c2-1)/REAL(NCell)) 
                  
                              neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
                  
                              IF (maincellindx==neighcellindx) THEN
                                 CALL VLAddFullAAShortCutOff()
                              ELSEIF (maincellindx<neighcellindx) THEN
                                 CALL VLAddFullABShortCutOff()
                              END IF
                              
                           END IF
                           
                        END DO
                     END DO
                  END DO
            
               END DO
            
            ELSE !no cell can show up as double and interaction is short ranged
! write(*,*) "large potential no double -- full list"
               DO maincellindx=1,LL%NumberCell
      
                  c1(3)=c1(3)+1
                  c1(2)=c1(2)+(c1(3)-1)/NCell(3)
                  c1(1)=c1(1)+(c1(2)-1)/NCell(2)
                  c1=MOD(c1-1,NCell)+1
         
                  IF (LLListRange(maincellindx)==0) CYCLE !no atoms present
         
                  DO j3=-nz,nz
                     DO j2=-ny,ny
                        DO j1=-nx,nx
               
                           c2=c1+(/j1,j2,j3/) !present cell + displacement
                           !BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
                  
                           !wrap the box using periodic BC
                           c2=c2-NCell*FLOOR(REAL(c2-1)/REAL(NCell)) 
                  
                           neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
                  
                           IF (maincellindx==neighcellindx) THEN
                              CALL VLAddFullAAShortCutOff()
                           ELSEIF (maincellindx<neighcellindx) THEN
                              CALL VLAddFullABShortCutOff()
                           END IF
                        END DO
                     END DO
                  END DO
            
               END DO
            
            END IF
            
         END IF
         
      END IF
      
      
      !This is the general case ... (specialized cases are uncommented above)
      !DO maincellindx=1,LL%NumberCell
      
      !   c1(3)=c1(3)+1
      !   c1(2)=c1(2)+(c1(3)-1)/NCell(3)
      !   c1(1)=c1(1)+(c1(2)-1)/NCell(2)
      !   c1=MOD(c1-1,NCell)+1
         
      !   IF (LLListRange(maincellindx)==0) CYCLE !no atoms present
         
      !   DO j3=-nz,nz
      !      DO j2=-ny,ny
      !         DO j1=-nx,nx
               
               !BoxPosition needs to be half box size based
                
                  !BoxPosition=(/j1,j2,j3/)/NCell !relative position of box containing the neighbor cell
      !            c2=c1+(/j1,j2,j3/) !present cell + displacement
      !            BoxPosition=0
      !            BoxPosition1=FLOOR(REAL(c2-1)/REAL(NCell))
      !            AllowPBC=1
      !            WHERE (LargePotentialCutOff) 
      !               BoxPosition=BoxPosition1
      !               AllowPBC=0
      !            END WHERE
                  
                  !wrap the box using periodic BC
      !            c2=c2-NCell*BoxPosition1 !FLOOR(REAL(c2-1)/REAL(NCell)) 
                  
      !            neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !is this inefficient??
                  
      !            CellOutsideBox=ANY(BoxPosition/=0)
      !            IF (maincellindx<=neighcellindx .OR. CellOutsideBox) CALL VLAdd() !consider a pair only once unless the neighbor cell lies outside the box with long interaction and small box size
      !         END DO
      !      END DO
      !   END DO
         
      !END DO
   
   CONTAINS
   
   INCLUDE "neighborlist_VLAdd.f90"
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   SUBROUTINE VLAdd()
   !this is the general code from which the file neighborlist_VLAdd.f90 has been generated
   !for special cases

      !cutoff used is set at the top
      !rules used:
      !Full list: any atom within cutoff radius of i should be added in the VL of i
      !Half list: any atom j within cutoff radius of i should be added in the VL of i provided i<j
      !  exception for half list: if j lies outside 1 box distance it has to be added even if j<i
      !Notice that the atom j that lies outside 1 box distance then it will not have its list, however,
      !that is perfectly fine. In Ewald summation, the box distance can be greater than one. In such
      !case the neighbors of i are completely available -- hence the pair interactions 
      !(for effectively more than N atoms in the system) is correct
      !cases considered:
      !i<j and a<b, A is main box and B is neighbor box, x means it works
      !Box1   Cell1   Atom1  Box2   Cell2   Atom2  Full Half    Comments
      ! A      a        i     A       a      i      x    x since dr=0 this will not go through
      ! A      a        i     A       a      j      x    x both half and full will be filled
      ! A      a        j     A       a      i      x    x since j>i only full list will be filled
      !----------------------------------------------------
      ! A      a        i     A       b      j      x    x both half and full will be filled
      ! A      a        j     A       b      i      x    x only full will be filled
      !----------------------------------------------------
      ! A      b        i     A       a      j      x    x since b>a this will not go through
      !----------------------------------------------------
      ! A      a        i     B       a      i      x    x both half and full will be filled
      ! A      a        i     B       a      j      x    x both half and full will be filled
      ! A      a        j     B       a      i      x    x both half and full will be filled
      ! A      a        i     B       b      j      x    x both half and full will be filled
      ! A      a        j     B       b      i      x    x both half and full will be filled
      ! A      b        i     B       a      j      x    x both half and full will be filled
      !----------------------------------------------------
      !The best way to compare the VL is to print them for periodic crystals of different sizes
      IMPLICIT NONE
      INTEGER :: i,j,istart,iend,jstart,jend
      INTEGER :: position,atomindx1,atomindx2
      INTEGER :: atomnumberneigh
      REAL(dp) :: dr(3),dr2mag,drmag
      
      istart=(maincellindx-1)*MaxAtomPerCell+1
      iend=istart+LLListRange(maincellindx)-1

      jstart=(neighcellindx-1)*MaxAtomPerCell+1
      jend=jstart+LLListRange(neighcellindx)-1
      
      
      DO i=istart,iend
         
         !atomindx1=LL%List(i) !atom index from first cell
         atomindx1=LLList(i) !atom index from first cell
         
         DO j=jstart,jend

            !atomindx2=LL%List(j) !atom index from other cell
            atomindx2=LLList(j) !atom index from other cell

            !dr=PBCdistance(AL,atomindx2,atomindx1)
            dr=AtomCoord(3*atomindx2-2:3*atomindx2)-AtomCoord(3*atomindx1-2:3*atomindx1)
            !this defn of dr is required in the force calc. 
            !to get the derivative dr/dx
            !WHERE (.NOT. LargePotentialCutOff) dr=dr-BoxSize*NINT(dr/BoxSize) !nearest image spacing using periodic BC
            dr=dr-BoxSize*NINT(dr/BoxSize)*AllowPBC !nearest image spacing using periodic BC, AllowPBC will be 0 if LargePotentialCutOff and 1 otherwise
            dr=dr+REAL(BoxPosition)*BoxSize !gets spacing between the two boxes
            dr2mag=DOT_PRODUCT(dr,dr)
            
            !dr2mag can be 0 when maincellindx==neighcellindx and atom is overlapping
            IF (dr2mag<=CutOff2 .AND. dr2mag>0.0001_dp) THEN

               drmag=SQRT(dr2mag)
               
               !add the atom only if full list or when atomindx1<atomindx2 with half list
               
               IF (.NOT. CreateHalfList .OR. &
                 (CreateHalfList .AND. (atomindx1<atomindx2 .OR. CellOutsideBox))) THEN
                  atomnumberneigh=VLListRange(atomindx1)+1 !VL%ListRange()=0 to begin with

                  IF (atomnumberneigh>MaxAtomPerAtom) THEN
                     WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                     WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                     WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                     STOP
                  END IF

                  position=(atomindx1-1)*MaxAtomPerAtom+atomnumberneigh
                  VLList(position)=atomindx2 !add atom2
                  VLdr(3*position-2:3*position)=dr
                  VLdrmag(position)=drmag
                  VLListRange(atomindx1)=atomnumberneigh 
                     !update count in VL for atom 1
                  IF (CellOutsideBox) VLListDomainAtom(position)=.FALSE. !atomindx2 does not belong to the domain
               END IF
               
               !fill this for other neighbor
               IF (maincellindx/=neighcellindx .AND. .NOT. CellOutsideBox) THEN  !avoid double
                  IF (.NOT. CreateHalfList .OR. (CreateHalfList .AND. atomindx2<atomindx1)) THEN

                     !atomnumberneigh=VL%ListRange(atomindx2)+1
                     atomnumberneigh=VLListRange(atomindx2)+1

                     IF (atomnumberneigh>MaxAtomPerAtom) THEN
                        WRITE(*,*) "$Err>> Max atoms per atom in VL is small"
                        WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                        WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                        STOP
                     END IF

                     position=(atomindx2-1)*MaxAtomPerAtom+atomnumberneigh
                     VLList(position)=atomindx1 !add atom2
                     VLdr(3*position-2:3*position)=-dr
                     VLdrmag(position)=drmag
                     VLListRange(atomindx2)=atomnumberneigh
                  END IF
               END IF

            END IF

         END DO

      END DO

   END SUBROUTINE VLAdd
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      
   END SUBROUTINE VLFill

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE UpdateOldVLFull(AL,VL,NAtoms)
   !small changes would have occurred in the atom coordinates that 
   !need to be incorporated in the VL

      IMPLICIT NONE
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: i,j,r1,r2,indx,MaxAtomPerAtom,NAtoms
      REAL(dp) :: dr(3),drmag,dr2mag
      INTEGER, DIMENSION(:), POINTER :: VLListRange,VLList
      REAL(dp), DIMENSION(:), POINTER :: VLdr,VLdrmag
      
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      VLList=>VL%List
      VLListRange=>VL%ListRange
      VLdr=>VL%dr
      VLdrmag=>VL%drmag

      DO i=1,NAtoms

         r1=(i-1)*MaxAtomPerAtom+1
         r2=r1-1+VLListRange(i)

         DO indx=r1,r2

            !j=VL%List(indx)
            j=VLList(indx)
            !dr=PBCdistance(AL,j,i)
            dr=AtomCoord(3*j-2:3*j)-AtomCoord(3*i-2:3*i)
            dr=dr-BoxSize*NINT(dr/BoxSize)
            dr2mag=DOT_PRODUCT(dr,dr)
            drmag=SQRT(dr2mag)
            !VL%dr(3*indx-2:3*indx)=dr
            !VL%drmag(indx)=drmag
            VLdr(3*indx-2:3*indx)=dr
            VLdrmag(indx)=drmag

         END DO

      END DO

   END SUBROUTINE UpdateOldVLFull
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE UpdateOldVLLocal(AL,NAtoms,iatom,errorstatus)
   !used when only handful atoms are displaced and their neighborlist can be updated locally
   !assumes that the box size is large enough so that multiple images of neighbors dont show-up
      IMPLICIT NONE
      INTEGER, PARAMETER :: maxneigh=300 !maximum number of neighors
      INTEGER :: iatom,NAtoms,nx,ny,nz,c1(3),c2(3),NCell(3),cellindx,neighcellindx
      INTEGER :: ix,iy,iz,jatom,j,r0,MaxAtomPerCell,MaxAtomPerAtom,r1,r2
      TYPE(SystemContainer), POINTER :: AL
      TYPE(LinkedListContainer), POINTER :: LL
      TYPE(VerletListContainer), POINTER :: VL
      INTEGER, DIMENSION(:), POINTER :: List,ListRange
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,VLdr,VLdrmag
      REAL(dp) :: CellSize(3),CutOff,dr(3),drmag,BoxSize(3)
      INTEGER, DIMENSION(maxneigh) :: NeighborIndex
      REAL(dp), DIMENSION(maxneigh) :: Neighbordrmag
      REAL(dp), DIMENSION(3*maxneigh) :: Neighbordr
      INTEGER, DIMENSION(:), POINTER :: VLList,VLListRange
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom
      LOGICAL :: checking
      INTEGER :: ineigh,errorstatus !,nNeigh
      
      WRITE(6,*) "THIS SUBROUTINE UpdateOldVLLocal HAS SOME ISSUES. IT IS NOT MATCHING WITH ADDVL"
      STOP
      
      !initialize
      LL=>AL%LL
      VL=>AL%VL
      VLList=>VL%List
      VLListRange=>VL%ListRange
      VLListDomainAtom=>VL%ListDomainAtom
      VLdr=>VL%dr
      VLdrmag=>VL%drmag
      
      !check if variables are well defined
      IF (.NOT. ASSOCIATED(LL) .OR. .NOT. ASSOCIATED(VL)) THEN
         WRITE(6,*) "Err>> Linked and/or verlet list is not associated in UpdateVLLocally"
         STOP
      END IF
      
      !initialize and obtain the potential information
      CellSize=LL%CellSize
      NCell=LL%NumberCell3D
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      CutOff=VL%CutOff+VL%Buffer  !cutoff is a local variable in this subroutine
      nx=CEILING(CutOff/CellSize(1)) !number of cells to be included
      ny=CEILING(CutOff/CellSize(2)) !number of cells to be included
      nz=CEILING(CutOff/CellSize(3)) !number of cells to be included
      
      ListRange=>LL%ListRange
      List=>LL%List
      AtomCoord=>AL%AtomCoord
      MaxAtomPerCell=LL%MaxAtomPerCell
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      
      !move iatom to the correct cell location in the linked list
      c1=CEILING((AtomCoord(3*iatom-2:3*iatom)-BoxSize)/CellSize)
      c1=MAX(c1,1)
      cellindx=((c1(1)-1)*NCell(2)+c1(2)-1)*NCell(3)+c1(3)
      
      !WE ASSUME THAT THE LINKED LIST IS ROUGHLY CORRECT -- THIS IS BECAUSE WE
      !ARE USING THE OLD LINKEDLIST WITHOUT UPDATING IT
      
      !find all atom neighbors and make their list
      !assumes that the box size is large enough so that multiple images of neighbors dont show-up
      !check if all atoms present are within domain
      
      r1=(iatom-1)*MaxAtomPerAtom
      errorstatus=0
      IF (ANY(.NOT. VLListDomainAtom(r1+1:r1+VLListRange(iatom)))) THEN !atom lies outside the domain
         errorstatus=1
         RETURN
      END IF
      
      !here we assume that we are working with a half-list -- this makes the code bit more complex &
      ! we are forced to use linked lists
      !WARNING: Some of the atoms might be missing from the VL but these atoms will be close to cutoff
      !nNeigh=0
      DO ix=-nx,nx
         DO iy=-ny,ny
            DO iz=-nz,nz
               c2=c1+(/ix,iy,iz/)
               c2=c2-NCell*FLOOR(REAL(c2-1)/REAL(NCell))
               neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3)
               r0=(neighcellindx-1)*MaxAtomPerCell !contents of cells
               DO j=r0+1,r0+ListRange(neighcellindx)
                  jatom=List(j)
                  IF (iatom==jatom) CYCLE
                  dr=AtomCoord(3*jatom-2:3*jatom)-AtomCoord(3*iatom-2:3*iatom)
                  dr=dr-BoxSize*NINT(dr/BoxSize)
                  drmag=SQRT(dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3))
                  IF (drmag<=cutoff) CALL UpdateOldVLLocalDriver()
               END DO
            END DO
         END DO
      END DO
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE UpdateOldVLLocalDriver()
         IMPLICIT NONE
         INTEGER :: r1,i1
         
         !update for iatoms list -- note VL can half/full list
         r1=(iatom-1)*MaxAtomPerAtom
         DO i1=r1+1,r1+VLListRange(iatom)
            IF (VLList(i1)==jatom) THEN
               VLdr(3*i1-2:3*i1)=dr
               VLdrmag(i1)=drmag
            END IF
         END DO
         
         !update for jatoms list -- note VL can half/full list
         r1=(jatom-1)*MaxAtomPerAtom
         DO i1=r1+1,r1+VLListRange(iatom)
            IF (VLList(i1)==iatom) THEN
               VLdr(3*i1-2:3*i1)=-dr
               VLdrmag(i1)=drmag
            END IF
         END DO
      END SUBROUTINE UpdateOldVLLocalDriver
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE UpdateOldVLLocal
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE NeighborList
