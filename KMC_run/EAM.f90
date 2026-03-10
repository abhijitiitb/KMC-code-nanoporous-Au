!LL=Linked list
!VL=Verlet list
!AL=Atom list

MODULE EAMSubroutines

   USE VARIABLE_TYPE
   IMPLICIT NONE
   INTEGER, PARAMETER :: maxatom=200000
   REAL(dp), DIMENSION(maxatom) :: DensityArray,EnergyAtomArr
   REAL(dp), DIMENSION(3*maxatom) :: DensityArray1,DensityArray2,DensityArray3
   LOGICAL :: CatchOutOfBound

   CONTAINS
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  Only energy evaluation xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EAMEnTaylor(AL,NAtoms) !EAM Energy Generic (n species) Taylor
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(EAMPotential), POINTER :: EAM
      INTEGER, DIMENSION(:), POINTER :: SizePP,SizeDT,SizeEE,VLListRange,VLList,AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom
      INTEGER :: ra,r0,r1,r2,iatom,jatom,j,section,NSpecies,NAtoms
      INTEGER :: s1,s2,smin,smax,pos,k,MaxAtomPerAtom
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag,VLdr,AtomCoord
      REAL(dp), DIMENSION(:), POINTER :: PP,DT,EE,RangePP,RangeDT,RangeEE,IncrPPinv,IncrDTinv,IncrEEinv
      REAL(dp) :: drmag,tmp,p,mincutoff,Density,dxinv
      REAL(dp) :: EnergyEmbedding,EnergyPP,EnergyPPTerm

      !NAtoms=AL%NAtoms
      NSpecies=NSpecies_global
      VL=>AL%VL
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      EAM=>AL%Potential%EAM

      EnergyPP=0._dp
      EnergyEmbedding=0._dp

      AtomSpecies=>AL%AtomSpecies
      AtomCoord=>AL%AtomCoord
      VLListRange=>VL%ListRange
      VLListDomainAtom=>VL%ListDomainAtom
      VLList=>VL%List
      VLdrmag=>VL%drmag
      VLdr=>VL%dr
      PP=>EAM%PP
      DT=>EAM%DT
      EE=>EAM%EE
      SizePP=>EAM%SizePP
      SizeDT=>EAM%SizeDT
      SizeEE=>EAM%SizeEE
      RangePP=>EAM%RangePP
      RangeDT=>EAM%RangeDT
      RangeEE=>EAM%RangeEE
      IncrPPinv=>EAM%IncrPPinv
      IncrDTinv=>EAM%IncrDTinv
      IncrEEinv=>EAM%IncrEEinv

      !AL%AtomForce=0._dp !never do this
      
      IF (NAtoms>maxatom) THEN
         WRITE(6,*) "$Err>> Increase value of maxatom in EAM"
         STOP
      END IF
      
      IF (AL%VL%ListType==2) THEN !use half-list
         DensityArray(1:NAtoms)=0._dp
        
         DO iatom=1,NAtoms !first deal with all pair parts
            INCLUDE "EAMEnHalfTaylorStep1.f90"
         END DO
         
         DO iatom=1,NAtoms !next deal with terms that use the pair parts
            INCLUDE "EAMEnHalfTaylorStep2.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergyEmbedding+EnergyPP

      ELSEIF (AL%VL%ListType==1) THEN !use full list

         DO iatom=1,NAtoms
            INCLUDE "EAMEnTaylorFull.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergyEmbedding+0.5_dp*EnergyPP

      ELSE

         WRITE(6,*) "$Err>> VL type could not be recognized"
         STOP

      END IF

   END SUBROUTINE EAMEnTaylor
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Energy + Force evaluation xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EAMForceLinearInterp(AL,NAtoms,errorstatus)
   !errorstatus=(0, no errors), (1, spacing less than lower PP cutoff), (2, spacing less than lower PP cutoff)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(EAMPotential), POINTER :: EAM
      INTEGER, DIMENSION(:), POINTER :: SizePP,SizeDT,SizeEE,VLListRange,VLList,AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom,AtomIsMoving
      INTEGER :: ra,r0,r1,r2,iatom,jatom,j,section,NSpecies,NAtoms,errorstatus
      INTEGER :: s1,s2,smin,smax,pos,k,k0,MaxAtomPerAtom
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag,VLdr,AtomCoord,AtomForce
      REAL(dp), DIMENSION(:), POINTER :: PP,DT,EE,RangePP,RangeDT,RangeEE,IncrPPinv,IncrDTinv,IncrEEinv
      REAL(dp) :: drmag,tmp,p,gi,gi1,slope,tgphi,tgrho,tgrad,grad(3),ro1,mincutoff,Density,dxinv
      REAL(dp) :: EnergyEmbedding,EnergyPP,EnergyPPTerm

      !NAtoms=AL%NAtoms
      NSpecies=NSpecies_global
      VL=>AL%VL
      IF (.NOT. ASSOCIATED(AL%VL%List)) THEN
      !system size must be quite large so we have transitioned to a LL formalism
         CALL EAMForceLinearInterpLL(AL,NAtoms,errorstatus)
         RETURN
      END IF
      
      CatchOutOfBound=.FALSE.
      
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      EAM=>AL%Potential%EAM

      EnergyPP=0._dp
      EnergyEmbedding=0._dp

      AtomSpecies=>AL%AtomSpecies
      AtomCoord=>AL%AtomCoord
      AtomForce=>AL%AtomForce
      AtomIsMoving=>AL%AtomIsMoving
      VLListRange=>VL%ListRange
      VLListDomainAtom=>VL%ListDomainAtom
      VLList=>VL%List
      VLdrmag=>VL%drmag
      VLdr=>VL%dr
      PP=>EAM%PP
      DT=>EAM%DT
      EE=>EAM%EE
      SizePP=>EAM%SizePP
      SizeDT=>EAM%SizeDT
      SizeEE=>EAM%SizeEE
      RangePP=>EAM%RangePP
      RangeDT=>EAM%RangeDT
      RangeEE=>EAM%RangeEE
      IncrPPinv=>EAM%IncrPPinv
      IncrDTinv=>EAM%IncrDTinv
      IncrEEinv=>EAM%IncrEEinv
      
      IF (NAtoms>maxatom) THEN
         WRITE(6,*) "$Err>> Increase value of maxatom in EAM"
         STOP
      END IF

      errorstatus=0
      IF (VL%ListType==2) THEN !use half-list

         DensityArray(1:NAtoms)=0._dp
         DensityArray1(1:3*NAtoms)=0._dp
      
         DO iatom=1,NAtoms
            INCLUDE "EAMForceHalfLinearInterpStep1.f90"
         END DO

         DO iatom=1,NAtoms 

            INCLUDE "EAMForceHalfLinearInterpStep2.f90"
                       
!           If(iatom==1) WRITE(6,*) DensityArray(iatom),iatom
!           If(iatom==1) WRITE(6,*) DensityArray1(iatom),iatom  

         Atom_Density(iatom)=DensityArray(iatom)
         END DO

         DO iatom=1,NAtoms
            INCLUDE "EAMForceHalfLinearInterpStep3.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergyEmbedding+EnergyPP

      ELSEIF (VL%ListType==1) THEN !use full list

         DO iatom=1,NAtoms
            INCLUDE "EAMForceFullLinearInterp.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergyEmbedding+0.5_dp*EnergyPP

      ELSE

         WRITE(6,*) "$Err>> VL type could not be recognized"
         STOP

      END IF
      
   END SUBROUTINE EAMForceLinearInterp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EAMForceLinearInterpLL(AL,NAtoms,errorstatus)
   !this subroutine employ LL neighborlist - typically used when system size is large
   !so that VL generation is expensive, especially if this is done without spatial 
   !parallelization
      IMPLICIT NONE
      INTEGER :: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(LinkedListContainer), POINTER :: LL
      TYPE(EAMPotential), POINTER :: EAM
      INTEGER, DIMENSION(:), POINTER :: SizePP,SizeDT,SizeEE,AtomSpecies,LLList,LLListRange
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      INTEGER :: ra,r0,r1,r2,iatom,jatom,j,section,NSpecies,NAtoms,errorstatus
      INTEGER :: c1(3),c2(3),NCell(3),maincellindx,neighcellindx,nx,ny,nz,j1,j2,j3
      INTEGER :: istart,iend,jstart,jend,tid
      INTEGER :: s1,s2,smin,smax,pos,k,k0,MaxAtomPerCell
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomForce
      REAL(dp), DIMENSION(:), POINTER :: PP,DT,EE,RangePP,RangeDT,RangeEE,IncrPPinv,IncrDTinv,IncrEEinv
      REAL(dp) :: CellSize(3),xc(3),dr(3),drmag,tmp,p,gi,gi1,slope,tgphi,tgrho,tgrad,grad(3)
      REAL(dp) :: ro1,mincutoff,Density,dxinv,cutoff
      REAL(dp) :: EnergyEmbedding,EnergyPP,EnergyPPTerm,BoxSize(3),PotentialEnergy
      INTEGER, PARAMETER :: chunksize=1000 !take 1000 atom at a time --
      !if there are million atoms and 100 processor - then each thread is called approximately 10 times
      
      !NAtoms=AL%NAtoms
      NSpecies=NSpecies_global
      NAtoms=AL%NAtoms
      VL=>AL%VL
      LL=>AL%LL
      EAM=>AL%Potential%EAM
      cutoff=VL%CutOff+VL%Buffer

      EnergyPP=0._dp
      EnergyEmbedding=0._dp
      
      LLList=>LL%List
      LLListRange=>LL%ListRange
      MaxAtomPerCell=LL%MaxAtomPerCell

      AtomSpecies=>AL%AtomSpecies
      AtomCoord=>AL%AtomCoord
      AtomForce=>AL%AtomForce
      AtomIsMoving=>AL%AtomIsMoving
      PP=>EAM%PP
      DT=>EAM%DT
      EE=>EAM%EE
      SizePP=>EAM%SizePP
      SizeDT=>EAM%SizeDT
      SizeEE=>EAM%SizeEE
      RangePP=>EAM%RangePP
      RangeDT=>EAM%RangeDT
      RangeEE=>EAM%RangeEE
      IncrPPinv=>EAM%IncrPPinv
      IncrDTinv=>EAM%IncrDTinv
      IncrEEinv=>EAM%IncrEEinv
      
      errorstatus=0
      
      CatchOutOfBound=.FALSE.

      nx=CEILING(VL%CutOff/LL%CellSize(1))
      ny=CEILING(VL%CutOff/LL%CellSize(2))
      nz=CEILING(VL%CutOff/LL%CellSize(3))
      NCell=LL%NumberCell3D
      
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      xc=AL%BoxSize(4:6)
      CellSize=LL%CellSize

      DensityArray(1:NAtoms)=0._dp
      DensityArray1(1:3*NAtoms)=0._dp
      
      EnergyAtomArr(1:NAtoms)=0._dp
      
      IF (NAtoms>maxatom) THEN
         WRITE(6,*) "$Err>> Increase value of maxatom in EAM"
         STOP
      END IF
!OMP parallelization at level of iatom 
!iatom is PRIVATE for each thread - so each thread maintains its values

!$OMP PARALLEL 
!!!$OMP& PRIVATE(EnergyAtom)

!-------------------------------------------------------------------------
!$OMP DO SCHEDULE(DYNAMIC,chunksize)
      DO iatom=1,NAtoms !this has been partially borrowed from VLFill
         
         CALL Step1(iatom,AtomCoord,AtomSpecies,LLList,LLListRange,xc,CellSize,BoxSize, &
            nx,ny,nz,NCell,MaxAtomPerCell,cutoff,IncrDTinv, &
            DT,RangeDT,SizeDT,errorstatus)
         
      END DO
!$OMP END DO NOWAIT

!$OMP BARRIER
!-------------------------------------------------------------------------

!$OMP DO SCHEDULE(DYNAMIC,chunksize)
      DO iatom=1,NAtoms !this has been partially borrowed from VLFill
         
         CALL Step2(iatom,AtomSpecies,EE,IncrEEinv,RangeEE,SizeEE,errorstatus)
         
      END DO
!$OMP END DO NOWAIT

!$OMP BARRIER
!-------------------------------------------------------------------------

!$OMP DO SCHEDULE(DYNAMIC,chunksize)

      DO iatom=1,NAtoms !this has been partially borrowed from VLFill
         
         !EnergyAtom=0._dp !initialize
         CALL Step3(iatom,AtomCoord,AtomSpecies,AtomForce, &
            xc,CellSize,BoxSize,MaxAtomPerCell, &
            LLList,LLListRange,NCell,NSpecies,cutoff,DT,SizeDT,RangeDT,IncrDTinv, &
            PP,SizePP,RangePP,IncrPPinv,nx,ny,nz,errorstatus)
         
      END DO
   write(*,*) "oops bug here??? why no omp end do"
   STOP
!$OMP BARRIER
!-------------------------------------------------------------------------
!$OMP END PARALLEL
!-------------------------------------------------------------------------
PotentialEnergy=0._dp
!!$OMP PARALLEL DO
!!$OMP& DEFAULT(SHARED) PRIVATE(iatom)
!!$OMP& SCHEDULE(STATIC,chunksize)
!!$OMP& REDUCTION(+:,PotentialEnergy)
   DO iatom=1,NAtoms
      PotentialEnergy=PotentialEnergy+EnergyAtomArr(iatom)
   END DO
!!$OMP END PARALLEL DO
!-------------------------------------------------------------------------
      
      !DO iatom=1,NAtoms
      !   PotentialEnergy=PotentialEnergy+EnergyAtomArr(iatom)
      !END DO
      
      AL%PotentialEnergy=AL%PotentialEnergy+PotentialEnergy
      
   END SUBROUTINE EAMForceLinearInterpLL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Step1(iatom,AtomCoord,AtomSpecies,LLList,LLListRange,xc,CellSize,BoxSize, &
      nx,ny,nz,NCell,MaxAtomPerCell,cutoff,IncrDTinv, &
      DT,RangeDT,SizeDT,errorstatus)
   !performs step1 of EAMForceLinearInterpLL
      IMPLICIT NONE
      INTEGER, DIMENSION(:), POINTER :: LLList,LLListRange,AtomSpecies,SizeDT
      REAL(dp), DIMENSION(:), POINTER :: DT
      REAL(dp), DIMENSION(:), POINTER :: IncrDTinv,RangeDT
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord
      INTEGER :: iatom,c1(3),c2(3),s1,s2,maincellindx,neighcellindx,k0,k,errorstatus
      INTEGER :: j1,j2,j3,nx,ny,nz,NCell(3),jstart,jend,j,jatom,MaxAtomPerCell,section
      REAL(dp) :: CellSize(3),xc(3),dr(3),drmag,cutoff,tmp,mincutoff,dxinv,p,BoxSize(3)
      
      c1=CEILING((AtomCoord(3*iatom-2:3*iatom)-xc)/CellSize)
      c1=MAX(c1,1)

      !Convert to scalar number
      maincellindx=((c1(1)-1)*NCell(2) + (c1(2)-1))*NCell(3) + c1(3)
      s1=AtomSpecies(iatom)
         
      !loop over neighbor cells
      DO j3=-nz,nz
         DO j2=-ny,ny
            DO j1=-nx,nx
               
               !Step 1. Obtain the neighbor cell
               c2(1)=c1(1)+j1
               c2(2)=c1(2)+j2
               c2(3)=c1(3)+j3
               c2=c2-NCell*FLOOR(REAL(c2-1)/REAL(NCell)) !wrap the box using periodic BC
            
               neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3)
               jstart=(neighcellindx-1)*MaxAtomPerCell+1
               jend=jstart+LLListRange(neighcellindx)-1
               
               !Step 2. For each atom in maincell perform force calculation
               DO j=jstart,jend
                  jatom=LLList(j)
                  IF (iatom/=jatom) THEN
                     INCLUDE "EAMForceLLLinearInterpStep1.f90" !find the density from atom j to atom i
                  END IF
               END DO
                  
            END DO
         END DO
      END DO
   END SUBROUTINE Step1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Step2(iatom,AtomSpecies,EE,IncrEEinv,RangeEE,SizeEE,errorstatus)
      IMPLICIT NONE
      INTEGER :: iatom,s1,section,k0,k,errorstatus
      INTEGER, DIMENSION(:), POINTER :: SizeEE,AtomSpecies
      REAL(dp), DIMENSION(:), POINTER :: EE,IncrEEinv,RangeEE
      REAL(dp) :: dxinv,mincutoff,tmp,p,EnergyEmbedding,gi,gi1
      
      s1=AtomSpecies(iatom)
      EnergyEmbedding=0._dp
      INCLUDE "EAMForceLLLinearInterpStep2.f90"
      EnergyAtomArr(iatom)=EnergyAtomArr(iatom)+EnergyEmbedding
   END SUBROUTINE Step2
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Step3(iatom,AtomCoord,AtomSpecies,AtomForce, &
      xc,CellSize,BoxSize,MaxAtomPerCell, &
      LLList,LLListRange,NCell,NSpecies,cutoff,DT,SizeDT,RangeDT,IncrDTinv, &
      PP,SizePP,RangePP,IncrPPinv,nx,ny,nz,errorstatus)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomForce
      REAL(dp), DIMENSION(:), POINTER :: DT,RangeDT,IncrDTinv,PP,RangePP,IncrPPinv
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies,SizeDT,SizePP,LLList,LLListRange
      INTEGER :: iatom,jatom,c1(3),c2(3),maincellindx,neighcellindx,s1,s2,j1,j2,j3,nx,ny,nz
      INTEGER :: j,jstart,jend,k,k0,section,smin,smax
      INTEGER :: NSpecies,NCell(3),errorstatus,MaxAtomPerCell
      REAL(dp) :: p,drmag,dr(3),cutoff,dxinv,mincutoff,tmp,gi,gi1,slope
      REAL(dp) :: tgphi,tgrho,tgrad,grad(3),CellSize(3),BoxSize(3),xc(3),EnergyPP
      
      c1=CEILING((AtomCoord(3*iatom-2:3*iatom)-xc)/CellSize)
      c1=MAX(c1,1)

      !Convert to scalar number
      maincellindx=((c1(1)-1)*NCell(2) + (c1(2)-1))*NCell(3) + c1(3)
      s1=AtomSpecies(iatom)
      
      EnergyPP=0._dp
      
      !loop over neighbor cells
      DO j3=-nz,nz
         DO j2=-ny,ny
            DO j1=-nx,nx
            
               !Step 1. Obtain the neighbor cell
               c2=c1+(/j1,j2,j3/)
               c2=c2-NCell*FLOOR(REAL(c2-1)/REAL(NCell)) !wrap the box using periodic BC
               neighcellindx=((c2(1)-1)*NCell(2)+c2(2)-1)*NCell(3)+c2(3) !inefficient??
               jstart=(neighcellindx-1)*MaxAtomPerCell+1
               jend=jstart+LLListRange(neighcellindx)-1
               
               !Step 3. Compute the forces and pair terms
               DO j=jstart,jend
                  jatom=LLList(j)
                  IF (iatom/=jatom) THEN
                     INCLUDE "EAMForceLLLinearInterpStep3.f90" !find the force from pair i-j
                  END IF
               END DO
            END DO
         END DO
      END DO
      
      EnergyAtomArr(iatom)=EnergyAtomArr(iatom)+EnergyPP/2._dp
   END SUBROUTINE Step3
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE EAMSubroutines
