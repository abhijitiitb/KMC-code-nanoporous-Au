MODULE SpatialAveragesPackage
   USE VARIABLE_TYPE
   USE nrprocedures
   IMPLICIT NONE 
   CONTAINS
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetSpeciesDistribution(AL,NDomains,ispecies,Nxyz,Nxx,Nxy,Nxz,Nyx,Nyy,Nyz,Nzx,Nzy,Nzz)
   !find the composition in a multicomponent system
   !Nxyz - find composition in (x,y,z) format, i.e., domains in 3D
   !Nxx - find composition averaged over y-z plane in x-direction
   !Nxy - find composition averaged over y-z plane in x-direction
   !Nxz - find composition averaged over y-z plane in x-direction
   !Nyx - find composition averaged over y-z plane in y-direction
   !Nyy - find composition averaged over y-z plane in y-direction
   !Nyz - find composition averaged over y-z plane in y-direction
   !Nzx - find composition averaged over y-z plane in z-direction
   !Nzy - find composition averaged over y-z plane in z-direction
   !Nzz - find composition averaged over y-z plane in z-direction
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(LinkedListContainer), POINTER :: LL
      INTEGER, DIMENSION(:), POINTER :: List,ListRange,AtomSpecies
      REAL(dp), DIMENSION(3) :: BoxSize
      INTEGER :: ispecies !species for which you need the composition
      INTEGER :: NDomains(3),ix,iy,iz,LLFactor(3),MaxAtomPerCell
      INTEGER :: MaxNcellx,MaxNcelly,MaxNcellz,icellx,icelly,icellz
      INTEGER :: Ncellx1,Ncellx2,Ncelly1,Ncelly2,Ncellz1,Ncellz2
      INTEGER :: Ncellx,Ncelly,Ncellz,r1,r2,j,CellIndex,jatom
      INTEGER, DIMENSION(:,:,:), POINTER :: Nxyz !count for number of atoms of each species type present in my domain
      INTEGER, DIMENSION(:), POINTER, OPTIONAL :: Nxx,Nxy,Nxz,Nyx,Nyy,Nyz,Nzx,Nzy,Nzz
      
      !Check whether LL is present
      IF (.NOT. ASSOCIATED(AL%LL)) THEN
         WRITE(6,*) "$Err>> Linked list should have been associated in GetConcentrationProfile"
         STOP
      END IF
      
      LL=>AL%LL
      ListRange=>LL%ListRange
      List=>LL%List
      AtomSpecies=>AL%AtomSpecies
      
      MaxNcellx=LL%NumberCell3D(1)
      MaxNcelly=LL%NumberCell3D(2)
      MaxNcellz=LL%NumberCell3D(3)
      
      Ncellx=MaxNcellx/NDomains(1)
      Ncelly=MaxNcelly/NDomains(2)
      Ncellz=MaxNcellz/NDomains(3)
      
      !get the cell index
      LLFactor(1)=LL%NumberCell3D(2)*LL%NumberCell3D(3) !x
      LLFactor(2)=LL%NumberCell3D(3) !y
      LLFactor(3)=1 !z
      MaxAtomPerCell=AL%LL%MaxAtomPerCell
      
      !find Nxyz
      DO ix=1,NDomains(1)
         DO iy=1,NDomains(2)
            DO iz=1,NDomains(3)
                
               Nxyz(ix,iy,iz)=0
               Ncellx1=(ix-1)*Ncellx+1
               Ncellx2=MIN(ix*Ncellx,MaxNcellx)
               Ncelly1=(iy-1)*Ncelly+1
               Ncelly2=MIN(iy*Ncelly,MaxNcelly)
               Ncellz1=(iz-1)*Ncellz+1
               Ncellz2=MIN(iz*Ncellz,MaxNcellz)
                
               DO icellx=Ncellx1,Ncellx2
                  DO icelly=Ncelly1,Ncelly2
                     DO icellz=Ncellz1,Ncellz2
                        CellIndex=(icellx-1)*LLFactor(1)+(icelly-1)*LLFactor(2)+icellz
                       
                        r1=(CellIndex-1)*MaxAtomPerCell+1
                        r2=(CellIndex-1)*MaxAtomPerCell+ListRange(cellIndex)
                        DO j=r1,r2
                           jatom=List(j)    ! this is the atom index
                           IF (AtomSpecies(jatom)==ispecies) THEN !from atom index atomspecies are getting
                              Nxyz(ix,iy,iz)=Nxyz(ix,iy,iz)+1
                           ENDIF
                        END DO
                  
                     END DO
                  END DO
               END DO
               
            END DO   
         END DO
      END DO  
      
      IF (PRESENT(Nxx)) THEN
         
      END IF
      
      IF (PRESENT(Nxy)) THEN
      END IF
      
      IF (PRESENT(Nxz)) THEN
      END IF
      
      IF (PRESENT(Nyx)) THEN
      END IF
      
      IF (PRESENT(Nyy)) THEN
      END IF
      
      IF (PRESENT(Nyz)) THEN
      END IF
      
      IF (PRESENT(Nzx)) THEN
      END IF
      
      IF (PRESENT(Nzy)) THEN
      END IF
      
      IF (PRESENT(Nzz)) THEN
      END IF
      
    END SUBROUTINE GetSpeciesDistribution
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    SUBROUTINE CalculateStrain(AL)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      
      IF (.NOT. ASSOCIATED(AL%VL)) THEN !need VL for the strain calculation
      END IF
    END SUBROUTINE CalculateStrain
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    SUBROUTINE CalculateAtomStrain(AL,iatom,strain,composition)
    !this has to be created by the user for the specific system
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: iatom
      REAL(dp), DIMENSION(3,3) :: strain
      REAL(dp), DIMENSION(:) :: composition
      REAL(dp) :: LatticeConstant,l
      REAL(dp), DIMENSION(3), TARGET :: rt1,rt2,rt3,rt4 !for Si-Ge we have four unstrained vectors
      REAL(dp), DIMENSION(3) :: r1p,r2p,r3p,r4p !strained system vectors
      REAL(dp), DIMENSION(:), POINTER :: r1,r2,r3,r4 !unstrained vectors
      
      !Calculate lattice constant
      LatticeConstant= Vegard(composition) !call this for unstrained case
      l=4._dp*LatticeConstant/SQRT(3._dp) !interatomic spacing in the unstrained case
      
      
      !Find the NN vectors in the unstrained case
      CALL UnstrainedSystemVectors(l,rt1,rt2,rt3,rt4)
         
      !Find the four nearest neighbors in the actual system (strained case)
      CALL StrainedSystemVectors(AL%VL,iatom,r1p,r2p,r3p,r4p)
      
      !Perform rotation to reorient the unstrained vectors
      CALL ReassignUnstrainedVectors(rt1,rt2,rt3,rt4,r1p,r2p,r3p,r4p,r1,r2,r3,r4)
      
      !Choose three vectors and create two matrices X' and X s.t. X'=E*X and E is the strain matrix
      !CALL
      
      !Compute the strain E in the system
      !CALL
         
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      FUNCTION Vegard(composition)
         IMPLICIT NONE
         REAL(dp), DIMENSION(:) :: composition
         REAL(dp) :: Vegard,x
         REAL, PARAMETER:: a_Si=5.43_dp
         REAL, PARAMETER:: a_Ge=5.64_dp
            
         x=composition(1)
         Vegard=x*a_Si+(1-x)*a_Ge
      END FUNCTION Vegard
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE CalculateAtomStrain
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE UnstrainedSystemVectors(l,r1,r2,r3,r4)
      IMPLICIT NONE
      REAL(dp) :: l
      REAL(dp), DIMENSION(:) ::r1,r2,r3,r4 ! also can give input vectors as allocatable
      
      r1=(/-l/4._dp,-l/4._dp,-l/4._dp/)
      r2=(/-l/4._dp,l/4._dp,l/4._dp/)
      r3=(/l/4._dp,-l/4._dp,l/4._dp/)
      r4=(/l/4._dp,l/4._dp,-l/4._dp/)
      
   END SUBROUTINE UnstrainedSystemVectors
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE StrainedSystemVectors(VL,iatom,r1p,r2p,r3p,r4p)
      IMPLICIT NONE
      INTEGER, PARAMETER :: MaxNNeigh=100 !max 100 neighbors possible
      TYPE(VerletListContainer), POINTER :: VL
      INTEGER :: iatom,index(MaxNNeigh),MaxAtomPerAtom
      INTEGER :: ra,r0,r1,r2,j
      REAL(dp) :: drmag(MaxNNeigh),r1p(3),r2p(3),r3p(3),r4p(3)
      
      !find the nearest neighbors using the verlet list (we are assuming that the VL has been created)
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      ra=(iatom-1)*MaxAtomPerAtom
      r0=VL%ListRange(iatom)
      
      r1=ra+1
      r2=ra+r0
      
      DO j=r1,r2
         index(j)=j
      END DO
      drmag(1:r0)=VL%drmag(r1:r2)
      
      CALL sort_heap_abhijit(drmag(1:r0),index(1:r0))
      r1p=VL%dr(3*index(1)-2:3*index(1))
      r2p=VL%dr(3*index(2)-2:3*index(2))
      r3p=VL%dr(3*index(3)-2:3*index(3))
      r4p=VL%dr(3*index(4)-2:3*index(4))
      
   END SUBROUTINE StrainedSystemVectors
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReassignUnstrainedVectors(rt1,rt2,rt3,rt4,r1p,r2p,r3p,r4p,r1,r2,r3,r4)
   !reassign the vectors belonging to the unstrained case so that r1 is closest to r1p,
   !r2 is closest to r2p and so on
      IMPLICIT NONE
      REAL(dp), DIMENSION(3), TARGET :: rt1,rt2,rt3,rt4,r1p,r2p,r3p,r4p
      REAL(dp), DIMENSION(:), POINTER :: r1,r2,r3,r4
      REAL(dp) :: drmag(4,4),drmin,spacing(24)
      INTEGER :: i,imin
      LOGICAL :: IsRemaining(4) !tells us whether it is remaining to be assigned
      
      IsRemaining=.TRUE.
      
      !following arrangements are possible --
      !r1, r2, r3 and r4 are presently pointers not pointing to anyone
      !rt1,rt2,rt3,rt4,r1p,r2p,r3p,r4p are vectors already defined
      !we wish to make r1 point to one of rt1,rt2,rt3,rt4 (simiilarly r2, r3 and r4 to point)
      !s.t. the r1 (r2,r3 and r4) so that it is closest to r1p (r2p,r3p,r4p)
      
      !to do this we create a matrix
      !drmag=
      !  _                                                    _
      ! |  |rt1-r1p]^2  |rt2-r1p]^2  |rt3-r1p]^2  |rt4-r1p]^2  |
      ! |  |rt1-r2p]^2  |rt2-r2p]^2  |rt3-r2p]^2  |rt4-r2p]^2  |
      ! |  |rt1-r3p]^2  |rt2-r3p]^2  |rt3-r3p]^2  |rt4-r3p]^2  |
      ! |_ |rt1-r4p]^2  |rt2-r4p]^2  |rt3-r4p]^2  |rt4-r4p]^2 _|
      !
      drmag(1,1)=GetSpacing2(rt1,r1p)
      drmag(2,1)=GetSpacing2(rt1,r2p)
      drmag(3,1)=GetSpacing2(rt1,r3p)
      drmag(4,1)=GetSpacing2(rt1,r4p)
      
      drmag(1,2)=GetSpacing2(rt2,r1p)
      drmag(2,2)=GetSpacing2(rt2,r2p)
      drmag(3,2)=GetSpacing2(rt2,r3p)
      drmag(4,2)=GetSpacing2(rt2,r4p)
      
      drmag(1,3)=GetSpacing2(rt3,r1p)
      drmag(2,3)=GetSpacing2(rt3,r2p)
      drmag(3,3)=GetSpacing2(rt3,r3p)
      drmag(4,3)=GetSpacing2(rt3,r4p)
      
      drmag(1,4)=GetSpacing2(rt4,r1p)
      drmag(2,4)=GetSpacing2(rt4,r2p)
      drmag(3,4)=GetSpacing2(rt4,r3p)
      drmag(4,4)=GetSpacing2(rt4,r4p)
      
      !following permutations are possible for [r1, r2, r3, r4]
      !(lets denote rt1=1, rt2=2, rt3=3 and rt4=4) 
      !for e.g., [2 1 3 4] it means that rt2 will be closest to r1p, rt1 will closest r2p, ..
      spacing(1)=drmag(1,1)+drmag(2,2)+drmag(3,3)+drmag(4,4) ! [1 2 3 4]
      spacing(2)=drmag(1,1)+drmag(2,2)+drmag(3,4)+drmag(4,3) ! [1 2 4 3]
      spacing(3)=drmag(1,1)+drmag(2,3)+drmag(3,2)+drmag(4,4) ! [1 3 2 4]
      spacing(4)=drmag(1,1)+drmag(2,3)+drmag(3,4)+drmag(4,2) ! [1 3 4 2]
      spacing(5)=drmag(1,1)+drmag(2,4)+drmag(3,3)+drmag(4,2) ! [1 4 3 2]
      spacing(6)=drmag(1,1)+drmag(2,4)+drmag(3,2)+drmag(4,3) ! [1 4 2 3]
      
      spacing(7)=drmag(1,2)+drmag(2,1)+drmag(3,3)+drmag(4,4) ! [2 1 3 4]
      spacing(8)=drmag(1,2)+drmag(2,1)+drmag(3,4)+drmag(4,3)! [2 1 4 3]
      spacing(9)=drmag(1,2)+drmag(2,3)+drmag(3,1)+drmag(4,4)! [2 3 1 4]
      spacing(10)=drmag(1,2)+drmag(2,3)+drmag(3,4)+drmag(4,1)! [2 3 4 1]
      spacing(11)=drmag(1,2)+drmag(2,4)+drmag(3,1)+drmag(4,3) ! [2 4 1 3]
      spacing(12)=drmag(1,2)+drmag(2,4)+drmag(3,3)+drmag(4,1)! [2 4 3 1]
      
      spacing(13)=drmag(1,3)+drmag(2,1)+drmag(3,2)+drmag(4,4) ! [3 1 2 4]
      spacing(14)=drmag(1,3)+drmag(2,1)+drmag(3,4)+drmag(4,2) ! [3 1 4 2]
      spacing(15)=drmag(1,3)+drmag(2,2)+drmag(3,1)+drmag(4,4) ! [3 2 1 4]
      spacing(16)=drmag(1,3)+drmag(2,2)+drmag(3,4)+drmag(4,1) ! [3 2 4 1]
      spacing(17)=drmag(1,3)+drmag(2,4)+drmag(3,1)+drmag(4,2) ! [3 4 1 2]
      spacing(18)=drmag(1,3)+drmag(2,4)+drmag(3,2)+drmag(4,1) ! [3 4 2 1]
      
      spacing(19)=drmag(1,4)+drmag(2,1)+drmag(3,2)+drmag(4,3) ! [4 1 2 3]
      spacing(20)=drmag(1,4)+drmag(2,1)+drmag(3,3)+drmag(4,2) ! [4 1 3 2]
      spacing(21)=drmag(1,4)+drmag(2,2)+drmag(3,1)+drmag(4,3) ! [4 2 1 3]
      spacing(22)=drmag(1,4)+drmag(2,2)+drmag(3,3)+drmag(4,1) ! [4 2 3 1]
      spacing(23)=drmag(1,4)+drmag(2,3)+drmag(3,2)+drmag(4,1) ! [4 3 2 1]
      spacing(24)=drmag(1,4)+drmag(2,3)+drmag(3,1)+drmag(4,2) ! [4 3 1 2]
      
      
      
      !get r2 so that it is closest to r2p -------------------------------
      
   END SUBROUTINE ReassignUnstrainedVectors
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetSpacing2(v1,v2)
      IMPLICIT NONE
      REAL(dp) :: GetSpacing2,v1(3),v2(3),diff(3)
      
      diff=v1-v2
      GetSpacing2=SUM(diff*diff)
   END FUNCTION GetSpacing2
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE SpatialAveragesPackage
