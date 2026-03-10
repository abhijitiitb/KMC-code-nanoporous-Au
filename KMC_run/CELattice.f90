
MODULE CEModel_Lattice

   USE VARIABLE_TYPE
   USE db_manipulate
   USE utilities
   USE io
   USE PotentialPackage
   USE NeighborList
   USE OptimizationAL
   
   
   IMPLICIT NONE
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!contains routines used by mainoption3.f90 for finding the barrier within the lattice approximation
   SUBROUTINE LatticeSetupArrays(AL,Atom_Index,AtomicSymbol)
!based on the AL provided the species array has to be setup
   IMPLICIT NONE
!   REAL(dp), DIMENSION(:), POINTER :: AtomCoord
   TYPE(SystemContainer), POINTER :: AL
   CHARACTER(len=2),INTENT(IN):: AtomicSymbol
   INTEGER :: NAtoms,No_Atoms
   INTEGER,OPTIONAL ::Atom_Index
!   INTEGER, DIMENSION(:), POINTER :: species_index,LatticeSpeciesArray, Env_Atom_index,LatticeEnvArray
   INTEGER, DIMENSION(:), POINTER :: species_index,Env_Atom_index
   REAL(dp), DIMENSION(3) :: sclfac
 
   sclfac=(/4.0801_dp/4.089_dp,4.0801_dp/4.089_dp,4.0801_dp/4.089_dp/) !Au  !IMPORTANT

   NAtoms=AL%NAtoms 
   
   NULLIFY(species_index)
!   NULLIFY(LatticeSpeciesArray)
   
   ALLOCATE(species_index(NAtoms))
!   ALLOCATE(LatticeSpeciesArray(NAtoms))
   
   !LatticeSpeciesArray
   !Species indexing
   !1: Cu
   !2: Ag
   !3: Au
   !4: Ni
   !5: Pd
   !6: Pt

   SELECT CASE (AtomicSymbol(1:2))
   
   CASE("Cu")
    species_index(Atom_Index)=1
   CASE("Ag")
    species_index(Atom_Index)=2
   CASE("Au")
    species_index(Atom_Index)=3
   CASE("Ni")
    species_index(Atom_Index)=4
   CASE("Pd")
    species_index(Atom_Index)=5
   CASE("Pt")
    species_index(Atom_Index)=6
    
  END SELECT
  
  LatticeSpeciesArray(Atom_Index)=species_index(Atom_Index)
!  Write(6,*) species_index(Atom_Index)
!  WRITE(6,*)LatticeSpeciesArray(Atom_Index),Atom_Index


   !LatticeSpeciesArray
   !Species indexing
   !1: Adatom(100)
   !2: Surface(100)
   !3: Subsurface(100)
   !4: Bulk
   
  ALLOCATE(Env_Atom_index(NAtoms)) 
!  ALLOCATE(LatticeEnvArray(NAtoms))
  
  IF(AL%AtomCoord(3*Atom_Index)>=(24.219_dp*sclfac(3)+1.0_dp)) THEN
!	WRITE(6,*)"Atom is Adatom"
	Env_Atom_index(Atom_Index)=1
  ELSEIF(AL%AtomCoord(3*Atom_Index)==(22.4895_dp*sclfac(3)+1.0_dp)) THEN
!    WRITE(6,*)"Atom is Surface Atom"
    Env_Atom_index(Atom_Index)=2
  ELSEIF(AL%AtomCoord(3*Atom_Index)==(20.445_dp*sclfac(3)+1.0_dp)) THEN
!    WRITE(6,*)"Atom is Subsurface Atom"
    Env_Atom_index(Atom_Index)=3
  ELSEIF(AL%AtomCoord(3*Atom_Index)<(20.445_dp*sclfac(3)+1.0_dp)) THEN
!    WRITE(6,*)"Atom is Bulk Atom"
    Env_Atom_index(Atom_Index)=4
  END IF 
  
  LatticeEnvArray(Atom_Index)=Env_Atom_index(Atom_Index)
!  WRITE(6,*)LatticeEnvArray(Atom_Index),Atom_Index
  
  
END SUBROUTINE LatticeSetupArrays
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE LatticeSetupPotentials()
!reads tables containing potential information
END SUBROUTINE LatticeSetupPotentials
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE LatticeSelectAtoms(AL,nmoving,rcut)
!from 1:nmoving are the atoms moving
!atoms within rcut distance are selected for activation barrier calculation
   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: AL
   INTEGER :: nmoving,NAtoms,j
   REAL(dp) :: rcut
   LOGICAL :: LatticeAtomsSelected(200000)
   
   INTEGER, DIMENSION(:), POINTER :: AtomSpecies
   REAL(dp), DIMENSION(:), POINTER :: AtomCoord=>NULL(),AtomVelocity=>NULL()
   TYPE(VerletListContainer), POINTER :: VL=>NULL()
   REAL(dp), DIMENSION(:), POINTER :: VLdrmag=>NULL(),VLdr=>NULL()
   INTEGER, DIMENSION(:), POINTER :: VLListRange=>NULL(),VLList=>NULL()
   INTEGER :: ra=0,r0=0,r_1=0,r_2=0,iatom=0,jatom=0,MaxAtomPerAtom=0,katom=0,s_1,s_2
   REAL(dp) :: drmagij,drij(3)
   LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
   
   
   NAtoms=AL%NAtoms  
   AtomSpecies=>AL%AtomSpecies
   AtomCoord=>AL%AtomCoord  
   
   LatticeAtomsSelected=.FALSE.
  
   IF (ASSOCIATED(AL%AtomIsMoving)) THEN
    AtomIsMoving=>AL%AtomIsMoving
   END IF  

   CALL PotentialInitialize(AL)
   CALL AddVerletList(AL,ListType=1,NAtoms=AL%NAtoms)  
  
   VL=>AL%VL
   VLList=>VL%List
   VLListRange=>VL%ListRange
   MaxAtomPerAtom=VL%MaxAtomPerAtom     
   VLdrmag=>VL%drmag
   VLdr=>VL%dr
  
   DO iatom=1,nmoving+1
        ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
        r0=VLListRange(iatom) !range for ith atom
        r_1=ra+1
        r_2=ra+r0
        s_1=AtomSpecies(iatom)
          DO j=r_1,r_2 !this will take care of the pair and many-body terms in one go
             jatom=VLList(j)
             s_2=AtomSpecies(jatom)
             drij=VLdr(3*j-2:3*j)  !vector between atoms i and j
             drmagij=VLdrmag(j)   ! spacing between atoms i and j
          !   WRITE(6,*) iatom,"-",jatom,":",drmagij
          IF(drmagij<=rcut) THEN 
			LatticeAtomsSelected(jatom)=.TRUE.
	     	WRITE(6,*) "Match found" ,jatom
            WRITE(6,*) LatticeAtomsSelected(jatom)
          END IF
          END DO
  END DO 

END SUBROUTINE LatticeSelectAtoms
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE LatticeSetupReference(state) !initial state
   IMPLICIT NONE
   CHARACTER :: state
END SUBROUTINE LatticeSetupReference
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE LatticeEnergy(E)
   IMPLICIT NONE
   REAL(dp) :: E
END SUBROUTINE LatticeEnergy
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

END MODULE CEModel_Lattice
