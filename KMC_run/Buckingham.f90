!LL=Linked list
!VL=Verlet list
!AL=Atom list

MODULE BuckinghamSubroutines

   USE VARIABLE_TYPE
   IMPLICIT NONE
   INTEGER, PARAMETER :: maxatom=20000
   
   CONTAINS
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  Only energy evaluation xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE BuckinghamEnTaylor(AL,NAtoms) !Buckingham Energy Generic (n species) Taylor
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(BuckinghamPotential), POINTER :: Buckingham
      INTEGER, DIMENSION(:), POINTER :: SizePP,VLListRange,VLList,AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom
      INTEGER :: ra,r0,r1,r2,iatom,jatom,j,section,NSpecies,NAtoms,MaxAtomPerAtom
      INTEGER :: s1,s2,smin,smax,pos,k,MaxPairPotentialTableSize
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag,VLdr,AtomCoord,AtomForce
      REAL(dp), DIMENSION(:), POINTER :: PP,RangePP,IncrPPinv
      REAL(dp) :: drmag,tmp,p
      REAL(dp) :: EnergyPP,mincutoff,dxinv,EnergyPPTerm

      !NAtoms=AL%NAtoms
      NSpecies=NSpecies_global
      MaxPairPotentialTableSize=MaxBuckTableSize
      
      VL=>AL%VL
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      Buckingham=>AL%Potential%Buckingham

      EnergyPP=0._dp

      AtomSpecies=>AL%AtomSpecies
      AtomCoord=>AL%AtomCoord
      AtomForce=>AL%AtomForce
      VLListRange=>VL%ListRange
      VLList=>VL%List
      VLdrmag=>VL%drmag
      VLListDomainAtom=>VL%ListDomainAtom
      VLdr=>VL%dr
      PP=>Buckingham%PP
      SizePP=>Buckingham%SizePP
      RangePP=>Buckingham%RangePP
      IncrPPinv=>Buckingham%IncrPPinv

      !AL%AtomForce=0._dp !never do this
      
      IF (AL%VL%ListType==2) THEN !use half-list
        
         DO iatom=1,NAtoms
            INCLUDE "PairPotentialEnTaylorHalf.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergyPP

      ELSEIF (AL%VL%ListType==1) THEN !use full list

         DO iatom=1,NAtoms
            INCLUDE "PairPotentialEnTaylorFull.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+0.5_dp*EnergyPP

      ELSE

         WRITE(6,*) "$Err>> VL type could not be recognized"
         STOP

      END IF

   END SUBROUTINE BuckinghamEnTaylor
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Energy + Force evaluation xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE BuckinghamForceLinearInterp(AL,NAtoms,errorstatus)
   !errorstatus=(0, no errors), (1, spacing less than lower PP cutoff)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(BuckinghamPotential), POINTER :: Buckingham
      INTEGER, DIMENSION(:), POINTER :: SizePP,VLListRange,VLList,AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom
      INTEGER :: ra,r0,r1,r2,iatom,jatom,j,section,NSpecies,NAtoms,MaxAtomPerAtom
      INTEGER :: s1,s2,smin,smax,pos,k,k0,MaxPairPotentialTableSize,errorstatus
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag,VLdr,AtomCoord,AtomForce
      REAL(dp), DIMENSION(:), POINTER :: PP,RangePP,IncrPPinv
      REAL(dp) :: drmag,tmp,p,gi,gi1,slope,tgphi,tgrho,tgrad,grad(3),ro1
      REAL(dp) :: EnergyPP,mincutoff,dxinv,EnergyPPTerm

      !NAtoms=AL%NAtoms
      NSpecies=NSpecies_global
      MaxPairPotentialTableSize=MaxBuckTableSize
      
      VL=>AL%VL
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      Buckingham=>AL%Potential%Buckingham

      EnergyPP=0._dp

      AtomSpecies=>AL%AtomSpecies
      AtomCoord=>AL%AtomCoord
      AtomForce=>AL%AtomForce
      VLListRange=>VL%ListRange
      VLList=>VL%List
      VLListDomainAtom=>VL%ListDomainAtom
      VLdrmag=>VL%drmag
      VLdr=>VL%dr
      PP=>Buckingham%PP
      SizePP=>Buckingham%SizePP
      RangePP=>Buckingham%RangePP
      IncrPPinv=>Buckingham%IncrPPinv
      
      errorstatus=0

      IF (VL%ListType==2) THEN !use half-list
         
         WRITE(6,*) "$Err>> Half list for Buckinghman not implemented"
         STOP
         DO iatom=1,NAtoms
            INCLUDE "PairPotentialForceHalfLinearInterp.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergyPP

      ELSEIF (VL%ListType==1) THEN !use full list

         DO iatom=1,NAtoms
            INCLUDE "PairPotentialForceFullLinearInterp.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+0.5_dp*EnergyPP

      ELSE

         WRITE(6,*) "$Err>> VL type could not be recognized"
         STOP

      END IF

   END SUBROUTINE BuckinghamForceLinearInterp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  Only Hessian evaluation xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE BuckinghamHessianGenLinInterp(AL,NAtoms)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: NAtoms
      
   END SUBROUTINE BuckinghamHessianGenLinInterp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
END MODULE BuckinghamSubroutines
