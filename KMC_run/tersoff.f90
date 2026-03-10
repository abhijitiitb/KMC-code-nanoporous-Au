MODULE TersoffSubroutines

   USE VARIABLE_TYPE
   IMPLICIT NONE
   INTEGER, PARAMETER :: maxatom=20000
  ! INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.0d0)
   
   CONTAINS
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  Only energy evaluation xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE TersoffEnTaylor(AL,NAtoms) !Tersoff Energy Generic (n species) Taylor
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(StillingerWeberPotential), POINTER :: SW
      INTEGER, DIMENSION(:), POINTER :: SizePP,SizeF3,VLListRange,VLList,AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom
      INTEGER :: ra,r0,r1,r2,iatom,jatom,katom,j,j1,section,section1,NSpecies,NAtoms,MaxAtomPerAtom
      INTEGER :: s1,s2,s3,smin,smax,smin1,smax1,pos,k,k0,MaxSWTableSize
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag,VLdr,AtomCoord,AtomForce
      REAL(dp), DIMENSION(:), POINTER :: PP,F3,RangePP,RangeF3,IncrPPinv,IncrF3inv,lam,epssqrt
      REAL(dp) :: drmagij,drmagik,drmagji,drmagjk,drmagki,drmagkj,phijik,phiijk,phiikj
      REAL(dp):: tmp,p,costerm,epsik,epsij,epsji,epsjk,epski,epskj,F3ij,F3ik,F3ji,F3jk,F3ki,F3kj,Lfn
      REAL(dp) :: EnergyPP,mincutoff,dxinv,EnergyPPTerm,drij(3),drik(3),drji(3),drjk(3),drki(3),drkj(3)
      

   END SUBROUTINE TersoffEnTaylor
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Energy + Force evaluation xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE TersoffForceLinearInterp(AL,NAtoms,errorstatus)
   !errorstatus=(0, no errors), (1, spacing less than lower PP cutoff)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(TersoffPotential), POINTER :: Tersoff
      INTEGER, DIMENSION(:), POINTER :: VLListRange,VLList,AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom
      INTEGER :: ra,r0,r1,r2,iatom,jatom,katom
      INTEGER :: raj,r0j,r1j,r2j
      INTEGER :: j,j1,sectionij,sectionik,sectionjk
      INTEGER :: NSpecies,NAtoms,MaxAtomPerAtom
      INTEGER :: s1,s2,s3,smin,smax,smin1,smax1,pos,k,k0,errorstatus
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag,VLdr,AtomCoord,AtomForce
      REAL(dp) :: EnergyTersoff,EnergySWTerm,mincutoff,drmagij,drmagik,drmagjk
      REAL(dp) :: dxinv,drij(3),drik(3),drjk(3)
      REAL(dp), DIMENSION(:), POINTER :: A,B,lam1,lam2,lam3,beta,n,c,d,h,R,RD,Rangefc
      REAL(dp) :: Aij,Bij,lam1ij,lam2ij,lam3ij,betaij,Rij,RDij,cij,dij,hij,nij
      REAL(dp) :: frepl,fattr,drt,fcut,b_ij,b_ji,Rangefcij
      REAL(dp) :: dfcut,dfcutik,dfcutjk,dfattr,dfrepl
      REAL(dp) :: drij_rixi(3),drij_rjxi(3),drik_rixi(3),drjk_rjxi(3)
      REAL(dp) :: zetaij,zetaji,drtik,fcutik,costheta,gtheta,lterm,expr,drtjk,fcutjk
      REAL(dp) :: gterm,zterm,bzn_ij,bzn_ji,bterm_ij,bterm_ji
      REAL(dp) :: dcostheta_i(3),dcostheta_j(3),dg_i(3),dg_j(3),dzetaij_i(3),dzetaij_j(3),dzetaji_i(3),dzetaji_j(3)
      REAL(dp) :: db_ij_drixi(3),db_ij_drjxi(3),db_ji_drixi(3),db_ji_drjxi(3),fterm(3),ForceTersoff_i(3),ForceTersoff_j(3)
      
      !NAtoms=AL%NAtoms
      NSpecies=NSpecies_global
      NAtoms=AL%NAtoms
      
      VL=>AL%VL
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      Tersoff=>AL%Potential%Tersoff

      EnergyTersoff=0._dp

      AtomSpecies=>AL%AtomSpecies
      AtomCoord=>AL%AtomCoord
      AtomForce=>AL%AtomForce
      VLListRange=>VL%ListRange
      VLList=>VL%List
      VLListDomainAtom=>VL%ListDomainAtom
      VLdrmag=>VL%drmag
      VLdr=>VL%dr
      A=>Tersoff%A
      B=>Tersoff%B
      lam1=>Tersoff%lam1
      lam2=>Tersoff%lam2
      lam3=>Tersoff%lam3
      beta=>Tersoff%beta
      n=>Tersoff%n
      c=>Tersoff%c
      d=>Tersoff%d
      h=>Tersoff%h
      R=>Tersoff%R
      RD=>Tersoff%RD
      Rangefc=>Tersoff%Rangefc
      
      errorstatus=0
       
      IF (VL%ListType==2) THEN !use half-list
      
         DO iatom=1,NAtoms
            !INCLUDE "stillingerweberForceEntaylorHalf.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergyTersoff
         
      ELSEIF (VL%ListType==1) THEN !use full list
         
         DO iatom=1,NAtoms
            INCLUDE "tersoffForceEntaylorFull.f90"
            Stop
         END DO
         write(*,*) "# atoms" , NAtoms
         AL%PotentialEnergy=AL%PotentialEnergy+EnergyTersoff

      ELSE

         WRITE(6,*) "$Err>> VL type could not be recognized"
         STOP

      END IF

   END SUBROUTINE TersoffForceLinearInterp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE TersoffSubroutines
