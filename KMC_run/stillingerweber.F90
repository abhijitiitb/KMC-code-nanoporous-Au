MODULE StillingerWeberSubroutines

   USE VARIABLE_TYPE
   IMPLICIT NONE
   INTEGER, PARAMETER :: maxatom=20000
  ! INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.0d0)
   
   CONTAINS
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  Only energy evaluation xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SWEnTaylor(AL,NAtoms) !SW Energy Generic (n species) Taylor
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
      

      !NAtoms=AL%NAtoms
      NSpecies=NSpecies_global
      
      VL=>AL%VL
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      SW=>AL%Potential%SW

      EnergyPP=0._dp

      AtomSpecies=>AL%AtomSpecies
      AtomCoord=>AL%AtomCoord
      AtomForce=>AL%AtomForce
      VLListRange=>VL%ListRange
      VLList=>VL%List
      VLdrmag=>VL%drmag
      VLListDomainAtom=>VL%ListDomainAtom
      VLdr=>VL%dr
      PP=>SW%PP
      F3=>SW%F3
      SizePP=>SW%SizePP
      SizeF3=>SW%SizeF3
      RangePP=>SW%RangePP
      RangeF3=>SW%RangeF3
      IncrPPinv=>SW%IncrPPinv
      IncrF3inv=>SW%IncrF3inv
      lam=>SW%lam
      epssqrt=>SW%epssqrt

      !AL%AtomForce=0._dp !never do this
      
      IF (AL%VL%ListType==2) THEN !use half-list
         
         WRITE(*,*) "SW>> Half list has not been implemented"
         STOP
         DO iatom=1,NAtoms
            !INCLUDE "stillingerweberEnTaylorHalf.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergyPP

      ELSEIF (AL%VL%ListType==1) THEN !use full list

         DO iatom=1,NAtoms
            !INCLUDE "stillingerweberEnTaylorFull.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+0.5_dp*EnergyPP

      ELSE

         WRITE(6,*) "$Err>> VL type could not be recognized"
         STOP

      END IF

   END SUBROUTINE SWEnTaylor
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Energy + Force evaluation xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SWForceLinearInterp(AL,NAtoms,errorstatus)
   !errorstatus=(0, no errors), (1, spacing less than lower PP cutoff)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(StillingerWeberPotential), POINTER :: SW
      INTEGER, DIMENSION(:), POINTER :: SizePP,SizeF3,VLListRange,VLList,AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom,AtomIsMoving
      INTEGER :: ra,r0,r1,r2,iatom,jatom,katom,j,j1,sectionij,sectionik,sectionjk
      INTEGER :: NSpecies,NAtoms,MaxAtomPerAtom
      INTEGER :: s1,s2,s3,smin,smax,smin1,smax1,pos,k,k0,errorstatus
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag,VLdr,AtomCoord,AtomForce
      REAL(dp), DIMENSION(:), POINTER :: PP,F3,RangePP,RangeF3,IncrPPinv,IncrF3inv,lam,epssqrt
      REAL(dp) :: drmagij,drmagik,drmagjk,F3ij,F3ik,F3jk,phi_ijk,phi_jik,phi_ikj
      REAL(dp):: tmp,p,gi,gi1,tgrad,F3ijgrad,F3ikgrad,F3jkgrad,costerm
      REAL(dp) :: epsik,epsij,epsjk,Lfnik,Lfnij,Lfnjk ! slope,tgphi,tgrho,grad(3),ro1
      REAL(dp) :: EnergySW,EnergySWTerm,mincutoff,dxinv,drij(3),drik(3),drjk(3)
      REAL(dp) :: costheta_ijk,costheta_jik,costheta_ikj
      REAL(dp) :: D1,D2,D3,D4,D5,D6,D7
      !REAL(dp) :: dhdr(3),drij1(3),drik1(3),drjk1(3)

      REAL(dp) :: forceoniatom, h_jik, h_ijk, h_ikj,EnergySWv3
      REAL(dp), DIMENSION(:), POINTER :: epsln,A,B,g,q,sigma,ac,cuttoff,lmd,gamma
      REAL(dp):: epslnsqrt
      REAL(dp)::  rij, rik, rjk, exporij, expogmrij, expogmrik, expogmrjk, rij4inv, rij2
      REAL(dp) :: dhdr(3),drij1(3),drik1(3),drjk1(3)

      !NAtoms=AL%NAtoms
      NSpecies=NSpecies_global
      NAtoms=AL%NAtoms
      
      VL=>AL%VL
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      SW=>AL%Potential%SW

      EnergySW=0._dp

      AtomSpecies=>AL%AtomSpecies
      AtomCoord=>AL%AtomCoord
      AtomForce=>AL%AtomForce
      AtomIsMoving=>AL%AtomIsMoving
      
      VLListRange=>VL%ListRange
      VLList=>VL%List
      VLListDomainAtom=>VL%ListDomainAtom
      VLdrmag=>VL%drmag
      VLdr=>VL%dr
      PP=>SW%PP
      F3=>SW%F3
      SizePP=>SW%SizePP
      SizeF3=>SW%SizeF3
      RangePP=>SW%RangePP
      RangeF3=>SW%RangeF3
      IncrPPinv=>SW%IncrPPinv
      IncrF3inv=>SW%IncrF3inv
      lam=>SW%lam
      epssqrt=>SW%epssqrt
      
      epsln=>SW%epsln                          ! 2.1678_dp
      A=>SW%A                                  ! 7.049556277_dp
      B=>SW%B                                  ! 0.6022245584_dp
      g=>SW%g                                  ! 4._dp
      q=>SW%q                                  ! 0._dp
      sigma=>SW%sigma                          ! 2.0951_dp
      ac=>SW%ac                                ! 1.8_dp
      lmd=>SW%lmd                              ! 21.0_dp
      gamma=>SW%gamma                          ! 1.2_dp
      cuttoff=>SW%cutoff                       ! 3.77_dp

      errorstatus=0
      IF (VL%ListType==2) THEN !use half-list
         
         WRITE(6,*) "Err>> SW-halflist is not working"
         STOP
         DO iatom=1,NAtoms
            INCLUDE "stillingerweberForceEntaylorHalf.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergySW
         
      ELSEIF (VL%ListType==1) THEN !use full list
        
         DO iatom=1,NAtoms
            !IF (ANY(AtomIsMoving(3*iatom-2:3*iatom))) THEN
               !INCLUDE "stillingerweberForceEntaylorFull.f90"
               INCLUDE "stillingerweberForceEntaylorFull-analytical.f90"
               !INCLUDE "stillingerweberForceEntaylorFull-table.f90"
            !END IF
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergySW

      ELSE

         WRITE(6,*) "$Err>> VL type could not be recognized"
         STOP

      END IF

   END SUBROUTINE SWForceLinearInterp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SWForceLinearInterp_OpenMP(AL,NAtoms,errorstatus)
   !errorstatus=(0, no errors), (1, spacing less than lower PP cutoff)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      TYPE(StillingerWeberPotential), POINTER :: SW
      INTEGER, DIMENSION(:), POINTER :: SizePP,SizeF3,VLListRange,VLList,AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom,AtomIsMoving
      INTEGER :: ra,r0,r1,r2,iatom,jatom,katom,j,j1,sectionij,sectionik,sectionjk
      INTEGER :: NSpecies,NAtoms,MaxAtomPerAtom
      INTEGER :: s1,s2,s3,smin,smax,smin1,smax1,pos,k,k0,errorstatus
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag,VLdr,AtomCoord,AtomForce
      REAL(dp), DIMENSION(:), POINTER :: PP,F3,RangePP,RangeF3,IncrPPinv,IncrF3inv,lam,epssqrt
      REAL(dp) :: drmagij,drmagik,drmagjk,F3ij,F3ik,F3jk,phi_ijk,phi_jik,phi_ikj
      REAL(dp):: tmp,p,gi,gi1,tgrad,F3ijgrad,F3ikgrad,F3jkgrad,costerm
      REAL(dp) :: epsik,epsij,epsjk,Lfnik,Lfnij,Lfnjk ! slope,tgphi,tgrho,grad(3),ro1
      REAL(dp) :: EnergySW,EnergySWTerm,mincutoff,dxinv,drij(3),drik(3),drjk(3)
      REAL(dp) :: costheta_ijk,costheta_jik,costheta_ikj
      REAL(dp) :: D1,D2,D3,D4,D5,D6,D7
      !REAL(dp) :: dhdr(3),drij1(3),drik1(3),drjk1(3)

      REAL(dp) :: forceoniatom, h_jik, h_ijk, h_ikj,EnergySWv3
      REAL(dp), DIMENSION(:), POINTER :: epsln,A,B,g,q,sigma,ac,cuttoff,lmd,gamma
      REAL(dp):: epslnsqrt
      REAL(dp)::  rij, rik, rjk, exporij, expogmrij, expogmrik, expogmrjk, rij4inv, rij2
      REAL(dp) :: dhdr(3),drij1(3),drik1(3),drjk1(3)
      INTEGER, PARAMETER :: chunksize=100

      !NAtoms=AL%NAtoms
      NSpecies=NSpecies_global
      NAtoms=AL%NAtoms
      
      VL=>AL%VL
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      SW=>AL%Potential%SW

      EnergySW=0._dp

      AtomSpecies=>AL%AtomSpecies
      AtomCoord=>AL%AtomCoord
      !AtomForce=>AL%AtomForce
      AtomIsMoving=>AL%AtomIsMoving
      
      VLListRange=>VL%ListRange
      VLList=>VL%List
      VLListDomainAtom=>VL%ListDomainAtom
      VLdrmag=>VL%drmag
      VLdr=>VL%dr
      PP=>SW%PP
      F3=>SW%F3
      SizePP=>SW%SizePP
      SizeF3=>SW%SizeF3
      RangePP=>SW%RangePP
      RangeF3=>SW%RangeF3
      IncrPPinv=>SW%IncrPPinv
      IncrF3inv=>SW%IncrF3inv
      lam=>SW%lam
      epssqrt=>SW%epssqrt
      
      epsln=>SW%epsln                          ! 2.1678_dp
      A=>SW%A                                  ! 7.049556277_dp
      B=>SW%B                                  ! 0.6022245584_dp
      g=>SW%g                                  ! 4._dp
      q=>SW%q                                  ! 0._dp
      sigma=>SW%sigma                          ! 2.0951_dp
      ac=>SW%ac                                ! 1.8_dp
      lmd=>SW%lmd                              ! 21.0_dp
      gamma=>SW%gamma                          ! 1.2_dp
      cuttoff=>SW%cutoff                       ! 3.77_dp

      errorstatus=0
      IF (VL%ListType==2) THEN !use half-list
         
         WRITE(6,*) "Err>> SW-halflist is not working"
         STOP
         DO iatom=1,NAtoms
            INCLUDE "stillingerweberForceEntaylorHalf.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergySW
         
      ELSEIF (VL%ListType==1) THEN !use full list

!$OMP PARALLEL PRIVATE(AtomForce,EnergySW) 
         ALLOCATE(AtomForce(3*NAtoms))
         AtomForce=0._dp
         EnergySW=0._dp
!$OMP DO SCHEDULE(DYNAMIC,chunksize)
         DO iatom=1,NAtoms
            CALL SW_analytical(iatom,EnergySW,AtomForce,AtomSpecies, &
     VLdr,VLdrmag,VLList,VLListRange,MaxAtomPerAtom,NSpecies, &
     epsln,A,B,g,q,sigma,ac,cuttoff,lmd,gamma) 
         END DO
!$OMP END DO
!$OMP CRITICAL
         AL%PotentialEnergy=AL%PotentialEnergy+EnergySW
         AL%AtomForce=AL%AtomForce+AtomForce
!$OMP END CRITICAL
         DEALLOCATE(AtomForce)
!$OMP END PARALLEL

      ELSE

         WRITE(6,*) "$Err>> VL type could not be recognized"
         STOP

      END IF

   END SUBROUTINE SWForceLinearInterp_OpenMP
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SW_analytical(iatom,EnergySW,AtomForce,AtomSpecies, &
      VLdr,VLdrmag,VLList,VLListRange,MaxAtomPerAtom,NSpecies, &
      epsln,A,B,g,q,sigma,ac,cuttoff,lmd,gamma)
      IMPLICIT NONE         

      INTEGER :: ra,r0,r1,r2,iatom,jatom,katom,j,j1,sectionij,sectionik,sectionjk
      INTEGER :: MaxAtomPerAtom,NSpecies
      INTEGER :: s1,s2,s3,smin,smax,smin1,smax1,pos,k,k0,errorstatus
      REAL(dp) :: drmagij,drmagik,drmagjk,F3ij,F3ik,F3jk,phi_ijk,phi_jik,phi_ikj
      REAL(dp):: tmp,p,gi,gi1,tgrad,F3ijgrad,F3ikgrad,F3jkgrad,costerm
      REAL(dp) :: epsik,epsij,epsjk,Lfnik,Lfnij,Lfnjk ! slope,tgphi,tgrho,grad(3),ro1
      REAL(dp) :: EnergySW,EnergySWTerm,mincutoff,dxinv,drij(3),drik(3),drjk(3)
      REAL(dp) :: costheta_ijk,costheta_jik,costheta_ikj
      REAL(dp) :: D1,D2,D3,D4,D5,D6,D7
      !REAL(dp) :: dhdr(3),drij1(3),drik1(3),drjk1(3)

      REAL(dp), DIMENSION(:), POINTER :: AtomForce,VLdr,VLdrmag
      REAL(dp), DIMENSION(:), POINTER :: epsln,A,B,g,q,sigma,ac,cuttoff,lmd,gamma
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies,VLListRange,VLList

      REAL(dp) :: forceoniatom, h_jik, h_ijk, h_ikj,EnergySWv3
      REAL(dp):: epslnsqrt
      REAL(dp)::  rij, rik, rjk, exporij, expogmrij, expogmrik, expogmrjk, rij4inv, rij2     
      REAL(dp) :: dhdr(3),drij1(3),drik1(3),drjk1(3)
  
      INCLUDE "stillingerweberForceEntaylorFull-analytical.f90"
   END SUBROUTINE SW_analytical
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   !SUBROUTINE SWForceLinearInterpLL(AL,NAtoms,errorstatus)
   !uses LL to find the energy and forces
   !   IMPLICIT NONE
      
      
   !   c1i=CEILING((AtomCoord(3*iatom-2:3*iatom)-xc)/CellSize) !location of atom i
   !   c1i=MAX(c1i,1)
      
   !END SUBROUTINE SWForceLinearInterpLL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE StillingerWeberSubroutines
