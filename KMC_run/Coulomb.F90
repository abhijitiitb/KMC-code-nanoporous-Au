!AL=Atom list
MODULE CoulombSubroutines
   !will ADD in energy, forces and hessian values for the input coordinates and Verlet list
   !MAKE SURE THAT AL energy, forces etc. have been properly initialized ...
   !the Ewald sum uses a VL to compute the real part - the cutoff used is 
   !alpha=sqrt(PI)*(tR/tF N/V^2)^(1/6)
   !rc=sqrt(p)
   USE VARIABLE_TYPE
   USE CoulombInitialize
   IMPLICIT NONE
         
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CoulombEnWolf(AL,NAtoms)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      INTEGER :: iatom,jatom,j,NAtoms,ra,r0,r1,r2,MaxAtomPerAtom
      REAL(dp) :: Rc,alfa,earc,q2sum,EnergyPP,z1,z2,drmag
      INTEGER, DIMENSION(:), POINTER :: VLList,VLListRange
      REAL(dp), DIMENSION(:), POINTER :: AtomCharge,VLdrmag
      
      !NAtoms=AL%NAtoms
      VL=>AL%VL
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      VLListRange=>VL%ListRange
      VLList=>VL%List
      VLdrmag=>VL%drmag
      AtomCharge=>AL%AtomCharge
      
      Rc=AL%Potential%Coulomb%Wolf%ReCutOff
      alfa=AL%Potential%Coulomb%Wolf%ALPHA
      
      earc=ERFC(alfa*Rc)/Rc
      
      q2sum=0._dp
      DO iatom=1,NAtoms
         q2sum=q2sum+AtomCharge(iatom)*AtomCharge(iatom)
      END DO
      
      IF (VL%ListType==2) THEN !half list
         DO iatom=1,NAtoms
            INCLUDE "WolfEnTaylorHalf.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergyPP-(earc/2._dp+alfa/SQRT(PI))*q2sum
      ELSE
         DO iatom=1,NAtoms
            INCLUDE "WolfEnTaylorFull.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+0.5_dp*EnergyPP-(0.5_dp*earc+alfa/SQRT(PI))*q2sum
      END IF
   END SUBROUTINE CoulombEnWolf
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CoulombEnEwald(AL,NAtoms)
      !compute coulomb energy terms using Ewald summation
      !Linked (Verlet) list is used when ListType=1 (2)
      IMPLICIT NONE
      REAL(dp), PARAMETER :: re_fac=14.399644549472361_dp !real term - conversion factor to get qi qj/(4pi eps0) r in eV [obtained as 1.6021765314e-19/4/pi/8.854187817e-12/1e-10]
      REAL(dp), PARAMETER :: im_fac=1.809512701237067e2_dp !fourier term - conversion factor [obtained as 1.6021765314e-19/1e-10/8.854187817e-12]
      REAL(dp), PARAMETER :: se_fac=-8.124129461602530_dp !self term - conversion factor [obtained as -1.6021765314e-19/8.854187817e-12/4/pi^1.5/1e-10]
      TYPE(SystemContainer), POINTER :: AL
      TYPE(EwaldParameters), POINTER :: Ew
      TYPE(LinkedListContainer), POINTER :: LL
      TYPE(VerletListContainer), POINTER :: VL
      INTEGER :: ListType
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomCharge,KVEC,VLdrmag
      INTEGER, DIMENSION(:), POINTER :: RVEC,LLList,LLListRange,VLList,VLListRange,KVector
      COMPLEX(dp), DIMENSION(:,:), POINTER :: EIKX,EIKY,EIKZ
      COMPLEX(dp) :: EIKR_SUM
      INTEGER :: ra,r0 !used for the Verlet list part
      INTEGER :: NAtoms,TOTr,TOTk,i,j,k,KX,KY,KZ,KSQ,KMAX,maincellindx,neighcellindx
      INTEGER :: c1(3),c2(3),NCell(3),c1start,c1end,c2start,c2end,ivec(3) !,nx,ny,nz
      INTEGER :: MaxAtomPerCell,MaxAtomPerAtom,atomindx1,atomindx2,c1pos,c2pos
      REAL(dp) :: z1,z2,r1(3),r2(3),dr(3),drbox(3),Volume
      REAL(dp) :: iLx2pi,iLy2pi,iLz2pi
      REAL(dp) :: vij,vd,vs,vk,vr,RX,RY,RZ,drmag,fac,ReCutOff
      REAL(dp) :: ALPHA,B,EIKR,BoxSize(3),iBoxSize(3)
      
      AtomCharge=>AL%AtomCharge
      AtomCoord=>AL%AtomCoord
      Ew=>AL%Potential%Coulomb%Ewald
      LL=>AL%LL
      LLList=>LL%List
      LLListRange=>LL%ListRange
      MaxAtomPerCell=LL%MaxAtomPerCell
      VL=>AL%VL
      VLList=>VL%List
      VLListRange=>VL%ListRange
      VLdrmag=>VL%drmag
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      
      !NAtoms=AL%NAtoms
      IF (NAtoms/=Ew%NAtoms) CALL EwaldSetUp(AL)
      ALPHA=Ew%ALPHA
      ReCutOff=Ew%ReCutOff
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Real space part
      vr=0._dp
      
      DO atomindx1=1,NAtoms !find the Ewald real space part for energy using a half list -- in case a full list is used, remember to set vr=vr/2.
         
         z1=AtomCharge(atomindx1)
         ra=(atomindx1-1)*MaxAtomPerAtom !start index minus 1
         r0=VLListRange(atomindx1) !range for ith atom
         !WRITE(*,*) "Main atom:",atomindx1,z1,ra,r0
         
         vij=0._dp
         DO j=ra+1,ra+r0
            atomindx2=VLList(j)
            drmag=VLdrmag(j)
            z2=AtomCharge(atomindx2)
            vij=vij+z2*ERFC(ALPHA*drmag)/drmag
            !WRITE(*,*) "Neig atom:",atomindx2,drmag,z2,vij
         END DO
         vr=vr+vij*z1
         
      END DO
      IF (AL%VL%ListType==1) vr=vr*0.5_dp !for full VL neighbor list
      vr=vr*re_fac
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Reciprocal space part
      
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      Volume=PRODUCT(BoxSize) !in Angstrom cube
      iBoxSize=1._dp/BoxSize
      !EIKX, EIKY, EIKZ no longer exist (see force calculation)
      !EIKX=>Ew%EIKX
      !EIKY=>Ew%EIKY
      !EIKZ=>Ew%EIKZ
      KMAX=Ew%ImCutOff
      
      iLx2pi=TWOPI*iBoxSize(1)
      iLy2pi=TWOPI*iBoxSize(2)
      iLz2pi=TWOPI*iBoxSize(3)
      
      EIKX(:,0)=CMPLX(1._dp,0._dp)
      EIKY(:,0)=CMPLX(1._dp,0._dp)
      EIKZ(:,0)=CMPLX(1._dp,0._dp)
      
      DO atomindx1=1,NAtoms
         EIKX(atomindx1,0)=CMPLX(1.0_dp,0.0_dp)
         EIKY(atomindx1,0)=CMPLX(1.0_dp,0.0_dp)
         EIKZ(atomindx1,0)=CMPLX(1.0_dp,0.0_dp)
         
         j=3*atomindx1-2; RX=AtomCoord(j)*iLx2pi
         EIKX(atomindx1,1)=CMPLX(COS(RX),SIN(RX))
         j=j+1; RY=AtomCoord(j)*iLy2pi
         EIKY(atomindx1,1)=CMPLX(COS(RY),SIN(RY))
         j=j+1; RZ=AtomCoord(j)*iLz2pi
         EIKZ(atomindx1,1)=CMPLX(COS(RZ),SIN(RZ))
      
         EIKY(atomindx1,-1)=CONJG(EIKY(atomindx1,1))
         EIKZ(atomindx1,-1)=CONJG(EIKZ(atomindx1,1))
      END DO
      
      !set up EIKX(1:NAtoms,k) etc.
      DO k=2,KMAX
         EIKX(:,k)=EIKX(:,k-1)*EIKX(:,1)
         EIKY(:,k)=EIKY(:,k-1)*EIKY(:,1)
         EIKY(:,-k)=CONJG(EIKY(:,k))
         EIKZ(:,k)=EIKZ(:,k-1)*EIKZ(:,1)
         EIKZ(:,-k)=CONJG(EIKZ(:,k))
      END DO
      
      vd=0._dp
      TOTk=Ew%TOTk
      KVector=>Ew%KVector
      KVEC=>Ew%KVEC
      DO k=1,TOTk
         KX=KVector(3*k-2)
         KY=KVector(3*k-1)
         KZ=KVector(3*k  )
         fac=2._dp
         IF (KX==0) fac=1._dp
         EIKR_SUM=SUM(AtomCharge*EIKX(:,KX)*EIKY(:,KY)*EIKZ(:,KZ))
         vd=vd+fac*KVEC(k)*CONJG(EIKR_SUM)*EIKR_SUM
      END DO
      
      vd=vd*im_fac/Volume
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Self-part
      vs=SUM(AtomCharge*AtomCharge)
      vs=ALPHA*vs*se_fac
      
      AL%PotentialEnergy=AL%PotentialEnergy+vr+vk+vs
      
   END SUBROUTINE CoulombEnEwald
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CoulombEnPPPM
   END SUBROUTINE CoulombEnPPPM
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CoulombEnFastMultipole
   END SUBROUTINE CoulombEnFastMultipole
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CoulombForceWolf(AL,NAtoms)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      INTEGER :: iatom,jatom,j,NAtoms,ra,r0,r1,r2,MaxAtomPerAtom
      REAL(dp) :: Rc,alfa,ear,earc,q2sum,EnergyPP,z1,z2,drmag,slope,grad(3),tgrad
      INTEGER, DIMENSION(:), POINTER :: VLList,VLListRange
      REAL(dp), DIMENSION(:), POINTER :: AtomCharge,AtomForce,VLdrmag,VLdr
      
      !NAtoms=AL%NAtoms
      VL=>AL%VL
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      VLListRange=>VL%ListRange
      VLList=>VL%List
      VLdrmag=>VL%drmag
      VLdr=>VL%dr
      AtomCharge=>AL%AtomCharge
      AtomForce=>AL%AtomForce
      
      Rc=AL%Potential%Coulomb%Wolf%ReCutOff
      alfa=AL%Potential%Coulomb%Wolf%ALPHA
      
      earc=ERFC(alfa*Rc)/Rc
      
      q2sum=0._dp
      DO iatom=1,NAtoms
         q2sum=q2sum+AtomCharge(iatom)*AtomCharge(iatom)
      END DO
      
      IF (VL%ListType==2) THEN !half list
         DO iatom=1,NAtoms
            INCLUDE "WolfForceHalfLinearInterp.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+EnergyPP-(earc/2._dp+alfa/SQRT(PI))*q2sum
      ELSE
         DO iatom=1,NAtoms
            INCLUDE "WolfForceFullLinearInterp.f90"
         END DO
         AL%PotentialEnergy=AL%PotentialEnergy+0.5_dp*EnergyPP-(0.5_dp*earc+alfa/SQRT(PI))*q2sum
      END IF
      
   END SUBROUTINE CoulombForceWolf
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CoulombForceEwald(AL,NAtoms,errorstatus)
      !compute coulomb energy and force terms using Ewald summation
      IMPLICIT NONE
      
      REAL(dp) :: theta0=1.91061193_dp,KSpring=1.322626_dp
      
      REAL(dp), PARAMETER :: re_fac=14.399644549472361_dp !real term - conversion factor to get qi qj/(4pi eps0) r in eV [obtained as 1.6021765314e-19/4/pi/8.854187817e-12/1e-10]
      REAL(dp), PARAMETER :: im_fac=1.809512701237067e2_dp !fourier term - conversion factor [obtained as 1.6021765314e-19/1e-10/8.854187817e-12]
      REAL(dp), PARAMETER :: se_fac=-8.124129461602530_dp !self term - conversion factor [obtained as -1.6021765314e-19/8.854187817e-12/4/pi^1.5/1e-10]
      INTEGER :: ListType,errorstatus
      TYPE(SystemContainer), POINTER :: AL
      TYPE(EwaldParameters), POINTER :: Ew
      TYPE(LinkedListContainer), POINTER :: LL
      TYPE(VerletListContainer), POINTER :: VL
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomCharge,KVEC,force,VLdrmag,VLdr
      INTEGER, DIMENSION(:), POINTER :: LLList,LLListRange,VLList,VLListRange,KVector,AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom
      COMPLEX(dp), DIMENSION(:,:), POINTER :: EIKX,EIKY,EIKZ
      COMPLEX(dp), DIMENSION(:), POINTER :: QEIKR
      COMPLEX(dp) :: EIKR_SUM
      INTEGER :: r0,ra !used with VL
      INTEGER :: NAtoms,TOTk,i,j,k,KX,KY,KZ,KSQ,KMAX,maincellindx,neighcellindx !,TOTr
      INTEGER :: c1(3),c2(3),NCell(3),nx,ny,nz,c1start,c1end,c2start,c2end,ivec(3)
      INTEGER :: MaxAtomPerCell,MaxAtomPerAtom,atomindx1,atomindx2,c1pos,c2pos
      REAL(dp) :: z1,z2,r1(3),r2(3),dr(3),drbox(3),Volume,erfckrij
      REAL(dp) :: iLx2pi,iLy2pi,iLz2pi,krij,kmag2,KSCALED(3),QCOSKR,QSINKR
      REAL(dp) :: vij,fij(3),vd,vs,vk,vr,RX,RY,RZ,drmag,fac,ReCutOff,v3
      REAL(dp) :: ALPHA,B,EIKR,BoxSize(3),iBoxSize(3),b1,fmax,fbu(3),BuckCutOff
      REAL(dp) :: vbuck
      INTEGER :: s1,s2,s3
      REAL(dp) :: eb,fb,qsum
      
      !Defined for LiFePO4
      INTEGER :: atomindx3
      REAL(dp) :: dr2(3),dr3(3),drmag2,ESpring,FSpring1(3),FSpring2(3),FSpring3(3)
      
      errorstatus=0
      
      AtomCharge=>AL%AtomCharge
      AtomCoord=>AL%AtomCoord
      AtomSpecies=>AL%AtomSpecies
      Ew=>AL%Potential%Coulomb%Ewald
      LL=>AL%LL
      LLList=>LL%List
      LLListRange=>LL%ListRange
      MaxAtomPerCell=LL%MaxAtomPerCell
      VL=>AL%VL
      VLList=>VL%List
      VLListRange=>VL%ListRange
      VLdrmag=>VL%drmag
      VLdr=>VL%dr
      VLListDomainAtom=>VL%ListDomainAtom
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      
      !NAtoms=AL%NAtoms
      IF (NAtoms/=Ew%NAtoms) THEN
         WRITE(6,*) "Setting up Ewald ..."
         CALL EwaldSetUp(AL)
      END IF
      ALPHA=Ew%ALPHA
      ReCutOff=Ew%ReCutOff
      !BuckCutOff=AL%Potential%MaxPotentialCutoff(4)
      !BuckCutOff=maxrcutbuck
      
      qsum=SUM(AtomCharge(1:NAtoms))
      IF (qsum/=0._dp) THEN
         WRITE(6,*) "$Err>> Net charge is ",qsum
         STOP
      END IF
      
      ALLOCATE(force(3*NAtoms))
      force=0._dp
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Real space part
      vr=0._dp
      vbuck=0._dp
      b1=2._dp*ALPHA/SQRT(PI)
      
      fmax=0._dp
      IF (AL%VL%ListType==1) THEN !full list
      
         DO atomindx1=1,NAtoms
         
         fbu=0._dp
         
            z1=AtomCharge(atomindx1)
            s1=AtomSpecies(atomindx1)
            ra=(atomindx1-1)*MaxAtomPerAtom !start index minus 1
            r0=VLListRange(atomindx1) !range for ith atom
            
            vij=0._dp
            v3=0._dp
            DO j=ra+1,ra+r0
               atomindx2=VLList(j)
               
               !According to MOLDY atomindx1 /= atomindx2
               !We have not implemented this condition
               
               
               dr=VLdr(3*j-2:3*j)
               drmag=VLdrmag(j)
               z2=AtomCharge(atomindx2)
               s2=AtomSpecies(atomindx2)

               !spring forces specifically for O-P-O in LiFePO4
               IF (s1==3 .AND. s2==4 .AND. drmag<1.8_dp) THEN
               !for LiFePO4 the average spacing between O and P is 1.6 Ang
                  DO k=j+1,ra+r0 !that should give unique triplets involving at least 1 O and 1 P
                     atomindx3=VLList(k)
                     s3=AtomSpecies(atomindx3)
                     dr3=VLdr(3*k-2:3*k)
                     drmag2=VLdrmag(k)
                     IF (s2==4 .AND. drmag2<1.8_dp) THEN
                        CALL SPRING(dr,dr3,drmag,drmag2,KSpring,theta0,ESpring,FSpring1,FSpring2,FSpring3)
                        v3=v3+ESpring
                        force(3*atomindx2-2:3*atomindx2)=force(3*atomindx2-2:3*atomindx2)+FSpring2 !atom j
                        force(3*atomindx1-2:3*atomindx1)=force(3*atomindx1-2:3*atomindx1)+FSpring1 !center atom
                        force(3*atomindx3-2:3*atomindx3)=force(3*atomindx3-2:3*atomindx3)+FSpring3 !atiom k
                     END IF
                  END DO
               END IF
               
               BuckCutOff=rcutbuck(s1,s2)
               !write(*,*) ".....",s1,s2,ABuck(s1,s2),rhoBuck(s1,s2),CBuck(s1,s2),ReCutOff
               eb=EBUCK(drmag,ABuck(s1,s2),rhoBuck(s1,s2),CBuck(s1,s2),BuckCutOff)
               fb=FBUCK(drmag,ABuck(s1,s2),rhoBuck(s1,s2),CBuck(s1,s2),BuckCutOff)

               IF (drmag<=ReCutOff) THEN
                  !write(*,*) atomindx2,dr,drmag
                  krij=ALPHA*drmag
                  erfckrij=ERFC(krij)/drmag
                  vij=z1*z2*erfckrij*re_fac+eb 
                  fij=z1*z2*(erfckrij+b1*EXP(-krij*krij))*dr/drmag/drmag*re_fac+fb*dr
               ELSE
                  vij=eb
                  fij=-fb*dr
               END IF
               vbuck=vbuck+eb
               !write(*,*) drmag,z1*z2*(erfckrij+b1*EXP(-krij*krij))*dr/drmag/drmag*re_fac
               fbu=fbu-fb*dr
               force(3*atomindx1-2:3*atomindx1)=force(3*atomindx1-2:3*atomindx1)+fij
            !write(*,*) "Ag ",dr
            !write(UNIT=6,FMT='(2I2,2I5,6ES15.5)') s1,s2,atomindx1,atomindx2,dr,-fb*dr
               vr=vr+vij
         !IF (atomindx1==1) THEN
         !   write(*,*)atomindx2,z2,krij,vij
         !END IF
            END DO
         END DO
         
         vr=vr*0.5_dp+v3 !for full VL neighbor list
         vbuck=vbuck*0.5_dp
         
      ELSE !VL is half list
         
         WRITE(6,*) "$Err>> Half-list has not been checked"
         STOP
         DO atomindx1=1,NAtoms

            z1=AtomCharge(atomindx1)
            s1=AtomSpecies(atomindx1)
            ra=(atomindx1-1)*MaxAtomPerAtom !start index minus 1
            r0=VLListRange(atomindx1) !range for ith atom
            
            DO j=ra+1,ra+r0
               atomindx2=VLList(j)
               dr=VLdr(3*j-2:3*j)
               drmag=VLdrmag(j)
               z2=AtomCharge(atomindx2)
               s2=AtomSpecies(atomindx2)

               !spring forces specifically for O-P-O in LiFePO4
               IF (s1==3 .AND. s2==4 .AND. drmag<1.8_dp) THEN
               !for LiFePO4 the average spacing between O and P is 1.6 Ang
                  DO k=j+1,ra+r0 !that should give unique triplets involving at least 1 O and 1 P
                     atomindx3=VLList(k)
                     s2=AtomSpecies(atomindx3)
                     dr2=VLdr(3*k-2:3*k)
                     drmag2=VLdrmag(k)
                     IF (s2==4 .AND. drmag2<1.8_dp) THEN
                        CALL SPRING(dr,dr2,drmag,drmag2,KSpring,theta0,ESpring,FSpring1,FSpring2,FSpring3)
                        vij=vij+ESpring
                        force(3*atomindx2-2:3*atomindx2)=force(3*atomindx2-2:3*atomindx2)+FSpring1
                        force(3*atomindx1-2:3*atomindx1)=force(3*atomindx1-2:3*atomindx1)+FSpring2
                        force(3*atomindx3-2:3*atomindx3)=force(3*atomindx3-2:3*atomindx3)+FSpring3
                     END IF
                  END DO
               END IF
               
               !buckingham forces
               eb=EBUCK(drmag,ABuck(s1,s2),rhoBuck(s1,s2),CBuck(s1,s2),ReCutOff)
               fb=FBUCK(drmag,ABuck(s1,s2),rhoBuck(s1,s2),CBuck(s1,s2),ReCutOff)
               
               !coulomb forces
               krij=ALPHA*drmag
               erfckrij=ERFC(krij)/drmag
               vij=vij+z1*z2*erfckrij*re_fac+eb
               fij=z1*z2*(erfckrij+b1*EXP(-krij*krij))*dr/drmag/drmag*re_fac-fb*dr
               force(3*atomindx1-2:3*atomindx1)=force(3*atomindx1-2:3*atomindx1)+fij
               IF (VLListDomainAtom(j)) &
                  force(3*atomindx2-2:3*atomindx2)=force(3*atomindx2-2:3*atomindx2)-fij
            END DO
            vr=vr+vij
      
         END DO
      END IF
      
      AL%AtomForce(1:3*NAtoms)=AL%AtomForce(1:3*NAtoms)-force(1:3*NAtoms) !eV/Angstrom
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Reciprocal space part -- atomic coordinates have to be in periodic box
      !CALL PBC(AL)
      
      force=0._dp
      ALLOCATE(QEIKR(NAtoms))
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      Volume=PRODUCT(BoxSize) !in Angstrom cube
      iBoxSize=1._dp/BoxSize
      !EIKX=>Ew%EIKX
      !EIKY=>Ew%EIKY
      !EIKZ=>Ew%EIKZ
      KMAX=Ew%ImCutOff
      ALLOCATE(EIKX(1:NAtoms,0:KMAX))
      ALLOCATE(EIKY(1:NAtoms,-KMAX:KMAX))
      ALLOCATE(EIKZ(1:NAtoms,-KMAX:KMAX))
      
      
      iLx2pi=TWOPI*iBoxSize(1)
      iLy2pi=TWOPI*iBoxSize(2)
      iLz2pi=TWOPI*iBoxSize(3)
      
      !WARNING
      !23/12/2012: gfortran appears to be buggy
      !when EIKX is being initialized, it made AL%Potential%SW to be associated
      !this did not happen with intel compiler
      DO atomindx1=1,NAtoms !effect of bug was witnessed here
         EIKX(atomindx1,0)=CMPLX(1.d0,0.d0)
         EIKY(atomindx1,0)=CMPLX(1.d0,0.d0)
         EIKZ(atomindx1,0)=CMPLX(1.d0,0.d0)
         
         j=3*atomindx1-2; RX=AtomCoord(j)*iLx2pi
         EIKX(atomindx1,1)=CMPLX(COS(RX),SIN(RX))
         j=j+1; RY=AtomCoord(j)*iLy2pi
         EIKY(atomindx1,1)=CMPLX(COS(RY),SIN(RY))
         j=j+1; RZ=AtomCoord(j)*iLz2pi
         EIKZ(atomindx1,1)=CMPLX(COS(RZ),SIN(RZ))
      
         EIKY(atomindx1,-1)=CONJG(EIKY(atomindx1,1))
         EIKZ(atomindx1,-1)=CONJG(EIKZ(atomindx1,1))
      END DO
      
      !set up EIKX(1:NAtoms,k) etc.
      DO k=2,KMAX
         EIKX(:,k)=EIKX(:,k-1)*EIKX(:,1)
         EIKY(:,k)=EIKY(:,k-1)*EIKY(:,1)
         EIKY(:,-k)=CONJG(EIKY(:,k))
         EIKZ(:,k)=EIKZ(:,k-1)*EIKZ(:,1)
         EIKZ(:,-k)=CONJG(EIKZ(:,k))
      END DO
      
      vd=0._dp
      TOTk=Ew%TOTk
      KVector=>Ew%KVector
      KVEC=>Ew%KVEC
      DO k=1,TOTk
         KX=KVector(3*k-2)
         KY=KVector(3*k-1)
         KZ=KVector(3*k  )
         fac=2._dp
         IF (KX==0) fac=1._dp
         QEIKR(1:NAtoms)=AtomCharge(1:NAtoms)*EIKX(1:NAtoms,KX)*EIKY(1:NAtoms,KY)*EIKZ(1:NAtoms,KZ) !dimension is natoms
         EIKR_SUM=SUM(QEIKR(1:NAtoms))
         QCOSKR=REAL(EIKR_SUM)
         QSINKR=AIMAG(EIKR_SUM)
         vd=vd+KVEC(k)*(CONJG(EIKR_SUM)*EIKR_SUM) !*fac
         KSCALED=REAL(KVECTOR(3*k-2:3*k),dp)*KVEC(k)*iBoxSize*TWOPI
         DO j=1,NAtoms
            force(3*j-2:3*j)=force(3*j-2:3*j)+ & 
             (AIMAG(QEIKR(j))*QCOSKR-REAL(QEIKR(j))*QSINKR)*KSCALED
            !alternatively, 
            !force(3*j-2:3*j)=force(3*j-2:3*j)+AIMAG(CONJ(EIKR_SUM)*QEIKR(j))*KSCALED
!IF (j==201) THEN
!IF(kx==3 .AND. ky==4 .and. kz==5 .AND. j==201) THEN
!write(*,*) (AIMAG(QEIKR(j))*QCOSKR-REAL(QEIKR(j))*QSINKR)*KSCALED*2._dp*im_fac/Volume/27.21_dp
!write(*,*) (AIMAG(QEIKR(j))*QCOSKR-REAL(QEIKR(j))*QSINKR)*KVEC(k)*2._dp*im_fac/Volume/27.21_dp
!write(*,*) KVECTOR(3*k-2:3*k),BoxSize
!write(*,*) REAL(KVECTOR(3*k-2:3*k),dp)*iBoxSize*TWOPI
!STOP
!END IF
         END DO
      END DO
!STOP
      vd=vd*im_fac/Volume
      force=2._dp*force*im_fac/Volume !eV/Angstrom
      AL%AtomForce=AL%AtomForce+force !eV/Angstrom 
!     write(*,*) AL%AtomForce(301:303)/27.21_dp
!     write(*,*) force(301:303)/27.21_dp
!     write(*,*) AL%AtomForce(601:603)/27.21_dp
!     write(*,*) force(601:603)/27.21_dp
     !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Self-part
      vs=SUM(AtomCharge(1:NAtoms)*AtomCharge(1:NAtoms))
      vs=ALPHA*vs*se_fac
      
      AL%PotentialEnergy=AL%PotentialEnergy+vr+vd+vs
      
!CALL FLUSH(6)
      DEALLOCATE(force)
      DEALLOCATE(QEIKR)
      DEALLOCATE(EIKX,EIKY,EIKZ)
      
   END SUBROUTINE CoulombForceEwald
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ERFC(x)
      IMPLICIT NONE
      !*******************************************************************
      !** APPROXIMATION TO THE COMPLEMENTARY ERROR FUNCTION             **
      !**                                                               **
      !** REFERENCE:                                                    **
      !**                                                               **
      !** ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS,    **
      !**    NATIONAL BUREAU OF STANDARDS, FORMULA 7.1.26               **
      !*******************************************************************

      REAL(dp), PARAMETER :: A1 = 0.254829592, A2 = -0.284496736
      REAL(dp), PARAMETER :: A3 = 1.421413741, A4 = -1.453152027
      REAL(dp), PARAMETER :: A5 = 1.061405429, P  =  0.3275911
      REAL(dp) :: T, X, XSQ, TP, ERFC

!*******************************************************************
      T  = 1.0 / ( 1.0 + P * X )
      XSQ = X * X
      TP = T * ( A1 + T * ( A2 + T * ( A3 + T * ( A4 + T * A5 ) ) ) )
      ERFC = TP * EXP ( -XSQ )
   END FUNCTION ERFC
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION EBUCK(r,A,rho,C,Rc)
      IMPLICIT NONE
      REAL(dp) :: r,A,rho,C,Rc,EBUCK

      EBUCK=0._dp
      IF (rho>1.e-5) EBUCK=A*EXP(-r/rho)-C/r**6 +C/Rc**6  -A*EXP(-Rc/rho) + &
         (Rc/20._dp)*(1._dp-(r/Rc)**20)*(6._dp*C/Rc**7-A*EXP(-Rc/rho)/rho)
         
   END FUNCTION EBUCK
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION FBUCK(r,A,rho,C,Rc)
      IMPLICIT NONE
      REAL(dp) :: r,A,rho,C,Rc,FBUCK

      FBUCK=0._dp
      IF (rho>1.e-5) &
        FBUCK=A*EXP(-r/rho)/(r*rho)-6._dp*C/r**8 !-((r/Rc)**19)*(A*EXP(-Rc/rho)/rho-6._dp*C/Rc**7)/r
   END FUNCTION FBUCK
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SPRING(r12,r13,drmag12,drmag13,K,theta0,ESpring,FSPRING1,FSPRING2,FSPRING3)
   !Subroutine finds the spring energy and the forces (3-body term)
   !FSPRING1 is the side atom in r1
   !FSPRING3 is the side atom in r2
   !FSPRING2 is the center atom
      IMPLICIT NONE
      REAL(dp), DIMENSION(3) :: r12,r13,FSPRING1,FSPRING2,FSPRING3
      REAL(dp) :: ESPRING,drmag12,drmag13, N, D, A, K, theta0, theta, var, D3, D2
      REAL(dp) :: dot_prod

      
      IF (K==0._dp) RETURN

      dot_prod = r13(1)*r12(1) + r13(2)*r12(2) + r13(3)*r12(3)
      D = drmag12*drmag13
      D2 = drmag12*drmag12*D
      D3 = drmag13*drmag13*D
      A = dot_prod/D
      theta = ACOS(A)
      var = -K*(theta - theta0)/(SQRT(1._dp-(dot_prod*dot_prod/D/D)))

      ESPRING = 0.5_dp*K*((theta - theta0)*(theta - theta0)) 
! for atom i (center atom)
      FSPRING1=var*( (r12+r13)/D + (r12/D2 + r13/D3)*dot_prod )
! for atom k    
      FSPRING3=var*( r13/D - r12*dot_prod/D2 )
! for atom j     
      FSPRING2=var*( r12/D - r13*dot_prod/D3 )
   
   END SUBROUTINE SPRING
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   !OpenMP implementation
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CoulombForceEwald_OpenMP(AL,NAtoms,errorstatus)
      !compute coulomb energy and force terms using Ewald summation
      IMPLICIT NONE
      
      REAL(dp) :: theta0=1.91061193_dp,KSpring=1.322626_dp
      
      REAL(dp), PARAMETER :: re_fac=14.399644549472361_dp !real term - conversion factor to get qi qj/(4pi eps0) r in eV [obtained as 1.6021765314e-19/4/pi/8.854187817e-12/1e-10]
      REAL(dp), PARAMETER :: im_fac=1.809512701237067e2_dp !fourier term - conversion factor [obtained as 1.6021765314e-19/1e-10/8.854187817e-12]
      REAL(dp), PARAMETER :: se_fac=-8.124129461602530_dp !self term - conversion factor [obtained as -1.6021765314e-19/8.854187817e-12/4/pi^1.5/1e-10]
      INTEGER :: ListType,errorstatus
      TYPE(SystemContainer), POINTER :: AL
      TYPE(EwaldParameters), POINTER :: Ew
      TYPE(LinkedListContainer), POINTER :: LL
      TYPE(VerletListContainer), POINTER :: VL
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomCharge,KVEC,force,force1,energy,VLdrmag,VLdr
      INTEGER, DIMENSION(:), POINTER :: LLList,LLListRange,VLList,VLListRange,KVector,AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: VLListDomainAtom
      COMPLEX(dp), DIMENSION(:,:), POINTER :: EIKX,EIKY,EIKZ
      COMPLEX(dp), DIMENSION(:), POINTER :: QEIKR
      COMPLEX(dp)  :: EIKR_SUM
      INTEGER :: r0,ra !used with VL
      INTEGER :: NAtoms,TOTk,i,j,k,KX,KY,KZ,KSQ,KMAX,maincellindx,neighcellindx !,TOTr
      INTEGER :: c1(3),c2(3),NCell(3),nx,ny,nz,c1start,c1end,c2start,c2end,ivec(3)
      INTEGER :: MaxAtomPerCell,MaxAtomPerAtom,atomindx1,atomindx2,c1pos,c2pos
      REAL(dp) :: z1,z2,r1(3),r2(3),dr(3),drbox(3),Volume,erfckrij
      REAL(dp) :: iLx2pi,iLy2pi,iLz2pi,krij,kmag2,KSCALED(3),QCOSKR,QSINKR
      REAL(dp) :: vij,fij(3),vd,vs,vk,vr,RX,RY,RZ,drmag,fac,ReCutOff,v3
      REAL(dp) :: ALPHA,B,EIKR,BoxSize(3),iBoxSize(3),b1,fmax,fbu(3),BuckCutOff
      REAL(dp) :: vbuck
      INTEGER, PARAMETER :: chunksize=100 !take 100 at a time --
      INTEGER :: s1,s2,s3
      REAL(dp) :: eb,fb,qsum
      
      !Defined for LiFePO4
      INTEGER :: atomindx3
      REAL(dp) :: dr2(3),dr3(3),drmag2,ESpring,FSpring1(3),FSpring2(3),FSpring3(3)
      
      errorstatus=0
      
      AtomCharge=>AL%AtomCharge
      AtomCoord=>AL%AtomCoord
      AtomSpecies=>AL%AtomSpecies
      Ew=>AL%Potential%Coulomb%Ewald
      LL=>AL%LL
      LLList=>LL%List
      LLListRange=>LL%ListRange
      MaxAtomPerCell=LL%MaxAtomPerCell
      VL=>AL%VL
      VLList=>VL%List
      VLListRange=>VL%ListRange
      VLdrmag=>VL%drmag
      VLdr=>VL%dr
      VLListDomainAtom=>VL%ListDomainAtom
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      
      !NAtoms=AL%NAtoms
      IF (NAtoms/=Ew%NAtoms) THEN
         WRITE(6,*) "Setting up Ewald ..."
         CALL EwaldSetUp(AL)
      END IF
      ALPHA=Ew%ALPHA
      ReCutOff=Ew%ReCutOff
      !BuckCutOff=AL%Potential%MaxPotentialCutoff(4)
      !BuckCutOff=maxrcutbuck
      
      qsum=SUM(AtomCharge(1:NAtoms))
      IF (qsum/=0._dp) THEN
         WRITE(6,*) "$Err>> Net charge is ",qsum
         STOP
      END IF
      
      ALLOCATE(force(3*NAtoms))
      ALLOCATE(energy(NAtoms))
      force=0._dp
      energy=0._dp
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Real space part
      vr=0._dp
      vbuck=0._dp
      b1=2._dp*ALPHA/SQRT(PI)
      
      fmax=0._dp
      IF (AL%VL%ListType==1) THEN !full list

!OMP parallelization at level of atomindx1
!atomindx1 is PRIVATE for each thread

!$OMP PARALLEL

!$OMP DO SCHEDULE(DYNAMIC,chunksize)
         DO atomindx1=1,NAtoms
         
            CALL ForceReal(atomindx1,energy,force,AtomCoord,AtomCharge,AtomSpecies,VLdr, &
               VLdrmag,VLList,VLListRange,re_fac,b1,ALPHA,ReCutOff,MaxAtomPerAtom)
            
         END DO
!$OMP END DO NOWAIT

!$OMP BARRIER
!$OMP END PARALLEL

         !vr=vr*0.5_dp+v3 !for full VL neighbor list
         !vbuck=vbuck*0.5_dp
         
      ELSE !VL is half list
         
         WRITE(6,*) "$Err>> Half-list has not been checked"
         STOP
         
      END IF
      
      AL%AtomForce(1:3*NAtoms)=AL%AtomForce(1:3*NAtoms)-force(1:3*NAtoms) !eV/Angstrom
      vr=SUM(energy(1:NAtoms))
      DEALLOCATE(energy)
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Reciprocal space part -- atomic coordinates have to be in periodic box
      !CALL PBC(AL)
      
      force=0._dp
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      Volume=PRODUCT(BoxSize) !in Angstrom cube
      iBoxSize=1._dp/BoxSize
      !EIKX=>Ew%EIKX
      !EIKY=>Ew%EIKY
      !EIKZ=>Ew%EIKZ
      KMAX=Ew%ImCutOff
!      ALLOCATE(QEIKR(NAtoms))
      ALLOCATE(EIKX(1:NAtoms,0:KMAX))
      ALLOCATE(EIKY(1:NAtoms,-KMAX:KMAX))
      ALLOCATE(EIKZ(1:NAtoms,-KMAX:KMAX))

      iLx2pi=TWOPI*iBoxSize(1)
      iLy2pi=TWOPI*iBoxSize(2)
      iLz2pi=TWOPI*iBoxSize(3)
      
      !WARNING
      !23/12/2012: gfortran appears to be buggy
      !when EIKX is being initialized, it made AL%Potential%SW to be associated
      !this did not happen with intel compiler
      DO atomindx1=1,NAtoms !effect of bug was witnessed here
         EIKX(atomindx1,0)=CMPLX(1.d0,0.d0)
         EIKY(atomindx1,0)=CMPLX(1.d0,0.d0)
         EIKZ(atomindx1,0)=CMPLX(1.d0,0.d0)
         
         j=3*atomindx1-2; RX=AtomCoord(j)*iLx2pi
         EIKX(atomindx1,1)=CMPLX(COS(RX),SIN(RX))
         j=j+1; RY=AtomCoord(j)*iLy2pi
         EIKY(atomindx1,1)=CMPLX(COS(RY),SIN(RY))
         j=j+1; RZ=AtomCoord(j)*iLz2pi
         EIKZ(atomindx1,1)=CMPLX(COS(RZ),SIN(RZ))
      
         EIKY(atomindx1,-1)=CONJG(EIKY(atomindx1,1))
         EIKZ(atomindx1,-1)=CONJG(EIKZ(atomindx1,1))
      END DO
      
      !set up EIKX(1:NAtoms,k) etc.
      DO k=2,KMAX
         EIKX(:,k)=EIKX(:,k-1)*EIKX(:,1)
         EIKY(:,k)=EIKY(:,k-1)*EIKY(:,1)
         EIKY(:,-k)=CONJG(EIKY(:,k))
         EIKZ(:,k)=EIKZ(:,k-1)*EIKZ(:,1)
         EIKZ(:,-k)=CONJG(EIKZ(:,k))
      END DO

      vd=0._dp
      TOTk=Ew%TOTk
      KVector=>Ew%KVector
      KVEC=>Ew%KVEC

!      DO k=1,TOTk
!         KX=KVector(3*k-2)
!         KY=KVector(3*k-1)
!         KZ=KVector(3*k  )
!         fac=2._dp
!         IF (KX==0) fac=1._dp
!         QEIKR(1:NAtoms)=AtomCharge(1:NAtoms)*EIKX(1:NAtoms,KX)*EIKY(1:NAtoms,KY)*EIKZ(1:NAtoms,KZ) !dimension is natoms
!         EIKR_SUM=SUM(QEIKR(1:NAtoms))
!         QCOSKR=REAL(EIKR_SUM)
!         QSINKR=AIMAG(EIKR_SUM)
!         vd=vd+KVEC(k)*(CONJG(EIKR_SUM)*EIKR_SUM) !*fac
!         KSCALED=REAL(KVECTOR(3*k-2:3*k),dp)*KVEC(k)*iBoxSize*TWOPI
!         DO j=1,NAtoms
!            force(3*j-2:3*j)=force(3*j-2:3*j)+ & 
!             (AIMAG(QEIKR(j))*QCOSKR-REAL(QEIKR(j))*QSINKR)*KSCALED
!         END DO
!      END DO

!xxxxxxxxxxxxxOpenMP attempt begins here xxxxxx
      ALLOCATE(energy(TOTk))
      energy=0._dp
!$OMP PARALLEL PRIVATE(force1,QEIKR)
      ALLOCATE(force1(3*NAtoms))
      ALLOCATE(QEIKR(NAtoms))
      force1=0._dp
!$OMP DO SCHEDULE(DYNAMIC,chunksize)
      DO k=1,TOTk
         force1(1:3*NAtoms)=force1(1:3*NAtoms)+ ForceImag1(k,NAtoms,AtomCharge,KVEC, &
           KVector,energy,EIKX,EIKY,EIKZ,QEIKR,iBoxSize,TWOPI)
      END DO
!$OMP END DO
!$OMP CRITICAL
      force=force+force1
!$OMP END CRITICAL
!$OMP BARRIER
      DEALLOCATE(QEIKR)
      DEALLOCATE(force1)
!$OMP END PARALLEL

!!!$OMP BARRIER
      vd=SUM(energy(1:TOTk))
!      write(*,*) "a1"; CALL FLUSH(6)
!stop

      DEALLOCATE(energy)
!xxxxxx OpenMP attempt ends here xxxx

      vd=vd*im_fac/Volume
      force=2._dp*force*im_fac/Volume !eV/Angstrom
      AL%AtomForce=AL%AtomForce+force !eV/Angstrom 
     !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Self-part
      vs=SUM(AtomCharge(1:NAtoms)*AtomCharge(1:NAtoms))
      vs=ALPHA*vs*se_fac
      
      AL%PotentialEnergy=AL%PotentialEnergy+vr+vd+vs
      
!CALL FLUSH(6)
   
      DEALLOCATE(force)
      DEALLOCATE(EIKX,EIKY,EIKZ)
      !DEALLOCATE(QEIKR)  
   END SUBROUTINE CoulombForceEwald_OpenMP
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ForceReal(atomindx1,energy,force,AtomCoord,AtomCharge,AtomSpecies,VLdr, &
      VLdrmag,VLList,VLListRange,re_fac,b1,ALPHA,ReCutOff,MaxAtomPerAtom)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), POINTER :: energy,force,AtomCoord,AtomCharge,VLdr,VLdrmag
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies,VLList,VLListRange
      REAL(dp) :: z1,z2,vij,v3,eb,fb,vbuck,vr,re_fac,dr(3),drmag,ReCutOff,erfckrij,fij(3)
      REAL(dp) :: ALPHA,BuckCutOff,b1,krij
      INTEGER :: s1,s2,atomindx1,atomindx2,ra,r0,j,MaxAtomPerAtom
         
      !fbu=0._dp
      
      z1=AtomCharge(atomindx1)
      s1=AtomSpecies(atomindx1)
      ra=(atomindx1-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(atomindx1) !range for ith atom
      
      vij=0._dp
      v3=0._dp
      vr=0._dp
      
      DO j=ra+1,ra+r0
         atomindx2=VLList(j)
         
         !According to MOLDY atomindx1 /= atomindx2
         !We have not implemented this condition
         
         dr=VLdr(3*j-2:3*j)
         drmag=VLdrmag(j)
         z2=AtomCharge(atomindx2)
         s2=AtomSpecies(atomindx2)

         !spring forces specifically for O-P-O in LiFePO4
         !IF (s1==3 .AND. s2==4 .AND. drmag<1.8_dp) THEN
         !for LiFePO4 the average spacing between O and P is 1.6 Ang
         !   DO k=j+1,ra+r0 !that should give unique triplets involving at least 1 O and 1 P
         !      atomindx3=VLList(k)
         !      s3=AtomSpecies(atomindx3)
         !      dr3=VLdr(3*k-2:3*k)
         !      drmag2=VLdrmag(k)
         !      IF (s2==4 .AND. drmag2<1.8_dp) THEN
         !         CALL SPRING(dr,dr3,drmag,drmag2,KSpring,theta0,ESpring,FSpring1,FSpring2,FSpring3)
         !         v3=v3+ESpring
         !        force(3*atomindx2-2:3*atomindx2)=force(3*atomindx2-2:3*atomindx2)+FSpring2 !atom j
         !         force(3*atomindx1-2:3*atomindx1)=force(3*atomindx1-2:3*atomindx1)+FSpring1 !center atom
         !         force(3*atomindx3-2:3*atomindx3)=force(3*atomindx3-2:3*atomindx3)+FSpring3 !atiom k
         !      END IF
         !   END DO
         !END IF
               
         BuckCutOff=rcutbuck(s1,s2)
         !write(*,*) ".....",s1,s2,ABuck(s1,s2),rhoBuck(s1,s2),CBuck(s1,s2),ReCutOff
         eb=EBUCK(drmag,ABuck(s1,s2),rhoBuck(s1,s2),CBuck(s1,s2),BuckCutOff)
         fb=FBUCK(drmag,ABuck(s1,s2),rhoBuck(s1,s2),CBuck(s1,s2),BuckCutOff)
         IF (drmag<=ReCutOff) THEN
            !write(*,*) atomindx2,dr,drmag
            krij=ALPHA*drmag
            erfckrij=ERFC(krij)/drmag
            vij=z1*z2*erfckrij*re_fac+eb 
            fij=z1*z2*(erfckrij+b1*EXP(-krij*krij))*dr/drmag/drmag*re_fac+fb*dr
         ELSE
            vij=eb
            fij=-fb*dr
         END IF
         !vbuck=vbuck+eb
         !write(*,*) drmag,z1*z2*(erfckrij+b1*EXP(-krij*krij))*dr/drmag/drmag*re_fac
         !fbu=fbu-fb*dr
         force(3*atomindx1-2:3*atomindx1)=force(3*atomindx1-2:3*atomindx1)+fij
         !write(*,*) "Ag ",dr
         !write(UNIT=6,FMT='(2I2,2I5,6ES15.5)') s1,s2,atomindx1,atomindx2,dr,-fb*dr
         vr=vr+vij
         !IF (atomindx1==1) THEN
         !   write(*,*)atomindx2,z2,krij,vij
         !END IF
      END DO
      
      energy(atomindx1)=vr*0.5_dp+v3
   END SUBROUTINE ForceReal
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ForceImag1(k,NAtoms,AtomCharge,KVEC,KVector,energy, &
      EIKX,EIKY,EIKZ,QEIKR,iBoxSize,TWOPI)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), POINTER :: energy,AtomCharge,KVEC
      INTEGER, DIMENSION(:), POINTER :: KVector
      COMPLEX(dp), DIMENSION(:,:), POINTER :: EIKX,EIKY,EIKZ
      COMPLEX(dp), DIMENSION(:), POINTER :: QEIKR
      COMPLEX(dp) :: EIKR_SUM
      INTEGER :: KX,KY,KZ,NAtoms,j,k,atomindx1
      REAL(dp) :: vd,ForceImag1(3*NAtoms),iBoxSize(3),TWOPI,KSCALED(3),QCOSKR,QSINKR
      
      KX=KVector(3*k-2)
      KY=KVector(3*k-1)
      KZ=KVector(3*k  )
      !ALLOCATE(QEIKR(NAtoms))
      QEIKR(1:NAtoms)=AtomCharge(1:NAtoms)* &
         EIKX(1:NAtoms,KX)*EIKY(1:NAtoms,KY)*EIKZ(1:NAtoms,KZ) !dimension is natoms
      EIKR_SUM=SUM(QEIKR(1:NAtoms))
      vd=KVEC(k)*(CONJG(EIKR_SUM)*EIKR_SUM) !*fac
      energy(k)=vd
      QCOSKR=REAL(EIKR_SUM)
      QSINKR=AIMAG(EIKR_SUM)
      KSCALED=REAL(KVECTOR(3*k-2:3*k),dp)*KVEC(k)*iBoxSize*TWOPI
      DO atomindx1=1,NAtoms
         ForceImag1(3*atomindx1-2:3*atomindx1)= &
            (AIMAG(QEIKR(atomindx1))*QCOSKR-REAL(QEIKR(atomindx1))*QSINKR)*KSCALED
      END DO
      !DEALLOCATE(QEIKR)
   END FUNCTION ForceImag1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE CoulombSubroutines
