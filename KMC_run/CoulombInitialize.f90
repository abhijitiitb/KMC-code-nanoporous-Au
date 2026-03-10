MODULE CoulombInitialize
   !initialize parameters needed by Ewald and Coloumb interactions
   USE VARIABLE_TYPE
   USE NeighborList
   IMPLICIT NONE
   REAL(dp), DIMENSION(5,5) :: ABuck=0._dp,rhoBuck=0._dp,CBuck=0._dp,rcutbuck=0._dp
   REAL(dp) :: maxrcutbuck
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WolfSetUp(AL)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      
      AL%Potential%MaxPotentialCutOff(3)=MAX(AL%Potential%MaxPotentialCutOff(3),AL%Potential%Coulomb%Wolf%ReCutOff)
      AL%Potential%VLBuffer=MAXVAL(AL%Potential%MaxPotentialCutOff)/50._dp
      
      WRITE(6,*) "CoulombWolf>> Setting up Wolf potential ..."
   END SUBROUTINE WolfSetUp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PPPMSetUp(AL)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
   END SUBROUTINE PPPMSetUp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE FFMSetUp(AL)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
   END SUBROUTINE FFMSetUp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EwaldSetUp(AL)
      IMPLICIT NONE
      INTEGER, PARAMETER :: MAXR=2500
      TYPE(SystemContainer), POINTER :: AL
      TYPE(InteracPotential), POINTER :: Potential
      TYPE(EwaldParameters), POINTER :: Ew
      INTEGER, DIMENSION(:), POINTER :: RVEC
      INTEGER :: NAtoms,CellSize(3),TOTr,ix,iy,iz,NCellsx,NCellsy,NCellsz
      CHARACTER(len=3) :: SpeciesName
      REAL(dp) :: p,Volume,ReCutOff2,lsq,distancex,distancey,distancez
      LOGICAL :: accept
      
      WRITE(*,*) "Coulomb>> Setting up Ewald potential ..."
      
      Potential=>AL%Potential
      
      NAtoms=AL%NAtoms
      Volume=PRODUCT(AL%BoxSize(1:3)-AL%BoxSize(4:6)) !in Angstrom cube
      
      IF (.NOT. ASSOCIATED(Potential)) THEN
         WRITE(*,*) "$Err>> The variable Potential in not associated in EwaldSetup"
         STOP
      END IF
      
      !in case the user has not manually setup these variables
      IF (.NOT. ASSOCIATED(Potential%Coulomb)) ALLOCATE(Potential%Coulomb)
      IF (.NOT. ASSOCIATED(Potential%Coulomb%Ewald)) ALLOCATE(Potential%Coulomb%Ewald)
      
      Ew=>Potential%Coulomb%Ewald
      
      Ew%NAtoms=NAtoms  !this tells us for how many atoms Ew was initialized
      IF (Ew%TimeRatio<=0._dp) Ew%TimeRatio=5.5_dp !from Moldy user manual - ratio of tR/tF
      IF (Ew%ALPHA<=0._dp) Ew%ALPHA=SQRT(PI)*(Ew%TimeRatio*REAL(NAtoms,dp)/Volume/Volume)**(1._dp/6._dp)
      IF (Ew%Precision<=0._dp) Ew%Precision=5.2e-5_dp
      p=LOG(1./Ew%Precision)
      IF (Ew%ReCutOff<=0._dp) Ew%ReCutOff=SQRT(p)/Ew%alpha !cutoff for real space
      IF (Ew%ImCutOff<=0) Ew%ImCutOff=CEILING(2._dp*SQRT(p)*Ew%ALPHA) !imaginary part
      
      WRITE(6,*) "Ewald parameter values in use"
      WRITE(6,*) "============================="
      WRITE(6,*) "Alpha:",Ew%ALPHA
      WRITE(6,*) "Time ratio:",Ew%TimeRatio
      WRITE(6,*) "Real space cutoff:",Ew%ReCutOff
      WRITE(6,*) "Imaginary space cutoff:",Ew%ImCutOff
      WRITE(6,*) "Precision:",Ew%Precision
      
      AL%Potential%MaxPotentialCutOff(3)=MAX(AL%Potential%MaxPotentialCutOff(3),Ew%ReCutOff)
      AL%Potential%VLBuffer=MAXVAL(AL%Potential%MaxPotentialCutOff)/100._dp
      
      CALL EwImSetUp(AL) !setup interpolation arrays for Ewald
      
      !Now set-up list of cells to search for while evaluating real term -- 
      !will be required if linked list is used in place of VL
      !it is assumed that the linked-list has already been set up
      !IF (.NOT. ASSOCIATED(AL%LL)) CALL AddVerletList(AL,NAtoms=NAtoms)
      !CellSize=AL%LL%CellSize
      !NCellsx=CEILING(Ew%ReCutOff/CellSize(1))
      !NCellsy=CEILING(Ew%ReCutOff/CellSize(2))
      !NCellsz=CEILING(Ew%ReCutOff/CellSize(3))
      !ReCutOff2=Ew%ReCutOff*Ew%ReCutOff
      
      !the ix,iy,iz that are accepted lie within the cutoff
      !use them by adding to the central cell for finding the neighbor cell
      !IF (.NOT. ASSOCIATED(Ew%RVEC)) ALLOCATE(Ew%RVEC(3*MAXR))
      !RVEC=>Ew%RVEC
      !TOTr=0
      !DO ix=-NCellsx,NCellsx
      !   distancex=REAL(ix,dp)*CellSize(1)
      !   DO iy=-NCellsy,NCellsy
      !      distancey=REAL(iy,dp)*CellSize(2)
      !      DO iz=-NCellsz,NCellsz
      !         distancez=REAL(iz,dp)*CellSize(3)
      !         lsq=distancex*distancex+distancey*distancey+distancez*distancez
      !         IF (lsq<ReCutOff2) THEN !this neighbor cell should be included in the Re space calculation
      !            TOTr=TOTr+1
      !            IF (TOTr>MAXR) THEN
      !               WRITE(8,*) "$Err>> Increase the size of the array RxVector"
      !               STOP
      !            END IF
      !            RVEC(3*TOTr-2)=ix
      !            RVEC(3*TOTr-1)=iy
      !            RVEC(3*TOTr  )=iz
      !         END IF
      !      END DO
      !   END DO
      !END DO
      
      !Ew%TOTr=TOTr

      !Buckingham parameters for MgO
      accept=.FALSE.
      IF (accept) THEN
      
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(1))
         IF (speciesname(1:2)/="Mg") THEN
            WRITE(6,*) "Err>> Error in Buckingham potential parameters, expected species 1 to be Mg"
            STOP
         END IF
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(2))
         IF (speciesname(1:1)/="O") THEN
            WRITE(6,*) "Err>> Error in Buckingham potential parameters, expected species 2 to be O"
            STOP
         END IF
         
         ABuck=0._dp
         ABuck(1,2)=1428.5_dp
         ABuck(2,1)=ABuck(1,2)
         ABuck(2,2)=22764._dp
    
         rhoBuck=0._dp
         rhoBuck(1,2)=0.2945_dp
         rhoBuck(2,1)=rhoBuck(1,2)
         rhoBuck(2,2)=0.149_dp

         CBuck=0._dp
         CBuck(2,2)=27.88_dp
         
         write(*,*) "rcutbuck not set"
         stop
      END IF

      !Buckingham parameters for YSZ (Y=1, Zr=2, O=3)
      accept=.TRUE.
      IF (accept) THEN
      
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(1))
         IF (speciesname(1:1)/="Y") THEN
            WRITE(6,*) "Err>> Error in Buckingham potential parameters, expected species 1 to be Y"
            STOP
         END IF
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(2))
         IF (speciesname(1:2)/="Zr") THEN
            WRITE(6,*) "Err>> Error in Buckingham potential parameters, expected species 2 to be Zr"
            STOP
         END IF
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(3))
         IF (speciesname(1:1)/="O") THEN
            WRITE(6,*) "Err>> Error in Buckingham potential parameters, expected species 3 to be O"
            STOP
         END IF
         
         ABuck=0.
         !ABuck(1,3)=1325.6_dp !Y-O
         ABuck(1,3)=1345.1_dp !Y-O
         !ABuck(2,3)=1024.6_dp !Zr-O
         ABuck(2,3)=985.87_dp !Zr-O
         !ABuck(3,3)=22764.3_dp !O-O
         ABuck(3,3)=22764.0_dp !O-O
         ABuck(3,1)=ABuck(1,3)
         ABuck(3,2)=ABuck(2,3)
    
         rhoBuck=0._dp
         !rhoBuck(1,3)=0.3461_dp !Y-O
         rhoBuck(1,3)=0.3491_dp !Y-O
         !rhoBuck(2,3)=0.376_dp !Zr-O
         rhoBuck(2,3)=0.376_dp !Zr-O
         rhoBuck(3,3)=0.149_dp
         rhoBuck(3,1)=rhoBuck(1,3)
         rhoBuck(3,2)=rhoBuck(2,3)

         CBuck=0._dp
         !CBuck(3,3)=27.89_dp
         CBuck(3,3)=27.88_dp
         
         rcutbuck=0._dp
         rcutbuck(1,3)=12._dp !Y-O
         rcutbuck(3,1)=12._dp
         rcutbuck(2,3)=10._dp !Zr-O
         rcutbuck(3,2)=10._dp
         rcutbuck(3,3)=10._dp !O-O
         maxrcutbuck=12._dp
      END IF
      
      !Buckingham parameters for LiFePO4 (Li=1, Fe=2, P=3, O=4)
      accept=.FALSE.
      IF (accept) THEN
      
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(1))
         IF (speciesname(1:2)/="Li") THEN
            WRITE(6,*) "Err>> Error in Buckingham potential parameters, expected species 1 to be Li"
            STOP
         END IF
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(2))
         IF (speciesname(1:2)/="Fe") THEN
            WRITE(6,*) "Err>> Error in Buckingham potential parameters, expected species 2 to be Fe"
            STOP
         END IF
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(3))
         IF (speciesname(1:1)/="P") THEN
            WRITE(6,*) "Err>> Error in Buckingham potential parameters, expected species 3 to be P"
            STOP
         END IF
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(4))
         IF (speciesname(1:1)/="O") THEN
            WRITE(6,*) "Err>> Error in Buckingham potential parameters, expected species 1 to be O"
            STOP
         END IF
         
         ABuck=0._dp
         ABuck(1,4)=632.1018_dp
         ABuck(2,4)=1105.2409_dp
         ABuck(3,4)=897.2648_dp
         ABuck(4,4)=22764.3_dp
         ABuck(4,1)=ABuck(1,4) !Li-O
         ABuck(4,2)=ABuck(2,4) !Fe-O
         ABuck(4,3)=ABuck(3,4) !P-O

         rhoBuck=0._dp
         rhoBuck(1,4)=0.2906_dp !Li-O
         rhoBuck(2,4)=0.3106_dp !Fe-O
         rhoBuck(3,4)=0.3577_dp !P-O
         rhoBuck(4,4)=0.149_dp  !O-O
         rhoBuck(4,1)=rhoBuck(1,4)
         rhoBuck(4,2)=rhoBuck(2,4)
         rhoBuck(4,3)=rhoBuck(3,4)

         CBuck=0._dp
         CBuck(4,4)=44.53_dp
         
         write(*,*) "rcutbuck not set"
         stop
      END IF

      WRITE(6,*) "===============i=================="
      WRITE(6,*) "Buckingham interaction parameters"
      WRITE(6,*) "================================="
      WRITE(6,*) "A:"
      WRITE(6,*) ABuck(1,:)
      WRITE(6,*) ABuck(2,:)
      WRITE(6,*) ABuck(3,:)
      WRITE(6,*) ABuck(4,:)
      WRITE(6,*) ABuck(5,:)
      WRITE(6,*) "rho:"
      WRITE(6,*) rhoBuck(1,:)
      WRITE(6,*) rhoBuck(2,:)
      WRITE(6,*) rhoBuck(3,:)
      WRITE(6,*) rhoBuck(4,:)
      WRITE(6,*) rhoBuck(5,:)
      WRITE(6,*) "C:"
      WRITE(6,*) CBuck(1,:)
      WRITE(6,*) CBuck(2,:)
      WRITE(6,*) CBuck(3,:)
      WRITE(6,*) CBuck(4,:)
      WRITE(6,*) CBuck(5,:)
      WRITE(6,*) "================================="
   END SUBROUTINE EwaldSetUp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
   SUBROUTINE EwImSetUp(AL)
      IMPLICIT NONE
      INTEGER, PARAMETER :: MAXK=10000
      TYPE(SystemContainer), POINTER :: AL
      TYPE(EwaldParameters), POINTER :: Ew
      INTEGER :: KSQMAX=0 !KSQMAX=27 used in A&T
      REAL(dp), DIMENSION(:), POINTER :: r
      REAL(dp) :: ALPHA,B,ReKX,ReKY,ReKZ,RKSQ
      REAL(dp) :: Lx,Ly,Lz
      INTEGER, DIMENSION(:), POINTER :: KVector
      REAL(dp), DIMENSION(:), POINTER :: KVEC
      INTEGER :: NAtoms,KX,KY,KZ,KSQ,TOTk,KMAX,KYMIN,KZMIN
      
      Ew=>AL%Potential%Coulomb%Ewald
      ALPHA=Ew%ALPHA
      B=0.25_dp/(ALPHA*ALPHA) !used in imaginary part
      KMAX=Ew%ImCutOff !A&T suggests a value of 5
      KSQMAX=KMAX*KMAX+2 !I think this should work well
      NAtoms=Ew%NAtoms
      
      Lx=AL%BoxSize(1)-AL%BoxSize(4)
      Ly=AL%BoxSize(2)-AL%BoxSize(5)
      Lz=AL%BoxSize(3)-AL%BoxSize(6)
      
      
      IF (.NOT. ASSOCIATED(Ew%KVEC)) ALLOCATE(Ew%KVEC(MAXK)) !stores the k^2 values
      IF (.NOT. ASSOCIATED(Ew%KVector)) ALLOCATE(Ew%KVector(3*MAXK)) !stores the k vector
      KVEC=>Ew%KVEC
      KVector=>Ew%KVector
      
      !IF (ASSOCIATED(Ew%EIKX)) DEALLOCATE(Ew%EIKX)
      !ALLOCATE(Ew%EIKX(1:NAtoms,0:KMAX))
      !IF (ASSOCIATED(Ew%EIKY)) DEALLOCATE(Ew%EIKY)
      !ALLOCATE(Ew%EIKY(1:NAtoms,-KMAX:KMAX))
      !IF (ASSOCIATED(Ew%EIKZ)) DEALLOCATE(Ew%EIKZ)
      !ALLOCATE(Ew%EIKZ(1:NAtoms,-KMAX:KMAX))
      
      TOTk=0

      KYMIN=0
      KZMIN=1
      DO KX=0,KMAX
         ReKX=TWOPI*REAL(KX,dp)/Lx
         !DO KY=-KMAX,KMAX
         DO KY=KYMIN,KMAX
            ReKY=TWOPI*REAL(KY,dp)/Ly
            !DO KZ=-KMAX,KMAX
            DO KZ=KZMIN,KMAX
               ReKZ=TWOPI*REAL(KZ,dp)/Lz
               KSQ=KX*KX + KY*KY + KZ*KZ
               
               IF (KSQ<=KSQMAX .AND. KSQ/=0) THEN
                  TOTk=TOTk+1
                  IF (TOTk > MAXK) THEN
                     WRITE(6,*) "$Err>> Increase the maximum size of kvector in EwImSetup"
                     STOP
                  END IF
                  RKSQ= ReKX*ReKX + ReKY*ReKY + ReKZ*ReKZ
                  KVEC(TOTk)=EXP(-B*RKSQ) / RKSQ
                  KVector(3*TOTk-2)=KX
                  KVector(3*TOTk-1)=KY
                  KVector(3*TOTk  )=KZ
               END IF
            END DO
            KZMIN=-KMAX
         END DO
         KYMIN=-KMAX
      END DO
      Ew%TOTk=TOTk
   END SUBROUTINE EwImSetUp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE CoulombInitialize
