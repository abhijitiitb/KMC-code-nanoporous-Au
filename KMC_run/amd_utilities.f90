MODULE utilities
   USE VARIABLE_TYPE
   USE db_manipulate
   USE randomgen
   
   INTERFACE PBCdistance
      MODULE PROCEDURE PBCdistance1,PBCdistance2,PBCdistance3
   END INTERFACE
   
   INTERFACE Compare
      MODULE PROCEDURE CompareAL,CompareChOS
   END INTERFACE
   
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   PURE FUNCTION RotationMatrix(e) !e is euler angle
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: e(3)
      REAL(dp) :: c(3),s(3),RotationMatrix(3,3)
      
      c=COS(e)
      s=SIN(e)
      
      !x-y-z frame
      !z-x-z frame
      RotationMatrix(1,1)=c(1)*c(3)-s(1)*c(2)*s(3)
      RotationMatrix(2,1)=s(1)*c(3)+c(1)*c(2)*s(3)
      RotationMatrix(3,1)=s(2)*s(3)
      RotationMatrix(1,2)=-c(1)*s(3)-s(1)*c(2)*c(3)
      RotationMatrix(2,2)=-s(1)*s(3)+c(1)*c(2)*c(3)
      RotationMatrix(3,2)=s(2)*c(3)
      RotationMatrix(1,3)=s(2)*s(1)
      RotationMatrix(2,3)=-s(2)*c(1)
      RotationMatrix(3,3)=c(2)
   END FUNCTION RotationMatrix
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PBC(AL)
   
      !Arrange atoms into period
   
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: i,r1,r2
      REAL(dp) :: xc(3),L(3),c(3)
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord
      
      xc=AL%BoxSize(4:6) !smallest corner
      L=AL%BoxSize(1:3)-AL%BoxSize(4:6) !actual box size
      AtomCoord=>AL%AtomCoord
      r1=-2; r2=0
      
      DO i=1,AL%NAtoms
      
         r1=r1+3; r2=r2+3
         
         c=AtomCoord(r1:r2)-xc
         
         AtomCoord(r1:r2)=c-FLOOR(c/L)*L+xc
         
      END DO
      
   END SUBROUTINE PBC
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   PURE FUNCTION PBCdistance1(AL,atom1,atom2)
   
      !difference vector from atom2
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: atom1,atom2
      REAL(dp), DIMENSION(3) :: r1,r2,vec1,vec2,vec3,PBCdistance1,BoxSize
      TYPE(SystemContainer), POINTER :: AL
      
      IF (atom1==atom2) THEN
         PBCdistance1=0._dp
         RETURN
      END IF
      
      r1=AL%AtomCoord(3*atom1-2:3*atom1)
      r2=AL%AtomCoord(3*atom2-2:3*atom2)
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      
      PBCdistance1=r1-r2 !vector from atom2
      PBCdistance1=PBCdistance1-BoxSize*NINT(PBCdistance1/BoxSize)
      
   END FUNCTION PBCdistance1
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION PBCdistance2(AL1,AL2)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL1,AL2
      REAL(dp), DIMENSION(3*AL1%NAtoms) :: PBCdistance2
      REAL(dp), DIMENSION(3*MaxNAtoms) :: BoxSize
      INTEGER :: NAtoms
      
      IF (AL1%NAtoms/=AL2%NAtoms) THEN
         WRITE(6,*) "$Err>> Size of Atom lists do not match to get distance"
         STOP
      END IF
      IF (ANY(AL1%BoxSize/=AL2%BoxSize)) THEN
         WRITE(6,*) "$Err>> Size of Atom list boxes do not match to get distance"
         STOP
      END IF
      
      NAtoms=AL1%NAtoms
      BoxSize(1:3*NAtoms)=PACK(SPREAD(AL1%BoxSize(1:3)-AL1%BoxSize(4:6),2,NAtoms),.TRUE.)
      PBCDistance2(1:3*NAtoms)=AL1%AtomCoord(1:3*NAtoms)-AL2%AtomCoord(1:3*NAtoms)
      PBCDistance2(1:3*NAtoms)=PBCdistance2(1:3*NAtoms)-BoxSize(1:3*NAtoms)* &
         NINT(PBCdistance2(1:3*NAtoms)/BoxSize(1:3*NAtoms))
   END FUNCTION PBCdistance2
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   PURE FUNCTION PBCdistance3(AL,coord1,coord2)
   
      IMPLICIT NONE
      REAL(dp), DIMENSION(3) :: r1,r2,vec1,vec2,vec3,PBCdistance3,BoxSize
      REAL(dp), DIMENSION(3), INTENT(IN) :: coord1,coord2
      TYPE(SystemContainer), POINTER :: AL
   
      r1=coord1
      r2=coord2
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      
      vec1=r1-r2 !vector from atom2
      PBCdistance3=vec1-BoxSize*NINT(vec1/BoxSize)
      !vec2=BoxSize+vec1 !displace q by Box
      !vec3=vec1-BoxSize !displace q by -Box
      
      !PBCdistance3=vec1
      !WHERE(ABS(vec2)<ABS(PBCdistance3)) PBCdistance3=vec2
      !WHERE(ABS(vec3)<ABS(PBCdistance3)) PBCdistance3=vec3
      
   END FUNCTION PBCdistance3
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION iPBC(i,n)
   !returns a PBC index given original index i and number of cells n
   !this is useful when neighboring cells are searched and the index i
   !can be out of bound
      IMPLICIT NONE
      INTEGER :: i,n,ipbc
      
      iPBC=i-n*FLOOR(REAL(i-1)/REAL(n))
   END FUNCTION iPBC
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CompareAL(AL1,AL2,ByPass1,Tolerance) !true if match
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL1,AL2
      LOGICAL :: CompareAL,ByPass
      LOGICAL, OPTIONAL :: ByPass1
      REAL(dp) :: Maxdx,Tolerance1
      REAL(dp), OPTIONAL :: Tolerance
   
      IF (PRESENT(ByPass1)) THEN
         ByPass=ByPass1
      ELSE
         ByPass=.FALSE.
      END IF
      
      IF (PRESENT(Tolerance)) THEN
         Tolerance1=Tolerance
      ELSE
         Tolerance1=1.0_dp
      END IF
      
      IF (.NOT. ByPass) THEN
         CompareAL=(AL1%NAtoms==AL2%NAtoms); IF (.NOT. CompareAL) RETURN
         !CompareAL=(AL1%NSpecies==AL2%NSpecies); IF (.NOT. CompareAL) RETURN
         CompareAL=(ALL(AL1%AtomSpecies==AL2%AtomSpecies)); IF (.NOT. CompareAL) RETURN
         CompareAL=(ALL(AL1%AtomIsMoving.eqv.AL2%AtomIsMoving)); IF (.NOT. CompareAL) RETURN
         CompareAL=(ALL(AL1%AtomCharge==AL2%AtomCharge)); IF (.NOT. CompareAL) RETURN
      END IF
      
      Maxdx=MAXVAL(PBCDistance(AL1,AL2))
      CompareAL=Maxdx<Tolerance1
   END FUNCTION CompareAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CompareChOS(ChOS1,ChOS2) !true if match
      IMPLICIT NONE
      TYPE(ChOSContainer), POINTER :: ChOS1,ChOS2
      LOGICAL :: CompareChOS
      
      CompareChOS=.FALSE.
   END FUNCTION CompareChOS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   PURE FUNCTION NORM(x,n)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      INTEGER :: i
      REAL(dp), DIMENSION(:), POINTER :: x
      REAL(dp) :: NORM,a
      
      NORM=0._dp
      DO i=1,n
         a=x(i)
         NORM=NORM+a*a
      END DO
      NORM=SQRT(NORM)
   END FUNCTION NORM
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION DET(A)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:,:) :: A
      REAL(dp) :: DET
      
      SELECT CASE(SIZE(A))
      !CASE(1); DET=ABS(A)
      CASE(4); DET=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      CASE(9); DET=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))- &
         A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))+ &
         A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      !CASE DEFAULT; CALL LUDecompose(A)
      END SELECT
   END FUNCTION DET
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION Determinant3(A)
      IMPLICIT NONE
      REAL(dp), DIMENSION(3,3) :: A
      REAL(dp) :: Determinant3
      
      Determinant3=A(1,1)
   END FUNCTION Determinant3
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetMovingAtoms(AL)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: GetMovingAtoms,i,r1,r2
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      
      GetMovingAtoms=0
      r1=-2; r2=0
      AtomIsMoving=>AL%AtomIsMoving
      DO i=1,AL%NAtoms
         r1=r1+3; r2=r2+3
         IF (ANY(AtomIsMoving(r1:r2))) GetMovingAtoms=GetMovingAtoms+1 
      END DO
   END FUNCTION GetMovingAtoms
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetTemperature(AL,constrained)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp) :: GetTemperature
      LOGICAL, OPTIONAL :: constrained
      INTEGER :: NMove(3)
      
      IF (PRESENT(constrained)) THEN
         !Maxwell-Boltzmann velocities are assigned 
         !for N-1 atoms, s.t. the system velocity is
         !zero ... this gives NDOF=N-1 (deg of freedom)
         IF (constrained) THEN
            NMove=AL%NMove-1
         ELSE
            NMove=AL%NMove
         END IF
      ELSE
         NMove=AL%NMove-1
      END IF
      
      !kboltzmann=8.314472_dp/6.02214179e23_dp
      CALL GetKineticEnergy(AL)
      !GetTemperature=2._dp*AL%KineticEnergy/3._dp/ &
      !  kboltzmann/REAL(NAtoms,dp)
      GetTemperature=2._dp*AL%KineticEnergy/ &
        kboltzmann/REAL(SUM(NMove),dp)
      
   END FUNCTION GetTemperature
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetKineticEnergy(AL)
      !unit of kinetic energy is eV
      !unit of velocity is sqrt(eV/amu)
      !unit of mass is amu
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(:), POINTER :: AtomMass,AtomVelocity
      REAL(dp) :: KE,m,v
      INTEGER :: iatom,NAtoms
      
      NAtoms=AL%NAtoms
      AtomMass=>AL%AtomMass
      AtomVelocity=>AL%AtomVelocity
      KE=0._dp
      DO iatom=1,3*NAtoms
         m=AtomMass(iatom)
         v=AtomVelocity(iatom)
         KE=KE+m*v*v
      END DO
      KE=0.5_dp*KE
      
      AL%KineticEnergy=KE
      
   END SUBROUTINE GetKineticEnergy
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetHistogram(x,nbin,xstart,xend)
      !make sure histogram is deleted at the end
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), POINTER :: x
      TYPE(RDFHistogram), POINTER :: GetHistogram
      INTEGER :: nbin,i,j
      REAL(dp), OPTIONAL :: xstart,xend
      REAL(dp) :: dx
      
      IF (.NOT.(PRESENT(xstart) .AND. PRESENT(xend))) THEN
         xstart=MINVAL(x)
         xend=MAXVAL(x)
      END IF
      GetHistogram%xstart=xstart
      GetHistogram%xend=xend
      GetHistogram%nbin=nbin
      CALL MakeSize(GetHistogram%hist,nbin)
      
      dx=(xstart-xend)/REAL(nbin,dp)
      DO i=1,SIZE(x)
         j=CEILING(x(i)/dx)
         GetHistogram%hist(j)=GetHistogram%hist(j)+1
      END DO
   END FUNCTION GetHistogram
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetAtomListMass(AL) 
      !transfers atom info from specieslist array
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: NAtoms,i,indx,atomicno,r1,r2
      REAL(dp) :: mass
      
      NAtoms=AL%NAtoms
      IF (NAtoms==0) THEN
         WRITE(*,*) "$Err>> Atom list is empty - cant get mass"
         STOP
      END IF
      !CALL MakeSize(AL%AtomMass,NAtoms)
      
      r1=-2; r2=0
      DO i=1,NAtoms
         r1=r1+3; r2=r2+3
         indx=AL%AtomSpecies(i)
         IF (indx==0) THEN
            WRITE(6,*) "$Err>> AtomSpecies not initialized in AtomList"
            STOP
         END IF
         atomicno=SpeciesDirectory_global(indx)
         mass=SpeciesList%AtomicMass(atomicno)
         AL%AtomMass(r1:r2)=mass !in amu
      END DO
   END SUBROUTINE GetAtomListMass
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetSpeciesSymbol(AtomicSymbol)
      !returns the species index in SpeciesDirectory corresponding to atomic symbol
      IMPLICIT NONE
      CHARACTER(len=3):: AtomicSymbol
      INTEGER :: i,GetSpeciesSymbol
      
      IF (NSpecies_global==0) THEN
         WRITE(6,*) "$Err>> Number of species is zero"
         STOP
      END IF
      
      IF (.NOT. ASSOCIATED(SpeciesDirectory_global)) THEN
         WRITE(6,*) "$Err>> SpeciesDirectory_global is not initialized"
         STOP
      END IF
      
      IF (ALL(SpeciesDirectory_global==0)) THEN
         WRITE(6,*) "$Err>> SpeciesDirectory_global is not initialized"
         STOP
      END IF
      
      DO i=1,NSpecies_global
         IF (TRIM(SpeciesList%AtomicSymbol(SpeciesDirectory_global(i))) &
            ==TRIM(AtomicSymbol)) THEN
            GetSpeciesSymbol=i
            EXIT
         END IF
      END DO
   END FUNCTION GetSpeciesSymbol
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetMDHistoryStatistics(MDSiml,n1,n2)
      !get mean and std deviation for snapshots from n1 to n2
      IMPLICIT NONE
      INTEGER :: n1,n2,i
      TYPE(MDContainer), POINTER :: MDSiml
      TYPE(SystemContainer), POINTER :: AL,AL1
      REAL(dp), DIMENSION(4) :: rdata=0._dp,rdata2=0._dp !PE,KE,PE+KE,Temperature
      
      IF (MDSiml%HistorySize<n2) THEN
         WRITE(6,*) "$Err>> Number of MD snapshots smaller than range provided"
         WRITE(UNIT=6,FMT='("... # MD snaps present:",I6)') MDSiml%HistorySize
      END IF
      
      AL=>MDSiml%History
      DO i=1,n1-1
         AL=>AL%NextNeigh
      END DO
      
      IF (n1<1 .OR. n2<1) THEN
         WRITE(6,*) "$Err>> range provided to GetMDHistoryStatistics is not a positive number"
         WRITE(UNIT=6,FMT='("... min range:")') n1
         WRITE(UNIT=6,FMT='("... max range:")') n2
         STOP
      END IF
      
      rdata=0._dp
      AL1=>AL
      DO i=n1,n2
         rdata(1)=rdata(1)+AL%PotentialEnergy
         rdata(2)=rdata(2)+AL%KineticEnergy
         rdata(3)=rdata(3)+AL%PotentialEnergy+AL%KineticEnergy
         rdata(4)=rdata(4)+GetTemperature(AL)
         AL=>AL%NextNeigh
      END DO
      rdata=rdata/REAL(n2-n1+1,dp)
      
      rdata2=0._dp
      AL=>AL1
      DO i=n1,n2
         rdata2(1)=rdata2(1)+(AL%PotentialEnergy-rdata(1))*(AL%PotentialEnergy-rdata(1))
         rdata2(2)=rdata2(2)+(AL%KineticEnergy-rdata(2))*(AL%KineticEnergy-rdata(2))
         rdata2(3)=rdata2(3)+(AL%PotentialEnergy+AL%KineticEnergy-rdata(3))*(AL%PotentialEnergy+AL%KineticEnergy-rdata(3))
         rdata2(4)=rdata2(4)+(GetTemperature(AL)-rdata(4))*(GetTemperature(AL)-rdata(4))
         AL=>AL%NextNeigh
      END DO
      rdata2=SQRT(rdata2)
      
      WRITE(UNIT=6,FMT='(">>Number of snapshots (total):",I6)') MDSiml%HistorySize
      WRITE(UNIT=6,FMT='(">>Number of snapshots used for averaging:",I6)') n2-n1+1
      WRITE(UNIT=6,FMT='(">>Potential energy:",ES15.5,"+/-",ES15.5)') rdata(1),rdata2(1)
      WRITE(UNIT=6,FMT='(">>Kinetic energy  :",ES15.5,"+/-",ES15.5)') rdata(2),rdata2(2)
      WRITE(UNIT=6,FMT='(">>Total energy    :",ES15.5,"+/-",ES15.5)') rdata(3),rdata2(3)
      WRITE(UNIT=6,FMT='(">>Temperature     :",ES15.5,"+/-",ES15.5)') rdata(4),rdata2(4)
   END SUBROUTINE GetMDHistoryStatistics
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetMaxwellBoltzmann(AL,Temperature,IsPrint)
      !unit of velocity is sqrt(eV/ang)
   
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), INTENT(IN) :: Temperature
      REAL(sp), DIMENSION(:), ALLOCATABLE :: v
      REAL(sp) :: randno
      INTEGER :: NAtoms,i,indx
      REAL(dp) :: const,vcom(3)
      LOGICAL, OPTIONAL :: IsPrint
      LOGICAL :: IsPrint1
      
      IF (PRESENT(IsPrint)) THEN
		IsPrint1=IsPrint
      ELSE
		IsPrint1=.TRUE.
      END IF
      
      IF (IsPrint1) WRITE(UNIT=*,FMT='("Ini>> Setting up atom velocities at temperature ", &
          F8.3," K")') Temperature

      NAtoms=AL%NAtoms
      IF (NAtoms==0) THEN
         WRITE(*,*) "$Err>> Atom list is empty - cant get velocity`"
         STOP
      END IF
      AL%AtomVelocity=0._dp
      
      ALLOCATE(v(3*NAtoms)) !note v is in single precision
      CALL gasdev_v(v) !this is scaled momentum
      
      const=SQRT(kboltzmann*Temperature)
          !standard deviation without mass
      WHERE (AL%AtomIsMoving) &
        AL%AtomVelocity=const*REAL(v,dp)/SQRT(AL%AtomMass) !final velocity
      DEALLOCATE(v)
          !at this point velocity is in a.u.
          ! to convert to SI units multiply with fvel

      !set net velocity to zero
      vcom=(/SUM(AL%AtomVelocity(1:3*NAtoms-2:3)), &
          SUM(AL%AtomVelocity(2:3*NAtoms-1:3)), &
          SUM(AL%AtomVelocity(3:3*NAtoms:3))/)
    
      WHERE (AL%NMove>0) 
         vcom=vcom/REAL(AL%NMove) !note NMove gives # free dimensions
      ELSEWHERE
         vcom=0._dp
      END WHERE
      
      DO indx=1,NAtoms
         WHERE (AL%AtomIsMoving(3*indx-2:3*indx)) &
           AL%AtomVelocity(3*indx-2:3*indx)= &
           AL%AtomVelocity(3*indx-2:3*indx) - vcom
      END DO ! in reality there is AL%Move-1 degrees of freedom in each dimension
      
      IF (IsPrint1) WRITE(UNIT=6,FMT='("   System temperature set to ",F10.3," K")') GetTemperature(AL)
      
   END SUBROUTINE GetMaxwellBoltzmann
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!    FUNCTION GetUnitNotInUse()
!       IMPLICIT NONE
!       INTEGER :: GetUnitNotInUse
!       
!       DO GetUnitNotInUse=UnitTmpStart,UnitTmpEnd
!       END DO
!    END FUNCTION GetUnitNotInUse
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   FUNCTION GetCenterOfMass(AL)
      !get Center of mass of atoms in AL
      !apply PBC before calling this fn if periodicity present
      
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(3) :: GetCenterOfMass,Mass
      INTEGER :: i,NAtoms,r1,r2
      
      NAtoms=AL%NAtoms
      GetCenterOfMass=0._dp
      Mass=0._dp
      
      r1=-2; r2=0
      DO i=1,NAtoms
         r1=r1+3; r2=r2+3
         Mass=Mass+AL%AtomMass(r1:r2)
         GetCenterOfMass=GetCenterOfMass+AL%AtomMass(r1:r2)*AL%AtomCoord(r1:r2)
      END DO
      GetCenterOfMass=GetCenterOfMass/Mass
   
   END FUNCTION GetCenterOfMass
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !Commented because ETIME not contained in ifort
   !FUNCTION GetCurrentCPUTime(IsNotInitialized)
   !   IMPLICIT NONE
   !   REAL, SAVE :: time_init
   !   REAL :: CPUTime,GetCurrentCPUTime,TARRAY(2)
   !   LOGICAL, INTENT(IN) :: IsNotInitialized
   
!  !    CALL CPU_TIME(CPUTime)
   !   !WALL TIME
   !   CPUTime=ETIME(TARRAY)
   !   IF (IsNotInitialized) THEN
   !      time_init=CPUTime
   !   END IF
   !   GetCurrentCPUTime=(CPUTime-time_init)
   !   RETURN
   ! END FUNCTION GetCurrentCPUTime
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   FUNCTION GetBoundingBox(AL)
      !Returns the maximum and minimum coordinates
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: i
      REAL(dp), DIMENSION(6) :: BBox,GetBoundingBox
      
      BBox(4:6)=1.e10_dp !lower bound
      BBox(1:3)=-1.e10_dp !upper bound
      DO i=1,AL%NAtoms
         BBox(4:6)=MIN(BBox(4:6),AL%AtomCoord(3*i-2:3*i))
         BBox(1:3)=MAX(BBox(1:3),AL%AtomCoord(3*i-2:3*i))
      END DO
      GetBoundingBox=BBox
   END FUNCTION GetBoundingBox
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetCoordination(AL,atomindx,opt,rcut)
      !returns # atoms based on opt
      !opt=0 - count # atoms with cutoff rcut - rcut has to be provided
      !opt=1 - count # atoms with Verlet cutoff, rcut is not needed, LL and VL up to date
      IMPLICIT NONE
      INTEGER :: atomindx,opt
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp) :: GetCoordination,rcut1
      REAL(dp), OPTIONAL :: rcut
      
      GetCoordination=0._dp
      SELECT CASE(opt)
      CASE(0)
         IF (.NOT. PRESENT(rcut)) THEN
            WRITE(*,*) "$Err>> Cutoff should have been provided with option 0 in GetCoordination"
            STOP
         END IF
         !complete this part
      CASE(1) !VL should be up to date
         GetCoordination=AL%VL%ListRange(atomindx)
      CASE DEFAULT
         WRITE(*,*) "$Err>> Option not available for GetCoordination"
         STOP
      END SELECT
   END FUNCTION GetCoordination
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetImmobile(AL,mask)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      LOGICAL, DIMENSION(:), POINTER :: mask
      INTEGER :: i
      
      IF (SIZE(mask)==AL%NAtoms) THEN
         CALL SetImmobile1(AL,mask)
      ELSEIF (SIZE(mask)==3*AL%NAtoms) THEN
         CALL SetImmobile3(AL,mask)
      ELSE
         WRITE(*,*) "$Err>> mask size provided to setimmobile is wierd"
         STOP
      END IF
      
      AL%NMove=0
      DO i=1,AL%NAtoms
         WHERE (AL%AtomIsMoving(3*i-2:3*i)) AL%NMove=AL%NMove+1
      END DO
   END SUBROUTINE SetImmobile
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetImmobile1(AL,mask)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      LOGICAL, DIMENSION(:), POINTER :: mask,mask1
      INTEGER :: i
      
      !when mask is TRUE then atom is immobile
      CALL MakeSize(mask1,3*AL%NAtoms)
      mask1=PACK(SPREAD(.NOT. mask,1,3),.TRUE.)
      !now if mask1 is FALSE atom is immobile
      WHERE(.NOT. mask1) AL%AtomIsMoving=mask1
      DEALLOCATE(mask1)
      
   END SUBROUTINE SetImmobile1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetImmobile3(AL,mask)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      LOGICAL, DIMENSION(:), POINTER :: mask,mask1
      
      !when mask is TRUE then atom is immobile
      CALL MakeSize(mask1,3*AL%NAtoms)
      mask1=.NOT. mask
      WHERE(mask) AL%AtomIsMoving=mask1
   END SUBROUTINE SetImmobile3
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetDriftVelocity(AL,DriftVelocity,mask)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(3) :: DriftVelocity
      REAL(dp), DIMENSION(:), POINTER :: vexpand
      LOGICAL, DIMENSION(:), POINTER, OPTIONAL :: mask !mask is of size AL%NAtoms
      LOGICAL, DIMENSION(:), POINTER :: mask1 !mask is of size AL%NAtoms
      
      CALL MakeSize(mask1,3*AL%NAtoms)
      IF (PRESENT(mask)) THEN
         mask1=mask
      ELSE
         mask1=.TRUE.
      END IF
      
      CALL MakeSize(vexpand,3*AL%NAtoms)
      vexpand=PACK(SPREAD(DriftVelocity,2,AL%NAtoms),.TRUE.)
      
      WHERE (mask1) AL%AtomDriftVelocity=vexpand
      DEALLOCATE(mask1,vexpand)
   END SUBROUTINE SetDriftVelocity
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetMDDetails(MDSiml,HistoryOn,HistoryFreq, &
      HistorySizeMax,dt,Temperature, &
      LangevinCoeff,TemperatureRamp)
      
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      LOGICAL, OPTIONAL :: HistoryOn
      INTEGER, OPTIONAL :: HistoryFreq,HistorySizeMax
      REAL(dp), OPTIONAL :: dt,Temperature,LangevinCoeff,TemperatureRamp
      
      IF (PRESENT(HistoryOn)) THEN
         MDSiml%HistoryOn=HistoryOn
      END IF
      
      !IF (PRESENT(HistoryFreq)) THEN
      !   IF (HistoryFreq<0) THEN
      !      WRITE(6,*) "$Err>> HistoryFreq in MDSiml is being set to negative number"
      !      STOP
      !   END IF
      !   MDSiml%HistoryFreq=HistoryFreq
      !END IF
      
      IF (PRESENT(HistorySizeMax)) THEN
         IF (HistorySizeMax<0) THEN
            WRITE(6,*) "$Err>> HistorySizeMax in MDSiml is being set to negative number"
            STOP
         END IF
         MDSiml%HistorySizeMax=HistorySizeMax
      END IF
      
      IF (PRESENT(dt)) THEN
         IF (dt<=0._dp) THEN
            WRITE(6,*) "$Err>> dt in MDSiml is being set to negative number or zero"
            STOP
         END IF
         MDSiml%dt=dt
      END IF
      
      IF (PRESENT(Temperature)) THEN
         IF (Temperature<=0._dp) THEN
            WRITE(6,*) "$Err>> Temperature in MDSiml is being set to negative number or zero"
            STOP
         END IF
         MDSiml%Temperature=Temperature
      END IF
      
      IF (PRESENT(LangevinCoeff)) THEN
         IF (LangevinCoeff<=0._dp) THEN
            WRITE(6,*) "$Err>> LangevinCoeff in MDSiml is being set to negative number or zero"
            STOP
         END IF
         MDSiml%LangevinCoeff=LangevinCoeff
      END IF
      
      IF (PRESENT(TemperatureRamp)) THEN
         MDSiml%TemperatureRamp=TemperatureRamp
      END IF
      
   END SUBROUTINE SetMDDetails
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchSymbol(AL,AtomicSymbol)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      CHARACTER(len=3):: AtomicSymbol
      INTEGER :: i,SearchSymbol
      
      DO i=1,NSpecies_global
         IF (TRIM(SpeciesList%AtomicSymbol(SpeciesDirectory_global(i))) &
            ==TRIM(AtomicSymbol)) THEN
            SearchSymbol=i
            EXIT
         END IF
      END DO
   END FUNCTION SearchSymbol
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION INT2CHAR(Natural)
      IMPLICIT NONE
      INTEGER :: Natural,DecimalInFocus(NINT2CHAR),DecimalPos
      CHARACTER(len=NINT2CHAR)::INT2CHAR
      CHARACTER :: Ch(NINT2CHAR)
      INTEGER :: i

      IF (Natural<=0) THEN
         WRITE(*,*) "Err>> Only natural numbers allowed in INT2CHAR not ",Natural
         STOP
      END IF

      Ch=' '
      DecimalPos=INT(LOG10(REAL(Natural))) 
      IF (DecimalPos+1>NINT2CHAR) THEN
         WRITE(*,*) 'Exceeded maximum limit handled by INT2CHAR'
         STOP
      END IF

      DO i=1,NINT2CHAR
         DecimalInFocus(i)=MOD(Natural,10**(i+1))/10**i
      END DO

      DO i=1,DecimalPos
         Ch(i)=CHAR(48+DecimalInFocus(DecimalPos+1-i))
      END DO
      Ch(DecimalPos+1)=CHAR(48+MOD(Natural,10))

      INT2CHAR=Ch(1)
      DO i=1,DecimalPos
         INT2CHAR=INT2CHAR(1:len_trim(INT2CHAR))//Ch(i+1)
      END DO
      RETURN
   END FUNCTION INT2CHAR
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CheckPotentialInformation(Potential,species1,species2)
      !checks whether the interaction potential for species1 and species 2 has been specified already
      IMPLICIT NONE
      TYPE(InteracPotential), POINTER :: Potential
      INTEGER :: species1,species2
      LOGICAL :: CheckPotentialInformation,Found
      
      Found=.FALSE.
      Found=Found .OR. ASSOCIATED(Potential%EAM)
      Found=Found .OR. ASSOCIATED(Potential%LJ)
      Found=Found .OR. ASSOCIATED(Potential%Coulomb)
      CheckPotentialInformation=Found
   END FUNCTION CheckPotentialInformation
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CentroSymmetricParameter(AL,iatom,nmax,rcut,csp,coordination,IsPrint)
   !finds the CSP for atom iatom present in AL with Nmax being given
   !based on Ju Li
   !see Kelchner, Plimpton and Hamilton, PRB 58, 17, 11085 (1998)
   !http://lammps.sandia.gov/doc/compute_centro_atom.html
   !nmax: number of nearest neighbors, e.g., fcc = 12, bcc = 8
   !rcut: nearest neighbor distance
      IMPLICIT NONE
      INTEGER, PARAMETER :: maxcoordination=100
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      INTEGER :: iatom,nmax,count,r1,r2,i,j,k,loc,nmaxuse,MaxAtomPerAtom,orig_listtype
      INTEGER :: indx(maxcoordination),coordination
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag,VLdr
      REAL(dp) :: csp,ri(3),rij(3),rj(3),rk(3),dr(3),dcut,Dsum,d2sum,Dj,drmag2
      REAL(dp), DIMENSION(maxcoordination) :: d,dcopy,dx,dy,dz
      REAL :: rcut
      LOGICAL, DIMENSION(maxcoordination) :: considered
      LOGICAL :: IsPrint1
      LOGICAL, OPTIONAL :: IsPrint
      
      csp=0.
      
      IF (PRESENT(IsPrint)) THEN
         IsPrint1=IsPrint
      ELSE
         IsPrint1=.FALSE.
      END IF
      
      !find nmaxuse
      nmaxuse=MAX((nmax/2),1)*2
      
      !sort neighbors of atom i according to their distances in ascending order
      VL=>AL%VL
      orig_listtype=VL%ListType
      IF (orig_listtype/=1) THEN
         WRITE(6,*) "Err>> CentroSymmetricParameter can work only for full verlet list"
         STOP
      END IF
      
      VLdrmag=>VL%drmag
      VLdr=>VL%dr
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      r1=(iatom-1)*MaxAtomPerAtom+1 !starting range of neighbors
      r2=r1-1+VL%ListRange(iatom) !end range of neighbors
      
      !find distances
      count=0
      DO j=r1,r2-1
         IF (VLdrmag(j)>rcut) CYCLE
         rj=VLdr(3*j-2:3*j)
         count=count+1
         d(count)=VLdrmag(j)
         dx(count)=rj(1)
         dy(count)=rj(2)
         dz(count)=rj(3)
         indx(count)=j
         IF (IsPrint1) WRITE(UNIT=6,FMT='(I7,":"I7,4ES15.3)') count,j,d(count),rj
      END DO
      coordination=count
      
      IF (IsPrint1) WRITE(6,*) "sort..."
      
      !sort and reset d, dx, dy and dz
      CALL sort_pick_abhijit(d(1:count),indx(1:count))
      DO i=1,count
         j=indx(i)
         d(i)=VLdrmag(j)
         rj=VLdr(3*j-2:3*j)
         dx(i)=rj(1)
         dy(i)=rj(2)
         dz(i)=rj(3)
         IF (IsPrint1) WRITE(UNIT=6,FMT='(I7,":",I7,4ES15.3)') i,j,d(i),rj
      END DO
      
      !find pairs and add to CSP
      d2sum=0._dp
      DO i=1,MIN(nmaxuse,count)
         d2sum=d2sum+d(i)*d(i)
      END DO
      IF (IsPrint1) WRITE(6,*) "Found d2sum:",d2sum
      
      considered=.FALSE.
      Dsum=0._dp
      DO i=1,MIN(nmaxuse,count)
         IF (considered(i)) CYCLE
         IF (IsPrint1) WRITE(6,*) "Atom i is ",i
         IF (IsPrint1) WRITE(6,*) "---------------"
         !find an atom which gives smallest |di+dj|
         Dj=10000._dp
         ri=(/dx(i),dy(i),dz(i)/) !vector for i
         k=0
         DO j=i+1,MIN(nmaxuse,count)
            IF (considered(j)) CYCLE
            rij=ri+(/dx(j),dy(j),dz(j)/) !vector for j
            drmag2=rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
            IF (Dj>drmag2) THEN
               Dj=drmag2
               k=j
               IF (IsPrint1) WRITE(6,*) "Atom j is ",j,rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
            END IF
         END DO
         IF (k==0) THEN
            Dj=0._dp
         ELSE
            considered(k)=.TRUE.
         END IF
         Dsum=Dsum+Dj
         IF (IsPrint1) WRITE(6,*) "Dsum for atom",i," is ",Dj
      END DO
      IF (IsPrint1) WRITE(*,*) "Total Dsum:",Dsum
      !Find CSP
      csp=0.5_dp*Dsum/d2sum
      IF (IsPrint1) WRITE(6,*) "CSP:",csp
      
   END SUBROUTINE CentroSymmetricParameter
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Welcome(IsPrint)
      IMPLICIT NONE
      INTEGER :: iunit,nentries,atmindx1,atmindx2,atmno1,atmno2,s1,s2,i,errorstatus
      INTEGER :: ispecies,nspecies
      CHARACTER(len=10) :: PotentialName
      CHARACTER :: CrapChar
      LOGICAL :: lexist,IsPrint1
      LOGICAL, OPTIONAL :: IsPrint
      
      IF (PRESENT(IsPrint)) THEN
         IsPrint1=IsPrint
      ELSE
         IsPrint1=.TRUE.
      END IF
      
      IF (IsPrint1) THEN
         WRITE(*,*) "********************************************************************"
         WRITE(*,*) "                 Multiscale simulation toolbox"
         WRITE(*,*) "                       Abhijit Chatterjee"
         WRITE(*,*) "********************************************************************"
      END IF

      !Read input.pot to get an idea of the number of species present
      
      iunit=301
      !INQUIRE(FILE="/usr/bin/license.amd.lic",EXIST=lexist)
      !IF (.NOT. lexist) THEN
         !WRITE(6,*) "Fatal error encountered."
         !STOP
      !END IF
      
      !INQUIRE(FILE="/usr/include/license.tad.lic",EXIST=lexist)
      !IF (.NOT. lexist) THEN
         !WRITE(6,*) "Fatal error encountered."
         !STOP
      !END IF
      
      OPEN(UNIT=iunit,FILE="input.pot",IOSTAT=errorstatus)
      IF (errorstatus/=0) THEN
         WRITE(*,*) "$Err>> Unable to open input.pot"
         STOP
      END IF
      
      READ(UNIT=iunit,FMT='(2I4)',IOSTAT=errorstatus) nentries,NSpecies_global
      write(*,*) nentries,NSpecies_global
      IF (errorstatus/=0) THEN
         WRITE(*,*) "$Err>> Unable to read # entries in input.pot"
         STOP
      END IF
      CALL MakeSize(SpeciesDirectory_global,NSpecies_global)
      SpeciesDirectory_global=0
      
      DO i=1,nentries
      
         READ(UNIT=iunit,FMT='(a10)',IOSTAT=errorstatus) PotentialName
         IF (errorstatus/=0) THEN
            WRITE(*,*) "$Err>> Unable to read PotentialName in input.pot at entry ",i
            STOP
         END IF
         SELECT CASE(PotentialName(1:2))
         CASE("LJ","BU") !LJ and Buckingham
            READ(iunit,*) s1,s2,atmno1,atmno2
            IF (XOR(s1==s2,atmno1==atmno2)) THEN
               WRITE(UNIT=6,FMT='("$Err>> One of species & atomic numbers matches,")',ADVANCE='NO')
               WRITE(UNIT=6,FMT='(" other does not, in entry number ",I2)') i
            END IF
            READ(iunit,*) CrapChar
            READ(iunit,*) CrapChar
            IF (SpeciesDirectory_global(s1)==0) THEN
               SpeciesDirectory_global(s1)=atmno1
               WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s1,SpeciesList%AtomicSymbol(atmno1)
            ELSE
               IF (SpeciesDirectory_global(s1)/=atmno1) THEN
                  WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                     "the atomic numbers already stored in db."
                  STOP
               END IF
            END IF
            IF (SpeciesDirectory_global(s2)==0) THEN
               SpeciesDirectory_global(s2)=atmno2
               WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s2,SpeciesList%AtomicSymbol(atmno2)
            ELSE
               IF (SpeciesDirectory_global(s2)/=atmno2) THEN
                  WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                  "the atomic numbers already stored in db."
                  STOP
               END IF
            END IF
         CASE("EA") !EAM
            IF (PotentialName(1:5)=="EAMPP" .OR. PotentialName(1:5)=="EAM  ") THEN
               READ(iunit,*) s1,s2,atmno1,atmno2
               IF (XOR(s1==s2,atmno1==atmno2)) THEN
                  WRITE(UNIT=6,FMT='("$Err>> One of species & atomic numbers matches,")',ADVANCE='NO')
                  WRITE(UNIT=6,FMT='(" other does not, in entry number ",I2)') i
               END IF
               READ(iunit,*) CrapChar
               READ(iunit,*) CrapChar
               IF (SpeciesDirectory_global(s1)==0) THEN
                  SpeciesDirectory_global(s1)=atmno1
                  WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s1,SpeciesList%AtomicSymbol(atmno1)
               ELSE
                  IF (SpeciesDirectory_global(s1)/=atmno1) THEN
                     WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                        "the atomic numbers already stored in db.",SpeciesDirectory_global(s1),atmno1
                     STOP
                  END IF
               END IF
               IF (SpeciesDirectory_global(s2)==0) THEN
                  SpeciesDirectory_global(s2)=atmno2
                  WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s2,SpeciesList%AtomicSymbol(atmno2)
               ELSE
                  IF (SpeciesDirectory_global(s2)/=atmno2) THEN
                     WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                        "the atomic numbers already stored in db.",SpeciesDirectory_global(s2),atmno2
                     STOP
                  END IF
               END IF
            ELSEIF (PotentialName(1:5)=="EAMDT" .OR. PotentialName(1:5)=="EAMEE") THEN
               READ(iunit,*) s1,atmno1
               READ(iunit,*) CrapChar
               READ(iunit,*) CrapChar
               IF (SpeciesDirectory_global(s1)==0) THEN
                  SpeciesDirectory_global(s1)=atmno1
                  WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s1,SpeciesList%AtomicSymbol(atmno1)
               ELSE
                  IF (SpeciesDirectory_global(s1)/=atmno1) THEN
                     WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                        "the atomic numbers already stored in db."
                     STOP
                  END IF
               END IF
            END IF
         CASE("CO") !COULOMB
            READ(iunit,*) nspecies
            DO ispecies=1,nspecies
               READ(iunit,*) s1,atmno1
               IF (SpeciesDirectory_global(s1)==0) THEN
                  SpeciesDirectory_global(s1)=atmno1
               ELSE
                  IF (SpeciesDirectory_global(s1)/=atmno1) THEN
                     WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                           "the atomic numbers already stored in db."
                     STOP
                  END IF
               END IF
            END DO
            IF (PotentialName(1:9)=="COULOMBEW") THEN
               READ(iunit,*) CrapChar
            END IF
            IF (PotentialName(1:9)=="COULOMBWO") THEN
               READ(iunit,*) CrapChar
               READ(iunit,*) CrapChar
            END IF
         CASE("SW")
            READ(iunit,*) s1,s2,atmno1,atmno2
            IF (XOR(s1==s2,atmno1==atmno2)) THEN
               WRITE(UNIT=6,FMT='("$Err>> One of species & atomic numbers matches,")',ADVANCE='NO')
               WRITE(UNIT=6,FMT='(" other does not, in entry number ",I2)') i
            END IF
            READ(iunit,*) CrapChar
            IF (SpeciesDirectory_global(s1)==0) THEN
               SpeciesDirectory_global(s1)=atmno1
               WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s1,SpeciesList%AtomicSymbol(atmno1)
            ELSE
               IF (SpeciesDirectory_global(s1)/=atmno1) THEN
                  WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                     "the atomic numbers already stored in db."
                  STOP
               END IF
            END IF
            IF (SpeciesDirectory_global(s2)==0) THEN
               SpeciesDirectory_global(s2)=atmno2
               WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s2,SpeciesList%AtomicSymbol(atmno2)
            ELSE
               IF (SpeciesDirectory_global(s2)/=atmno2) THEN
                  WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                  "the atomic numbers already stored in db."
                  STOP
               END IF
            END IF
         CASE("TE")
            READ(iunit,*) s1,s2,atmno1,atmno2
            IF (XOR(s1==s2,atmno1==atmno2)) THEN
               WRITE(UNIT=6,FMT='("$Err>> One of species & atomic numbers matches,")',ADVANCE='NO')
               WRITE(UNIT=6,FMT='(" other does not, in entry number ",I2)') i
            END IF
            READ(iunit,*) CrapChar
            IF (SpeciesDirectory_global(s1)==0) THEN
               SpeciesDirectory_global(s1)=atmno1
               WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s1,SpeciesList%AtomicSymbol(atmno1)
            ELSE
               IF (SpeciesDirectory_global(s1)/=atmno1) THEN
                  WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                     "the atomic numbers already stored in db."
                  STOP
               END IF
            END IF
            IF (SpeciesDirectory_global(s2)==0) THEN
               SpeciesDirectory_global(s2)=atmno2
               WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s2,SpeciesList%AtomicSymbol(atmno2)
            ELSE
               IF (SpeciesDirectory_global(s2)/=atmno2) THEN
                  WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                  "the atomic numbers already stored in db."
                  STOP
               END IF
            END IF
         
         CASE("FE")
         
         CASE("AI")
         
         CASE("ME")
         
         CASE("TB")

         CASE("DF") !density functional theory
            READ(iunit,*) s1,s2,atmno1,atmno2
            IF (s1/=s2 .OR. atmno1/=atmno2) THEN
               WRITE(UNIT=6,FMT='("$Err>> Species and/or atomic numbers mismatch,")',ADVANCE='NO')
               WRITE(UNIT=6,FMT='(" at entry ",I2)') i
               STOP
            END IF
            IF (SpeciesDirectory_global(s1)==0) THEN
               SpeciesDirectory_global(s1)=atmno1
               WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s1,SpeciesList%AtomicSymbol(atmno1)
            ELSE
               IF (SpeciesDirectory_global(s1)/=atmno1) THEN
                  WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                  "the atomic numbers already stored in db.",SpeciesDirectory_global(s1),atmno1
                  STOP
               END IF
            END IF
         CASE("OP") !opls potential
            READ(iunit,*) s1,s2,atmno1,atmno2
            IF (XOR(s1==s2,atmno1==atmno2)) THEN
            WRITE(UNIT=6,FMT='("$Err>> One of species & atomic numbers matches,")',ADVANCE='NO')
            WRITE(UNIT=6,FMT='(" other does not, in entry number ",I2)') i
            END IF
            READ(iunit,*) CrapChar
            IF (SpeciesDirectory_global(s1)==0) THEN
               SpeciesDirectory_global(s1)=atmno1
               WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s1,SpeciesList%AtomicSymbol(atmno1)
            ELSE
               IF (SpeciesDirectory_global(s1)/=atmno1) THEN
                  WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                  "the atomic numbers already stored in db."
                  STOP
               END IF
            END IF
            IF (SpeciesDirectory_global(s2)==0) THEN
               SpeciesDirectory_global(s2)=atmno2
               WRITE(UNIT=6,FMT='("Species ",I3," is ", a3)') s2,SpeciesList%AtomicSymbol(atmno2)
            ELSE
               IF (SpeciesDirectory_global(s2)/=atmno2) THEN
                  WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
                   "the atomic numbers already stored in db."
                  STOP
               END IF
            END IF
         CASE DEFAULT
            WRITE(6,*) "Unknown potential type has been specified in input.pot for entry number ",i
            WRITE(6,*) "Allowed potentials are LJ, EAM (EAMPP, EAMDT, EAMEE), COULOMBEW, "
            WRITE(6,*) "     BUCKINGHAM, SW (Stillinger-Weber), TE (Tersoff), FE (FENE), "
            WRITE(6,*) "     AI (AIREBO), ME (MEAM) and TB (Tight-binding)"
            STOP
         END SELECT
      !   READ(UNIT=iunit,FMT='(2I4)',ADVANCE='NO') atmindx1,atmindx2 !species index
      !   READ(UNIT=iunit,FMT='(2I4)') atmno1,atmno2 !atomic #s
      !   IF (XOR(atmindx1==atmindx2,atmno1==atmno2)) THEN
      !      WRITE(UNIT=6,FMT='("$Err>> One of species & atomic numbers matches,")',ADVANCE='NO')
      !      WRITE(UNIT=6,FMT='(" other does not, in line ",I2)') i+1
      !   END IF 
      !   IF (SpeciesDirectory_global(atmindx1)==0) THEN
      !      SpeciesDirectory_global(atmindx1)=atmno1
      !   ELSE
      !      IF (SpeciesDirectory_global(atmindx1)/=atmno1) THEN
      !         WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
      !           "the atomic numbers already stored in db."
      !         STOP
      !      END IF
      !   END IF
      !   IF (SpeciesDirectory_global(atmindx2)==0) THEN
      !      SpeciesDirectory_global(atmindx2)=atmno2
      !   ELSE
      !      IF (SpeciesDirectory_global(atmindx2)/=atmno2) THEN
      !         WRITE(6,*) "$Err>> Mismatch in atomic numbers provided and ", &
      !           "the atomic numbers already stored in db."
      !         STOP
      !      END IF
      !   END IF
         
      END DO
      CLOSE(301)
   END SUBROUTINE Welcome
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE utilities
