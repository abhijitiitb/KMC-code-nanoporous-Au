MODULE Crystal 
   !create structure based on crystal unit cell
   !unit cell provided below
   USE VARIABLE_TYPE
   USE db_manipulate
   USE utilities
   USE io
   USE ran_state
   USE Ecuyer_random
   IMPLICIT NONE
   REAL(dp), DIMENSION(9) :: Basis=0
   
   INTERFACE Create
      MODULE PROCEDURE CreateSurf,CreateDislocationEdge
   END INTERFACE
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE CreateFreeSurface()
   END SUBROUTINE CreateFreeSurface
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateSurf(AL,NCopies,CrystalType, &
     LatticeConst,AtomicNo,dimn,boxsclfac)
      !Create an AL with certain crystal structure (perfect crystal)
      !NCopies gives the # unit cells copied in each dimn
      !CrystalType is the unit cell type
      !   if CrystalType is "None" then AL is untouched in UnitCell
      !LatticeConst
      !   if LatticeConst is 0 then AL is untouched in FillCrystalDetails
      !AtomicNo for each species in AL%SpeciesDirectory
      !   where AL%SpeciesDirectory (3) implies first element
      !   type is Lithium
      !dimn gives dimensions where surface is exists
      !boxsclfac helps expand perfect crystal box to larger size
      
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER, INTENT(IN) :: NCopies(3)
      INTEGER, DIMENSION(:), INTENT(IN) :: AtomicNo
      CHARACTER(len=10), INTENT(IN) :: CrystalType
      REAL(dp), DIMENSION(:), INTENT(IN) :: LatticeConst
      INTEGER, DIMENSION(3), OPTIONAL :: dimn
      REAL(dp), DIMENSION(3), OPTIONAL :: boxsclfac
      INTEGER :: i,d
      CHARACTER(len=30) :: filename
      
      !========================Preliminary output======================
      WRITE(6,*) "Crystal>> Creating crystal structure ..."
      WRITE(6,*) "... ",TRIM(CrystalType)," crystal"
      
      DO i=1,SIZE(AtomicNo)
         WRITE(UNIT=6,FMT='(" ... species[",i2,"]:")',ADVANCE="NO") i
         WRITE(6,*) TRIM(SpeciesList%AtomicSymbol(AtomicNo(i)))
      END DO
      
      IF (SIZE(LatticeConst)==1) THEN
         WRITE(UNIT=*,FMT='(" ... lattice constant [",F7.3,",",F7.3,",",F7.3,"]")') &
            LatticeConst(1),LatticeConst(1),LatticeConst(1)
      ELSEIF (SIZE(LatticeConst)==3) THEN
         WRITE(UNIT=*,FMT='(" ... lattice constant [",F7.3,",",F7.3,",",F7.3,"]")') &
            LatticeConst(1),LatticeConst(2),LatticeConst(3)
      ELSE
         WRITE(6,*) "$Err>> LatticeConst should be either dimn 1 or 3"
         STOP
      END IF
      
      WRITE(UNIT=*,FMT='(" ... copies of unit cell ",I5,"  x",I5,"  x",I5)') &
         NCopies(1),NCopies(2),NCopies(3)
      
      !===========================Building the crystal=====================
      CALL UnitCell(AL,CrystalType) !get the unit cell
      WRITE(6,*) "... assigned unit cell"
      CALL FillCrystalDetails(AL,LatticeConst,AtomicNo)
      WRITE(6,*) "... filled unit cell details"
      CALL Expand(AL,NCopies) !make NCopies
      WRITE(6,*) "... expanded unit cell to specified number of copies"
      
      !============================create surface=========================
      IF (PRESENT(dimn) .AND. PRESENT(boxsclfac)) &
         AL%BoxSize(1:3)=boxsclfac*(AL%BoxSize(1:3)-AL%BoxSize(4:6))+ &
         AL%BoxSize(4:6)
   END SUBROUTINE CreateSurf
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateDislocationScrew()
   END SUBROUTINE CreateDislocationScrew
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateDislocationEdge()
   END SUBROUTINE CreateDislocationEdge
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateGBTwist(AL,NCopies,CrystalType,LatticeConst,AtomicNo,dimn,euler) 
      !grain boundary - 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1,AL2
      INTEGER, INTENT(IN) :: NCopies(3),AtomicNo(3),dimn
      CHARACTER(len=10), INTENT(IN) :: CrystalType
      REAL(dp), DIMENSION(:), INTENT(IN) :: LatticeConst
      REAL(dp), DIMENSION(:) :: euler
      
      
   END SUBROUTINE CreateGBTwist
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateGlass()
      IMPLICIT NONE
   END SUBROUTINE CreateGlass
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE InsertVacancy(n,imode)
      !imode=1 random,2 random clusters,3 anistropic random
      IMPLICIT NONE
      INTEGER :: n,imode
   END SUBROUTINE InsertVacancy
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE InsertInterstitial(n,imode)
      !imode=1 random,2 random clusters,3 anistropic random
      IMPLICIT NONE
      INTEGER :: n,imode
   END SUBROUTINE InsertInterstitial
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Scl(AL,sclfac1) 
      !scale lattice by factor
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(:) :: sclfac1
      REAL(dp), DIMENSION(3) :: sclfac
      INTEGER :: nsize,i
      
      nsize=SIZE(sclfac1)
      IF (nsize==1) THEN
         sclfac=PACK(SPREAD(sclfac1,2,3),.TRUE.)
      ELSEIF (nsize==3) THEN
         sclfac=sclfac1
      ELSE
         WRITE(6,*) "$Err>> Scaling factor for crystal has wierd # of dimensions"
         STOP
      END IF
      
      AL%BoxSize(1:3)=(AL%BoxSize(1:3)-AL%BoxSize(4:6))*sclfac+ &
           AL%BoxSize(4:6)
      DO i=1,AL%NAtoms
         AL%AtomCoord(3*i-2:3*i)=(AL%AtomCoord(3*i-2:3*i)- &
              AL%BoxSize(4:6))*sclfac+AL%BoxSize(4:6)
      END DO
   END SUBROUTINE Scl
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE FillCrystalDetails(AL,LatticeConst,atomicno)
   !When AL contains only one unit cell - gives the correct length dimensions to the systems
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(:), INTENT(IN) :: LatticeConst
      REAL(dp) :: refcoord(3),abox(3)
      REAL(dp), DIMENSION(:), POINTER :: LatticeConstexp,refcoordexp
      INTEGER, DIMENSION(:) :: atomicno
      INTEGER :: NAtoms,atmno,i

      IF (ALL(LatticeConst==0.0_dp)) THEN
         WRITE(6,*) "Err>> Lattice constant should be non-zero"
         RETURN
      END IF
      
      IF (NSpecies_global<SIZE(atomicno)) THEN
      !IF (AL%NSpeciesType/=SIZE(atomicno)) THEN
         WRITE(6,*) "$Warning: species atomicnumber array of incorrect size"
      END IF
      
      !Get correct atomic # and mass for AL
      !AL%SpeciesDirectory=atomicno !order should be correct
      CALL GetAtomListMass(AL)
      
      !Setup correct unit cell size
      refcoord=AL%BoxSize(4:6) !gives min bounding box for UnitCell
      NAtoms=AL%NAtoms
      CALL MakeSize(refcoordexp,3*NAtoms)
      CALL MakeSize(LatticeConstexp,3*NAtoms)
      
      refcoordexp=PACK(SPREAD(refcoord,2,NAtoms),.TRUE.)
      IF (SIZE(LatticeConst)==1) THEN !all sides of UnitCell are of equal size
         !Check if LatticeConst is really equal in all dimn
         !IF (.NOT. (AL%BoxSize(1)==AL%BoxSize(2) .AND. &
         !   AL%BoxSize(1)==AL%BoxSize(3)) ) THEN
         !   WRITE(6,*) "$Err>> Original unit cell has sides of different length"
         !   STOP
         !END IF
         !Spread LatticeConst over a larger array
         !Note LatticeConst is defn such than when
         !multiplied to BoxSize of UnitCell the correct
         !size is obtained
         LatticeConstexp=PACK(SPREAD(LatticeConst,2,3*NAtoms),.TRUE.)
         abox=PACK(SPREAD(LatticeConst,2,3),.TRUE.)
      ELSEIF (SIZE(LatticeConst)==3) THEN !different sizes
         LatticeConstexp=PACK(SPREAD(LatticeConst,2,NAtoms),.TRUE.)
         abox=LatticeConst
      ELSE
         WRITE(6,*) "$Err: Unknown format of lattice scaling encountered"
         STOP
      END IF
      
      AL%AtomCoord=(AL%AtomCoord-refcoordexp)*LatticeConstexp+refcoordexp
      AL%BoxSize(1:3)=(AL%BoxSize(1:3)-refcoord)*abox+refcoord
      Basis(1:3)=Basis(1:3)*abox
      Basis(4:6)=Basis(4:6)*abox
      Basis(7:9)=Basis(7:9)*abox
      DEALLOCATE(refcoordexp,LatticeConstexp)
      
   END SUBROUTINE FillCrystalDetails
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Cleave(AL,plane,coord,mask,UseOldMask)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      LOGICAL, DIMENSION(:), POINTER :: mask
      REAL(dp), DIMENSION(3,3) :: plane
      REAL(dp), DIMENSION(3) :: coord
      LOGICAL, OPTIONAL :: UseOldMask
      LOGICAL :: UseOldMask1
      INTEGER :: i,count,dimn

      IF (PRESENT(UseOldMask)) THEN
         UseOldMask1=UseOldMask
      ELSE
         UseOldMask1=.FALSE.
      END IF
      
      count=0
      DO i=1,3
         IF (plane(1,i)==plane(2,i) .AND. plane(2,i)==plane(3,i)) THEN
            count=count+1
            dimn=i
         END IF
      END DO
      IF (count>1) THEN
         WRITE(6,*) "$Err>> Straight line provided instead of plane"
         STOP
      ELSEIF (count==0) THEN
         CALL Cleave3(AL,plane,coord,mask)
      ELSE
         CALL Cleave1(AL,plane,coord,mask,dimn,UseOldMask1)
      END IF
   END SUBROUTINE Cleave
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Cleave1(AL,plane,coord,mask,dimn,UseOldMask)
      !coord is a location which is not deleted
      !cut along x or y or z planes
      !mask is of size NAtoms 
      !Role of mask
      !  mask can be used to perform complicated cuts to the crystal
      !  e.g., carving a slice out of the crystal
      !  only atoms with mask true are checked for selection
      !  NOTE: in some cases it makes sense to save changes before further 
      !     cuts are made
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      LOGICAL, DIMENSION(:), POINTER :: mask
      REAL(dp), DIMENSION(3,3) :: plane
      REAL(dp), DIMENSION(3) :: coord
      INTEGER :: dimn,r1,MarkedAtoms,i,NAtoms
      LOGICAL, OPTIONAL :: UseOldMask

      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(6,*) "$Err>> Atom list provide to cleave is not associated"
         STOP
      END IF
      NAtoms=AL%NAtoms
      IF (NAtoms==0) THEN
         WRITE(6,*) "$Err>> No atoms present in Atom list provided to cleave"
         STOP
      END IF
      WRITE(UNIT=*,FMT='(" ...cleaving crystal along dimension ",I1)') dimn
      
      !if mask already provided then work with it
      IF (PRESENT(UseOldMask)) THEN
         IF (.NOT. UseOldMask) THEN
            CALL MakeSize(mask,NAtoms)
            mask=.TRUE.
         END IF
      ELSE !clear old mask
         IF (.NOT. ASSOCIATED(mask)) THEN
            CALL MakeSize(mask,NAtoms)
            mask=.TRUE.
         END IF
      END IF
      
      MarkedAtoms=0
      r1=-3+dimn
      IF (coord(dimn)-plane(1,dimn)==0._dp) THEN
         WRITE(6,*) "$Err>> side of plane to be selected is unclear"
         STOP
      END IF

      DO i=1,NAtoms
         IF (.NOT. mask(i)) CYCLE !ignore
         r1=r1+3
         !selected atoms should lie on same side as coord
         mask(i)=(AL%AtomCoord(r1)-plane(1,dimn))/(coord(dimn)-plane(1,dimn))>0._dp
         IF (mask(i)) MarkedAtoms=MarkedAtoms+1
      END DO
      
      WRITE(UNIT=*,FMT='(" ...marked ",I5," atoms out of ",I5," atoms")') &
         MarkedAtoms,NAtoms

   END SUBROUTINE Cleave1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Cleave3(AL,plane,coord,mask) !cut using a plane
      !coord is a location which is not deleted
      !mask is of size NAtoms
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      LOGICAL, DIMENSION(:), POINTER :: mask
      REAL(dp), INTENT(IN) :: plane(3,3),coord(3)
      REAL(dp) :: A(3),D,coord1(3),s,s1,mat(3,3)
      INTEGER :: i,r1,r2,NAtoms,MarkedAtoms
      
      !Eq. of plane - background
      !plane(1,:)=(x1,y1,z1) - point 1
      !plane(2,:)=(x2,y2,z2) - point 2
      !plane(3,:)=(x3,y3,z3) - point 3
      !these 3 points determine the plane
      
      !Eq. of plane Ax+By+Cz=D where normal to the plane is (A,B,C)
      !A = y1 (z2 - z3) + y2 (z3 - z1) + y3 (z1 - z2) 
      !B = z1 (x2 - x3) + z2 (x3 - x1) + z3 (x1 - x2) 
      !C = x1 (y2 - y3) + x2 (y3 - y1) + x3 (y1 - y2) 
      !D = x1 (y2 z3 - y3 z2) + x2 (y3 z1 - y1 z3) + x3 (y1 z2 - y2 z1)
      
      !The sign of s = Ax + By + Cz - D determines 
      !which side the point (x,y,z) lies with respect 
      !to the plane. If s > 0 then the point lies on 
      !the same side as the normal (A,B,C). If s < 0 
      !then it lies on the opposite side, if s = 0 then 
      !the point (x,y,z) lies on the plane. 
      
      WRITE(6,*) " Cleave3 has not been tested"
      STOP
      
      !Our A,B,C are represented by A(3)
      mat=plane; mat(:,1)=1; A(1)=DET(mat)
      mat=plane; mat(:,2)=1; A(2)=DET(mat)
      mat=plane; mat(:,3)=1; A(3)=DET(mat)
      mat=plane; D=DET(mat)
      IF (D==0._dp) THEN
         WRITE(6,*) "$Err>> Choose positions in plane differently"
         WRITE(6,*) "  so that determinant in cleave is not zero"
         STOP
      END IF
      s=DOT_PRODUCT(A,(/0._dp,1._dp,51.8_dp/))-D
      s=DOT_PRODUCT(A,coord)-D
      s=s/ABS(s)
      
      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(6,*) "$Err>> Atom list provide to cleave is not associated"
         STOP
      END IF
      NAtoms=AL%NAtoms
      IF (NAtoms==0) THEN
         WRITE(6,*) "$Err>> No atoms present in Atom list provided to cleave"
         STOP
      END IF
      WRITE(6,*) "...cleaving crystal"
      
      !if mask already provided then work with it
      IF (.NOT. ASSOCIATED(mask)) THEN
         CALL MakeSize(mask,NAtoms)
         mask=.TRUE.
      END IF
      
      MarkedAtoms=0
      r1=-2; r2=0
      DO i=1,NAtoms
         IF (.NOT. mask(i)) CYCLE !ignore
         r1=r1+3; r2=r2+3
         coord1=AL%AtomCoord(r1:r2)
         s1=DOT_PRODUCT(A,coord1)-D
         s1=s1/ABS(s1)
         mask(i)=(s1==s)
         IF (mask(i)) MarkedAtoms=MarkedAtoms+1
      END DO
      
      WRITE(UNIT=*,FMT='(" ... marked ",I5," atoms out of ",I5," atoms")') &
         MarkedAtoms,NAtoms
      !make sure mask is deleted at the end
   END SUBROUTINE Cleave3
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE CarveCircle(AL,centercoord,radius,mask,normaldirection)
   !selects all atoms outside a circle of radius "radius" with center "centercoord" that need to be deleted
   !normaldirection is the direction in which the normal to the circle plane is pointing - it can be 1, 2 or 3. For eg., normaldirection=3 means CarveCirle in the xy plane
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp) :: radius,centercoord(3),coord(3),dcoord(3),rad,rad2,radius2
      LOGICAL, DIMENSION(:), POINTER :: mask
      INTEGER :: i,NAtoms,normaldirection,r1,r2
      
      NAtoms=AL%NAtoms
      IF (NAtoms==0) THEN
         WRITE(6,*) "$Err>> No atoms present in Atom list"
         STOP
      END IF
      WRITE(6,*) "...selecting atoms inside circle"
      
      !if mask already provided then work with it
      IF (.NOT. ASSOCIATED(mask)) THEN
         CALL MakeSize(mask,NAtoms)
         mask=.TRUE.
      END IF
      
      r1=-2; r2=0
      radius2=radius*radius
      DO i=1,NAtoms
         r1=r1+3
         r2=r2+3
         coord=AL%AtomCoord(r1:r2)
         dcoord=coord-centercoord
         dcoord(normaldirection)=0._dp
         rad2=dcoord(1)*dcoord(1)+dcoord(2)*dcoord(2)+dcoord(3)*dcoord(3)
         mask(i)= rad2>radius2
      END DO
      
   END SUBROUTINE CarveCircle
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE CarveSphere(AL,centercoord,radius,mask)
   !selects all atoms outside a sphere of radius "radius" with center "centercoord" that need to be deleted
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp) :: radius,centercoord(3),coord(3),dcoord(3),rad,rad2,radius2
      LOGICAL, DIMENSION(:), POINTER :: mask
      INTEGER :: i,NAtoms,r1,r2
      
      NAtoms=AL%NAtoms
      IF (NAtoms==0) THEN
         WRITE(6,*) "$Err>> No atoms present in Atom list"
         STOP
      END IF
      WRITE(6,*) "...selecting atoms inside circle"
      
      !if mask already provided then work with it
      IF (.NOT. ASSOCIATED(mask)) THEN
         CALL MakeSize(mask,NAtoms)
         mask=.TRUE.
      END IF
      
      r1=-2; r2=0
      radius2=radius*radius
      DO i=1,NAtoms
         r1=r1+3
         r2=r2+3
         coord=AL%AtomCoord(r1:r2)
         dcoord=coord-centercoord
         rad2=dcoord(1)*dcoord(1)+dcoord(2)*dcoord(2)+dcoord(3)*dcoord(3)
         mask(i)= rad2>radius2
      END DO
      
   END SUBROUTINE CarveSphere
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE SelectLowCoordinationAtoms(AL,coordinationdistance,mincoordination,probability,mask)
   !selects the atoms that have a coordination number lower than mincoordination
   !if probability is present then atoms are removed according to the probability
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), OPTIONAL :: probability
      REAL(dp) :: mincoordination,probability1,coordinationdistance,rnd
      REAL(dp), DIMENSION(:), POINTER :: Coordination
      LOGICAL, DIMENSION(:), POINTER :: mask
      INTEGER :: NAtoms,i
      
      IF (PRESENT(probability)) THEN
         probability1=probability
      ELSE
         probability1=1.
      END IF
      
      NAtoms=AL%NAtoms
      !if mask already provided then work with it
      IF (.NOT. ASSOCIATED(mask)) THEN
         CALL MakeSize(mask,NAtoms)
         mask=.FALSE.
      END IF
      
      CALL FindCoordination(AL,coordinationdistance,Coordination)
      
      DO i=1,NAtoms
         IF (Coordination(i)<mincoordination) THEN
            rnd=taus88()
            mask(i)= rnd<probability1
         END IF
      END DO
      
   END SUBROUTINE SelectLowCoordinationAtoms
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE RandomizeAtoms(AL,displacement)
   !randomly displace atoms
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord
      REAL(dp) :: displacement
      INTEGER :: icoord
      
      AtomCoord=>AL%AtomCoord
      DO icoord=1,3*AL%NAtoms
         AtomCoord(icoord)=AtomCoord(icoord)+(1._dp-2.*taus88())*displacement
      END DO
      CALL PBC(AL)
   END SUBROUTINE RandomizeAtoms
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE FindCoordination(AL,coordinationdistance,Coordination)
   
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp) :: coordinationdistance,coordinationdistance2,coord1(3),coord2(3),dcoord(3),distance2
      REAL(dp), DIMENSION(:), POINTER :: Coordination
      INTEGER :: i,j,NAtoms
      
      !if Coordination already provided then work with it
      NAtoms=AL%NAtoms
      IF (.NOT. ASSOCIATED(Coordination)) CALL MakeSize(Coordination,NAtoms)
      
      coordinationdistance2=coordinationdistance*coordinationdistance
      DO i=1,NAtoms-1
         coord1=AL%AtomCoord(3*i-2:3*i)
         DO j=i+1,NAtoms
            coord2=AL%AtomCoord(3*j-2:3*j)
            dcoord=coord1-coord2
            distance2=dcoord(1)*dcoord(1)+dcoord(2)*dcoord(2)+dcoord(3)*dcoord(3)
            IF (distance2<coordinationdistance2) THEN
               Coordination(i)=Coordination(i)+1
               Coordination(j)=Coordination(j)+1
            END IF
         END DO
      END DO
      
   END SUBROUTINE FindCoordination
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE AlignCOMWithBoxCenter(AL)
      !aligns center of mass with box center
      !can be used with PBC to get better crystal shapes
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(3) :: COM,BoxCenter
      
      COM=GetCenterOfMass(AL)
      BoxCenter=(AL%BoxSize(4:6)+AL%BoxSize(1:3))/2._dp
      CALL Translate(AL,BoxCenter-COM,translatebox=.FALSE.)
      
   END SUBROUTINE AlignCOMWithBoxCenter
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Translate(AL,t,mask,translatebox)
      !mask specifies which coordinates are allowed to be translated
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp) :: t(3)
      REAL(dp), DIMENSION(:), POINTER :: texpand
      LOGICAL, DIMENSION(:), POINTER, OPTIONAL :: mask
      LOGICAL, OPTIONAL :: translatebox
      LOGICAL, DIMENSION(:), POINTER :: mask1
      LOGICAL :: translatebox1
      INTEGER :: i
      
      CALL MakeSize(mask1,3*AL%NAtoms)
      IF (PRESENT(mask)) THEN
         mask1=mask
      ELSE
         mask1=.TRUE.
      END IF
      
      IF (PRESENT(translatebox)) THEN
         translatebox1=translatebox
      ELSE
         translatebox1=.TRUE.
      END IF
      
      CALL MakeSize(texpand,3*AL%NAtoms)
      texpand=PACK(SPREAD(t,2,AL%NAtoms),.TRUE.)
      WHERE (mask1) AL%AtomCoord=AL%AtomCoord+texpand
      DEALLOCATE(mask1,texpand)
      
      !move upper and lower bounds
      !IF (translatebox1) AL%BoxSize=AL%BoxSize+PACK(SPREAD(t,2,2),.TRUE.) 
   END SUBROUTINE Translate
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RotateSpherical(AL,angle,origin,mask,units)
      !Rotate crystal in spherical angles
      !angle is size 2 in degrees by default
      !angle(1) is psi coordinate
      !angle(2) is theta coordinate
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(2) :: angle,angle1
      REAL(dp), DIMENSION(3) :: origin
      REAL(dp), DIMENSION(:), POINTER :: coord
      REAL(dp) :: radius,radiusxy,theta,phi
      LOGICAL, DIMENSION(:), POINTER, OPTIONAL :: mask
      LOGICAL, DIMENSION(:), POINTER :: mask1
      CHARACTER(len=3), OPTIONAL :: units
      CHARACTER(len=3) :: units1
      INTEGER :: i,NAtoms,r1,r2
      
      NAtoms=AL%NAtoms
      CALL MakeSize(mask1,NAtoms) !mask1 should be of size NAtoms
      IF (PRESENT(mask)) THEN
         mask1=mask
      ELSE
         mask1=.TRUE.
      END IF
      IF (PRESENT(units)) THEN
         units1=units
      ELSE
         units1="deg"
      END IF
      WRITE(6,*) "... about to rotate crystal in spherical coord"
      WRITE(UNIT=*,FMT='("   ... about origin [",3f8.3,"]")') origin
      WRITE(UNIT=*,FMT='("   ... rotation [",2f8.3,"]")',ADVANCE="NO") angle
      SELECT CASE (units1(1:3))
      CASE("rad")
         angle1=angle
      CASE("deg")
         angle1=angle*PI/180._dp !convert to radians
      CASE DEFAULT
         WRITE(6,*) "$Err>> Unknown spherical angle units specified"
         STOP
      END SELECT
      WRITE(6,*) TRIM(units1)
      
      r1=-2; r2=0
      CALL MakeSize(coord,3)
      DO i=1,NAtoms
         r1=r1+3; r2=r2+3
         coord=AL%AtomCoord(r1:r2)-origin !before rotation
         radius=NORM(coord,3)
         
         IF (radius<1.e-10_dp) THEN
            WRITE(6,*) "$Err>> Radius is too small, select different origin"
            STOP
         END IF
         !get z angle theta
         theta=ACOS(coord(3)/radius)
         !get xy angle phi
         radiusxy=ABS(radius*SIN(theta))
         
         IF (radiusxy==0._dp .AND. angle(2)/=0._dp) THEN
            WRITE(6,*) "$Err>> Undefined rotation of pole requested"
            STOP
         END IF
         
         IF (radiusxy<1.e-10_dp) THEN
            WRITE(6,*) "$Err>> Radius is too small, select different origin"
            STOP
         END IF
         
         phi=ACOS(ABS(coord(1))/radiusxy) !this is in the first quadrant
         IF (coord(1)>0._dp) THEN
            IF (coord(2)>0._dp) THEN !1st quadrant
            ELSE !4th quadrant
               phi=2*PI-phi
            END IF
         ELSE
            IF (coord(2)>0._dp) THEN !2nd quadrant
               phi=PI-phi
            ELSE !3rd quadrant
               phi=PI+phi
            END IF
         END IF
         !now check if phi is correct
         IF (ABS(radiusxy*COS(phi)-coord(1))>0.0001_dp .OR.  &
            ABS(radiusxy*SIN(phi)-coord(2))>0.0001_dp) THEN
            WRITE(6,*) "$Err>> Error in calculating phi"
            STOP
         END IF
         
         !add the angle to rotate
         phi=phi+angle1(1) !phi & theta can be rotated as much as you want
         theta=theta-angle1(2)
         radiusxy=radius*SIN(theta)
         coord(1)=radiusxy*COS(phi)
         coord(2)=radiusxy*SIN(phi)
         coord(3)=radius*COS(theta)
         AL%AtomCoord(r1:r2)=coord+origin
      END DO
      DEALLOCATE(coord)
      
      WRITE(6,*) "...BoxSize of Atom list will be incorrect due to rotation"
      WRITE(6,*) "...rotation completed"
   END SUBROUTINE RotateSpherical
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RotateEuler(AL,euler,origin,mask,units)
      !rotate about origin with angle euler
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(3) :: euler,euler1,origin
      REAL(dp), DIMENSION(:), POINTER :: originexpand
      REAL(dp) :: rotmatrix(3,3)
      LOGICAL, DIMENSION(:), POINTER, OPTIONAL :: mask
      LOGICAL, DIMENSION(:), POINTER :: mask1
      CHARACTER(len=3), OPTIONAL :: units
      CHARACTER(len=3) :: units1
      INTEGER :: i,NAtoms,r1,r2
      
      WRITE(6,*) "... about to rotate crystal "
      WRITE(UNIT=*,FMT='("   ... about origin [",3f8.3,"]")') origin
      NAtoms=AL%NAtoms
      rotmatrix=RotationMatrix(euler)
      CALL MakeSize(originexpand,3*NAtoms)
      CALL MakeSize(mask1,3*NAtoms)
      IF (PRESENT(mask)) THEN
         mask1=mask
      ELSE
         mask1=.TRUE.
      END IF
      IF (PRESENT(units)) THEN
         units1=units
      ELSE
         units1="deg"
      END IF
      
      SELECT CASE (units(1:3))
      CASE("rad")
      CASE("deg")
      CASE DEFAULT
      END SELECT
      originexpand=PACK(SPREAD(origin,2,AL%NAtoms),.TRUE.)
      WHERE(mask1) AL%AtomCoord=AL%AtomCoord-originexpand
      r1=-2; r2=0
      DO i=1,NAtoms
         r1=r1+3; r2=r2+3
         AL%AtomCoord(r1:r2)=MATMUL(rotmatrix,AL%AtomCoord(r1:r2))
      END DO
      WHERE(mask1) AL%AtomCoord=AL%AtomCoord+originexpand
      
      DEALLOCATE(originexpand,mask1)
   END SUBROUTINE RotateEuler
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ChangeSpeciesAtomicNumber(AL,atomicnumber,mask)
      !mask corresponds to atoms that are marked (see Mark)
      !for marked atoms change their species type
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      LOGICAL, DIMENSION(:), POINTER, OPTIONAL :: mask
      LOGICAL, DIMENSION(:), POINTER :: mask1
      INTEGER, INTENT(IN) :: atomicnumber
      INTEGER :: speciesindx
      
      CALL MakeSize(mask1,AL%NAtoms)
      IF (PRESENT(mask)) THEN
         mask1=mask
      ELSE
         mask1=.TRUE.
      END IF
      
      !speciesindx=GetSpecies(AL,atomicnumber)
      speciesindx=GetSpecies(atomicnumber)
      WHERE(mask1) AL%AtomSpecies=speciesindx
      DEALLOCATE(mask1)
   END SUBROUTINE ChangeSpeciesAtomicNumber
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ChangeSpeciesIndex(AL,newspeciesindx,mask)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: newspeciesindx,NAtoms
      LOGICAL, DIMENSION(:), POINTER, OPTIONAL :: mask !mask is of size AL%NAtoms
      LOGICAL, DIMENSION(:), POINTER :: mask1 !mask is of size AL%NAtoms
         
      NAtoms=AL%NAtoms
      CALL MakeSize(mask1,3*AL%NAtoms)
      IF (PRESENT(mask)) THEN
         mask1=mask
      ELSE
         mask1=.TRUE.
      END IF
      
      WHERE (mask1) AL%AtomSpecies=newspeciesindx
      DEALLOCATE(mask1)
      WRITE(6,*) "Change atom species to ",newspeciesindx
   END SUBROUTINE ChangeSpeciesIndex
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE FreezeAtoms(AL,mask,iprint)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      LOGICAL, DIMENSION(:), INTENT(IN) :: mask
      INTEGER, OPTIONAL :: iprint
      INTEGER :: i,iprint1
      
      IF (PRESENT(iprint)) THEN
         iprint1=iprint
      ELSE
         iprint1=0
      END IF
      
      DO i=1,AL%NAtoms
         IF (mask(i)) THEN
            AL%AtomVelocity(3*i-2:3*i)=0._dp
            AL%AtomIsMoving(3*i-2:3*i)=.FALSE.
            AL%NMove=AL%NMove-1
         END IF
         IF (iprint1>0 .AND. mask(i)) WRITE(UNIT=6,FMT='("   ..Frozen atom: ",I4)') i
      END DO
   END SUBROUTINE FreezeAtoms
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReplaceAtomSpeciesRandomly(AL,NAtoms,NReplace,species1,species2)
   !Replaces NReplace atoms of type species1 with atoms of type species2
   !AtomSpecies: array holding the species types of the atoms
   !NAtoms: total number of atoms present
   !NReplace: number of atoms to be replaced
      IMPLICIT NONE
      INTEGER :: NAtoms,NReplace,species1,species2,NAtomSpecies,iatom
      TYPE(SystemContainer), POINTER :: AL
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      REAL(dp), DIMENSION(:), POINTER :: AtomMass
      LOGICAL :: NotSatisfied
      REAL(dp) :: r
      
      AtomSpecies=>AL%AtomSpecies
      AtomMass=>AL%AtomMass
      
      NAtomSpecies=0 !number of atoms of type RemoveSpecies
      DO iatom=1,NAtoms
         IF (AtomSpecies(iatom)==species1) NAtomSpecies=NAtomSpecies+1
      END DO
      IF (NAtomSpecies<NReplace) THEN
         WRITE(6,*) "Err>> Number of atoms to be replaced is more than number of atoms present"
         STOP
      END IF
      
      NotSatisfied=.TRUE.
      NAtomSpecies=0
      DO WHILE (NotSatisfied)
         r=taus88()
         iatom=CEILING(r*REAL(NAtoms))
         IF (AtomSpecies(iatom)==species1 .AND. NotSatisfied) THEN
            r=taus88()
            IF (r>0.5_dp) THEN
               NAtomSpecies=NAtomSpecies+1
               !WRITE(UNIT=6,FMT='(I3,". Atom ",i6," was species ",i2)',ADVANCE="NO") &
               !   NAtomSpecies,iatom,AtomSpecies(iatom)
               AtomSpecies(iatom)=species2 !this atom should be removed
               !WRITE(UNIT=6,FMT='(" and now is species ",i2)') AtomSpecies(iatom)
            END IF
         END IF
         NotSatisfied=NAtomSpecies<NReplace
      END DO
      
      CALL GetAtomListMass(AL)
      
   END SUBROUTINE ReplaceAtomSpeciesRandomly
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReplaceAtomSpeciesCompositionProfile(AL,NAtoms,dimn,species1,species2)
   !Replaces NReplace atoms of type species1 with atoms of type species2
   !AtomSpecies: array holding the species types of the atoms
   !NAtoms: total number of atoms present
   !NReplace: number of atoms to be replaced
   !dimn: x=1, y=2, z=3
      IMPLICIT NONE
      INTEGER :: NAtoms,species1,species2,NAtomSpecies,iatom
      TYPE(SystemContainer), POINTER :: AL
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      REAL(dp), DIMENSION(:), POINTER :: AtomMass,AtomCoord
      INTEGER :: dimn
      REAL(dp) :: r,m,position
      
      AtomCoord=>AL%AtomCoord
      AtomSpecies=>AL%AtomSpecies
      
      NAtomSpecies=0
      
      DO iatom=1,NAtoms
         position=AtomCoord(3*(iatom-1)+dimn)
         
         IF (AtomSpecies(iatom)==species1) THEN
            r=CompositionProfile(position)
            m=taus88()
            IF (m>r) THEN !we need to change the species type
               NAtomSpecies=NAtomSpecies+1
               AtomSpecies(iatom)=species2 !this atom should be removed
               !WRITE(UNIT=6,FMT='(I3,". Atom ",i6," was species ",i2)',ADVANCE="NO") &
               !   NAtomSpecies,iatom,AtomSpecies(iatom)
            END IF
         END IF
      END DO
      
      CALL GetAtomListMass(AL)
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      FUNCTION CompositionProfile(pos)
      !to be created by the user -- this will be specific to the system in hand
         IMPLICIT NONE
         REAL(dp) :: pos
         REAL(dp) :: CompositionProfile
         REAL, PARAMETER:: R=29.865_dp, D=13.575_dp
         
         IF (pos<=16.29_dp) THEN
            CompositionProfile=1.0_dp ! pure species1 is here
         ELSEIF (pos>=43.44_dp) THEN
            CompositionProfile=0.0_dp  ! pure Germanium is here
         ELSE
            CompositionProfile=0.5_dp-0.5_dp*SIN(1.57_dp*((pos-R)/D))
         END IF
      END FUNCTION CompositionProfile
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      
   END SUBROUTINE ReplaceAtomSpeciesCompositionProfile
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteAtomsRandomly(AtomSpecies,NAtoms,NRemove,RemoveSpecies,mask)
   !removes NRemove atoms from the list
   !NAtoms: total number of atoms present
   !NRemove: number of atoms to be remove of speciestype RemoveSpecies
   !RemoveSpecies: is species type for e.g., 1, 2 ...
   !mask is true if species is to be retained
      IMPLICIT NONE
      INTEGER :: NAtoms,NRemove,RemoveSpecies,iatom,NAtomSpecies
      INTEGER, DIMENSION(:) :: AtomSpecies
      LOGICAL, DIMENSION(:) :: mask
      LOGICAL :: NotSatisfied
      REAL(dp) :: r
      
      NAtomSpecies=0 !number of atoms of type RemoveSpecies
      DO iatom=1,NAtoms
         IF (AtomSpecies(iatom)==RemoveSpecies) NAtomSpecies=NAtomSpecies+1
      END DO
      IF (NAtomSpecies<NRemove) THEN
         WRITE(6,*) "Err>> Number of atoms to be removed is more than number of atoms present"
         STOP
      END IF
      
      NotSatisfied=.TRUE.
      NAtomSpecies=0
      mask(1:NAtoms)=.TRUE.
      DO WHILE (NotSatisfied)
         r=taus88()
         iatom=CEILING(r*REAL(NAtoms))
         IF (AtomSpecies(iatom)==RemoveSpecies .AND. mask(iatom) .AND. NotSatisfied) THEN
            r=taus88()
            IF (r>0.5_dp) THEN
               NAtomSpecies=NAtomSpecies+1
               mask(iatom)=.FALSE. !this atom should be removed
               WRITE(UNIT=6,FMT='(I3,". Atom ",i6," of species ",i2," is removed")') &
                  NAtomSpecies,iatom,AtomSpecies(iatom)
            END IF
         END IF
         NotSatisfied=NAtomSpecies<NRemove
      END DO
   END SUBROUTINE DeleteAtomsRandomly
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE InsertAtomsRandomly(AtomSpecies,NAtoms,NRemove,RemoveSpecies,mask)
   !removes NRemove atoms from the list
   !NAtoms: total number of atoms present
   !NRemove: number of atoms to be remove of speciestype RemoveSpecies
   !RemoveSpecies: is species type for e.g., 1, 2 ...
   !mask is true if species is to be retained
      IMPLICIT NONE
      INTEGER :: NAtoms,NRemove,RemoveSpecies,iatom,NAtomSpecies
      INTEGER, DIMENSION(:) :: AtomSpecies
      LOGICAL, DIMENSION(:) :: mask
      LOGICAL :: NotSatisfied
      REAL(dp) :: r
      CHARACTER(len=100) :: filename1
      
!      NAtomSpecies=0 !number of atoms of type RemoveSpecies
!      DO iatom=1,NAtoms
!         IF (AtomSpecies(iatom)==RemoveSpecies) NAtomSpecies=NAtomSpecies+1
!      END DO
!      IF (NAtomSpecies<NRemove) THEN
!         WRITE(6,*) "Err>> Number of atoms to be removed is more than number of atoms present"
!         STOP
!      END IF
      
!      NotSatisfied=.TRUE.
!      NAtomSpecies=0
!      mask(1:NAtoms)=.TRUE.
      
!     filename1="Ag.xyz" 
!     OPEN(UNIT=101,FILE=filename1) 
      
!      DO WHILE (NotSatisfied)
!         r=taus88()
!         iatom=CEILING(r*REAL(NAtoms))
!         IF (AtomSpecies(iatom)==RemoveSpecies .AND. mask(iatom) .AND. NotSatisfied) THEN
!            r=taus88()
!            IF (r>0.5_dp) THEN
!               NAtomSpecies=NAtomSpecies+1
!               mask(iatom)=.FALSE. !this atom should be removed
!               WRITE(UNIT=6,FMT='(I3,". Atom ",i6," of species ",i2," is removed")') &
!                  NAtomSpecies,iatom,AtomSpecies(iatom)
!            END IF
!         END IF
!         NotSatisfied=NAtomSpecies<NRemove
!      END DO
!      
!    CLOSE(101)  
   END SUBROUTINE InsertAtomsRandomly
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   
   
   SUBROUTINE MarkAtomsInSpace(AL,Box,mask)
      !all atoms lying inside the box get marked
      ! assumed that periodic operations have been done before hand
      ! box is dimension 6 written as (xmax,ymax,zmax,xmin,ymin,zmin)
      IMPLICIT NONE
      REAL(dp), DIMENSION(6) :: Box
      REAL(dp), DIMENSION(3) :: coord
      TYPE(SystemContainer), POINTER :: AL
      LOGICAL, DIMENSION(:), POINTER :: mask
      INTEGER :: i
      
      CALL MakeSize(mask,AL%NAtoms)
      
      DO i=1,AL%NAtoms
         coord=AL%AtomCoord(3*i-2:3*i)
         mask(i)= (coord(1)<=box(1) .AND. coord(1)>=box(4) .AND. &
            coord(2)<=box(2) .AND. coord(2)>=box(5) .AND. &
            coord(3)<=box(3) .AND. coord(3)>=box(6))
      END DO
   END SUBROUTINE MarkAtomsInSpace
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Add(AL1,AL) 
      !insert AL1 into AL - used by sbrtn Expand
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL1,AL2,AL
      INTEGER :: r1,r2
      
      IF (AL%NAtoms==0 .OR. AL1%NAtoms==0) THEN
         WRITE(6,*) "$Err: One or more of atomic lists to be joined is empty"
         STOP
      END IF
      
      NULLIFY(AL2)
      CALL Merge(AL1,AL,AL2)
      CALL Delete(AL)
      AL=>AL2
   END SUBROUTINE Add
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Del(AL,mask)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      LOGICAL, DIMENSION(:), POINTER :: mask

      IF (SIZE(mask)==AL%NAtoms) THEN
         CALL Del1(AL,mask)
      ELSEIF (SIZE(mask)==3*AL%NAtoms) THEN
         WRITE(6,*) "$Err>> this del has not been implemented"
         STOP
      ELSE
         WRITE(6,*) "...number of atoms in Atom list while deleting is ",AL%NAtoms
         WRITE(6,*) "$Err>> size of mask passed to del is wierd ",SIZE(mask)
         STOP
      END IF
   END SUBROUTINE Del
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Del1(AL,mask)
      !all atoms with mask .TRUE. are to be deleted
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1
      LOGICAL, DIMENSION(:), POINTER :: mask !mask is of size AL%NAtoms
      LOGICAL, DIMENSION(:), POINTER :: mask1 !mask is of size AL%NAtoms
      INTEGER, DIMENSION(:), POINTER :: ivec
      INTEGER :: NAtoms,NAtoms1,i
      
      NAtoms=AL%NAtoms
      mask=.NOT. mask

      !count # of atoms
      CALL MakeSize(ivec,NAtoms)
      WHERE (mask) 
         ivec=1
      ELSEWHERE
         ivec=0
      END WHERE

      NAtoms1=SUM(ivec)
      DEALLOCATE(ivec)

      WRITE(UNIT=6,FMT='(" ...trimming atom list, keeping ", &
           I5," atoms out of ",I5," atoms")') NAtoms1,NAtoms

      CALL MakeSize(mask1,3*NAtoms)
      mask1=PACK(SPREAD(mask,1,3),.TRUE.) !for use with coordinates etc.
      NULLIFY(AL1)
      CALL MakeSize(AL1,NAtoms1)
      
      AL1%PotentialEnergy=0._dp
      AL1%KineticEnergy=0._dp
      AL1%NAtoms=NAtoms1
      !AL1%NSpeciesType=AL%NSpeciesType
      AL1%BoundaryCondX=AL%BoundaryCondX
      AL1%BoundaryCondY=AL%BoundaryCondY
      AL1%BoundaryCondZ=AL%BoundaryCondZ
      AL1%BoxSize=AL%BoxSize

      !CALL MakeSize(AL1%SpeciesDirectory,SIZE(AL%SpeciesDirectory))
      !AL1%SpeciesDirectory=AL%SpeciesDirectory
      
      AL1%AtomSpecies=PACK(AL%AtomSpecies,Mask)
      AL1%AtomCharge=PACK(AL%AtomCharge,Mask)

      AL1%AtomCoord=PACK(AL%AtomCoord,Mask1)
      AL1%AtomVelocity=PACK(AL%AtomVelocity,Mask1)
      AL1%AtomDriftVelocity=PACK(AL%AtomDriftVelocity,Mask1)
      AL1%AtomForce=PACK(AL%AtomForce,Mask1)
      AL1%AtomMass=PACK(AL%AtomMass,Mask1)
      AL1%AtomShiftCoord=PACK(AL%AtomShiftCoord,Mask1)
      AL1%AtomIsMoving=PACK(AL%AtomIsMoving,Mask1)
      DEALLOCATE(mask,mask1)
      
      CALL Delete(AL)
      AL=>AL1
      NULLIFY(AL1)
   END SUBROUTINE Del1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Expand(AL,n)
      !assuming UnitCell can be expanded in Cartesian coordinates
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: n(3)
      
      IF (n(1)>1) CALL Expand0(AL,n(1),1) !expand along x
      IF (n(2)>1) CALL Expand0(AL,n(2),2) !expand along y
      IF (n(3)>1) CALL Expand0(AL,n(3),3) !expand along z
      
   END SUBROUTINE Expand
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Expand0(AL,n,dimn)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1
      REAL(dp) :: t(3)
      INTEGER :: n,dimn,i

      NULLIFY(AL1)
      !t=0._dp !translation
      !t(dimn)=AL%BoxSize(dimn)-AL%BoxSize(dimn+3) !only along dimn of interest
      t=Basis(3*dimn-2:3*dimn)

      CALL MakeSize(AL1,AL%NAtoms)
      CALL Copy(AL,AL1) !AL1 is the copy which is translated & added

      !CALL MakeSize(AL1%SpeciesDirectory,SIZE(AL%SpeciesDirectory))
      !AL1%SpeciesDirectory=AL%SpeciesDirectory
      DO i=2,n
         CALL Translate(AL1,t) !translate all atoms by t*i
         CALL Add(AL1,AL) !merge the two AL
         AL%BoxSize(dimn)=AL%BoxSize(dimn)+AL1%BoxSize(dimn)-AL1%BoxSize(dimn+3)
      END DO
      
      CALL Delete(AL1)
   END SUBROUTINE Expand0
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE FlipCrystal(AL,plane) !reflection about a plane
      IMPLICIT NONE
      !see Cleave for background on equation of plane
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), INTENT(IN) :: plane(3,3)
      REAL(dp), DIMENSION(:), POINTER :: A
      REAL(dp) :: D,coord(3),mat(3,3),Anorm,Aunit(3),distance
      INTEGER :: i,r1,r2
      
      !A is normal
      CALL MakeSize(A,3)
      mat=plane; mat(:,1)=1; A(1)=DET(mat)
      mat=plane; mat(:,2)=1; A(2)=DET(mat)
      mat=plane; mat(:,3)=1; A(3)=DET(mat)
      mat=plane; D=DET(mat)
      
      Anorm=NORM(A,3)
      Aunit=A/Anorm
      r1=-2; r2=0
      DO i=1,AL%NAtoms
         r1=r1+3; r2=r2+3
         coord=AL%AtomCoord(r1:r2)
         distance=(DOT_PRODUCT(coord,A)-D)/Anorm ; !signed distance from the plane
         AL%AtomCoord(r1:r2)=coord-distance*Aunit !reflect
      END DO
      DEALLOCATE(A)
   END SUBROUTINE FlipCrystal
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MoveBoxRange(AL,dimn,maxbox,minbox,maxboxincr,minboxincr)
   !Modify the box size along a certain cartesian dimension 
   
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: dimn
      REAL(dp), OPTIONAL :: maxbox,minbox,maxboxincr,minboxincr
      REAL(dp) :: maxbox1,minbox1,maxboxincr1,minboxincr1
      LOGICAL :: imx1,imx2,imxr1,imxr2
      
      imx1=PRESENT(maxbox)
      imx2=PRESENT(minbox)
      imxr1=PRESENT(maxboxincr)
      imxr2=PRESENT(minboxincr)
       
      IF ((imx1 .OR. imx2) .AND. (imxr1 .OR. imxr2)) THEN
         WRITE(6,*) "$Err>> Both new box size and box increment provided to MoveBoxRange"
         STOP
      END IF
      
      IF (imx1) THEN
         AL%BoxSize(dimn)=maxbox
      END IF
      IF (imx2) THEN
         AL%BoxSize(dimn+3)=minbox
      END IF
      IF (imxr1) THEN
         AL%BoxSize(dimn)=AL%BoxSize(dimn)+maxboxincr
      END IF
      IF (imxr2) THEN
         AL%BoxSize(dimn+3)=AL%BoxSize(dimn+3)+minboxincr
      END IF
      
      WRITE(6,*) "New box size:"
      WRITE(UNIT=*,FMT='(" Maximum...",3F10.3)') AL%BoxSize(1:3)
      WRITE(UNIT=*,FMT='(" Minimum...",3F10.3)') AL%BoxSize(4:6)
   END SUBROUTINE MoveBoxRange
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE UnitCell(UC,CrystalType) !return unitcell information as desired
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      CHARACTER(len=10) :: CrystalType
      
      SELECT CASE (CrystalType(1:3))
      CASE("scu") !simple cubic
         CALL sc(UC)
      CASE("fcc")
         !CALL fcc(UC)
         SELECT CASE (CrystalType(4:6))
         CASE ("100"); CALL fcc100(UC)
         CASE ("110"); CALL fcc110(UC)
         CASE ("111"); CALL fcc111(UC)
         CASE ("210"); CALL fcc210(UC)
         CASE ("910"); CALL fcc910(UC)
         CASE DEFAULT; CALL fcc100(UC) !bulk
         END SELECT
      CASE("bcc")
         CALL bcc(UC)
         !SELECT CASE (CrystalType(4:6))
         !CASE ("100"); CALL bcc100(UC)
         !CASE ("110"); CALL bcc110(UC)
         !CASE ("111"); CALL bcc111(UC)
         !CASE DEFAULT; CALL bcc(UC) !bulk
         !END SELECT
      CASE("hcp")
         CALL hcp(UC)
      CASE("hex")
         CALL hex(UC)
      CASE("dia")
         CALL diamond(UC)
      CASE("gra")
         CALL graphite(UC)
      CASE("per")
         CALL Perovskite(UC)
      CASE("nac") !NaCl
         CALL NaCl(UC)
      CASE("csc") !CsCl
         CALL CsCl(UC)
      CASE("caf") !CaF2
         CALL CaF2(UC)
      CASE("rut") !Rutile
         CALL Rutile(UC)
      CASE("Oli") !Olivine
         CALL Olivine(UC)
      CASE("Cha") !creates a chain of alternating positive and negative ions
         CALL ChainIons(UC)
      CASE("LF1") !LiFePO4
         CALL LiFePO4(UC)
      CASE("LF2") !Li0.5FePO4
         CALL LiHalfFePO4(UC)
      CASE("LF3") !FePO4
         CALL FePO4(UC)
      CASE("Te1")
         CALL ZrO2Tetragonal(UC)
      CASE("Si1") !fcc+half interstial si
         CALL Silicon(UC)
      CASE("wur")
         CALL Wurtzite(UC)
      CASE("ZnS")
         CALL ZnS(UC)
      CASE("None")
         RETURN
      CASE DEFAULT
         WRITE(6,*) "Err:Unknown crystal structure passed"
         STOP
      END SELECT
      
   END SUBROUTINE UnitCell
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE sc(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 1
      !Examples: 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,1) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=1
      UC%NMove(1:3)=1
      !UC%NSpeciesType=1
      !CALL MakeSize(UC%SpeciesDirectory,1)
      !UC%SpeciesDirectory=1 !replace with atomic # later

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/)
      
      UC%AtomSpecies=1
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      
      Basis(1:3)=(/1.0_dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1.0_dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1.0_dp/)
   END SUBROUTINE sc
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE fcc100(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 1 
      !Examples: 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,4) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=4
      UC%NMove(1:3)=4
      !UC%NSpeciesType=1
      !CALL MakeSize(UC%SpeciesDirectory,1)
      !UC%SpeciesDirectory=1

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/)
      UC%AtomCoord(4:6)=(/0.5_dp,0.0_dp,0.5_dp/)
      UC%AtomCoord(7:9)=(/0.5_dp,0.5_dp,0.0_dp/)
      UC%AtomCoord(10:12)=(/0.0_dp,0.5_dp,0.5_dp/)
      
      UC%AtomSpecies=1
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      
      Basis(1:3)=(/1.0_dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1.0_dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1.0_dp/)
   END SUBROUTINE fcc100
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE fcc110(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 1
      !Examples: 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      REAL(dp) :: sqrt2inv
      
      CALL MakeSize(UC,2) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=2
      UC%NMove(1:3)=2
      !UC%NSpeciesType=1
      !CALL MakeSize(UC%SpeciesDirectory,1)
      !UC%SpeciesDirectory=1 !replace with atomic # later

      sqrt2inv=1._dp/SQRT(2._dp)
      UC%BoxSize=(/sqrt2inv,1.0_dp,sqrt2inv,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/)
      UC%AtomCoord(4:6)=0.5_dp*(/sqrt2inv,1.0_dp,sqrt2inv/)
      
      UC%AtomSpecies=1
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      
      Basis(1:3)=(/sqrt2inv,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1.0_dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,sqrt2inv/)
   END SUBROUTINE fcc110
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE fcc111(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 1
      !Examples: 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,6) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=6
      UC%NMove(1:3)=6
      
      UC%BoxSize=0._dp
      UC%BoxSize(1)=0.707106781186547_dp
      UC%BoxSize(2)=1.22474487139159_dp
      UC%BoxSize(3)=1.73205080756888_dp
      UC%AtomCoord(1:3)=(/0.35355339_dp,0.20412415_dp,0.00000000_dp/)
      UC%AtomCoord(4:6)=(/0.00000000_dp,0.81649658_dp,0.00000000_dp/)
      UC%AtomCoord(7:9)=(/0.00000000_dp,0.00000000_dp,-0.57735027_dp/)
      UC%AtomCoord(10:12)=(/0.35355339_dp,0.61237244_dp,-0.57735027_dp/)
      UC%AtomCoord(13:15)=(/0.00000000_dp,0.40824829_dp,-1.15470054_dp/)
      UC%AtomCoord(16:18)=(/.35355339_dp,1.02062073_dp,-1.15470054_dp/)
      !UC%AtomCoord(1:3)=(/0._dp,0._dp,0._dp/)
      !UC%AtomCoord(4:6)=(/0.000000000000000_dp,0.408248290463863_dp,1.154700538379252_dp/)
      !UC%AtomCoord(7:9)=(/0.000000000000000_dp,0.204124145231932_dp,0.577350269189626_dp/)
      !UC%AtomCoord(10:12)=(/0.353553390593274_dp,0.816496580927726_dp,0.577350269189626_dp/)
      !UC%AtomCoord(13:15)=(/0.353553390593274_dp,0.612372435695794_dp,0.000000000000000_dp/)
      !UC%AtomCoord(16:18)=(/0.353553390593274_dp,1.020620726159658_dp,1.154700538379252_dp/)

      !UC%BoxSize=0._dp
      !UC%BoxSize(1)=0.707106781186547_dp
      !UC%BoxSize(2)=0.612372435695794_dp
      !UC%BoxSize(3)=0.577350269189626_dp
      !UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/)
      
      UC%AtomSpecies=1
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      
      Basis(1:3)=(/0.707106781186547_dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1.22474487139159_dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1.73205080756888_dp/)
      !Basis(1:3)=(/0.707106781186547_dp,0.0_dp,0.0_dp/) !repeat along x-direction
      !Basis(4:6)=(/0.353553390593274_dp,0.612372435695794_dp,0.0_dp/) !repeat along y-direction
      !Basis(7:9)=(/0.353553390593274_dp,0.204124145231932_dp,0.577350269189626_dp/)
   END SUBROUTINE fcc111
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE fcc210(UC)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,20) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=20
      UC%NMove(1:3)=20
      
      UC%BoxSize=0._dp
      UC%BoxSize(1)=2.23606798_dp
      UC%BoxSize(2)=1.00000000_dp
      UC%BoxSize(3)=2.23606798_dp
      UC%AtomCoord(1:3)  =(/0.89442719_dp,0.00000000_dp,0.00000000_dp/)
      UC%AtomCoord(4:6)  =(/2.01246118_dp,0.50000000_dp,0.00000000_dp/)
      UC%AtomCoord(7:9)  =(/1.56524758_dp,0.00000000_dp,-0.22360680_dp/)
      UC%AtomCoord(10:12)=(/0.44721360_dp,0.50000000_dp,-0.22360680_dp/)
      UC%AtomCoord(13:15)=(/0.00000000_dp,0.00000000_dp,-0.44721360_dp/)
      UC%AtomCoord(16:18)=(/1.11803399_dp,0.50000000_dp,-0.44721360_dp/)
      UC%AtomCoord(19:21)=(/0.67082039_dp,0.00000000_dp,-0.67082039_dp/)
      UC%AtomCoord(22:24)=(/1.78885438_dp,0.50000000_dp,-0.67082039_dp/)
      UC%AtomCoord(25:27)=(/1.34164079_dp,0.00000000_dp,-0.89442719_dp/)
      UC%AtomCoord(28:30)=(/0.22360680_dp,0.50000000_dp,-0.89442719_dp/)
      UC%AtomCoord(31:33)=(/2.01246118_dp,0.00000000_dp,-1.11803399_dp/)
      UC%AtomCoord(34:36)=(/0.89442719_dp,0.50000000_dp,-1.11803399_dp/)
      UC%AtomCoord(37:39)=(/0.44721360_dp,0.00000000_dp,-1.34164079_dp/)
      UC%AtomCoord(40:42)=(/1.56524758_dp,0.50000000_dp,-1.34164079_dp/)
      UC%AtomCoord(43:45)=(/1.11803399_dp,0.00000000_dp,-1.56524758_dp/)
      UC%AtomCoord(46:48)=(/0.00000000_dp,0.50000000_dp,-1.56524758_dp/)
      UC%AtomCoord(49:51)=(/1.78885438_dp,0.00000000_dp,-1.78885438_dp/)
      UC%AtomCoord(52:54)=(/0.67082039_dp,0.50000000_dp,-1.78885438_dp/)
      UC%AtomCoord(55:57)=(/0.22360680_dp,0.00000000_dp,-2.01246118_dp/)
      UC%AtomCoord(58:60)=(/1.34164079_dp,0.50000000_dp,-2.01246118_dp/)
      
      UC%AtomSpecies=1
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      
      Basis(1:3)=(/2.23606798_dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1.0_dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,2.23606798_dp/)
   END SUBROUTINE fcc210
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE fcc910(UC)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,82) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=82
      UC%NMove(1:3)=82
      
      UC%BoxSize=0._dp
      UC%BoxSize(1)=4.52769257_dp
      UC%BoxSize(2)=1.00000000_dp
      UC%BoxSize(3)=4.52769257_dp
      
      UC%AtomCoord=(/3.03686697_dp,0.00000000_dp, 0.00000000_dp, &
         2.53992510_dp,0.50000000_dp,-0.05521576_dp,2.04298323_dp,0.00000000_dp,-0.11043153_dp, &
         1.54604137_dp,0.50000000_dp,-0.16564729_dp,1.04909950_dp,0.00000000_dp,-0.22086305_dp, &
         0.55215763_dp,0.50000000_dp,-0.27607882_dp,0.05521576_dp,0.00000000_dp,-0.33129458_dp, &
         4.08596646_dp,0.50000000_dp,-0.38651034_dp,3.58902460_dp,0.00000000_dp,-0.44172610_dp, &
         3.09208273_dp,0.50000000_dp,-0.49694187_dp,2.59514086_dp,0.00000000_dp,-0.55215763_dp, &
         2.09819900_dp,0.50000000_dp,-0.60737339_dp,1.60125713_dp,0.00000000_dp,-0.66258916_dp, &
         1.10431526_dp,0.50000000_dp,-0.71780492_dp,0.60737339_dp,0.00000000_dp,-0.77302068_dp, &
         0.11043153_dp,0.50000000_dp,-0.82823645_dp,4.14118223_dp,0.00000000_dp,-0.88345221_dp, &
         3.64424036_dp,0.50000000_dp,-0.93866797_dp,3.14729849_dp,0.00000000_dp,-0.99388373_dp, &
         2.65035663_dp,0.50000000_dp,-1.04909950_dp,2.15341476_dp,0.00000000_dp,-1.10431526_dp, &
         1.65647289_dp,0.50000000_dp,-1.15953102_dp,1.15953102_dp,0.00000000_dp,-1.21474679_dp, &
         0.66258916_dp,0.50000000_dp,-1.26996255_dp,0.16564729_dp,0.00000000_dp,-1.32517831_dp, &
         4.19639799_dp,0.50000000_dp,-1.38039408_dp,3.69945612_dp,0.00000000_dp,-1.43560984_dp, &
         3.20251426_dp,0.50000000_dp,-1.49082560_dp,2.70557239_dp,0.00000000_dp,-1.54604137_dp, &
         2.20863052_dp,0.50000000_dp,-1.60125713_dp,1.71168865_dp,0.00000000_dp,-1.65647289_dp, &
         1.21474679_dp,0.50000000_dp,-1.71168865_dp,0.71780492_dp,0.00000000_dp,-1.76690442_dp, &
         0.22086305_dp,0.50000000_dp,-1.82212018_dp,4.25161375_dp,0.00000000_dp,-1.87733594_dp, &
         3.75467189_dp,0.50000000_dp,-1.93255171_dp,3.25773002_dp,0.00000000_dp,-1.98776747_dp, &
         2.76078815_dp,0.50000000_dp,-2.04298323_dp,2.26384628_dp,0.00000000_dp,-2.09819900_dp, &
         1.76690442_dp,0.50000000_dp,-2.15341476_dp,1.26996255_dp,0.00000000_dp,-2.20863052_dp, &
         0.77302068_dp,0.50000000_dp,-2.26384628_dp,0.27607882_dp,0.00000000_dp,-2.31906205_dp, &
         4.30682952_dp,0.50000000_dp,-2.37427781_dp,3.80988765_dp,0.00000000_dp,-2.42949357_dp, &
         3.31294578_dp,0.50000000_dp,-2.48470934_dp,2.81600391_dp,0.00000000_dp,-2.53992510_dp, &
         2.31906205_dp,0.50000000_dp,-2.59514086_dp,1.82212018_dp,0.00000000_dp,-2.65035663_dp, &
         1.32517831_dp,0.50000000_dp,-2.70557239_dp,0.82823645_dp,0.00000000_dp,-2.76078815_dp, &
         0.33129458_dp,0.50000000_dp,-2.81600391_dp,4.36204528_dp,0.00000000_dp,-2.87121968_dp, &
         3.86510341_dp,0.50000000_dp,-2.92643544_dp,3.36816155_dp,0.00000000_dp,-2.98165120_dp, &
         2.87121968_dp,0.50000000_dp,-3.03686697_dp,2.37427781_dp,0.00000000_dp,-3.09208273_dp, &
         1.87733594_dp,0.50000000_dp,-3.14729849_dp,1.38039408_dp,0.00000000_dp,-3.20251426_dp, &
         0.88345221_dp,0.50000000_dp,-3.25773002_dp,0.38651034_dp,0.00000000_dp,-3.31294578_dp, &
         4.41726104_dp,0.50000000_dp,-3.36816155_dp,3.92031918_dp,0.00000000_dp,-3.42337731_dp, &
         3.42337731_dp,0.50000000_dp,-3.47859307_dp,2.92643544_dp,0.00000000_dp,-3.53380883_dp, &
         2.42949357_dp,0.50000000_dp,-3.58902460_dp,1.93255171_dp,0.00000000_dp,-3.64424036_dp, &
         1.43560984_dp,0.50000000_dp,-3.69945612_dp,0.93866797_dp,0.00000000_dp,-3.75467189_dp, &
         0.44172610_dp,0.50000000_dp,-3.80988765_dp,4.47247681_dp,0.00000000_dp,-3.86510341_dp, &
         3.97553494_dp,0.50000000_dp,-3.92031918_dp,3.47859307_dp,0.00000000_dp,-3.97553494_dp, &
         2.98165120_dp,0.50000000_dp,-4.03075070_dp,2.48470934_dp,0.00000000_dp,-4.08596646_dp, &
         1.98776747_dp,0.50000000_dp,-4.14118223_dp,1.49082560_dp,0.00000000_dp,-4.19639799_dp, &
         0.99388373_dp,0.50000000_dp,-4.25161375_dp,0.49694187_dp,0.00000000_dp,-4.30682952_dp, &
         0.00000000_dp,0.50000000_dp,-4.36204528_dp,4.03075070_dp,0.00000000_dp,-4.41726104_dp, &
         3.53380883_dp,0.50000000_dp,-4.47247681_dp/)
      
      UC%AtomSpecies=1
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      
      Basis(1:3)=(/4.52769257_dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1.0_dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,4.52769257_dp/)
   END SUBROUTINE fcc910
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE bcc(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 1
      !Examples:
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating bcc-bulk crystal"
      
      CALL MakeSize(UC,2) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=2
      UC%NMove(1:3)=2
      !UC%NSpeciesType=1
      !CALL MakeSize(UC%SpeciesDirectory,1)
      !UC%SpeciesDirectory=1

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/)
      UC%AtomCoord(4:6)=(/0.5_dp,0.5_dp,0.5_dp/)
      
      UC%AtomSpecies=1
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      
      Basis(1:3)=(/1.0_dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1.0_dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1.0_dp/)
   END SUBROUTINE bcc
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE hcp(UC) !A3 crystal structure
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1, i.e., a=c=1
      !# atom species present = 1
      !Examples: 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating hcp-bulk crystal"
      
      CALL MakeSize(UC,4) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=4
      UC%NMove(1:3)=4
      !UC%NSpeciesType=1
      !CALL MakeSize(UC%SpeciesDirectory,1)
      !UC%SpeciesDirectory=1 !replace with atomic # later

      UC%BoxSize=0._dp
      UC%BoxSize(1)=0.707106781186547_dp
      UC%BoxSize(2)=1.22474487139159_dp
      UC%BoxSize(3)=1.15470053837925_dp
      UC%AtomCoord(1:3)=(/0.353553390593274_dp,0.0_dp,0.0_dp/)
      UC%AtomCoord(4:6)=(/0._dp,0.612372435695794_dp,0._dp/)
      UC%AtomCoord(7:9)=(/0._dp,0.204124145231932_dp,-0.577350269189626_dp/)
      UC%AtomCoord(10:12)=(/0.353553390593274_dp,0.816496580927726_dp,-0.577350269189626_dp/)
      
      UC%AtomSpecies=1
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      
      Basis(1:3)=(/0.707106781186547_dp,0.0_dp,0.0_dp/) !repeat along x-direction
      Basis(4:6)=(/0._dp,1.22474487139159_dp,0.0_dp/) !repeat along y-direction
      Basis(7:9)=(/0.0_dp,0.0_dp,1.15470053837925_dp/)
   END SUBROUTINE hcp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE hex(UC)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
   END SUBROUTINE hex
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Wurtzite(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 2
      !Examples: CdS,GaN,SiC,ZnS
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating Wurtzite bulk crystal"

      CALL MakeSize(UC,4)
      UC%NAtoms=4
      UC%NMove(1:3)=4
      !UC%NSpeciesType=2
      !CALL MakeSize(UC%SpeciesDirectory,2)
      !UC%SpeciesDirectory=(/1,2/) !replace with atomic # later

      !considering parallelogram parallelopipe a cuboid structure to generate unit cell
      UC%BoxSize=(/1._dp,SQRT(3._dp)/2._dp,2._dp*SQRT(2._dp/3._dp),0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.00000_dp,0.00000_dp,0.00000_dp/) ! base atom along x Cd atoms
      UC%AtomCoord(4:6)=(/0.50000_dp,1/2._dp/SQRT(3._dp),1/2._dp/SQRT(6._dp)/) !tetrahedral atom second plane S
      UC%AtomCoord(7:9)=(/0.50000_dp,1/2._dp/SQRT(3._dp),SQRT(2._dp/3._dp)/) ! tetrahedral atom third plane Cd
      UC%AtomCoord(10:12)=(/0.00000_dp,0.00000_dp,SQRT(2._dp/3._dp)+1/2._dp/SQRT(6._dp)/) ! fourth plane S atoms
      
      
      UC%AtomSpecies(1)=1
      UC%AtomSpecies(2)=2
      UC%AtomSpecies(3)=1
      UC%AtomSpecies(4)=2
      
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.5_dp,SQRT(3._dp)/2._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,2._dp*SQRT(2._dp/3._dp)/)
   END SUBROUTINE Wurtzite
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ZnS(UC)
      !Zinc Blende structure is similar to silicon diamond cubice structure (tetrahedral positions are occupied by sulphur)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 2 Zinc and Sulphur
      !Examples: 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating ZnS-bulk crystal"
      
      CALL MakeSize(UC,8) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=8
      UC%NMove=8
      !UC%NSpeciesType=2
      !CALL MakeSize(UC%SpeciesDirectory,2)
      !UC%SpeciesDirectory=(/1,2/) !replace with atomic # later

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/) !Zn
      UC%AtomCoord(4:6)=(/0.0_dp,0.5_dp,0.5_dp/) !Zn
      UC%AtomCoord(7:9)=(/0.5_dp,0.5_dp,0.0_dp/) !Zn
      UC%AtomCoord(10:12)=(/0.5_dp,0.0_dp,0.5_dp/) !Zn
      UC%AtomCoord(13:15)=(/0.75_dp,0.25_dp,0.75_dp/) !int S
      UC%AtomCoord(16:18)=(/0.25_dp,0.75_dp,0.75_dp/) !int S
      UC%AtomCoord(19:21)=(/0.25_dp,0.25_dp,0.25_dp/) !int S
      UC%AtomCoord(22:24)=(/0.75_dp,0.75_dp,0.25_dp/) !int S

      UC%AtomSpecies(1:4)=1
      UC%AtomSpecies(5:8)=2
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE ZnS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CsCl(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 1
      !Examples:
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating CsCl-bulk crystal"

      CALL MakeSize(UC,2) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=2
      UC%NMove(1:3)=2
      !UC%NSpeciesType=2
      !CALL MakeSize(UC%SpeciesDirectory,2)
      !UC%SpeciesDirectory=(/1,2/)

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/) !cation
      UC%AtomCoord(4:6)=(/0.5_dp,0.5_dp,0.5_dp/) !anion
      
      UC%AtomSpecies=(/1,2/)
      UC%AtomCharge=(/1.0_dp,-1.0_dp/)
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/) !repeat along x-direction
      Basis(4:6)=(/0._dp,1._dp,0.0_dp/) !repeat along y-direction
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE CsCl
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ChainIons(UC)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating a 1D chain of ions"
      
      CALL MakeSize(UC,2) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=2
      UC%NMove(1:3)=2
      !UC%NSpeciesType=2
      !CALL MakeSize(UC%SpeciesDirectory,2)
      !UC%SpeciesDirectory=(/1,2/)

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/) !cation
      UC%AtomCoord(4:6)=(/0.5_dp,0.0_dp,0.0_dp/) !anion
      
      UC%AtomSpecies=(/1,2/)
      UC%AtomCharge=(/1.0_dp,-1.0_dp/)
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE ChainIons
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE NaCl(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 2
      !Examples: 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating NaCl-bulk crystal"

      CALL MakeSize(UC,8)
      UC%NAtoms=8
      UC%NMove(1:3)=8
      !UC%NSpeciesType=2
      !CALL MakeSize(UC%SpeciesDirectory,2)
      !UC%SpeciesDirectory=(/1,2/) !replace with atomic # later

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)  =(/0.0_dp,0.0_dp,0.0_dp/) !cation
      UC%AtomCoord(4:6)  =(/0.5_dp,0.5_dp,0.0_dp/) !cation
      UC%AtomCoord(7:9)  =(/0.0_dp,0.5_dp,0.5_dp/) !cation
      UC%AtomCoord(10:12)=(/0.5_dp,0.0_dp,0.5_dp/) !cation
      UC%AtomCoord(13:15)=(/0.5_dp,0.0_dp,0.0_dp/) ! anion
      UC%AtomCoord(16:18)=(/0.0_dp,0.5_dp,0.0_dp/) ! anion
      UC%AtomCoord(19:21)=(/0.0_dp,0.0_dp,0.5_dp/) ! anion
      UC%AtomCoord(22:24)=(/0.5_dp,0.5_dp,0.5_dp/) ! anion

      UC%AtomCharge(1:4)=1._dp
      UC%AtomCharge(5:8)=-1._dp
      UC%AtomSpecies(1:4)=1
      UC%AtomSpecies(5:8)=2
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE NaCl
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SiC(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 1
      !Examples: 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating SiC-bulk crystal"
      
      CALL MakeSize(UC,8) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=8
      UC%NMove(1:3)=8
      !UC%NSpeciesType=2
      !CALL MakeSize(UC%SpeciesDirectory,2)
      !UC%SpeciesDirectory=(/1,2/) !replace with atomic # later

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/) !Si
      UC%AtomCoord(4:6)=(/0.75_dp,0.25_dp,0.75_dp/) !C
      UC%AtomCoord(7:9)=(/0.25_dp,0.75_dp,0.75_dp/) !C
      UC%AtomCoord(10:12)=(/0.25_dp,0.25_dp,0.25_dp/) !C
      UC%AtomCoord(13:15)=(/0.75_dp,0.75_dp,0.25_dp/) !C
      UC%AtomCoord(16:18)=(/0.0_dp,0.5_dp,0.5_dp/) !Si
      UC%AtomCoord(19:21)=(/0.5_dp,0.5_dp,0.0_dp/) !Si
      UC%AtomCoord(22:24)=(/0.5_dp,0.0_dp,0.5_dp/) !Si

      UC%AtomSpecies=(/1,2,2,2,2,1,1,1/)
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE SiC
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CaF2(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !cation fcc sites. anion tetrahedral sites
      !Lattice constant = 1
      !# atom species present = 1
      !Examples:
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,12) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=12
      UC%NMove=12
      !UC%NSpeciesType=2
      !CALL MakeSize(UC%SpeciesDirectory,2)
      !UC%SpeciesDirectory=(/1,2/)

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)  =(/0.0_dp,0.0_dp,0.0_dp/) !cation
      UC%AtomCoord(4:6)  =(/0.5_dp,0.5_dp,0.0_dp/) !cation
      UC%AtomCoord(7:9)  =(/0.0_dp,0.5_dp,0.5_dp/) !cation
      UC%AtomCoord(10:12)=(/0.5_dp,0.0_dp,0.5_dp/) !cation
      UC%AtomCoord(13:15)=(/0.25_dp,0.25_dp,0.25_dp/) ! anion
      UC%AtomCoord(16:18)=(/0.25_dp,0.75_dp,0.25_dp/) ! anion
      UC%AtomCoord(19:21)=(/0.75_dp,0.25_dp,0.25_dp/) ! anion
      UC%AtomCoord(22:24)=(/0.75_dp,0.75_dp,0.25_dp/) ! anion
      UC%AtomCoord(25:27)=(/0.25_dp,0.25_dp,0.75_dp/) ! anion
      UC%AtomCoord(28:30)=(/0.25_dp,0.75_dp,0.75_dp/) ! anion
      UC%AtomCoord(31:33)=(/0.75_dp,0.25_dp,0.75_dp/) ! anion
      UC%AtomCoord(34:36)=(/0.75_dp,0.75_dp,0.75_dp/) ! anion

      UC%AtomSpecies(1:4)=1
      UC%AtomSpecies(5:12)=2
      UC%AtomCharge(1:4)=2._dp
      UC%AtomCharge(5:12)=-1._dp
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE CaF2
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CaFx(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !cation fcc sites. anion tetrahedral sites
      !Lattice constant = 1
      !# atom species present = 1
      !Examples:
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,12) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=12
      UC%NMove=12
      !UC%NSpeciesType=2
      !CALL MakeSize(UC%SpeciesDirectory,2)
      !UC%SpeciesDirectory=(/1,2/)

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)  =(/0.0_dp,0.0_dp,0.0_dp/) !cation
      UC%AtomCoord(4:6)  =(/0.5_dp,0.5_dp,0.0_dp/) !cation
      UC%AtomCoord(7:9)  =(/0.0_dp,0.5_dp,0.5_dp/) !cation
      UC%AtomCoord(10:12)=(/0.5_dp,0.0_dp,0.5_dp/) !cation
      UC%AtomCoord(13:15)=(/0.25_dp,0.25_dp,0.25_dp/) ! anion
      UC%AtomCoord(16:18)=(/0.25_dp,0.75_dp,0.25_dp/) ! anion
      UC%AtomCoord(19:21)=(/0.75_dp,0.25_dp,0.25_dp/) ! anion
      UC%AtomCoord(22:24)=(/0.75_dp,0.75_dp,0.25_dp/) ! anion
      UC%AtomCoord(25:27)=(/0.25_dp,0.25_dp,0.75_dp/) ! anion
      UC%AtomCoord(28:30)=(/0.25_dp,0.75_dp,0.75_dp/) ! anion
      UC%AtomCoord(31:33)=(/0.75_dp,0.25_dp,0.75_dp/) ! anion
      UC%AtomCoord(34:36)=(/0.75_dp,0.75_dp,0.75_dp/) ! anion

      UC%AtomSpecies(1:4)=1
      UC%AtomSpecies(5:12)=2
      UC%AtomCharge(1:4)=2._dp
      UC%AtomCharge(5:12)=-1._dp
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE CaFx
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Olivine(UC) 
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !e.g., LiFePO4
      !Lattice constant = 1
      !# atom species present = 1
      !Examples:
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,6) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=6
      UC%NMove=6
      !UC%NSpeciesType=4
      !CALL MakeSize(UC%SpeciesDirectory,4)
      !UC%SpeciesDirectory=(/1,2,3,4/)

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)  =(/0.0_dp,0.0_dp,0.0_dp/) !cation, e.g., Li+
      UC%AtomCoord(4:6)  =(/0.29_dp,0.25_dp,0.43_dp/) !cation, e.g. Fe2+
      UC%AtomCoord(7:9)  =(/0.09_dp,0.25_dp,0.75_dp/) !cation, e.g., P5+
      UC%AtomCoord(10:12)=(/0.09_dp,0.25_dp,0.19_dp/) !anion, e.g., O2-
      UC%AtomCoord(13:15)=(/0.17_dp,0.05_dp,0.30_dp/) !anion, e.g., O2-
      UC%AtomCoord(16:18)=(/0.45_dp,0.25_dp,0.19_dp/) !anion, e.g., O2-

      UC%AtomSpecies(1)=1 
      UC%AtomSpecies(2)=2
      UC%AtomSpecies(3)=3 
      UC%AtomSpecies(4:6)=4
      
      UC%AtomCharge(1)=1._dp
      UC%AtomCharge(2)=2._dp
      UC%AtomCharge(3)=5._dp
      UC%AtomCharge(4:6)=-2._dp
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE Olivine
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!   SUBROUTINE LiFePO4(UC)
!   !Unit cell for LiFePO4 obtained from Phys. Rev. B 68, 165107 (2003) - originally! from Chem. Mater. 13 1570 (2001)
!      IMPLICIT NONE
!      TYPE(SystemContainer), POINTER :: UC
      
!      CALL MakeSize(UC,6) !Elements not allocated @ this point - SpeciesDirectory
!      UC%NAtoms=6
!      UC%NMove=6
!      !UC%BoxSize=(/10.06_dp,5.89_dp,4.65_dp,0._dp,0._dp,0._dp/)
!      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      
!      UC%AtomCoord(1:3)  =(/0.0_dp,0.0_dp,0.0_dp/) !cation, e.g., Li+
!      UC%AtomCoord(4:6)  =(/0.28_dp,0.25_dp,0.98_dp/) !cation, e.g. Fe2+
!      UC%AtomCoord(7:9)  =(/0.09_dp,0.25_dp,0.42_dp/) !cation, e.g., P5+
!      UC%AtomCoord(10:12)=(/0.09_dp,0.25_dp,0.75_dp/) !anion, e.g., O2-
!      UC%AtomCoord(13:15)=(/0.45_dp,0.25_dp,0.21_dp/) !anion, e.g., O2-
!      UC%AtomCoord(16:18)=(/0.17_dp,0.04_dp,0.29_dp/) !anion, e.g., O2-

!      UC%AtomSpecies(1)=1 
!      UC%AtomSpecies(2)=2
!      UC%AtomSpecies(3)=3 
!      UC%AtomSpecies(4:6)=4
      
!      UC%AtomCharge(1)=1._dp
!      UC%AtomCharge(2)=2._dp
!      UC%AtomCharge(3)=5._dp
!      UC%AtomCharge(4:6)=-2._dp
!      UC%AtomIsMoving=.TRUE.
!      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
!      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
!      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
!   END SUBROUTINE LiFePO4
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    SUBROUTINE ZrO2Tetragonal(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !cation fcc sites. anion tetrahedral sites
      !Lattice constant = 1
      !# atom species present = 1
      !Examples:
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,12) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=12
      UC%NMove=12
      !UC%NSpeciesType=2
      !CALL MakeSize(UC%SpeciesDirectory,2)
      !UC%SpeciesDirectory=(/1,2/)

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/) !cation
      UC%AtomCoord(4:6)=(/0.5_dp,0.0_dp,0.5_dp/) !cation
      UC%AtomCoord(7:9)=(/0.5_dp,0.5_dp,0.0_dp/) !cation
      UC%AtomCoord(10:12)=(/0.0_dp,0.5_dp,0.5_dp/) !cation
      UC%AtomCoord(13:15)=(/0.75_dp,0.25_dp,0.3074_dp/) ! anion
      UC%AtomCoord(16:18)=(/0.25_dp,0.75_dp,0.3047_dp/) ! anion
      UC%AtomCoord(19:21)=(/0.75_dp,0.75_dp,0.1926_dp/) ! anion
      UC%AtomCoord(22:24)=(/0.25_dp,0.25_dp,0.1926_dp/) ! anion
      UC%AtomCoord(25:27)=(/0.75_dp,0.25_dp,0.8074_dp/) ! anion
      UC%AtomCoord(28:30)=(/0.25_dp,0.75_dp,0.8074_dp/) ! anion
      UC%AtomCoord(31:33)=(/0.75_dp,0.75_dp,0.6926_dp/) ! anion
      UC%AtomCoord(34:36)=(/0.25_dp,0.25_dp,0.6926_dp/) ! anion

      UC%AtomSpecies(1:4)=1
      UC%AtomSpecies(5:12)=2
      UC%AtomCharge(1:4)=4._dp
      UC%AtomCharge(5:12)=-2._dp
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE ZrO2Tetragonal
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Silicon(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 1
      !Examples: 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating Si-bulk crystal"
      
      CALL MakeSize(UC,8) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=8
      UC%NMove=8
      !UC%NSpeciesType=2
      !CALL MakeSize(UC%SpeciesDirectory,2)
      !UC%SpeciesDirectory=(/1,2/) !replace with atomic # later

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/) !Si
      UC%AtomCoord(4:6)=(/0.0_dp,0.5_dp,0.5_dp/) !Si
      UC%AtomCoord(7:9)=(/0.5_dp,0.5_dp,0.0_dp/) !Si
      UC%AtomCoord(10:12)=(/0.5_dp,0.0_dp,0.5_dp/) !Si
      UC%AtomCoord(13:15)=(/0.75_dp,0.25_dp,0.75_dp/) !int si
      UC%AtomCoord(16:18)=(/0.25_dp,0.75_dp,0.75_dp/) !int si
      UC%AtomCoord(19:21)=(/0.25_dp,0.25_dp,0.25_dp/) !int si
      UC%AtomCoord(22:24)=(/0.75_dp,0.75_dp,0.25_dp/) !int si

      UC%AtomSpecies=1
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE Silicon
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LiFePO4(UC)
   !Unit cell for LiFePO4 obtained from Phys. Rev. B 68, 165107 (2003) - originally from Chem. Mater. 13 1570 (2001)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,28) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=28
      UC%NMove=28
      !UC%BoxSize=(/10.06_dp,5.89_dp,4.65_dp,0._dp,0._dp,0._dp/)
      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      
      UC%AtomCoord(1:84)  =(/0.0_dp,0.0_dp,0.0_dp, & !cation, e.g., Li+
        0.5_dp,0.5_dp,0.5_dp,  & !Li
        0._dp,0.5_dp,0._dp, & !Li
        0.5_dp,0._dp,0.5_dp, & !Li
      
        7._dp/25._dp,0.25_dp,0.98_dp, & !Fe
        18._dp/25._dp,0.75_dp,0.2_dp, &
        0.78_dp,0.25_dp,13._dp/25._dp, &
        0.22_dp,0.75_dp,12._dp/25._dp, &
      
        0.09_dp,0.25_dp,0.42_dp, &
        0.91_dp,0.75_dp,.58_dp, &
        0.59_dp,.25_dp,2._dp/25._dp, &
        0.41_dp,0.75_dp,23._dp/25._dp, &
      
        0.09_dp,0.25_dp,0.75_dp, &
        9._dp/20._dp,0.25_dp,0.21_dp, &
        0.17_dp,1._dp/25._dp,.29_dp, &
        0.91_dp,0.75_dp,0.25_dp, &
        11._dp/20._dp,.75_dp,.79_dp, &
        .83_dp,24._dp/25._dp,.71_dp, &
        .59_dp,0.25_dp,0.75_dp, &
        .41_dp,.75_dp,.25_dp, &
        19._dp/20._dp,.25_dp,.29_dp, &
        1._dp/20._dp,.75_dp,.71_dp, &
        .67_dp,.46_dp,.21_dp, &
        .83_dp,.54_dp,.71_dp, &
        .33_dp,24._dp/25._dp,.79_dp, &
        .33_dp,.54_dp,.79_dp, &
        .17_dp,.46_dp,.29_dp, &
        .67_dp,1._dp/25._dp,.21_dp /)

      UC%AtomSpecies(1:4)=1 
      UC%AtomSpecies(5:8)=2
      UC%AtomSpecies(9:12)=3 
      UC%AtomSpecies(13:28)=4
      
      UC%AtomCharge(1:4)=1._dp
      UC%AtomCharge(5:8)=2._dp
      UC%AtomCharge(9:12)=5._dp
      UC%AtomCharge(13:28)=-2._dp
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE LiFePO4
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LiHalfFePO4(UC)
   !Unit cell for Li0.5FePO4 obtained from Phys. Rev. B 68, 165107 (2003) - originally from Chem. Mater. 13 1570 (2001)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,7) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=7
      UC%NMove=7
      !UC%BoxSize=(/9.96_dp,5.83_dp,4.7_dp,0._dp,0._dp,0._dp/)
      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      
      UC%AtomCoord(1:3)  =(/0.0_dp,0.0_dp,0.0_dp/) !cation, e.g., Li+
      UC%AtomCoord(4:6)  =(/0.28_dp,0.26_dp,0.98_dp/) !cation, e.g. Fe2+
      UC%AtomCoord(7:9)  =(/0.10_dp,0.24_dp,0.42_dp/) !cation, e.g., P5+
      UC%AtomCoord(10:12)=(/0.11_dp,0.23_dp,0.73_dp/) !anion, e.g., O2-
      UC%AtomCoord(13:15)=(/0.45_dp,0.24_dp,0.18_dp/) !anion, e.g., O2-
      UC%AtomCoord(16:18)=(/0.17_dp,0.04_dp,0.27_dp/) !anion, e.g., O2-
      UC%AtomCoord(19:21)=(/0.32_dp,0.55_dp,0.80_dp/) !anion, e.g., O2-

      UC%AtomSpecies(1)=1 
      UC%AtomSpecies(2)=2
      UC%AtomSpecies(3)=3 
      UC%AtomSpecies(4:7)=4
      
      UC%AtomCharge(1)=1._dp
      UC%AtomCharge(2)=2._dp
      UC%AtomCharge(3)=5._dp
      UC%AtomCharge(4:7)=-2._dp
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE LiHalfFePO4
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE FePO4(UC)
   !Unit cell for FePO4 obtained from Phys. Rev. B 68, 165107 (2003) - originally from Chem. Mater. 13 1570 (2001)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      CALL MakeSize(UC,5) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=5
      UC%NMove=5
      !UC%BoxSize=(/9.81_dp,5.79_dp,4.79_dp,0._dp,0._dp,0._dp/)
      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      
      UC%AtomCoord(1:3)  =(/0.27_dp,0.25_dp,0.95_dp/) !cation, e.g. Fe2+
      UC%AtomCoord(4:6)  =(/0.09_dp,0.25_dp,0.40_dp/) !cation, e.g., P5+
      UC%AtomCoord(7:9)=(/0.12_dp,0.25_dp,0.71_dp/) !anion, e.g., O2-
      UC%AtomCoord(10:12)=(/0.44_dp,0.25_dp,0.15_dp/) !anion, e.g., O2-
      UC%AtomCoord(13:15)=(/0.17_dp,0.04_dp,0.25_dp/) !anion, e.g., O2-

      UC%AtomSpecies(1)=1 
      UC%AtomSpecies(2)=2
      UC%AtomSpecies(3:5)=3 
      
      UC%AtomCharge(1)=2._dp
      UC%AtomCharge(2)=5._dp
      UC%AtomCharge(3:5)=-2._dp
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE FePO4
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Rutile(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !cation1 sc sites. cation2 face centered. anion center
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
!       CALL MakeSize(UC%Coord,5)
!       UC%NAtoms=5
!       UC%NSpeciesType=3
!       !UC%SpeciesDirectory=(/atindx(1),atindx(2),atindx(3)/)
!       UC%AtomCharge(1)=4 (1:4)=2._dp; UC%AtomCharge(5:12)=-1._dp
!       UC%AtomSpecies(1:4)=1; UC%AtomSpecies(5:12)=2
!       UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
!       UC%AtomCoord(1:3)  =(0.0_dp,0.0_dp,0.0_dp) !cation1
!       UC%AtomCoord(3:6)  =(0.5_dp,0.5_dp,0.0_dp) !cation2
!       UC%AtomCoord(7:9)  =(0.0_dp,0.5_dp,0.5_dp) !cation2
!       UC%AtomCoord(10:12)=(0.5_dp,0.0_dp,0.5_dp) !cation2
!       UC%AtomCoord(13:15)=(0.5_dp,0.5_dp,0.5_dp) ! anion
   END SUBROUTINE Rutile
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Perovskite(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 1
      !Examples:
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating ABO3-bulk crystal"

      CALL MakeSize(UC,5) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=5
      UC%NMove=5
      !UC%NSpeciesType=3
      !CALL MakeSize(UC%SpeciesDirectory,3)
      !UC%SpeciesDirectory=(/1,2,3/)

      UC%BoxSize=(/1._dp,1._dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.0_dp,0.0_dp/) !cation1(+2 charge)
      UC%AtomCoord(4:6)=(/0.5_dp,0.5_dp,0.5_dp/) !cation2(+4 charge)
      UC%AtomCoord(7:9)=(/0.5_dp,0.5_dp,0.0_dp/) !anion
      UC%AtomCoord(10:12)=(/0.5_dp,0.0_dp,0.5_dp/) !anion
      UC%AtomCoord(13:15)=(/0.0_dp,0.5_dp,0.5_dp/) !anion
      
      UC%AtomSpecies=(/1,2,3,3,3/)
      UC%AtomCharge=(/2.0_dp,4.0_dp,-2.0_dp,-2.0_dp,-2.0_dp/)
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/1._dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,1._dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,1._dp/)
   END SUBROUTINE Perovskite
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE graphite(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating graphite-bulk crystal"

      CALL MakeSize(UC,32) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=32
      UC%NMove(1:3)=32
      !UC%NSpeciesType=1
      !CALL MakeSize(UC%SpeciesDirectory,1)
      !UC%SpeciesDirectory=1

      UC%BoxSize=(/9.26400000_dp,5.34857289_dp,8.41040650_dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord=(/ .00000000_dp, .00000000_dp, .00000000_dp, &
         .77200000_dp,4.01142967_dp, .00000000_dp, &
         7.72000000_dp, .00000000_dp, .00000000_dp, &
         .77200000_dp,1.33714322_dp, .00000000_dp, &
         .00000000_dp,2.67428645_dp, .00000000_dp, &
         6.94800000_dp,4.01142967_dp, .00000000_dp, &
         2.31600000_dp,4.01142967_dp, .00000000_dp, &
         6.94800000_dp,1.33714322_dp, .00000000_dp, &
         2.31600000_dp,1.33714322_dp, .00000000_dp, &
         7.72000000_dp,2.67428645_dp, .00000000_dp, &
         3.08800000_dp, .00000000_dp, .00000000_dp, &
         3.08800000_dp,2.67428645_dp, .00000000_dp, &
         5.40400000_dp,4.01142967_dp, .00000000_dp, &
         5.40400000_dp,1.33714322_dp, .00000000_dp, &
         .00000000_dp, .00000000_dp,-4.20520325_dp, &
         .77200000_dp,4.01142967_dp,-4.20520325_dp, &
         7.72000000_dp, .00000000_dp,-4.20520325_dp, &
         .77200000_dp,1.33714322_dp,-4.20520325_dp, &
         4.63200000_dp, .00000000_dp, .00000000_dp, &
         6.94800000_dp,4.01142967_dp,-4.20520325_dp, &
         2.31600000_dp,4.01142967_dp,-4.20520325_dp, &
         .00000000_dp,2.67428645_dp,-4.20520325_dp, &
         6.94800000_dp,1.33714322_dp,-4.20520325_dp, &
         2.31600000_dp,1.33714322_dp,-4.20520325_dp, &
         7.72000000_dp,2.67428645_dp,-4.20520325_dp, &
         3.08800000_dp, .00000000_dp,-4.20520325_dp, &
         4.63200000_dp,2.67428645_dp, .00000000_dp, &
         3.08800000_dp,2.67428645_dp,-4.20520325_dp, &
         5.40400000_dp,4.01142967_dp,-4.20520325_dp, &
         5.40400000_dp,1.33714322_dp,-4.20520325_dp, &
         4.63200000_dp, .00000000_dp,-4.20520325_dp, &
         4.63200000_dp,2.67428645_dp,-4.20520325_dp/)
      
      UC%AtomSpecies=1
      UC%AtomCharge=(/1.0_dp,-1.0_dp/)
      UC%AtomIsMoving=.TRUE.
      Basis(1:3)=(/9.26400000_dp,0.0_dp,0.0_dp/)
      Basis(4:6)=(/0.0_dp,5.34857289_dp,0.0_dp/)
      Basis(7:9)=(/0.0_dp,0.0_dp,8.41040650_dp/)
   END SUBROUTINE graphite
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE alumina(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 1
      !Examples: 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating alumina-bulk crystal"
      
      CALL MakeSize(UC,1) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=30
      UC%NMove(1:3)=30
      !UC%NSpeciesType=2
      !CALL MakeSize(UC%SpeciesDirectory,2)
      !UC%SpeciesDirectory=2 !replace with atomic # later

      UC%BoxSize=(/0.70710678_dp,0.70710678_dp,1._dp,0._dp,0._dp,0._dp/)
      
      UC%AtomCoord(1:3)=(/0._dp,0._dp,0.355_dp/) !Al
      UC%AtomCoord(4:6)=(/0._dp,0._dp,0.645_dp/)
      UC%AtomCoord(7:9)=(/0._dp,0._dp,0.855_dp/)
      UC%AtomCoord(10:12)=(/0._dp,0._dp,0.145_dp/)
      UC%AtomCoord(13:15)=(/0.33333333_dp,0.66666666_dp,0.02167_dp/)
      UC%AtomCoord(16:18)=(/0.33333333_dp,0.66666666_dp,0.52170_dp/)
      UC%AtomCoord(19:21)=(/0.33333333_dp,0.66666666_dp,0.31167_dp/)
      UC%AtomCoord(22:24)=(/0.33333333_dp,0.66666666_dp,0.81167_dp/)
      UC%AtomCoord(25:27)=(/0.66666666_dp,0.33333333_dp,0.68833_dp/)
      UC%AtomCoord(28:30)=(/0.66666666_dp,0.33333333_dp,0.18833_dp/)
      UC%AtomCoord(31:33)=(/0.66666666_dp,0.33333333_dp,0.97833_dp/)
      UC%AtomCoord(34:36)=(/0.66666666_dp,0.33333333_dp,0.47883_dp/)
      
      UC%AtomCoord(37:39)=(/0.3030_dp,0._dp,0.25_dp/) !O
      UC%AtomCoord(40:42)=(/0.6970_dp,0._dp,0.75_dp/)
      UC%AtomCoord(43:45)=(/0._dp,0.303_dp,0.25_dp/)
      UC%AtomCoord(46:48)=(/0.6970_dp,0.697_dp,0.25_dp/)
      UC%AtomCoord(49:51)=(/0._dp,0.697_dp,0.75_dp/)
      UC%AtomCoord(52:54)=(/0.3030_dp,0.3030_dp,0.75_dp/)
      UC%AtomCoord(55:57)=(/0.6360_dp,0.666_dp,0.9166_dp/)
      UC%AtomCoord(58:60)=(/0.3330_dp,0.96970_dp,0.9166_dp/)
      UC%AtomCoord(61:63)=(/0.03033_dp,0.3636_dp,0.9166_dp/)
      UC%AtomCoord(64:66)=(/0.3330_dp,0.3636_dp,0.4166_dp/)
      UC%AtomCoord(67:69)=(/0.03033_dp,0.666_dp,0.4166_dp/)
      UC%AtomCoord(70:72)=(/0.63630_dp,0.9696_dp,0.4166_dp/)
      UC%AtomCoord(73:75)=(/0.9696_dp,0.333_dp,0.5833_dp/)
      UC%AtomCoord(76:78)=(/0.66_dp,0.6363_dp,0.5833_dp/)
      UC%AtomCoord(79:81)=(/0.3636_dp,0.03033_dp,0.5833_dp/)
      UC%AtomCoord(82:84)=(/0.66_dp,0.03033_dp,0.0833_dp/)
      UC%AtomCoord(85:87)=(/0.3636_dp,0.33_dp,0.0833_dp/)
      UC%AtomCoord(88:90)=(/0.96967_dp,0.636_dp,0.0833_dp/)
      
      UC%AtomSpecies=2
      UC%AtomCharge(1:12)=3._dp
      UC%AtomCharge(13:30)=-2._dp
      UC%AtomIsMoving=.TRUE.
   END SUBROUTINE alumina
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE diamond(UC)
      !Add NAtoms,NMove,NSpeciesType,AtomCharge,AtomCoord,AtomIsMoving,AtomSpecies,BoxSize,Species
      !Lattice constant = 1
      !# atom species present = 1
      !Examples: 
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: UC
      
      WRITE(6,*) "... Creating diamond-bulk crystal"
      
      CALL MakeSize(UC,1) !Elements not allocated @ this point - SpeciesDirectory
      UC%NAtoms=1
      UC%NMove(1:3)=1
      !UC%NSpeciesType=1
      !CALL MakeSize(UC%SpeciesDirectory,1)
      !UC%SpeciesDirectory=1 !replace with atomic # later

      UC%BoxSize=(/0.70710678_dp,0.70710678_dp,1._dp,0._dp,0._dp,0._dp/)
      UC%AtomCoord(1:3)=(/0.0_dp,0.35355339_dp,0.75000000_dp/)
      UC%AtomCoord(4:6)=(/0.0_dp,0.0_dp,0.5_dp/)
      UC%AtomCoord(7:9)=(/0.35355339_dp,0.0_dp,0.25_dp/)
      UC%AtomCoord(10:12)=(/0.35355339_dp,0.35355339_dp,0.0_dp/)
      
      UC%AtomSpecies=1
      UC%AtomCharge=0._dp
      UC%AtomIsMoving=.TRUE.
   END SUBROUTINE diamond
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddGaussianNoise(AL,disp,mask)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(sp), DIMENSION(:), ALLOCATABLE :: v
      LOGICAL, DIMENSION(:), POINTER, OPTIONAL :: mask !mask is of size AL%NAtoms
      LOGICAL, DIMENSION(:), POINTER :: mask1 !mask is of size AL%NAtoms
      REAL(dp) :: disp
      INTEGER :: NAtoms,i
      
      NAtoms=AL%NAtoms
      CALL MakeSize(mask1,3*AL%NAtoms)
      IF (PRESENT(mask)) THEN
         mask1=mask
      ELSE
         mask1=.TRUE.
      END IF
      
      ALLOCATE(v(3*NAtoms)) !note v is in single precision
      IF (3*NAtoms>600000) THEN
      !intel fortran compiler gives problems when the size is greater than 7800000
         DO i=1,3*NAtoms
            CALL gasdev_s(v(i))
         END DO
      ELSE
         CALL gasdev_v(v)
      END IF
      WHERE (ABS(v)>1._dp) v=1._dp
      v=v*disp
      WHERE (mask1) AL%AtomCoord=AL%AtomCoord+v
      DEALLOCATE(v)
      DEALLOCATE(mask1)
   END SUBROUTINE AddGaussianNoise
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RadialDistributionFn(AL,rcut,nbins,atom,species)
      !finds the radial distribution function for all atoms or only atom selected
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER, OPTIONAL :: atom,species
      REAL(dp) :: rcut,binsize,density,radius
      INTEGER :: nbins,ibin,iatom,ispecies
      REAL(dp), DIMENSION(:), POINTER :: rdf !this contains the information

      ALLOCATE(rdf(nbins))
      rdf=0
      binsize=rcut/REAL(nbins)
      WRITE(6,*) "Binsize for radial distribution function:",binsize

      IF (PRESENT(species)) THEN
         ispecies=species
      ELSE
         ispecies=-999
      END IF
      
      IF (PRESENT(atom)) THEN
         iatom=atom
         CALL UpdateRDF(iatom)
      ELSE
         DO iatom=1,AL%NAtoms
            IF (AL%AtomSpecies(iatom)==ispecies .OR. ispecies==-999) CALL UpdateRDF(iatom)
         END DO
      END IF

      density=REAL(AL%NAtoms)/PRODUCT(AL%BoxSize(1:3)-AL%BoxSize(4:6))
      rdf=rdf/density/REAL(AL%NAtoms)
      DO ibin=1,nbins
         radius=binsize*(REAL(ibin)-0.5)
         WRITE(6,*) radius,rdf(ibin)/4._dp/PI/radius/radius/binsize
      END DO

      DEALLOCATE(rdf)

      CONTAINS
      SUBROUTINE UpdateRDF(iatom)
         INTEGER :: jatom,ibin,iatom
         REAL(dp) :: dr(3),drmag

         DO jatom=1,AL%NAtoms
            IF (iatom==jatom) CYCLE
            dr=PBCdistance1(AL,iatom,jatom)
            drmag=SQRT(DOT_PRODUCT(dr,dr))
            ibin=CEILING(drmag/binsize)
            IF (ibin<=nbins) rdf(ibin)=rdf(ibin)+1._dp
         END DO
      END SUBROUTINE UpdateRDF
   END SUBROUTINE RadialDistributionFn
   !  Recognize the crystal structure from the RadialDistributionFn
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AngleDistributionFn
   END SUBROUTINE AngleDistributionFn
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE FilterAtoms(AL,filtermode,FilteredAtomCoord,FilteredAtomSpecies, &
      CrystalType,LatticeConst,NFilteredAtom)
   !mask contains information regarding whether the atom has to be retained
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,UC
      TYPE(LinkedListContainer), POINTER :: LL
      INTEGER :: filtermode,NFilteredAtom
      REAL(dp), DIMENSION(3*MaxNFilteredAtoms) :: FilteredAtomCoord
      INTEGER, DIMENSION(MaxNFilteredAtoms) :: FilteredAtomSpecies
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomCoordu
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      REAL(dp) :: tol,tolsq,drmag,dr(3),LatticeConst(3),coord(3),coord1(3),coordu(3),v(3)
      REAL(dp) :: BoxSize(3),xc(3)
      CHARACTER(len=10), INTENT(IN) :: CrystalType
      INTEGER :: iatom,ispecies,NAtoms,iunitatom,NUnitCells(3),a,b,c
      LOGICAL :: Found
      
      SELECT CASE (filtermode)
      CASE(1) !sharad maheshwari's filter for radiation damage
         NAtoms=AL%NAtoms
         NFilteredAtom=0
         AtomCoord=>AL%AtomCoord
         AtomSpecies=>AL%AtomSpecies
         !obtain the unit cell from crystal.f90
         CALL UnitCell(UC,CrystalType)
         tol=0.05 !dimensionless quantity
         tolsq=tol*tol
         !Step 1. find the interstitial atoms
         DO iatom=1,NAtoms
            coord=AtomCoord(3*iatom-2:3*iatom)
            coord1=MOD(coord,LatticeConst)/LatticeConst !wrap the coordinate & make dimensionless
            ispecies=AtomSpecies(iatom)
            !now compare positions with unit cell atoms
            AtomCoordu=>UC%AtomCoord
            Found=.FALSE.
            DO iunitatom=1,UC%NAtoms
               coordu=AtomCoordu(3*iunitatom-2:3*iunitatom)
               dr=coord1-coordu
               Found=dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3) < tolsq
               IF (Found) EXIT
            END DO
            IF (.NOT. Found) THEN
               NFilteredAtom=NFilteredAtom+1
               FilteredAtomCoord(3*NFilteredAtom-2:3*NFilteredAtom)=coord
               FilteredAtomSpecies(NFilteredAtom)=ispecies
            END IF
         END DO
         !Step 2. find the vacancies - assumes that the linked list is correctly setup
         BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
         xc=AL%BoxSize(4:6)
         LL=>AL%LL
         NUnitCells=FLOOR(BoxSize/LatticeConst)
         DO a=1,NUnitCells(1)
            DO b=1,NUnitCells(2)
               DO c=1,NUnitCells(3)
                  v=a*Basis(1:3)+b*Basis(4:6)+c*Basis(7:9)
                  v=v*LatticeConst !from dimensionless get units of Angstrom
                  DO iunitatom=1,UC%NAtoms
                     coordu=AtomCoordu(3*iunitatom-2:3*iunitatom)+v !translate unit cell to correct position
                     Found=LocateAtomInCrystal()
                     IF (Found) EXIT
                  END DO
                  IF (.NOT. Found) THEN
                     NFilteredAtom=NFilteredAtom+1
                     FilteredAtomCoord(3*NFilteredAtom-2:3*NFilteredAtom)=coord
                     FilteredAtomSpecies(NFilteredAtom)=ispecies
                  END IF
               END DO
            END DO
         END DO
      CASE DEFAULT
         WRITE(6,*) "$Err>> Filtermode has not been implemented yet in FilterAtomsForPrinting"
         STOP
      END SELECT
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      FUNCTION LocateAtomInCrystal()
         IMPLICIT NONE
         LOGICAL :: LocateAtomInCrystal
         
         LocateAtomInCrystal=.FALSE.
      END FUNCTION LocateAtomInCrystal
   END SUBROUTINE FilterAtoms
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE Crystal
