!xxxxxxxxxxxxxxxxxxxxxxxxxxxx  
SUBROUTINE SetupSpeciesType(isite)
!setup the species index
   IMPLICIT NONE
   INTEGER :: isite
   CHARACTER(len=3) :: Symbol
   
   Symbol=SpeciesList%AtomicSymbol(SpeciesDirectory_global(AL%AtomSpecies(isite)))
   IF(Symbol=="Fe") THEN
      SpeciesType(isite)=1
   ELSEIF (Symbol=="Pt") THEN
      SpeciesType(isite)=2 
   ELSE  
      SpeciesType(isite)=0
   END IF
END SUBROUTINE SetupSpeciesType
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE NearestNeighborSiteList(isite,isprint,coreradius,center)
!find list of index which are neighbors
   IMPLICIT NONE
   INTEGER :: isite,j,jsite,n,ra,r0,r1,r2,s1,s2,m,p
   REAL(dp) ::drmagij,coord(3)
   REAL(dp), OPTIONAL :: coreradius,center(3)
   LOGICAL, OPTIONAL :: isprint
   LOGICAL :: isprint1,IsSiteCore

   n=0 !number of sites
   
   ra=(isite-1)*MaxAtomPerAtom !start index minus 1
   r0=VLListRange(isite) !range for ith atom
   r1=ra+1
   r2=ra+r0
   s1=AtomSpecies(isite)
   coord=AtomCoord(3*isite-2:3*isite)

   isprint1=.FALSE.
   IF (PRESENT(isprint)) isprint1=isprint

   DO j=r1,r2 !this will take care of the pair and many-body terms in one go
      jsite=VLList(j)
      IF (isprint1) THEN
         WRITE(6,*) jsite,VLdrmag(j)
      END IF
      s2=AtomSpecies(jsite)
      drmagij=VLdrmag(j)   ! spacing between atoms i and j
      IF(drmagij<=3.55_dp) THEN
            n=n+1
            NeighborSiteList((isite-1)*12+n)=jsite
      END IF 
!      IF (isite==347805 .AND. drmagij<5.1_dp) WRITE(6,*) jsite, s2,   drmagij 
   END DO

   NumberNeighborSites(isite)=n

!   IF (PRESENT(coreradius) .AND. PRESENT(center)) THEN
!      IsSiteCore=(SQRT(DOT_PRODUCT(coord-center,coord-center))<coreradius)
!      IF (NumberNeighborSites(isite)/=12 .AND. IsSiteCore) THEN
!         DO j=r1,r2
!            jsite=VLList(j)
!            WRITE(6,*) jsite,VLdrmag(j)
!         END DO
!         STOP
!      END IF
!   END IF

   m=0
   n=0
   p=0
   DO j=r1,r2 !this will take care of the pair and many-body terms in one go
      jsite=VLList(j)
      s2=AtomSpecies(jsite)
      drmagij=VLdrmag(j)   ! spacing between atoms i and j
      IF(s2>0 .AND. drmagij<3.55_dp) m=m+1
      IF(s2>0 .AND. drmagij<5.35_dp) THEN
        n=n+1
        IF(s2==2) p=p+1
      END IF
   END DO

   FirstNearestNeighborSites(isite)=m
   ThirdNearestNeighborSites(isite)=n
   Local_composition(isite)=0._dp
   IF (n>0) Local_composition(isite)=DBLE(p)/DBLE(n)

END SUBROUTINE NearestNeighborSiteList
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE InitializeRates()
   IMPLICIT NONE
   REAL(dp), PARAMETER :: kboltzmann=8.62e-5
   INTEGER :: NAtm
   
   RateDiffusion=0.d0
   RateDissolution=0.d0
   DO NAtm=0,11
      RateDiffusion(NAtm)=prefac*(EXP(-(NAtm*xi+e0)/(kboltzmann*Temperature1)))
      RateDissolution(NAtm)=pnu_E*(exp(-((NAtm*xi)-phi_E)/(kboltzmann*Temperature1)))
   END DO

   !RateDiffusion(0)=0._dp !cant allow this to move

   WRITE(6,*) "NAtm   RateDiffusion RateDissolution"
   DO NAtm=0,12
      WRITE(UNIT=6,FMT='(i3,2es20.8)') NAtm, &
        RateDiffusion(NAtm),RateDissolution(NAtm)
   END DO
END SUBROUTINE InitializeRates
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx 
SUBROUTINE UpdateRateArray(isite)
!update the rate at site i
   IMPLICIT NONE
   INTEGER :: isite,NAtm,k,jsite,n
   
   Rate((isite-1)*13+1:isite*13)=0._dp !initialize
   IF (speciestype(isite)==0) RETURN !site is already vacant
   
   !Find number of atoms present in the neighbor site
   NAtm=0
   n=NumberNeighborSites(isite)
   DO k=1,n
      jsite=NeighborSiteList((isite-1)*12+k)
      !write(6,*) "Site #:",k,jsite
      !call flush(6)
      IF (jsite>0) THEN
         IF(speciestype(jsite)>0) NAtm=NAtm+1
      END IF
   END DO
   
   DO k=1,n  !n is fixed for a particular site
      jsite=NeighborSiteList((isite-1)*12+k)
      IF (jsite>0) THEN !neighbor site is present
         IF (speciestype(jsite)==0) THEN !vacancy present
            Rate((isite-1)*13+k)=RateDiffusion(NAtm)
         END IF
      END IF
   END DO
   
   Coordination(isite)=NAtm

   IF (speciestype(isite)==1) Rate(isite*13)=RateDissolution(NAtm)
   
END SUBROUTINE UpdateRateArray
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE UpdatesAffectedSites(isite,jmove)
!updates the species at the sites participating in the move
!as well as the list of sites whose rates are affected by the move
   IMPLICIT NONE
   INTEGER :: isite,jmove,jsite,ksite,k,m
   
   IF(jmove==13) THEN !dissolution
      speciestype(isite)=0
   ELSE !diffusion event
      jsite=NeighborSiteList((isite-1)*12+jmove)
      IF (jsite==0) THEN
         WRITE(6,*) "Site not found!"
         STOP
      END IF
      IF (speciestype(jsite)==0) THEN
         speciestype(jsite)=speciestype(isite)
         speciestype(isite)=0
      ELSE
         WRITE(6,*) "Cannot hop to neighbor site which is occupied"
         WRITE(6,*) "isite,sp,coord:",isite,speciestype(isite),AtomCoord(3*isite-2:3*isite)
         WRITE(6,*) "jsite:",jsite
         WRITE(6,*) "jmove:",jmove
         CALL PrintListofNNSites(isite)
         STOP
      END IF
   END IF
   
   NAffectedSites=1
   AffectedSiteList(1)=isite
   DO k=1,NumberNeighborSites(isite)
      NAffectedSites=NAffectedSites+1
      AffectedSiteList(NAffectedSites)=NeighborSiteList((isite-1)*12+k)
   END DO
   IF (VERBOSE) WRITE(6,*) "# atoms added as affected sites:",NAffectedSites
   
   IF(jmove<13) THEN !add neighbors of jatom
      DO k=1,NumberNeighborSites(jsite)
        ksite=NeighborSiteList((jsite-1)*12+k)
        IF (.NOT. ANY(AffectedSiteList(1:NAffectedSites)==ksite)) THEN !add this atom
           NAffectedSites=NAffectedSites+1
           AffectedSiteList(NAffectedSites)=ksite
        END IF
      END DO
   END IF
!   WRITE(*,*)AffectedSiteList

   !update rates
   NUpdatedPositions=0
   DO k=1,NAffectedSites
      ksite=AffectedSiteList(k)
      CALL UpdateRateArray(ksite)
      DO m=1,13
         NUpdatedPositions=NUpdatedPositions+1
         UpdatedPositions(NUpdatedPositions)=(ksite-1)*13+m
      END DO
   END DO
   IF (NUpdatedPositions>24*13) THEN
      WRITE(6,*) "Unexpected number of rates that were updated"
      STOP
   END IF
   IF (VERBOSE) WRITE(6,*) "# rates to be updated in BT:",NUpdatedPositions
   CALL BinaryTreeLocalUpdate(UpdatedPositions(1:NUpdatedPositions), &
      NUpdatedPositions)
END SUBROUTINE UpdatesAffectedSites
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE PrintFrame()
   IMPLICIT NONE
   !LOGICAL :: Satisfied
   INTEGER, SAVE :: iframe=0
   REAL(dp) :: EOD
   REAL(dp), SAVE ::TargetEOD=.20

   !Satisfied=.FALSE.
   !IF (.NOT. Satisfied) RETURN
   EOD=DBLE(NDissolved)/DBLE(NElectroactive_initial)
   IF (EOD<TargetEOD) RETURN

   iframe=iframe+1
   !IF (time>5.e9 .aND. time<1.e10) THEN
   IF (EOD>TargetEOD) THEN
      WRITE(6,*) "Printing xyz for nanowire ..."
      CALL GetPrintableAtoms(.FALSE.)
      OPEN(UNIT=302,FILE="NW."//TRIM(INT2CHAR(iframe))//".xyz") !for analysis
      CALL PrintFrame1(302,NPrintableAtoms,PrintableAtom,0)
      CLOSE(302)
   END IF
   CALL GetPrintableAtoms(.TRUE.)
   CALL CharacterizeSurface() !identifies (111), (100) and other surface atoms
     !among atoms which are printable
   CALL CharacterizeFacets()
   CALL PrintCoordination()
!   CALL CharacterizeSubsurface()
   !IF (time>5.e9 .aND. time<1.e10) THEN
!   IF (EOD>TargetEOD) THEN
!      CALL PrintSurface(111,iframe) !
!      CALL PrintSurface(100,iframe)
!      CALL PrintSurface(0,iframe) !others
!!      CALL PrintSubSurface(111,iframe)
!!      CALL PrintSubSurface(100,iframe)
!   END IF
   IF (EOD>TargetEOD) THEN
      TargetEOD=TargetEOD+.20_dp
      IF (TargetEOD>=.20_dp) TargetEOD=.20
   END IF


   WRITE(6,*) "Printing xyz for nanowire shell...",iframe," EOD:",EOD," time:",time
!   OPEN(UNIT=302,FILE="NWShell."//TRIM(INT2CHAR(iframe))//".xyz") !for plotting
!   CALL PrintFrame1(302,NPrintableAtoms,PrintableAtom,0)
!   CLOSE(302)

   WRITE(548,*) kmcx,kmcd,time,fracdiffusion
IF (EOD>.20) STOP
!IF (iframe>250000000000000) THEN
!WRITE(6,*) "iframe more than 250"
!STOP
!END IF
END SUBROUTINE PrintFrame
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE PrintCoordination()
   IMPLICIT NONE
   INTEGER :: i,n1(0:12),n2(0:12),x

   n1=0; n2=0
   DO isite=1,NSites
      IF (speciestype(isite)==1) THEN
         x=coordination(isite)
         n1(x)=n1(x)+1
      END IF
      IF (speciestype(isite)==2) THEN
         x=coordination(isite)
         n2(x)=n2(x)+1
      END IF
   END DO
   WRITE(UNIT=543,FMT='(14ES12.5)') kmcx,kmcd,time,REAL(n1)/REAL(NSites)
   WRITE(UNIT=544,FMT='(14ES12.5)') kmcx,kmcd,time,REAL(n2)/REAL(NSites)
   WRITE(UNIT=545,FMT='(14ES12.5)') kmcx,kmcd,time,REAL(n1+n2)/REAL(NSites)
   WRITE(549,*) kmcx,kmcd,time,diffncoordination(0:12)/SUM(diffncoordination)
   WRITE(550,*) kmcx,kmcd,time,dislncoordination(0:12)/SUM(dislncoordination)

END SUBROUTINE PrintCoordination
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE PrintFrame1(iunit,num,mask,printtype)
   IMPLICIT NONE
   INTEGER :: iunit,num,printtype
   LOGICAL, DIMENSION(:), POINTER :: mask

   WRITE(UNIT=iunit,FMT='(4i10," NAtomsInFrame, NSites, NElectroactive, Nonactive")')  &
     num,NSites, & !NPrintableAtoms,NSites, &
     NElectroactive-NDissolved,Nonactive   ! changed here
   WRITE(UNIT=iunit,FMT='(4es15.5)') AL%BoxSize(1:3),time
   DO isite=1,NSites
      IF (mask(isite)) & !IF (PrintableAtom(isite)) &
         WRITE(UNIT=iunit,FMT='(a4,3f20.8,i8,"   ",i8,i8)') &
          PrintSpeciesType(isite,printtype),AtomCoord(3*isite-2:3*isite), &
            isite,FirstNearestNeighborSites(isite),ThirdNearestNeighborSites(isite)
   END DO
END SUBROUTINE PrintFrame1
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE PrintSurface(isurf,iframe)
   IMPLICIT NONE
   INTEGER :: isurf,iframe,n,i
   INTEGER, DIMENSION(:), POINTER :: arr

   SELECT CASE(isurf)
   CASE(111) !111
      OPEN(UNIT=302,FILE="NWSurf111."//TRIM(INT2CHAR(iframe))//".xyz") 
      arr=>IsSurf111atom
   CASE(100) !100
      OPEN(UNIT=302,FILE="NWSurf100."//TRIM(INT2CHAR(iframe))//".xyz") 
      arr=>IsSurf100atom
   CASE(0) !others
      OPEN(UNIT=302,FILE="NWSurfOthers."//TRIM(INT2CHAR(iframe))//".xyz") 
      arr=>IsSurfOtherAtom
   END SELECT
   
   mask=.FALSE.
   n=0
   DO i=1,NSites
      IF (arr(i)>0 .AND. coordination(i)<11) THEN
         n=n+1
         mask(i)=.TRUE.
      END IF
   END DO
   IF (isurf==111) WRITE(6,*) "# sites with 111:",n
   IF (isurf==100) WRITE(6,*) "# sites with 100:",n
   !CALL PrintFrame1(302,n,mask,isurf)
   CALL PrintFrame1(302,n,mask,0)
   CLOSE(302)
END SUBROUTINE PrintSurface
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE PrintSubSurface(isurf,iframe)
   IMPLICIT NONE
   INTEGER :: isurf,iframe,n,i
   INTEGER, DIMENSION(:), POINTER :: arr

   SELECT CASE(isurf)
   CASE(111) !111
      OPEN(UNIT=302,FILE="NWSubSurf111."//TRIM(INT2CHAR(iframe))//".xyz") 
      arr=>IsSubSurf111atom
   CASE(100) !100
      OPEN(UNIT=302,FILE="NWSubSurf100."//TRIM(INT2CHAR(iframe))//".xyz") 
      arr=>IsSubSurf100atom
   END SELECT
   
   mask=.FALSE.
   n=0
   DO i=1,NSites
      IF (arr(i)>0) THEN
         n=n+1
         mask(i)=.TRUE.
      END IF
   END DO
   IF (isurf==111) WRITE(6,*) "# subsurf sites with 111:",n
   IF (isurf==100) WRITE(6,*) "# subsurf sites with 100:",n
   CALL PrintFrame1(302,n,mask,0)
   CLOSE(302)
END SUBROUTINE PrintSubSurface
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
FUNCTION PrintSpeciesType(isite,printtype)
   IMPLICIT NONE
   CHARACTER(len=4) :: PrintSpeciesType
   INTEGER :: isite,printtype

   IF (printtype==0) THEN
      SELECT CASE (speciestype(isite))
      CASE(1); PrintSpeciesType="Fe  "
      CASE(2); PrintSpeciesType="Pt  "
      CASE DEFAULT
      WRITE(6,*) "Unexpected vacancy encountered"
      STOP
      END SELECT
   ELSEIF (printtype==111) THEN
      PrintSpeciesType="111"
   ELSEIF (printtype==100) THEN
      PrintSpeciesType="100"
   END IF
END FUNCTION PrintSpeciesType
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE PrintListofNNSites(isite)
   IMPLICIT NONE
   INTEGER :: isite,i,jsite

   DO i=1,NumberNeighborSites(isite)
      jsite=NeighborSiteList((isite-1)*12+i)
      WRITE(UNIT=6,FMT='("NN ",i2,":",i8," sp: ",i2," coord:",3f10.3," Rate:",3f10.3)') &
        i,jsite,speciestype(jsite),AtomCoord(3*jsite-2:3*jsite), &
        Rate((isite-1)*13+i)
   END DO
   DO i=NumberNeighborSites(isite)+1,12
      WRITE(UNIT=6,FMT='("NN ",i2," Rate:",f10.3)') i,Rate((isite-1)*13+i)
   END DO
   WRITE(UNIT=6,FMT='("NN ",i2," DissolRate:",f10.3)') 13,Rate((isite-1)*13+13)
END SUBROUTINE PrintListofNNSites
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE GetPrintableAtoms(PrintShell)
   IMPLICIT NONE
   LOGICAL :: PrintShell,UnderCoordNeigh
   INTEGER :: isite,n

   PrintableAtom=.FALSE.
   IF (PrintShell) THEN
      NPrintableAtoms=0
      DO isite=1,NSites
         IF (speciestype(isite)>0) THEN
            IF (Coordination(isite)<12) THEN
               PrintableAtom(isite)=.TRUE.
            ELSE
               CALL NeighborIsUndercoordinated(isite,UnderCoordNeigh)
               PrintableAtom(isite)=UnderCoordNeigh
            END IF
            IF (PrintableAtom(isite)) NPrintableAtoms=NPrintableAtoms+1
         END IF
      END DO
   ELSE
      NPrintableAtoms=NParticle-NDissolved !total number of atoms present
!n=0  ! changed here
      DO isite=1,NSites
         PrintableAtom(isite)=speciestype(isite)>0
!if (PrintableAtom(isite)) n=n+1   ! changed here
      END DO
!write(6,*) n   ! changed here
   END IF
END SUBROUTINE GetPrintableAtoms
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE NeighborIsUndercoordinated(isite,result)
   IMPLICIT NONE
   INTEGER :: isite,k,ksite
   LOGICAL :: result

   result=.FALSE.
   DO k=1,NumberNeighborSites(isite)
      ksite=NeighborSiteList((isite-1)*12+k)
      IF (speciestype(ksite)>0) THEN
         result=result .OR. Coordination(ksite)<12
      END IF
   END DO
END SUBROUTINE NeighborIsUndercoordinated
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE CharacterizeSurface()
   IMPLICIT NONE
   INTEGER :: isite

   IsSurf111atom=0 !all sites bulk/surface initialized
   IsSurf100atom=0
   IsSurfOtherAtom=0
   DO isite=1,NSites
      SurfaceInfo(isite)="                 "
      IF (PrintableAtom(isite)) THEN !lets look at its neighbor
         IsSurf111atom(isite)=CharacterizeSurface111(isite) !check if this is 111 surface
         IsSurf100atom(isite)=CharacterizeSurface100(isite) !check if this is 100 surface
         !IF (IsSurf111atom(isite)>0 .AND. IsSurf100atom(isite)>0) THEN
         !   WRITE(6,*) "Atom found to be both 111 and 100 surface atoms!"
         !END IF
         IF (IsSurf111atom(isite)>0 .AND. IsSurf100atom(isite)>0) THEN
            IsSurf111atom(isite)=0
            IsSurf100atom(isite)=0
         END IF
         IF (IsSurf111atom(isite)==0 .AND. IsSurf100atom(isite)==0) THEN
            IsSurfOtherAtom(isite)=1 !some other type of surface
            SurfaceInfo(isite)="Other"
         END IF
      END IF
   END DO
END SUBROUTINE CharacterizeSurface
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
FUNCTION CharacterizeSurface111(isite)
!NOTE: if atom lies at intersection of two 111 planes then only one will be assigned
!0 implies does not lie on 111 plane
   IMPLICIT NONE
   INTEGER :: isite,iplane,jsite,CharacterizeSurface111,ivec,n
   REAL :: vec(18),c0(3),cs(3)
   REAL, DIMENSION(3), PARAMETER :: L=(/244.806000,244.806000,734.418000/)

   DO iplane=1,4
     CALL GetFacetCoordinates(vec,111,iplane,18)

      !search coordinate
      c0=AtomCoord(3*isite-2:3*isite) !coordinate of site
      n=0
      !write(6,*) c0,"<-"
      DO ivec=1,6
         cs=c0+vec(3*ivec-2:3*ivec)
         !write(6,*) cs,".."
         cs=cs-L*FLOOR(cs/L) !periodic boundary condition
         !write(6,*) cs,"..pbc"
         jsite=LocateNeighSite(isite,cs)
         !write(6,*) "jsite",jsite
         !stop
         IF (jsite>0) THEN !jsite=0 implies site not found
            !write(*,*) jsite,":jsite"
            IF (speciestype(jsite)>0) THEN !actual atom present
               IF (PrintableAtom(jsite)) THEN !site lying in plane is also surface atom
                  n=n+1
                  !write(*,*) n
               END IF
            END IF
         END IF
      END DO
      
      !finalize
      IF (n>5) THEN !at least 3 NN should be surface atoms in 111 plane
         CharacterizeSurface111=iplane
         !stop
         RETURN
      END IF
   END DO
   !stop
   CharacterizeSurface111=0
END FUNCTION CharacterizeSurface111
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
FUNCTION CharacterizeSurface100(isite)
!NOTE: if atom lies at intersection of two 100 planes then only one will be assigned
!0 implies does not lie on 100 plane
   IMPLICIT NONE
   INTEGER :: isite,iplane,jsite,CharacterizeSurface100,ivec,n
   REAL :: z,a,vec(12),c0(3),cs(3)
   REAL, DIMENSION(3), PARAMETER :: L=(/244.806000,244.806000,734.418000/)

   DO iplane=1,3
      SELECT CASE(iplane) !select plane
      CASE(1); vec=(/-a,z,-a, a,z,a, -a,z,a, a,z,-a/)
      CASE(2); vec=(/-a,a,z, a,-a,z, -a,-a,z, a,a,z/)
      CASE(3); vec=(/z,-a,a, z,a,-a, z,-a,-a, z,a,a/)
      END SELECT
      CALL GetFacetCoordinates(vec,100,iplane,12)

      !search coordinate
      c0=AtomCoord(3*isite-2:3*isite) !coordinate of site
      n=0
      !write(6,*) c0,"<-"
      DO ivec=1,4
         cs=c0+vec(3*ivec-2:3*ivec)
         !write(6,*) cs,".."
         cs=cs-L*FLOOR(cs/L) !periodic boundary condition
         !write(6,*) cs,"..pbc"
         jsite=LocateNeighSite(isite,cs)
         !write(6,*) "jsite",jsite
         !stop
         IF (jsite>0) THEN !jsite=0 implies site not found
            !write(*,*) jsite,":jsite"
            IF (speciestype(jsite)>0) THEN !actual atom present
               IF (PrintableAtom(jsite)) THEN !site lying in plane is also surface atom
                  n=n+1
                  !write(*,*) n
               END IF
            END IF
         END IF
      END DO
      
      !finalize
      IF (n>3) THEN !at least 3 NN should be surface atoms in 111 plane
         CharacterizeSurface100=iplane
         !stop
         RETURN
      END IF
   END DO
   !stop

   CharacterizeSurface100=0
END FUNCTION CharacterizeSurface100
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
FUNCTION LocateNeighSite(isite,coord) !locate neighbor of isite which has
   !coordinate coord
   IMPLICIT NONE
   INTEGER :: isite,LocateNeighSite,k,ksite
   REAL :: coord(3),c(3),dcoord(3)
   REAL, PARAMETER :: tol=0.100000
   
   DO k=1,NumberNeighborSites(isite)
      ksite=NeighborSiteList((isite-1)*12+k)
      c=AtomCoord(3*ksite-2:3*ksite)
      dcoord=ABS(coord-c)
      !write(6,*) c,"_"
     ! IF (ANY(dcoord>15.0)) THEN
         !WRITE(6,*) "something wierd ..."
         !STOP
     ! END IF
      IF (ALL(dcoord<tol)) THEN
         LocateNeighSite=ksite
         RETURN
      END IF
   END DO
   LocateNeighSite=0
END FUNCTION LocateNeighSite
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE CharacterizeSubsurface()
   IMPLICIT NONE
   CALL CharacterizeSubsurface111()
   CALL CharacterizeSubsurface100()
END SUBROUTINE CharacterizeSubsurface
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE CharacterizeSubsurface111()
   IMPLICIT NONE
   INTEGER :: isite,k,ksite

   IsSubSurf111atom=0
   DO isite=1,NSites
      IF (coordination(isite)<12) CYCLE !this forces us to look at bulk atoms only
      DO k=1,NumberNeighborSites(isite)
         ksite=NeighborSiteList((isite-1)*12+k)
         IF (IsSurf111atom(ksite)>0) THEN !neighbor is 111 atom
            IsSubSurf111atom(isite)=1
         END IF
      END DO
   END DO
END SUBROUTINE CharacterizeSubsurface111
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE CharacterizeSubsurface100()
   IMPLICIT NONE
   INTEGER :: isite,k,ksite

   IsSubSurf100atom=0
   DO isite=1,NSites
      IF (coordination(isite)<12) CYCLE !this forces us to look at bulk atoms only
      DO k=1,NumberNeighborSites(isite)
         ksite=NeighborSiteList((isite-1)*12+k)
         IF (IsSurf100atom(ksite)>0) THEN !neighbor is 111 atom
            IsSubSurf100atom(isite)=1
         END IF
      END DO
   END DO
END SUBROUTINE CharacterizeSubsurface100
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
FUNCTION NumberFe(isite)
   IMPLICIT NONE
   INTEGER :: NumberFe,isite,k,jsite,n,NFe

   IF (speciestype(isite)==0) RETURN !site is already vacant

   !Find number of atoms present in the neighbor site
   NumberFe=0
   NFe=0 !added by aditya
   n=NumberNeighborSites(isite)
   DO k=1,n
      jsite=NeighborSiteList((isite-1)*12+k)
      !write(6,*) "Site #:",k,jsite
      !call flush(6)
      IF (jsite>0) THEN
         IF(speciestype(jsite)==1) NFe=NFe+1  !added by aditya
      END IF
   END DO

   NumberFe=NFe

END FUNCTION NumberFe
  !xxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Store_local_composition(AL1)

   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: AL1
   INTEGER :: i

   OPEN(UNIT=234,FILE="Local_composition")
   DO i=1,NSites
      WRITE(234,*) Local_composition(i)
   END DO
   CLOSE(234)
END SUBROUTINE Store_local_composition
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Read_local_composition(AL1,AL)

   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: AL1,AL
   INTEGER :: i,pos,ncells(3)
   REAL(dp) :: coord(3),cellsize(3)
   REAL(dp),ALLOCATABLE,DIMENSION(:) :: Local_composition0

   ALLOCATE(Local_composition0(AL1%NAtoms))
   OPEN(UNIT=234,FILE="Local_composition")
   DO i=1,AL1%NAtoms
      READ(234,*) Local_composition0(i)
   END DO
   CLOSE(234)

   cellsize=AL%LL%CellSize
   ncells=AL%LL%NumberCell3D

   Local_composition=0
   DO i=1,AL1%NAtoms
      coord=AL1%AtomCoord(3*i-2:3*i)
      pos=GetIndex(AL,coord,cellsize,ncells)
      Local_composition(pos)=Local_composition0(i)
   END DO

END SUBROUTINE Read_local_composition
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Fill_defect_sites
   IMPLICIT NONE
   INTEGER :: i,j,k,nEN

!   nEN=0
!   DO i=1,NSites
!      IF (AL%AtomSpecies(i)==2) nEN=nEN+1
!   END DO
!
!WRITE(6,*) nEN

!   OPEN(UNIT=236,FILE="tauss-values")
!   OPEN(UNIT=237,FILE="Local_composition1")

   DO j=1,11
      k=42-j
      DO i=1,NSites
         IF (AL%AtomSpecies(i)==0) THEN
            IF (ThirdNearestNeighborSites(i)>k) THEN
               IF (taus88()>Local_composition(i)) THEN
!                   WRITE(236,*) taus88()
!                   WRITE(237,*) i
                   AL%AtomSpecies(i)=1  !electroactive species inserted
               ELSE
!                   IF (nEN<52367) THEN
                   AL%AtomSpecies(i)=2
!                      nEN=nEN+1
!                   ELSE
!                      AL%AtomSpecies(i)=1
!                   END IF
               END IF
            END IF
         END IF
      END DO
   END DO
!   CLOSE(236)
!   CLOSE(237)

END SUBROUTINE Fill_defect_sites
   !xxxxxxxxxxxxxxxxxxxxxxxxxx

