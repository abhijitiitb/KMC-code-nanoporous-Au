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
SUBROUTINE NearestNeighborSiteList(isite)
!find list of index which are neighbors
   IMPLICIT NONE
   INTEGER :: isite,j,jsite,n,ra,r0,r1,r2,s1,s2
   REAL(dp) ::drmagij

   n=0 !number of sites
   
   ra=(isite-1)*MaxAtomPerAtom !start index minus 1
   r0=VLListRange(isite) !range for ith atom
   r1=ra+1
   r2=ra+r0
   s1=AtomSpecies(isite)

   DO j=r1,r2 !this will take care of the pair and many-body terms in one go
      jsite=VLList(j)
      s2=AtomSpecies(jsite)
      drmagij=VLdrmag(j)   ! spacing between atoms i and j
      IF(drmagij<=3.01_dp) THEN
         n=n+1
         NeighborSiteList((isite-1)*12+n)=jsite
      END IF  
   END DO
   NumberNeighborSites(isite)=n
END SUBROUTINE NearestNeighborSiteList
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE InitializeRates()
   IMPLICIT NONE
   REAL(dp), PARAMETER :: kboltzmann=8.62e-5
   INTEGER :: NAtm
   
   RateDiffusion=0.d0
   RateDissolution=0.d0
   DO NAtm=0,11
      RateDiffusion(NAtm)=prefac*(EXP(-(NAtm*xi)/(kboltzmann*Temperature1)))
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

   !Satisfied=.FALSE.
   !IF (.NOT. Satisfied) RETURN
   iframe=iframe+1
   !OPEN(UNIT=302,FILE="NW."//TRIM(INT2CHAR(iframe))//".xyz") !for analysis
   !CALL GetPrintableAtoms(.FALSE.)
   !WRITE(6,*) "Printing xyz for nanowire ..."
   !CALL PrintFrame1()
   !CLOSE(302)

   OPEN(UNIT=302,FILE="NWShell."//TRIM(INT2CHAR(iframe))//".xyz") !for plotting
   CALL GetPrintableAtoms(.TRUE.)
   WRITE(6,*) "Printing xyz for nanowire shell...",iframe
   CALL PrintFrame1()
   CALL CharacterizeSurface() !identifies (111), (100) and other surface atoms
     !among atoms which are printable
   CALL PrintSurface(111,466) !
   CALL PrintSurface(100,467)
   CALL PrintSurface(0,468) !others
   CALL PrintFacetAreaDistribution(111)
   CALL PrintFacetAreaDistribution(100)
   CLOSE(302)
IF (iframe>250) THEN
WRITE(6,*) "iframe more than 250"
STOP
END IF
END SUBROUTINE PrintFrame
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE PrintFrame1()
   IMPLICIT NONE

   WRITE(UNIT=302,FMT='(4i10," NAtomsInFrame, NSites, NElectroactive, Nonactive")')  &
     NPrintableAtoms,NSites, &
     NElectroactive-NDissolved,NSites-NElectroactive
   WRITE(UNIT=302,FMT='(4es15.5)') AL%BoxSize(1:3),time
   DO isite=1,NSites
      IF (PrintableAtom(isite)) &
         WRITE(UNIT=302,FMT='(a4,3f20.8)') &
          PrintSpeciesType(isite),AtomCoord(3*isite-2:3*isite)
   END DO
END SUBROUTINE PrintFrame1
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
FUNCTION PrintSpeciesType(isite)
   IMPLICIT NONE
   CHARACTER(len=4) :: PrintSpeciesType
   INTEGER :: isite

   SELECT CASE (speciestype(isite))
   CASE(1); PrintSpeciesType="Fe  "
   CASE(2); PrintSpeciesType="Pt  "
   CASE DEFAULT
      WRITE(6,*) "Unexpected vacancy encountered"
      STOP
   END SELECT
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
   INTEGER :: isite

   PrintableAtom=.FALSE.
   IF (PrintShell) THEN
      NPrintableAtoms=0
      DO isite=1,NSites
         IF (speciestype(isite)>0) THEN
         CALL NeighborIsUndercoordinated(isite,UnderCoordNeigh)
         PrintableAtom(isite)=(Coordination(isite)<12 .OR. UnderCoordNeigh)
         IF (PrintableAtom(isite)) NPrintableAtoms=NPrintableAtoms+1
         END IF
      END DO
   ELSE
      NPrintableAtoms=NSites-NDissolved !total number of atoms present
      DO isite=1,NSites
         PrintableAtom(isite)=speciestype(isite)>0
      END DO
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
      IF (PrintableAtom(isite)) THEN !lets look at its neighbors
         IsSurf111atom(isite)=CharacterizeSurface111(isite) !check if this is 111 surface
         IsSurf100atom(isite)=CharacterizeSurface100(isite) !check if this is 111 surface
         IF (IsSurf111atom(isite)>0 .AND. IsSurf100atom(isite)>0) THEN
            WRITE(6,*) "Atom found to be both 111 and 100 surface atoms!"
         END IF
         IF (.NOT. (IsSurf111atom(isite)>0 .OR. IsSurf100atom(isite)>0)) THEN
            IsSurfOtherAtom(isite)=1 !some other type of surface
         END IF
      END IF
   END DO
END SUBROUTINE CharacterizeSurface
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
FUNCTION CharacterizeSurface111(isite)
!NOTE: if atom lies at intersection of two 111 planes then only one will be assigned
!0 implies does not lie on 111 plane
   IMPLICIT NONE
   INTEGER :: isite,iplane,jsite,CharacterizeSurface111,ivec
   REAL :: z,a,vec(18),c0(3),cs(3)
   REAL, DIMENSION(3), PARAMETER :: L=(/244.806000,244.806000,734.418000/)
   LOGICAL :: IsMatch

   a=2.8600000000
   z=0.00000000000
   DO iplane=1,4
      SELECT CASE(iplane) !select plane
      CASE(1); vec=(/-a,z,-a, a,z,a, -a,-a,z, a,a,z, z,-a,a, z,a,-a/)
      CASE(2); vec=(/-a,z,a, a,z,-a, a,-a,z, -a,a,z, z,-a,a, z,a,-a/)
      CASE(3); vec=(/-a,z,-a, a,z,a, -a,a,z, a,-a,z, z,a,a, z,-a,-a/)
      CASE(4); vec=(/-a,z,a, a,z,-a, -a,-a,z, a,a,z, z,a,a, z,-a,-a/)
      END SELECT

      !search coordinate
      c0=AtomCoord(3*isite-2:3*isite) !coordinate of site
      IsMatch=.TRUE.
      DO ivec=1,6
         cs=c0+vec(3*ivec-2:3*ivec)
         cs=cs-L*NINT(cs/L) !periodic boundary condition
         jsite=LocateNeighSite(isite,cs)
         IF (jsite>0) THEN !jsite=0 implies site not found
            IsMatch=IsMatch .AND. PrintableAtom(jsite) !site lying in plane is also surface atom
         END IF
      END DO
      
      !finalize
      IF (IsMatch) THEN
         CharacterizeSurface111=iplane
         RETURN
      END IF
   END DO
   CharacterizeSurface111=0
END FUNCTION CharacterizeSurface111
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
FUNCTION CharacterizeSurface100(isite)
   IMPLICIT NONE
   INTEGER :: isite,CharacterizeSurface100

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
      IF (ALL(dcoord<tol)) THEN
         LocateNeighSite=ksite
         RETURN
      END IF
   END DO
   LocateNeighSite=0
END FUNCTION LocateNeighSite
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
