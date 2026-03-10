!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE CharacterizeFacets()
   IMPLICIT NONE
   INTEGER :: isurf,nsp111_1,nsp111_2,nsp100_1,nsp100_2,nspother_1,nspother_2
   INTEGER :: nspbulk_1,nspbulk_2
   INTEGER :: nfacets,nfacetatoms(100000)
   REAL :: avfacets,sb,s111,s100,soth,s

   nfacetatoms=0; nfacets=0
   CALL DetectFacet(111,nfacets,nfacetatoms)
   avfacets=0
   nsp111_1=0 !#atoms species1
   nsp111_2=0 !#atoms species2
   DO isite=1,NSites
      IF (IsSurf111atom(isite)>0) CALL CountSpecies(isite,nsp111_1,nsp111_2)
   END DO
   IF (nfacets>0) avfacets=SUM(nfacetatoms(1:nfacets))/REAL(nfacets)
   WRITE(6,*) "111 Facet info:",nfacets,avfacets,MAXVAL(nfacetatoms),"<<<<<"
   WRITE(546,*)   kmcx,kmcd,time,nfacets,avfacets,MAXVAL(nfacetatoms), &
     SUM(nfacetatoms(1:nfacets))

   nfacetatoms=0; nfacets=0
   CALL DetectFacet(100,nfacets,nfacetatoms)
   avfacets=0
   nsp100_1=0 !#atoms species1
   nsp100_2=0 !#atoms species2
   DO isite=1,NSites
      IF (IsSurf100atom(isite)>0) CALL CountSpecies(isite,nsp100_1,nsp100_2)
   END DO
   IF (nfacets>0) avfacets=SUM(nfacetatoms(1:nfacets))/REAL(nfacets)
   WRITE(6,*) "100 Facet info:",nfacets,avfacets,MAXVAL(nfacetatoms),"<<<<<"
   WRITE(547,*) kmcx,kmcd,time,nfacets,avfacets,MAXVAL(nfacetatoms), &
      SUM(nfacetatoms(1:nfacets))

   !also lets add bulk information
   nspbulk_1=0 !#atoms species1
   nspbulk_2=0 !#atoms species2
   DO isite=1,NSites
      IF (Coordination(isite)==12) CALL CountSpecies(isite,nspbulk_1,nspbulk_2)
   END DO

   !also lets add other surface atom information
   !edge/corner atoms ... subsurface (CN=11) not included
   nspother_1=0 !#atoms species1
   nspother_2=0 !#atoms species2
   DO isite=1,NSites
      IF (Coordination(isite)<11 .AND. .NOT.(IsSurf111atom(isite)>0 .OR. IsSurf100atom(isite)>0)) &
        CALL CountSpecies(isite,nspother_1,nspother_2)
   END DO

   s=nspbulk_1+nspbulk_2+nsp111_1+nsp111_2+nsp100_1+nsp100_2+nspother_1+nspother_2
   sb=MAX(nspbulk_1+nspbulk_2,1) !so that division by 0 is not performed
   s111=MAX(nsp111_1+nsp111_2,1)
   s100=MAX(nsp100_1+nsp100_2,1)
   soth=MAX(nspother_1+nspother_2,1)
   WRITE(551,*) kmcx,kmcd,time,&
      nspbulk_1,nspbulk_2,nsp111_1,nsp111_2, &
      nsp100_1,nsp100_2,nspother_1,nspother_2,s,sb,s111,s100,soth, &
      REAL(nspbulk_1)/sb,REAL(nspbulk_2)/sb,REAL(nsp111_1)/s111,REAL(nsp111_2)/s111, &
      REAL(nsp100_1)/s100,REAL(nsp100_2)/s100,REAL(nspother_1)/soth,REAL(nspother_2)/soth

END SUBROUTINE CharacterizeFacets
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CountSpecies(isite,nsp1,nsp2)
      INTEGER :: nsp1,nsp2,isite
      SELECT CASE(speciestype(isite))
      CASE(1); nsp1=nsp1+1
      CASE(2); nsp2=nsp2+1
      END SELECT
   END SUBROUTINE CountSpecies
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE DetectFacet(isurf,nfacets,nfacetatoms)
!loop through all sites & locate atoms belonging to each 111 facet
   IMPLICIT NONE
   INTEGER :: isite,n,isurf,iplane,nf
   INTEGER, DIMENSION(:), POINTER :: arr
   LOGICAL, DIMENSION(:), POINTER :: mask
   REAL, DIMENSION(:), POINTER :: vec
   INTEGER :: nfacets,nfacetatoms(100000)
   CHARACTER(len=20) :: cplane,cplane0

   ALLOCATE(mask(NSites)) !whether site has to be analyzed
   mask=.TRUE. !all have to be analyzed
   SELECT CASE(isurf)
   CASE (111)
      arr=>IsSurf111atom; n=18 !contains surface plane info
      cplane0="[111]"
   CASE (100)
      arr=>IsSurf100atom; n=12 !n is # directions for plane
      cplane0="[100]"
   END SELECT
   ALLOCATE(vec(n))

   nfacets=0
   DO isite=1,NSites
      IF (mask(isite) .AND. arr(isite)>0 .AND. coordination(isite)<11) THEN !site has not been considered
         iplane=arr(isite) !tells which plane is being considered
         CALL GetFacetCoordinates(vec,isurf,iplane,n) !vectors for iplane
         nf=0 !# of atoms in facet containing isite
         nfacets=nfacets+1
         cplane=TRIM(cplane0)//TRIM(INT2CHAR(iplane))
         CALL TagFacetAtoms(isite,iplane,vec,mask,arr,n,nf,nfacets,cplane)
         nfacetatoms(nfacets)=nf
      END IF
   END DO
   DEALLOCATE(mask,vec)
END SUBROUTINE DetectFacet
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx`
SUBROUTINE GetFacetCoordinates(vec,isurf,iplane,n)
!finds the NN coordinates for a given facet
   IMPLICIT NONE
   INTEGER :: isurf,iplane,n
   REAL, DIMENSION(n) :: vec
   REAL :: z,a,c1,c2,c3,c4
   
   a=2.0400000000
   z=0.00000000000
   IF (NWOrientation==1) THEN  !(100) oriented NW
      IF (isurf==111) THEN !(100) oriented NW
         SELECT CASE(iplane) !select plane
         CASE(1); vec=(/-a,z,-a, a,z,a, -a,-a,z, a,a,z, z,-a,a, z,a,-a/)
         CASE(2); vec=(/-a,z,a, a,z,-a, a,-a,z, -a,a,z, z,-a,a, z,a,-a/)
         CASE(3); vec=(/-a,z,-a, a,z,a, -a,a,z, a,-a,z, z,a,a, z,-a,-a/)
         CASE(4); vec=(/-a,z,a, a,z,-a, -a,-a,z, a,a,z, z,a,a, z,-a,-a/)
         CASE DEFAULT; WRITE(6,*) "Unexpected err1"; STOP 
         END SELECT
      ELSEIF (isurf==100) THEN !(100) oriented NW
         SELECT CASE(iplane) !select plane
         CASE(1); vec=(/-a,z,-a, a,z,a, -a,z,a, a,z,-a/)
         CASE(2); vec=(/-a,a,z, a,-a,z, -a,-a,z, a,a,z/)
         CASE(3); vec=(/z,-a,a, z,a,-a, z,-a,-a, z,a,a/)
         CASE DEFAULT; WRITE(6,*) "Unexpected err2"; STOP 
         END SELECT
      END IF
   
   ELSEIF (NWOrientation==2) THEN !(110) oriented NW
      c1=1.443
      c2=2.041
      c3=2.885
      c4=0.
      IF (isurf==111) THEN 
         SELECT CASE(iplane)
         CASE(1); vec=(/-c1,-c2,-c1, c4,c4,-c3, c1,c2,-c1, c1,c2,c1, c4,c4,c3, -c1,-c2,c1/)
         CASE(2); vec=(/c4,c4,-c3, c1,-c2,-c1, c1,-c2,c1, c4,c4,c3, -c1,c2,c1, -c1,c2,-c1/)
         CASE(3); vec=(/c3,c4,c4, c1,-c2,-c1, -c1,-c2,-c1, -c3,c4,c4, -c1,c2,c1, c1,c2,c1/)
         CASE(4); vec=(/c3,c4,c4, c1,-c2,c1, -c1,-c2,c1, -c3,c4,c4, -c1,c2,-c1, c1,c2,-c1/)
         CASE DEFAULT; WRITE(6,*) "Unexpected err1"; STOP 
         END SELECT
      ELSEIF (isurf==100) THEN 
         SELECT CASE(iplane)
         CASE(1); vec=(/c4,c4,-c3, -c3,c4,c4, c4,c4,c3, c3,c4,c4/)
         CASE(2); vec=(/-c1,c2,c1, c1,c2,-c1, c1,-c2,-c1, -c1,-c2,c1/)
         CASE(3); vec=(/c1,c2,c1, -c1,c2,-c1, -c1,-c2,-c1, c1,-c2,c1/)
         CASE DEFAULT; WRITE(6,*) "Unexpected err1"; STOP 
         END SELECT
      END IF
   END IF


END SUBROUTINE GetFacetCoordinates
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx`
RECURSIVE SUBROUTINE TagFacetAtoms(isite,iplane,vec,mask,arr,nvec,nf, &
  nfacets,cplane)
!isite - site that needs to be tagged
!iplane - selected plane from selected surface isurf
!isurf - selected surface (111) or (100)
!mask - list of sites that need to be analyzed
!arr - plane type for each site
   IMPLICIT NONE
   INTEGER :: isite,jsite,iplane,nvec,nf,ivec,nfacets
   REAL, DIMENSION(nvec) :: vec
   REAL :: c0(3),cs(3)
   REAL, DIMENSION(3), PARAMETER :: L=(/244.806000,244.806000,734.418000/)
   LOGICAL, DIMENSION(:), POINTER :: mask
   INTEGER, DIMENSION(:), POINTER :: arr
   CHARACTER(len=20) :: cplane

   mask(isite)=.FALSE.
   nf=nf+1 !# facet atoms
   c0=AtomCoord(3*isite-2:3*isite)
   SurfaceInfo(isite)=TRIM(SurfaceInfo(isite))//TRIM(cplane)//"-"//TRIM(INT2CHAR(nfacets))
   !WRITE(6,*) SurfaceInfo(isite)
   DO ivec=1,nvec/3
      cs=c0+vec(3*ivec-2:3*ivec)
      cs=cs-L*FLOOR(cs/L) !periodic boundary condition
      jsite=LocateNeighSite(isite,cs)
      IF (mask(jsite) .AND. arr(jsite)==iplane .AND. coordination(isite)<11) &
         CALL TagFacetAtoms(jsite,iplane,vec,mask,arr,nvec,nf,nfacets,cplane)
   END DO
END SUBROUTINE TagFacetAtoms
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE AnalyzeNWFrame()
   IMPLICIT NONE
   INTEGER :: nligament,iligament
   
   WRITE(6,*) "Analyzing ligaments ..."
   ALLOCATE(LigamentIndex(NSites))
   CALL FindLigaments(nligament)
   !nligament=1
   !LigamentIndex(1:NSites)=1
   
   !STOP

   WRITE(6,*) "Analyzing size distribution of ligaments..."
   DO iligament=1,nligament
      CALL GetSizeDistribution(iligament) 
   END DO

   DEALLOCATE(LigamentIndex)
END SUBROUTINE AnalyzeNWFrame
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE FindLigaments(nligament)
   IMPLICIT NONE
   INTEGER :: isite,iligament,taglvl,i,n,nligament
   INTEGER, DIMENSION(:), ALLOCATABLE :: LigamentAtomCount

   LigamentIndex=0
   iligament=0
   DO isite=1,NSites
      IF (LigamentIndex(isite)==0 .AND. SpeciesType(isite)>0) THEN
         taglvl=0
         iligament=iligament+1
         CALL TagLigamentSite(isite,iligament,taglvl)
      END IF
   END DO

   !retagging
   ALLOCATE(LigamentAtomCount(iligament))
   LigamentAtomCount=0
   DO isite=1,NSites
      i=LigamentIndex(isite)
      LigamentAtomCount(i)=LigamentAtomCount(i)+1
   END DO
   nligament=0
   DO i=1,iligament
      IF (LigamentAtomCount(i)>0) THEN !atoms are present in ligament
         nligament=nligament+1
         IF (nligament<iligament) THEN
            CALL ResetLigamentIndex(i,nligament)
         ELSEIF (nligament>iligament) THEN
         END IF
      END IF
   END DO

   WRITE(6,*) "Number of ligaments found:",nligament

!   WRITE(6,*) "Printing ligament files"
!   LigamentAtomCount=0
!   DO isite=1,NSites
!      i=LigamentIndex(isite)
!      LigamentAtomCount(i)=LigamentAtomCount(i)+1
!  END DO
!   DO i=1,nligament
!      WRITE(6,*) "Ligament:",i," #sites:",LigamentAtomCount(i)
!      CALL PrintLigament(i,LigamentAtomCount(i))
!   END DO
!   
!   DEALLOCATE(LigamentAtomCount)
END SUBROUTINE FindLigaments
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
RECURSIVE SUBROUTINE TagLigamentSite(isite,iligament,taglvl)
   INTEGER :: isite,iligament,taglvl,jligament
   INTEGER :: j,jsite

   LigamentIndex(isite)=iligament

   taglvl=taglvl+1
   IF (taglvl>5) RETURN
   
   !check neighbors
   DO j=1,NumberNeighborSites(isite)
      jsite=NeighborSiteList((isite-1)*12+j)
      IF (LigamentIndex(jsite)>0) THEN !this site has already been considered
         jligament=LigamentIndex(jsite)
         CALL ResetLigamentIndex(jligament,iligament) !since we are searching iligament
         !recursively it is difficult to change iligament on the fly, we can
         !change jligament to iligament more easily
      ELSEIF (LigamentIndex(jsite)==0  .AND. SpeciesType(isite)>0) THEN
         CALL TagLigamentSite(jsite,iligament,taglvl)
      END IF
   END DO

END SUBROUTINE TagLigamentSite
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE ResetLigamentIndex(i,j)
!replace i with j
   IMPLICIT NONE
   INTEGER :: i,j

   WHERE (LigamentIndex==i) LigamentIndex=j
END SUBROUTINE ResetLigamentIndex
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE PrintLigament(i,n)
   IMPLICIT NONE
   INTEGER :: i,n
   CHARACTER(len=100) :: filename1
   CHARACTER(len=NINT2CHAR) :: txt

   txt=INT2CHAR(i)
   filename1="Ligament."//TRIM(txt)//".xyz"
   OPEN(UNIT=435,FILE=TRIM(filename1))
   WRITE(435,*) n
   WRITE(435,*) "Ligaments for file:",TRIM(filename)
   DO isite=1,NSites
      IF (LigamentIndex(isite)==i) THEN
         WRITE(UNIT=435,FMT='(a4,3f20.8,i8,"   ",a20)') &
            PrintSpeciesType(isite,0),AtomCoord(3*isite-2:3*isite)
      END IF
   END DO
   CLOSE(435)
END SUBROUTINE PrintLigament
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE GetSizeDistribution(iligament)
!finds the size distribution for a ligament
   IMPLICIT NONE
   INTEGER, PARAMETER :: nrot=4000
   INTEGER :: iligament,iattempt,nlocations,n,nsel,isite,c
   INTEGER, DIMENSION(1000000) :: LigamentSiteList
   REAL(dp), DIMENSION(3,3*nrot) :: rmatrices
   REAL(dp), DIMENSION(3,nrot) :: eulerangles

   c=0
   DO isite=1,NSites
      IF (LigamentIndex(isite)==iligament) c=c+1
   END DO
   WRITE(6,*) "Number of sites belonging to the ligament:",c
   
!   CALL SetupLigamentSites(iligament,LigamentSiteList,n)!Setup list of sites belonging to a ligament
!   CALL SetupRotationMatrices(rmatrices,nrot,eulerangles)
!   IF (c<100) RETURN
!
!   nlocations=100
!   nsel=0 !#sites studied/selected
!   DO iattempt=1,nlocations !check for size at nlocations
!      IF (n==nsel) THEN
!         WRITE(6,*) "All sites have been studied"
!         EXIT
!      END IF
!      WRITE(6,*) "xxxxxxxxxxxxxxxxAttempt ",iattempt,"xxxxxxxxxxxxxxxxxxxxxx"
!      WRITE(6,*) "# sites studied so far:",nsel,"/",n
!      CALL SelectLigamentSite(LigamentSiteList,n,nsel,isite) !select isite randomly from ligament
!        !nsel gives the number of sites that have already been studied
!        !1:nsel are sites in LigamentSiteList that have been studied
!      WRITE(6,*) "Site,ligament:",isite,iligament
!      CALL GetCrossSection(isite,LigamentSiteList,n,nsel,rmatrices,nrot,eulerangles)!find the plane which has the fewest atoms in the cross-section
!      !mark all atoms within cross-section
!      !note down the diameter
!   END DO
END SUBROUTINE GetSizeDistribution
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE SetupLigamentSites(i,LigamentSiteList,n)
!n denotes the number of sites belonging to the ligament i
   IMPLICIT NONE
   INTEGER :: i,n,s
   INTEGER, DIMENSION(1000000), INTENT(INOUT) :: LigamentSiteList

   n=0
   s=SIZE(LigamentSiteList)
   DO isite=1,NSites
      IF (LigamentIndex(isite)==i) THEN
         n=n+1
         IF (n>s) THEN
            WRITE(6,*) "Array size exceeded in SUBROUTINE SetupLigamentSites"
            STOP
         END IF
         LigamentSiteList(n)=isite
      END IF
   END DO
END SUBROUTINE SetupLigamentSites
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE SelectLigamentSite(arr,n,ns,i)
   IMPLICIT NONE
   INTEGER :: n,ns,i,si,sj
   INTEGER, DIMENSION(:) :: arr
   REAL :: r

   r=taus88()
   i=INT(r*REAL(n-ns)+ns) !position of site selected
   !next we swap positions and bring site i to position ns+1
   !site at ns+1 will be moved to the original position of i
   si=arr(i)
   ns=ns+1 !number of selcted sites incremented by 1
   sj=arr(ns)
   arr(i)=sj
   arr(ns)=si
   i=si

END SUBROUTINE SelectLigamentSite
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE GetCrossSection(isite,arr,n,ns,rmatrices,nrot,e)
   IMPLICIT NONE
   INTEGER :: n,ns,isite,jsite,ksite,irot,area,minarea,j,ix,iy,nx
   INTEGER :: srot,nrot
   REAL(dp), INTENT(IN) :: rmatrices(3,3*nrot)
   REAL(dp), INTENT(IN) :: e(3,nrot)
   INTEGER, DIMENSION(:) :: arr
   LOGICAL :: Pixels(-50:50,-50:50),PixelsOfInterest(-50:50,-50:50),success
   REAL(dp) :: ci(3),cj(3),rmat(3,3),d,BoxSize(3),com(3)

   ci=AtomCoord(3*isite-2:3*isite)
   WRITE(6,*) "ci :",ci
   minarea=1e8 !initialize
   srot=0
   BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
   WRITE(6,*) "BoxSize :",BoxSize
   WRITE(6,*) "nrot :",nrot

   d=3.2 !thickness of disk
   DO irot=1,nrot !loop over all rotation matrices
      rmat=GetRotationMatrix(rmatrices,irot,nrot)
      Pixels=.FALSE. !reset all pixels to zero value
      PixelsOfInterest=.FALSE.
      success=.TRUE.
      DO j=1,n !loop over all sites
         jsite=arr(j)
         cj=AtomCoord(3*jsite-2:3*jsite)-ci !set ci at
         cj=cj-BoxSize*NINT(cj/BoxSize)
         cj=MATMUL(rmat,cj)
         IF (ABS(cj(3))<d) THEN !this site lies within a disk of 2d thickness
            !note that site present in the pixel
            ix=FLOOR(cj(1)/3.2)
            iy=FLOOR(cj(2)/3.2)
            IF (ABS(ix)<=50 .AND. ABS(iy)<=50) Pixels(ix,iy)=.TRUE.
         END IF
      END DO
      area=0
      CALL AnalyzePixels(0,0,Pixels,PixelsOfInterest,area)

      IF (area<minarea) THEN
         srot=irot !selected rotation
         minarea=area
      END IF
   END DO

   !print cross-sections and find whether we should proceed
   DO irot=srot,srot
      rmat=GetRotationMatrix(rmatrices,irot,nrot)
      Pixels=.FALSE. !reset all pixels to zero value
      PixelsOfInterest=.FALSE.
      success=.TRUE.
      DO j=1,n !loop over all sites
         jsite=arr(j)
         cj=AtomCoord(3*jsite-2:3*jsite)-ci !set ci at
         cj=cj-BoxSize*NINT(cj/BoxSize)
         cj=MATMUL(rmat,cj)
         IF (ABS(cj(3))<d) THEN !this site lies within a disk of 2d thickness
            !note that site present in the pixel
            ix=FLOOR(cj(1)/3.2)
            iy=FLOOR(cj(2)/3.2)
            IF (ABS(ix)<=50 .AND. ABS(iy)<=50) Pixels(ix,iy)=.TRUE.
         END IF
      END DO
!~       DO ix=-30,30
!~          DO iy=-30,30
!~             IF (Pixels(ix,iy)) THEN
!~                                 WRITE(UNIT=6,FMT='(A2)',ADVANCE="NO") " 1"
!~             ELSE
!~                                 WRITE(UNIT=6,FMT='(A2)',ADVANCE="NO") " 0"
!~                         END IF
!~          END DO
!~          WRITE(6,*) ""
!~       END DO
      area=0
      CALL AnalyzePixels(0,0,Pixels,PixelsOfInterest,area)
      
      DO j=1,n
         jsite=arr(j)
         cj=AtomCoord(3*jsite-2:3*jsite)-ci
         cj=MATMUL(rmat,cj)
         IF (ABS(cj(3))<d .AND. j>ns) THEN !mark it
            ix=FLOOR(cj(1)/2.5)
            iy=FLOOR(cj(2)/2.5)
            IF (ix<=50 .AND. iy<=50) THEN
               IF (PixelsOfInterest(ix,iy) .AND. j<=ns .AND. jsite/=isite) THEN
                  success=.FALSE. !this site has already been selected
                  WRITE(6,*) "oops ..site",jsite," is at position ",j
               END IF
            END IF
         END IF
      END DO
      IF (.NOT. success) THEN
         WRITE(6,*) "Atoms belonging to cross-section already analyzed .. search aborted"
         RETURN
      END IF
      
      WRITE(6,*) "-----After processing-----"
!~       DO ix=-30,30
!~          DO iy=-30,30
!~             IF (PixelsOfInterest(ix,iy)) THEN
!~                                 WRITE(UNIT=6,FMT='(A2)',ADVANCE="NO") " 1"
!~             ELSE
!~                                 WRITE(UNIT=6,FMT='(A2)',ADVANCE="NO") " 0"
!~                         END IF
!~          END DO
!~          WRITE(6,*) ""
!~       END DO
   END DO
   WRITE(6,*) "minimum area:",minarea

   !selected sites should be tagged as selected
   IF (srot>0) THEN
      rmat=GetRotationMatrix(rmatrices,srot,nrot)
      WRITE(6,*) "phi:",e(1,srot)
      WRITE(6,*) "theta:",e(2,srot)
      WRITE(6,*) "psi:",e(3,srot)
      WRITE(6,*) "Area :",minarea
      WRITE(UNIT=6,FMT='("Atoms marked ...")',ADVANCE="no")
      WRITE(UNIT=6,FMT='(i8)',ADVANCE="no") isite-1
      !mark atoms and find center of mass
      nx=1 !isite already included
      com=ci
      DO j=1,n
         jsite=arr(j)
         cj=AtomCoord(3*jsite-2:3*jsite)-ci
         cj=MATMUL(rmat,cj)
         IF (ABS(cj(3))<d .AND. j>ns) THEN !mark it
            ix=FLOOR(cj(1)/2.5)
            iy=FLOOR(cj(2)/2.5)
            IF (ix<=50 .AND. iy<=50) THEN
               IF (PixelsOfInterest(ix,iy)) THEN
                  ns=ns+1
                  ksite=arr(ns)
                  arr(ns)=jsite
                  arr(j)=ksite
                  nx=nx+1
                  com=com+AtomCoord(3*jsite-2:3*jsite)
                  WRITE(UNIT=6,FMT='(i8)',ADVANCE="no") jsite-1
               END IF
            END IF
         END IF
      END DO
      WRITE(6,*) ""
      com=com/REAL(nx) !center of mass
      WRITE(6,*) "Center of mass:",com
   END IF
END SUBROUTINE GetCrossSection
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
!!!!SUBROUTINE GetCrossSection(isite,arr,n,ns,rmatrices,nrot,e)
!!!!   IMPLICIT NONE
!!!!   INTEGER :: n,ns,isite,jsite,ksite,irot,area,minarea,j,ix,iy,nx
!!!!   INTEGER :: srot,nrot
!!!!   REAL(dp), INTENT(IN) :: rmatrices(3,3*nrot)
!!!!   REAL(dp), INTENT(IN) :: e(3,nrot)
!!!!   INTEGER, DIMENSION(:) :: arr
!!!!   LOGICAL :: Pixels(-50:50,-50:50),PixelsOfInterest(-50:50,-50:50)
!!!!   LOGICAL :: success,AcceptableCrossSection
!!!!   REAL(dp) :: ci(3),cj(3),rmat(3,3),d,BoxSize(3),com(3)
!!!!
!!!!   ci=AtomCoord(3*isite-2:3*isite)
!!!!   minarea=1e8 !initialize
!!!!   srot=0
!!!!   BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
!!!!   d=3. !thickness of disk
!!!!   DO irot=1,nrot !loop over all rotation matrices
!!!!      rmat=GetRotationMatrix(rmatrices,irot,nrot)
!!!!      Pixels=.FALSE. !reset all pixels to zero value
!!!!      PixelsOfInterest=.FALSE.
!!!!      success=.TRUE.
!!!!      DO j=1,n !loop over all sites
!!!!         jsite=arr(j)
!!!!         cj=AtomCoord(3*jsite-2:3*jsite)-ci !set ci at
!!!!         cj=cj-BoxSize*NINT(cj/BoxSize)
!!!!         cj=MATMUL(rmat,cj)
!!!!         IF (ABS(cj(3))<d) THEN !this site lies within a disk of 2d thickness 
!!!!            IF (j<=ns .AND. jsite/=isite) success=.FALSE. !this site has already been selected 
!!!!            !note that site present in the pixel
!!!!            ix=FLOOR(cj(1)/2.5)
!!!!            iy=FLOOR(cj(2)/2.5)
!!!!            IF (ABS(ix)<=50 .AND. ABS(iy)<=50) Pixels(ix,iy)=.TRUE.
!!!!         END IF
!!!!      END DO
!!!!      area=0
!!!!      CALL AnalyzePixels(0,0,Pixels,PixelsOfInterest,area)
!!!!!      WRITE(6,*) irot,area,success
!!!!      IF (area<minarea) THEN
!!!!         srot=irot !selected rotation
!!!!         AcceptableCrossSection=success !the cross-section with lowest area
!!!!         !dont contain sites that have been already studied
!!!!         minarea=area
!!!!      END IF
!!!!   END DO
!!!!
!!!!   !selected sites should be tagged as selected 
!!!!   IF (srot>0 .AND. AcceptableCrossSection) THEN
!!!!      rmat=GetRotationMatrix(rmatrices,srot,nrot)
!!!!      WRITE(6,*) "phi:",e(1,srot)
!!!!      WRITE(6,*) "theta:",e(2,srot) 
!!!!      WRITE(6,*) "psi:",e(3,srot)
!!!!      WRITE(6,*) "Area :",minarea
!!!!      WRITE(6,*) "Atoms marked ..."
!!!!      WRITE(UNIT=6,FMT='(i8)',ADVANCE="no") isite-1
!!!!      !mark atoms and find center of mass
!!!!      nx=1 !isite already included
!!!!      com=ci
!!!!      DO j=1,n
!!!!         jsite=arr(j)
!!!!         cj=AtomCoord(3*jsite-2:3*jsite)-ci
!!!!         cj=MATMUL(rmat,cj)
!!!!         IF (ABS(cj(3))<d .AND. j>ns) THEN !mark it
!!!!            ix=FLOOR(cj(1)/2.5)
!!!!            iy=FLOOR(cj(2)/2.5)
!!!!            IF (ix<=50 .AND. iy<=50) THEN
!!!!               IF (PixelsOfInterest(ix,iy)) THEN
!!!!                  ns=ns+1
!!!!                  ksite=arr(ns)
!!!!                  arr(ns)=jsite
!!!!                  arr(j)=ksite
!!!!                  nx=nx+1
!!!!                  com=com+AtomCoord(3*jsite-2:3*jsite)
!!!!      WRITE(UNIT=6,FMT='(i8)',ADVANCE="no") jsite-1
!!!!               END IF
!!!!            END IF
!!!!         END IF
!!!!      END DO
!!!!      WRITE(6,*) ""
!!!!      com=com/REAL(nx) !center of mass
!!!!      WRITE(6,*) "Center of mass:",com
!!!!   END IF
!!!!END SUBROUTINE GetCrossSection
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE SetupRotationMatrices(r,n,e)
   IMPLICIT NONE
   INTEGER :: n,n1
   REAL(dp), DIMENSION(3,3*n), INTENT(INOUT) :: r
   REAL(dp), DIMENSION(3,n), INTENT(INOUT) :: e 
   REAL(dp) :: pi,phi,theta,psi,r1(3,3),r2(3,3)
   REAL(dp) :: rphi(3,3),cphi,sphi
   REAL (dp) :: rtheta(3,3),ctheta,stheta
   REAL(dp) :: rpsi(3,3),cpsi,spsi
   INTEGER :: iphi,itheta,ipsi,nphi,ntheta,npsi
   
   pi=4._dp*ATAN(1._dp)
   nphi=20
   ntheta=nphi/2
   npsi=nphi

   r=0._dp

   n1=0
   DO iphi=1,nphi
      phi=REAL(iphi-1)/REAL(nphi)*2.*pi
      rphi=0._dp; cphi=cos(phi); sphi=sin(phi); rphi(1,1)=cphi; rphi(2,2)=cphi;
      rphi(1,2)=sphi; rphi(2,1)=-sphi; rphi(3,3)=1._dp
      DO itheta=1,ntheta
         theta=REAL(itheta-1)/REAL(ntheta)*pi
         rtheta=0._dp; ctheta=cos(theta); stheta=sin(theta)
         rtheta(1,1)=1._dp; rtheta(2,2)=ctheta; rtheta(3,3)=ctheta;
         rtheta(2,3)=stheta; rtheta(3,2)=-stheta
         r1=MATMUL(rtheta,rphi)
         DO ipsi=1,npsi
            n1=n1+1
            IF (n1>n) THEN
               WRITE(6,*) "Rotation matrices are more than space allocated"
               STOP
            END IF
            psi=REAL(ipsi-1)/REAL(npsi)*2.*pi
            rpsi=0._dp; cpsi=cos(psi); spsi=sin(psi); rpsi(1,1)=cpsi;
            rpsi(2,2)=cpsi; rpsi(1,2)=spsi; rpsi(2,1)=-spsi; rpsi(3,3)=1._dp
            r2=MATMUL(rpsi,r1)
            r(1:3,3*n1-2)=r2(1:3,1)
            r(1:3,3*n1-1)=r2(1:3,2)
            r(1:3,3*n1)=r2(1:3,3)
            e(1:3,n1)=(/phi,theta,psi/)
         END DO
      END DO
   END DO
   
   IF (n1/=n) THEN
      WRITE(6,*) "Space allocated does not match size of rotation matrices"
      STOP
   END IF
   WRITE(6,*) "Rotation matrices are initialized .. #matrices ..",n
END SUBROUTINE SetupRotationMatrices
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE GetEulerAngle(srot,phi,theta,psi)
!definitions should match ones given in SetupRotationMatrices
   IMPLICIT NONE
   REAL :: phi,theta,psi,pi
   INTEGER :: srot,nphi,ntheta,npsi,ipsix,ithetax,iphix

   pi=4._dp*ATAN(1._dp)
   nphi=20
   ntheta=nphi/2
   npsi=nphi

   ipsix=MOD(srot-1,npsi)
   ithetax=MOD((srot-1-ipsix)/npsi,ntheta)
   iphix=(srot-1-ipsix-ithetax*npsi)/ntheta/npsi

   psi=REAL(ipsix)/REAL(npsi)*2.*pi
   theta=REAL(ithetax)/REAL(ntheta)*pi
   phi=REAL(iphix)/REAL(nphi)*2.*pi

END SUBROUTINE GetEulerAngle 
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
FUNCTION GetRotationMatrix(r,i,n)
!get ith of n rotation matrices
   IMPLICIT NONE
   INTEGER :: i,n
   REAL(dp) :: r(3,3*n),GetRotationMatrix(3,3)

   GetRotationMatrix(1:3,1)=r(1:3,3*i-2)
   GetRotationMatrix(1:3,2)=r(1:3,3*i-1)
   GetRotationMatrix(1:3,3)=r(1:3,3*i)
END FUNCTION GetRotationMatrix
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
RECURSIVE SUBROUTINE AnalyzePixels(ix,iy,parr,parr1,area)
   IMPLICIT NONE
   INTEGER :: ix,iy,area,ix1,iy1
   LOGICAL :: parr(-50:50,-50:50),parr1(-50:50,-50:50)

   IF (ABS(ix)>50 .OR. ABS(iy)>50.) RETURN
   IF (parr(ix,iy)) THEN
      parr(ix,iy)=.FALSE. !this pixel is analyzed
      parr1(ix,iy)=.TRUE. !this pixel will be of interest
      area=area+1
      ix1=ix+1; iy1=iy; CALL AnalyzePixels(ix1,iy1,parr,parr1,area)
      ix1=ix-1; iy1=iy; CALL AnalyzePixels(ix1,iy1,parr,parr1,area)
      ix1=ix; iy1=iy+1; CALL AnalyzePixels(ix1,iy1,parr,parr1,area)
      ix1=ix; iy1=iy-1; CALL AnalyzePixels(ix1,iy1,parr,parr1,area)
   END IF
END SUBROUTINE AnalyzePixels
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
