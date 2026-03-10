MODULE AlignNanoparticle
   USE VARIABLE_TYPE
   USE utilities
   USE OptimPackage
   TYPE(SystemContainer), POINTER :: gloALref,gloAL
   REAL(dp) :: cob(3)
   INTEGER :: giatom,gjatom
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   !Part 1. Code for use with MD trajectory
   SUBROUTINE Align(AL,ALref,iatom,jatom)
   !this will align frame AL to frame ALref
   !using iatom and jatom as fixed positions
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,Alref
      INTEGER :: iatom,jatom,xatom,iter,ITMAX,errorstatus
      REAL(dp) :: dr(3),theta(3),R(3,3),cob(3)
      REAL(dp) :: xtol,ftol,gtol,fret

      !translate AL to match iatom position
      CALL CenterOfMassCorrect(AL)
      cob=AL%BoxSize(1:3)/2. !center of box
!~       dr=ALref%AtomCoord(3*iatom-2:3*iatom)-AL%AtomCoord(3*iatom-2:3*iatom)
!~     !write(6,*) AL%AtomCoord(3*xatom-2:3*xatom)
!~     !write(6,*) ALref%AtomCoord(3*xatom-2:3*xatom)
!~     !  DO xatom=1,AL%NAtoms
!~     !     AL%AtomCoord(3*xatom-2:3*xatom)=AL%AtomCoord(3*xatom-2:3*xatom)+dr
!~     !  END DO

      !rotate
      theta=0.3_dp ! (angles)
      ITMAX=100  ! (maximum iterations)
      gloALref=>ALref
      gloAL=>AL
      giatom=iatom
      gjatom=jatom
      xtol=1.e-10; ftol=1.e-10; gtol=1.e-10  ! tolerance of x=distance, f=function, g=gradient
      write(6,*) "Performing optimzation.."
      CALL ConjugateGradient(p=theta,ftol=ftol,gtol=gtol,xtol=xtol,iter=iter, &
         fret=fret,func=Sep,dfunc=dSep,ITMAX=ITMAX,errorstatus=errorstatus)
      WRITE(6,*) "theta:",theta,Sep(theta,errorstatus)

      R=RotationMatrix(theta)
      DO xatom=1,AL%NAtoms
         dr=AL%AtomCoord(3*xatom-2:3*xatom)-cob
         dr=MATMUL(R,dr)
         AL%AtomCoord(3*xatom-2:3*xatom)=cob+dr
      END DO

   END SUBROUTINE Align
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION Sep(t,errorstatus)
   !Part 1. Code for use with MD trajectory
      USE VARIABLE_TYPE
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(dp), DIMENSION(:), INTENT(IN) :: t
      REAL(dp) :: Sep,R(3,3),d(3),dref(3)

      R=RotationMatrix(t)
      Sep=0._dp
      errorstatus=0
      
      !iatom
      dref=gloALref%AtomCoord(3*giatom-2:3*giatom)-cob
      dref=dref/SQRT(DOT_PRODUCT(dref,dref))
      d=gloAL%AtomCoord(3*giatom-2:3*giatom)-cob
      d=d/SQRT(DOT_PRODUCT(d,d))
      d=MATMUL(R,d)
      Sep=Sep+DOT_PRODUCT(d-dref,d-dref)
      !jatom
      dref=gloALref%AtomCoord(3*gjatom-2:3*gjatom)-cob
      dref=dref/SQRT(DOT_PRODUCT(dref,dref))
      d=gloAL%AtomCoord(3*gjatom-2:3*gjatom)-cob
      d=d/SQRT(DOT_PRODUCT(d,d))
      d=MATMUL(R,d)
      Sep=Sep+DOT_PRODUCT(d-dref,d-dref)
   END FUNCTION Sep
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION dSep(t,errorstatus)
   !Part 1. Code for use with MD trajectory
      USE VARIABLE_TYPE
      IMPLICIT NONE
      INTEGER :: i,errorstatus
      REAL(dp), DIMENSION(:), INTENT(IN) :: t
      REAL(dp), DIMENSION(size(t)) :: dSep
      REAL(dp), DIMENSION(3) :: tx
      REAL(dp) :: fp,fn

      errorstatus=0
      DO i=1,3
         tx=t; tx(i)=tx(i)+0.01; fp=Sep(tx,errorstatus)
         tx=t; tx(i)=tx(i)-0.01; fn=Sep(tx,errorstatus)
         dSep(i)=(fp-fn)/0.02
      END DO
   END FUNCTION 
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CenterOfMassCorrect(AL)
   !Part 1. Code for use with MD trajectory
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp) :: com(3),dr(3)
      INTEGER :: iatom

      com=0._dp
      DO iatom=1,AL%NAtoms
         com=com+AL%AtomCoord(3*iatom-2:3*iatom)
      END DO
      com=com/REAL(AL%NAtoms)
      cob=AL%BoxSize(1:3)/2. !cneter of box
      dr=cob-com
      DO iatom=1,AL%NAtoms
         AL%AtomCoord(3*iatom-2:3*iatom)=AL%AtomCoord(3*iatom-2:3*iatom)+dr
      END DO
   END SUBROUTINE CenterOfMassCorrect 
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckForSubstantialMovement(AL,ALref,SubstantialMovementFound)
   !Part 1. Code for use with MD trajectory
   !both AL and ALref have been aligned already
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,ALref
      INTEGER :: iatom
      LOGICAL, DIMENSION(:) :: SubstantialMovementFound
      REAL(dp) :: tol=1.,dr(3)
      
      DO iatom=1,AL%NAtoms
         dr=AL%AtomCoord(3*iatom-2:3*iatom)-ALref%AtomCoord(3*iatom-2:3*iatom) 
         IF (SQRT(DOT_PRODUCT(dr,dr))>tol) THEN
            SubstantialMovementFound(iatom)=.TRUE.
         END IF
      END DO
   END SUBROUTINE CheckForSubstantialMovement
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   

   !Part 2. Code for latching NP with large slab
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LatchParticle(AL,ALref,boxsize_x,boxsize_y,boxsize_z)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,ALref
      INTEGER :: i,j,pos,ncells(3),ncollisions
      REAL(dp) :: CellSize(3),dr(3),cm(3), boxsize_x,boxsize_y,boxsize_z
      INTEGER, PARAMETER :: n=14 !number of lattice sites
      REAL(dp), PARAMETER :: a=4.0801
      REAL(dp), DIMENSION(3*n) :: unitcellpositions=a*(/ &
      0.,0.,0.,   1.,0.,0.,   0.,1.,0.,   0.,0.,1.,   1.,1.,0.,  1.,0.,1., 0.,1.,1., &
      1.,1.,1.,   0.5,0.5,0.,   0.5,0.,0.5,   0.,0.5,0.5,    1.,0.5,0.5, 0.5,1.,0.5, &
      0.5,0.5,1./)
      REAL(dp) :: cellcoord(3),pos1(3)
      REAL(dp) :: NearestSlabCoord(3),coord(3),dist(n),cb(3)
      INTEGER :: cell(3),species(1),overlap_fnn,overlap_snn,overlap_tnn,overlap_frnn

      gloALref=>ALref
      gloAL=>AL

      cm=COM(AL)
      !cb=(/126.48310000,126.48310000,126.48310000/)   !COM(ALref)
      !cb=(/200.0000,200.0000,200.0000/)
	cb(1)=(boxsize_x*4.0801)/2.0_dp
	cb(2)=(boxsize_y*4.0801)/2.0_dp
	cb(3)=(boxsize_z*4.0801)/2.0_dp
    WRITE(6,*) "cb3:",cb
      coord=cm
      cell=MAX(1,CEILING(coord/a)) !cell indices
      cm=REAL(cell-1)*a !location of cell corner
      dr=cb-cm

      !make AL coordinates lattice values
      DO i=1,AL%NAtoms
         coord=AL%AtomCoord(3*i-2:3*i)
         coord=coord+dr
         cell=MAX(1,CEILING(coord/a)) !cell indices
         cellcoord=REAL(cell-1)*a !location of cell corner
         DO j=1,n
            pos1=unitcellpositions(3*j-2:3*j)+cellcoord
            pos1=coord-pos1
            dist(j)=DOT_PRODUCT(pos1,pos1)
         END DO
         j=MINLOC(dist,1)
         NearestSlabCoord=unitcellpositions(3*j-2:3*j)+cellcoord
         AL%AtomCoord(3*i-2:3*i)=NearestSlabCoord !closest match to slab coordinates
      END DO
      WRITE(6,*) "Updated NP coordinate positions to best match slab"

      !read species into slab
      ncollisions=0
      ALref%AtomSpecies=0
      cellsize=ALref%LL%CellSize
      ncells=ALref%LL%NumberCell3D

      OPEN(UNIT=236,FILE="overlaps")
      OPEN(UNIT=237,FILE="overlap_species")

      DO i=1,AL%NAtoms
         coord=AL%AtomCoord(3*i-2:3*i)
         species=AL%AtomSpecies(i)
        !c=NearestSlabCoord(c)
         pos=GetIndex(ALref,coord,cellsize,ncells) !identify site with coordinate coord
         IF (ALref%AtomSpecies(pos)/=0) THEN
             ncollisions=ncollisions+1
             WRITE(236,*) coord
             WRITE(237,*) species
         ELSE 
             ALref%AtomSpecies(pos)=AL%AtomSpecies(i)
         END IF
      END DO
      WRITE(6,*) "Number of collisions while matching:",ncollisions
      CLOSE(236)
      CLOSE(237)

      CALL Latch_to_firstnearestneighbours(ALref,ncollisions,overlap_fnn)
      CALL Latch_to_secondnearestneighbours(ALref,overlap_fnn,overlap_snn)
!      CALL Latch_to_thirdnearestneighbours(ALref,overlap_snn,overlap_tnn)
!      CALL Latch_to_fourthnearestneighbours(ALref,overlap_tnn,overlap_frnn)

!      CALL Fill_vacancies_with_CN_12(ALref)
!
!      CALL Fill_vacancies_with_CN_11(ALref)
!      CALL Fill_vacancies_with_CN_12(ALref)
!
!      CALL Fill_vacancies_with_CN_10(ALref)
!      CALL Fill_vacancies_with_CN_11(ALref)
!      CALL Fill_vacancies_with_CN_12(ALref)
!
!      CALL Fill_vacancies_with_CN_9(ALref)
!      CALL Fill_vacancies_with_CN_10(ALref)
!      CALL Fill_vacancies_with_CN_11(ALref)
!      CALL Fill_vacancies_with_CN_12(ALref)

!      CALL Fill_vacancies_with_CN_8(ALref)
!      CALL Fill_vacancies_with_CN_9(ALref)
!      CALL Fill_vacancies_with_CN_10(ALref)
!      CALL Fill_vacancies_with_CN_11(ALref)
!      CALL Fill_vacancies_with_CN_12(ALref)

   END SUBROUTINE LatchParticle
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetIndex(ALref,coord,cellsize,ncells)
   !retruns the atom index in AL with coordinate closest to coord
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: ALref
      REAL(dp) :: coord(3),cellsize(3),dmin
      INTEGER :: icell(3)
      INTEGER :: GetIndex,ncells(3),v(3),ix,iy,iz,LLFactor(3)
      INTEGER :: MaxAtomPerCell,xatom
     
      GetIndex=0
      icell=MAX(1,CEILING(coord/cellsize))
      LLFactor(1)=ncells(2)*ncells(3) !x
      LLFactor(2)=ncells(3) !y
      LLFactor(3)=1 !z
      MaxAtomPerCell=ALref%LL%MaxAtomPerCell
      dmin=1000._dp
      xatom=0
      DO ix=-1,1
         DO iy=-1,1
            DO iz=-1,1
               v=(/ix,iy,iz/); CALL GI(v)
            END DO
         END DO
      END DO
      GetIndex=xatom
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE GI(v)
      !search in a cell with location icell+v
         IMPLICIT NONE
         INTEGER :: v(3),jcell(3),cellindex,i,pos,jatom
         REAL(dp) :: coord1(3),distance,coordx(3)

         jcell=icell+v
         WHERE (jcell==0) jcell=ncells
         WHERE (jcell==ncells+1) jcell=1
         jcell(1:2)=jcell(1:2)-1 !convert to scalar number
         cellindex=SUM(jcell*LLFactor)

         !search cell
         DO i=1,ALref%LL%ListRange(cellindex)
            pos=(cellindex-1)*MaxAtomPerCell+i
            jatom=ALref%LL%List(pos)
            coord1=ALref%AtomCoord(3*jatom-2:3*jatom)
            coordx=PBCdistance3(ALref,coord1,coord) !get relative distance
            distance=SQRT(DOT_PRODUCT(coordx,coordx))
            IF (distance<dmin) THEN
               xatom=jatom
               dmin=distance
            END IF
         END DO
      END SUBROUTINE GI
      !xxxxxxxxxxxxxxxxxxxxxxx
   END FUNCTION GetIndex
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION LatchResidual(t,errorstatus)
      USE VARIABLE_TYPE
      IMPLICIT NONE
      INTEGER :: errorstatus,i
      REAL(dp), DIMENSION(:), INTENT(IN) :: t
      REAL(dp) :: LatchResidual,R(3,3),d(3),cm(3)
      REAL(dp) :: translate(3),theta(3),c(3)

      errorstatus=0
      LatchResidual=0._dp 

      translate=t(1:3)
      theta=t(4:6)
      cm=COM(gloAL) !center of mass without translation

      !rotate coordinates about com, also apply translation
      R=RotationMatrix(theta)
      DO i=1,gloAL%NAtoms
         c=gloAL%AtomCoord(3*i-2:3*i)-cm
         c=MATMUL(R,c)+cm !rotate the coord
         c=c+translate !add translation vector
         d=NearestSlabCoord(c) !find nearest slab coord
         d=d-c
         LatchResidual=LatchResidual+DOT_PRODUCT(d,d)
      END DO
   END FUNCTION LatchResidual
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION dLatchResidual(t,errorstatus)
      USE VARIABLE_TYPE
      IMPLICIT NONE
      INTEGER :: errorstatus,i
      REAL(dp), DIMENSION(:), INTENT(IN) :: t
      REAL(dp) :: dLatchResidual(6),tx(6),fp,fn

      dLatchResidual=0._dp
      errorstatus=0
      DO i=1,6
         tx=t; tx(i)=tx(i)+0.01; fp=LatchResidual(tx,errorstatus)
         tx=t; tx(i)=tx(i)-0.01; fn=LatchResidual(tx,errorstatus)
         dLatchResidual(i)=(fp-fn)/0.02
      END DO
   END FUNCTION dLatchResidual
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION COM(AL)
   !returns center of mass of particle
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp) :: com(3),dr(3)
      INTEGER :: iatom
      
      com=0._dp
      DO iatom=1,AL%NAtoms
         com=com+AL%AtomCoord(3*iatom-2:3*iatom)
      END DO
      com=com/REAL(AL%NAtoms)
   END FUNCTION COM
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION NearestSlabCoord(coord)
      IMPLICIT NONE
      INTEGER, PARAMETER :: n=14 !number of lattice sites
      REAL(dp), PARAMETER :: a=4.0801
      REAL(dp), DIMENSION(3*n) :: unitcellpositions=a*(/ &
      0.,0.,0.,   1.,0.,0.,   0.,1.,0.,   0.,0.,1.,   1.,1.,0.,  1.,0.,1.,   0.,1.,1., &
      1.,1.,1.,   0.5,0.5,0.,   0.5,0.,0.5,   0.,0.5,0.5,    1.,0.5,0.5,    0.5,1.,0.5, &
      0.5,0.5,1./)
      REAL(dp) :: cellcoord(3),pos(3)
      REAL(dp) :: NearestSlabCoord(3),coord(3),dist(n)
      INTEGER :: cell(3),i
      
      cell=MAX(1,CEILING(coord/a)) !cell indices
      cellcoord=REAL(cell-1)*a !location of cell corner
      DO i=1,n
         pos=unitcellpositions(3*i-2:3*i)+cellcoord
         pos=coord-pos
         dist(i)=DOT_PRODUCT(pos,pos)
      END DO
      i=MINLOC(dist,1)
      NearestSlabCoord=unitcellpositions(3*i-2:3*i)+cellcoord
   END FUNCTION NearestSlabCoord
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Fill_vacancies_with_CN_12(ALref)
   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: ALref
   INTEGER :: i,k,pos,nPt,nFe,ncells(3),occupied_sites,coordination
   REAL(dp) :: CellSize(3),coord(3),fnn(3)
   INTEGER, PARAMETER :: nn=12 !first nearest neighbours
   REAL(dp), PARAMETER :: a=4.0801
   REAL(dp), DIMENSION(3*nn) :: firstnearestneighbours=a*(/ &
   0.5,0.5,0.,  -0.5,-0.5,0.,  -0.5,0.5,0.,  0.5,-0.5,0., &
   0.,0.5,0.5,  0.,-0.5,-0.5,  0.,-0.5,0.5,  0.,0.5,-0.5, &
   0.5,0.,0.5,  -0.5,0.,-0.5,  -0.5,0.,0.5, 0.5,0.,-0.5/)

   cellsize=ALref%LL%CellSize
   ncells=ALref%LL%NumberCell3D

   occupied_sites=0
   DO i=1,ALref%NAtoms
      IF (ALref%AtomSpecies(i)==0) THEN
         coord=ALref%AtomCoord(3*i-2:3*i)
         IF (coord(1)>10. .AND. coord(1)<215. .AND. coord(2)>10. .AND. coord(2)<215. .AND. coord(3)>10. .AND. coord(3)<215.) THEN
             coordination=0
             nPt=0
             nFe=0
             DO k=1,12
                fnn=firstnearestneighbours(3*k-2:3*k)+coord
                pos=GetIndex(ALref,fnn,cellsize,ncells)
                IF (ALref%AtomSpecies(pos)/=0) coordination=coordination+1
                IF (ALref%AtomSpecies(pos)==1) nFe=nFe+1
                IF (ALref%AtomSpecies(pos)==2) nPt=nPt+1
             END DO
             IF (coordination==12 .AND. nFe<=10) ALref%AtomSpecies(i)=1
             IF (coordination==12 .AND. nFe>10) ALref%AtomSpecies(i)=2
             IF (coordination==12) occupied_sites=occupied_sites+1
         END IF
      END IF
   END DO
   WRITE(6,*) occupied_sites

   END SUBROUTINE
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Fill_vacancies_with_CN_11(ALref)
   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: ALref
   INTEGER :: i,k,pos,nPt,nFe,ncells(3),occupied_sites,coordination
   REAL(dp) :: CellSize(3),coord(3),fnn(3)
   INTEGER, PARAMETER :: nn=12 !first nearest neighbours
   REAL(dp), PARAMETER :: a=4.0801
   REAL(dp), DIMENSION(3*nn) :: firstnearestneighbours=a*(/ &
   0.5,0.5,0.,  -0.5,-0.5,0.,  -0.5,0.5,0.,  0.5,-0.5,0., &
   0.,0.5,0.5,  0.,-0.5,-0.5,  0.,-0.5,0.5,  0.,0.5,-0.5, &
   0.5,0.,0.5,  -0.5,0.,-0.5,  -0.5,0.,0.5, 0.5,0.,-0.5/)

   cellsize=ALref%LL%CellSize
   ncells=ALref%LL%NumberCell3D

   occupied_sites=0
   DO i=1,ALref%NAtoms
      IF (ALref%AtomSpecies(i)==0) THEN
         coord=ALref%AtomCoord(3*i-2:3*i)
         IF (coord(1)>10. .AND. coord(1)<215. .AND. coord(2)>10. .AND. coord(2)<215. .AND. coord(3)>10. .AND. coord(3)<215.) THEN
             coordination=0
             nPt=0
             nFe=0
             DO k=1,12
                fnn=firstnearestneighbours(3*k-2:3*k)+coord
                pos=GetIndex(ALref,fnn,cellsize,ncells)
                IF (ALref%AtomSpecies(pos)/=0) coordination=coordination+1
                IF (ALref%AtomSpecies(pos)==1) nFe=nFe+1
                IF (ALref%AtomSpecies(pos)==2) nPt=nPt+1
             END DO
             IF (coordination==11 .AND. nFe<=9) ALref%AtomSpecies(i)=1
             IF (coordination==11 .AND. nFe>9) ALref%AtomSpecies(i)=2
             IF (coordination==11) occupied_sites=occupied_sites+1
         END IF
      END IF
   END DO
   WRITE(6,*) occupied_sites

   END SUBROUTINE
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Fill_vacancies_with_CN_10(ALref)
   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: ALref
   INTEGER :: i,k,pos,nPt,nFe,ncells(3),occupied_sites,coordination
   REAL(dp) :: CellSize(3),coord(3),fnn(3)
   INTEGER, PARAMETER :: nn=12 !first nearest neighbours
   REAL(dp), PARAMETER :: a=4.0801
   REAL(dp), DIMENSION(3*nn) :: firstnearestneighbours=a*(/ &
   0.5,0.5,0.,  -0.5,-0.5,0.,  -0.5,0.5,0.,  0.5,-0.5,0., &
   0.,0.5,0.5,  0.,-0.5,-0.5,  0.,-0.5,0.5,  0.,0.5,-0.5, &
   0.5,0.,0.5,  -0.5,0.,-0.5,  -0.5,0.,0.5, 0.5,0.,-0.5/)

   cellsize=ALref%LL%CellSize
   ncells=ALref%LL%NumberCell3D

   occupied_sites=0
   DO i=1,ALref%NAtoms
      IF (ALref%AtomSpecies(i)==0) THEN
         coord=ALref%AtomCoord(3*i-2:3*i)
         IF (coord(1)>10. .AND. coord(1)<215. .AND. coord(2)>10. .AND. coord(2)<215. .AND. coord(3)>10. .AND. coord(3)<215.) THEN
             coordination=0
             nPt=0
             nFe=0
             DO k=1,12
                fnn=firstnearestneighbours(3*k-2:3*k)+coord
                pos=GetIndex(ALref,fnn,cellsize,ncells)
                IF (ALref%AtomSpecies(pos)/=0) coordination=coordination+1
                IF (ALref%AtomSpecies(pos)==1) nFe=nFe+1
                IF (ALref%AtomSpecies(pos)==2) nPt=nPt+1
             END DO
             IF (coordination==10 .AND. nFe<=8) ALref%AtomSpecies(i)=1
             IF (coordination==10 .AND. nFe>8) ALref%AtomSpecies(i)=2
             IF (coordination==10) occupied_sites=occupied_sites+1
         END IF
      END IF
   END DO
   WRITE(6,*) occupied_sites

   END SUBROUTINE
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Fill_vacancies_with_CN_9(ALref)
   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: ALref
   INTEGER :: i,k,pos,nPt,nFe,ncells(3),occupied_sites,coordination
   REAL(dp) :: CellSize(3),coord(3),fnn(3)
   INTEGER, PARAMETER :: nn=12 !first nearest neighbours
   REAL(dp), PARAMETER :: a=4.0801
   REAL(dp), DIMENSION(3*nn) :: firstnearestneighbours=a*(/ &
   0.5,0.5,0.,  -0.5,-0.5,0.,  -0.5,0.5,0.,  0.5,-0.5,0., &
   0.,0.5,0.5,  0.,-0.5,-0.5,  0.,-0.5,0.5,  0.,0.5,-0.5, &
   0.5,0.,0.5,  -0.5,0.,-0.5,  -0.5,0.,0.5, 0.5,0.,-0.5/)

   cellsize=ALref%LL%CellSize
   ncells=ALref%LL%NumberCell3D

   occupied_sites=0
   DO i=1,ALref%NAtoms
      IF (ALref%AtomSpecies(i)==0) THEN
         coord=ALref%AtomCoord(3*i-2:3*i)
         IF (coord(1)>10. .AND. coord(1)<215. .AND. coord(2)>10. .AND. coord(2)<215. .AND. coord(3)>10. .AND. coord(3)<215.) THEN
             coordination=0
             nPt=0
             nFe=0
             DO k=1,12
                fnn=firstnearestneighbours(3*k-2:3*k)+coord
                pos=GetIndex(ALref,fnn,cellsize,ncells)
                IF (ALref%AtomSpecies(pos)/=0) coordination=coordination+1
                IF (ALref%AtomSpecies(pos)==1) nFe=nFe+1
                IF (ALref%AtomSpecies(pos)==2) nPt=nPt+1
             END DO
             IF (coordination==9 .AND. nFe<=7) ALref%AtomSpecies(i)=1
             IF (coordination==9 .AND. nFe>7) ALref%AtomSpecies(i)=2
             IF (coordination==9) occupied_sites=occupied_sites+1
         END IF
      END IF
   END DO
   WRITE(6,*) occupied_sites

   END SUBROUTINE
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Fill_vacancies_with_CN_8(ALref)
   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: ALref
   INTEGER :: i,k,pos,nPt,nFe,ncells(3),occupied_sites,coordination
   REAL(dp) :: CellSize(3),coord(3),fnn(3)
   INTEGER, PARAMETER :: nn=12 !first nearest neighbours
   REAL(dp), PARAMETER :: a=4.0801
   REAL(dp), DIMENSION(3*nn) :: firstnearestneighbours=a*(/ &
   0.5,0.5,0.,  -0.5,-0.5,0.,  -0.5,0.5,0.,  0.5,-0.5,0., &
   0.,0.5,0.5,  0.,-0.5,-0.5,  0.,-0.5,0.5,  0.,0.5,-0.5, &
   0.5,0.,0.5,  -0.5,0.,-0.5,  -0.5,0.,0.5, 0.5,0.,-0.5/)

   cellsize=ALref%LL%CellSize
   ncells=ALref%LL%NumberCell3D

   occupied_sites=0
   DO i=1,ALref%NAtoms
      IF (ALref%AtomSpecies(i)==0) THEN
         coord=ALref%AtomCoord(3*i-2:3*i)
         IF (coord(1)>10. .AND. coord(1)<215. .AND. coord(2)>10. .AND. coord(2)<215. .AND. coord(3)>10. .AND. coord(3)<215.) THEN
             coordination=0
             nPt=0
             nFe=0
             DO k=1,12
                fnn=firstnearestneighbours(3*k-2:3*k)+coord
                pos=GetIndex(ALref,fnn,cellsize,ncells)
                IF (ALref%AtomSpecies(pos)/=0) coordination=coordination+1
                IF (ALref%AtomSpecies(pos)==1) nFe=nFe+1
                IF (ALref%AtomSpecies(pos)==2) nPt=nPt+1
             END DO
             IF (coordination==8 .AND. nFe<=6) ALref%AtomSpecies(i)=1
             IF (coordination==8 .AND. nFe>6) ALref%AtomSpecies(i)=2
             IF (coordination==8) occupied_sites=occupied_sites+1
         END IF
      END IF
   END DO
   WRITE(6,*) occupied_sites

   END SUBROUTINE
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Latch_to_firstnearestneighbours(ALref,ncollisions,overlap_fnn)
   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: ALref
   INTEGER :: i,j,k,pos,ncollisions,species,cn(12),ncells(3),coordination,overlap_fnn
   REAL(dp) :: CellSize(3),coord(3),fnn(3),fnnv(3)
   INTEGER, PARAMETER :: nn=12 !first nearest neighbours
   REAL(dp), PARAMETER :: a=4.0801
   REAL(dp), DIMENSION(3*nn) :: firstnearestneighbours=a*(/ &
   0.5,0.5,0.,  -0.5,-0.5,0.,  -0.5,0.5,0.,  0.5,-0.5,0., &
   0.,0.5,0.5,  0.,-0.5,-0.5,  0.,-0.5,0.5,  0.,0.5,-0.5, &
   0.5,0.,0.5,  -0.5,0.,-0.5,  -0.5,0.,0.5, 0.5,0.,-0.5/)
   REAL(dp),ALLOCATABLE,DIMENSION(:,:) :: overlap_coord
   REAL(dp),ALLOCATABLE,DIMENSION(:) :: overlap_coords,overlap_species


   cellsize=ALref%LL%CellSize
   ncells=ALref%LL%NumberCell3D

   ALLOCATE(overlap_coord(ncollisions,3))
   ALLOCATE(overlap_species(ncollisions))

   OPEN(UNIT=236,FILE="overlaps")
   OPEN(UNIT=237,FILE="overlap_species")

   DO i=1,ncollisions
      READ(236,*) (overlap_coord(i,j),j =1,3)   !these are AL coords
      READ(237,*) overlap_species(i)
   END DO
   CLOSE(236)
   CLOSE(237)

   ALLOCATE(overlap_coords(3*ncollisions))
   DO i=1,ncollisions
      DO j=1,3
         overlap_coords(3*i-3+j)=overlap_coord(i,j) !these are AL coords
      END DO
   END DO

   OPEN(UNIT=239,FILE="fnn_overlaps")
   OPEN(UNIT=240,FILE="fnn_overlap_species")

   overlap_fnn=0
   DO i=1,ncollisions
      coord=overlap_coords(3*i-2:3*i)  !these are AL coords
      species=overlap_species(i)
      !WRITE(6,*) coord
      DO j=1,12
         fnn=firstnearestneighbours(3*j-2:3*j)+coord
         pos=GetIndex(ALref,fnn,cellsize,ncells) !identify site withcoordinate fnn
         IF (ALref%AtomSpecies(pos)==0) THEN ! find coordination number of this site
            coordination=0
            DO k=1,12
               fnnv=firstnearestneighbours(3*k-2:3*k)+fnn
               pos=GetIndex(ALref,fnnv,cellsize,ncells)
               IF (ALref%AtomSpecies(pos)/=0) THEN
                  coordination=coordination+1
                  !WRITE(6,*) k
               END IF
            END DO
            cn(j)=coordination
            !WRITE(6,*) j
            !WRITE(6,*) coordination
         ELSE
            cn(j)=0
         END IF
      END DO
      j=MAXLOC(cn,1)
      !WRITE(6,*) j
      fnn=firstnearestneighbours(3*j-2:3*j)+coord
      pos=GetIndex(ALref,fnn,cellsize,ncells) !identify site with coordinate fnn
      IF (ALref%AtomSpecies(pos)/=0) THEN
          overlap_fnn=overlap_fnn+1
          WRITE(239,*) coord  ! this should be coord, earlier it was fnn
          WRITE(240,*) species
      ELSE
          ALref%AtomSpecies(pos)=overlap_species(i)
      END IF
   END DO
   WRITE(6,*) "Number of collisions among fnn while matching:",overlap_fnn
   CLOSE(239)
   CLOSE(240)

   END SUBROUTINE
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Latch_to_secondnearestneighbours(ALref,overlap_fnn,overlap_snn) 
   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: ALref
   INTEGER :: i,j,k,l,pos,ncollisions,species,cn(12),ncells(3),scn(6)
   INTEGER :: coordination,overlap_fnn,overlap_snn
   REAL(dp) :: CellSize(3),coord(3),fnn(3),fnnv(3),snn(3),snnv(3)
   INTEGER, PARAMETER :: nn=12 !first nearest neighbours
   REAL(dp), PARAMETER :: a=4.0801
   REAL(dp), DIMENSION(3*nn) :: firstnearestneighbours=a*(/ &
   0.5,0.5,0.,  -0.5,-0.5,0.,  -0.5,0.5,0.,  0.5,-0.5,0., &
   0.,0.5,0.5,  0.,-0.5,-0.5,  0.,-0.5,0.5,  0.,0.5,-0.5, &
   0.5,0.,0.5,  -0.5,0.,-0.5,  -0.5,0.,0.5, 0.5,0.,-0.5/)
   REAL(dp), DIMENSION(3*6) :: secondnearestneighbours=a*(/ &
   1.,0.,0.,  -1.,0.,0.,  0.,1.,0.,  0.,-1.,0.,  0.,0.,1., 0.,0.,-1./)
   REAL(dp),ALLOCATABLE,DIMENSION(:,:) :: overlap_coord,fnn_overlap_coord
   REAL(dp),ALLOCATABLE,DIMENSION(:) :: overlap_coords,overlap_species,fnn_overlap_coords,fnn_overlap_species


   cellsize=ALref%LL%CellSize
   ncells=ALref%LL%NumberCell3D

   ALLOCATE(fnn_overlap_coord(overlap_fnn,3))
   ALLOCATE(fnn_overlap_species(overlap_fnn))

   OPEN(UNIT=239,FILE="fnn_overlaps")
   OPEN(UNIT=240,FILE="fnn_overlap_species")

   DO i=1,overlap_fnn
      READ(239,*) (fnn_overlap_coord(i,j),j =1,3)
      READ(240,*) fnn_overlap_species(i)
   END DO
   CLOSE(239)
   CLOSE(240)

   ALLOCATE(fnn_overlap_coords(3*overlap_fnn))
   DO i=1,overlap_fnn
      DO j=1,3
         fnn_overlap_coords(3*i-3+j)=fnn_overlap_coord(i,j)
      END DO
   END DO

   OPEN(UNIT=242,FILE="snn_overlaps")
   OPEN(UNIT=243,FILE="snn_overlap_species")

   overlap_snn=0
   DO i=1,overlap_fnn
      coord=fnn_overlap_coords(3*i-2:3*i)
      species=fnn_overlap_species(i)
      !WRITE(6,*) coord
      DO l=1,6
         snn=secondnearestneighbours(3*l-2:3*l)+coord
         pos=GetIndex(ALref,snn,cellsize,ncells) !identify site with coordinate fnn
         IF (ALref%AtomSpecies(pos)==0) THEN ! find coordination number of this site
            coordination=0
            DO k=1,12
               snnv=firstnearestneighbours(3*k-2:3*k)+snn
               pos=GetIndex(ALref,snnv,cellsize,ncells)
               IF (ALref%AtomSpecies(pos)/=0) THEN
                  coordination=coordination+1
                  !WRITE(6,*) k
               END IF
            END DO
            scn(l)=coordination
            !WRITE(6,*) l
            !WRITE(6,*) coordination
         ELSE
            scn(l)=0
         END IF
      END DO
      l=MAXLOC(scn,1)
      !WRITE(6,*) l
      snn=secondnearestneighbours(3*l-2:3*l)+coord
      pos=GetIndex(ALref,snn,cellsize,ncells) !identify site with coordinate snn
      IF (ALref%AtomSpecies(pos)/=0) THEN
          overlap_snn=overlap_snn+1
          WRITE(242,*) coord ! earlier this was snn
          WRITE(243,*) species
      ELSE
          ALref%AtomSpecies(pos)=fnn_overlap_species(i)
      END IF
   END DO
   WRITE(6,*) "Number of collisions among snn while matching:",overlap_snn
   CLOSE(242)
   CLOSE(243)
                                                                
   END SUBROUTINE
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Latch_to_thirdnearestneighbours(ALref,overlap_snn,overlap_tnn)
   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: ALref
   INTEGER :: i,j,k,l,pos,ncollisions,species,ncells(3),tcn(24)
   INTEGER :: coordination,overlap_snn,overlap_tnn
   REAL(dp) :: CellSize(3),coord(3),tnn(3),tnnv(3)
   INTEGER, PARAMETER :: nn=12 !first nearest neighbours
   REAL(dp), PARAMETER :: a=4.0801
   REAL(dp), DIMENSION(3*nn) :: firstnearestneighbours=a*(/ &
   0.5,0.5,0.,  -0.5,-0.5,0.,  -0.5,0.5,0.,  0.5,-0.5,0., &
   0.,0.5,0.5,  0.,-0.5,-0.5,  0.,-0.5,0.5,  0.,0.5,-0.5, &
   0.5,0.,0.5,  -0.5,0.,-0.5,  -0.5,0.,0.5, 0.5,0.,-0.5/)
   REAL(dp), DIMENSION(3*24) :: thirdnearestneighbours=a*(/ &
   0.5,0.5,1., 0.5,-0.5,1., -0.5,0.5,1., -0.5,-0.5,1., 0.5,0.5,-1.,0.5,-0.5,-1., -0.5,0.5,-1., -0.5,-0.5,-1., &
   0.5,1.,0.5, 0.5,1.,-0.5, -0.5,1.,0.5, -0.5,1.,-0.5, 0.5,-1.,0.5,0.5,-1.,-0.5, -0.5,-1.,0.5, -0.5,-1.,-0.5, &
   1.,0.5,0.5, 1.,0.5,-0.5, 1.,-0.5,0.5, 1.,-0.5,-0.5, -1.,0.5,0.5,-1.,0.5,-0.5, -1.,-0.5,0.5, -1.,-0.5,-0.5/)
   REAL(dp),ALLOCATABLE,DIMENSION(:,:) :: snn_overlap_coord
   REAL(dp),ALLOCATABLE,DIMENSION(:) :: snn_overlap_coords,snn_overlap_species

   cellsize=ALref%LL%CellSize
   ncells=ALref%LL%NumberCell3D

   ALLOCATE(snn_overlap_coord(overlap_snn,3))
   ALLOCATE(snn_overlap_species(overlap_snn))

   OPEN(UNIT=240,FILE="snn_overlaps")
   OPEN(UNIT=241,FILE="snn_overlap_species")

   DO i=1,overlap_snn
      READ(240,*) (snn_overlap_coord(i,j),j =1,3)
      READ(241,*) snn_overlap_species(i)
   END DO
   CLOSE(240)
   CLOSE(241)

   ALLOCATE(snn_overlap_coords(3*overlap_snn))
   DO i=1,overlap_snn
      DO j=1,3
         snn_overlap_coords(3*i-3+j)=snn_overlap_coord(i,j)
      END DO
   END DO

   OPEN(UNIT=242,FILE="tnn_overlaps")
   OPEN(UNIT=243,FILE="tnn_overlap_species")

   overlap_tnn=0
   DO i=1,overlap_snn
      coord=snn_overlap_coords(3*i-2:3*i)
      species=snn_overlap_species(i)
      !WRITE(6,*) coord
      !pos=GetIndex(ALref,coord,cellsize,ncells) !identify site with coordinate coord
      DO l=1,24
         tnn=thirdnearestneighbours(3*l-2:3*l)+coord
         pos=GetIndex(ALref,tnn,cellsize,ncells) !identify site with coordinate fnn
         IF (ALref%AtomSpecies(pos)==0) THEN ! find coordination number of this site
            coordination=0
            DO k=1,12
               tnnv=firstnearestneighbours(3*k-2:3*k)+tnn
               pos=GetIndex(ALref,tnnv,cellsize,ncells)
               IF (ALref%AtomSpecies(pos)/=0) THEN
                  coordination=coordination+1
                  !WRITE(6,*) k
               END IF
            END DO
            tcn(l)=coordination
            !WRITE(6,*) l
            !WRITE(6,*) coordination
         ELSE
            tcn(l)=0
         END IF
      END DO
      l=MAXLOC(tcn,1)
      !WRITE(6,*) l
      tnn=thirdnearestneighbours(3*l-2:3*l)+coord
      pos=GetIndex(ALref,tnn,cellsize,ncells) !identify site with coordinate snn
      IF (ALref%AtomSpecies(pos)/=0) THEN
          overlap_tnn=overlap_tnn+1
          WRITE(242,*) coord !earlier this was tnn
          WRITE(243,*) species
      ELSE
          ALref%AtomSpecies(pos)=snn_overlap_species(i)
      END IF
   END DO
   WRITE(6,*) "Number of collisions among tnn while matching:",overlap_tnn
   CLOSE(242)
   CLOSE(243)

   END SUBROUTINE
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Latch_to_fourthnearestneighbours(ALref,overlap_tnn,overlap_frnn)
   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: ALref
   INTEGER :: i,j,k,l,pos,ncollisions,species,ncells(3),frcn(12)
   INTEGER :: coordination,overlap_tnn,overlap_frnn
   REAL(dp) :: CellSize(3),coord(3),frnn(3),frnnv(3)
   INTEGER, PARAMETER :: nn=12 !first nearest neighbours
   REAL(dp), PARAMETER :: a=4.0801
   REAL(dp), DIMENSION(3*nn) :: firstnearestneighbours=a*(/ &
   0.5,0.5,0.,  -0.5,-0.5,0.,  -0.5,0.5,0.,  0.5,-0.5,0., &
   0.,0.5,0.5,  0.,-0.5,-0.5,  0.,-0.5,0.5,  0.,0.5,-0.5, &
   0.5,0.,0.5,  -0.5,0.,-0.5,  -0.5,0.,0.5, 0.5,0.,-0.5/)
   REAL(dp), DIMENSION(3*12) :: fourthnearestneighbours=a*(/ &
   1.,1.,0., 1.,0.,1., 0.,1.,1., -1.,-1.,0., -1.,0.,-1., 0.,-1.,-1., &
   1.,-1.,0., 1.,0.,-1., 0.,1.,-1., -1.,1.,0., -1.,0.,1., 0.,-1.,1./)
   REAL(dp),ALLOCATABLE,DIMENSION(:,:) :: tnn_overlap_coord
   REAL(dp),ALLOCATABLE,DIMENSION(:) :: tnn_overlap_coords,tnn_overlap_species

   cellsize=ALref%LL%CellSize
   ncells=ALref%LL%NumberCell3D

   ALLOCATE(tnn_overlap_coord(overlap_tnn,3))
   ALLOCATE(tnn_overlap_species(overlap_tnn))

   OPEN(UNIT=242,FILE="tnn_overlaps")
   OPEN(UNIT=243,FILE="tnn_overlap_species")

   DO i=1,overlap_tnn
      READ(242,*) (tnn_overlap_coord(i,j),j =1,3)
      READ(243,*) tnn_overlap_species(i)
   END DO
   CLOSE(242)
   CLOSE(243)

   ALLOCATE(tnn_overlap_coords(3*overlap_tnn))
   DO i=1,overlap_tnn
      DO j=1,3
         tnn_overlap_coords(3*i-3+j)=tnn_overlap_coord(i,j)
      END DO
   END DO

   OPEN(UNIT=244,FILE="frnn_overlaps")
   OPEN(UNIT=245,FILE="frnn_overlap_species")

   overlap_frnn=0
   DO i=1,overlap_tnn
      coord=tnn_overlap_coords(3*i-2:3*i)
      species=tnn_overlap_species(i)
      DO l=1,12
         frnn=fourthnearestneighbours(3*l-2:3*l)+coord
         pos=GetIndex(ALref,frnn,cellsize,ncells) !identify site with coordinate fnn
         IF (ALref%AtomSpecies(pos)==0) THEN ! find coordination number of this site
            coordination=0
            DO k=1,12
               frnnv=firstnearestneighbours(3*k-2:3*k)+frnn
               pos=GetIndex(ALref,frnnv,cellsize,ncells)
               IF (ALref%AtomSpecies(pos)/=0) THEN
                  coordination=coordination+1
               END IF
            END DO
            frcn(l)=coordination
            !WRITE(6,*) l
            !WRITE(6,*) coordination
         ELSE
            frcn(l)=0
         END IF
      END DO
      l=MAXLOC(frcn,1)
      !WRITE(6,*) l
      frnn=fourthnearestneighbours(3*l-2:3*l)+coord
      pos=GetIndex(ALref,frnn,cellsize,ncells) !identify site with coordinate snn
      IF (ALref%AtomSpecies(pos)/=0) THEN
          overlap_frnn=overlap_frnn+1
          WRITE(244,*) coord !earlier this was tnn
          WRITE(245,*) species
      ELSE
          ALref%AtomSpecies(pos)=tnn_overlap_species(i)
      END IF
   END DO
   WRITE(6,*) "Number of collisions among tnn while matching:",overlap_frnn
   CLOSE(244)
   CLOSE(245)

   END SUBROUTINE
   !xxxxxxxxxxxxxxxxxxxxxxxxxx

END MODULE AlignNanoparticle
