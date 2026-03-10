!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE MPI_check

   IMPLICIT NONE

   CALL SetUpSpecies()
   CALL Welcome()

   CALL MPI_INIT(ierr)
   IF (ierr/=0) THEN
      WRITE(6,*) "Err>> Error in initializing MPI"
      STOP
   END IF

   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)
   IF (ierr/=0) THEN
      WRITE(6,*) "Err>> Error in initializing MPI_COMM_WORLD"
      STOP
   END IF
   num_procs=num_procs-1 !# slaves -- this is required by mpi_manager.f90

   CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
   IF (ierr/=0) THEN
      WRITE(6,*) "Err>> Error in obtaining rank"
      STOP
   END IF
   irank=irank+1
   charint2char=INT2CHAR(irank)
   TxtHeader="n"//TRIM(charint2char)//":"

END SUBROUTINE MPI_check
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Loadfile(AL,filename,NSites)

   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: AL
   CHARACTER(len=100) :: filename
   INTEGER,PARAMETER :: int64=selected_int_kind(15)
   INTEGER :: isite,NSites,i!,ispecies
   INTEGER(kind=int64) :: ispecies
   REAL(dp) :: coord(3),cb(3)

   OPEN(UNIT=102,FILE=filename)
   CALL MakeSize(AL,NSites)

   Do isite=1,NSites
     READ(102,*) i,ispecies,coord
     AL%AtomCoord(3*isite-2:3*isite)=coord
     AL%AtomSpecies(isite)=ispecies
   END DO
   CLOSE(102)

   AL%BoxSize(1)=boxsize_x*4.0801_dp
   AL%BoxSize(2)=boxsize_y*4.0801_dp
   AL%BoxSize(3)=boxsize_z*4.0801_dp
   cb=AL%BoxSize(1:3)/2.0  !center of box
   WRITE(6,*) "cb1:",cb

   AL%AtomSpecies=1
   CALL PotentialInitialize(AL)
   CALL AddVerletList(AL,ListType=1,NAtoms=NSites)

   OPEN(UNIT=102,FILE=filename)
   Do isite=1,NSites
     READ(102,*) i,ispecies,coord
     AL%AtomCoord(3*isite-2:3*isite)=coord
     AL%AtomSpecies(isite)=ispecies
   END DO
   CLOSE(102)

   WRITE(6,*) "Read file ", TRIM(filename)
END SUBROUTINE Loadfile
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Read_input_parameters 
   
   IMPLICIT NONE

   !KMC simulation
   !!!!!Prepare fcc(111) Structure of Pt and Replace Au Atoms randomly with
   !Fe Atoms
   OPEN(UNIT=442,FILE="kmc.parameters")
   READ(442,*) filename
   READ(442,*) CalculationType
   READ(442,*) Temperature1
   READ(442,*) composition
   READ(442,*) e0                             !new addition
   READ(442,*) prefac
   READ(442,*) xi
   READ(442,*) pnu_E
   READ(442,*) phi_E
   READ(442,*) iseed
   READ(442,*) tagx
   CALL InitializeTaus88(iseed)
   READ(442,*) PrintFreq
   READ(442,*) NWOrientation
   READ(442,*) IsReplaceAtoms
   READ(442,*) non_lattice 
   READ(442,*) NAtoms_input
   READ(442,*) NElectroactive_initial
   CLOSE(442)

   WRITE(6,*) "Using coordinates from file:",TRIM(filename)
   IF (CalculationType==1) THEN
      WRITE(6,*) "Temperature:",Temperature1
      WRITE(6,*) "Composition:",composition
      WRITE(6,*) "Diffusion prefactor:",prefac
      WRITE(6,*) "xi:",xi
      WRITE(6,*) "pnu_E:",pnu_E
      WRITE(6,*) "phi_E:",phi_E
      IF (NWOrientation==1) THEN
         WRITE(6,*) "(100)-oriented NW"
      ELSEIF (NWOrientation==2) THEN
         WRITE(6,*) "(110)-oriented NW"
      ENDIF
   ELSE
      WRITE(6,*) "ANALYSIS OF FRAME ..."
   END IF

END SUBROUTINE Read_input_parameters  
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Read_input_coordinates

   IMPLICIT NONE
   REAL(dp) :: coord(3),dr(3),cb(3),cm(3)
   INTEGER :: i
! Preparing AL1 the coordinates of particle

   IF (non_lattice) THEN   !load file
   filename="coordinates"
   CALL Loadfile(AL1,filename,NAtoms_input)   !nanoparticle
   !crystaltype="fcc"
   !filename="np.xyz"
   !CALL WriteXYZ(AL1,filename)
   ELSE
     filename="latchednp.xyz"
     CALL ReadCoord(AL1,filename)
   END IF

   cm=0._dp
   DO i=1,AL1%NAtoms
      cm=cm+AL1%AtomCoord(3*i-2:3*i)
   END DO
   cm=cm/REAL(AL1%NAtoms)

   cb(1)=(boxsize_x*4.0801)/2.0_dp
   cb(2)=(boxsize_y*4.0801)/2.0_dp
   cb(3)=(boxsize_z*4.0801)/2.0_dp
    WRITE(6,*) "cb2:",cb
    
   dr=cb-cm

   DO i=1,AL1%NAtoms
      coord=AL1%AtomCoord(3*i-2:3*i)
      coord=coord+dr
      AL1%AtomCoord(3*i-2:3*i)=coord
   END DO

! Preparing AL the coordinates of lattice points

   crystaltype="fcc"
   CALL Create(AL,(/int(boxsize_x),int(boxsize_y),int(boxsize_z)/),crystaltype,(/4.0801_dp/),(/79/))
   NSites=AL%NAtoms
   filename="box_AL.xyz"
   CALL WriteXYZ(AL,filename)

END SUBROUTINE Read_input_coordinates
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Initialize_variables(AL)

   IMPLICIT NONE

   TYPE(SystemContainer), POINTER :: AL

   !initialize few variables
   AtomSpecies=>AL%AtomSpecies
   AtomCoord=>AL%AtomCoord

   CALL PotentialInitialize(AL)
   CALL AddVerletList(AL,ListType=1,NAtoms=AL%NAtoms) !add verlet list to reference

   VL=>AL%VL
   VLList=>VL%List
   VLListRange=>VL%ListRange
   MaxAtomPerAtom=VL%MaxAtomPerAtom
   VLdrmag=>VL%drmag
   NSites=AL%NAtoms

   ALLOCATE(SpeciesType(NSites)) !species type: electroactive/noble/vacancy
   ALLOCATE(Coordination(NSites)) !coordination number
   ALLOCATE(NeighborSiteList(12*NSites)) !site index that is neighnor
   ALLOCATE(NumberNeighborSites(NSites)) !number of site not occupied by vacancy
   ALLOCATE(FirstNearestNeighborSites(NSites))
   ALLOCATE(ThirdNearestNeighborSites(NSites))
   ALLOCATE(Local_composition(NSites))
   ALLOCATE(Rate(13*NSites)) !rate @ each lattice site
   ALLOCATE(PrintableAtom(NSites)) !
   ALLOCATE(DiffusionFlux(NSites))
   ALLOCATE(SurfaceInfo(NSites))
   SurfaceInfo=" "
   ALLOCATE(IsSurf111atom(NSites),IsSurf100atom(NSites),IsSurfOtherAtom(NSites),mask(NSites))
   ALLOCATE(IsSubSurf111atom(NSites),IsSubSurf100atom(NSites))

END SUBROUTINE Initialize_variables
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Setup_Replaceatoms(AL)

   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: AL
   !replace atoms
   !IsReplaceAtoms=.TRUE. !to replace atoms value should be .TRUE.
   IF (CalculationType==1 .AND. IsReplaceAtoms) THEN
      DO isite=1,AL%NAtoms
         IF (AL%AtomSpecies(isite)==0) CYCLE

         IF (taus88()>composition) THEN
            AL%AtomSpecies(isite)=1  !electroactive species inserted
         END IF
      END DO
   END IF

END SUBROUTINE Setup_Replaceatoms
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx  
SUBROUTINE Print_info_of_LatchedParticle(AL)

   IMPLICIT NONE
   TYPE(SystemContainer), POINTER :: AL
   INTEGER :: i

   !electroactive and vacant sites
   NVacancies=0 !number of vacant sites
   NElectroactive=0
   Nonactive=0
   DO i=1,AL%NAtoms
      SELECT CASE(AL%AtomSpecies(i))
      CASE(0); NVacancies=NVacancies+1
      CASE(1); NElectroactive=NElectroactive+1
      CASE(2); Nonactive=Nonactive+1
      END SELECT
   END DO
!   NSites=AL%NAtoms
   NParticle=NElectroactive+Nonactive

   WRITE(6,*) " After latching the nanoparticle and taking collisions into account:"
   WRITE(6,*) " # electroactive species:",NElectroactive
   WRITE(6,*) " # vacant sites:",NVacancies
   WRITE(6,*) " # nonactive species:",Nonactive
   WRITE(6,*) " # number of atoms in the particle",NParticle

END SUBROUTINE Print_info_of_LatchedParticle
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Initialize_Species_and_Neighbours

   IMPLICIT NONE

   !initialize few variables
   SpeciesType=0
   NeighborSiteList=0
   Coordination=0
   NumberNeighborSites=0

   !setup species
   DO isite=1,NSites
      CALL SetupSpeciesType(isite)
   END DO
   !setup neighbor list
   DO isite=1,NSites
      CALL NearestNeighborSiteList(isite,isprint=.FALSE.,coreradius=coreradius,center=center)
   END DO

END SUBROUTINE Initialize_Species_and_Neighbours
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Print_LatchedParticle

   IMPLICIT NONE

   NDissolved=0
   CALL GetPrintableAtoms(.FALSE.)
   filename="latchedslab.xyz"
   OPEN(UNIT=302,FILE=filename)
   CALL PrintFrame1(302,NPrintableAtoms,PrintableAtom,0)
   CLOSE(302)

END SUBROUTINE Print_LatchedParticle
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Count_surface_atoms

   IMPLICIT NONE

   !Print coordination statistics
   cntr1=0; cntr2=0; cntr3=0
   DO isite=1,NSites
      IF (Coordination(isite)<12 .AND. SpeciesType(isite)==1) cntr1=cntr1+1
      IF (Coordination(isite)<12 .AND. SpeciesType(isite)==2) cntr2=cntr2+1
   END DO
   WRITE(6,*) "# surface Fe atoms:",cntr1
   WRITE(6,*) "# surface Pt atoms:",cntr2
   WRITE(6,*) "# surface atoms:",cntr1+cntr2

END SUBROUTINE Count_surface_atoms
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE Initialize_KMC

   IMPLICIT NONE
   !initialize rate array
   time=0._dp
   kmcx=0. !# kmc moves
   kmcd=0. !# kmc dissolution events
   Rate=0._dp

   CALL InitializeRates()
   DO isite=1,NSites
      CALL UpdateRateArray(isite)
   END DO
   write(*,*) ANY(Rate>0.)
   CALL BinaryTreeInitialize(Rate,13*NSites)
   write(*,*) ANY(Rate>0.)
   BT=>BinaryTree
   CALL PrintFrame()

END SUBROUTINE Initialize_KMC
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE translate_particle_to_corner

   IMPLICIT NONE

   REAL(dp) :: c(3),a,b,d
   REAL(dp), DIMENSION(193530) :: X_AtomCoord,Y_AtomCoord,Z_AtomCoord
   REAL(dp), DIMENSION(193530*3) :: UpdatedCoord

   AtomCoord=>AL1%AtomCoord
   DO isite=1,193530
      X_AtomCoord(isite)=AtomCoord(3*isite-2)
      Y_AtomCoord(isite)=AtomCoord(3*isite-1)
      Z_AtomCoord(isite)=AtomCoord(3*isite)   
   END DO

   !finding xmin,ymin,zmin
   c(1)=MINVAL(X_AtomCoord)
   c(2)=MINVAL(Y_AtomCoord)
   c(3)=MINVAL(Z_AtomCoord)
   WRITE(6,*) c

   !translating the particle to the corner
   DO isite=1,193530
      UpdatedCoord(3*isite-2:3*isite)=AtomCoord(3*isite-2:3*isite)-c
   END DO

   !moving atoms to the grid points
   UpdatedCoord=NINT(UpdatedCoord/2.04005)
   UpdatedCoord=UpdatedCoord*2.04005_dp   

   !moving non-lattice atoms to the lattice points
   DO isite=1,193530
      a=MOD(UpdatedCoord(3*isite-2),4.0801)
      b=MOD(UpdatedCoord(3*isite-1),4.0801)
      d=MOD(UpdatedCoord(3*isite),4.0801)

      IF (a .EQ. 0.0 .AND. b .EQ. 0.0 .AND. d .NE. 0.0) UpdatedCoord(3*isite)=UpdatedCoord(3*isite)+2.04005
      IF (a .NE. 0.0 .AND. b .NE. 0.0 .AND. d .NE. 0.0) UpdatedCoord(3*isite)=UpdatedCoord(3*isite)+2.04005
      IF (a .NE. 0.0 .AND. b .EQ. 0.0 .AND. d .EQ. 0.0) UpdatedCoord(3*isite)=UpdatedCoord(3*isite)+2.04005
      IF (a .EQ. 0.0 .AND. b .NE. 0.0 .AND. d .EQ. 0.0) UpdatedCoord(3*isite)=UpdatedCoord(3*isite)+2.04005
   END DO

   OPEN(UNIT=302,FILE="Axes")
      WRITE(302,*) UpdatedCoord
   CLOSE(302)

END SUBROUTINE translate_particle_to_corner
!xxxxxxxxxxxxxxxxxxxxxxxxxxxx
