MODULE IO_PDB
!contains subroutines for reading and printing pdb files
   USE VARIABLE_TYPE
   USE db_manipulate
   USE OptimPackage
   IMPLICIT NONE
   INTEGER, PARAMETER :: MaxAtom=10000
   
   !variables used by code
   
   CHARACTER(len=100), PRIVATE :: Header
   CHARACTER(len=6), DIMENSION(MaxAtom), PRIVATE :: recordtype
   INTEGER, DIMENSION(MaxAtom), PRIVATE :: atomID
   CHARACTER(len=4), DIMENSION(MaxAtom), PRIVATE :: atomname
   CHARACTER(len=3), DIMENSION(MaxAtom), PRIVATE :: residuename
   CHARACTER(len=3), DIMENSION(MaxAtom), PRIVATE :: segmentname,charge
   CHARACTER(len=1), DIMENSION(MaxAtom), PRIVATE :: chainID,altLoc
   INTEGER, DIMENSION(MaxAtom), PRIVATE :: residueID
   REAL, DIMENSION(3*MaxAtom), PRIVATE :: coord
   REAL, DIMENSION(MaxAtom), PRIVATE :: occupancy,temperaturefactor
   
   REAL, DIMENSION(3*MaxAtom), PRIVATE :: Reducedcoord
   REAL, DIMENSION(3*MaxAtom), PRIVATE, TARGET :: pdbsys1,pdbsys2
   REAL, DIMENSION(3*MaxAtom), PRIVATE :: pdbsys10,pdbsys20
   
   INTEGER, PRIVATE :: NAtoms,ReducedNAtoms
   
   TYPE PDBContainer
      CHARACTER(len=100) :: Header
      CHARACTER(len=6), DIMENSION(:), POINTER :: recordtype
      INTEGER, DIMENSION(:), POINTER :: atomID
      CHARACTER(len=4), DIMENSION(:), POINTER :: atomname
      CHARACTER(len=3), DIMENSION(:), POINTER :: residuename
      CHARACTER(len=3), DIMENSION(:), POINTER :: segmentname,charge
      CHARACTER(len=1), DIMENSION(:), POINTER :: chainID,altLoc
      INTEGER, DIMENSION(:), POINTER :: residueID
      REAL, DIMENSION(:), POINTER :: coord
      REAL, DIMENSION(:), POINTER :: occupancy,temperaturefactor
      
      REAL, DIMENSION(:), POINTER :: Reducedcoord
      
      INTEGER, PRIVATE :: NAtoms,ReducedNAtoms
   END TYPE PDBContainer
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReadPDB(filename)
   !iopt=0  read l
   !iopt=1  ignore water
      IMPLICIT NONE
      CHARACTER(len=6) :: recordtype1
      INTEGER :: atomID1
      CHARACTER(len=4) :: atomname1
      CHARACTER(len=3) :: residuename1
      CHARACTER(len=3) :: segmentname1,charge1
      CHARACTER(len=1) :: buffer1,tloc1,chainID1
      INTEGER :: residueID1
      REAL :: x,y,z,occupancy1,temperaturefactor1
      CHARACTER(len=100) :: filename,buffer100
      CHARACTER(len=6) :: buffer10
      LOGICAL :: NotSatisfied
      INTEGER :: iunit,errorstatus,iatom
      
      NotSatisfied=.TRUE.
      iunit=424
      iatom=0
      OPEN(UNIT=iunit,FILE=TRIM(filename))
      DO WHILE (NotSatisfied)
         READ(UNIT=iunit,FMT='(a6)',ADVANCE="NO",IOSTAT=errorstatus) recordtype1
         SELECT CASE(TRIM(recordtype1))
         CASE("REMARK")
            READ(UNIT=iunit,FMT='(a100)',IOSTAT=errorstatus) buffer100
            Header=recordtype1//buffer100
         CASE("CRYST1")
            READ(UNIT=iunit,FMT='(a100)',IOSTAT=errorstatus) buffer100
            Header=recordtype1//buffer100
         CASE("ATOM  ")
            READ(UNIT=iunit,FMT='(i5,a1,a4,a1)',ADVANCE="NO") atomID1,buffer1,atomname1,tloc1
            READ(UNIT=iunit,FMT='(a3,2a1)',ADVANCE="NO") residuename1,buffer1,chainID1
            READ(UNIT=iunit,FMT='(i4,4a1)',ADVANCE="NO") residueID1,buffer1,buffer1,buffer1,buffer1
            READ(UNIT=iunit,FMT='(3f8.3,2f6.2)',ADVANCE="NO") x,y,z,occupancy1,temperaturefactor1
            READ(UNIT=iunit,FMT='(a6)',ADVANCE="NO") buffer10
            READ(UNIT=iunit,FMT='(2a3)') segmentname1,charge1
         CASE("END")
            EXIT
         END SELECT
         IF (recordtype1=="ATOM" .OR. recordtype1=="HETATM") THEN
            iatom=iatom+1
            recordtype(iatom)=recordtype1
            atomID(iatom)=atomID1
            atomname(iatom)=atomname1
            altloc(iatom)=tloc1
            residuename(iatom)=residuename1
            chainID(iatom)=chainID1
            residueID(iatom)=residueID1
            coord(3*iatom-2:3*iatom)=(/x,y,z/)
            occupancy(iatom)=occupancy1
            temperaturefactor(iatom)=temperaturefactor1
            segmentname(iatom)=segmentname1
            charge(iatom)=charge1
         END IF
      END DO
      CLOSE(iunit)
      NAtoms=iatom
      WRITE(6,*) "Read PDB file from ...",TRIM(filename)
   END SUBROUTINE ReadPDB
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WritePDB(filename)
      IMPLICIT NONE
      CHARACTER(len=100) :: filename
      INTEGER :: iatom,iunit
      CHARACTER :: buffer1
      
      iunit=424
      iatom=0
      buffer1=" "
      OPEN(UNIT=iunit,FILE=TRIM(filename))
      WRITE(UNIT=iunit,FMT='(a100)') Header(1:100)
      DO iatom=1,NAtoms
         WRITE(UNIT=iunit,FMT='(a6,i5,a1,a4,a1)',ADVANCE="NO") recordtype(iatom),atomID(iatom),buffer1, &
            atomname(iatom),altLoc(iatom)
         WRITE(UNIT=iunit,FMT='(a3,2a1)',ADVANCE="NO") residuename(iatom),buffer1,chainID(iatom)
         WRITE(UNIT=iunit,FMT='(i4,4a1)',ADVANCE="NO") residueID(iatom),buffer1,buffer1,buffer1,buffer1
         WRITE(UNIT=iunit,FMT='(3f8.3,2f6.2)',ADVANCE="NO") coord(3*iatom-2:3*iatom),occupancy(iatom), &
           temperaturefactor(iatom)
         WRITE(UNIT=iunit,FMT='("      ",2a3)') segmentname(iatom),charge(iatom)
      END DO
      WRITE(iunit,FMT='("END")')
      CLOSE(iunit)
      WRITE(6,*) "Wrote PDB file to ...",TRIM(filename)
   END SUBROUTINE WritePDB
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SelectWaterInPDB(mask)
      IMPLICIT NONE
      LOGICAL, DIMENSION(:), POINTER :: mask
      INTEGER :: i
      
      IF (ASSOCIATED(mask)) THEN
         WRITE(6,*) "mask in SelectWaterInPDB cannot be associated"
         STOP
      END IF
      
      ALLOCATE(mask(NAtoms))
      mask=.FALSE.
      DO i=1,NAtoms
         IF (residuename(i)=="TIP") mask(i)=.TRUE.
      END DO
   END SUBROUTINE SelectWaterInPDB
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ConvertPDB2PDBAL(PDBAL,mask)
   !pdb data stored in global arrays will be transferred to PDBAL
      IMPLICIT NONE
      TYPE(PDBContainer), POINTER :: PDBAL
      LOGICAL, DIMENSION(:), OPTIONAL :: mask
      LOGICAL, DIMENSION(:), ALLOCATABLE :: mask1
      INTEGER :: i,nenteries
      
      IF (PRESENT(mask)) THEN
         ALLOCATE(mask1(SIZE(mask)))
         mask1=mask
         nenteries=0
         DO i=1,SIZE(mask)
            IF (mask(i)) nenteries=nenteries+1
         END DO
      ELSE
         ALLOCATE(mask1(NAtoms))
         mask1=.TRUE.
         nenteries=NAtoms
      END IF
      
      IF (ASSOCIATED(PDBAL)) CALL DeletePDBAL(PDBAL)
      CALL MakePDBAL(PDBAL,nenteries)
      PDBAL%NAtoms=nenteries
      PDBAL%ReducedNAtoms=ReducedNAtoms
      DO i=1,nenteries
         IF (mask1(i)) THEN
            PDBAL%recordtype(i)=recordtype(i)
            PDBAL%atomID(i)=atomID(i)
            PDBAL%atomname(i)=atomname(i)
            PDBAL%residuename(i)=residuename(i)
            PDBAL%segmentname(i)=segmentname(i)
            PDBAL%charge(i)=charge(i)
            PDBAL%chainID(i)=chainID(i)
            PDBAL%altLoc(i)=altLoc(i)
            PDBAL%residueID(i)=residueID(i)
            PDBAL%coord(3*i-2:3*i)=coord(3*i-2:3*i)
            PDBAL%occupancy(i)=occupancy(i)
            PDBAL%temperaturefactor(i)=temperaturefactor(i)
            PDBAL%Reducedcoord(3*i-2:3*i)=Reducedcoord(3*i-2:3*i)
         END IF
      END DO
      
      DEALLOCATE(mask1)
   END SUBROUTINE ConvertPDB2PDBAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ConvertPDBAL2PDB(PDBAL)
      IMPLICIT NONE
      TYPE(PDBContainer), POINTER :: PDBAL
      INTEGER :: i
      
      NAtoms=PDBAL%NAtoms
      ReducedNAtoms=PDBAL%ReducedNAtoms
      DO i=1,NAtoms
         recordtype(i)=PDBAL%recordtype(i)
         atomID(i)=PDBAL%atomID(i)
         atomname(i)=PDBAL%atomname(i)
         residuename(i)=PDBAL%residuename(i)
         segmentname(i)=PDBAL%segmentname(i)
         charge(i)=PDBAL%charge(i)
         chainID(i)=PDBAL%chainID(i)
         altLoc(i)=PDBAL%altLoc(i)
         residueID(i)=PDBAL%residueID(i)
         coord(3*i-2:3*i)=PDBAL%coord(3*i-2:3*i)
         occupancy(i)=PDBAL%occupancy(i)
         temperaturefactor(i)=PDBAL%temperaturefactor(i)
         Reducedcoord(3*i-2:3*i)=PDBAL%Reducedcoord(3*i-2:3*i)
      END DO
   END SUBROUTINE ConvertPDBAL2PDB
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeletePDBAL(PDBAL)
      IMPLICIT NONE
      TYPE(PDBContainer), POINTER :: PDBAL
      
      DEALLOCATE(PDBAL)
      
   END SUBROUTINE DeletePDBAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MakePDBAL(PDBAL,NAtoms)
      IMPLICIT NONE
      TYPE(PDBContainer), POINTER :: PDBAL
      INTEGER :: NAtoms
      
      ALLOCATE(PDBAL)
      PDBAL%NAtoms=NAtoms
      ALLOCATE(PDBAL%recordtype(NAtoms))
      ALLOCATE(PDBAL%atomID(NAtoms))
      ALLOCATE(PDBAL%atomname(NAtoms))
      ALLOCATE(PDBAL%residuename(NAtoms))
      ALLOCATE(PDBAL%segmentname(NAtoms))
      ALLOCATE(PDBAL%charge(NAtoms))
      ALLOCATE(PDBAL%chainID(NAtoms))
      ALLOCATE(PDBAL%altLoc(NAtoms))
      ALLOCATE(PDBAL%residueID(NAtoms))
      ALLOCATE(PDBAL%coord(3*NAtoms))
      ALLOCATE(PDBAL%occupancy(NAtoms))
      ALLOCATE(PDBAL%temperaturefactor(NAtoms))
      ALLOCATE(PDBAL%Reducedcoord(3*NAtoms))
   END SUBROUTINE MakePDBAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateNAMD_MD(filename,psffile,pdbfile)
      IMPLICIT NONE
      CHARACTER(len=100) :: filename,psffile,pdbfile
      INTEGER :: iunit
      
      iunit=424
      OPEN(UNIT=iunit,FILE=filename)
      WRITE(UNIT=iunit,FMT='("set outName	        output/NP_300")')
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# MD SETUP")')
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("timestep                2.0")') 
      WRITE(UNIT=iunit,FMT='("numSteps                5000000")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# INPUT")') 
      WRITE(UNIT=iunit,FMT='("structure   	           ")',ADVANCE="NO") 
      WRITE(iunit,*) TRIM(psffile)
      WRITE(UNIT=iunit,FMT='("coordinates	          ",a100)',ADVANCE="NO") 
      WRITE(iunit,*) TRIM(pdbfile)
      WRITE(UNIT=iunit,FMT='("binCoordinates  	   output/equilsolv.coor.old")') 
      WRITE(UNIT=iunit,FMT='("binVelocities              output/equilsolv.vel.old")') 
      WRITE(UNIT=iunit,FMT='("extendedSystem             output/equilsolv.xsc")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# OUTPUT")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("outputenergies          2500")') 
      WRITE(UNIT=iunit,FMT='("outputtiming            2500")') 
      WRITE(UNIT=iunit,FMT='("restartfreq             2500")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("binaryOutput            no")') 
      WRITE(UNIT=iunit,FMT='("binaryRestart           yes")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("outputname              ",a22)')  '[format "%s" $outName]'
      WRITE(UNIT=iunit,FMT='("restartname             $outName")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# DCD")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("DCDfile                 ${outName}.dcd")') 
      WRITE(UNIT=iunit,FMT='("DCDfreq                 1000")') 
      WRITE(UNIT=iunit,FMT='("DCDUnitCell             yes")') 
      WRITE(UNIT=iunit,FMT='("xstFreq            	2500")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("parameters              ../../common/par_all36_prot.prm")') 
      WRITE(UNIT=iunit,FMT='("paraTypeCharmm          on")')
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# PME")')
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("PME                     yes")')
      WRITE(UNIT=iunit,FMT='("PMETolerance            10e-6")')
      WRITE(UNIT=iunit,FMT='("PMEInterpOrder          4")')
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("PMEGridSpacing          1.0")')
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# WRAP ATOMS FOR OUTPUT")')
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("wrapAll                 on")')
      WRITE(UNIT=iunit,FMT='("wrapNearest             yes")')
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# CONSTANT-T")')
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("langevin                on")')
      WRITE(UNIT=iunit,FMT='("langevinTemp            300.0")')
      WRITE(UNIT=iunit,FMT='("langevinDamping         1.0")')
      WRITE(UNIT=iunit,FMT='("langevinHydrogen         off")')
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# SPACE PARTITIONING")')
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("splitpatch              hydrogen")') 
      WRITE(UNIT=iunit,FMT='("hgroupcutoff            2.8")') 
      WRITE(UNIT=iunit,FMT='("stepspercycle           20")') 
      WRITE(UNIT=iunit,FMT='("margin                  1.0")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# CUT-OFFS")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("switching               on")') 
      WRITE(UNIT=iunit,FMT='("switchdist               9.0")') 
      WRITE(UNIT=iunit,FMT='("cutoff                  11.0")') 
      WRITE(UNIT=iunit,FMT='("pairlistdist            13.0")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# RESPA ")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("fullElectFrequency      2")') 
      WRITE(UNIT=iunit,FMT='("nonbondedFreq           1")') 
      WRITE(UNIT=iunit,FMT='("# 1-4 NON-BONDED")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("exclude                 scaled1-4")') 
      WRITE(UNIT=iunit,FMT='("1-4scaling              1.0")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# COM")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("commotion               no")') 
      WRITE(UNIT=iunit,FMT='("zeroMomentum            yes")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("# SHAKE")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT='("rigidbonds              l")') 
      WRITE(UNIT=iunit,FMT=*) 
      WRITE(UNIT=iunit,FMT=*) 
      CLOSE(iunit)
   END SUBROUTINE CreateNAMD_MD
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateNAMD_Optimize(filename)
   !optimize a pdb file
      IMPLICIT NONE
      CHARACTER(len=100) :: filename,psffile,pdbfile
   END SUBROUTINE CreateNAMD_Optimize
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EliminateRedundantDegrees()
      IMPLICIT NONE
      INTEGER :: iatom,jatom
      
      jatom=0
      DO iatom=1,NAtoms
         IF (atomname(iatom)==" OH2") THEN
            !WRITE(6,*) "Encountered water"
         ELSE
            jatom=jatom+1
            Reducedcoord(3*jatom-2:3*jatom)=coord(3*iatom-2:3*iatom)
            WRITE(*,*) atomname(iatom),Reducedcoord(3*jatom-2:3*jatom)
         END IF
      END DO
      ReducedNAtoms=jatom
      WRITE(6,*) "Number of molecules in reduced system:",ReducedNAtoms
      WRITE(6,*) "Number of molecules in full system:",NAtoms
   END SUBROUTINE EliminateRedundantDegrees
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION IsPDBALMatch(PDBAL1,PDBAL2)
   !determines whether PDBAL1 and PDBAL2 match
      IMPLICIT NONE
      TYPE(PDBContainer), POINTER :: PDBAL1,PDBAL2
      REAL, DIMENSION(:), POINTER :: x
      LOGICAL :: IsPDBALMatch
      INTEGER :: NAtoms,ITMAX,iter,errorstatus1,i
      REAL :: ftol,xtol,gtol,fret,com1(3),com2(3),tra(3)
      
      NAtoms=PDBAL1%NAtoms
      IsPDBALMatch=.FALSE.
      IF (PDBAL2%NAtoms/=NAtoms) THEN
         IsPDBALMatch=.FALSE.
         RETURN
      END IF
      
      pdbsys10(1:3*NAtoms)=PDBAL1%coord(1:3*NAtoms) !reference
      pdbsys20(1:3*NAtoms)=PDBAL2%coord(1:3*NAtoms) !reference
      pdbsys2(1:3*NAtoms)=PDBAL2%coord(1:3*NAtoms) !to be optimized
      
      !com
      com1=0.
      com2=0.
      DO i=1,NAtoms
         com1=com1+pdbsys10(3*i-2:3*i)
         com2=com2+pdbsys2(3*i-2:3*i)
      END DO
      com1=com1/REAL(NAtoms)
      com2=com2/REAL(NAtoms)
      tra=com1-com2
      !DO i=1,NAtoms
      !   pdbsys2(3*i-2:3*i)=pdbsys2(3*i-2:3*i)+tra
      !END DO
      
      x=>pdbsys2(1:3*NAtoms)
      ftol=1.e-4
      gtol=1.e-4
      xtol=1.e-4
      ITMAX=300
      CALL SteepestDescent(x,ftol,gtol,xtol,iter,fret,pdbdistance,dpdbdistance, &
         dx=0.0001,ITMAX=ITMAX,IsPrint1=.FALSE.,errorstatus=errorstatus1)
      
      CALL Lbfgs(x,ftol,gtol,xtol,iter,fret,pdbdistance,dpdbdistance,IsPrint1=.FALSE., &
         ITMAX=ITMAX,errorstatus=errorstatus1)
      STOP
   END FUNCTION IsPDBALMatch
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION pdbdistance(x,errorstatus)
      USE VARIABLE_TYPE
      IMPLICIT NONE
      INTEGER :: errorstatus,NAtoms,i
      REAL, DIMENSION(:), INTENT(IN) :: x
      REAL :: pdbdistance,coord(3),length2,length20
      
      pdbdistance=0.
      errorstatus=0
      !use stiff springs within pdb2 and weak springs between pdb1 and pdb2
      NAtoms=SIZE(x)/3
      DO i=1,NAtoms-1
         coord=(x(3*i-2:3*i)-x(3*i+1:3*i+3))
         length2=DOT_PRODUCT(coord,coord)
         coord=pdbsys20(3*i-2:3*i)-pdbsys20(3*i+1:3*i+3)
         length20=DOT_PRODUCT(coord,coord)
         
         pdbdistance=pdbdistance+10.*ABS(length2-length20) !spring constant 100
      END DO
      DO i=1,NAtoms
         coord=pdbsys2(3*i-2:3*i)-pdbsys10(3*i+1:3*i+3)
         length2=DOT_PRODUCT(coord,coord) !spring constant 1.
         pdbdistance=pdbdistance+length2
      END DO
      write(*,*) pdbdistance
   END FUNCTION pdbdistance
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION dpdbdistance(x,errorstatus)
      USE VARIABLE_TYPE
      IMPLICIT NONE
      INTEGER :: errorstatus,i
      REAL, DIMENSION(:), INTENT(IN) :: x
      REAL, DIMENSION(SIZE(x)) :: dpdbdistance,x1
      REAL :: fp,fn
      
      dpdbdistance=0.
      errorstatus=0
      NAtoms=SIZE(x)/3
      x1=x
      DO i=1,3*NAtoms
         x1(i)=x1(i)+0.01
         fp=pdbdistance(x1,errorstatus)
         x1(i)=x1(i)-0.02
         fn=pdbdistance(x1,errorstatus)
         x1(i)=x1(i)+0.01
         dpdbdistance(i)=(fp-fn)/0.02
      END DO
   END FUNCTION dpdbdistance
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION rotatepdb(x,errorstatus)
   !We require pdbsys2 to be rotated about com1 so that it matches with pdbsys1
   !suppose r1 is vector in pdbsys1 and r2 is corresponding vector in pdbsys2
   !R is the rotation matrix
   !Ideally, R*r2=r1
   !here R=[v1 v2 v3; v4 v5 v6; v7 v8 v9] is 3x3 matrix
   !x is [v1 v2 v3 v4 v5 v6 v7 v8 v9], we want to optimize x to get best match
   !rotatepdb finds the sum of squared distances between the vectors in pdbsys1 and pdbsys2
   ! after rotating by the matrix x
      REAL, DIMENSION(:), INTENT(IN) :: x
      REAL :: coord1(3),coord2(3),rotatepdb
      INTEGER :: errorstatus,i,n
      
      rotatepdb=0.
      DO i=1,n !for l vectors/molecules in pdb
         coord1=pdbsys1(3*i-2:3*i)
         coord2=pdbsys2(3*i-2:3*i)
         rotatepdb=rotatepdb+(DOT_PRODUCT(x(1:3),coord2)-coord1(1))**2
         rotatepdb=rotatepdb+(DOT_PRODUCT(x(4:6),coord2)-coord1(2))**2
         rotatepdb=rotatepdb+(DOT_PRODUCT(x(7:9),coord2)-coord1(3))**2
      END DO
   END FUNCTION rotatepdb
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION drotatepdb(x,errorstatus)
   !find gradient of rotatepdb
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN) :: x
      REAL, DIMENSION(SIZE(x)) :: drotatepdb,xuse
      REAL :: fp,fn
      INTEGER :: errorstatus,i
      
      drotatepdb=0.
      xuse=x
      DO i=1,9
         xuse(i)=xuse(i)+0.001
         fp=rotatepdb(xuse,errorstatus)
         xuse(i)=xuse(i)-0.002
         fn=rotatepdb(xuse,errorstatus)
         xuse(i)=xuse(i)+0.001
         drotatepdb(i)=(fp-fn)/0.002
      END DO
   END FUNCTION drotatepdb
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE IO_PDB
