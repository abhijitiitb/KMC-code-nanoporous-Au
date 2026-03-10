MODULE MDProcess
!Principal Author: Abhijit Chatterjee
!Department of Chemical Engineering, IIT Kanpur
!Provided by user find the processes possible in the current system
   USE KMC_VARIABLE_TYPE
   USE KMCDbManipulate
   USE AMDProcess
   USE utilities
   USE KMCConfigurations
   USE KMCRelax
   
   IMPLICIT NONE
   INTEGER :: NDomainX,NDomainY,NDomainZ
   REAL(dp) :: LDomainx,LDomainy,LDomainz
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateProcessListFullMD()
   !used for create process list of full system
      IMPLICIT NONE
      INTEGER :: ix,iy,iz,ifile
      INTEGER :: NFilesMD
      
      !do in three stages
      !Stage 1 : set up xyz file domainwise
      ifile=0
      NFilesMD=0
      DO iz=1,NDomainZ
         DO iy=1,NDomainY
            DO ix=1,NDomainX
               ifile=ix+((iy-1)*NdomainX)+((iz-1)*NdomainY*NdomainX)
               NFilesMD=NFilesMD+1
               CALL LEKMCSetUpXYZ(ix,iy,iz,ifile)
            END DO
         END DO
      END DO
      
      !Stage 2 : Perform MD 
      CALL LEKMCProcessSearch(NFilesMD)
      
      !Stage 3 : Generate Local Processes
      CALL GenerateLEKMCProcessListMD()
   END SUBROUTINE GenerateProcessListFullMD
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LEKMCSetUpXYZ(ix,iy,iz,ifile)
   !sets up the files for doing MD
      IMPLICIT NONE
      INTEGER :: ix,iy,iz,ifile
      INTEGER :: icx,icy,icz,iunit,icmin(3),icmax(3)
      INTEGER :: NAtoms,NMove(3),BoxSize(6),NAtomsInRegion(4),NMoveInRegion(12)
      REAL(dp) :: rangecoord(6) !contains the min and max atom coordinates along x, y and z
      REAL(dp) :: LBuffer,rc,vacuum
      CHARACTER(len=100) :: filename
      CHARACTER(len=NINT2CHAR) :: charifile
      
      NAtoms=0
      NMove=0
      NAtomsInRegion=0
      NMoveInRegion=0
      
      
      charifile=INT2CHAR(ifile)
      filename="LEKMCMD."//TRIM(charifile)//".xyz"
      OPEN(UNIT=351,FILE=TRIM(filename))
      
      !find the regions and the cells in the regions
      !region 1
      rangecoord=0._dp
      rangecoord(1)=REAL(ix-1,dp)*LDomainx; rangecoord(2)=REAL(ix,dp)*LDomainx
      rangecoord(3)=REAL(iy-1,dp)*LDomainy; rangecoord(4)=REAL(iy,dp)*LDomainy
      rangecoord(5)=REAL(iz-1,dp)*LDomainz; rangecoord(6)=REAL(iz,dp)*LDomainz
      CALL AddAffectedCellsFromLEKMCMD(rangecoord,value=1)
      
      !region 2
      rangecoord(1)=rangecoord(1)-LBuffer
      rangecoord(2)=rangecoord(2)+LBuffer
      rangecoord(3)=rangecoord(3)-LBuffer
      rangecoord(4)=rangecoord(4)+LBuffer
      rangecoord(5)=rangecoord(5)-LBuffer
      rangecoord(6)=rangecoord(6)+LBuffer
      CALL AddAffectedCellsFromLEKMCMD(rangecoord,value=2)
      
      !region 3
      rangecoord(1)=rangecoord(1)-rc
      rangecoord(2)=rangecoord(2)+rc
      rangecoord(3)=rangecoord(3)-rc
      rangecoord(4)=rangecoord(4)+rc
      rangecoord(5)=rangecoord(5)-rc
      rangecoord(6)=rangecoord(6)+rc
      CALL AddAffectedCellsFromLEKMCMD(rangecoord,value=3)
      
      !region 4
      rangecoord(1)=rangecoord(1)-rc
      rangecoord(2)=rangecoord(2)+rc
      rangecoord(3)=rangecoord(3)-rc
      rangecoord(4)=rangecoord(4)+rc
      rangecoord(5)=rangecoord(5)-rc
      rangecoord(6)=rangecoord(6)+rc
      CALL AddAffectedCellsFromLEKMCMD(rangecoord,value=4)
      
      !write number of atoms to file
      WRITE(UNIT=351,FMT='(20i8)') SUM(NAtomsInRegion),SUM(NMoveInRegion),NAtomsInRegion,NMoveInRegion
      BoxSize=0._dp
      BoxSize(4:6)=(/rangecoord(2)-rangecoord(1),rangecoord(4)-rangecoord(3), &
         rangecoord(6)-rangecoord(5)/)
      vacuum=20._dp
      BoxSize(4:6)=BoxSize(4:6)+vacuum
      WRITE(UNIT=351,FMT='(6F15.8,3I5," LEKMC-MD Domain Index")') BoxSize,ix,iy,iz
      
      CALL AddCellContent(value=1,iunit=351) !add atoms in region 1 to file filename
      CALL AddCellContent(value=2,iunit=351) !add atoms in region 2 to file filename
      CALL AddCellContent(value=3,iunit=351) !add atoms in region 3 to file filename
      CALL AddCellContent(value=4,iunit=351) !add atoms in region 4 to file filename
      
      !reset the affected cells to zero
      CALL ResetAffectedCellInfo()
      
      CLOSE(351)
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE AddAffectedCellsFromLEKMCMD(rangecoord,value)
      !add the periodic boundary condition part
         IMPLICIT NONE
         INTEGER :: value
         REAL(dp) :: rangecoord(6)
         INTEGER :: icx,icy,icz,icmin(3),icmax(3)
         INTEGER :: icx1,icy1,icz1
         
         icmin=(/CEILING(rangecoord(1)/LDomainx),CEILING(rangecoord(3)/LDomainy), &
            CEILING(rangecoord(5)/LDomainz)/)
         icmax=(/CEILING(rangecoord(2)/LDomainx),CEILING(rangecoord(4)/LDomainy), &
            CEILING(rangecoord(6)/LDomainz)/)
         DO icz=icmin(3),icmax(3)
            DO icy=icmin(2),icmax(2)
               DO icx=icmin(1),icmax(1)
                  icx1=iPBC(icx,NDomainX)
                  icy1=iPBC(icy,NDomainY)
                  icz1=iPBC(icz,NDomainZ)
                  NAtomsInRegion(value)=NAtomsInRegion(value)+Cell(icx1,icy1,icz1)%NAtomsInCell
                  NMoveInRegion(3*value-2:3*value)=NMoveInRegion(3*value-2:3*value)+ &
                     Cell(icx1,icy1,icz1)%NMoveInCell
                  CALL AddCellToAffectedList((/icx1,icy1,icz1/),value=value)
               END DO
            END DO
         END DO
      END SUBROUTINE AddAffectedCellsFromLEKMCMD
   END SUBROUTINE LEKMCSetUpXYZ
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddCellContent(value,iunit)
   !adds the atom coordinates to a file
      IMPLICIT NONE
      INTEGER :: value,icx,icy,icz,iunit
      
      DO icz=1,NCells(3)
         DO icy=1,NCells(2)
            DO icx=1,NCells(1)
               IF (IsCellAffected(icx,icy,icz)==value) CALL AddCellContent1(icx,icy,icz,iunit)
            END DO
         END DO
      END DO
   END SUBROUTINE AddCellContent
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddCellContent1(icx,icy,icz,iunit)
      IMPLICIT NONE
      INTEGER :: icx,icy,icz,iunit
      TYPE(KMCAtomList), POINTER :: AL
      TYPE(KMCAtom), POINTER :: Atom
      CHARACTER(len=3) :: speciesname
      
      AL=>Cell(icx,icy,icz)%AL
      DO WHILE (ASSOCIATED(AL))
         Atom=>AL%Atom
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(Atom%Species))
         WRITE(UNIT=iunit,FMT='(a3,3F15.8,3L)') speciesname, &
            Atom%Coord,Atom%IsMoving,Atom%IsMoving,Atom%IsMoving
         AL=>AL%NextNeigh
      END DO
      
   END SUBROUTINE AddCellContent1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE LEKMCProcessSearch(nfiles)
   !searches for MD processes using the files created in LEKMCSetUpXYZ
      IMPLICIT NONE
      INTEGER :: nfiles
      
      
   END SUBROUTINE LEKMCProcessSearch
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateLEKMCProcessListMD()
   !Implicit converts the MD information to LEKMC processes
      IMPLICIT NONE
   END SUBROUTINE GenerateLEKMCProcessListMD
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   

END MODULE MDProcess
