MODULE ModFile !Level 2
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
   USE KMC_VARIABLE_TYPE
   IMPLICIT NONE

   CONTAINS

!===================================================================
SUBROUTINE OpenFile(iUnit,OptFileName,OptOperation)
   IMPLICIT NONE
   INTEGER :: iUnit
   CHARACTER(len=30), OPTIONAL :: OptFileName
   CHARACTER, OPTIONAL :: OptOperation
   INTEGER :: errorstatus=0
   CHARACTER(len=30) :: FileName

   IF (PRESENT(OptFileName)) THEN
      FileName=OptFileName
   ELSE
      SELECT CASE (iUnit)
      CASE(UnitScrap); FileName=FileScrap
      CASE(UnitKMCInputConfig); FileName=FileKMCInputConfig
      CASE(UnitTADStart); FileName=FileTADStart
      CASE(UnitTADInput); FileName=FileTADInput
      CASE(UnitProcessList); FileName=FileProcessList
      CASE(UnitInitialize); FileName=FileInitialize
      CASE(UnitRelaxCoord); FileName=FileRelaxCoord
      CASE(UnitRelaxInput); FileName=FileRelaxInput
      CASE(UnitLatticeSnap); FileName=FileLatticeSnap
      CASE(UnitTADEnd); FileName=FileTADEnd
      CASE(UnitNEBResults); FileName=FileNEBResults
      CASE(UnitProcessListRec); FileName=FileProcessListRec
      CASE(UnitTADAddn); FileName=FileTADAddn
      CASE(UnitTimeSnaps); FileName=FileTimeSnaps
      CASE DEFAULT
         WRITE(*,*) 'Open>>iUnit=',iUnit
         WRITE(*,*) 'Open>>Incorrect iUnit passed to OpenFile.'
      END SELECT
   END IF
   
   IF (PRESENT(OptOperation)) THEN
      SELECT CASE (OptOperation)
      CASE ('r'); OPEN(UNIT=iUnit,FILE=FileName,STATUS='OLD',ACTION='READ',IOSTAT=errorstatus)
      CASE ('w'); OPEN(UNIT=iUnit,FILE=FileName,STATUS='REPLACE',ACTION='WRITE',IOSTAT=errorstatus)
      CASE ('a'); OPEN(UNIT=iUnit,FILE=FileName,STATUS='UNKNOWN',ACTION='READWRITE',IOSTAT=errorstatus)
      CASE DEFAULT
      END SELECT
   ELSE
      SELECT CASE (iUnit)
      CASE(UnitScrap,UnitTADStart,UnitTADInput) !to make compatible with prev version
         OPEN(UNIT=iUnit,FILE=FileName,STATUS='REPLACE',IOSTAT=errorstatus)
      CASE(UnitKMCInputConfig,UnitProcessList,UnitInitialize)
         OPEN(UNIT=iUnit,FILE=FileName,STATUS='OLD',IOSTAT=errorstatus)
      CASE DEFAULT
         OPEN(UNIT=iUnit,FILE=FileName,IOSTAT=errorstatus)
      END SELECT
   END IF
  
   IF (iUnit/=UnitScrap .AND. errorstatus/=0) THEN 
      WRITE(*,*) 'Open>>Error opening file: '//FileName
      STOP
   END IF
  ! WRITE(*,*) 'Open>> File '//FileName
END SUBROUTINE OpenFile
!===========================================================

END MODULE ModFile
