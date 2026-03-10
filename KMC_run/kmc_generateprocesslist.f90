MODULE GenerateProcess
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
   !Module is used to generate process list to be used by KMC
   USE KMC_VARIABLE_TYPE
   USE AMDProcess
   USE LEKMCMD
   USE DimerProcess
   USE UserProcess
   
   IMPLICIT NONE
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateProcessListFull()
   !Last checked: Jan 09, 2011
      IMPLICIT NONE
      
      IF (EnabledTAD) CALL GenerateTADProcessListFull()
      IF (EnabledDimer) CALL GenerateDimerProcessListFull()
      IF (EnabledUserKMCProcess) CALL GenerateUserProcessListFull()
   END SUBROUTINE GenerateProcessListFull
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateProcessListCell(CellIndex)
   !Last checked: Jan 09, 2011
      IMPLICIT NONE
      INTEGER :: CellIndex(3)
      
      !IF (EnabledTAD) CALL GenerateTADProcessListCell(CellIndex)
      IF (EnabledDimer) CALL GenerateDimerProcessListCell(CellIndex)
      IF (EnabledUserKMCProcess) CALL GenerateUserProcessListCell(CellIndex)
   END SUBROUTINE GenerateProcessListCell
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE GenerateProcess
