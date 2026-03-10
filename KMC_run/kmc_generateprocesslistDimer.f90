MODULE DimerProcess
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
   USE KMC_VARIABLE_TYPE
   IMPLICIT NONE
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateDimerProcessListFull()
   END SUBROUTINE GenerateDimerProcessListFull
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateDimerProcessListCell(CellIndex)
      IMPLICIT NONE
      INTEGER :: CellIndex(3)
   END SUBROUTINE GenerateDimerProcessListCell
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE DimerProcess
