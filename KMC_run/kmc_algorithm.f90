MODULE KMCAlgorithm
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
   !USE KMCUpdateSystemCarve
   USE KMCUpdateSystemDomainDecomposition
   IMPLICIT NONE
   !INTEGER :: AffectedCells1(3*SizeAffectedCells)=0
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE KMC()
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCProcess), POINTER :: prc
      
      !CALL Check() !checks if the current KMC database is self-consistent
      WRITE(6,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
      CALL SelectProcess(prc,atom)
      
      KMCStep=KMCStep+1
      IF (KMCStep==NKMCStepsPerBlock) THEN
         KMCBlock=KMCBlock+1
         KMCStep=0
      END IF
      
      CALL UpdateSystem(prc,atom)
      WRITE(6,*) "Reached here..."
      CALL FLUSH(6)
      STOP
   END SUBROUTINE KMC
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SelectProcess(prc,atom)
   !Selects the process from the list
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessList), POINTER :: ActiveProcess
      TYPE(KMCAtomList), POINTER :: AL1
      REAL(dp) :: KMCProcessRateSelected,prcrate
      REAL(dp) :: random_no
      LOGICAL :: Found
      INTEGER :: position,i
      
      !Find total rate
      CALL GetKMCSumProcessRates()
      !IF (KMCSumProcessRates<6.8e10 .OR. KMCSumProcessRates>7.0e10) THEN
         !WRITE(6,*) "KMCSumProcessRates are incorrect ..."
         WRITE(6,*) "... present value of KMCSumProcessRates is ",KMCSumProcessRates
         !WRITE(6,*) "...number of processes is",KMCSumProcessRates/34705467380.*4.
         !CALL PrintAtoms(PrintAtomsFormat,.FALSE.)
         !CALL GenSystemRecord()
         !STOP
      !END IF
      
      !Search for process
      random_no=taus88()
      KMCProcessRateSelected=KMCSumProcessRates*random_no
      !WRITE(6,*) "KMCSumProcessRates:",KMCSumProcessRates
      !WRITE(6,*) "KMCProcessRateSelected:",KMCProcessRateSelected
      ActiveProcess=>ActiveProcessList
      DO WHILE (ASSOCIATED(ActiveProcess))
         prc=>ActiveProcess%Process
         prcrate=prc%Rate*REAL(prc%NumberSubscribers,dp)
         !write(*,*) prc%CfgType%ConfigurationIndex,prc%Rate,prc%NumberSubscribers
         Found=prcrate>=KMCProcessRateSelected
         !write(*,*) "Found process:",Found
         IF (Found) EXIT
         KMCProcessRateSelected=KMCProcessRateSelected-prcrate
         ActiveProcess=>ActiveProcess%NextNeigh
      END DO
      
      IF (.NOT. ASSOCIATED(prc)) THEN
         WRITE(6,*) "Process could not be found"
CALL CheckKMC()
WRITE(6,*) "No issues found with the KMC databases -- yet prc is not associated"
         STOP
      END IF
      
      IF (KMCSumProcessRates<=0._dp) THEN
         WRITE(6,*) "Sum of process rates is less than or equal to zero"
CALL CheckKMC()
WRITE(6,*) "No issues found with the KMC databases -- yet sum of rates is zero"
         STOP
      ELSE
         KMCTime=KMCTime+taus88()/KMCSumProcessRates
      END IF
      
      position=CEILING(KMCProcessRateSelected/prc%Rate)
      !write(6,*) "Cfg of Process selected:",prc%CfgType%ConfigurationIndex
      !write(6,*) "Process rate:",prc%Rate
      !write(6,*) "Number of subscribers:",prc%NumberSubscribers
      !WRite(6,*) "Position selected:",position
      IF (prc%NumberSubscribers<position) THEN
         WRITE(6,*) "Number of subscribers is smaller than expected"
         WRITE(6,*) "Searching for atom position:",position
         WRITE(6,*) "Number of subscribers:",prc%NumberSubscribers
         STOP
      END IF
      
      AL1=>prc%AL
      DO i=1,position-1
         AL1=>AL1%NextNeigh
         IF (.NOT. ASSOCIATED(AL1)) THEN
            WRITE(6,*) "Unable to find atom for executing process"
            STOP
         END IF
      END DO 
      atom=>AL1%Atom
      !WRITE(6,*) "Done ... selected atom:",atom%Index,atom%Coord
   END SUBROUTINE SelectProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SelectProcess1(prc,atom)
   !Selects the process from the list
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCAtomList), POINTER :: AL1
      REAL(dp) :: KMCProcessRateSelected,prcrate
      REAL(dp) :: random_no
      LOGICAL :: Found
      INTEGER :: position,i
      
      !Find total rate
      CALL GetKMCSumProcessRates1()
      
      !Search for process
      random_no=taus88()
      KMCProcessRateSelected=KMCSumProcessRates*random_no
      prc=>KMC_PL
      DO WHILE (ASSOCIATED(prc))
         prcrate=prc%Rate*REAL(prc%NumberSubscribers,dp)
         Found=prcrate>=KMCProcessRateSelected
         IF (Found) EXIT
         KMCProcessRateSelected=KMCProcessRateSelected-prcrate
         prc=>prc%NextNeigh
      END DO
      
      IF (.NOT. ASSOCIATED(prc)) THEN
         WRITE(6,*) "Process could not be found"
         STOP
      END IF
      
      IF (KMCSumProcessRates<=0._dp) THEN
         WRITE(6,*) "Sum of process rates is less than or equal to zero"
         STOP
      ELSE
         KMCTime=KMCTime+taus88()/KMCSumProcessRates
      END IF
      
      position=CEILING(KMCProcessRateSelected/prc%Rate)
      IF (prc%NumberSubscribers<position) THEN
         WRITE(6,*) "Number of subscribers is smaller than expected"
         WRITE(6,*) "Searching for atom position:",position
         WRITE(6,*) "Number of subscribers:",prc%NumberSubscribers
         STOP
      END IF
      
      AL1=>prc%AL
      DO i=1,position-1
         AL1=>AL1%NextNeigh
         IF (.NOT. ASSOCIATED(AL1)) THEN
            WRITE(6,*) "Unable to find atom for executing process"
            STOP
         END IF
      END DO 
      atom=>AL1%Atom
   END SUBROUTINE SelectProcess1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !SUBROUTINE AssignAllProcessInAffectedCells()
   !   IMPLICIT NONE
   !   TYPE(KMCProcess), POINTER :: prc
      
   !   prc=>KMC_PL
   !   DO WHILE (ASSOCIATED(prc))
   !      CALL AssignProcessInAffectedCells(prc)
   !      prc=>prc%NextNeigh
   !   END DO
   !END SUBROUTINE AssignAllProcessInAffectedCells
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE KMCAlgorithm
