MODULE ModError !Level 6
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
   USE KMC_VARIABLE_TYPE
   USE KMCUtilities, ONLY : GetCurrentCPUTime,PrintAtoms
   IMPLICIT NONE
   CHARACTER(len=150) :: ErrorMessage

   CONTAINS

!===============================================
SUBROUTINE ReportError()
   IMPLICIT NONE
   REAL :: CPUTimeRead
   
   WRITE(6,*) 'Err>>reporting error/....'
   WRITE(6,*) ErrorMessage
   WRITE(UnitScrap,*) ErrorMessage
   WRITE(6,*) 'Err>>Error reported in scrap file'
   
   WRITE(UnitScrap,*) 'Err>>TAD calls statistics'
   WRITE(UnitScrap,*) 'Err>>CPU time:',CPUTime(2)
   WRITE(UnitScrap,*) 'Err>>#calls:',TotalCalls(2)
   WRITE(UnitScrap,*) 'Err>>Total CPU statistics'
   CPUTimeRead=GetCurrentCPUTime(.FALSE.)
   CPUTime(6)=CPUTimeRead
   WRITE(UnitScrap,*) 'Err>>CPU time:',CPUTime(6)
   CLOSE(UnitScrap)
      
   WRITE(UnitLatticeSnap,*) 'Err>>Saving current lattice snap before closing'
   IF (KMCStep>1) CALL PrintAtoms(PrintAtomsFormat)
   CLOSE(UnitLatticeSnap)
   !CALL GenSystemRecord !prints the config and process list
   STOP
END SUBROUTINE ReportError
!===============================================
SUBROUTINE CheckForCorrectness(errorstatus,val)
   IMPLICIT NONE
   INTEGER :: errorstatus,val
   
   SELECT CASE (val)
   CASE(1)
      IF (errorstatus/=0) THEN 
         ErrorMessage="Ini>>Unable to read the variable NAtomKMC from file "//FileInitialize
         CALL ReportError()
      END IF
      WRITE(UnitScrap,FMT='("Ini>>Total number of atoms present: ",i5)') NAtomKMC
      IF (NAtomKMC<=0) THEN 
         ErrorMessage="Ini>>Total number of atoms in the KMC lattice has to be greater than zero"
         CALL ReportError()
      END IF
   CASE(2)
      IF (errorstatus/=0) THEN 
         ErrorMessage='Ini>>Unable to read the variable NMoveKMC from file '//FileInitialize
         CALL ReportError()
      END IF
      WRITE(UnitScrap,FMT='("Ini>>Total number of moving atoms present: ",i5)') NMoveKMC
      IF (NMoveKMC<=0) THEN 
         ErrorMessage='Ini>>Total number of atoms in the KMC lattice that can move has to be greater than zero'
         CALL ReportError()
      END IF
   CASE(3)
      IF (errorstatus/=0) THEN 
         ErrorMessage='Ini>>Unable to read the variable NKMBlock from file '//FileInitialize
         CALL ReportError()
      END IF
      WRITE(UnitScrap,FMT='("Ini>>Total number of KMC events: ",i10)') NKMCBlock
      IF (NKMCBlock<=0) THEN 
         ErrorMessage='Ini>>Total number of KMC events has to be greater than zero'
         CALL ReportError()
      END IF
   CASE(4)
      IF (errorstatus/=0) THEN 
         ErrorMessage='Ini>>Unable to read the variable MaxKMCTime from file '//FileInitialize
         CALL ReportError()
      END IF
      WRITE(UnitScrap,FMT='("Ini>>Total KMC time accessed:",f10.3)') MaxKMCTime
      IF (MaxKMCTime<=0._dp) THEN 
         ErrorMessage='Ini>>Total KMC time to be accessed has to be greater than zero'
         CALL ReportError()
      END IF
   CASE(5)
      IF (errorstatus/=0) THEN 
         ErrorMessage='Ini>>Unable to read the variable EnvCutoffRadius from file '//FileInitialize
         CALL ReportError()
      END IF
      WRITE(UnitScrap,FMT='("Ini>>Cut off radius for radial bins:",f10.3)') EnvCutoffRadius
      IF (EnvCutoffRadius<=0.) THEN 
         ErrorMessage='Ini>>Cutoff radius for radial bins has to be greater than zero'
         CALL ReportError()
      END IF
   CASE(6)
      IF (errorstatus/=0) THEN 
         ErrorMessage='Ini>>Unable to read the variable BinSizeInEnv from file '//FileInitialize
         CALL ReportError()
      END IF
      WRITE(UNIT=UnitScrap,FMT='("Ini>>Total number of bins along radial direction:",i5)') NumberBinsInEnv 
      WRITE(UNIT=UnitScrap,FMT='("Ini>>Radial bin size forced by the code:",f10.3)') BinSizeInEnv
      IF (NumberBinsInEnv<=0) THEN 
         ErrorMessage='Ini>>Total number of bins along radial direction has to be greater than zero'
         CALL ReportError()
      END IF
   CASE(7)
      IF (errorstatus/=0) THEN 
         ErrorMessage='Ini>>Unable to read the variable AtomPositionUncertainity from file '//FileInitialize
         CALL ReportError()
      END IF
      WRITE(UnitScrap,FMT='("Ini>>Tolerance in radial positions (error bar):",f10.3)') RadialTolerance
      IF (RadialTolerance<=0._dp) THEN 
         ErrorMessage='Ini>>Tolerance in radial positions has to be greater than zero'
         CALL ReportError()
      END IF
   CASE(8)
      IF (errorstatus/=0) THEN 
         ErrorMessage='Ini>>Unable to read the variable SubcellFactor from file '//FileInitialize
         CALL ReportError()
      END IF
   CASE(9)
      IF (errorstatus/=0) THEN 
         ErrorMessage='Ini>>Unable to read the variable Temperature from file '//FileInitialize
         CALL ReportError()
      END IF
      IF (errorstatus/=0) THEN 
         ErrorMessage='Ini>>Unable to read the variable TemperatureTADHigh from file '//FileInitialize
         CALL ReportError()
      END IF
      WRITE(UnitScrap,FMT='("Ini>>Simulation temperature:",i5)') TADLoTemperature
    WRITE(UnitScrap,FMT='("Ini>>TAD high temperature:",i5)') TADHiTemperature
    IF (TADLoTemperature<=0. .OR. TADHiTemperature<=0.) THEN 
       ErrorMessage='Ini>>Both simulation temperature and high temperature used for TAD should be greater than zero'
       CALL ReportError()
    END IF
    IF (TADLoTemperature>TADHiTemperature) THEN 
       ErrorMessage='Ini>>Simulation temperature should be smaller than TAD high temperature'
       CALL ReportError()
    END IF
 CASE(10)
    IF (KMCMinSumProcessRates<1.e-3) THEN 
       ErrorMessage='Ini>>Select a value of KMCMinSumProcessRates>0'
       CALL ReportError()
    END IF
 CASE(11)
    IF (TADMinPrefac<=1.e4) THEN
       ErrorMessage='Ini>>Min prefactor for TAD has to be set postive, but > 1e4'
       CALL ReportError()
    END IF
 CASE(12)
    IF (TADDelta<=1.e-7) THEN
       ErrorMessage='Ini>>Delta for TAD has to be set to positive, but > 1e-7'
       CALL ReportError()
    END IF
 CASE(13)
    IF (TimeLowStop<=0._dp) THEN
       ErrorMessage='Ini>>Low T stop time for TAD < 0.'
       CALL ReportError()
    END IF
 CASE(14)
    IF (errorstatus/=0) THEN 
       ErrorMessage='Ini>>Unable to read the variable CellSize(1:3) from file '//FileKMCInputConfig
       CALL ReportError()
    END IF
    WRITE(UnitScrap,FMT='("Ini>>KMC lattice size:",3f10.3)') KMCCellSize
    IF (ANY(KMCCellSize<=0._dp)) THEN 
       ErrorMessage='Ini>>KMC lattice dimension(s) has to be greater than zero'
       CALL ReportError()
    END IF
 CASE(15)
    IF (errorstatus/=0) THEN
       WRITE(UnitScrap,*) 'Ini>>Reading details of atom ',errorstatus
       ErrorMessage='Ini>>Error in reading line of '//FileKMCInputConfig
       CALL ReportError()
    END IF
 CASE(16)
    IF (NAtomKMC/=errorstatus) THEN
       WRITE(UnitScrap,*) 'Ini>>Number of atoms coordinates read: ',errorstatus
       ErrorMessage='Ini>>Number of atoms in input file does not match number suggested by user'
       CALL ReportError()
    END IF
 END SELECT
END SUBROUTINE CheckForCorrectness
!*********************************************************
SUBROUTINE CheckAllSitesAssociated()
  IMPLICIT NONE
  TYPE(KMCAtom), POINTER :: atom
  LOGICAL :: IsNotAssociated
  
  IsNotAssociated=.FALSE.
  atom=>KMC_AL
  DO WHILE (ASSOCIATED(atom))
     IF (.NOT. ASSOCIATED(atom%Cfg)) THEN
        WRITE(6,*) 'Err>>Site not pointing to any cfg:',atom%Index
        IsNotAssociated=.TRUE.
     END IF
     atom=>atom%NextNeigh
   END DO
   IF (IsNotAssociated) STOP
 END SUBROUTINE CheckAllSitesAssociated
 !*********************************************************
END MODULE ModError
