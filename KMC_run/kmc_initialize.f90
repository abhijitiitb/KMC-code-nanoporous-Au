MODULE KMCInitialize !Level 7
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur

   USE KMC_VARIABLE_TYPE
   USE KMCUtilities
   USE KMCDbManipulate
   USE ModError
   USE KMCRecords
   USE ModFile
   USE KMCRelax
   USE GenerateProcess
   USE Ecuyer_random
   IMPLICIT NONE

   REAL :: CPUTimeRead
   TYPE(KMCProcess), POINTER :: prc

   CONTAINS

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE InitializeKMC()
      IMPLICIT NONE
      LOGICAL :: IsSnapOn

      CALL SYSTEM('mkdir ./StoreTADKMC')
      CALL SYSTEM('mkdir ./StoreTADKMC/SystemConfigs')
      CALL SYSTEM('mkdir ./StoreTADKMC/SystemProcesses')
      CALL SYSTEM('mkdir ./LatticeSnaps')
      CALL SYSTEM('mkdir ./CfgImg')
      CALL SYSTEM('mv Scrap.t Scrap.old.t')
      
      CPUTime=GetCurrentCPUTime(.TRUE.) !initialize CPU time variables for all processors
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Step 1: Initial steps...
      CALL ReadInput()
      WRITE(6,*) 'Ini>>Initialized KMC random seed :',idum
      CALL InitializeTaus88(idum)
      CALL ReadDataFromPreviousRun()
      CALL ReadCoordinates()
      CALL PrintSpeciesNotation()
      !CALL SetOrigin()
      !IsSnapOn=CheckIfSnapOn()
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Step 2: Build rest of the database
      CALL SetupCells()
      CALL RelaxFullKMCSystem()
      
      !IF (IsSnapOn) THEN
      !   WRITE(UNIT=*,FMT='("Ini>> Snapping atoms to lattice ...")',ADVANCE="NO")
      !ELSE
      !   WRITE(UNIT=*,FMT='("Ini>> Minimizing energy of system ...")',ADVANCE="NO")
      !END IF
      !CALL RelaxSystem(IsSnapOn)
      WRITE(6,*) " [DONE]"
      CALL PrintAtoms(PrintAtomsFormat,.FALSE.)
      !CALL SetUpFuzzy()
      !WARNING: make sure that if previous results read then
      !new subcells or official not created

      CALL SetupCfg()
      !CALL PrintCfgClasses()
      CALL GenerateInitialParametersReport()
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Step 3: Finish up with initializing
      CPUTime(6)=GetCurrentCPUTime(.FALSE.)
      WRITE(UnitScrap,*) ' '
      WRITE(UnitScrap,*) 'CPU time (min) after initialization=',CPUTime(6)
      WRITE(6,*) 'Ini>>Total # of configs:',NumberCfgTypes

      !this is for processes that have already been observed
      !CALL GenSystemRecord()
      !CALL PrintRecordedCfgSequence()

      prc=>KMC_PL
      DO WHILE (ASSOCIATED(prc))
         CALL AssignProcessFull(prc)
         prc=>prc%NextNeigh
      END DO
      
      CALL PrintProcessAssignments()
      CALL PrintAtomProcessAssignments()
      
      CALL GenerateProcessListFull() !generates the possible processes for the current system
      CALL GetKMCSumProcessRates()
   
      !TotRate=0._dp
      CALL GenSystemRecord()
      WRITE(6,*) 'Ini>>Completed initialization'
      WRITE(UnitTimeSnaps,*) "   KMCSteps       KMCTime(s)    CPUTime(s)   NumberCfgTypes  NumberProcesses  NAtomKMC"
      !WRITE(UnitTimeSnaps,*) "         KMC event |        KMC  time |"// &
      !   "        NOffic   |  NConfig  |   NAtoms |  NTADSimlPerEvent |   CPU time KMC"//&
      !   "   |    CPU time TAD   |     CPUTime All  |    #TADCall  |    Total Rate   |"//&
      !   " Rate Prefactor | Rate Barrier"

   END SUBROUTINE InitializeKMC
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReadInput()
      IMPLICIT NONE
      INTEGER :: j
      REAL :: PositionUncertainity
   
      WRITE(UNIT=6,FMT='("Ini>>Reading input file ...")',ADVANCE="NO") 
      SpeciesName=(/'Ag','Cu','Au','mo','pt','zn','pb','sn','fe','ru'/)

      CALL OpenFile(UnitScrap)
      CALL OpenFile(UnitTimeSnaps)
      CALL GetWelcomeMessage()
      WRITE(UnitScrap,*) '   Ini>>Reading file ',TRIM(FileKMCInputConfig)
  
      !Initialize time
      KMCTime=0._dp
      CPUTime=0.
      CPUStartValue=0.
      TotalCalls=0

      NULLIFY(RecordedCfg)
      CALL CreateRedundantCfg()

      CALL OpenFile(UnitInitialize)
      READ(UnitInitialize,*) NAtomKMC
      READ(UnitInitialize,*) NMoveKMC
      READ(UnitInitialize,*) NAtomPrint
      READ(UnitInitialize,*) NKMCBlock
      READ(UnitInitialize,*) NKMCStepsPerBlock
      READ(UnitInitialize,*) MaxKMCTime
      READ(UnitInitialize,*) KMCMinSumProcessRates
      READ(UnitInitialize,*) Temperature
      kBT=GetkBT(Temperature)
      
      READ(UnitInitialize,*) EnvCutOffRadius
      READ(UnitInitialize,*) BinSizeInEnv
      NumberBinsInEnv=CEILING(EnvCutOffRadius/BinSizeInEnv) !this is for the core
      EnvCutOffRadius=REAL(NumberBinsInEnv)*BinSizeInEnv
      EnvCutOffRadius2=EnvCutOffRadius**2

      READ(UnitInitialize,*) PositionUncertainity
      RadialTolerance=2.*PositionUncertainity

      READ(UnitInitialize,*) (KMCCellSize(j),j=1,3)
      NumberNeighborCellsInEnv=CEILING(EnvCutOffRadius/KMCCellSize(1))
      IF (NumberNeighborCellsInEnv<1 .OR. ANY(NCells/2>NumberNeighborCellsInEnv)) THEN
         WRITE(UnitScrap,*) "Ini>> NumberNeighborCellsInEnv=",NumberNeighborCellsInEnv
         ErrorMessage="Ini>>NumberNeighborCellsInEnv too small or too large"
         CALL ReportError()
      END IF
      Rad2Correction=(EnvCutOffRadius+RadialTolerance)**2 !we add a RadialTolerance to EnvCutOffRadius
      READ(UnitInitialize,*) BondOrderFlag
      WRITE(UnitScrap,*) "BondOrderFlag:",BondOrderFlag
   

      READ(UnitInitialize,*) EnabledTAD
      READ(UnitInitialize,*) EnabledDimer
      READ(UnitInitialize,*) EnabledUserKMCProcess
      READ(UnitInitialize,*) EnabledNEBWithUserKMCProcess
      IF (.NOT. EnabledTAD .AND. .NOT. EnabledDimer .AND. .NOT. EnabledUserKMCProcess) THEN
         WRITE(6,*) "Input file does not provide the KMC processes"
         WRITE(UnitScrap,*) "Input file does not provide the KMC processes"
         STOP
      END IF
      
      !IF (EnabledTAD) CALL ReadTADInput()
      IF (EnabledDimer) CALL ReadDimerInput()
      
      READ(UnitInitialize,*) idum

      READ(UnitInitialize,*) iFreqSnapshot
      READ(UnitInitialize,*) TimeFreqSnapshot
      READ(UnitInitialize,*) SnapshotType
      READ(UnitInitialize,*) KMCRelaxOption
      READ(UnitInitialize,*) UnitCellType,LatticeConst
      
      CALL InitializeUnitCell()
      
      READ(UnitInitialize,*) KMCTime
      READ(UnitInitialize,*) PrintAtomsFormat
      READ(UnitInitialize,*) IsConfigPrint
      READ(UnitInitialize,*) Restart
      READ(UnitInitialize,*) ReadProcessConfigList
      READ(UnitInitialize,*) ResetNoFirings
      READ(UnitInitialize,*) ShellBased
      READ(UnitInitialize,*) MakeDir
      READ(UnitInitialize,*) DontByPass
      READ(UnitInitialize,*) SlowDownFastProcess

      !READ(UnitInitialize,*) FuzzyOption
      !READ(UnitInitialize,*) FuzzyDelta
      !READ(UnitInitialize,*) FuzzyAlpha
      READ(UnitInitialize,*) NClusterProcs
      READ(UnitInitialize,*) IsDeposition
      IF (IsDeposition) THEN
         WRITE(6,*) " "
         WRITE(6,*) "Ini>> Deposition on"
      END IF
      WRITE(6,*) "NSpeciesType=",NSpeciesType
      READ(UnitInitialize,*) (DepositionRate(j), j=1,NSpeciesType)
      READ(UnitInitialize,*) RandomizeTADRun
      READ(UnitInitialize,*) IsIgnoreTAD
      READ(UnitInitialize,*) ForceNoPrint
      READ(UnitInitialize,*) IgnoreFailedTAD
      !READ(UnitInitialize,*) IgnoreRotationMatch
      READ(UnitInitialize,*) FillUpFirings
   
      IF (FillUpFirings .AND. ResetNoFirings) THEN
         WRITE(6,*) "Err>>  Cannot have FillUpFirings and ResetNoFirings @ same time"
         STOP
      END IF

      READ(UnitInitialize,*) BOCutoffForCfg
      READ(UnitInitialize,*) MaxCfgInDatabase
      READ(UnitInitialize,*) MaxPrcInDatabase

      CALL InitializeUtilities()
      IF (.NOT. IsDeposition) THEN
         DepositionRate=0._dp
         KMCTotDepositionRate=0._dp
      ELSE
         KMCTotDepositionRate=SUM(DepositionRate)
      END IF
      KMCSumProcessRates=KMCTotDepositionRate
      
      WRITE(UnitScrap,*) ' '
      CLOSE(UNIT=UnitInitialize)
   
      WRITE(UNIT=UnitScrap,FMT='("Ini>>Time initialized to: ",f10.3)') KMCTime
      KMCBlock=0
      KMCStep=0
      NumberCfgTypes=0
      WRITE(6,*) " [DONE]"
   END SUBROUTINE ReadInput
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReadTADInput()
      IMPLICIT NONE
      
      READ(UnitInitialize,*) TADLoTemperature
      READ(UnitInitialize,*) TADHiTemperature
      kBTLo=GetkBT(REAL(TADLoTemperature))
      kBTHi=GetkBT(REAL(TADHiTemperature))

      READ(UnitInitialize,*) TADMinPrefac
      READ(UnitInitialize,*) TADDelta
      CALL GetPTADCoeff()

      READ(UnitInitialize,*) TimeLowStop
      READ(UnitInitialize,*) TADmin_dTime
      READ(UnitInitialize,*) TADToDoOverDrive
      READ(UnitInitialize,*) TADProcSelecFlag
      READ(UnitInitialize,*) LazyTAD
      READ(UnitInitialize,*) TimeLazyTAD
      
      WRITE(UnitScrap,*) "TAD minprefac,delta,lowstop,dTime:", &
         TADMinPrefac,TADDelta,TimeLowStop,TADmin_dTime
   END SUBROUTINE ReadTADInput
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReadDimerInput()
      IMPLICIT NONE
   END SUBROUTINE ReadDimerInput
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateInitialParametersReport()
      IMPLICIT NONE
      
      WRITE(UNIT=*,FMT='("Ini>> Printing input values ..")',ADVANCE="NO")
      WRITE(UNIT=UnitScrap,FMT='("# atoms:",I5)') NAtomKMC
      WRITE(UNIT=UnitScrap,FMT='("# moving atoms:",I5)') NMoveKMC
      WRITE(UNIT=UnitScrap,FMT='("kmc blocks:",I10)') NKMCBlock
      WRITE(UNIT=UnitScrap,FMT='("kmc steps per block:",I10)') NKMCStepsPerBlock
      WRITE(UNIT=UnitScrap,FMT='("maximum time:",F10.3)') MaxKMCTime
      WRITE(UNIT=UnitScrap,FMT='("radial cutoff:",F10.3)') EnvCutOffRadius
      WRITE(UNIT=UnitScrap,FMT='("radial bin size:",F10.3)') BinSizeInEnv
      WRITE(UNIT=UnitScrap,FMT='("radial position tol:",F10.3)') RadialTolerance
      WRITE(UNIT=UnitScrap,FMT='("subcell size:",F10.3,F10.3,F10.3)') KMCCellSize
      WRITE(UNIT=UnitScrap,FMT='("defn of bond order used:",I5)') BondOrderFlag
      !IF (EnabledTAD) THEN
      !   WRITE(UNIT=UnitScrap,FMT='("system temperature:",I5)') TADLoTemperature
      !   WRITE(UNIT=UnitScrap,FMT='("TAD high temperature:",I5)') TADHiTemperature
      !   WRITE(UNIT=UnitScrap,FMT='("smallest KMC rate possible:",F20.8)') KMCMinSumProcessRates
      !   WRITE(UNIT=UnitScrap,FMT='("tad stop time:",F20.8)') TimeLowStop
      !   WRITE(UNIT=UnitScrap,FMT='("shortest tad time:",F20.8)') TADmin_dTime
      !END IF
      WRITE(UNIT=UnitScrap,FMT='("KMC relaxation type:",I5)') KMCRelaxOption
      WRITE(UNIT=UnitScrap,FMT='("unit cell type:",I5)') UnitCellType
      WRITE(UNIT=UnitScrap,FMT='("lattice const:",F10.3,F10.3,F10.3)') LatticeConst
      WRITE(UNIT=UnitScrap,FMT='("current kmc time:",F10.3)') KMCTime
      WRITE(UNIT=UnitScrap,FMT='("number of subcells:",3I5)') NCells
      WRITE(6,*) " [DONE]"
   END SUBROUTINE GenerateInitialParametersReport
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReadCoordinates() !reads the XYZ format
      IMPLICIT NONE
      INTEGER :: i=0,NAtomKMCcounter=0,NMoveKMCcounter=0
      INTEGER :: errorstatus,j,pos
      TYPE(KMCAtom), POINTER :: atom=>NULL()
      CHARACTER :: CrapChar
      CHARACTER(len=2) :: AtomicSymbol
      
      WRITE(UNIT=6,FMT='("Ini>> Reading input coordinates ...")')
      WRITE(UnitScrap,*) '   Ini>>Reading file ',TRIM(FileKMCInputConfig),' for input data'
      IF (Restart) THEN
         WRITE(6,*) "Ini>> Restarting old run"
         OPEN(UNIT=UnitKMCInputConfig,FILE=TRIM(FileStoreState)//".old")
      ELSE
         WRITE(6,*) "Ini>> Performing new run"
         CALL OpenFile(UnitKMCInputConfig)
      END IF
      
      READ(UnitKMCInputConfig,*) NAtomKMCcounter,NMoveKMCcounter
      IF (NAtomKMC/=NAtomKMCCounter .OR. NMoveKMC/=NMoveKMCCounter) THEN
         WRITE(6,*) "Error>> NAtomKMC or NMoveKMC values do not match in input files"
         WRITE(6,*) "NAtomKMC read:",NAtomKMC
         WRITE(6,*) "NMoveKMC read:",NMoveKMC
         WRITE(6,*) "NAtomKMC expected:",NAtomKMCCounter
         WRITE(6,*) "NMoveKMC expected:",NMoveKMCCounter
         STOP
      END IF
      NAtomKMC=0
      NMoveKMC=NMoveKMCCounter
      IF (NAtomKMCcounter<NMoveKMCcounter) THEN
         ErrorMessage="Ini>>Number moving atoms greater than number atoms"
         CALL ReportError()
      END IF

      !READ(UnitKMCInputConfig,*) CrapChar
      READ(UnitKMCInputConfig,*) (KMCBoxSize(j),j=1,3)
   
      NULLIFY(atom)
      NULLIFY(KMC_AL)
      DO i=1,NAtomKMCCounter
         CALL AddAtom(atom)
         IF (.NOT. ASSOCIATED(KMC_AL)) KMC_AL=>atom
         atom%Index=i
         READ(UnitKMCInputConfig,*) AtomicSymbol,atom%Coord,atom%IsMoving
         
         
         
         atom%Coord(1)=atom%Coord(1)+1.
         atom%Coord(2)=atom%Coord(2)+1.
         
         
         
         
         DO pos=1,10
            IF (TRIM(AtomicSymbol)==TRIM(SpeciesName(pos))) THEN
               atom%Species=pos
               EXIT
            END IF
         END DO
         !READ(UNIT=UnitKMCInputConfig,FMT='(3f20.8,i5)',IOSTAT=errorstatus) atom%Coord,atom%Species
         !IF (errorstatus>0) WRITE(6,*) 'Ini>>Unable to read coord for atom index:',i
         !IF (i>NMoveKMC) THEN
         !   atom%IsMoving=.FALSE.
         !ELSE
         !   atom%IsMoving=.TRUE.
         !END IF
      END DO

      IF (NAtomKMCcounter /= NAtomKMC) THEN
         ErrorMessage="Ini>>Incorrect # atoms in input file"
         CALL ReportError()
      END IF

      CLOSE(UNIT=UnitKMCInputConfig)
      
      WRITE(6,*) " [DONE]"
   END SUBROUTINE ReadCoordinates
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReadCoordinates1() !reads the clsman format
      IMPLICIT NONE
      INTEGER :: i=0,NAtomKMCcounter=0,NMoveKMCcounter=0
      INTEGER :: errorstatus,j
      TYPE(KMCAtom), POINTER :: atom=>NULL()
      CHARACTER :: CrapChar
      
      WRITE(UNIT=6,FMT='("Ini>> Reading input coordinates ...")')
      WRITE(UnitScrap,*) '   Ini>>Reading file ',TRIM(FileKMCInputConfig),' for input data'
      IF (Restart) THEN
         WRITE(6,*) "Ini>> Restarting old run"
         OPEN(UNIT=UnitKMCInputConfig,FILE=TRIM(FileStoreState)//".old")
      ELSE
         WRITE(6,*) "Ini>> Performing new run"
         CALL OpenFile(UnitKMCInputConfig)
      END IF
      
      READ(UnitKMCInputConfig,*) NAtomKMCcounter,NMoveKMCcounter
      IF (NAtomKMC/=NAtomKMCCounter .OR. NMoveKMC/=NMoveKMCCounter) THEN
         WRITE(6,*) "Error>> NAtomKMC or NMoveKMC values do not match in input files"
         WRITE(6,*) "NAtomKMC read:",NAtomKMC
         WRITE(6,*) "NMoveKMC read:",NMoveKMC
         WRITE(6,*) "NAtomKMC expected:",NAtomKMCCounter
         WRITE(6,*) "NMoveKMC expected:",NMoveKMCCounter
         STOP
      END IF
      NAtomKMC=0
      NMoveKMC=NMoveKMCCounter
      IF (NAtomKMCcounter<NMoveKMCcounter) THEN
         ErrorMessage="Ini>>Number moving atoms greater than number atoms"
         CALL ReportError()
      END IF

      READ(UnitKMCInputConfig,*) CrapChar
      READ(UNIT=UnitKMCInputConfig,FMT='(3f20.8)') (KMCBoxSize(j),j=1,3)
   
      NULLIFY(atom)
      NULLIFY(KMC_AL)
      DO i=1,NAtomKMCCounter
         CALL AddAtom(atom)
         IF (.NOT. ASSOCIATED(KMC_AL)) KMC_AL=>atom
         atom%Index=i
         READ(UNIT=UnitKMCInputConfig,FMT='(3f20.8,i5)',IOSTAT=errorstatus) atom%Coord,atom%Species
         IF (errorstatus>0) WRITE(6,*) 'Ini>>Unable to read coord for atom index:',i
         IF (i>NMoveKMC) THEN
            atom%IsMoving=.FALSE.
         ELSE
            atom%IsMoving=.TRUE.
         END IF
      END DO

      IF (NAtomKMCcounter /= NAtomKMC) THEN
         ErrorMessage="Ini>>Incorrect # atoms in input file"
         CALL ReportError()
      END IF

      CLOSE(UNIT=UnitKMCInputConfig)
      
      WRITE(6,*) " [DONE]"
   END SUBROUTINE ReadCoordinates1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetupCfg()
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom=>NULL()

      WRITE(UNIT=6,FMT='("Ini>> Setting neighbor lists ..")',ADVANCE="NO") 
      IF (ANY(KMCBoxSize/2.<EnvCutOffRadius)) THEN
         WRITE(6,*) "Cutoff for environment should be smaller than one-half boxsize"
         WRITE(UnitScrap,*) "Cutoff for environment should be smaller than one-half boxsize"
         STOP
      END IF
      atom=>KMC_AL
      !WRITE(6,*) "Printing configs ..."
      !WRITE(UnitScrap,*) "Printing configs ..."
      DO WHILE (ASSOCIATED(atom))
         IF (atom%IsMoving) CALL AddLocalCfg(atom)
         CALL AssignConfig(atom)
         !WRITE(UNIT=UnitScrap,FMT='(I4,3ES10.3,I6)') atom%Index,atom%Coord,atom%Cfg%ConfigurationIndex
         atom=>atom%NextNeigh
      END DO
      WRITE(6,*) " [DONE]"
   END SUBROUTINE SetupCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReadDataFromPreviousRun()
      IMPLICIT NONE

      IF (ReadProcessConfigList) THEN
         WRITE(6,*) "Ini>> Reading previously recorded configs and processes ..."
         CALL RecvrSystemRecord()
         CALL SYSTEM('cp '//TRIM(FileProcessListRec)//' '//TRIM(FileProcessListRec)//'.old')
         WRITE(6,*) " [DONE]"
      END IF
   END SUBROUTINE ReadDataFromPreviousRun
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetupCells()
      !size subcell is (CellSize(1)*CellSize(2)*CellSize(3))
      IMPLICIT NONE
      INTEGER :: indx
      TYPE (KMCAtom), POINTER :: atom=>NULL()
   
      WRITE(UNIT=*,FMT='("Ini>> Setting up subcells ..")',ADVANCE="NO")
      NCells=NINT(KMCBoxSize/KMCCellSize)
      ALLOCATE(Cell(NCells(1),NCells(2),NCells(3))); NAllottedCellAtom=NAllottedCellAtom+PRODUCT(NCells)
      ALLOCATE(IsCellAffected(NCells(1),NCells(2),NCells(3))); NAllottedIsCellAffected=&
        NAllottedIsCellAffected+PRODUCT(NCells)
      IsCellAffected=0
      
      WRITE(6,FMT='("Ini>>Allocated space to "// &
         "Subcell of size ",3i5)') (NCells(indx),indx=1,3)

      IF (ANY(NCells<8) .AND. (EnabledTAD .OR. EnabledDimer)) THEN
         ErrorMessage="Number of subcells in each direction " // &
            "should be >=8, reduce EnvCutOffRadius or increase lattice"
         CALL ReportError()
      END IF
      
      IF (ANY(CEILING(REAL(NCells)/2)<=NumberNeighborCellsInEnv)) THEN
         WRITE(6,*) "NumberNeighborCellsInEnv is too large"
         STOP
      END IF
   
      IF (ANY(ABS(KMCBoxSize-KMCCellSize*REAL(NCells))>0.01)) THEN
         WRITE(6,*) 'Ini>>CellSize:',KMCCellSize
         WRITE(6,*) "Actual system boxsize:",KMCBoxSize
         WRITE(6,*) 'Ini>>Projected boxsize:',KMCCellSize*REAL(NCells)
         WRITE(6,*) "Error>> Box size should be integer multiple of cell size"
         WRITE(UnitScrap,*) "Box size should be integer multiple of cell size"
         STOP
      END IF

      atom=>KMC_AL
      DO WHILE (ASSOCIATED(atom))
         CALL RefreshCell(atom)
         atom=>atom%NextNeigh
      END DO 
      WRITE(UnitScrap,*) '         ========================='
      WRITE(6,*) " [DONE]"
   END SUBROUTINE SetupCells
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE KMCInitialize

