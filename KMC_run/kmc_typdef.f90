MODULE KMC_VARIABLE_TYPE
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
   USE VARIABLE_TYPE, ONLY : I4B,I2B,I1B,sp,dp,PI
   IMPLICIT NONE
   !INTEGER, PARAMETER :: i4b = SELECTED_INT_KIND(9)
   !INTEGER, PARAMETER :: i2b = SELECTED_INT_KIND(4)
   !INTEGER, PARAMETER :: i1b = SELECTED_INT_KIND(2)
   !INTEGER, PARAMETER :: sp = KIND(1.0)
   !INTEGER, PARAMETER :: dp = KIND(1.0d0)

   INTEGER, PARAMETER :: NSpeciesType=2 !# atomic species present (excluding vacancies)

   ! File unit numbers and names
   INTEGER, PARAMETER :: UnitScrap=101 !for reporting current status of the code (o/p)
   INTEGER, PARAMETER :: UnitKMCInputConfig=102 !for reading initial config (i/p)
   INTEGER, PARAMETER :: UnitTADStart=103
   INTEGER, PARAMETER :: UnitTADInput=104
   INTEGER, PARAMETER :: UnitInitialize=105
   INTEGER, PARAMETER :: UnitLatticeSnap=106
   INTEGER, PARAMETER :: UnitRelaxCoord=107
   INTEGER, PARAMETER :: UnitRelaxInput=108
   INTEGER, PARAMETER :: UnitTADEnd=109
   !INTEGER, PARAMETER :: UnitTADEnd=110 (reserved)
   INTEGER, PARAMETER :: UnitNEBResults=111
   INTEGER, PARAMETER :: UnitProcessList=112
   INTEGER, PARAMETER :: UnitProcessListRec=114
   INTEGER, PARAMETER :: UnitTADAddn=115
   INTEGER, PARAMETER :: UnitBlok=116
   INTEGER, PARAMETER :: UnitTimeSnaps=117
   INTEGER, PARAMETER :: UnitStoreState=118
   INTEGER, PARAMETER :: SizeAffectedCells=4000 !maximum # atoms that can get affected in adaptive KMC

   CHARACTER(len=30), PARAMETER :: FileScrap="Scrap.t"
   CHARACTER(len=30), PARAMETER :: FileKMCInputConfig="InputCoord.xyz"  !"StepAgInput" !
   CHARACTER(len=30), PARAMETER :: FileInitialize="input.initial"
   CHARACTER(len=30), PARAMETER :: FileTADStart="TAD.start"
   CHARACTER(len=30), PARAMETER :: FileTADInput="TAD.input"
   CHARACTER(len=3), PARAMETER :: Casename="TAD" !storage depends on the official #
   !CHARACTER(len=30), PARAMETER :: FileLatticeSnap="Lattice.xyz"
   CHARACTER(len=30), PARAMETER :: FileProcessList="TAD.output"
   CHARACTER(len=30), PARAMETER :: FileRelaxCoord="Rlx.start"
   CHARACTER(len=30), PARAMETER :: FileRelaxInput="TAD.input"
   CHARACTER(len=30), PARAMETER :: FileRelaxedPositions="TAD.output"
   CHARACTER(len=30), PARAMETER :: FileTADEnd="TAD.end"
   CHARACTER(len=30), PARAMETER :: FileNEBResults="TAD.neb" 
   CHARACTER(len=30), PARAMETER :: FileProcessListRec="./StoreTADKMC/SystemDetails"
   CHARACTER(len=30), PARAMETER :: FileTADAddn="TAD.addn"
   CHARACTER(len=30), PARAMETER :: FileStoreState="./StoreTADKMC/LastSavedState"
   CHARACTER(len=30), PARAMETER :: FileTimeSnaps="Snapshots.t"
   
   !REAL(dp), PARAMETER :: PI=3.141592654_dp

   INTEGER, PARAMETER :: MaxKMCProcessAtoms=10000
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE CellAtom
      INTEGER :: NAtomsInCell=0,NMoveInCell(3)=0
      TYPE (KMCAtomList), POINTER :: AL=>NULL()
   END TYPE CellAtom
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE (CellAtom), ALLOCATABLE :: Cell(:,:,:)
   INTEGER, DIMENSION(:,:,:), POINTER :: IsCellAffected !TRUE when the cell is affected due to a KMC process being selected
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE KMCAtom
      TYPE (KMCAtom), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL()
      REAL, DIMENSION(3) :: Coord=0.
      INTEGER :: Species=0,Index=0
      INTEGER :: CellIndex(3)=0
      REAL :: BondOrder=0.
      REAL(dp) :: CatalogTimeAvailable=0._dp
      !Why need CatalogTimeAvailable?
      !One of the major issues with self-learning KMC methods is to know when it is time to search
      !for new processes. If all processes are known for a material then each time the system 
      !moves to a new state one can identify the processes that get more/less activated, and the atoms 
      !that activate the process. Each atom can be linked to a process -- hence if an atom participates
      !in a process, we make the processes involving this atom to be less active by 1 - this has to be 
      !done carefully since many atoms are participating in a process. After figuring the new atom 
      !configurations (locally), we search for processes that are activated (again locally).
      !
      TYPE(Environment), POINTER :: Env=>NULL()
      TYPE(Configuration), POINTER :: Cfg=>NULL()
      TYPE(KMCProcessSubscriptionInfo), POINTER :: PrcSubscriberInfo=>NULL()
      LOGICAL :: IsMoving=.TRUE.
   END TYPE KMCAtom
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE(KMCAtom), POINTER :: KMC_AL=>NULL() !this contains the KMC atom
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE KMCProcessSubscriptionInfo
      TYPE(KMCProcessSubscriptionInfo), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL()
      TYPE(KMCProcess), POINTER :: Process=>NULL()
      LOGICAL :: IsActive=.FALSE.
      TYPE(KMCAtom), POINTER :: ActiveAtom=>NULL()
   END TYPE KMCProcessSubscriptionInfo
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE KMCProcess
      TYPE(KMCProcess), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL() !Linked list
      TYPE(KMCProcessAtom), POINTER :: ReactantAL=>NULL() !participating atoms
      TYPE(Configuration), POINTER :: CfgType=>NULL() !which configuration process belongs
      TYPE(KMCProcessList), POINTER :: APLPosition=>NULL()
      TYPE(KMCAtomList), POINTER :: AL=>NULL() !subscribers to this process
      INTEGER :: NumberProcessAtoms=0,NumberSubscribers=0,Index=0
      INTEGER :: FiringCount=0 !# times fired so far
      INTEGER :: MaxFiringAllowed=0 !max # times firing is allowed
      REAL(dp) :: Rate=0._dp,Frequency=0._dp,FwdBarrier=0._dp
      !REAL :: MaxPrcAtomDisplacement=0.
      INTEGER :: RecordNumber=0 !permanent record number
      INTEGER :: LastUpdate=0 !last KMC iteration where the process was accessed
      INTEGER :: FirstAdded=0 !KMC iteration where the process was added
   END TYPE KMCProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE(KMCProcess), POINTER :: KMC_PL=>NULL() !this contains the process list
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE KMCProcessAtom
      TYPE(KMCProcessAtom), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL()
      INTEGER :: CfgIndex=0,Species=0
      REAL, DIMENSION(3) :: RelInitialPos=0.,RelSaddlePos=0.,RelFinalPos=0.
   END TYPE KMCProcessAtom
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE Environment
      TYPE(Histogram), POINTER :: ShortRangeConfig=>NULL()
      !TYPE(Long), POINTER :: LongRangeConfig=NULL()
   END TYPE Environment
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE Histogram
      TYPE(Bin), POINTER :: BinInfo=>NULL()
      !TYPE(BinAtom), POINTER :: RedundantNeighbors=>NULL()
      INTEGER, DIMENSION(NSpeciesType) :: Population=0
   END TYPE Histogram
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE Bin
      TYPE (Bin), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL()
      TYPE (BinAtom), POINTER :: AL=>NULL()
      INTEGER :: BinNumber=0
      INTEGER, DIMENSION(NSpeciesType) :: InBin=0,BinEdgeP=0,BinEdgeN=0
   END TYPE Bin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE BinAtom
      TYPE(BinAtom), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL()
      TYPE(KMCAtom), POINTER :: Atom=>NULL()
      INTEGER :: Species=0
      REAL :: distance=0.,RelCoord(3)=0.
   END TYPE BinAtom
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE KMCAtomList
      TYPE(KMCAtomList), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL()
      TYPE(KMCAtom), POINTER :: Atom=>NULL()
   END TYPE KMCAtomList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE Configuration
      TYPE(Configuration), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL()
      TYPE(KMCAtomList), POINTER :: AtomsWithCfg=>NULL()
      TYPE(Environment), POINTER :: Env=>NULL()
      TYPE(KMCProcessList), POINTER :: PL=>NULL() !processes associated with this cfg
      INTEGER :: ConfigurationIndex=0,Species=0
      INTEGER :: NAtomsWithCfg=0
      INTEGER :: RecFileNumber=0
      REAL :: BondOrder=0.
   END TYPE Configuration
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE(Configuration), POINTER :: RecordedCfg=>NULL(), RedundantCfg=>NULL()
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE ConfigurationList
      TYPE(ConfigurationList), POINTER :: PrevNeigh=>NULL(),NextNeigh=>NULL()
      TYPE(Configuration), POINTER :: Cfg=>NULL()
   END TYPE ConfigurationList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE(ConfigurationList), POINTER :: NewCfgList=>NULL()
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE KMCProcessList
      TYPE(KMCProcessList), POINTER :: PrevNeigh=>NULL(),NextNeigh=>NULL()
      TYPE(KMCProcess), POINTER :: Process=>NULL()
   END TYPE KMCProcessList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE(KMCProcessList), POINTER :: NewProcessList=>NULL() !newly created processes
   TYPE(KMCProcessList), POINTER :: ActiveProcessList=>NULL() !active processes
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   !Variables used in simulation
   
   INTEGER :: NAtomKMC=0,NMoveKMC=0,NCells(3)=0,NAtomPrint=0
   INTEGER :: NAffectedCells,AffectedCells(3*SizeAffectedCells) !upto 100 affected cells during a KMC process can be stored
      
   INTEGER :: KMCMode=0 !0 is off-lattice with no relax, 1 is off-lattice with relax, 2 is on-lattice
   INTEGER :: KMCBlock=0,KMCStep=0 !tracks number of iterations for KMC
   INTEGER :: NKMCBlock=0,NKMCStepsPerBlock=0 !maximum iterations for KMC
   
   CHARACTER(len=3) :: SpeciesName(10)
   REAL :: Temperature=0. !system temperature
      
   REAL(dp) :: KMCTime=0._dp,MaxKMCTime=0._dp,KMCSumProcessRates=0._dp,KMCTotDepositionRate=0._dp
   REAL(dp) :: KMCMinSumProcessRates=0._dp,TotalRateSnapshot=0._dp
   
   LOGICAL :: IsDeposition=.FALSE.,EnabledTAD=.FALSE.,EnabledDimer=.FALSE.,EnabledUserKMCProcess=.FALSE.
   LOGICAL :: EnabledNEBWithUserKMCProcess=.FALSE.
   
   REAL :: KMCCellSize(3)=0.,KMCBoxSize(3)=0.
   REAL :: EnvCutOffRadius=0.,EnvCutOffRadius2=0. !radius for local environment
   REAL :: BinSizeInEnv=0. !bin width
   REAL :: PositionUncertainty=0. !uncertainty in the position of an atom
   REAL :: RadialTolerance=0.,Rad2Correction=0.
      
   INTEGER :: NumberBinsInEnv=0
   INTEGER :: NumberNeighborCellsInEnv=0
   INTEGER :: NumberCfgTypes=0,NumberProcesses=0,NumberActiveProcesses=0,NFailedKMCMoves=0
   INTEGER :: MaxCfgInDatabase=0,MaxPrcInDatabase=0
   INTEGER :: SnapshotType=1,BondOrderFlag=2,KMCRelaxOption=2,UnitCellType=0
   INTEGER :: iseed=2222,idum=2222
   REAL :: LatticeConst(3)=0.,KMCMaxLatticeSpacing=0.
   REAL :: BOCutoffForCfg=1000.
   REAL :: DepositionRate(NSpeciesType)
   
   !TAD parameters
   INTEGER :: NClusterProcs,NTADSimlPerEvent
   REAL :: TADLoTemperature=0.,TADHiTemperature=0.
   REAL :: TADpmi=0.,TADpma=0.,TADMinPrefac=0._dp,TADDelta=0.,TADAlpha=0.,TADNuMin=0.
   INTEGER :: TADProcSelecFlag=1
   REAL(dp) :: TimeLowStop=0._dp,TimeLazyTAD=10._dp,TADmin_dTime=0._dp
   REAL :: kBT=0.,kBThi=0.,kBTlo=0.
   LOGICAL :: TADToDoOverDrive=.FALSE.,ShellBased=.TRUE.,DontBypass=.FALSE.
   LOGICAL :: LazyTAD=.FALSE.,RandomizeTADRun=.FALSE.,IsIgnoreTAD=.FALSE.
   LOGICAL :: SlowDownFastProcess=.FALSE.,IgnoreFailedTAD=.FALSE.
   LOGICAL :: FillUpFirings=.TRUE.

   !Print variables
   INTEGER :: SnapCounter=0
   INTEGER :: iFreqSnapshot=100000
   REAL(dp) :: TimeFreqSnapshot=1.e10,PrnRatePrefac=-10._dp,PrnRateBarrier=-10._dp
   LOGICAL :: IsConfigPrint=.TRUE.,ForceNoPrint=.FALSE.
   CHARACTER :: PrintAtomsFormat='x'
   CHARACTER(len=30) :: FileLatticeSnap="Lattice"
   
   !CPU and counts
   INTEGER :: TotalCalls(7)=0
   REAL :: CPUTime(7),CPUStartValue(7)=0.,CPUAnalysis(6)=0.
   REAL :: ProcessorTime(100) !assuming 100 procs
   !Indices:     1. RecordedConfigsFile, 2. TADRunCalls, 3. RelaxCalls
   !4. NEB calls, 5. KMC calls, 6. Entire siml, 7. Print lattice calls
   
   !Restart variables
   LOGICAL :: ReadProcessConfigList=.TRUE.,ResetNoFirings=.TRUE.,MakeDir=.FALSE.
   LOGICAL :: Restart=.FALSE.
   
   !variables for storing KMC memory requirements
   INTEGER :: NAllottedCellAtom=0 
   INTEGER :: NAllottedIsCellAffected=0
   INTEGER :: NAllottedKMCAtom=0
   INTEGER :: NAllottedKMCProcessSubscriptionInfo=0
   INTEGER :: NAllottedKMCProcess=0
   INTEGER :: NAllottedKMCProcessAtom=0
   INTEGER :: NAllottedKMCProcessList=0
   INTEGER :: NAllottedEnvironment=0
   INTEGER :: NAllottedHistogram=0
   INTEGER :: NAllottedBin=0
   INTEGER :: NAllottedBinAtom=0
   INTEGER :: NAllottedKMCAtomList=0
   INTEGER :: NAllottedConfiguration=0
   INTEGER :: NAllottedConfigurationList=0
   INTEGER :: NAllottedOthers=0
   
END MODULE KMC_VARIABLE_TYPE
!If restarting then
! Manually modify Process.rec according to your needs
! Make sure input.coord is correct
! Check if # firings need to be reset and if old process file needs to be read
! Make copies of previous files such as snaps, process and snapshots

