MODULE AMD_VARIABLE_TYPE
   USE VARIABLE_TYPE
   USE KMC_VARIABLE_TYPE
   IMPLICIT NONE
   
   INTEGER, PARAMETER :: MaxSizeBEArray=100
   INTEGER, DIMENSION(MaxSizeBEArray) :: TargetNBasinEscapes=0
   REAL(dp), DIMENSION(MaxSizeBEArray) :: TargetStopTime=0._dp
   INTEGER :: NTarget=0

   LOGICAL :: TPMD=.FALSE. !whether TPMD method is running
   LOGICAL :: TPMDSuperbasin=.FALSE. !whether TPMD superbasin method is running
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE AMDStateContainer !used for global TAD simulations
   
      TYPE(AMDStateContainer), POINTER :: PrevNeigh=>NULL()
      TYPE(AMDStateContainer), POINTER :: NextNeigh=>NULL()
      TYPE(SuperbasinContainer), POINTER :: superbasin=>NULL()
      INTEGER :: Index=0,NVisits=0 !index for the state, # visits
      INTEGER :: NProcessType=0 !number of types of processes from the state
      INTEGER :: NMDTransitions=0 !number of transitions from the state
      INTEGER :: NMDFailed=0 !number of failed transitions from the state
      INTEGER :: WindowSize=5 !number of transitions that should have been seen
      !before any decision about rates estimates and validity time can be made
      INTEGER :: PreviousWindowEnd=0 !last transition # where a window was completed
      TYPE(AMDProcessContainer), POINTER :: Process=>NULL()
      TYPE(KMCProcess), POINTER :: LEKMCProcess=>NULL()
      REAL(dp) :: TotalProcessRate=0._dp
      REAL(dp) :: delta=0.1_dp !error associated with observing an unknown process
      REAL(dp) :: theta=0.5_dp !probability that at least one transition has been unknown
      REAL(dp) :: KMCValidTime=-1._dp !time for which KMC is valid
      REAL(dp) :: KMCTimeExpired=0._dp !time used up in KMC so far
      REAL(dp) :: RateEstimateUnknown=0._dp !upper bound for the sum of rates of unknown processes
      REAL(dp) :: MDTime=0._dp
      INTEGER, DIMENSION(MaxSizeBEArray) :: NBasinEscapes=0
      !REAL(dp), DIMENSION(MaxSizeBEArray) :: BasinEscapeTime=0._dp
      REAL(dp) :: BasinEscapeTime=0._dp
      REAL(dp) :: PotentialEnergy=0._dp
      
   END TYPE AMDStateContainer
   TYPE (AMDStateContainer), POINTER :: ListOfGlobalAMDStates=>NULL()
   INTEGER :: NGlobalAMDStates=0
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE AMDProcessContainer
      
      TYPE(AMDProcessContainer), POINTER :: PrevNeigh=>NULL()
      TYPE(AMDProcessContainer), POINTER :: NextNeigh=>NULL()
      
      TYPE(AMDStateContainer), POINTER :: InitialState=>NULL(),FinalState=>NULL()
      TYPE(AMDStateContainer), POINTER ::  TPMDActualInitial=>NULL() !in TPMD superbasin, the state from
         !where the escape occurs could be different from the pointer InitialState
      !TYPE(ListOfScalars), POINTER :: dTime=>NULL() !time required to MD to observe this transition
      REAL(dp), DIMENSION(MaxNTicks) :: Ticks=1.e6_dp !time for the next escape
      REAL(dp) :: Eact=0._dp,prefactor=0._dp,Rate=0._dp
      INTEGER :: Index=0,NTicks=0,AccruedTicks=0,Count=0 !accrued ticks count all ticks attempted while NTicks is stored number
      
      !TAD variables
      REAL(dp) :: TimeAccrued=0._dp !this is the time elapsed when the process was selected previous times
      REAL(dp) :: MDTimeAccrued=0._dp !this is total time accumulated doing pure MD
   END TYPE AMDProcessContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE SuperbasinstateContainer
      TYPE(SuperbasinstateContainer), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL() !contains the list of ALs in superbasin
      TYPE(AMDStateContainer), POINTER ::  AMDstate=>NULL()
   END TYPE SuperbasinstateContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE SuperbasinContainer !holds the superbasin
      TYPE(SuperbasinContainer), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL()
      TYPE(SuperbasinstateContainer), POINTER :: state=>NULL()
      INTEGER :: NStates=0,index=0
   END TYPE SuperbasinContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE(SuperbasinContainer), POINTER :: ListSuperbasins=>NULL()
   INTEGER :: NAMDSuperbasins=0,RunningIndexSuperbasins=0
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE TADParameters !parameters used for TAD calculation
      REAL(dp) :: TemperatureHigh=0._dp
      REAL(dp) :: dTemperature=0._dp,dTime=0._dp !for ramping up
      REAL(dp) :: minprefactor=0._dp,delta=0._dp
   END TYPE TADParameters
   !xxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE AMD_VARIABLE_TYPE
