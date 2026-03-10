MODULE AMDProcess
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
   USE KMC_VARIABLE_TYPE
   USE KMCProcessSearchParallel
   USE AMD_VARIABLE_TYPE
   USE TPMDModule 
   USE Ecuyer_random
   USE amd_db_manipulate
   IMPLICIT NONE
   
   CONTAINS
   INCLUDE "kmc_generateprocesslistAMD.MLE-RateEstimation.f90" !subroutines used for Vijesh's BES with MLE for rates
   !estimation of unknown rates is now being performed using timescale spectrum
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ParallelReplicaMDKMC(state,NKMCMove,StopTime,Time,nsearch,Temperature, &
      imode,StoreLEKMCProcess,DetectTransitionFrequency,TADParam,WriteTrajectoryFile)
      !Performs a parallel replica + KMC
      !state is the current state of the system -- at the end the final state is returned
      !NKMCMove - number of KMC moves to be performed
      !Time - total KMC time elapsed
      !nsearch - number of search trajectories attempted for each state
      !DetectTransitionFrequency - frequency of checking for a process in MD
      !WriteTrajectoryFile - writes the states visited by BCMD (it is a log file)
      !imode - running mode for ParallelReplicaMDKMC
      !imode 1: Perform nsearches for a single state and report back w/o a validity time
      !imode 2: Generate validity time
      !imode 3: parallel replica MD
      !imode 4: temperature accelerated dynamics - Sorenson
      !imode 5: temperature accelerated dynamics - eTAD - Montalenti
      !imode 6: temperature accelerated dynamics - pTAD - Chatterjee
      !imode 7: hyperdynamics - bond-boost method
      
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state,statefinal
      TYPE(AMDProcessContainer), POINTER :: process,process1
      TYPE(TADParameters), POINTER, OPTIONAL :: TADParam
      TYPE(SystemContainer), POINTER :: AL1
      INTEGER :: NKMCMove,iKMCMove,nsearch,DetectTransitionFrequency,StoppingOption,nsuccess
      INTEGER :: KMCIterations,imode,isearch,ObtainRate,i
      REAL(dp) :: Time,Time0,Temperature,MDStopTime,dt,rand1,rand2,ProcessRate,TotalProcessRate,RandProcessRate
      REAL(dp) :: TotalMDTime,CollectedMDTime,dtparrep,cputime
      REAL(dp), OPTIONAL :: StopTime
      REAL(dp) :: TemperatureHigh,delta,minprefactor,tLshort,tHStop,dTime,cputime1
      CHARACTER(len=100), OPTIONAL :: WriteTrajectoryFile
      CHARACTER(len=100) :: filename
      CHARACTER(len=NINT2CHAR) :: charstateindx,charisearch
      LOGICAL :: NotSatisfied,KMCIsValid,Found,generate_rand1,StoreLEKMCProcess
      LOGICAL :: GetValidityTime,lexist,MDRequired
      
      IF (PRESENT(WriteTrajectoryFile)) THEN
         OPEN(UNIT=378,FILE=TRIM(WriteTrajectoryFile))
         WRITE(378,*) "#iterations  StateIndex          kTot      ProcessSelected    dt"
      END IF
      
      NotSatisfied=.TRUE.
      KMCIterations=0
      generate_rand1=.TRUE.
      Time0=Time
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      IF (imode==1) THEN
      !Perform nsearches for a single state and report back with list of processes
      !(not more than one state is allowed)
      !if TPMD option is selected that TPMD is performed
         
         Time0=0._dp
         DO WHILE (KMCIterations<NKMCMove)
            StoppingOption=0
            MDStopTime=1.e9_dp
            ObtainRate=0
            NotSatisfied=.TRUE.
            MDRequired=.FALSE.
            DO WHILE (NotSatisfied)
               TotalMDTime=0._dp
               WRITE(6,*) ">>Current state is:",state%Index
               CALL FLUSH(6)
               IF (MDRequired) THEN
                  WRITE(6,*) "MD is required ..."
                  CALL AMDProcessSearch(state%Index,nsearch,Temperature,DetectTransitionFrequency,MDStopTime, &
                     TotalMDTime,nsuccess,StoppingOption,ObtainRate)
                  !update the state information
                  state%NMDTransitions=state%NMDTransitions+nsuccess
                  state%NMDFailed=state%NMDFailed+nsearch-nsuccess
                  state%MDTime=state%MDTime+TotalMDTime !this includes success and failures
         
                  WRITE(6,*) "Number success is ",nsuccess
                  WRITE(6,*) "MD time total is",state%MDTime
         
                  !TPMD steps are perfomed to collect escape times at Temperature
                  CALL AddProcessToAMDGlobalState(state,nsearch,CollectedMDTime,FindRate=.TRUE., &
                     StoreLEKMCProcess=StoreLEKMCProcess,Temperature=Temperature,RemoveFiles=.FALSE.) !CollectedMDTime is total time from successful trajectories
                  WRITE(6,*) "CollectedMDTime,TotalMDTime",CollectedMDTime,TotalMDTime
               END IF

               IF (TPMD) THEN !select escape with shortest time
                  CALL TPMDSelectProcess(state,Temperature,dt,TADParam,MDRequired) !NewBasinTime is the updated time in basin
                  IF (.NOT. MDRequired) THEN
                     IF (dt<0._dp) THEN !this can happen rarely when a new process with a short tick is found
                        WRITE(6,*) "oops ... dt found to be negative"
                        STOP
                     END IF
                     Time0=Time0+dt !time increment
                     Time=Time0
                     CALL CPU_TIME(cputime1)
                     WRITE(6,*) ">>state after move:",state%Index,"@ cpu time",cputime1/60._dp," min"
                  END IF
                  NotSatisfied=MDRequired
               END IF
            END DO
            KMCIterations=KMCIterations+1
            WRITE(6,*) "Current KMC time and move:",Time,KMCIterations
         END DO
      END IF
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      
      IF (imode==2) THEN !Generate validity time for catalog
      
      state%NVisits=state%NVisits+1 !number of visits to the state
      ObtainRate=1 !obtain both activation barrier and prefactor
      
      DO WHILE (NotSatisfied)
         KMCIsValid=state%KMCValidTime>state%KMCTimeExpired
         IF (KMCIsValid) THEN !catalog from current state could be valid
            !IF (generate_rand1) THEN !generate a new value of rand1
               rand1=taus88()
               generate_rand1=.FALSE.
            !END IF
            TotalProcessRate=state%TotalProcessRate
            dt=-log(rand1)/TotalProcessRate
            KMCIsValid=state%KMCValidTime>state%KMCTimeExpired+dt
            WRITE(378,*) ">>dt,tv,tkmcnew,isvalid>>>",dt,state%KMCValidTime, &
                 state%KMCTimeExpired+dt,KMCIsValid
            !WRITE(6,*) ">>>>>>>>>>>>>>>>",dt,state%KMCValidTime,state%KMCTimeExpired+dt,KMCIsValid
            IF (KMCIsValid) THEN !perform a KMC step
               state%KMCTimeExpired=state%KMCTimeExpired+dt
               Time=Time+dt
               CALL WriteGlobalAMDStateInfo(state) !update the validity time expired
               generate_rand1=.TRUE. !we have accepted the rand1 value stored earlier -- so now we need to generate a new one
               rand2=taus88()
               !select a process
               RandProcessRate=TotalProcessRate*rand2
               process=>state%Process
               DO WHILE (ASSOCIATED(process))
                  ProcessRate=process%Rate
                  Found= RandProcessRate<ProcessRate
                  IF (Found) THEN
                     state=>process%FinalState
                     state%NVisits=state%NVisits+1 !number of visits to the state
                     EXIT
                  ELSE
                     RandProcessRate=RandProcessRate-ProcessRate
                  END IF
                  process=>process%NextNeigh
               END DO
               CALL CPU_TIME(cputime)
               IF (PRESENT(WriteTrajectoryFile)) WRITE(378,*) "Accepted the process .. current state is ",state%Index
               CALL FLUSH(378)
               KMCIterations=KMCIterations+1
               NotSatisfied=KMCIterations<NKMCMove
               IF (PRESENT(WriteTrajectoryFile)) &
                  WRITE(UNIT=378,FMT='(2I5,ES15.5,I5,3E15.5)') KMCIterations,state%Index, &
                  TotalProcessRate,process%Index,dt,Time,cputime/60.
            END IF
         END IF
         IF (.NOT. KMCIsValid) THEN !catalog from current state is not valid
            StoppingOption=0
            GetValidityTime=.TRUE.
            MDStopTime=1.e3_dp
            IF (PRESENT(WriteTrajectoryFile)) THEN
               WRITE(378,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
               WRITE(378,*) "MDStopTimePerProcessor:",MDStopTime
               CALL FLUSH(378)
            END IF
            CALL AMDProcessSearch(state%Index,nsearch,Temperature, &
               DetectTransitionFrequency,MDStopTime, &
               TotalMDTime,nsuccess,StoppingOption,ObtainRate)
            !update the state information
            state%NMDTransitions=state%NMDTransitions+nsuccess
            state%NMDFailed=state%NMDFailed+nsearch-nsuccess
            !state%MDTime=state%MDTime+TotalMDTime !this includes success and failures
            
            CALL UpdateSuperBasin(state,nsearch) !add newly found superbasin states
            
            CALL AddProcessToAMDGlobalState(state,nsearch,CollectedMDTime, &
               FindRate=.TRUE.,StoreLEKMCProcess=StoreLEKMCProcess, &
               Temperature=Temperature,RemoveFiles=.TRUE.) !CollectedMDTime is total time from successful trajectories
            state%MDTime=state%MDTime+CollectedMDTime
            CALL GetCatalogValidityTime(state) !obtain validity time
            
            IF (PRESENT(WriteTrajectoryFile)) THEN
               WRITE(378,*) "State index:",state%Index
               WRITE(378,*) "KMC validity:",state%KMCValidTime
               WRITE(378,*) "Successful transitions so far:",state%NMDTransitions
               WRITE(378,*) "Failed transitions so far:",state%NMDFailed
               WRITE(378,*) "Number of processes known so far:",state%NProcessType
               WRITE(378,*) "KMCRate:",state%TotalProcessRate
               !WRITE(378,*) "TargetNBasinEscapes:",TargetNBasinEscapes
               !WRITE(378,*) "NBasinEscapes:",state%NBasinEscapes(1:NTarget)
               !WRITE(378,*) "BasinEscapeTime:",state%BasinEscapeTime(1:NTarget)
               WRITE(378,*) "Total residence time:",state%MDTime
               CALL FLUSH(378)
            END IF
         END IF
         NotSatisfied=KMCIterations<NKMCMove
      END DO
      
      END IF
      
      IF (imode==3) THEN !parallel replica MD
      
      DO WHILE (NotSatisfied)
         StoppingOption=1
         MDStopTime=1.e3_dp
         ObtainRate=0 !dont get the rates
WRITE(6,*) "Begin MD simulations ..."; CALL FLUSH(6)
WRITE(378,*) "Begin MD simulations ..."; CALL FLUSH(378)
         CALL AMDProcessSearch(state%Index,nsearch,Temperature,DetectTransitionFrequency,MDStopTime, &
            TotalMDTime,nsuccess,StoppingOption,ObtainRate)
         charstateindx=INT2CHAR(state%Index)
WRITE(378,*) "Completed MD simulations ...",nsearch
WRITE(378,*) "TotalMDTime",TotalMDTime; CALL FLUSH(378)
         
         dtparrep=1.e8 !it is possible for two escapes to have occurred so we select the shortest one
         
         DO isearch=1,nsearch
            !WRITE(378,*) "isearch:",isearch
            charisearch=INT2CHAR(isearch)
            filename="./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".finalstate.xyz"
         
            !check if file exists
            INQUIRE(FILE=filename,EXIST=lexist)
            IF (.NOT. lexist) THEN
               WRITE(6,*) "Warning>> Unable to find ",TRIM(filename)
               CYCLE
            END IF
            
            !Check if the final state is already recorded
            NULLIFY(AL1)
            CALL ReadCoord(AL1,filename)
            filename="./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".param"
            
            OPEN(UNIT=385,FILE=filename)
            READ(385,*) dt
            CLOSE(385)

!WRITE(378,*) "dt,dtparrep:",dt,dtparrep
            IF (dt<dtparrep) THEN
               dtparrep=dt !reset the value of dtparrep so that finally we select the smallest escape time
               statefinal=> SearchGlobalAMDState(AL1)
               IF (ASSOCIATED(statefinal)) THEN
                  IF (state%Index==statefinal%Index) THEN
WRITE(378,*) "Initial and final states are identical while adding process to AMD state"
WRITE(378,*) "Initial state index:",state%Index
WRITE(378,*) "Final state index:",statefinal%Index
CALL FLUSH(6)
!STOP
                  END IF
               ELSE !could not find final state in the list
                  statefinal=> AddOptimizedALToGlobalAMDStates(AL1)
               END IF
WRITE(378,*) " ... final state is ",statefinal%Index, " after time interval:",dt,TotalMDTime
            END IF
            
            IF (lexist) THEN !delete files
               CALL  SYSTEM("rm "//"./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".param")
               CALL SYSTEM("rm "//"./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".finalstate.xyz")
            END IF
            
         END DO
         
         Time=Time+TotalMDTime
         KMCIterations=KMCIterations+1
         state=>statefinal !this is new state
WRITE(378,*) "Move to new state ..."; CALL FLUSH(378)
WRITE(6,*) "Move to new state ..."; CALL FLUSH(6)
         IF (PRESENT(WriteTrajectoryFile)) THEN
            WRITE(378,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
            WRITE(378,*) "Current state index:",state%Index
            WRITE(378,*) "KMC time:",Time
            WRITE(378,*) "KMC iteration:",KMCIterations
            CALL FLUSH(378)
         END IF
         
         NotSatisfied=KMCIterations<NKMCMove
         IF (PRESENT(StopTime) .AND. NotSatisfied) NotSatisfied= (Time-Time0)<StopTime
         
WRITE(378,*) "reached end of the KMC move"
      END DO
      END IF
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      IF (imode==4) THEN !TAD method
         
         !get the parameters
         TemperatureHigh=TADParam%TemperatureHigh
         IF (TemperatureHigh<Temperature) THEN
            WRITE(6,*) "$Err> TAD temperature is lower than actual system temperature"
            STOP
         END IF
         minprefactor=TADParam%minprefactor
         IF (minprefactor<=0._dp) THEN
            WRITE(6,*) "$Err> Minimum prefactor is zero"
            STOP
         END IF
         delta=TADParam%delta
         IF (delta<=0._dp) THEN
            WRITE(6,*) "$Err> TAD delta is zero"
            STOP
         END IF
         
         ObtainRate=1 !obtain both activation barrier and prefactor
         DO WHILE (NotSatisfied)
            !Step 1. Get t_{L,short} and check if MD searches required
            process=>state%Process
            tLshort=1.e3

!NOTE TIME ACCRUED IS PROBABLY mESSEDUP
            DO WHILE (ASSOCIATED(process))
               IF (process%Ticks(1)+process%TimeAccrued<tLshort) THEN
                  tLshort=process%Ticks(1)+process%TimeAccrued
                  process1=>process !this is the process with the shortest time
               END IF
            END DO
            tHStop=log(1._dp/delta)/minprefactor* &
               (minprefactor*tLshort/log(1._dp/delta))**(TemperatureHigh/Temperature)
            !state%MDTime stores the high temperature time
            MDRequired=tHStop>state%MDTime
            
            !Step 2. if MD search required then perform more TAD
            IF (MDRequired) THEN
               StoppingOption=0
               MDStopTime=1.e3_dp
               IF (PRESENT(WriteTrajectoryFile)) THEN
                  WRITE(378,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
                  WRITE(378,*) "MDStopTimePerProcessor:",MDStopTime
                  CALL FLUSH(378)
               END IF
               CALL AMDProcessSearch(state%Index,nsearch,TemperatureHigh, &
                  DetectTransitionFrequency,MDStopTime, &
                  TotalMDTime,nsuccess,StoppingOption,ObtainRate)
               !update the state information
               state%NMDTransitions=state%NMDTransitions+nsuccess
               state%NMDFailed=state%NMDFailed+nsearch-nsuccess
               !state%MDTime=state%MDTime+TotalMDTime !this includes success and failures
               
               CALL UpdateSuperBasin(state,nsearch) !add newly found superbasin states
            
               CALL AddProcessToAMDGlobalState(state,nsearch,CollectedMDTime, &
                  FindRate=.TRUE.,StoreLEKMCProcess=StoreLEKMCProcess, &
                  Temperature=Temperature,RemoveFiles=.TRUE.) !CollectedMDTime is total time from successful trajectories
               
               state%MDTime=state%MDTime+CollectedMDTime !high temperature time so far
            
               IF (PRESENT(WriteTrajectoryFile)) THEN
                  WRITE(378,*) "State index:",state%Index
                  WRITE(378,*) "Successful transitions so far:",state%NMDTransitions
                  WRITE(378,*) "Failed transitions so far:",state%NMDFailed
                  WRITE(378,*) "Number of processes known so far:",state%NProcessType
                  WRITE(378,*) "KMCRate:",state%TotalProcessRate
                  WRITE(378,*) "Total high temperature residence time:",state%MDTime
                  CALL FLUSH(378)
               END IF
               
            ELSE !Step 3. if stop criterion satisfied then select process and advance
               !process1 is the process with the shortest time
               
               dTime=process1%Ticks(1)* & !time advance for this process
                  EXP(process1%Eact*(1._dp/Temperature))
               Time=Time+process1%Ticks(1)
               process%TimeAccrued=process%TimeAccrued+process1%Ticks(2)
               KMCIterations=KMCIterations+1
               
               process1%NTicks=process1%NTicks-1
               DO i=1,process1%NTicks
                  process1%Ticks(i)=process1%Ticks(i+1)
               END DO
               
               NotSatisfied=KMCIterations<NKMCMove
            END IF
         END DO
      END IF
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      IF (PRESENT(WriteTrajectoryFile)) CLOSE(378)
   END SUBROUTINE ParallelReplicaMDKMC
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!   FUNCTION GetKMCValidityTime(state,method)
!   !finds the time for which the KMC catalog is valid
!   !old method - not in use
!      IMPLICIT NONE
!      TYPE(AMDStateContainer), POINTER :: state
!      INTEGER :: method
!      REAL(dp) :: GetKMCValidityTime,RateEstimateUnknown
!      INTEGER :: NProcessType,NMDTransitions,NMDFailed,WindowSize,PreviousWindowEnd
!      INTEGER :: NKnownProcessType
!      
!      WRITE(6,*) "WARNING >>> FUNCTION GetKMCValidityTime IS NO LONGER IN USE"
!      STOP
!      GetKMCValidityTime=-1._dp
!      SELECT CASE(method)
!      CASE(1) !Bhute-Chatterjee approach (Version 1)
!         IF (state%TotalProcessRate>1.e6_dp) &
!            GetKMCValidityTime=-log(1._dp-state%delta)/ &
!            ((1._dp-state%theta)**(-1./REAL(state%NMDTransitions))-1._dp)/ &
!            state%TotalProcessRate
!      CASE(2) !Bhute-Chatterjee approach (Version 2)
!         NProcessType=state%NProcessType !#processes for state
!         NMDTransitions=state%NMDTransitions !#searches resulting in a transition
!         NMDFailed=state%NMDFailed !#searches that failed in finding a transition
!         WindowSize=state%WindowSize !#transitions that form a window
!         PreviousWindowEnd=state%PreviousWindowEnd !last count where window was completed
!         NKnownProcessType=state%NKnownProcessType  
!         RateEstimateUnknown=state%RateEstimateUnknown
!      !Overview of the method
!      !This method works in terms of windows of successful transitions. 
!      !Initially no processes are known. 
!      !Comment 1. If NMDTransitions<WindowSize and PreviousWindowEnd=0
!      !then MD will search for more processes. 
!      !Comment 2. If NMDTransitions>=WindowSize and PreviousWindowEnd=0 
!      !set PreviousWindowEnd to NMDTransitions 
!      !and NKnownProcessType to NProcessType. 
!      !Comment 3a. no new processes are observed during the next window,
!      !i.e., NProcessType=NKnownProcessType and 
!      !NMDTransitions>=PreviousWindowEnd+WindowSize, 
!      !then we can estimate the rate. 
!      !Comment 3b. Otherwise, a new process has been observed during the next window
!      !i.e., NProcessType>NKnownProcessType
!      !NMDTransitions>=PreviousWindowEnd+WindowSize, 
!      !WindowSize becomes NProcessType
!         !check if there are enough processes that have been observed so far
!         IF (NMDTransitions<WindowSize .AND. PreviousWindowEnd==0) RETURN !comment 1
!         IF (NMDTransitions>=WindowSize .AND. PreviousWindowEnd==0) THEN
!            PreviousWindowEnd=NMDTransitions
!            state%PreviousWindowEnd=NMDTransitions !comment 2
!            NKnownProcessType=NProcessType
!            state%NKnownProcessType=NProcessType
!            RETURN
!         END IF
!         IF (NMDTransitions>=PreviousWindowEnd+WindowSize) THEN !we have completed next window
!            IF (NProcessType==NKnownProcessType) THEN !comment 3a
!               RateEstimateUnknown=10._dp*state%TotalProcessRate/REAL(NMDTransitions,dp)
!               state%RateEstimateUnknown=RateEstimateUnknown
!            ELSEIF (NProcessType>NKnownProcessType) THEN !comment 3b
!               NKnownProcessType=NProcessType
!               state%NKnownProcessType=NProcessType
!               WindowSize=NProcessType
!               state%WindowSize=NProcessType !MAX(NProcessType,WindowSize)
!            ELSE
!               WRITE(6,*) "Err>> NProcessType is found to be less than NKnownProcessType in AMD-KMC"
!               STOP
!            END IF
!            PreviousWindowEnd=NMDTransitions !move into the next window
!            state%PreviousWindowEnd=NMDTransitions
!         END IF
!         IF (RateEstimateUnknown>0._dp) GetKMCValidityTime=-log(1._dp-state%delta)/RateEstimateUnknown
!      END SELECT
!   END FUNCTION GetKMCValidityTime
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetMDStopTime(state,tstopdimensionless)
   !finds the MD stop time for different methods
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state
      REAL(dp) :: ProcessRate,GetMDStopTime,tstopdimensionless
      
      ProcessRate=state%TotalProcessRate
      IF (ProcessRate>0._dp) THEN
         GetMDStopTime=tstopdimensionless/ProcessRate
      ELSE
         GetMDStopTime=1.e3_dp !very long time
      END IF
   END FUNCTION GetMDStopTime
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AMDProcessSearch(stateindex,nsearch,Temperature,DetectTransitionFrequency, &
      MDTimePerBCMD,TotalMDTime,nsuccess,StoppingOption,ObtainRate)
   !search processes from the current state
      IMPLICIT NONE
      INTEGER :: stateindex,nsearch,DetectTransitionFrequency
      REAL(dp) :: Temperature
      REAL(dp) :: MDTimePerBCMD,TotalMDTime1
      REAL(dp), OPTIONAL :: TotalMDTime
      INTEGER :: nsuccess,StoppingOption,ObtainRate
      
      TotalMDTime1=0._dp
      nsuccess=0
      CALL PerformMDProcessSearchParallel(stateindex,nsearch,Temperature,DetectTransitionFrequency, &
         MDTimePerBCMD,TotalMDTime1,nsuccess,StoppingOption,ObtainRate)
      IF (PRESENT(TotalMDTime)) TotalMDTime=TotalMDTime1
   END SUBROUTINE AMDProcessSearch
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !SUBROUTINE ConvertFullKMCToAMDState()
   !   IMPLICIT NONE
      
   !END SUBROUTINE ConvertFullKMCToAMDState
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !FUNCTION CheckLocalProcessAbsent()
   !   IMPLICIT NONE
   !   LOGICAL :: CheckLocalProcessAbsent
      
   !   CheckLocalProcessAbsent=.TRUE.
      
   !END FUNCTION CheckLocalProcessAbsent
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE AMDProcess
