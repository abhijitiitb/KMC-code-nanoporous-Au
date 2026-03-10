MODULE amd_db_manipulate
   USE VARIABLE_TYPE
   USE AMD_VARIABLE_TYPE
   USE utilities
   USE IO
   USE NEB_package
   USE randomgen
   USE ran_state
   USE KMCDbManipulate
   IMPLICIT NONE
   
   REAL(dp), DIMENSION(3*MaxNAtoms), PRIVATE ::  AtomCoordStored,AtomMassStored
   REAL(dp), DIMENSION(MaxNAtoms), PRIVATE ::  AtomChargeStored
   INTEGER, DIMENSION(MaxNAtoms), PRIVATE ::  AtomSpeciesStored
   LOGICAL, DIMENSION(3*MaxNAtoms), PRIVATE ::  AtomIsMovingStored
   
   INTERFACE SearchGlobalAMDState
      MODULE PROCEDURE SearchGlobalAMDStateAL,SearchGlobalAMDStateIndex
   END INTERFACE SearchGlobalAMDState
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION AddOptimizedALToGlobalAMDStates(AL)
   !AL contains the optimized coordinates, global AMD refers to the full system state
   !The state is new, i.e., it has not been visited before, hence the state is added
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(AMDStateContainer), POINTER :: state,AddOptimizedALToGlobalAMDStates
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomCharge
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      INTEGER :: iatom,NAtoms,i1,i3,indx,errorstatus
      CHARACTER(len=100) :: filename,filename2
      CHARACTER(len=NINT2CHAR) :: charindex
      CHARACTER(len=3) :: speciesname
      
      ALLOCATE(state)
      IF (ASSOCIATED(ListOfGlobalAMDStates)) THEN
         state%NextNeigh=>ListOfGlobalAMDStates
         ListOfGlobalAMDStates%PrevNeigh=>state
      ELSE
         CALL SYSTEM("mkdir ./GlobalAMD")
      END IF
      ListOfGlobalAMDStates=>state
      
      NGlobalAMDStates=NGlobalAMDStates+1
      
      WRITE(6,*) ">>>o2>>> Number of AMD states:",NGlobalAMDStates
      
      state%Index=NGlobalAMDStates
      charindex=INT2CHAR(NGlobalAMDStates)
      filename="./GlobalAMD/state."//TRIM(charindex)//".xyz"
      
      filename2="./GlobalAMD/state."//TRIM(charindex)//".tpmdstages"
      CALL SYSTEM("rm "//TRIM(filename2))

      OPEN(UNIT=305,FILE=filename,IOSTAT=errorstatus)
      NAtoms=AL%NAtoms
      AtomCoord=>AL%AtomCoord(1:3*NAtoms)
      AtomSpecies=>AL%AtomSpecies(1:NAtoms)
      AtomCharge=>AL%AtomCharge(1:NAtoms)
      AtomIsMoving=>AL%AtomIsMoving(1:3*NAtoms)
      state%PotentialEnergy=AL%PotentialEnergy
      WRITE(UNIT=305,FMT='(4i8,ES20.8,"  NAtoms Nmove PotentialEnergy")') &
        NAtoms,AL%NMove,AL%PotentialEnergy
      WRITE(UNIT=305,FMT='(6f14.5," StateIndex:",I6)') AL%BoxSize,NGlobalAMDStates
      
      i3=0
      DO iatom=1,NAtoms
         i1=i3+1
         i3=i1+2
         indx=AtomSpecies(iatom)
         IF (indx==0) THEN
            WRITE(6,*) "$Err>> Atom species was found to be zero in AddStateToListOfGlobalStates"
            STOP
         END IF
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(indx))
         WRITE(UNIT=305,FMT='(a3,3f14.5,"  ",3L2,f20.8,I3)') speciesname,AtomCoord(i1:i3), &
            AtomIsMoving(i1:i3),AtomCharge(iatom),indx
         !WRITE(UNIT=305,FMT='(a3,3f14.5,3L,i5,3f8.3)') speciesname,AtomCoord(i1:i3),AtomIsMoving(i1:i3),AtomSpecies(iatom),AtomCharge(iatom)
      END DO
      CLOSE(305)
      
      CALL WriteGlobalAMDStateInfo(state) !create the info file
      AddOptimizedALToGlobalAMDStates=>state
      
   END FUNCTION AddOptimizedALToGlobalAMDStates
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE UpdateSuperBasin(state,nsearch)
   !checks whether files with superbasin state information have been created
   !this state information is present in ./GlobalAMD/state.##.*.newbasinstates.xyz
   !Note: This is currently being used while generating validity of kmc catalog
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: ALtmp,ALtmp1
      TYPE(AMDStateContainer), POINTER :: state,state1
      INTEGER :: nsearch,isearch
      CHARACTER(len=NINT2CHAR) :: charstateindx,charisearch
      CHARACTER(len=100) :: filename
      LOGICAL :: lexist
      
      DO isearch=1,nsearch
      
         !Step 1. Check whether file exists
         charisearch=INT2CHAR(isearch)
         charstateindx=INT2CHAR(state%Index)
         filename="./GlobalAMD/state."//TRIM(charstateindx)//"."//TRIM(charisearch)//".newbasinstates.xyz"
         INQUIRE(FILE=filename,EXIST=lexist)
         IF (.NOT. lexist) CYCLE
      
         !Step 2. Read the entire file into ALtmp -- these could be new superbasin states
         CALL ReadCoord(ALtmp,filename)
         
         !Step 3. Go through each superbasin state
         ALtmp1=>ALtmp
         DO WHILE (ASSOCIATED(ALtmp1))
            !Step 3a. check whether current SB state is new to the SB stored in AL
            state1=>SearchGlobalAMDState(ALtmp1)
            
            !Step 3b. make adjustments to the superbasin datastructure
            IF (ASSOCIATED(state1)) THEN
               IF (state1%Index==state%Index) THEN
                  WRITE(6,*) "$Err>> State index match in UpdateSuperBasin"
                  STOP
               END IF
               CALL MergeSuperbasins(state,state1)
            ELSE
               state1=> AddOptimizedALToGlobalAMDStates(ALtmp1)
               CALL MergeSuperbasins(state,state1)
            END IF
            !Step 3c. Make rearrangements for the states
            ALtmp1=>ALtmp1%NextNeigh
         END DO
         !Step 4. Remove the file
         CALL SYSTEM("rm "//TRIM(filename))
      END DO
   END SUBROUTINE UpdateSuperBasin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MergeSuperbasins(state1,state2)
   !merges the superbasin of the two states
   !also takes care of merging the timeaxis for the superbasins
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state1,state2
      TYPE(SuperbasinContainer), POINTER :: sb1,sb2,sbc,sbn
      TYPE(SuperbasinstateContainer), POINTER :: sbstate
      
      sb1=>state1%superbasin
      sb2=>state2%superbasin
      
      IF (ASSOCIATED(sb1)) THEN
         IF (ASSOCIATED(sb2)) THEN !both sb1 and sb2 are associated
            !choose one of the two superbasins
            IF (sb1%Index/=sb2%Index) THEN !two different superbasins
               WRITE(6,*) "1."; CALL FLUSH(6)
               sbc=>sb1 !superbasin chosen as merged superbasin
               sbn=>sb2 !superbasin not chosen
               IF (sb1%NStates>sb2%NStates) THEN
                  sbc=>sb2 !this will be merged superbasin
                  sbn=>sb1
               END IF
               WRITE(6,*) "2."; CALL FLUSH(6)
               sbstate=>sbn%state
               DO WHILE (ASSOCIATED(sbstate))
                  CALL AddSuperBasinState(sbstate%AMDstate,sbc)
                  sbstate=>sbstate%NextNeigh
               END DO
               WRITE(6,*) "3."; CALL FLUSH(6)
               CALL DeleteSuperbasin(sbn)
            END IF
         ELSE !add state2 to sb1
            CALL AddSuperbasinState(state2,sb1)
            CALL WriteSuperbasinInfo(sb1)
         END IF
      ELSE
         IF (ASSOCIATED(sb2)) THEN !add state1 to sb2
            CALL AddSuperbasinState(state1,sb2)
            CALL WriteSuperbasinInfo(sb2)
         ELSE !both sb1 and sb2 are not associated, create a new superbasin and add these states to the superbasin
            ALLOCATE(sb1)
            CALL AddSuperBasin(sb1)
            CALL AddSuperbasinState(state1,sb1)
            CALL AddSuperbasinState(state2,sb1)
            CALL WriteSuperbasinInfo(sb1)
         END IF
      END IF
   END SUBROUTINE MergeSuperbasins
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddSuperBasin(sb)
   !adds a superbasin to the database
      IMPLICIT NONE
      TYPE(SuperbasinContainer), POINTER :: sb
      
      !add sb to the data structure
      IF (ASSOCIATED(ListSuperbasins)) THEN
         sb%NextNeigh=>ListSuperbasins
         ListSuperbasins%PrevNeigh=>sb
      END IF
      ListSuperbasins=>sb
      
      NAMDSuperbasins=NAMDSuperbasins+1
      RunningIndexSuperbasins=RunningIndexSuperbasins+1
      WRITE(6,*) "Current running index for SB:",RunningIndexSuperbasins
      CALL FLUSH(6)
      sb%Index=RunningIndexSuperbasins
      
   END SUBROUTINE AddSuperBasin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddSuperBasinState(state,sb,IsUpdate)
   !adds state to sb
      IMPLICIT NONE
      TYPE(SuperbasinContainer), POINTER :: sb
      TYPE(AMDStateContainer), POINTER :: state
      TYPE(SuperbasinstateContainer), POINTER :: sbstate
      TYPE(AMDProcessContainer), POINTER :: process
      REAL(dp) :: SB_BasinEscapeTime,BasinEscapeTime
      LOGICAL :: IsUpdate1
      LOGICAL, OPTIONAL :: IsUpdate

      IF (PRESENT(IsUpdate)) THEN
         IsUpdate1=IsUpdate
      ELSE
         IsUpdate1=.TRUE.
      END IF
      
      !in case the state is in correct superbasin already
      !then there is no need to go further
      IF (ASSOCIATED(state%superbasin)) THEN
         IF (state%superbasin%Index==sb%Index) RETURN
      END IF
      
      !get the process times of state on the same time axis as that of sb
      IF (ASSOCIATED(sb%state)) THEN
         IF (ASSOCIATED(sb%state%AMDstate)) THEN
            SB_BasinEscapeTime=sb%state%AMDstate%BasinEscapeTime
            BasinEscapeTime=state%BasinEscapeTime
            IF (BasinEscapeTime/=SB_BasinEscapeTime) THEN
               process=>state%Process
               DO WHILE (ASSOCIATED(process))
                  CALL ShiftTimeAxis(process=process,ref_time_old=BasinEscapeTime,ref_time_new=SB_BasinEscapeTime)
                  CALL PrintAMDProcess(process)
                  process=>process%NextNeigh
               END DO
               state%BasinEscapeTime=SB_BasinEscapeTime
            END IF
         END IF
      END IF
      
      !add the state to the superbasin
      state%superbasin=>sb
      NULLIFY(sbstate)
      ALLOCATE(sbstate)
      sbstate%AMDstate=>state
      sbstate%NextNeigh=>sb%state
      IF (ASSOCIATED(sb%state)) sb%state%PrevNeigh=>sbstate
      sb%state=>sbstate
      sb%NStates=sb%NStates+1
      
      !printing
      IF (IsUpdate1) CALL WriteGlobalAMDStateInfo(state)
   END SUBROUTINE AddSuperBasinState
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ShiftTimeAxis(process,ref_time_old,ref_time_new)
   !shift the time axis so that instead of the reference time ref_time_old
   !that was used earlier, now we use ref_time_new as the reference time for the process
      IMPLICIT NONE
      TYPE(AMDProcessContainer), POINTER :: process
      REAL(dp) :: ref_time_old,ref_time_new,dt
      INTEGER :: nticks
      
      dt=ref_time_new-ref_time_old !shift in the time axis
      nticks=process%NTicks
      process%Ticks(1:nticks)=process%Ticks(1:nticks)+dt
   END SUBROUTINE ShiftTimeAxis
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteSuperbasin(sb)
      IMPLICIT NONE
      TYPE(SuperbasinContainer), POINTER :: sb
      TYPE(SuperbasinstateContainer), POINTER :: sbstate,sbstate1
      CHARACTER(len=100) :: filename
      CHARACTER(len=NINT2CHAR) :: charisuperbasin
      INTEGER :: index
      
      !remove the superbasin states
      sbstate=>sb%state
      NULLIFY(sb%state)
      sb%NStates=0
      DO WHILE (ASSOCIATED(sbstate))
         sbstate1=>sbstate
         sbstate=>sbstate%NextNeigh
         NULLIFY(sbstate1%AMDstate)
         NULLIFY(sbstate1%NextNeigh)
         NULLIFY(sbstate1%PrevNeigh)
         DEALLOCATE(sbstate1)
      END DO
      
      !remove the superbasin
      index=sb%Index
      IF (ASSOCIATED(sb%PrevNeigh)) THEN
         IF (ASSOCIATED(sb%NextNeigh)) THEN
            sb%NextNeigh%PrevNeigh=>sb%PrevNeigh
            sb%PrevNeigh%NextNeigh=>sb%NextNeigh
         ELSE
            NULLIFY(sb%PrevNeigh%NextNeigh)
         END IF
      ELSE
         IF (ASSOCIATED(sb%NextNeigh)) THEN
            NULLIFY(sb%NextNeigh%PrevNeigh)
            ListSuperbasins=>sb%NextNeigh
         ELSE
            NULLIFY(ListSuperbasins)
         END IF
      END IF
      DEALLOCATE(sb)
      NAMDSuperbasins=NAMDSuperbasins-1
      
      !remove the file for superbasin
      charisuperbasin=INT2CHAR(index)
      filename="./GlobalAMD/superbasin."//TRIM(charisuperbasin)//".info"
      CALL SYSTEM("rm "//TRIM(filename))
   END SUBROUTINE DeleteSuperbasin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddProcessToAMDGlobalState(stateinitial,nsearch,CollectedMDTime, &
      FindRate,StoreLEKMCProcess,Temperature,RemoveFiles, &
      AtomCfgIndex,AtomSpecies,NAtoms)
      !adds a process to a state if the process has not been seen before
      !if it has been seen before then the tick mark is updated
      !can add validity time if GetValidityTime is .TRUE.
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: stateinitial,statefinal,TPMDActualInitial
      TYPE(AMDProcessContainer), POINTER :: process
      TYPE(SuperbasinContainer), POINTER :: SB
      TYPE(SuperbasinStateContainer), POINTER :: SBstate
      TYPE(SystemContainer), POINTER :: AL0,AL1,AL2,AL3
      TYPE(ChOSContainer), POINTER :: chos
      REAL(dp) :: SpringConstant,TSEnergy,CollectedMDTime,MDTime,Temperature,TotalProcessRate,EscapeTime
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomCoord1
      REAL(sp) :: distance,distancev(3),maxdistance(MaxNAtoms),maxoveralldisplacement
      INTEGER :: isearch,nsearch,NTicks,i,nimages,ntransitionstates,errorstatus,image,iatom,tmp2
      LOGICAL :: SearchProcess,IsMatch,lexist,FindRate,RemoveFiles,StoreLEKMCProcess
      LOGICAL :: FindRate1,NotSatisfied
      REAL(dp) :: Rate,prefactor,tmp1,PotentialEnergy
      CHARACTER(len=100) :: filename,filename2
      CHARACTER(len=NINT2CHAR) :: charstateindx,charisearch,charprocessindx
      CHARACTER :: CrapChar
      
      INTEGER, OPTIONAL :: NAtoms !number of atoms in the state files - required with LEKMC
      INTEGER, DIMENSION(:), OPTIONAL :: AtomCfgIndex,AtomSpecies !this contains the list of atomic environments -- needed for LEKMC
      TYPE(KMCProcess), POINTER :: prc
      INTEGER :: ProcessMasterAtomEnvIndex,ProcessMasterAtomIndex
      REAL, DIMENSION(3) :: initialcoord,saddlecoord,finalcoord,ReferenceCoord

      REAL(dp), DIMENSION(100000) :: TPMDStageTemperature,TPMDStageTime
      INTEGER, DIMENSION(100000) :: TPMDStage
      INTEGER :: TPMDnstages,TPMDnstages1
      REAL(dp) :: k0,kstage,stateresidenceprobability,TPMDTemperature,PartitionFn,RefEnergy
      INTEGER :: InitialStateIndex,EscapeStateIndex,TPMD_SB_initialstate
      REAL(dp), DIMENSION(:), POINTER :: stateprobability
      REAL(dp), DIMENSION(:,:), POINTER :: TransitionMatrix
      INTEGER, DIMENSION(:), POINTER :: StateMap
      INTEGER :: errorstatus1
      LOGICAL :: EquilibriumSB
      
      CollectedMDTime=0._dp !this is the time only from the successful events
      TPMDnstages=0
      
      charstateindx=INT2CHAR(stateinitial%Index)
      filename2="./GlobalAMD/state."//TRIM(charstateindx)//".tpmdstages"
     
      IF (TPMD) THEN !read .tpmdstage file
         OPEN(UNIT=540,FILE=TRIM(filename2))
         NotSatisfied=.TRUE.
         DO WHILE (NotSatisfied)
            TPMDnstages=TPMDnstages+1
            READ(UNIT=540,FMT=*,IOSTAT=errorstatus) TPMDStageTemperature(TPMDnstages), &
               TPMDStageTime(TPMDnstages),TPMDStage(TPMDnstages)
            IF (errorstatus/=0) THEN
               TPMDnstages=TPMDnstages-1
               NotSatisfied=.FALSE.
            END IF
         END DO
         WRITE(6,*) "Number of TPMD stages from past:",TPMDnstages
         CALL FLUSH(6)
      END IF !note unit 540 needs to be kept open as info will be written
      
      filename="./GlobalAMD/state."//TRIM(charstateindx)//".xyz"
      NULLIFY(AL0)
      CALL ReadCoord(AL0,filename)
      CALL PotentialInitialize(AL0) 
      CALL AddVerletList(AL0,NAtoms=AL0%NAtoms)
      IF (.NOT. ASSOCIATED(stateinitial)) THEN
         WRITE(6,*) "State should have been associated in subroutine AddAMDStateProcess"
         STOP
      END IF
      
      DO isearch=1,nsearch
      
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
         IF (TPMDSuperbasin) THEN !we need potential energy
            OPEN(UNIT=534,FILE=filename) !record potential energy of the system
            READ(UNIT=534,FMT=*,IOSTAT=errorstatus) tmp1,tmp1,tmp1,tmp1,AL1%PotentialEnergy
            IF (errorstatus/=0) THEN
               WRITE(6,*) "Format of xyz file does not provide potential energy"
               STOP
            END IF
            CLOSE(534)
         END IF
         
         !WRITE(6,*) "------------------------------"
         !WRITE(6,*) "Attempt for search:",isearch
         statefinal=> SearchGlobalAMDState(AL1)
         
         SearchProcess=.TRUE.
         IF (ASSOCIATED(statefinal)) THEN
            IF (stateinitial%Index==statefinal%Index) THEN
               WRITE(6,*) "Initial and final states are identical while adding process to AMD state"
               CALL Delete(AL1)
               CYCLE
               !STOP
            END IF
         ELSE !could not find final state in the list
            !WRITE(6,*) "Search failed"
            statefinal=> AddOptimizedALToGlobalAMDStates(AL1)
            SearchProcess=.FALSE. !this is a new state so the process is also new
            !WRITE(6,*) "Number of global states:",NGlobalAMDStates
         END IF
         
         !Search for process in the initial state
         NULLIFY(process)
         IF (SearchProcess) process=>SearchGlobalAMDProcess(stateinitial,statefinal)
         
         !read .param file
         ntransitionstates=0
         filename="./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".param"
         OPEN(UNIT=385,FILE=filename)
         IF (TPMD) THEN
            READ(385,*) EscapeTime,MDTime !escape time will be computed next
         ELSE
            READ(385,*) MDTime
            EscapeTime=MDTime
         END IF
         READ(385,*) PotentialEnergy !we dont need to read temperature
         stateinitial%PotentialEnergy=PotentialEnergy !this is useful for TPMD
         
         READ(385,*) tmp2
         READ(385,*) Rate,prefactor,TSEnergy,ntransitionstates,TPMD_SB_initialstate
         !TPMD_SB_initialstate is the state from whcih the escape occurred
         IF (TPMDSuperbasin) THEN
            TPMDActualInitial=>SearchGlobalAMDStateIndex(TPMD_SB_initialstate)
         ELSE
            TPMDActualInitial=>stateinitial
         END IF
         CALL FLUSH(6)
         
         IF (TPMD .AND. Rate>0._dp) THEN
            !state.??.tpmdstage is already positioned correctly
            READ(385,*) TPMDnstages1 !number of stages
            DO i=0,TPMDnstages1
               READ(385,*) TPMDStageTemperature(TPMDnstages+1+i),TPMDStageTime(TPMDnstages+1+i), &
                 TPMDStage(TPMDnstages+1+i)
               WRITE(540,*) TPMDStageTemperature(TPMDnstages+1+i),TPMDStageTime(TPMDnstages+1+i), &
                 TPMDStage(TPMDnstages+1+i)
            END DO
            TPMDnstages=TPMDnstages+1+TPMDnstages1
            WRITE(6,*) "... number of TPMD stages now becomes:",TPMDnstages
            !obtain the correct escape time by accounting for all stages so far
            EscapeTime=0._dp
            SB=>stateinitial%superbasin
            stateresidenceprobability=1._dp
            
            !IF (TPMDSuperbasin .AND. ASSOCIATED(SB)) THEN
            !   CALL InitializeTransitionMatrix(SB,TransitionMatrix, & !TransitionMatrix is for SB
            !     stateprobability,StateMap,errorstatus1) !stateprobability,StateMap are arrays, StateMap holds the order in which states are stored
            !   CALL PrintTransitionMatrix(TransitionMatrix)
            !   InitialStateIndex=LocateSBstateInStateMap(StateMap,stateinitial%Index)
            !   EscapeStateIndex=LocateSBstateInStateMap(StateMap,TPMD_SB_initialstate) !we need to track probability of this state
            !   stateprobability=0._dp
            !   stateprobability(InitialStateIndex)=1._dp !initially system was in this state
             !  CALL PrintStateProbability(stateprobability)
            !END IF
            
            k0=prefactor*EXP(-TSEnergy/kboltzmann/Temperature)
            IF (TPMDSuperbasin .AND. ASSOCIATED(SB)) THEN
               WRITE(6,*) "Equilibrium SB is being implemented ..."
               WRITE(6,*) "SB index:",SB%index
               WRITE(6,*) "# states in SB:",SB%NStates
               RefEnergy=stateinitial%PotentialEnergy
               SBstate=>SB%state
               PartitionFn=0._dp
               DO WHILE (ASSOCIATED(SBstate))
                  WRITE(6,*) "state, energy:",SBstate%AMDstate%Index, &
                    SBstate%AMDstate%PotentialEnergy
                  PartitionFn=PartitionFn+EXP(-(SBstate%AMDstate%PotentialEnergy -RefEnergy)/ &
                    kboltzmann/Temperature)
                  SBstate=>SBstate%NextNeigh
               END DO
               stateresidenceprobability=1./PartitionFn
               k0=k0*stateresidenceprobability
               WRITE(6,*) "Partition function at system temperature",PartitionFn
               WRITE(6,*) "State probabilty at system temperature:",stateresidenceprobability
            END IF
            
            EquilibriumSB=.TRUE. !assume that equilibrium attained in SB
            
            DO i=1,TPMDnstages
               IF (TPMDSuperbasin .AND. ASSOCIATED(SB)) THEN
                  !Equilibrium approach
             !     IF (EquilibriumSB) THEN
                     TPMDTemperature=TPMDStageTemperature(i)
                     RefEnergy=stateinitial%PotentialEnergy
                     SBstate=>SB%state
                     PartitionFn=0._dp
                     DO WHILE (ASSOCIATED(SBstate))
                        PartitionFn=PartitionFn+EXP(-(SBstate%AMDstate%PotentialEnergy -RefEnergy)/ &
                          kboltzmann/TPMDTemperature)
                        SBstate=>SBstate%NextNeigh
                     END DO
                     stateresidenceprobability=1./PartitionFn

              !    ELSE
                  !Non-equilibrium approach
                  !remember to deallocate the matrices created
                  !presently non-eq option not used so dealloc is disabled at
                  !end of the subroutine
              !       CALL CreateTransitionMatrix(SB,TransitionMatrix,StateMap,TPMDStageTemperature(i))
              !       CALL PrintTransitionMatrix(TransitionMatrix)
              !       IF (TPMDStage(i)==0) THEN
              !          WRITE(6,*) "Using dt:",TPMDStageTime(i)
              !          TransitionMatrix=-TransitionMatrix*TPMDStageTime(i)
              !       ELSE
              !          WRITE(6,*) "Using dt:",TPMDStageTime(i)-TPMDStageTime(i-1)
              !          TransitionMatrix=-TransitionMatrix*(TPMDStageTime(i)-TPMDStageTime(i-1))
              !       END IF
              !       CALL PrintTransitionMatrix(TransitionMatrix)
              !       CALL EXPM(TransitionMatrix)
              !       CALL PrintTransitionMatrix(TransitionMatrix)
              !       stateprobability=MATMUL(TransitionMatrix,stateprobability)
              !       CALL PrintStateProbability(stateprobability)
              !       stateresidenceprobability=stateprobability(EscapeStateIndex)
              !       WRITE(6,*) "TPMD stateresidenceprobability:",stateresidenceprobability
              !    END IF
               END IF
               
               kstage=prefactor*EXP(-TSEnergy/kboltzmann/TPMDStageTemperature(i))*stateresidenceprobability
               
               IF (TPMDStage(i)==0) THEN !EscapeTime is initialized earlier
                  EscapeTime=EscapeTime+kstage*TPMDStageTime(i)/k0
               ELSE
                  EscapeTime=EscapeTime+kstage*(TPMDStageTime(i)-TPMDStageTime(i-1))/k0
               END IF
             !  IF (TPMDSuperbasin .AND. ASSOCIATED(SB)) THEN
             !     CALL CreateTransitionMatrix(SB,TransitionMatrix,StateMap,Temperature)
             !     stateprobability=0._dp
             !     stateprobability(InitialStateIndex)=1._dp
             !     CALL CorrectedEscapeTime(EscapeTime,stateprobability,TransitionMatrix,EscapeStateIndex)
             !  END IF
            END DO
            
            IF (ASSOCIATED(process)) THEN
               WRITE(6,*) "... escape time for process index:",EscapeTime,process%Index
            ELSE
               WRITE(6,*) "... escape time for process index:",EscapeTime," process does not exist"
            END IF
            IF (TPMDSuperbasin .AND. ASSOCIATED(SB)) THEN
               !DEALLOCATE(TransitionMatrix,StateMap,stateprobability)
            END IF
         END IF
         CLOSE(385)
         
         IF (.NOT. ASSOCIATED(process)) THEN !process has not been found
            !IsMatch=CompareALToAL(AL0,AL1,xtol=0.5_dp)
            !IF (.NOT. IsMatch) THEN !initial and final states dont match so we can do NEB
               ntransitionstates=1
               
               FindRate1=FindRate .AND. ntransitionstates==-1 !-1 indicates that prefactor and activation barrier were never computed
               IF (FindRate1 .OR. StoreLEKMCProcess) THEN
                  !NEB calculation
                  TSEnergy=0._dp
                  nimages=11
                  springconstant=0.1_dp
                  NULLIFY(AL2)
                  CALL PotentialInitialize(AL1)
                  !CALL AddVerletList(AL0,NAtoms=AL1%NAtoms)
                  NULLIFY(chos)
                  CALL NEB(AL0,AL1,nimages=nimages,imode=1101,omode=13,interpmode=1, &
                     springconstant=springconstant,ITMAX=100,TSEnergy=TSEnergy, &
                     ntransitionstates=ntransitionstates,TS=AL2,NEBChos=chos,ReuseVL=.FALSE., &
                     errorstatus=errorstatus)
                  IF (errorstatus/=0) THEN !maybe linear interpolation did not work so we use mixed linear-nonlinear interpolation
                     WRITE(6,*) "$Err>> Linear interpolation did not work ..."
                     CALL DeleteChos(chos)
                     CALL Delete(AL2)
                     NULLIFY(chos)
                     NULLIFY(AL2)
                     CALL NEB(AL0,AL1,nimages=nimages,imode=1101,omode=13,interpmode=3, &
                        springconstant=springconstant,ITMAX=100,TSEnergy=TSEnergy, &
                        ntransitionstates=ntransitionstates,TS=AL2,NEBChos=chos,ReuseVL=.FALSE., &
                        errorstatus=errorstatus)
                  END IF
                  prefactor=1.e13_dp
                  Rate=prefactor*EXP(-TSEnergy/kboltzmann/Temperature)
               END IF
               
               IF (ntransitionstates==1 .AND. StoreLEKMCProcess) THEN !create the LEKMC process
                  NAtoms=AL0%NAtoms
                  WRITE(6,*) "LEKMC with MD is being tested right now + there might be memory leaks"
                  STOP
                  
                  !check if all the required information is present
                  IF (.NOT. (PRESENT(AtomCfgIndex) .AND. PRESENT(AtomSpecies) .AND. PRESENT(NAtoms))) THEN
                     WRITE(6,*) "Err>> Additional information required to generated LEKMC process with MD"
                     STOP
                  END IF
                  
                  !now store the relative coordinates
                  AL3=>chos%AL
                  AtomCoord=>AL0%AtomCoord
                  maxdistance=0.
                  maxoveralldisplacement=0.
                  DO image=2,nimages
                     AL3=>AL3%NextNeigh
                     AtomCoord1=>AL3%AtomCoord
                     DO iatom=1,NAtoms
                        distancev=PBCdistance3(AL0,AtomCoord1(3*iatom-2:3*iatom),AtomCoord(3*iatom-2:3*iatom))!PBCdistance3 required updation
                        distance=SQRT(distancev(1)*distancev(1)+distancev(2)*distancev(2)+distancev(3)*distancev(3))
                        maxdistance(iatom)=MAX(maxdistance(iatom),distance) !this is the distance by which the atom moves
                     END DO
                     maxoveralldisplacement=MAX(maxoveralldisplacement,maxdistance(iatom))
                  END DO
                  
                  !find atom to which the process belongs - this is done on the basis of the atom that moves the most and a low environment index
                  ProcessMasterAtomEnvIndex=99999999
                  DO iatom=1,NAtoms
                     IF (AtomCfgIndex(iatom)<0) CYCLE
                     IF (maxdistance(iatom)>0.6*maxoveralldisplacement .AND. &
                        AtomCfgIndex(iatom)<ProcessMasterAtomEnvIndex) THEN
                        
                        ProcessMasterAtomEnvIndex=AtomCfgIndex(iatom)
                        ProcessMasterAtomIndex=iatom
                     END IF
                  END DO
                  
                  !create the LEKMC process
                  ReferenceCoord=AtomCoord(ProcessMasterAtomIndex)
                  ALLOCATE(prc)
                  DO iatom=1,NAtoms
                     distance=maxdistance(iatom)
                     IF (distance>0.1) THEN !it is a moving atom
                        initialcoord=AL0%AtomCoord(3*iatom-2:3*iatom) !relative coord
                        saddlecoord=AL2%AtomCoord(3*iatom-2:3*iatom)
                        finalcoord=AL1%AtomCoord(3*iatom-2:3*iatom)
                        CALL AddAtomToProcess(prc,initialcoord,finalcoord,saddlecoord,AtomCfgIndex(iatom), &
                           AtomSpecies(iatom))
                     END IF
                  END DO
                  
                  !add the rates here
                  prc%FwdBarrier=REAL(TSEnergy)
                  prc%Frequency=1.e13_dp
                  prc%Rate=prc%Frequency*EXP(-TSEnergy/kboltzmann/Temperature)
                  
                  !add process to list of LEKMC processes in state
                  
               END IF
               
               !delete temporary variables that have been created
               IF (FindRate1) THEN
                  CALL DeleteChOS(chos,.TRUE.) !delete entire chos
                  CALL Delete(AL2)
               END IF
                  
               IF (ntransitionstates==1) THEN !NEB makes sense, now create the new AMD process
                  ALLOCATE(process)
                  IF (ASSOCIATED(stateinitial%Process)) THEN
                     stateinitial%Process%PrevNeigh=>process
                     process%NextNeigh=>stateinitial%process
                  END IF
                  stateinitial%Process=>process
                  stateinitial%NProcessType=stateinitial%NProcessType+1
                  process%Index=stateinitial%NProcessType
            
                  
                  process%InitialState=>stateinitial !in TPMD we regard stateinitial as state from which the process was sighted
                  IF (TPMD) THEN
                     IF (TPMDSuperbasin) THEN !the escape could have actually occurred from another state TPMD_SB_initialstate
                        !TPMD_SB_initialstate is integer value for state
                        process%TPMDActualInitial=>TPMDActualInitial !this is the state from which the escape actually occurred
                     ELSE
                        process%TPMDActualInitial=>stateinitial
                     END IF
                  END IF
                  process%FinalState=>statefinal
            
                  process%Eact=TSEnergy  !save the activation energy
                  process%Rate=Rate
                  process%Prefactor=prefactor
                  stateinitial%TotalProcessRate=stateinitial%TotalProcessRate+process%Rate
                  
                  !IF (FindRate) THEN
                  !   charprocessindx=INT2CHAR(process%Index)
                  !   filename="./GlobalAMD/ProcessChOS."//TRIM(charstateindx)//"."//TRIM(charprocessindx)//".xyz"
                  !   CALL SYSTEM("mv -v chos2.xyz "//TRIM(filename))
                  !END IF
                  WRITE(6,*) "Noted information about the process"
               ELSE
                  WRITE(6,*) "Rejecting the processes .. "
                  WRITE(6,*) "... since number of transition states found to be ",ntransitionstates
                  WRITE(6,*) "... and activation energy found to be ",TSEnergy
               END IF
            !END IF
         END IF
         
         IF (ASSOCIATED(process)) THEN !update the tickmarks and print output
            
            NTicks=process%NTicks
            NTicks=NTicks+1
            
            CollectedMDTime=CollectedMDTime+MDTime
            process%MDTimeAccrued=process%MDTimeAccrued+MDTime
            
            IF (NTicks<=MaxNTicks) THEN
               IF (TPMD) THEN !we need absolute time
                  process%Ticks(NTicks)=EscapeTime
                  process%TimeAccrued=EscapeTime !this will be useful while artificially adding ticks
               ELSE
                  process%Ticks(NTicks)=EscapeTime
               END IF
            ELSE
               WRITE(6,FMT='("Warning>> Number of stored ticks has exceeded MaxNTicks for process",I4,"."I4)') &
                  stateinitial%Index,process%Index
            END IF
            
            process%AccruedTicks= process%AccruedTicks+1 !total that should have been there
            NTicks=MIN(NTicks,MaxNTicks)
            process%NTicks=NTicks !actual stored
            
            !write the process information to file
            charprocessindx=INT2CHAR(process%Index)
            filename="./GlobalAMD/Process."//TRIM(charstateindx)//"."//TRIM(charprocessindx)
         
            OPEN(UNIT=385,FILE=filename)
            CALL PrintAMDProcess(process)
            CLOSE(385)
            
         END IF
         
         CALL Delete(AL1) !delete final state
      
      END DO
      
      CALL WriteGlobalAMDStateInfo(stateinitial)
      
      TotalProcessRate=stateinitial%TotalProcessRate

      IF (TPMD) CLOSE(540)
      
      !Validity time based on m/tB calculation
      DO isearch=1,nsearch
         charisearch=INT2CHAR(isearch)
      !   IF (GetValidityTime) THEN
      !      filename="./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".param"
      !      INQUIRE(FILE=filename,EXIST=lexist)
      !      IF (.NOT. lexist) CYCLE
      !      OPEN(UNIT=385,FILE=filename)
      !      READ(385,*) MDTime
      !      CLOSE(385)
      !      DO i=1,NTarget
      !         IF (stateinitial%NBasinEscapes(i)==TargetNBasinEscapes(i)) CYCLE
      !         !Note: TargetStopTime is dimensionless quantity
      !         IF (MDTime<TargetStopTime(i)/TotalProcessRate) THEN
      !            stateinitial%NBasinEscapes(i)=stateinitial%NBasinEscapes(i)+1
      !            stateinitial%BasinEscapeTime(i)=stateinitial%BasinEscapeTime(i)+MDTime
      !         END IF
      !         IF (stateinitial%NBasinEscapes(i)==TargetNBasinEscapes(i)) THEN
      !            WRITE(6,*) "Calculating the validity time ..."
      !            WRITE(6,*) "Current validity time:",stateinitial%KMCValidTime
      !            WRITE(6,*) "m:",stateinitial%NBasinEscapes(i)
      !            WRITE(6,*) "tB:",stateinitial%BasinEscapeTime(i)
      !            WRITE(6,*) "kK:",TotalProcessRate
      !            stateinitial%KMCValidTime=MAX(stateinitial%KMCValidTime, &
      !            -LOG(1-stateinitial%delta)/ &
      !            (REAL(stateinitial%NBasinEscapes(i),dp)/stateinitial%BasinEscapeTime(i)- &
      !            TotalProcessRate))
      !            WRITE(6,*) "New validity time:",stateinitial%KMCValidTime
      !            CALL FLUSH(6)
      !         END IF
      !      END DO
      !   END IF
         IF (RemoveFiles) THEN
            CALL SYSTEM("rm "//"./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".param")
            CALL SYSTEM("rm "//"./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".finalstate.xyz")
          END IF
      END DO
      CALL Delete(AL0) !delete the initial state
      
   END SUBROUTINE AddProcessToAMDGlobalState
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetCatalogValidityTime(state)
   !finds the validity time for the catalog belonging to a state
      IMPLICIT NONE
      INTEGER, PARAMETER :: MaxBand=1000 !maximum number of bands
      TYPE(AMDStateContainer), POINTER :: state
      INTEGER :: nband,iband
      REAL(dp) :: band(2*MaxBand),ku,spectral_band_width
      !The array band contains the rate information for each band - odd element is the largest and even is smallest rate in the band 
      
      spectral_band_width=5._dp
      CALL AssignBands(state,nband,band,MaxBand,spectral_band_width) !create bands
      WRITE(378,*) Nband
      WRITE(378,*) band(1:2*Nband)
      IF (nband==0) RETURN
      
      ku=0._dp
      DO iband=1,nband
         ku=ku+kuBand(state,iband,band,MaxBand) !from accessible bands
      END DO
      ku=ku+10._dp/state%MDTime
      WRITE(378,*) "Total ku:",ku
      
      IF (state%KMCValidTime<0._dp .OR. ku<state%RateEstimateUnknown) &
         state%RateEstimateUnknown=ku
      
      IF (state%RateEstimateUnknown>0._dp) &
         state%KMCValidTime=-LOG(1._dp-state%delta)/state%RateEstimateUnknown
      WRITE(378,*) "Validity time:",state%KMCValidTime
      
      CALL WriteGlobalAMDStateInfo(state)
   END SUBROUTINE GetCatalogValidityTime
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AssignBands(state,nband,band,MaxBand,spectral_band_width)
      IMPLICIT NONE
      INTEGER :: nband,MaxBand
      INTEGER :: NProcessesConsidered,NProcessType
      TYPE(AMDStateContainer), POINTER :: state
      TYPE(AMDProcessContainer), POINTER :: process
      REAL(dp), DIMENSION(2*MaxBand) :: band
      REAL(dp) :: spectral_band_width,kmax,kmin,ku
      LOGICAL :: NotSatisfied
      
      band=0._dp
      nband=0
      
      NProcessType=state%NProcessType
      IF (NProcessType==0) RETURN
      
      !create the first band
      kmax=0._dp
      process=>state%Process
      DO WHILE (ASSOCIATED(process))
         kmax=MAX(kmax,process%Rate)
         write(378,*) "Process rate :",process%Rate
         process=>process%NextNeigh
      END DO
      kmin=kmax/spectral_band_width
      nband=1
      band(1)=kmax
      band(2)=kmin
      
      !create the remaining bands
      NotSatisfied=.TRUE.
      DO WHILE (NotSatisfied)
         kmax=0._dp
         process=>state%Process
         DO WHILE (ASSOCIATED(process))
            IF (process%Rate<kmin) kmax=MAX(kmax,process%Rate)
            process=>process%NextNeigh
         END DO
         IF (kmax==0._dp) THEN
            NotSatisfied=.FALSE.
         ELSE
            kmin=kmax/spectral_band_width
            nband=nband+1
            band(2*nband-1)=kmax
            band(2*nband)=kmin
         END IF
      END DO
   END SUBROUTINE AssignBands
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION kuBand(state,iband,band,MaxBand)
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state
      TYPE(AMDProcessContainer), POINTER :: prc
      REAL(dp) :: kuBand,kb
      REAL(dp), SAVE :: kb1
      INTEGER :: iband,MaxBand
      INTEGER :: Np,Nk,Nt
      INTEGER, SAVE :: Np1,Nt1,Nk1
      REAL(dp), DIMENSION(2*MaxBand) :: band
      
      kuBand=0._dp
      
      !get the number of processes from the band
      !Np1 is estimated number of processes from band 1
      !kb1 is the average rate in band 1
      !Nt1 is the number of transitions from 1
      WRITE(378,*) "##############Band #",iband
      IF (iband==1) THEN
         CALL MLEBinomialDistribution(state,iband,band,MaxBand,kb1,Np1,Nk1,Nt1)
         kb=kb1
         Np=Np1
         Nk=Nk1 !number of known processes
      ELSE
         !CALL MLEUniformDistribution(state,iband,band,MaxBand,kb,Np,Nk,kb1,Np1,Nt1) !kb1,Np1 and Nt1 are inputs
         CALL MLEBinomialDistribution(state,iband,band,MaxBand,kb,Np,Nk,Nt)
      END IF
      
      kuBand=kb*REAL(MAX(Np,Nk)-Nk,dp) !if a band cannot contribute then kuBand is zero
      WRITE(378,*) "Found Np,Nk,kb:",Np,Nk,kb
      CALL FLUSH(378)
      !kuBand is zero when
      ! - not enough datapoints to perform MLEBinomialDistribution
      ! - tb>10/kb for the band
      ! - if one of the later bands is using band 1 and not no estimate for band 1 is available, for e.g., the first band has only one process
      
   END FUNCTION kuBand
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MLEBinomialDistribution(state,iband,band,MaxBand,kb,Np,Nk,Nt)
   !kb- average rate for the band
   !Np- estimated number of processes from the band
   !Nk- number of known processes
   !Nt- number of transitions for the band
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state
      TYPE(AMDProcessContainer), POINTER :: process
      INTEGER :: iband,MaxBand,Np,Nk,Nt,idata,ndatapoints
      INTEGER :: Count(10),n
      REAL(dp), DIMENSION(2*MaxBand) :: band
      LOGICAL :: IsFitPossible,NotConverged
      REAL(dp) :: errleft,errcenter,errright
      REAL(dp) :: kb,tb,kmax,kmin
      
      kmax=band(2*iband-1) !maximum rate in the selected band
      kmin=band(2*iband) !minimum rate in the selected band
      Count=0 !Count(n) store number of times a process has been observed n times
      
      !create the histogram and obtain the average rate
      kb=0._dp !this is number averaged rate
      Nk=0 !this is the number of known processes in the band
      Nt=0 !number of transitions observed for the band
      Np=0
      process=>state%Process
      DO WHILE (ASSOCIATED(process))
         IF (process%Rate<=kmax .AND. process%Rate>kmin) THEN
            n=process%NTicks
            IF (n<=10) Count(n)=Count(n)+1 !IMPORTANT MaxNTicks HAS TO BE greater than 10
            WRITE(378,*) "Number of ticks for process in BD fit:",n
            kb=kb+process%Rate*REAL(n,dp)
            Nk=Nk+1 !number known process
            nt=nt+n
         END IF
         process=>process%NextNeigh
      END DO
      kb=kb/REAL(nt)
      WRITE(378,*) "kb:",kb," nt:",nt
      
      !if the residence time is greater than 10./kb return
      !state%MDTime<=10._dp/kwtav implies that definitely we should have scanned through the band
      tb=state%MDTime
      WRITE(378,*) "tb, 10/kb",tb,10._dp/kb
      IF (tb>10._dp/kb) THEN
         Np=Nk
         RETURN
      END IF
      
      !check if we can perform fitting
      ndatapoints=0
      DO idata=1,10
         IF (Count(idata)>0) ndatapoints=ndatapoints+1
      END DO
      
      CALL FLUSH(378)
      
      IsFitPossible=ndatapoints>1
      !ndatapoints>3 implies we have enough data points to pin down Np
      
      IF (.NOT. IsFitPossible) RETURN
      
      !we have enough data points to do a fit
      Np=Nk !total number of processes predicted from fitting
      NotConverged=.TRUE.
      OPEN(UNIT=344,FILE="CHECK-ERRORCONVERGENCE")
      WRITE(344,*) Np-1,Np,Np+1; CALL FLUSH(344)

      errleft=NpSquaredError(Count,Np-1,Nt)
      errcenter=NpSquaredError(Count,Np,Nt)
      errright=NpSquaredError(Count,Np+1,Nt)
      
      DO WHILE (NotConverged)
         NotConverged=.NOT. (errcenter<errleft .AND. errcenter<errright)
         
         IF (Np==0) EXIT
         WRITE(344,*) Np; CALL FLUSH(344)
         IF (errcenter>errleft) THEN !go towards the left
            Np=Np-1
            WRITE(344,*) Np; CALL FLUSH(344)
            errright=errcenter
            errcenter=errleft
            errleft=NpSquaredError(Count,Np-1,Nt)
         ELSEIF (errcenter>errright) THEN !go towards the right
            Np=Np+1
            WRITE(344,*) Np; CALL FLUSH(344)
            errleft=errcenter
            errcenter=errright
            errright=NpSquaredError(Count,Np+1,Nt)
         END IF
      END DO
      WRITE(344,*) "Completed"
      CLOSE(344)
   END SUBROUTINE MLEBinomialDistribution
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MLEUniformDistribution(state,iband,band,MaxBand,kb,Np,Nk,kb1,Nk1,Nt1)
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state
      TYPE(AMDProcessContainer), POINTER :: process
      INTEGER :: iband,MaxBand,Np,NK,Nt1,Nk1,Ntb,n
      REAL(dp), DIMENSION(2*MaxBand) :: band
      REAL(dp) :: kb,kb1,kmax,kmin,tb
      
      kmax=band(2*iband-1) !maximum rate in the selected band
      kmin=band(2*iband) !minimum rate in the selected band
      Ntb=0 !number of tranistions from b
      kb=0 !average rate from b
      Np=0 !estimated number of process from b
      Nk=0 !known processes from b
      process=>state%Process
      DO WHILE (ASSOCIATED(process))
         IF (process%Rate<=kmax .AND. process%Rate>kmin) THEN
            n=process%NTicks
            WRITE(378,*) "Number of ticks for process in UD fit:",n
            kb=kb+process%Rate*REAL(n,dp)
            Nk=Nk+1 !number known process
            Ntb=Ntb+n
         END IF
         process=>process%NextNeigh
      END DO
      kb=kb/REAL(Ntb)
      
      !if the residence time is greater than 10./kb return
      !state%MDTime<=10._dp/kwtav implies that definitely we should have scanned through the band
      tb=state%MDTime
      WRITE(378,*) "tb, 10/kb",tb,10._dp/kb
      IF (tb>10._dp/kb) THEN
         Np=Nk
         RETURN
      END IF
      
      WRITE(378,*) "kb:",kb," nt:",ntb
      Np=Nk1*(kb1/kb)*Ntb/Nt1
   END SUBROUTINE MLEUniformDistribution
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION NpSquaredError(Count,Np,m)
   !m is the number of transitions from the band
   !Np is number of processes from the band
      IMPLICIT NONE
      INTEGER :: Np,m,Count(10),i
      REAL(dp) :: NpSquaredError,err
      
      NpSquaredError=0._dp
      DO i=1,10
         err= Count(i) -  Np*Combination(m,i)* &
            (1._dp/REAL(Np,dp))**i  *   (1._dp-1._dp/REAL(Np,dp))**(m-i)
         NpSquaredError=NpSquaredError+err*err
      END DO
      
   END FUNCTION NpSquaredError
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION Combination(m,i)
      IMPLICIT NONE
      INTEGER :: i,j,k,m
      REAL(dp) :: Combination
      
      Combination=1._dp
      k=MIN(i,m-i)
      DO j=1,k
         Combination=Combination*REAL(m-k+j)/REAL(j)
      END DO
   END FUNCTION Combination
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION AddAMDGlobalProcessWithDetails(finalstate,prefactor,activationbarrier,Temperature)
   !creates a process that stores the relavent process information
      IMPLICIT NONE
      TYPE(AMDProcessContainer), POINTER :: AddAMDGlobalProcessWithDetails
      TYPE(AMDStateContainer), POINTER :: finalstate
      REAL(dp) :: activationbarrier,Temperature,rate,kBT,prefactor
      
      NULLIFY(AddAMDGlobalProcessWithDetails)
      ALLOCATE(AddAMDGlobalProcessWithDetails)
      AddAMDGlobalProcessWithDetails%FinalState=>finalstate
      AddAMDGlobalProcessWithDetails%prefactor=prefactor
      AddAMDGlobalProcessWithDetails%Eact=activationbarrier
      
      kBT=1.3806503e-23_dp*Temperature/1.60217646e-19
      AddAMDGlobalProcessWithDetails%Rate=prefactor*EXP(-activationbarrier/kBT)
   END FUNCTION AddAMDGlobalProcessWithDetails
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddAMDGlobalProcessTick(process,tick,exceededarraysize)
   !adds a tick to the process "process"
      IMPLICIT NONE
      TYPE(AMDProcessContainer), POINTER :: process
      REAL(dp) :: tick
      LOGICAL :: exceededarraysize
      INTEGER :: NTicks
      
      NTicks=process%NTicks+1
      IF (NTicks<MaxNTicks) THEN
         process%Ticks(NTicks)=tick
         exceededarraysize=.FALSE.
      ELSE
         exceededarraysize=.TRUE.
      END IF
   END SUBROUTINE AddAMDGlobalProcessTick
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchGlobalAMDProcess(stateinitial,statefinal)
   !finds the process from the stateinitial to statefinal
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: stateinitial,statefinal 
      TYPE(AMDProcessContainer), POINTER :: SearchGlobalAMDProcess
      
      SearchGlobalAMDProcess=>stateinitial%Process
      DO WHILE (ASSOCIATED(SearchGlobalAMDProcess))
         IF (SearchGlobalAMDProcess%finalstate%Index==statefinal%Index) RETURN
         SearchGlobalAMDProcess=>SearchGlobalAMDProcess%NextNeigh
    END DO
   END FUNCTION SearchGlobalAMDProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchAMDSuperbasinIndex(SBIndex,IsForceCreate)
   !searches for the superbasin with index SBIndex
      TYPE(SuperbasinContainer), POINTER :: SB,SearchAMDSuperbasinIndex
      INTEGER :: SBIndex
      LOGICAL, OPTIONAL :: IsForceCreate

      SearchAMDSuperbasinIndex=>NULL()
      IF (SBIndex<0) THEN
         WRITE(6,*) "Err>> Negative SBindex passed to SearchAMDSuperbasinIndex"
         STOP
      END IF

      SB=>ListSuperbasins
      DO WHILE (ASSOCIATED(SB))
         IF (SB%index==SBIndex) THEN
            SearchAMDSuperbasinIndex=>SB
            RETURN
         END IF
         SB=>SB%NextNeigh
      END DO

      !SB was not found
      IF (PRESENT(IsForceCreate)) THEN
         IF (IsForceCreate) THEN
            NULLIFY(SB)
            ALLOCATE(SB)
            CALL AddSuperBasin(SB)
            SB%Index=SBIndex
            RunningIndexSuperbasins=MAX(SBIndex,RunningIndexSuperbasins) !useful for assigning SB index for future SB
            SearchAMDSuperbasinIndex=>SB
         END IF
      END IF
   END FUNCTION SearchAMDSuperbasinIndex
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchGlobalAMDStateAL(AL)
   !finds the AMD state that matches with AL
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(AMDStateContainer), POINTER :: SearchGlobalAMDStateAL
      
      SearchGlobalAMDStateAL=>ListOfGlobalAMDStates
      DO WHILE (ASSOCIATED(SearchGlobalAMDStateAL))
         !WRITE(6,*) "Comparing to state ",SearchGlobalAMDState%Index
         IF (CompareGlobalAMDState(SearchGlobalAMDStateAL,AL,tolerance=1._dp)) RETURN
         SearchGlobalAMDStateAL=>SearchGlobalAMDStateAL%NextNeigh
      END DO
   END FUNCTION SearchGlobalAMDStateAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchGlobalAMDStateIndex(index)
   !finds the AMD state with matching index
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: SearchGlobalAMDStateIndex
      INTEGER :: Index
      
      SearchGlobalAMDStateIndex=>ListOfGlobalAMDStates
      DO WHILE (ASSOCIATED(SearchGlobalAMDStateIndex))
         IF (SearchGlobalAMDStateIndex%Index==Index) RETURN
         SearchGlobalAMDStateIndex=>SearchGlobalAMDStateIndex%NextNeigh
      END DO
   END FUNCTION SearchGlobalAMDStateIndex
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CompareGlobalAMDState(state,AL,tolerance)
   !true if the state matches with atom list AL
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(AMDStateContainer), POINTER :: state
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomCharge
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      REAL(dp) :: BoxSize(3),tolerance,coord(3),coord1(3),Charge
      LOGICAL, SAVE :: NotInitialized=.TRUE.
      LOGICAL :: CompareGlobalAMDState,Found,IsMoving(3)
      CHARACTER(len=100) :: filename
      CHARACTER(len=NINT2CHAR) :: charindex
      CHARACTER(len=3) :: CrapChar,species1,species2
      CHARACTER :: CrapChar1
      INTEGER :: i1,i3,indx,iatom,NAtoms,NMove(3)
      
      !read in coordinates from files
      indx=state%Index
      charindex=INT2CHAR(indx)
      !filename="./GlobalAMD/"//TRIM(charindex)
      filename="./GlobalAMD/state."//TRIM(charindex)//".xyz"
      
      AtomCoord=>AL%AtomCoord
      AtomSpecies=>AL%AtomSpecies
      AtomCharge=>AL%AtomCharge
      AtomIsMoving=>AL%AtomIsMoving
      BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
      
      CompareGlobalAMDState=.FALSE.
      OPEN(UNIT=305,FILE=filename)
      READ(305,*) NAtoms,NMove
      IF (AL%NAtoms/=NAtoms .OR. ANY(AL%NMove/=NMove)) THEN
         !WRITE(6,*) "Number of atoms does not match"
         CLOSE(305)
         RETURN
      END IF
      READ(305,*) CrapChar1
      i3=0
      CompareGlobalAMDState=.FALSE.
      DO iatom=1,NAtoms
         i1=i3+1
         i3=i1+2
         species1=SpeciesList%AtomicSymbol(SpeciesDirectory_global(AtomSpecies(iatom)))
         READ(305,*) species2,Coord,IsMoving !,Charge
         coord1=AtomCoord(i1:i3)-coord
         coord1=coord1-BoxSize*NINT(coord1/BoxSize) !relative distance
         Found= MAXVAL(ABS(coord1))<tolerance .AND. species1==species2 .AND. &
            ALL(AtomIsMoving(i1:i3) .EQV. IsMoving) !.AND. Charge==AtomCharge(iatom)
         IF (.NOT. Found) THEN
            !WRITE(6,*) "Atom coordinates are different"
            CLOSE(305)
            RETURN
         END IF
      END DO
      CLOSE(305)
      !WRITE(6,*) "Match found"
      
      CompareGlobalAMDState=.TRUE.
   END FUNCTION CompareGlobalAMDState
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteGlobalAMDState(state,deletefile)
   !deletes state from the global AMD state list
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state
      LOGICAL :: deletefile
      CHARACTER(len=100) :: filename
      CHARACTER(len=NINT2CHAR) :: charindex
      INTEGER :: indx
      
      IF (ASSOCIATED(state%NextNeigh)) THEN
         IF (ASSOCIATED(state%PrevNeigh)) THEN
            state%NextNeigh%PrevNeigh=>state%PrevNeigh
            state%PrevNeigh%NextNeigh=>state%NextNeigh
         ELSE
            ListOfGlobalAMDStates=>state%NextNeigh
            NULLIFY(ListOfGlobalAMDStates%PrevNeigh)
         END IF
      ELSE
         IF (ASSOCIATED(state%PrevNeigh)) THEN
            NULLIFY(state%PrevNeigh%NextNeigh)
         ELSE
            NULLIFY(ListOfGlobalAMDStates)
         END IF
      END IF
      NULLIFY(state%PrevNeigh)
      NULLIFY(state%NextNeigh)
      CALL DeleteAMDGlobalStateAllProcesses(state)
      NGlobalAMDStates=NGlobalAMDStates-1
      
      IF (deletefile) THEN
         indx=state%Index
         charindex=INT2CHAR(indx)
         !filename="./GlobalAMD/"//TRIM(charindex)
         filename="./GlobalAMD/state."//TRIM(charindex)//".xyz"
         CALL SYSTEM("rm -v "//TRIM(filename))
      END IF
      DEALLOCATE(state)
   END SUBROUTINE DeleteGlobalAMDState
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteAMDGlobalStateAllProcesses(state)
   !deletes all process from the state "state"
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state
      TYPE(AMDProcessContainer), POINTER :: process,process1
      
      process=>state%Process
      DO WHILE (ASSOCIATED(process))
         process1=>process%NextNeigh
         CALL DeleteAMDGlobalStateProcess(state,process)
         process=>process1
      END DO
   END SUBROUTINE DeleteAMDGlobalStateAllProcesses
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteAMDGlobalStateProcess(state,process)
   !Removes a specific process from the list of processes of a state
      IMPLICIT NONE
      TYPE(AMDProcessContainer), POINTER :: process
      TYPE(AMDStateContainer), POINTER :: state
      
      IF (.NOT. ASSOCIATED(process)) THEN
         WRITE(6,*) "Err>> Process is not associated in DeleteAMDGlobalStateProcess"
         STOP
      END IF
      
      IF (ASSOCIATED(process%PrevNeigh)) THEN
         IF (ASSOCIATED(process%NextNeigh)) THEN
            process%PrevNeigh%NextNeigh=>process%NextNeigh
            process%NextNeigh%PrevNeigh=>process%PrevNeigh
         ELSE
            NULLIFY(process%PrevNeigh%NextNeigh)
         END IF
      ELSE
         IF (ASSOCIATED(process%NextNeigh)) THEN
            state%Process=>process%NextNeigh
            NULLIFY(state%Process%PrevNeigh)
         ELSE
            NULLIFY(state%Process)
         END IF
      END IF
      NULLIFY(process%PrevNeigh)
      NULLIFY(process%NextNeigh)
      
      state%NProcessType=state%NProcessType-1
      IF (state%NProcessType<0) THEN
         WRITE(6,*) "Err>> Number processes belonging to state less than zero in DeleteAMDGlobalStateProcess"
         STOP
      END IF
      
      NULLIFY(process%InitialState)
      NULLIFY(process%FinalState)
      DEALLOCATE(process)
      
   END SUBROUTINE DeleteAMDGlobalStateProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ConvertToAL(NAtoms,Temperature,iseed)
   !converts array data into AL information, sets up the velocities, linked list and verlet list
   !used by slave processors at the beginning of parallel MD runs
      IMPLICIT NONE
      LOGICAL :: NotInitialized=.TRUE.
      TYPE(SystemContainer), POINTER :: ALSbrtn,ConvertToAL !AL is the pointer to the system container ALSbrtn
      INTEGER :: NAtoms,iseed
      REAL(dp) :: Temperature
      REAL(dp), DIMENSION(NAtoms) :: AtomCharge
      REAL(dp), DIMENSION(3*NAtoms) :: AtomCoord,AtomMass
      INTEGER, DIMENSION(NAtoms) :: AtomSpecies
      
      IF (NotInitialized) THEN
         NotInitialized=.FALSE.
         CALL MakeSize(ALSbrtn,NAtoms)
         !read potential information
      END IF
      ALSbrtn%AtomCoord(1:3*NAtoms)=AtomCoordStored(1:3*NAtoms)
      ALSbrtn%AtomMass(1:3*NAtoms)=AtomMassStored(1:3*NAtoms)
      ALSbrtn%AtomIsMoving(1:3*NAtoms)=AtomIsMovingStored(1:3*NAtoms)
      ALSbrtn%AtomSpecies(1:NAtoms)=AtomSpeciesStored(1:NAtoms)
      ALSbrtn%AtomCharge(1:NAtoms)=AtomChargeStored(1:NAtoms)
      ConvertToAL=>ALSbrtn
   END FUNCTION ConvertToAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RecoverStateInformation(IsPrint,NoTicks)
   !used while restarting any AMD-KMC run
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state
      TYPE(AMDProcessContainer), POINTER :: process
      TYPE(SuperbasinContainer), POINTER :: SB
      LOGICAL :: NotSatisfied,lexist,lexist1,IsPrint1,NoTicks1
      INTEGER :: istate,nprocess,iprocess,NTicks
      INTEGER :: stateinitial_index,statefinal_index,SBIndex,apparent_initial_index
      CHARACTER(len=NINT2CHAR) :: charstateindx,charprocessindx
      CHARACTER(len=100) :: filename
      LOGICAL, OPTIONAL :: IsPrint,NoTicks
      
      IF (PRESENT(IsPrint)) THEN
         IsPrint1=IsPrint
      ELSE
         IsPrint1=.FALSE.
      END IF
      IF (PRESENT(NoTicks)) THEN
         NoTicks1=NoTicks
      ELSE
         NoTicks1=.FALSE.
      END IF

      !read all state files -- state.xx.xyz
      istate=0
      !iprocess=0
      NULLIFY(ListOfGlobalAMDStates)
      NotSatisfied=.TRUE.
      DO WHILE (NotSatisfied)
         istate=istate+1
         charstateindx=INT2CHAR(istate)
         filename="./GlobalAMD/state."//TRIM(charstateindx)//".xyz"
         INQUIRE(FILE=filename,EXIST=lexist1)
         
         IF (.NOT. lexist1) EXIT !did not find the file -- end of list of states
         
         ALLOCATE(state)
         NGlobalAMDStates=NGlobalAMDStates+1
         state%index=istate
         IF (ASSOCIATED(ListOfGlobalAMDStates)) THEN
            ListOfGlobalAMDStates%PrevNeigh=>state
            state%NextNeigh=>ListOfGlobalAMDStates
            ListOfGlobalAMDStates=>state
         ELSE
            ListOfGlobalAMDStates=>state
         END IF
         
         filename="./GlobalAMD/state."//TRIM(charstateindx)//".info"
         INQUIRE(FILE=filename,EXIST=lexist)
         IF (lexist) THEN
            !read info file
            IF (IsPrint1) WRITE(6,*) "Reading ",TRIM(filename)
            OPEN(UNIT=385,FILE=filename)
            READ(385,*) state%index,state%NVisits,SBIndex,state%PotentialEnergy
            IF (SBIndex>0) THEN !for TPMD method
               SB=>SearchAMDSuperbasinIndex(SBindex,.TRUE.)
               CALL AddSuperBasinState(state,SB,.FALSE.)
            END IF
            READ(385,*) state%NProcessType
            READ(385,*) state%NMDTransitions
            READ(385,*) state%NMDFailed
            READ(385,*) state%TotalProcessRate
            READ(385,*) state%delta
            IF (TPMD) THEN
               READ(385,*) state%BasinEscapeTime !time elapsed in basin
            ELSE !read validity time
               READ(385,*) state%KMCValidTime
            END IF
            READ(385,*) state%KMCTimeExpired
            READ(385,*) state%RateEstimateUnknown
            READ(385,*) state%MDTime
         END IF
         CLOSE(385)
         
         WRITE(6,*) "Added state #",istate
         
      END DO
      
      !read all process files -- process.yy.xx
      DO istate=1,NGlobalAMDStates
         state=>SearchGlobalAMDStateIndex(istate)
         IF (state%Index/= istate) THEN
            WRITE(6,*) "Err>> State index does not match with istate in RecoverStateInformation"
            STOP
         END IF
         charstateindx=INT2CHAR(istate)
         nprocess=state%NProcessType
         !WRITE(6,*) "Number of Processes:",nprocess
         DO iprocess=1,nprocess
            ALLOCATE(process)
          !  WRITE(6,*) "Process Index:",iprocess
            charprocessindx=INT2CHAR(iprocess)
            filename="./GlobalAMD/Process."//TRIM(charstateindx)//"."//TRIM(charprocessindx)
            INQUIRE(FILE=filename,EXIST=lexist)
            IF (.NOT. lexist) THEN
               WRITE(6,*) "Err>> Process file not found in RecoverStateInformation:",TRIM(filename)
              ! WRITE(6,*) "State index:",istate
              ! WRITE(6,*) "Number of processes:",nprocess
               STOP
            END IF
            WRITE(6,*) ".... accessing ..",TRIM(filename)
            OPEN(UNIT=385,FILE=filename)
            IF (TPMD) THEN
               READ(385,*) process%NTicks,process%AccruedTicks,process%MDTimeAccrued
            ELSE
               READ(385,*) process%NTicks,process%AccruedTicks
            END IF

            READ(385,*) process%Ticks(1:process%NTicks)
            IF (NoTicks1) process%NTicks=0 !after reading ticks, force to no ticks
            READ(385,*) process%Eact
            READ(385,*) apparent_initial_index,stateinitial_index !state from
              !escape actually occurred, state to which the process belongs
            READ(385,*) statefinal_index
            process%InitialState=>SearchGlobalAMDStateIndex(stateinitial_index)
            process%FinalState=>SearchGlobalAMDStateIndex(statefinal_index)
            process%TPMDActualInitial=>SearchGlobalAMDStateIndex(apparent_initial_index)
            READ(385,*) process%Prefactor
            READ(385,*) process%Rate
            READ(385,*) process%Count
            READ(385,*) process%Index
            READ(385,*) process%TimeAccrued
            CLOSE(385)
            IF (ASSOCIATED(state%process)) THEN
               state%process%PrevNeigh=>process
               process%NextNeigh=>state%process
               state%process=>process
            ELSE
               state%process=>process
            END IF
         END DO
         IF (nprocess>0) WRITE(UNIT=6,FMT='(" Added ",I5," processes to state",I5)') nprocess,istate 
      END DO

      IF (IsPrint1) CALL PrintAMDSummary()
   END SUBROUTINE RecoverStateInformation
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintAMDProcess(process)
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: stateinitial,statefinal
      TYPE(AMDProcessContainer), POINTER :: process
      INTEGER :: NTicks,i

      IF (TPMD) THEN
         WRITE(385,*) process%NTicks,process%AccruedTicks,process%MDTimeAccrued," # stored & overall ticks, MD time accrued"
      ELSE
         WRITE(385,*) process%NTicks,process%AccruedTicks," # stored & overall ticks"
      END IF
      NTicks=process%NTicks
      stateinitial=>process%InitialState
      statefinal=>process%FinalState
      DO i=1,NTicks
         WRITE(UNIT=385,FMT='(ES15.5," ")',ADVANCE="NO") process%Ticks(i)
      END DO
      IF (TPMD) THEN
         WRITE(385,*) "    .... list of tpmd projected tickmarks"
      ELSE
         WRITE(385,*) "    .... list of tickmarks"
      END IF
      WRITE(UNIT=385,FMT='(ES15.5," Activation Barrier")') process%Eact
      IF (TPMD) THEN
         WRITE(UNIT=385,FMT='(2I4," apparent initial state, actual initial state")') &
            stateinitial%Index,process%TPMDActualInitial%Index
      ELSE
         WRITE(UNIT=385,FMT='(I4," Initial state")') stateinitial%Index
      END IF
      WRITE(UNIT=385,FMT='(I4," Final state")') statefinal%Index
      WRITE(UNIT=385,FMT='(ES15.5," Prefactor")') process%Prefactor
      WRITE(UNIT=385,FMT='(ES15.5," Process rate")') process%Rate
      WRITE(UNIT=385,FMT='(I4," Count")') process%Count
      WRITE(385,*) process%Index," process index"
      WRITE(UNIT=385,FMT='(ES15.5," accrued time")') process%TimeAccrued
   END SUBROUTINE PrintAMDProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintAMDSummary()
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state
      TYPE(AMDProcessContainer), POINTER :: process
      TYPE(SuperbasinContainer), POINTER :: sb
      INTEGER :: sbindex

      !print states
      WRITE(6,*) "Number of AMD states:",NGlobalAMDStates
      state=>ListOfGlobalAMDStates
      !WRITE(6,*) "    Index  #processes  superbasin  energy(eV)"
      DO WHILE (ASSOCIATED(state)) 
         sb=>state%superbasin
         sbindex=0
         IF (ASSOCIATED(sb)) sbindex=sb%Index
         WRITE(UNIT=6,FMT='("State, #processes, index, energy (eV)",I5,2I10,E20.8)') &
            state%Index,state%NProcessType,sbindex,state%PotentialEnergy
         IF (state%NProcessType>0) THEN
            WRITE(6,*) "       Process      Ea          nu   ", &
               "        Initialstate    Finalstate"
            process=>state%process
            DO WHILE (ASSOCIATED(process))
               WRITE(UNIT=6,FMT='("      ",i4,2E15.5,"    ",2I10)') &
                 process%Index,process%Eact,process%Prefactor, &
                 process%InitialState%Index,process%FinalState%Index
               process=>process%NextNeigh
            END DO
         END IF
         state=>state%NextNeigh
      END DO

      sb=>ListSuperbasins
      IF (ASSOCIATED(sb)) WRITE(6,*) "SuperbasinIndex  #states"
      DO WHILE (ASSOCIATED(sb))
         WRITE(UNIT=6,FMT='("    ",i4,"       ",i6)') sb%Index,sb%NStates
         sb=>sb%NextNeigh
      END DO
   END SUBROUTINE PrintAMDSummary
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION LocateSBstateInStateMap(StateMap,stateindex)
      IMPLICIT NONE
      INTEGER :: LocateSBstateInStateMap,stateindex
      INTEGER, DIMENSION(:) :: StateMap
      LOGICAL :: NotFound

      LocateSBstateInStateMap=0
      NotFound=.TRUE.
      IF (NotFound) THEN
      END IF
   END FUNCTION LocateSBstateInStateMap
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE InitializeTransitionMatrix(SB,TransitionMatrix,stateprobability,StateMap,errorstatus)
   !sets up the arrays and matrix required to solve the master equation
   !TransitionMatrix,stateprobability,StateMap will be allocated
      IMPLICIT NONE
      TYPE(SuperbasinContainer), POINTER :: SB
      TYPE(SuperbasinStateContainer), POINTER :: state
      REAL(dp), DIMENSION(:,:), POINTER :: TransitionMatrix
      REAL(dp), DIMENSION(:), POINTER :: stateprobability
      INTEGER, DIMENSION(:), POINTER :: StateMap
      INTEGER :: NStates,errorstatus,istate

      state=>SB%state
      NStates=SB%NStates
      errorstatus=0
      IF (NStates>50) THEN !too big SB cannot handle
         errorstatus=1
         RETURN
      END IF

      ALLOCATE(TransitionMatrix(NStates,NStates))
      ALLOCATE(stateprobability(NStates),StateMap(NStates))
      DO WHILE (ASSOCIATED(state))
         istate=istate+1
         StateMap(istate)=state%AMDstate%Index
         state=>state%NextNeigh
      END DO
   END SUBROUTINE InitializeTransitionMatrix
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateTransitionMatrix(SB,TransitionMatrix,StateMap,Temperature)
      IMPLICIT NONE
      TYPE(SuperbasinContainer), POINTER :: SB
      TYPE(SuperbasinStateContainer), POINTER :: state
      TYPE(AMDStateContainer), POINTER :: amdstate
      TYPE(AMDProcessContainer), POINTER :: process
      REAL(dp), DIMENSION(:,:), POINTER :: TransitionMatrix
      INTEGER, DIMENSION(:), POINTER :: StateMap
      REAL(dp) :: Temperature,rate
      INTEGER :: IState,FState,i,j

      TransitionMatrix=0._dp
      DO WHILE (ASSOCIATED(state))
         amdstate=>state%AMDstate
         process=>amdstate%process
         DO WHILE (ASSOCIATED(process))
            IState=process%InitialState%Index
            FState=process%FinalState%Index
            IF (process%Prefactor>0._dp) THEN
               i=LocateSBstateInStateMap(StateMap,IState)
               j=LocateSBstateInStateMap(StateMap,FState)
               rate=process%Prefactor*EXP(-process%Eact/kboltzmann/Temperature)
               TransitionMatrix(i,j)=TransitionMatrix(i,j)+rate
               TransitionMatrix(j,i)=TransitionMatrix(j,i)-rate
            END IF
            process=>process%NextNeigh
         END DO
         state=>state%NextNeigh
      END DO
   END SUBROUTINE CreateTransitionMatrix
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EXPM(matrix)
   !Taylor series based approach to find EXPM, matrix will be replaced with EXP(matrix)
   !matrix is transition matrix
      IMPLICIT NONE
      REAL(dp), DIMENSION(:,:), POINTER :: matrix,matrixpow,expmatrix
      INTEGER :: iter,i,n
      LOGICAL :: NotSatisfied

      n=SIZE(matrix,1)
      ALLOCATE(matrixpow(n,n),expmatrix(n,n))
      DO i=1,n
         expmatrix(i,i)=1._dp !identity matrix
      END DO

      matrixpow=matrix
      NotSatisfied=.TRUE.
      iter=1
      DO WHILE (NotSatisfied)
         expmatrix=expmatrix+matrixpow
         iter=iter+1  !for next iteration
         matrixpow=MATMUL(matrix,matrixpow)/REAL(iter) !for next iteration
         NotSatisfied=ANY(ABS(matrixpow)>0.0001)
      END DO
      matrix=expmatrix
      DEALLOCATE(expmatrix,matrixpow)
      WRITE(6,*) "Exponentiation of matrix converged, # iterations:",iter-1
   END SUBROUTINE EXPM
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CorrectedEscapeTime(EscapeTime,stateprobability,TransitionMatrix,EscapeStateIndex)
   !presently we assume that the equilibrium probility is attained
   !in such a case one just needs to find exp(-TransitionMatrix*EscapeTime)
      IMPLICIT NONE
      REAL(dp) :: EscapeTime
      REAL(dp), DIMENSION(:,:), POINTER :: TransitionMatrix,matrix
      REAL(dp), DIMENSION(:), POINTER :: stateprobability
      INTEGER :: EscapeStateIndex

      TransitionMatrix=-TransitionMatrix*EscapeTime
      CALL EXPM(TransitionMatrix)
      stateprobability=MATMUL(TransitionMatrix,stateprobability)
   END SUBROUTINE CorrectedEscapeTime
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintTransitionMatrix(TransitionMatrix)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:,:), POINTER :: TransitionMatrix
      INTEGER :: n,i,j

      n=SIZE(TransitionMatrix,1)
      WRITE(6,*) "Size of transition matrix is",n
      WRITE(6,*) "Transition matrix is:"
      DO i=1,n
         DO j=1,n
            WRITE(UNIT=6,FMT='(ES10.4)',ADVANCE="NO") TransitionMatrix(i,j)
         END DO
         WRITE(6,*) ""
      END DO
   END SUBROUTINE PrintTransitionMatrix
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintStateProbability(stateprobability)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), POINTER :: stateprobability
      INTEGER :: i

      WRITE(6,*) "State probability:"
      DO i=1,SIZE(stateprobability)
         WRITE(UNIT=6,FMT='(ES15.5)',ADVANCE="NO") stateprobability(i)
      END DO
      WRITE(6,*) ""
      
   END SUBROUTINE PrintStateProbability
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE amd_db_manipulate

