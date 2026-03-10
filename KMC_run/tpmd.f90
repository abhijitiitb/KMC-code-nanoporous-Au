MODULE TPMDModule
!perform temperature programmed molecular dynamics calculations
   USE AMD_VARIABLE_TYPE
   USE amd_db_manipulate
   USE Ecuyer_random
   USE MDPackage
   IMPLICIT NONE
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE TPMDSelectProcess(state,SystemTemperature,dt,TADParam,MDRequired)
      IMPLICIT NONE
      INTEGER, PARAMETER :: MaxCountSBProcess=4
      TYPE(AMDStateContainer), POINTER :: state,finalstate
      TYPE(AMDProcessContainer), POINTER :: process,process_selected,process1
      TYPE(SuperbasinContainer), POINTER :: sb,sb1,sb2
      TYPE(SuperbasinstateContainer), POINTER :: SBstate
      TYPE(TADParameters), POINTER, OPTIONAL :: TADParam
      REAL(dp) :: SystemTemperature,NewBasinTime,timemin,Ew,nu_min,alfa,dt
      REAL(dp) :: StageTemperature,StageTime,k0,timeprojected,kstage
      CHARACTER(len=100) :: filename
      CHARACTER(len=NINT2CHAR) :: charstateindx,charprocessindx
      INTEGER :: errorstatus,NTicks,stateindx,SBindx,finalSBindx
      LOGICAL :: MDRequired,IsSBProcess,IsFinalStateNotSBstate

      !find shortest time
      !in case of superbasin ... situation is bit more complicated
      ! escapes between superbasin states are to be ignored
      ! the initial state matters, because the superbasin might be large
      ! hence the probability evolution will depend on initial condition
      process=>state%Process
      NULLIFY(process_selected)
      timemin=9999._dp
      DO WHILE (ASSOCIATED(process))
         IsFinalStateNotSBstate=.TRUE.
         IF (TPMDSuperbasin) THEN !find the superbasin state of initial and final state
            sb1=>state%superbasin
            sb2=>process%FinalState%superbasin
            IF (ASSOCIATED(sb1) .AND. ASSOCIATED(sb2)) IsFinalStateNotSBstate=sb1%Index/=sb2%Index
            !this way ticks between superbasin states can be ignored using IsFinalStateNotSBstate
         END IF
         IF (IsFinalStateNotSBstate) THEN
            !check if the activation barrier is large
            !small activation barrier could result in small time step
            !this becomes a problem since KMC time step will then be 
            !negative
            IF (process%Eact>0.1_dp) THEN !allow process
               IF (process%NTicks>0) THEN
                  IF (timemin>process%Ticks(1)) THEN
                     process_selected=>process
                     timemin=process%Ticks(1)
                  END IF
               ELSEIF (process%NTicks==0) THEN
                  MDRequired=.TRUE.
                  dt=0._dp
                  RETURN
               ELSE
                  WRITE(6,*) "Err>> ... nticks is less than zero"
                  STOP
               END IF
            ELSE
               WRITE(6,*) "Err>> Activation barrier of process is too small"
               WRITE(6,*) " ... ignoring process"
            END IF
         END IF
         process=>process%NextNeigh
      END DO

      dt=- state%BasinEscapeTime !time elapsed is Final-Initial, we are initializing dt = -Initial
      !in TPMD the BasinEscapeTime denotes the reference time
      !all processes tickmarks are provided with respect to this reference time
      !for instance if process tick mark is 1.53e-6 s and BasinEscapeTime is 3.2e-9 s
      ! then the process will happen after 1.53e-6 - 3.2e-9 s
      ! as the dynamics is advanced, the BasinEscapeTime is advanced to the new time
      ! Comment 1. in essence this tells us the duration for which we have seen a basin
      !Comment 1 is not valid when a superbasin is formed
      ! when two superbasins are merged, their process times need to be placed on the same time axis
      ! first all superbasin states are required to have the same BasinEscapeTime
      ! this helps while reusing the ticks as superbasins are created by adjusting the process times
      ! to the common time axis
      !Comment 2. The KMC dynamics is more interested in the changes in the BasinEscapeTime 
      ! and not the actual process times
      
      !find the Ew -- worst case barrier
      alfa=0.1_dp
      nu_min=1.e12
      Ew= - kboltzmann*SystemTemperature* LOG( -LOG(alfa)/nu_min/timemin )
      WRITE(6,*) "Worst case barrier:",Ew
      CALL FLUSH(6)

      charstateindx=INT2CHAR(state%Index)

      !printing info
      WRITE(6,*) "Printing process files ..."
      CALL FLUSH(6)
      CALL SYSTEM("cat ./GlobalAMD/Process."//TRIM(charstateindx)//".*")
      CALL FLUSH(6)
      !WRITE(6,*) "Printing MD param files ..."
      !CALL FLUSH(6)
      !CALL SYSTEM("cat ./GlobalAMD/MD.1*.param")
      !CALL FLUSH(6)
      !WRITE(6,*) "Priting MD state file ..."
      !CALL FLUSH(6)
      !CALL SYSTEM("cat ./GlobalAMD/state.1.info")
      !CALL FLUSH(6)
      !WRITE(6,*) "Priting tpmd stages ..."
      !CALL FLUSH(6)
      !CALL SYSTEM("cat ./GlobalAMD/state.1.tpmdstages")
      !CALL FLUSH(6)

      !read tpmd stages and find projected time
      filename="./GlobalAMD/state."//TRIM(charstateindx)//".tpmdstages"
      errorstatus=0
      k0=nu_min*EXP(-Ew/kboltzmann/SystemTemperature)
      timeprojected=0._dp
      OPEN(UNIT=540,FILE=TRIM(filename))
      WRITE(*,*) TRIM(filename)
      DO WHILE (errorstatus==0) 
         READ(UNIT=540,FMT=*,IOSTAT=errorstatus) StageTemperature,StageTime
         IF (errorstatus==0) THEN
      !      WRITE(*,*) StageTemperature,StageTime
            kstage=nu_min*EXP(-Ew/kboltzmann/StageTemperature)
            timeprojected=timeprojected+kstage*StageTime/k0 !time when the worst case process is to be observed
         END IF
      END DO
      CLOSE(540)
      WRITE(6,*) "Projected worst case time, min time:",timeprojected,timemin
      
      !determine if additional MD is required
      MDRequired=.TRUE.
      IF (timeprojected>timemin) THEN !can perform the move
         MDRequired=.FALSE.
         
         !update time
         NewBasinTime=timemin
         state%BasinEscapeTime=NewBasinTime !updated basin time - when last escape occurred
         dt=dt+NewBasinTime !time increment for this basin
         
         !update process count & check if state part of new SB
         process_selected%Count=process_selected%Count+1
         process_selected%TimeAccrued=timemin
         NTicks=process_selected%NTicks
         WRITE(6,*) "Earlier ticks:",process_selected%Ticks(1:process_selected%NTicks)
         process_selected%Ticks(1:NTicks-1)=process_selected%Ticks(2:NTicks)
         process_selected%Ticks(NTicks:MaxNTicks)=1000000._dp
         process_selected%NTicks=MAX(0,NTicks-1) !number of ticks currently available
         !process_selected%AccruedTicks=MAX(0,process_selected%AccruedTicks-1)  -- AccruedTicks is total ticks collected by MD
         WRITE(6,*) "New ticks:",process_selected%Ticks(1:process_selected%NTicks)
         CALL FLUSH(6)
         !MDRequired= process_selected%NTicks==0 !since the selected process occurs at these times 

         !check if selected process is part of superbasin
         IF (TPMDSuperbasin .AND. process_selected%Count>MaxCountSBProcess) THEN
            WRITE(6,*) "Checking whether process belongs to superbasin ..."
            
            !determine existing SB of starting state
            stateindx=state%Index
            SBindx=0
            IF (ASSOCIATED(state%superbasin)) SBindx=state%superbasin%Index

            !determine existing SB of final state
            finalstate=>process_selected%FinalState
            finalSBindx=0
            IF (ASSOCIATED(finalstate%superbasin)) finalSBindx=finalstate%superbasin%Index
            WRITE(6,*) "State & superbasin:",stateindx,SBindx
            WRITE(6,*) "State & superbasin:",finalstate%Index,finalSBindx
            CALL FLUSH(6)

            IF (SBindx==0 .OR. finalSBindx==0 .OR. SBindx/=finalSBindx) THEN !current perception is the two states dont belong to same SB
               !determine whether a SB exist between the two processes
               process1=>finalstate%Process
               IsSBProcess=.FALSE.
               DO WHILE (ASSOCIATED(process1))
                  IF (process1%FinalState%Index==stateindx .AND. process1%Count>MaxCountSBProcess) THEN
                     IsSBProcess=.TRUE.
                     EXIT
                  END IF
                  process1=>process1%NextNeigh
               END DO

               !if superbasin is present make updates
               IF (IsSBProcess) THEN
                  WRITE(6,*) "Following states are part of superbasin:",stateindx,finalstate%Index
                  IF (SBindx==0) THEN
                     IF (finalSBindx==0) THEN
                        NULLIFY(sb)
                        ALLOCATE(sb)
                        CALL SYSTEM("cat ./GlobalAMD/state."//TRIM(charstateindx)//".info")
                        CALL SYSTEM("cat ./GlobalAMD/Process."//TRIM(charstateindx)//".*")
                        WRITE(6,*) "Adding two states to a new SB"
                        CALL AddSuperBasinState(state,sb)
                        CALL AddSuperBasinState(finalstate,sb)
                        CALL AddSuperBasin(sb)
                        CALL WriteSuperbasinInfo(sb)
                        CALL SYSTEM("cat ./GlobalAMD/state."//TRIM(charstateindx)//".info")
                        CALL SYSTEM("cat ./GlobalAMD/Process."//TRIM(charstateindx)//".*")
                        CALL FLUSH(6)
                     ELSE !finalSBindx>0
                        CALL SYSTEM("cat ./GlobalAMD/state."//TRIM(charstateindx)//".info")
                        CALL SYSTEM("cat ./GlobalAMD/Process."//TRIM(charstateindx)//".*")
                        WRITE(6,*) "Adding initial state to SB of final state"
                        sb=>finalstate%superbasin
                        CALL AddSuperBasinState(state,sb)
                        CALL WriteSuperbasinInfo(sb)
                        CALL SYSTEM("cat ./GlobalAMD/state."//TRIM(charstateindx)//".info")
                        CALL SYSTEM("cat ./GlobalAMD/Process."//TRIM(charstateindx)//".*")
                        CALL FLUSH(6)
                     END IF
                  ELSE !SBindx>0
                     IF (finalSBindx==0) THEN
                        CALL SYSTEM("cat ./GlobalAMD/state."//TRIM(charstateindx)//".info")
                        CALL SYSTEM("cat ./GlobalAMD/Process."//TRIM(charstateindx)//".*")
                        WRITE(6,*) "Adding final state to SB of initial state"
                        sb=>state%superbasin
                        CALL AddSuperBasinState(finalstate,sb)
                        CALL WriteSuperbasinInfo(sb)
                        CALL SYSTEM("cat ./GlobalAMD/state."//TRIM(charstateindx)//".info")
                        CALL SYSTEM("cat ./GlobalAMD/Process."//TRIM(charstateindx)//".*")
                        CALL FLUSH(6)
                     ELSE !finalSBindx>0 and is different from SBindx
                        CALL SYSTEM("cat ./GlobalAMD/state."//TRIM(charstateindx)//".info")
                        CALL SYSTEM("cat ./GlobalAMD/Process."//TRIM(charstateindx)//".*")
                        WRITE(6,*) "Merging SB of initial & final state"
                        CALL MergeSuperbasins(state,finalstate)
                        sb=>state%superbasin
                        CALL WriteSuperbasinInfo(sb)
                        CALL SYSTEM("cat ./GlobalAMD/state."//TRIM(charstateindx)//".info")
                        CALL SYSTEM("cat ./GlobalAMD/Process."//TRIM(charstateindx)//".*")
                        CALL FLUSH(6)
                     END IF
                  END IF
                  WRITE(6,*) "Summary of superbasin info:"
                  WRITE(6,*) "Initial state belongs to: ",state%superbasin%Index
                  WRITE(6,*) "Final   state belongs to: ",finalstate%superbasin%Index
                  WRITE(6,*) "# Superbasin states:",state%superbasin%Nstates
                  WRITE(6,*) "Superbasin index:",state%superbasin%Index
                  WRITE(UNIT=6,FMT='("Following states are present in SB:")',ADVANCE="NO")
                  SBstate=>state%superbasin%state
                  DO WHILE (ASSOCIATED(SBstate))
                     WRITE(UNIT=6,FMT='(I4)',ADVANCE="NO") SBstate%AMDstate%Index
                     SBstate=>SBstate%NextNeigh
                  END DO
                  WRITE(6,*) ""
                  CALL WriteGlobalAMDStateInfo(state)
                  CALL WriteGlobalAMDStateInfo(finalstate)
                  CALL FLUSH(6)
               END IF
            END IF
         END IF

         !generating artifical ticks can be risky
         !IF (process_selected%NTicks==0) THEN
         !   process_selected%NTicks=1
         !   k0=process_selected%prefactor*EXP(-process_selected%Eact/kboltzmann/SystemTemperature)
         !   process_selected%Ticks(1)=-LOG(taus88())/k0
         !   WRITE(6,*) "Increment for artificial tick:",process_selected%Ticks(1)
         !   !generating ticks artificially is perfectly safe. the worst case time depends on tpmdstages
         !   !if the worst case time is less than the tick generated we need to rebuild the catalog
         !   !in some cases a new tick from the updated catalog could be less than the basin escape time
         !   process_selected%Ticks(1)=process_selected%TimeAccrued+process_selected%Ticks(1)
         !   process_selected%TimeAccrued=process_selected%Ticks(1)
         !   WRITE(6,*) "Generated a tickmark artificially @ time:", process_selected%TimeAccrued
         !END IF

         charprocessindx=INT2CHAR(process_selected%Index)
         filename="./GlobalAMD/Process."//TRIM(charstateindx)//"."//TRIM(charprocessindx)
         OPEN(UNIT=385,FILE=filename)
         CALL PrintAMDProcess(process_selected)
         CALL FLUSH(6)
         CLOSE(385)
         CALL FLUSH(6)
         !   process=>state%Process
         !   DO WHILE (ASSOCIATED(process))
         !      write(6,*) "Proces index, ticks:",process%Index,process%NTicks
         !      IF (process%NTicks<=0 .AND. process%Rate>0._dp) THEN
         !         WRITE(6,*) "... oops no ticks found for", process%Index
         !         STOP
         !      END IF
         !      process=>process%NextNeigh
         !   END DO
         WRITE(6,*) "Selected process index is:",process_selected%Index
         WRITE(6,*) "New basin time (this could be modified after merging with new SB):",state%BasinEscapeTime
         state=>process_selected%FinalState
         WRITE(6,*) "New state is index:",state%Index
      ELSE !cannot perform escape
         dt=0._dp
      END IF
      CALL FLUSH(6)

   END SUBROUTINE TPMDSelectProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE StateTPMD(state,SystemTemperature,Temperature,TemperatureProgram,DetectTransitionFrequency, &
      prefacmin,delta,imode)
   !written in 2012
   !Performs TPMD calculations for the state using one search at a time
   !SystemTemperature using the TPMD validity
   !imode=0 regular TPMD
   !imode=1 e-TPMD
   !Temperature = starting temperature for TPMD
   !TemperatureProgram = temperature programming to be used for TPMD
   !nsearch = number of MD searches performed
   !prefacmin = minimum prefactor in TPMD
   !delta = error in TPMD
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state
      REAL(dp) :: SystemTemperature,Temperature,prefacmin,delta,MDTime
      INTEGER :: DetectTransitionFrequency,imode,ObtainRate,isearch,TemperatureProgram
      LOGICAL :: success
      
      !step 1: Perform MD search
      isearch=1
      success=.FALSE.
      MDTime=0._dp
      ObtainRate=1
      CALL SearchProcessesFromState(state%index,Temperature,DetectTransitionFrequency,isearch, &
         success,MDTime,ObtainRate,TemperatureProgram)
      
      !step 2: Project escape times to escape times for the SystemTemperature
      
      !step 3: decide how many escapes are within the confidence level
      
   END SUBROUTINE StateTPMD
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE TPMDModule
