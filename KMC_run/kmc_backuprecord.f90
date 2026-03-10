MODULE KMCRecords !Level 5
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
   USE KMC_VARIABLE_TYPE
   USE KMCUtilities, ONLY : INT2CHAR
   USE ModFile, ONLY : OpenFile
   USE KMCDbManipulate
   USE KMCConfigurations
   IMPLICIT NONE

   LOGICAL, PARAMETER :: PrintCfgXYZ=.TRUE.,PrintPrcXYZ=.TRUE.
   INTEGER :: CounterCfgFileNumber=0,CounterProcessFileNumber=0
   SAVE

   CONTAINS

!--------------------------------------------------------------------------------
!              Recovering and recording of system status
!--------------------------------------------------------------------------------
SUBROUTINE GenSystemRecord()
   IMPLICIT NONE
   TYPE(Configuration), POINTER :: config
   TYPE(Bin), POINTER :: bin1
   TYPE(BinAtom), POINTER :: binatom1
   TYPE(KMCProcess), POINTER :: prc
   TYPE(KMCProcessAtom), POINTER :: processatom1
   TYPE(KMCProcessList), POINTER :: process1
   TYPE(KMCAtomList), POINTER :: AL
   INTEGER :: indx,indx1,RecFileNumber,iStore
   INTEGER :: n,NEnvironmentAtoms
   LOGICAL :: IsOld
   CHARACTER(len=NINT2CHAR) :: ChO
   CHARACTER(len=250) :: FileName

   CALL SYSTEM('cp '//TRIM(FileProcessListRec)//' '//TRIM(FileProcessListRec)//'.old')
   CALL SYSTEM('cp '//TRIM(FileStoreState)//' '//TRIM(FileStoreState)//'.old')

   OPEN(UNIT=UnitProcessListRec,FILE=TRIM(FileProcessListRec))
   WRITE(6,*) 'Recr>>Storing systems'
   WRITE(UnitProcessListRec,*) KMCTime,SnapCounter,KMCSumProcessRates,' [A1]'
   WRITE(UnitProcessListRec,*) CPUTime(2)+CPUStartValue(2),CPUTime(5)+CPUStartValue(5), &
     CPUTime(6)+CPUStartValue(6),' [A2]'
   WRITE(UnitProcessListRec,*) NumberCfgTypes,NumberProcesses,'  [A3]'
   WRITE(UnitProcessListRec,*) KMCBlock,KMCStep,idum,' [A4]'
   WRITE(UnitProcessListRec,*) TotalCalls(2),TotalCalls(5),TotalCalls(7),' [A5]'
   WRITE(UnitProcessListRec,*) Temperature,TADLoTemperature,TADHiTemperature,' [A6]'
   WRITE(UnitProcessListRec,*) SnapCounter,NAtomKMC,NMoveKMC,' [A7]'
   CLOSE(UnitProcessListRec)
   
   !----------------CONFIGURATIONS------------------
   !reach the end of NewCfgList
   IF (ASSOCIATED(NewCfgList)) THEN
      DO WHILE (ASSOCIATED(NewCfgList%NextNeigh))
         NewCfgList=>NewCfgList%NextNeigh
      END DO
   END IF
   
   iStore=0
   DO WHILE (ASSOCIATED(NewCfgList)) !need to store in ascending order, so also need to provide config before this config
      config=>NewCfgList%Cfg
      IF (config%RecFileNumber>0) THEN
         RecFileNumber=config%RecFileNumber
         IsOld=.TRUE.
      ELSE
         CounterCfgFileNumber=CounterCfgFileNumber+1
         RecFileNumber=CounterCfgFileNumber
         config%RecFileNumber=RecFileNumber
         IsOld=.FALSE.
      END IF
      ChO=INT2CHAR(RecFileNumber)
      FileName='./StoreTADKMC/SystemConfigs/C'//TRIM(ChO)
      IF (IsOld) CALL SYSTEM('cp '//TRIM(FileName)//' '//TRIM(FileName)//'.old')
      OPEN(UNIT=UnitProcessListRec,FILE=TRIM(FileName))
      WRITE(UNIT=UnitProcessListRec,FMT='(3I10,"  [C1]")') config%ConfigurationIndex,KMCBlock,KMCStep
      WRITE(UnitProcessListRec,*) config%Species,' [C2]'
      WRITE(UnitProcessListRec,*) config%Env%ShortRangeConfig%Population(1:NSpeciesType),' [C3]'
      WRITE(UnitProcessListRec,*) config%BondOrder,config%NAtomsWithCfg,' [C4]'
      
      bin1=>config%Env%ShortRangeConfig%BinInfo
      NEnvironmentAtoms=0
      DO WHILE (ASSOCIATED(bin1))
         WRITE(UNIT=UnitProcessListRec,FMT='(I10," [C5]")') bin1%BinNumber
         DO indx=1,NSpeciesType
            WRITE(UNIT=UnitProcessListRec,FMT='(3I10," [C6]")') bin1%InBin(indx), &
              bin1%BinEdgeP(indx),bin1%BinEdgeN(indx)
         END DO
         binatom1=>bin1%AL
         DO WHILE (ASSOCIATED(binatom1))
            NEnvironmentAtoms=NEnvironmentAtoms+1
            WRITE(UNIT=UnitProcessListRec,FMT='(4F18.4,I10," [C7]")') binatom1%distance, &
              binatom1%RelCoord,binatom1%Species
            binatom1=>binatom1%NextNeigh
         END DO
         WRITE(UNIT=UnitProcessListRec,FMT='(4F18.4,I10," [C7]")') 0.0,0.0,0.0,0.0,0
         bin1=>bin1%NextNeigh
      END DO
      WRITE(UnitProcessListRec,*) '0  [C5]'
      CLOSE(UnitProcessListRec)
      IF (IsOld) CALL SYSTEM('rm '//TRIM(FileName)//'.old')
      
      IF (PrintCfgXYZ) THEN
         NEnvironmentAtoms=NEnvironmentAtoms+1 !add the center atom as well
         OPEN(UNIT=UnitProcessListRec,FILE=TRIM(FileName)//'.xyz')
         WRITE(UnitProcessListRec,*) NEnvironmentAtoms
         WRITE(UnitProcessListRec,*) "BO:",config%BondOrder
         WRITE(UnitProcessListRec,* ) SpeciesName(config%Species+1),0.,0.,0.
         bin1=>config%Env%ShortRangeConfig%BinInfo
         DO WHILE (ASSOCIATED(bin1))
            binatom1=>bin1%AL
            DO WHILE (ASSOCIATED(binatom1))
               WRITE(UnitProcessListRec,*) SpeciesName(binatom1%Species),binatom1%RelCoord
               binatom1=>binatom1%NextNeigh
            END DO
            bin1=>bin1%NextNeigh
         END DO
         !Below the atoms that belong to this configuration are printed to the file. However, this is now commented out because in any case the cfg is not added to the new list when atoms join/leave this cfg. See subroutines AddToCfg and RemoveFromCfg where CALL AddCfgToNewList(config) is commented out 
         !WRITE(UnitProcessListRec,*) ""
         !WRITE(UnitProcessListRec,*) "Atoms with this cfg:"
         !AL=>config%AtomsWithCfg
         !DO WHILE (ASSOCIATED(AL))
         !   WRITE(UnitProcessListRec,*) AL%Atom%Index,AL%Atom%Coord
         !   AL=>AL%NextNeigh
         !END DO
         CLOSE(UnitProcessListRec)
      END IF
      
      IF (ASSOCIATED(NewCfgList%PrevNeigh)) THEN
         NULLIFY(NewCfgList%Cfg)
         NewCfgList=>NewCfgList%PrevNeigh
         NULLIFY(NewCfgList%NextNeigh%PrevNeigh)
         DEALLOCATE(NewCfgList%NextNeigh); NAllottedConfigurationList=NAllottedConfigurationList-1
      ELSE
         NULLIFY(NewCfgList%Cfg)
         DEALLOCATE(NewCfgList); NAllottedConfigurationList=NAllottedConfigurationList-1
      END IF
      iStore=iStore+1
   END DO  !ASSOCIATED(config)
   WRITE(UnitScrap,*) "Gen>># configs stored:",istore
   WRITE(6,*) "Gen>># configs stored:",istore
   !--------------END OF CONFIGURATIONS---------------------
   !--------------------PROCESSES---------------------------
   istore=0

   process1=>NewProcessList
   IF (ASSOCIATED(process1)) WRITE(UnitScrap,*) "FOLLOWING ARE THE PROCESSES THE NEED TO BE RECORDED"
   DO WHILE (ASSOCIATED(process1))
      WRITE(UNIT=UnitScrap,FMT='(I5,",")',ADVANCE="NO") process1%Process%Index
      process1=>process1%NextNeigh
   END DO

   DO WHILE (ASSOCIATED(NewProcessList))
      process1=>NewProcessList
      prc=>process1%Process
      RecFileNumber=prc%RecordNumber
      IF (RecFileNumber<=0 .OR. prc%Index<=0) THEN
         WRITE(6,*) "Err>> While writing processes, record number is zero "
         WRITE(6,*) "Process index:",prc%Index
         WRITE(6,*) "Process record number:",prc%RecordNumber
         CALL PrintProcess(prc,UnitScrap)
         WRITE(6,*) "Total number of processes so far:",NumberProcesses
         CLOSE(UnitScrap)
         STOP
      END IF 
      ChO=INT2CHAR(RecFileNumber)
      FileName='./StoreTADKMC/SystemProcesses/P'//TRIM(ChO)
      IF (IsOld) CALL SYSTEM('cp '//TRIM(FileName)//' '//TRIM(FileName)//'.old')
      OPEN(UNIT=UnitProcessListRec,FILE=TRIM(FileName))
      WRITE(UNIT=UnitProcessListRec,FMT='(I10," [P1]")') prc%NumberProcessAtoms
      WRITE(UNIT=UnitProcessListRec,FMT='(I10," [P1a]")') prc%CfgType%ConfigurationIndex
      WRITE(UNIT=UnitProcessListRec,FMT='(2I18," [P2]")') prc%FiringCount,prc%MaxFiringAllowed
      WRITE(UNIT=UnitProcessListRec,FMT='(2ES18.4,F18.4," [P3]")') prc%Rate,prc%Frequency, prc%FwdBarrier
      processatom1=>prc%ReactantAL
      DO WHILE (ASSOCIATED(processatom1))
         WRITE(UNIT=UnitProcessListRec,FMT='(3F18.4,"  [P4]")') processatom1%RelInitialPos
         WRITE(UNIT=UnitProcessListRec,FMT='(3F18.4,"  [P4a]")') processatom1%RelSaddlePos
         WRITE(UNIT=UnitProcessListRec,FMT='(3F18.4,"  [P5]")') processatom1%RelFinalPos
         WRITE(UNIT=UnitProcessListRec,FMT='(2I10,"  [P6]")') processatom1%CfgIndex,processatom1%Species
         processatom1=>processatom1%NextNeigh
      END DO
      WRITE(UnitProcessListRec,*) '0    [P1]'
      CLOSE(UnitProcessListRec)
      IF (IsOld) CALL SYSTEM('rm '//TRIM(FileName)//'.old')
      
      IF (PrintPrcXYZ) THEN
         OPEN(UNIT=UnitProcessListRec,FILE=TRIM(FileName)//".xyz")
         WRITE(UNIT=UnitProcessListRec,FMT='(I10," [P1]")') prc%NumberProcessAtoms
         WRITE(UnitProcessListRec,*) "Initial state"
         processatom1=>prc%ReactantAL
         DO WHILE (ASSOCIATED(processatom1))
            WRITE(UnitProcessListRec,*) SpeciesName(processatom1%Species),processatom1%RelInitialPos
            processatom1=>processatom1%NextNeigh
         END DO
         WRITE(UNIT=UnitProcessListRec,FMT='(I10," [P1]")') prc%NumberProcessAtoms
         WRITE(UnitProcessListRec,*) "Saddle state"
         processatom1=>prc%ReactantAL
         DO WHILE (ASSOCIATED(processatom1))
            WRITE(UnitProcessListRec,*) SpeciesName(processatom1%Species),processatom1%RelSaddlePos
            processatom1=>processatom1%NextNeigh
         END DO
         WRITE(UNIT=UnitProcessListRec,FMT='(I10," [P1]")') prc%NumberProcessAtoms
         WRITE(UnitProcessListRec,*) "Final state"
         processatom1=>prc%ReactantAL
         DO WHILE (ASSOCIATED(processatom1))
            WRITE(UnitProcessListRec,*) SpeciesName(processatom1%Species),processatom1%RelFinalPos
            processatom1=>processatom1%NextNeigh
         END DO
         CLOSE(UnitProcessListRec)
      END IF
      
      IF (ASSOCIATED(NewProcessList%NextNeigh)) THEN
         NULLIFY(NewProcessList%Process)
         NewProcessList=>NewProcessList%NextNeigh
         NULLIFY(NewProcessList%PrevNeigh%NextNeigh)
         DEALLOCATE(NewProcessList%PrevNeigh); NAllottedKMCProcessList=NAllottedKMCProcessList-1
      ELSE
         NULLIFY(NewProcessList%Process)
         DEALLOCATE(NewProcessList); NAllottedKMCProcessList=NAllottedKMCProcessList-1
      END IF
      istore=istore+1
   END DO
   
   WRITE(UnitScrap,*) "Gen>># processes stored:",istore
   WRITE(6,*) "Gen>># processes stored:",istore
   !--------------END OF OFFICIALS---------------------
   WRITE(6,*) "Gen>> Succesfully Stored system"
   
!ADVANTAGES OF DELETING CONFIGS AND PROCESSES
! 1) Relevant processes/cfgs are always stored in memory and can be used quickly
! 2) All processes/cfgs are stored
! 3) Additional cost required in storing process atoms, as more new processes are found the older processes
!   are observed less often and they end up wasting memory due to the stored process atoms

!   DETAILED COMPARISON OF MEMORY SAVED BY DELETING PROCESS/CFG INFORMATION
!   Memory requirement (approx. POINTER(1 byte), REAL(4 byte), INTEGER(2 byte))

!   ATOM PROCESS (with process not in active use, hence there are no subsrcibers to process)
!        Variable                        Deleting all info (bytes saved)        Deleting process atoms
!     NextNeigh, PrevNeigh                          x (2)
!     ReactantAL (assuming 20 process atoms)        x (400)                          x (400)
!     CfgType                                       x (1)
!     APLPosition                                   x (1)
!     NumberProcessAtoms                            x (2)
!     NumebrSubscribers                             x (2)
!     Index,FiringCount,MaxFiringAllowed            x (6)
!     Rate,Freq,Barrier                             x (12)
!     LastUpdate,FirstAdded                         x (4)
!     ---------------------------------------------------------------------------------------------
!     TOTAL SAVED PER PROCESS DELETE                   430 B                             400 B
!     TOTAL SAVED FOR 10000 PROCESSES                  4 MB
!     TOTAL SAVED FOR 10000000 PROCESSES               4 GB

!   ATOM CONFIGURATION (with cfgs not in use, AtomsWithCfg is NULL)
!        Variable                        Deleting all info (bytes saved)        Deleting process atoms
!     NextNeigh, PrevNeigh, AtomsWithCfg            x (3)
!     Env (assuming 30 atoms in env)                x (600)                          x (600)
!     Others                                        x  (20)
!     ----------------------------------------------------------------------------------------------
!     TOTAL SAVED PER CFG DELETED                   630 B                                600 B
!     TOTAL SAVED FOR 10000 CFGS                    6 MB
!     TOTAL SAVED FOR 10000000 CFGS                 6 GB


   !First delete processes
   !Next delete configurations
   !IF (NumberCfgTypes>MaxCfgInDatabase) THEN
   !   WRITE(6,*) "Gen>> Attempting to delete redundant configurations"
   !   NDeleted=0
   !   CurrKMCStep=KMCStep+KMCBlock*NKMCStepsPeBlock
   !   DeleteBeforeKMCStep=CurrKMCStep-NDeleteWindow

   !   cfg=>
   !   WRITE(6,*) "Gen>> Number of configurations deleted:",NDeleted
   !END IF
   !WRITE(6,*)  
END SUBROUTINE GenSystemRecord
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE RecvrSystemRecord !called at the begining of the run to read stored configs
   IMPLICIT NONE
   TYPE(Configuration), POINTER :: config,config1
   TYPE(Bin), POINTER :: bin1
   TYPE(BinAtom), POINTER :: binatom1
   TYPE(KMCProcess), POINTER :: prc
   TYPE(KMCProcessAtom), POINTER :: prcatm
   INTEGER :: indx,indx1,SizeEnvironment
   INTEGER :: errorstatus
   INTEGER :: TempLoRec,SnapCounterVal
   INTEGER :: ConfigurationIndex,NAtomKMCrec,NMoveKMCRec
   INTEGER :: BinNumber,SpeciesIndex,NoAtomsInvolved,NumberProcesses1
   LOGICAL :: NotSatisfied,NotSatisfied1 !,IsPresent
   REAL :: distance,Coord(3)
   CHARACTER(len=NINT2CHAR) :: ChO
   CHARACTER(len=250) :: FileName
   
   OPEN(UNIT=UnitProcessListRec,FILE=TRIM(FileProcessListRec))
   WRITE(6,*) 'Recv>>Reading list of configurations & processes encountered earlier'
   WRITE(UnitScrap,*) 'Recv>>Reading list of configurations & processes encountered earlier'

   READ(UnitProcessListRec,*) KMCTime,SnapCounterVal ![A1]
   READ(UnitProcessListRec,*) CPUStartValue(2),CPUStartValue(5),CPUStartValue(6) ![A2]
   READ(UnitProcessListRec,*) NumberCfgTypes,NumberProcesses1 ![A3]
   READ(UnitProcessListRec,*) KMCBlock,KMCStep,idum ![A4]
   READ(UnitProcessListRec,*) TotalCalls(2),TotalCalls(5),TotalCalls(7) ![A5]
   READ(UnitProcessListRec,*) Temperature  ![A6]
   READ(UnitProcessListRec,*) SnapCounter,NAtomKMCrec,NMoveKMCrec ![A7]
   CLOSE(UnitProcessListRec)
   CounterCfgFileNumber=NumberCfgTypes

   IF (ResetNoFirings) THEN
      KMCTime=0._dp
      CPUStartValue(5)=0.
      KMCBlock=0
      KMCStep=0
      TotalCalls(5)=0
      idum=ABS(idum)
      SnapCounter=0
   ELSE
      idum=ABS(idum)+46+KMCStep+KMCBlock*NKMCStepsPerBlock
      IF (Restart) THEN !start from where left previously
         NAtomKMC=NAtomKMCrec
         NMoveKMC=NMoveKMCrec
         SnapCounter=SnapCounterVal
      ELSE !reset values
         KMCTime=0._dp
         CPUStartValue(5)=0.
         KMCBlock=0
         KMCStep=0
         TotalCalls(5)=0
         idum=ABS(idum)
         SnapCounter=0
      END IF
   END IF
   
   !----------------------READ ATOM ENVIRONMENTS----------------------
   IF (.NOT. ASSOCIATED(RecordedCfg)) CALL CreateRedundantCfg()
   DO indx=1,NumberCfgTypes
      ChO=INT2CHAR(indx)
      FileName='./StoreTADKMC/SystemConfigs/C'//TRIM(ChO)
      OPEN(UNIT=UnitProcessListRec,FILE=TRIM(FileName))
      READ(UnitProcessListRec,*) ConfigurationIndex ![C1]
      ALLOCATE(config); NAllottedConfiguration=NAllottedConfiguration+1
      
      config%ConfigurationIndex=ConfigurationIndex
      config%RecFileNumber=indx
      READ(UnitProcessListRec,*) config%Species ![C2]
      ALLOCATE(config%Env); NAllottedEnvironment=NAllottedEnvironment+1
      ALLOCATE(config%Env%ShortRangeConfig); NAllottedHistogram=NAllottedHistogram+1
      READ(UnitProcessListRec,*) config%Env%ShortRangeConfig%Population(1:NSpeciesType) ![C3]
      READ(UnitProcessListRec,*) config%BondOrder ![C4]
      
      NotSatisfied=.TRUE.
      NULLIFY(bin1)
      DO WHILE (NotSatisfied) !start of bin
         READ(UnitProcessListRec,*) BinNumber ![C5]
         IF (BinNumber==0) EXIT !equivalent of NotSatisfied=.FALSE.
         IF (.NOT. ASSOCIATED(config%Env%ShortRangeConfig%BinInfo)) THEN
            ALLOCATE(config%Env%ShortRangeConfig%BinInfo); NAllottedBin=NAllottedBin+1
            bin1=>config%Env%ShortRangeConfig%BinInfo
            NULLIFY(bin1%NextNeigh)
            NULLIFY(bin1%PrevNeigh)
         ELSE
            ALLOCATE(bin1%NextNeigh); NAllottedBin=NAllottedBin+1
            bin1%NextNeigh%PrevNeigh=>bin1
            bin1=>bin1%NextNeigh
            NULLIFY(bin1%NextNeigh)
         END IF 
         bin1%BinNumber=BinNumber
         DO indx1=1,NSpeciesType
            READ(UnitProcessListRec,*) bin1%InBin(indx1),bin1%BinEdgeP(indx1),bin1%BinEdgeN(indx1) ![C6]
         END DO
         NotSatisfied1=.TRUE.
         NULLIFY(binatom1)
         DO WHILE (NotSatisfied1)
            READ(UnitProcessListRec,*) distance,Coord,SpeciesIndex ![C7]
            IF (distance==0.) EXIT
            IF (.NOT. ASSOCIATED(bin1%AL)) THEN
               ALLOCATE(bin1%AL); NAllottedBinAtom=NAllottedBinAtom+1
               binatom1=>bin1%AL
               NULLIFY(binatom1%NextNeigh)
               NULLIFY(binatom1%PrevNeigh)
            ELSE
               ALLOCATE(binatom1%NextNeigh); NAllottedBinAtom=NAllottedBinAtom+1
               binatom1%NextNeigh%PrevNeigh=>binatom1
               binatom1=>binatom1%NextNeigh
               NULLIFY(binatom1%NextNeigh)
            END IF !end of atm list
            binatom1%distance=distance
            binatom1%RelCoord=Coord
            binatom1%Species=SpeciesIndex
         END DO !end of atom list
      END DO !end of bin
      CLOSE(UNIT=UnitProcessListRec)
      CALL AddCfgToRecordedCfg(config)
   END DO !end of config
   
   WRITE(6,*) 'Recv>>Total # configs read from previous record:',NumberCfgTypes
   !----------------------END OF CONFIGURATIONS----------------------   

   !----------------------READ ATOM PROCESSES----------------------
   DO indx=1,NumberProcesses1
      ChO=INT2CHAR(indx)
      FileName='./StoreTADKMC/SystemProcesses/P'//TRIM(ChO)
      OPEN(UNIT=UnitProcessListRec,FILE=TRIM(FileName),IOSTAT=errorstatus)
      
      IF (errorstatus/=0) CYCLE !file not found -- process got corrupted, lets proceed further
      
      !IF (.NOT. ASSOCIATED(KMC_PL)) THEN
      !   ALLOCATE(KMC_PL)
      !   NULLIFY(KMC_PL%NextNeigh)
      !   NULLIFY(KMC_PL%PrevNeigh)
      !ELSE !at least one element exists
      !   ALLOCATE(KMC_PL%PrevNeigh)
      !   KMC_PL%PrevNeigh%NextNeigh=>KMC_PL
      !   KMC_PL=>KMC_PL%PrevNeigh
      !   NULLIFY(KMC_PL%PrevNeigh)
      !END IF
      
      ALLOCATE(prc); NAllottedKMCProcess=NAllottedKMCProcess+1

      !Add process information
      !prc=>KMC_PL !allocate a process
      prc%RecordNumber=indx
      prc%Index=indx

      READ(UnitProcessListRec,*) prc%NumberProcessAtoms ![P1]
      READ(UnitProcessListRec,*) ConfigurationIndex ![P1a]
      prc%CfgType=>SearchCfg(ConfigurationIndex)
      IF (.NOT. ASSOCIATED(prc%CfgType)) THEN
         WRITE(6,*) "Err>> Unable to find configuration while recovering process"
         WRITE(6,*) "Process index:",indx
         WRITE(6,*) "Configuration index being searched:",ConfigurationIndex
         STOP
      END IF
      READ(UNIT=UnitProcessListRec,FMT='(2I18)',IOSTAT=errorstatus) &
         prc%FiringCount,prc%MaxFiringAllowed ![P2]
      IF (errorstatus/=0) THEN
         WRITE(6,*) "Error in reading process counts, setting these to default values"
         prc%FiringCount=0
         prc%MaxFiringAllowed=1
      END IF
      IF (ResetNoFirings) prc%FiringCount=0
      IF (FillUpFirings) prc%FiringCount=prc%MaxFiringAllowed

      IF (prc%FiringCount>prc%MaxFiringAllowed) THEN
         !CALL AddToTADToDoList(ptr3)
      END IF

      READ(UnitProcessListRec,*) prc%Rate,prc%Frequency,prc%FwdBarrier ![P3]

      prc%Rate=prc%Frequency*EXP(-prc%FwdBarrier/kBT)
      
      NULLIFY(prcatm)
      DO indx1=1,prc%NumberProcessAtoms
         IF (.NOT. ASSOCIATED(prc%ReactantAL)) THEN
            ALLOCATE(prc%ReactantAL); NAllottedKMCProcessAtom=NAllottedKMCProcessAtom+1
            prcatm=>prc%ReactantAL
            NULLIFY(prcatm%NextNeigh)
            NULLIFY(prcatm%PrevNeigh)
         ELSE
            ALLOCATE(prcatm%NextNeigh); NAllottedKMCProcessAtom=NAllottedKMCProcessAtom+1
            prcatm%NextNeigh%PrevNeigh=>prcatm
            prcatm=>prcatm%NextNeigh
            NULLIFY(prcatm%NextNeigh)
         END IF
         READ(UnitProcessListRec,*) prcatm%RelInitialPos ![P4]
         READ(UnitProcessListRec,*) prcatm%RelSaddlePos ![P4a]
         READ(UnitProcessListRec,*) prcatm%RelFinalPos ![P5]
         READ(UnitProcessListRec,*) prcatm%CfgIndex,prcatm%Species ![P6]
      END DO !end of participants
      CLOSE(UNIT=UnitProcessListRec)
      CALL AddProcessToPL(prc)
      
      prc%FirstAdded=-1 !process added at KMCStep=-1 -- this undoes the FirstAdded being set to 0 by AddProcessToPL
   END DO !end of official
   
   WRITE(6,*) 'Recv>>Total # processes read from previous record:',NumberProcesses
   !---------------------END OF PROCESSES----------------------

   !CALL CheckConsistency
   CALL CheckKMC()
END SUBROUTINE RecvrSystemRecord
!--------------------------------
END MODULE KMCRecords
