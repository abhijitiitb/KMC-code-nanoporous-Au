MODULE LEKMCMD
!provides an interface for finding and transfering process information between MD and kmc
   USE VARIABLE_TYPE
   USE KMC_VARIABLE_TYPE
   USE amd_db_manipulate
   !USE KMCProcessSearchParallel
   USE AMDProcess
   IMPLICIT NONE
   
   TYPE(SystemContainer), POINTER, PRIVATE :: AL=>NULL()
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetupAL(NAtoms,NMove,AtomCoord,AtomSpecies,AtomIsMoving,VLType)
   !adds the atom information to AL while employing periodic boundary conditions to fold back the system
   !ReferenceCoord: coordinate around which the domain was created
      IMPLICIT NONE
      INTEGER :: NAtoms,NMove(3)
      INTEGER :: VLType
      REAL(dp), DIMENSION(3*NAtoms) :: AtomCoord
      INTEGER, DIMENSION(NAtoms) :: AtomSpecies
      LOGICAL, DIMENSION(3*NAtoms) :: AtomIsMoving
      CHARACTER(len=100) :: filename
      
      IF (ASSOCIATED(AL)) THEN
         WRITE(6,*) "Err>> AL is associated in SetupAL"
         STOP
      ELSE
         !IF (imode==1) THEN !this is typically used by full update
         !   CALL ReadCoord(AL,filename)
         !   CALL PotentialInitialize(AL)
         !   CALL AddVerletList(AL,ListType=2,NAtoms=AL%NAtoms)
         !ELSE!used when local domains have to be created
            CALL MakeSize(AL,NAtoms)
            AL%NAtoms=NAtoms
            AL%NMove=NMove
            AL%AtomCoord(1:3*NAtoms)=AtomCoord(1:3*NAtoms)
            AL%AtomSpecies(1:NAtoms)=AtomSpecies(1:NAtoms)
            AL%AtomIsMoving(1:3*NAtoms)=AtomIsMoving(1:3*NAtoms)
            AL%BoxSize(1:3)=KMCBoxSize
            AL%BoxSize(4:6)=0.
            CALL PotentialInitialize(AL)
            CALL AddVerletList(AL,ListType=VLType,NAtoms=AL%NAtoms)
         !END IF
      END IF
   END SUBROUTINE SetupAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE FindMDProcessForLEKMC(NAtoms,AtomEnvironmentIndex,Temperature, &
      nsearch,LEKMCProcess,NLEKMCProcess)
   !finds a list of processes for LEKMC given the starting state of the system
   !inputs: AL,AtomEnvironmentIndex
   !output: LEKMCProcess,NLEKMCProcess
   !Note that the AMD times need to be passed
      IMPLICIT NONE
      INTEGER :: NAtoms,NLEKMCProcess
      REAL(dp) :: Temperature
      INTEGER, DIMENSION(NAtoms) :: AtomEnvironmentIndex
      TYPE(KMCProcess), POINTER :: LEKMCProcess
      TYPE(AMDStateContainer), POINTER :: state
      !parameters passed by LEKMC
      INTEGER :: nsearch,DetectTransitionFrequency,nsuccess,StoppingOption
      REAL(dp) :: MDTimePerBCMD,TotalMDTime,CollectedMDTime
      
      IF (ASSOCIATED(LEKMCProcess)) THEN
         WRITE(6,*) "Err>> The list of processes provided to FindMDProcessForLEKMC should not be associated"
         STOP
      END IF
      NULLIFY(LEKMCProcess)
      
      !Step 1. Add the state to AMD state list, recover old AMD confidence measures
      CALL OptimizeAL(AL,3,NAtoms=AL%NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
         ITMAX=500,MaxDisplacement=0.05_dp)
      state=>AddOptimizedALToGlobalAMDStates(AL) !starting state for AMD
      
      !Step 2. Find processes for the state
      DetectTransitionFrequency=100
      CALL ParallelReplicaMDKMC(state=state,NKMCMove=1,Time=KMCTime,nsearch=nsearch, &
         Temperature=Temperature,imode=1,StoreLEKMCProcess=.TRUE., &
         DetectTransitionFrequency=DetectTransitionFrequency)
      !CALL AddProcessToAMDGlobalState(state,nsearch=4000,CollectedMDTime=CollectedMDTime, &
      !   FindRate=.FALSE.,StoreLEKMCProcess=.TRUE.,Temperature=Temperature,RemoveFiles=.FALSE., &
      !   AtomCfgIndex=AtomEnvironmentIndex,NAtoms=NAtoms)
      !CALL PerformMDProcessSearchParallel(stateindex=1,nsearch=nsearch,Temperature=Temperature, &
      !   DetectTransitionFrequency=DetectTransitionFrequency,MDTimePerBCMD=MDTimePerBCMD, &
      !   TotalMDTime=TotalMDTime,nsuccess=nsuccess,StoppingOption=StoppingOption)
      
      !Step 3. Find NEB activation barriers and the moving atoms for each process 
      !and create a process list for the state
      CALL RecoverProcessFromMDProcess(LEKMCProcess,NLEKMCProcess,state)
      
      !Step 4. Deallocate state information and AL
      CALL DeleteAMDGlobalStateAllProcesses(state)
      CALL DeleteGlobalAMDState(state,deletefile=.TRUE.)
      
   END SUBROUTINE FindMDProcessForLEKMC
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RecoverProcessFromMDProcess(LEKMCProcess,NLEKMCProcess,state)
   !recovers all process information from the MD data
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: LEKMCProcess
      TYPE(AMDStateContainer), POINTER :: state
      TYPE(AMDProcessContainer), POINTER :: amdprc
      INTEGER :: NLEKMCProcess
      
      NLEKMCProcess=state%NProcessType
      !LEKMCProcess=>state%LEKMCProcess
      IF (NLEKMCProcess==0) THEN
         WRITE(6,*) "Err>> Number of processes found using MD for LEKMC is zero"
         STOP
      END IF
   END SUBROUTINE RecoverProcessFromMDProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateTADProcessListFull()
   !finds processes for the full system
   !Used when 1) the system is being initialized or 2) the system size is small and there is no point in doing the local version
      IMPLICIT NONE
      INTEGER :: NAtoms,NMove(3),iatom
      REAL(dp), DIMENSION(3*MaxNAtoms) :: AtomCoord
      INTEGER, DIMENSION(MaxNAtoms) :: AtomSpecies,AtomEnvironmentIndex
      LOGICAL, DIMENSION(3*MaxNAtoms) :: AtomIsMoving
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCProcess), POINTER :: LEKMCProcess
      INTEGER :: nsearch,NLEKMCProcess
      
      !Step 1. Convert KMC_AL to MD file for which processes are found
      NAtoms=NAtomKMC
      NMove=NMoveKMC
      atom=>KMC_AL
      iatom=0
      DO WHILE (ASSOCIATED(atom))
         iatom=iatom+1
         AtomCoord(3*iatom-2:3*iatom)=atom%Coord
         AtomSpecies(iatom)=atom%Species
         AtomIsMoving(3*iatom-2:3*iatom)=atom%IsMoving
         AtomEnvironmentIndex(iatom)=atom%Cfg%ConfigurationIndex
         atom=>atom%NextNeigh
      END DO
      CALL SetupAL(NAtoms,NMove,AtomCoord=AtomCoord(1:3*NAtoms), &
         AtomSpecies=AtomSpecies(1:NAtoms),AtomIsMoving=AtomIsMoving(1:3*NAtoms),VLType=2)
      
      !Step 2. Search for processes in parallel
      nsearch=1000
      NULLIFY(LEKMCProcess)
      CALL FindMDProcessForLEKMC(NAtoms=NAtoms,AtomEnvironmentIndex=AtomEnvironmentIndex, &
         Temperature=REAL(Temperature,dp),nsearch=nsearch,LEKMCProcess=LEKMCProcess, &
         NLEKMCProcess=NLEKMCProcess) !searches for the processes and finds an estimate for the error
      
      !Step 3. Obtain the local processes (this involves finding the process atoms and the activation barrier)
      !Environment of atoms is already known, check if such a local process already exists
   !   LocalProcessAbsent= CheckLocalProcessAbsent()
      !If the process does not exist, find the barriers and rates and add it to the process list 
      !IF (LocalProcessAbsent) CALL AMDLocateProcessAtoms()
      
   END SUBROUTINE GenerateTADProcessListFull
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateTADProcessListCell(CellIndex)
   !local version based on domains
      IMPLICIT NONE
      INTEGER :: CellIndex(3)
   END SUBROUTINE GenerateTADProcessListCell
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE LEKMCMD
