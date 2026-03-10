MODULE MDPackage !NVE or NVT MD calculations
   USE VARIABLE_TYPE
   USE AMD_VARIABLE_TYPE
   USE PotentialPackage
   USE Utilities
   USE db_manipulate
   USE ran_state
   USE Ecuyer_random
   USE OptimizationAL
   USE io
   USE NeighborList
   USE NEB_package
   USE mpimanager
   USE mpi
   USE nmode
   USE Crystal
   USE PatternRecognition
   IMPLICIT NONE
   INTEGER, PRIVATE :: ObtainTSTRate,ObtainTSTRate_ntransitionstates
   REAL(dp), PRIVATE :: ObtainTSTRate_prefactor,ObtainTSTRate_TSEnergy
   LOGICAL, PRIVATE :: CheckMDExit
   
   CONTAINS 
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EnableMDHistory(MDSiml,HistorySizeMax,iprint)
   !enable storage of snapshots in the trajectory
   ! typically, the trajectory stored in MDContainer has the variable
   ! HistoryOn=.FALSE.
   !The past history of the trajectory is stored in a list of 
   ! SystemContainer %History%. Each element of the list contains the
   ! time snaps as well as other details that may be required
   !HistorySizeMax is the maximum number of snaps that can be stored
   !HistoryFreq is the frequency of storing snaps, expressed in terms
   ! of number of time steps
   
   !Last checked: July 03, 2009
   !Last modified: 
   !Last modified by: AC
   
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      INTEGER, OPTIONAL :: HistorySizeMax,iprint!,HistoryFreq
      
      CALL DeleteMDHistory(MDSiml) !enabling the MDHistory would delete any stored history
      MDSiml%HistoryOn=.TRUE.
      IF (PRESENT(HistorySizeMax)) MDSiml%HistorySizeMax=HistorySizeMax
      !IF (PRESENT(HistoryFreq)) MDSiml%HistoryFreq=HistoryFreq
      
      MDSiml%HistorySize=0
      
      IF (PRESENT(iprint)) THEN
         IF (iprint>1) THEN
            WRITE(6,*) " "
            WRITE(6,*) "Enabling history storage in MD trajectory"
            WRITE(UNIT=6,FMT='("     MD history size:",I4)') HistorySizeMax
            !WRITE(UNIT=6,FMT='("     MD history store frequency:",I4)') HistoryFreq
         END IF
      END IF
      
   END SUBROUTINE EnableMDHistory
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DisableMDHistory(MDSiml)   
   !Last checked: July 03, 2009
   !Last modified: 
   !Last modified by: AC
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      
      MDSiml%HistoryOn=.FALSE.
      CALL DeleteMDHistory(MDSiml)
      
   END SUBROUTINE DisableMDHistory
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SearchProcessesFromState(stateindex,Temperature,DetectTransitionFrequency,isearch, &
      success,MDTime,ObtainRate,TemperatureProgram)
   !searches for processes from a state using parallel processing
   !stateindex specifies the initial state -- this info is read from a file
   !the output is stored as process.<stateindex>.<isearch>.attempt where <isearch> is the search attempt (used to create a file with final state)
   !DetectTransitionFrequency is the frequency with which the process is searched
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(MDContainer), POINTER :: MDSiml=>NULL()
      TYPE(AMDStateContainer), POINTER :: state
      CHARACTER(len=100) :: filename,filename2
      REAL(dp) :: Temperature,Rate,MDTime
      INTEGER :: isearch,nsearches,stateindex,DetectTransitionFrequency
      CHARACTER(len=10) :: integrator
      CHARACTER(len=NINT2CHAR) :: charstateindx,charisearch
      INTEGER :: errorstatus
      LOGICAL :: NotSatisfied,success
      INTEGER :: ObtainRate
      INTEGER, OPTIONAL :: TemperatureProgram
      
      NotSatisfied=.TRUE.
      
      ObtainTSTRate=ObtainRate
      ObtainTSTRate_ntransitionstates=-1 !indicates that the numbers have not been calculated
      ObtainTSTRate_TSEnergy=-1._dp !indicates that the numbers have not been calculated
      ObtainTSTRate_prefactor=-1._dp !indicates that the numbers have not been calculated
      
      Rate=0._dp
      DO WHILE (NotSatisfied)
         charstateindx=INT2CHAR(stateindex)
         charisearch=INT2CHAR(isearch)
         filename="./GlobalAMD/state."//TRIM(charstateindx)//".xyz"
         NULLIFY(AL)
         CALL ReadCoord(AL,filename,iprint=2)
         CALL GetMaxwellBoltzmann(AL=AL,Temperature=Temperature,IsPrint=.FALSE.) !initialize velocities
         !CALL Analyze(AL)
         CALL PotentialInitialize(AL,IsPrint=.FALSE.)
         !CALL WritePotential(AL)
         CALL AddVerletList(AL,NAtoms=AL%NAtoms,ListType=1) !SPECIFY THE LIST TYPE
         
         integrator="langevin  "
         CALL AddALToMD(AL,MDSiml,integrator=integrator,coeff=1.e12_dp, &
            temperature=Temperature,TemperatureRamp=0._dp,dt=4.e-15_dp,IsPrint=.FALSE.)
         
         filename="./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".xyz"
         !CALL MD(MDSiml,DetectTransitionFrequency=DetectTransitionFrequency, &
            !WriteXYZTrajectoryFile=filename,errorstatus=errorstatus) !print MD trajectory
         !filename2="tee"//"."//TRIM(charisearch)//".xyz"
         
         CALL MD(MDSiml,DetectTransitionFrequency=DetectTransitionFrequency,errorstatus=errorstatus)! &
         !WriteXYZTrajectoryFile=filename2,errorstatus=errorstatus) !do not print MD trajectory
         
         !Now MDSiml contains the final state -- store the final state
         MDTime=MDSiml%AL%Time !MD time elapsed
         IF (errorstatus==0) THEN
            NotSatisfied=.FALSE.
            
            filename="./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".finalstate.xyz"
            CALL WriteXYZ(MDSiml%AL,filename)
            filename="./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".param"
            OPEN(UNIT=385,FILE=filename)
            WRITE(385,*) MDSiml%AL%Time,"  Time elapsed"
            WRITE(385,*) Temperature,"  Temperature"
            WRITE(385,*) DetectTransitionFrequency,"  DetectTransitionFrequency"
            IF (ObtainTSTRate==1) THEN
               IF (ObtainTSTRate_TSEnergy==-9999._dp) THEN
                  ObtainTSTRate_TSEnergy=0._dp
                  ObtainTSTRate_ntransitionstates=1 !this will be an issue in few cases
                  ObtainTSTRate_prefactor=0._dp
               END IF
               Rate=ObtainTSTRate_prefactor*EXP(-ObtainTSTRate_TSEnergy/kboltzmann/Temperature)
               WRITE(385,*) Rate,ObtainTSTRate_prefactor,ObtainTSTRate_TSEnergy, &
                  ObtainTSTRate_ntransitionstates, &
                  "   Process Rate, prefactor, activation barrier, ntransitionstates"
            ELSE
               WRITE(385,*) Rate,"   Process Rate"
            END IF
            CLOSE(385)
            success=.TRUE.
         ELSEIF (errorstatus==2) THEN !MD was forcibly stopped using file MDTransitionSearchStatus
            NotSatisfied=.FALSE.
            success=.FALSE.
         ELSEIF (errorstatus==3) THEN !MD lost the transition
            NotSatisfied=.TRUE.
            success=.FALSE.
         ELSE !could not find a transition
            NotSatisfied=.FALSE.
            success=.FALSE.
         END IF
         
         CALL Delete(MDSiml) !AL also gets deleted as a result of this
         
      END DO
      
   END SUBROUTINE SearchProcessesFromState
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DepositSurface(AL,Temperature,AtomicSymbol,AtomCharge,WriteChOSTrajectoryInfo, &
      WriteXYZTrajectoryFile,DepositCoord)
      !deposit 1 atom of AtomicSymbol at temperature
      !deposition on surface along negative z direction
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1
      TYPE(MDContainer), POINTER :: MDSiml
      REAL(dp) :: Temperature,coordination,const,BoundingBox(6)
      REAL(dp) :: AtmForce,zcurrent,vmag,h1
      REAL(dp), DIMENSION(3), OPTIONAL :: DepositCoord
      REAL(sp) :: h
      REAL(sp), DIMENSION(:), ALLOCATABLE :: v
      REAL(dp), OPTIONAL :: AtomCharge
      INTEGER :: NAtomsDeposit,i,r1,iteration,NAtoms,ListType,errorstatus,iatom
      CHARACTER(len=3), OPTIONAL :: AtomicSymbol !should this be optional?
      CHARACTER(len=10) :: integrator
      LOGICAL :: Found,WriteChOSTrajectory,WriteXYZTrajectory
      LOGICAL, OPTIONAL :: WriteChOSTrajectoryInfo
      CHARACTER(len=100), OPTIONAL :: WriteXYZTrajectoryFile
      CHARACTER(len=100) :: filename
      
      NAtoms=AL%NAtoms
      ListType=AL%VL%ListType
      IF (ListType/=1 .AND. ListType/=2) THEN
         WRITE(6,*) "$Err>> VL should be assigned to AL in DepositSurface"
         STOP
      END IF
      
      CALL ResizeToLargerAL(AL,NAtoms+1) !Add atom to AL - copy AL to a MDSiml - do MD
      AL%NMove=AL%NMove+1
      
      i=SearchSymbol(AL,AtomicSymbol)
      IF (i>0) THEN
         r1=AL%NAtoms !new atom added
         AL%AtomSpecies(r1)=i
         AL%AtomMass(3*r1-2:3*r1)=SpeciesList%AtomicMass(SpeciesDirectory_global(i))
      ELSE
         WRITE(*,*) "$Err>> Attempted to deposit unknown species ",AtomicSymbol
         STOP
      END IF
      
      IF (PRESENT(AtomCharge)) AL%AtomCharge(r1)=AtomCharge !desired charge
      AL%AtomIsMoving(3*r1-2:3*r1)=.TRUE.
      
      !Find location to deposit the atom from
      BoundingBox=GetBoundingBox(AL)
      IF (BoundingBox(3)-BoundingBox(6)>0.8_dp*(AL%BoxSize(3)-AL%BoxSize(6))) THEN
         WRITE(*,*) "$Err>> It appears no free surface is present"
         STOP
      END IF
      
      IF (PRESENT(DepositCoord)) THEN
         AL%AtomCoord(3*r1-2:3*r1-1)=DepositCoord(1:2)
      ELSE
         h1=taus88() !CALL ran0_s(h)
         AL%AtomCoord(3*r1-2)=h1*(AL%BoxSize(1)-AL%BoxSize(4))
         h1=taus88() !CALL ran0_s(h)
         AL%AtomCoord(3*r1-1)=h1*(AL%BoxSize(2)-AL%BoxSize(5))
      END IF
      AL%AtomCoord(3*r1)=0.1_dp*AL%BoxSize(3)+0.9_dp*BoundingBox(3)
      
      CALL Delete(AL%VL)
      CALL Delete(AL%LL)
      
      filename="deposit.tmp.xyz"
      CALL WriteXYZ(AL,filename)
      !CALL Delete(AL)
      !NULLIFY(AL)
      
      !CALL ReadCoord(AL,filename)
      !CALL PotentialInitialize(AL)
      !CALL GetMaxwellBoltzmann(AL,Temperature)
      CALL AddVerletList(AL,ListType=1,NAtoms=AL%NAtoms)
      
      ALLOCATE(v(3*r1))
      CALL gasdev_v(v) !this is scALed momentum
      const=SQRT(kboltzmann*Temperature)
      !standard deviation without mass
      WHERE (AL%AtomIsMoving) &
         AL%AtomVelocity=const*REAL(v,dp)/SQRT(AL%AtomMass(r1)) !finAL velocity
      vmag=SQRT(SUM(AL%AtomVelocity(3*r1-2:3*r1)*AL%AtomVelocity(3*r1-2:3*r1)))
      DEALLOCATE(v)
      
      CALL ran0_s(h)
      IF (h<0.5) THEN
         AL%AtomVelocity(3*r1-2)=-vmag*0.1
      ELSE
         AL%AtomVelocity(3*r1-2)=vmag*0.1
      END IF
      CALL ran0_s(h)
      IF (h<0.5) THEN
         AL%AtomVelocity(3*r1-1)=-vmag*0.1
      ELSE
         AL%AtomVelocity(3*r1-1)=vmag*0.1
      END IF
      AL%AtomVelocity(3*r1)=-vmag*0.98995
      !CALL Delete(AL%VL)
      !CALL Delete(AL%LL)
      !CALL AddVerletList(AL,ListType=1,NAtoms=NAtoms+1)
      
      CALL AnalyzeAL(AL)
      IF (PRESENT(WriteXYZTrajectoryFile)) THEN
         WriteXYZTrajectory=.TRUE.
         OPEN(UNIT=303,FILE=TRIM(WriteXYZTrajectoryFile))
      ELSE
         WriteXYZTrajectory=.FALSE.
      END IF
      IF (PRESENT(WriteChOSTrajectoryInfo)) THEN
         WriteChOSTrajectory=WriteChOSTrajectoryInfo
      ELSE
         WriteChOSTrajectory=.FALSE.
      END IF

      NULLIFY(MDSiml) 
      !integrator="verlet" 
      integrator="vrscale" !velocity-rescaling
      CALL AddALToMD(AL,MDSiml,dt=4.e-15_dp,integrator=integrator,Temperature=temperature) !we need to check each time because atom is moving down
      IF (WriteChOSTrajectory .OR. WriteXYZTrajectory) &
        CALL EnableMDHistory(MDSiml,HistorySizeMax=200,iprint=1)
      
      Found=.FALSE.
      iteration=0
      !open(Unit=305,FILE="agcc.txt")
      DO WHILE (.NOT. Found)
         iteration=iteration+1
         CALL MD(MDSiml,100,IsNeighborListInitialized=.TRUE.)
         IF (WriteXYZTrajectory) CALL WriteXYZ(AL,fileopen=.FALSE.,fileclose=.FALSE.,iunit=303,iprint=1)
         IF (WriteChOSTrajectory) CALL AddALToMDHistory(MDSiml)
         Found=SUM(ABS(AL%AtomForce(3*r1-2:3*r1)))>0.0_dp
         write(*,*) "DepositAtom:",AL%AtomCoord(3*r1-2:3*r1),AL%AtomVelocity(3*r1)
         !AL%AtomVelocity(3*r1)=-ABS(AL%AtomVelocity(3*r1))
      END DO

      Found=.FALSE.
      zcurrent=AL%AtomCoord(3*r1)
      DO WHILE (.NOT. Found)
         CALL MD(MDSiml,100,IsNeighborListInitialized=.TRUE.)
         IF (WriteXYZTrajectory) CALL WriteXYZ(AL,fileopen=.FALSE.,fileclose=.FALSE.,iunit=303,iprint=1)
         IF (WriteChOSTrajectory) CALL AddALToMDHistory(MDSiml)
         Found= ABS(AL%AtomCoord(3*r1)-zcurrent)<0.3_dp
         zcurrent=AL%AtomCoord(3*r1)
         !write(305,*) AL%AtomCoord(3*r1),AL%AtomVelocity(3*r1)
         !AL%AtomVelocity(3*r1)=-ABS(AL%AtomVelocity(3*r1))
      END DO
      close(305)
      WRITE(6,*) "Completed deposition"
      
      !Minimization
      CALL OptimizeAL(AL,3,NAtoms=AL%NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
         ITMAX=500,MaxDisplacement=0.05_dp)

      !Copy MDSiml%AL to AL
      CALL Delete(MDSiml,AL)
      IF (WriteXYZTrajectory) THEN
         CLOSE(303)
         CALL SYSTEM("cat "//TRIM(WriteXYZTrajectoryFile)//" "//TRIM(filename)//" > "//"tmp.00001")
         CALL SYSTEM("mv tmp.00001 "//TRIM(WriteXYZTrajectoryFile))
         WRITE(6,*) "Deposition trajectory has been written to file ",TRIM(WriteXYZTrajectoryFile)
      END IF
      
      IF (WriteChOSTrajectory) THEN
         WRITE(6,*) "Deposition trajectory has been added to ChOS "
      END IF
      
      NULLIFY(MDSiml)
   END SUBROUTINE DepositSurface
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MDTransitionSearch(MDSiml,ALInitial,ALFinal,AtomicProcess,NAtoms)
   !searches the time when the transition took place in an MD calculation, and find the
   !final state, activation barrier, prefactor
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      TYPE(MDProcessContainer), POINTER :: AtomicProcess
      TYPE(SystemContainer), POINTER :: ALtrack,ALInitial,ALFinal,AL_TS
      INTEGER :: imagebegin,imageend,image,imageno,NAtoms,ntransitionstates
      LOGICAL :: NotFound,matchfinal,matchinitial
      INTEGER :: errorstatus
      REAL(dp) :: TSEnergy

      !Find the correct time block where the process occured by doing a binary search
      !AL%Time will be alloted the "correct" time
      imagebegin=1
      imageend=MDSiml%HistorySize
      image=(imagebegin+imageend)/2
      image=MAX(2,image) !mid-point, e.g., 23/2=11
      NotFound=.TRUE.
      !find the image which matches the final state and the previous image matches initial state
      !the number of images stored in history exceeds the number of images stored in one search cycle in MD
      !consequently, the first image is going to be that of the initial state and the final image is the 
      !final state
      DO WHILE (NotFound)
         ALtrack=>MDSiml%History
         DO imageno=2,image
            ALtrack=>ALtrack%NextNeigh
         END DO
         !optimize followed by comparison
         CALL OptimizeAL(ALtrack,3,NAtoms=NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
            ITMAX=500,MaxDisplacement=0.05_dp)
         CALL OptimizeAL(ALtrack%PrevNeigh,3,NAtoms=NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
            ITMAX=500,MaxDisplacement=0.05_dp)
         matchfinal=CompareALToAL(ALtrack,ALFinal,xtol=1.0_dp)
         IF (matchfinal .AND. CompareALToAL(ALtrack%PrevNeigh,ALInitial,xtol=1._dp)) THEN !transition found
            NotFound=.FALSE.
            EXIT
         END IF
         matchinitial=CompareALToAL(ALtrack,ALInitial,xtol=1._dp)
         IF (.NOT. matchfinal .AND. .NOT. matchinitial) THEN !something is wrong
            WRITE(6,*) "Intermediate state encountered"
            NULLIFY(AtomicProcess)
            RETURN
         END IF 
         IF (matchfinal) imageend=MAX(image,imagebegin+1) !move towards first image
         IF (matchinitial) imagebegin=MIN(image,imageend-1)
         image=MAX(2,(imagebegin+imageend)/2)
      END DO
      ALFinal%Time=(ALtrack%Time+ALtrack%PrevNeigh%Time)/2._dp

      !Find the activation barrier using NEB 
        !For global state we do not need the local displacments of atoms

      CALL NEB(ALInitial,ALFinal,nimages=9,imode=1101,omode=13,interpmode=1, &
        springconstant=1._dp,ITMAX=100,TSEnergy=TSEnergy,ntransitionstates=ntransitionstates, &
        TS=AL_TS,errorstatus=errorstatus)
      CALL Delete(AL_TS)
      IF (ntransitionstates>1) THEN   
         CALL NEB(ALInitial,ALFinal,nimages=15,imode=1101,omode=13,interpmode=1, &
           springconstant=1._dp,ITMAX=100,TSEnergy=TSEnergy,ntransitionstates=ntransitionstates, &
           TS=AL_TS,errorstatus=errorstatus)
         WRITE(6,*) "Multiple transition states found in NEB"
         NULLIFY(AtomicProcess)
         CALL Delete(AL_TS)
         RETURN
      END IF
      !Find the prefactor using Vineyard
      ALLOCATE(AtomicProcess) !this will be added to the system at the end 
      AtomicProcess%prefactor=1.e13
      AtomicProcess%Eact=TSEnergy
   END SUBROUTINE MDTransitionSearch
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MD( &
      MDSiml, & !data structure containing the system
      nmd, &    !number of steps for which MD needs to be performed
      IsNeighborListInitialized, & !is the VL initialized for the system
      !AtomicProcess, &  !information about the atomic process - initially this should be null
      DetectTransitionFrequency, & !frequency of detecting atomic process (used as nmd in such a case)
      DeleteHistory, & !delete the history stored so far
      iprint, & !print level
      WriteXYZTrajectoryFile, & !if provided the sbrtn will open a new file with this name write the trajectory to it -- unit file number 303
      WriteChOSTrajectoryInfo, &  !store snapshots in ChOS
      WriteEnergyTrajectoryFile, &  !print energy information -- unit number 304
      WriteTrajectoryFrequency, & !frequency of printing information - default value is nmd1/10
      WriteNAtoms, & !(optional) number of atoms to be printed in XYZ trajectory file 
      WritePeriodic, & !wrap the coordinates using PBC while printing the trajectory
      errorstatus)
      
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      TYPE(SystemContainer), POINTER :: AL,ALreference,AL2 !,ALOptim
      !TYPE(MDProcessContainer), POINTER, OPTIONAL :: AtomicProcess
      INTEGER, OPTIONAL :: nmd,DetectTransitionFrequency,iprint,WriteTrajectoryFrequency,WriteNAtoms
      INTEGER, OPTIONAL :: errorstatus
      LOGICAL, OPTIONAL :: IsNeighborListInitialized,WriteChOSTrajectoryInfo,DeleteHistory,WritePeriodic
      INTEGER :: nmd1,NAtoms,i,WriteTrajectoryFrequency1,iprint1,WriteNAtoms1
      REAL(dp), DIMENSION(:), POINTER :: xp,xc,xf,m,xorig
      LOGICAL :: DetectTransition1,NoTransitionOccurred,DeleteHistory1,IsVLInitialized
      LOGICAL :: WriteXYZTrajectory,WriteEnergyTrajectory,WriteChOSTrajectory,WritePeriodic1
      REAL(dp) :: cputime1,cputime2,dt,dt2,fv,dTime
      REAL(dp) :: c0,c1,c2,a1,a2,ldt,ldtinv
      REAL(dp) :: TemperatureRamp !used by Temperature scheduler
      CHARACTER(len=100), OPTIONAL :: WriteXYZTrajectoryFile,WriteEnergyTrajectoryFile
      REAL(dp), DIMENSION(:), POINTER :: c,sim,sigr,sigv,vf,vc,ap,ac,crv,DriftDistance
      REAL(sp), DIMENSION(:), POINTER :: g1,g2
      CHARACTER(len=100) :: filename
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomVelocity,AtomForce
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      LOGICAL :: NotThermalized,MDContinue,MatchOccurred
      INTEGER :: rn,nattempts,errorstatusread,msgtype,status=0 !,ierr,status,TAG
      REAL(dp) :: MDStopTime
      REAL(dp) :: dx1
      CHARACTER(len=100) :: filna
      INTEGER :: ntransitionstates
      REAL(dp) :: TSEnergy,prefactor

      CALL CPU_TIME(cputime1)
      INCLUDE "md_setup.f90"
      
      NAtoms=AL%NAtoms
      rn=3*NAtoms
      dt=MDSiml%dt
      IF (dt==0._dp) THEN
         WRITE(*,*) TRIM(TxtHeader)//"$Err>> MD time step has not been setup"
         STOP
      END IF
      fv=SQRT(fev/famu)/fang !conversion factor to get 
      !Angstrom from time*velocity (time=sec, velocity=SQRT(eV/amu))
      ! as well as Velocity from time*accl
      
      NULLIFY(xp)
      NULLIFY(xc)
      NULLIFY(xf)
      NULLIFY(m)
      
      CALL MakeSize(xp,3*NAtoms)
      CALL MakeSize(xc,3*NAtoms)
      CALL MakeSize(xf,3*NAtoms)
      CALL MakeSize(m,3*NAtoms)

      AtomCoord=>AL%AtomCoord
      AtomVelocity=>AL%AtomVelocity
      AtomForce=>AL%AtomForce
      AtomIsMoving=>AL%AtomIsMoving
      
      IF (DeleteHistory1) CALL DeleteMDHistory(MDSiml)
      IF (WriteXYZTrajectory) CALL WriteXYZ(AL,fileopen=.FALSE.,fileclose=.FALSE.,iunit=303,iprint=1)
      IF (WriteChOSTrajectory) CALL AddALToMDHistory(MDSiml,iprint1)
      
      IF (PRESENT(errorstatus)) errorstatus=1 !there is error
      
      IF (MDSiml%Integrator(1:4)=="verl") THEN
         INCLUDE "md_verlet_setup.f90"
      END IF
      IF (MDSiml%Integrator(1:4)=="lang") THEN
         INCLUDE "md_langevin_setup.f90"
      END IF
      
      xc=AtomCoord(1:rn)
      ALLOCATE(DriftDistance(3*NAtoms))
      DriftDistance=AL%AtomDriftVelocity(1:rn)*1.e10*dt !AL%AtomDriftVelocity is in m/s
      
      NULLIFY(ALreference)
      NULLIFY(AL2)
      
      CheckMDExit=.FALSE.
      IF (DetectTransition1) THEN
         
         !CALL DeleteTransition(AtomicProcess) !see md_setup where this has been commented out
         !IF (ASSOCIATED(AtomicProcess)) THEN
         !   WRITE(6,*) "$Err>> AtomicProcess should have been null in subroutine MD"
         !   STOP
         !END IF
         
         CheckMDExit=.TRUE.
         
         CALL Delete(ALreference)
         CALL GeometricOptimization()
         CALL MakeSize(ALreference,NAtoms)
         CALL Copy(AL,ALreference) !ALreference contains the optimized initial state
         CALL MakeSize(AL2,NAtoms)
         CALL Copy(AL,AL2) !ALreference contains the optimized initial state

         !CALL Delete(ALOptim)
         !CALL MakeSize(ALOptim,NAtoms)
         !CALL Copy(AL,ALOptim) !created to analyze the error with optimization

         !thermalize
         !xxxxxxxxxxxxxxxxxxxx
         nmd1=50
         WriteTrajectoryFrequency1=1000 !if WriteTrajectoryFrequency1 is larger than nmd1 then 
         !nothing will be written to file even though the file might be open
         NotThermalized=.TRUE.
         nattempts=0
         
         DO WHILE (NotThermalized)
            IF (ANY(AL%NMove==0)) THEN
               WRITE(6,*) TRIM(TxtHeader)//"NMove is zero"
               STOP
            END IF
            CALL RunMD()
            CALL GeometricOptimization()
            
            IF (CompareALToAL(ALreference,AL,xtol=1.5_dp)) THEN
               WRITE(6,*) TRIM(TxtHeader)//"Thermalized system successfully"
               NotThermalized=.FALSE.
               AtomCoord(1:rn)=xc !current MD position
            ELSE
               WRITE(6,*) TRIM(TxtHeader)//"Attempt to thermalize failed ... attempting again (curr attempt:", &
                  nattempts,")"
               nattempts=nattempts+1
               IF (nattempts>5) THEN
                  filename="AMDProcessSearchAL.xyz"
                  CALL WriteXYZ(AL,filename)
                  filename="AMDProcessSearchALReference.xyz"
                  CALL WriteXYZ(ALreference,filename)
                  STOP
               END IF
               AL%AtomCoord(1:rn)=ALreference%AtomCoord(1:rn)
               xc=AtomCoord(1:rn)
               CALL GetMaxwellBoltzmann(AL,MDSiml%Temperature) !reinitialize velocities
               WRITE(6,*) AL%AtomVelocity(1501:1503)
            END IF
         END DO
         
         !transition search
         !xxxxxxxxxxxxxxxxxxxx
         !nmd1=DetectTransitionFrequency*20 !since we shall store we do not need to optimize frequently
         nmd1=DetectTransitionFrequency
         WriteTrajectoryFrequency1=DetectTransitionFrequency
         
         AtomCoord(1:rn)=xc !have the latest thermalized coordinates

         NoTransitionOccurred=.TRUE.
         AL%Time=0._dp
         DO WHILE (NoTransitionOccurred)
            
            !check if the master is telling the MD to continue or not
            !old way of doing this
            !nattempts=0
            !DO WHILE (nattempts<5)
            !   nattempts=nattempts+1
            !   OPEN(UNIT=318,FILE="MDTransitionSearchStatus")
            !   READ(UNIT=318,FMT=*,IOSTAT=errorstatusread) MDStopTime
            !   CLOSE(318)
            !   IF (errorstatusread==0) EXIT
            !END DO
            !new way
            CALL MDExitQuery()
            IF (.NOT. MDContinue) THEN !time to stop
               WRITE(6,*) TRIM(TxtHeader)//"Current MD time:",AL%Time
               WRITE(6,*) TRIM(TxtHeader)//"Found stop time:",MDStopTime
               IF (PRESENT(errorstatus)) errorstatus=2
               CALL Delete(ALreference)
               RETURN
            END IF
            
            IF (ANY(AL%NMove==0)) THEN
               WRITE(6,*) TRIM(TxtHeader)//"NMove is zero"
               STOP
            END IF
            
            
            CALL RunMD()
            
            !continue with transition search
            !ALOptim%AtomCoord(1:xc)=AtomCoord(1:xc) !created to analyze the optimization related bug

            !CALL OptimizeAL(AL,omode,NAtoms=NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
            !   ITMAX=500,MaxDisplacement=0.05_dp)
            !filename="./GlobalAMD/Optim.xyz"
            
            IF (.NOT. MDContinue) THEN !time to stop
               WRITE(6,*) TRIM(TxtHeader)//"Current MD time:",AL%Time
               WRITE(6,*) TRIM(TxtHeader)//"Found stop time:",MDStopTime
               IF (PRESENT(errorstatus)) errorstatus=2
               CALL Delete(ALreference)
               RETURN
            END IF
            
            CALL GeometricOptimization()
            NoTransitionOccurred=CompareALToAL(ALreference,AL,xtol=1.5_dp,dx1=dx1)
            
            !IF (dx1>0.4_dp .AND. dx1<1.0_dp) THEN
          
            !   write(6,*)"dx =",dx1
            !   AL%AtomCoord=xc
            !   filna=TRIM(TxtHeader)//"ALOptimiz.xyz"
            !   CALL WriteXYZ(AL,filna)
            !   filna=TRIM(TxtHeader)//"ALRefere.xyz"
            !   CALL WriteXYZ(ALreference,filna)
         
            !   CALL NEB(ALreference,AL,nimages=15,imode=1101,omode=13,interpmode=1, &
            !      springconstant=1._dp,ITMAX=50,ntransitionstates=ntransitionstates, &
            !      TSEnergy=TSEnergy,TS=AL2,errorstatus=errorstatus)
            !   WRITE(UNIT=6,FMT='("Transition state energy: ",F8.3," eV")') TSEnergy
            !   CALL FLUSH(6)
            !   STOP
          
            !END IF
            
            IF (.NOT. NoTransitionOccurred) THEN !transition occurred, another attempt at the optimization
               CALL GeometricOptimization()
               NoTransitionOccurred=CompareALToAL(ALreference,AL,xtol=1.5_dp)
            END IF
            
            IF (NoTransitionOccurred) THEN !transition has not occurred
               AL%AtomCoord=xc !important that xc stores current state of system
            ELSE
               !CALL AddProcessToState(AtomicProcess,AL) !this is the final state
               WRITE(UNIT=6,FMT='("Net MD time elapsed:",ES15.3)') AL%Time
            END IF
         END DO
         
         !Transition has occurred, now lets do some extra MD
         nmd1=100 !extra 1000 steps of MD
         dTime=AL%Time
         CALL Delete(ALreference)
         CALL MakeSize(ALreference,NAtoms)
         CALL Copy(AL,ALreference) !ALreference contains the optimized final state
         
         DeleteHistory1=.FALSE. !to store the extra part
         CALL RunMD()
         AL%Time=dTime
         CALL GeometricOptimization()
         MatchOccurred= CompareALToAL(ALreference,AL,xtol=1.5_dp)
         
         IF (.NOT. MatchOccurred) THEN !transition occurred, another attempt at the optimization
            CALL GeometricOptimization()
            MatchOccurred=CompareALToAL(ALreference,AL,xtol=1.5_dp)
         END IF
            
         IF (MatchOccurred) THEN
            WRITE(6,*) TRIM(TxtHeader)//"Successfully found transition"
            IF (PRESENT(errorstatus)) errorstatus=0
            IF (ObtainTSTRate==1) THEN
               CALL HTST(AL2,ALreference,ntransitionstates,TSEnergy,prefactor)
               ObtainTSTRate_ntransitionstates=ntransitionstates
               ObtainTSTRate_TSEnergy=TSEnergy
               ObtainTSTRate_prefactor=prefactor
            END IF
         ELSE
            WRITE(6,*) TRIM(TxtHeader)//"$Err>> Lost the transition that was discovered"
            IF (PRESENT(errorstatus)) THEN
               errorstatus=3 !failed to find a proper transition
            ELSE
               STOP
            END IF
         END IF
         !CALL MDTransitionSearch(MDSiml,ALreference,AL,AtomicProcess,NAtoms)
         
         CALL Delete(ALreference)
         CALL Delete(AL2)

      ELSE
         CALL RunMD()
         IF (PRESENT(errorstatus)) errorstatus=0 !there is error
      END IF
      
      DEALLOCATE(xp,xc,xf,m,DriftDistance)
      !IF (MDSiml%Integrator(1:4)=="verl") THEN
      !   DEALLOCATE(xp)
      !END IF
      IF (MDSiml%Integrator(1:4)=="lang") THEN
         DEALLOCATE(ap,ac,vf,vc,sim,sigr,sigv,g1,g2,crv,c)
      END IF

      IF (DetectTransition1) CALL Delete(ALreference)
      IF (WriteXYZTrajectory) CLOSE(303)
      IF (WriteEnergyTrajectory) CLOSE(304)
      CALL CPU_TIME(cputime2)
      IF (iprint1>0) WRITE(UNIT=6,FMT='(">> MD CPU time per atom per step",ES15.5)') &
         (cputime2-cputime1)/REAL(AL%NAtoms,dp)/REAL(nmd1,dp)
   CONTAINS
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GeometricOptimization()
   !performs energy optimization with AL
      IMPLICIT NONE
      INTEGER :: iattempt,errorstatus2
      
      !WRITE(6,*) "Optimization is being performed"
      DO iattempt=1,5
         errorstatus2=0
         !steepest descent
         !CALL OptimizeAL(AL,1,NAtoms=AL%NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
         !   ITMAX=4,MaxDisplacement=0.5_dp)
         !l-bfgs
         CALL OptimizeAL(AL,3,NAtoms=AL%NAtoms,ftol=1.e-5_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
            ITMAX=50,MaxDisplacement=0.05_dp,errorstatus=errorstatus2)
         IF (errorstatus2==0) EXIT
         WRITE(6,*) TRIM(TxtHeader)//"Attempt to optimize structure failed: attempt #",iattempt; CALL FLUSH(6)
      END DO
      IF (errorstatus2/=0) THEN
         WRITE(6,*) "Err>> Optimization could not be completed"
         STOP
      END IF
      !WRITE(6,*) "Optimization is completed"
   END SUBROUTINE GeometricOptimization
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RunMD()
      IMPLICIT NONE
      
      IF (MDSiml%Integrator(1:4)=="verl") CALL MDVerlet()
      IF (MDSiml%Integrator(1:4)=="velo") CALL MDVelocityVerlet()
      IF (MDSiml%Integrator(1:4)=="leap") CALL MDLeapFrog()
      IF (MDSiml%Integrator(1:4)=="lang") CALL MDLangevin()
      IF (MDSiml%Integrator(1:4)=="vrsc") CALL MDVelocityRescale()
   END SUBROUTINE RunMD
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MDVerlet(velocityrescalefreq)
      IMPLICIT NONE
      INTEGER :: i,j,iblock,errorstatus,ip
      REAL(dp) :: invtime,Time,Temperature
      INTEGER, OPTIONAL :: velocityrescalefreq
      INTEGER :: velocityrescalefreq1
      CHARACTER(len=100) :: filename
      LOGICAL :: velocityrescale
      
      INTEGER :: atmaxdisp
      REAL(dp) :: maxdisp
      
      IsVLInitialized=.NOT. AL%VL%UnInitialized
      
      velocityrescale=.FALSE.
      velocityrescalefreq1=10000
      IF (PRESENT(velocityrescalefreq)) THEN
         velocityrescale=.TRUE.
         velocityrescalefreq1=velocityrescalefreq
      END IF
      
      dt2=dt*dt
      m=AL%AtomMass/dt2/fv/fv !now force/m is Ang units
      
      DO i=1,rn
         IF (AtomIsMoving(i)) THEN
            xp(i)=xc(i)-AtomVelocity(i)*dt*fv
         ELSE
            xp(i)=xc(i)
            AtomVelocity(i)=0._dp
         END IF
      END DO
      !WHERE (AtomIsMoving(1:rn))
      !   xp=xc-AtomVelocity(1:rn)*dt*fv
      !ELSEWHERE
      !   xp=xc
      !   AtomVelocity(1:rn)=0._dp
      !END WHERE
      
      !WRITE(*,*) "Max vel contri",MAXVAL(ABS(AL%AtomVelocity*dt*fv))
      !CALL AddVerletList(AL,ListType=AL%VL%ListType, &
      !   IsUninitialized=.NOT. IsVLInitialized,NAtoms=NAtoms)
      CALL GetForcesVL(AL,NAtoms,errorstatus)
      IF (errorstatus/=0) STOP
      invtime=0.5_dp/dt/fv !from ang/s to SQRT(eV/amu)
      
      Time=AL%Time
      
      iblock=0
      xf=xc
      DO i=1,nmd !md loop
         iblock=iblock+1
         CALL AddVerletListQuick(AL,NAtoms)
         CALL GetForcesVL(AL,NAtoms,errorstatus)
         IF (errorstatus/=0) STOP
         
         !write(*,*) "Difference:",AtomCoord(1:6)-xp(1:6)
         DO ip=1,rn
            IF (AtomIsMoving(ip)) THEN
               xf(ip)=xc(ip)+xc(ip)-xp(ip)+AtomForce(ip)/m(ip)
            ELSE
               xf(ip)=xc(ip)
            END IF
         END DO
         !WHERE (AtomIsMoving(1:rn))
         !   xf=xc+xc-xp+AtomForce(1:rn)/m
         !   !xf=AtomCoord+AtomCoord-xp+AtomForce/m !role of PBC should be accounted for -- !in this form AtomCoord could be significantly different from xp because PBC was !used for AtomCoord
         !ELSEWHERE
         !   xf=xc
         !   !AtomForce=0._dp
         !END WHERE
         !WRITE(*,*) MAXVAL(ABS(xf-xc)),MAXVAL(ABS(AtomForce))
         !WRITE(UNIT=6,FMT='(" >> Is VL correct? ",L)') CheckVL(AL)
         xp=xc
         xc=xf
         AtomCoord(1:rn)=xc
         
         Time=Time+dt
         IF (iblock==WriteTrajectoryFrequency1) THEN
            iblock=0
            AL%Time=Time
            DO ip=1,rn
               IF (AtomIsMoving(ip)) THEN
                  AtomVelocity(ip)=2._dp*invtime*(xc(ip)-xp(ip))
               ELSE
                  xf(ip)=xc(ip)
                  AtomVelocity(ip)=0._dp
               END IF
            END DO
            !WHERE (AtomIsMoving(1:rn))
            !   !xf=xc+xc-xp+AtomForce(1:rn)/m -- why do we have this?? is this because of 
            !   AtomVelocity(1:rn)=invtime*(xf-xp)
            !ELSEWHERE
            !   xf=xc
            !   AtomVelocity=0._dp
            !END WHERE
            CALL PBC(AL)
            IF (WriteXYZTrajectory) CALL WriteXYZ(AL,fileopen=.FALSE.,fileclose=.FALSE.,iunit=303,iprint=1)
            !IF (WriteXYZTrajectory) CALL WriteXYZCSP(AL,fileopen=.FALSE.,fileclose=.FALSE.,iunit=303,iprint=1, &
            !   min_csp=0.5,nmax=8,rcut=2.6) !use CentroSymmetricParameter for printing with csp>0.5
               
            IF (WriteChOSTrajectory) CALL AddALToMDHistory(MDSiml,iprint=0)
            IF (WriteEnergyTrajectory) THEN
               CALL GetKineticEnergy(AL)
               WRITE(UNIT=304,FMT='(I8,6ES16.6)') i,AL%Time,MDSiml%Temperature, &
                 AL%PotentialEnergy,AL%KineticEnergy,AL%KineticEnergy+AL%PotentialEnergy,GetTemperature(AL)
               CALL FLUSH(304)
            END IF
         END IF
         
      END DO
      
      AL%Time=Time
      !another iteration to get the velocity
      DO ip=1,rn
         IF (AtomIsMoving(ip)) THEN
            xf(ip)=xc(ip)+xc(ip)-xp(ip)+AtomForce(ip)/m(ip)
            AtomVelocity(ip)=invtime*(xf(ip)-xp(ip))
         ELSE
            xf(ip)=xc(ip)
            AtomVelocity(ip)=0._dp
         END IF
      END DO
      !WHERE (AtomIsMoving)
      !   xf=xc+xc-xp+AtomForce(1:rn)/m
      !   AtomVelocity(1:rn)=invtime*(xf-xp)
      !ELSEWHERE
      !   xf=xc
      !   AtomVelocity(1:rn)=0._dp
      !END WHERE
      
      IF (velocityrescale .AND. MOD(i,velocityrescalefreq1)==0) THEN
         Temperature=GetTemperature(AL)
         DO ip=1,3*AL%NAtoms
            AtomVelocity(ip)=AtomVelocity(ip)*SQRT(MDSiml%Temperature/Temperature)
         END DO
      END IF
         
   END SUBROUTINE MDVerlet
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MDLeapFrog()
      IMPLICIT NONE
      
   END SUBROUTINE MDLeapFrog
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MDVelocityVerlet()
      IMPLICIT NONE
      INTEGER :: i,j,iblock,errorstatus,ip
      REAL(dp) :: invtime,Time
      
      WRITE(6,*) "MDVelocityVerlet has not been checked"
      STOP
      
      IsVLInitialized=.NOT. AL%VL%UnInitialized
      
      dt2=dt*dt
      m=AL%AtomMass/dt2/fv/fv !now force/m is Ang units
      
      CALL GetForcesVL(AL,NAtoms,errorstatus)
      IF (errorstatus/=0) STOP
      invtime=0.5_dp/dt/fv !from ang/s to SQRT(eV/amu)
      
      Time=AL%Time
      
      iblock=0
      xf=xc
      DO i=1,nmd !md loop
         iblock=iblock+1
         
         DO ip=1,rn
            IF (AtomIsMoving(ip)) THEN
               xc(ip)=xc(ip)+vc(ip)*dt*fv+AtomForce(ip)/m(ip)/2._dp !update the positions
               vc(ip)=vc(ip)+AtomForce(ip)/m(ip)/2._dp/fv/dt !partial update to the velocity
            ELSE
               xf(ip)=xc(ip)
            END IF
         END DO
         
         AtomCoord(1:rn)=xc
         CALL AddVerletListQuick(AL,NAtoms)
         CALL GetForcesVL(AL,NAtoms,errorstatus)
         IF (errorstatus/=0) STOP
         
         DO ip=1,rn
            IF (AtomIsMoving(ip)) THEN
               vc(ip)=AtomForce(ip)/m(ip)/2._dp/fv/dt !remaining update to the velocity
            ELSE
               xf(ip)=xc(ip)
            END IF
         END DO
         
         Time=Time+dt
         IF (iblock==WriteTrajectoryFrequency1) THEN
            iblock=0
            AL%Time=Time
            DO ip=1,rn
               IF (AtomIsMoving(ip)) THEN
                  AtomVelocity(ip)=2._dp*invtime*(xc(ip)-xp(ip))
               ELSE
                  xf(ip)=xc(ip)
                  AtomVelocity(ip)=0._dp
               END IF
            END DO
            !WHERE (AtomIsMoving(1:rn))
            !   !xf=xc+xc-xp+AtomForce(1:rn)/m -- why do we have this?? is this because of 
            !   AtomVelocity(1:rn)=invtime*(xf-xp)
            !ELSEWHERE
            !   xf=xc
            !   AtomVelocity=0._dp
            !END WHERE
            CALL PBC(AL)
            IF (WriteXYZTrajectory) CALL WriteXYZ(AL,fileopen=.FALSE.,fileclose=.FALSE.,iunit=303,iprint=1)
            !IF (WriteXYZTrajectory) CALL WriteXYZCSP(AL,fileopen=.FALSE.,fileclose=.FALSE.,iunit=303,iprint=1, &
            !   min_csp=0.5,nmax=8,rcut=2.6) !use CentroSymmetricParameter for printing with csp>0.5
               
            IF (WriteChOSTrajectory) CALL AddALToMDHistory(MDSiml,iprint=0)
            IF (WriteEnergyTrajectory) THEN
               CALL GetKineticEnergy(AL)
               WRITE(UNIT=304,FMT='(I8,6ES16.6)') i,AL%Time,MDSiml%Temperature, &
                 AL%PotentialEnergy,AL%KineticEnergy,AL%KineticEnergy+AL%PotentialEnergy,GetTemperature(AL)
               CALL FLUSH(304)
            END IF
         END IF
         
      END DO
      
      AL%Time=Time
      !another iteration to get the velocity
      DO ip=1,rn
         IF (AtomIsMoving(ip)) THEN
            xf(ip)=xc(ip)+xc(ip)-xp(ip)+AtomForce(ip)/m(ip)
            AtomVelocity(ip)=invtime*(xf(ip)-xp(ip))
         ELSE
            xf(ip)=xc(ip)
            AtomVelocity(ip)=0._dp
         END IF
      END DO
      !WHERE (AtomIsMoving)
      !   xf=xc+xc-xp+AtomForce(1:rn)/m
      !   AtomVelocity(1:rn)=invtime*(xf-xp)
      !ELSEWHERE
      !   xf=xc
      !   AtomVelocity(1:rn)=0._dp
      !END WHERE
      
   END SUBROUTINE MDVelocityVerlet
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MDLangevin()
      !Langevin dynamics
      
      IMPLICIT NONE
      REAL(dp) :: Time,c1u,c2u,TemperatureInitially,TemperatureRatioSqrt
      INTEGER :: i,j,iblock,errorstatus,CheckMDExitFreq
      CHARACTER(len=100) :: filename

      vc=AtomVelocity(1:rn) !based on current temperature (note temperature can change)
      
      IsVLInitialized= .NOT. AL%VL%UnInitialized
      
      CALL AddVerletList(AL,IsUninitialized=.NOT. IsVLInitialized,NAtoms=NAtoms)
      CALL GetForcesVL(AL,NAtoms,errorstatus)
      IF (errorstatus/=0) STOP
      ac=AtomForce(1:rn)/m

      c1u=c1*fv*dt !to reduce # multiplications
      c2u=c2*fv*dt
      
      Time=MDSiml%AL%Time
      
      !temperature scheduling
      TemperatureInitially=MDSiml%Temperature !initial temperature
      
      CheckMDExitFreq=100000000
      IF (CheckMDExit) CheckMDExitFreq=100
      
      iblock=0
      DO i=1,nmd1
         iblock=iblock+1
         CALL gasdev_v(g1) !added Dec. 08, 2009
         CALL gasdev_v(g2) !added Dec. 08, 2009
         
         TemperatureRatioSqrt=SQRT(1._dp+TemperatureRamp*REAL(i,dp)/TemperatureInitially)
         
         g2=sigv*(crv*g1+c*g2)*TemperatureRatioSqrt !units of velocity (TemperatureRatioSqrt is present to ramp the temperature)
         g1=sigr*g1*TemperatureRatioSqrt !units of distance
         WHERE (AtomIsMoving)
            !xc=xc + fv*dt*c1*vc + fv*dt*c2*fv*dt*ac + REAL(g1,dp)
            xc=xc + c1u*vc + c2u*fv*dt*ac + g1
            ap=ac
         END WHERE

         xc=xc+DriftDistance
         
         !write(*,*) MAXVAL(AtomCoord(1:3*NAtoms)-xc(1:3*NAtoms))
         !stop
         AtomCoord(1:rn)=xc
         CALL AddVerletListQuick(AL,NAtoms)
         !make sure that the linked list is updated at each iteration
         !this is going to be important for coulombic interactions which requires the linked list
         CALL GetForcesVL(AL,NAtoms,errorstatus)
         IF (errorstatus/=0) STOP
         ac=AtomForce(1:rn)/m
         !WHERE (AL%AtomIsMoving) vc=c0*vc + fv*dt*((c1-c2)*ap + c2*ac) + REAL(g2,dp)
         WHERE (AtomIsMoving) vc=c0*vc + (c1u-c2u)*ap + c2u*ac + g2
      
         Time=Time+dt
         
         IF (iblock==WriteTrajectoryFrequency1) THEN
            iblock=0
            WHERE (AtomIsMoving) AtomVelocity=vc

            !Print the trajectory file
            AtomCoord(1:rn)=xc(1:rn)
            IF (WritePeriodic1) CALL PBC(AL) !wrap around using periodic BC
            !CALL Optimize(AL,3)
            IF (WriteXYZTrajectory) CALL WriteXYZ(AL,fileopen=.FALSE.,fileclose=.FALSE.,iunit=303,iprint=1, &
               NAtoms=WriteNAtoms1,filtermode=1)
            !IF (WriteXYZTrajectory) CALL WriteXYZCSP(AL,fileopen=.FALSE.,fileclose=.FALSE.,iunit=303,iprint=1, &
            !   min_csp=0.5,nmax=8,rcut=2.6) !use CentroSymmetricParameter for printing with csp>0.5

            !Store in ChOS
            !CALL PBC(AL)
            IF (WriteChOSTrajectory) CALL AddALToMDHistory(MDSiml,iprint1)

            IF (WriteEnergyTrajectory) THEN
               CALL GetKineticEnergy(AL)
               WRITE(UNIT=304,FMT='(I8,6ES16.6)') i,AL%Time,MDSiml%Temperature, &
                 AL%PotentialEnergy,AL%KineticEnergy,AL%KineticEnergy+AL%PotentialEnergy,GetTemperature(AL)
            END IF
            
         END IF
         MDSiml%AL%Time=Time
         MDSiml%Temperature=TemperatureInitially+TemperatureRamp*REAL(i,dp) !update the temperature
         !INCLUDE "md_langevin_temperature_sched.f90"
         
         IF (MOD(i,CheckMDExitFreq)==0) THEN
            CALL MDExitQuery()
            IF (.NOT. MDContinue) EXIT
         END IF
            
      END DO
      WHERE (AtomIsMoving) AtomVelocity(1:rn)=vc
      
   END SUBROUTINE MDLangevin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MDAndersen()
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      
      
   END SUBROUTINE MDAndersen
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MDVelocityRescale()
   !uses velocity rescaling for constant temperature
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      
      CALL MDVerlet(velocityrescalefreq=10)
      
   END SUBROUTINE MDVelocityRescale
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MDNoseHoover()
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      
   END SUBROUTINE MDNoseHoover
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MDParrinelloRahmanConstPressure()
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      
   END SUBROUTINE MDParrinelloRahmanConstPressure
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MDParrinelloRahmanConstStress()
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      
   END SUBROUTINE MDParrinelloRahmanConstStress
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   INCLUDE "md.bondconstraints.f90"
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MDExitQuery()
      IMPLICIT NONE
      
      TAG=0
      status=0
      ierr=0
      msgtype=2
      CALL MPI_SEND(irank,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(msgtype,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_RECV(MDStopTime,1,MPI_DOUBLE_PRECISION,0,TAG,MPI_COMM_WORLD,status,ierr)
      
      MDContinue=AL%Time<MDStopTime
   END SUBROUTINE MDExitQuery
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE MD
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE MDPackage
