
      IF (.NOT. ASSOCIATED(MDSiml)) THEN
         WRITE(*,*) "$Err>> MD object has not been setup"
         STOP
      END IF
      
      IF (PRESENT(WriteChOSTrajectoryInfo)) THEN
         WriteChOSTrajectory=WriteChOSTrajectoryInfo
      ELSE
         WriteChOSTrajectory=.FALSE.
      END IF
      
      IF (PRESENT(DetectTransitionFrequency)) THEN
         DetectTransition1=.TRUE.
         nmd1=DetectTransitionFrequency
         CALL EnableMDHistory(MDSiml,HistorySizeMax=30,iprint=1) !chos to be created
         WriteChOSTrajectory=.TRUE.
         IF (PRESENT(nmd)) THEN
            WRITE(6,*) "$Err>> nmd should not have been provided during transition detection"
            STOP
         END IF
      ELSE
         DetectTransition1=.FALSE.
         IF (PRESENT(nmd)) THEN
            nmd1=nmd
         ELSE
            WRITE(6,*) "$Err>> Specify number of MD steps"
            STOP
         END IF
      END IF
      
      IF (WriteChOSTrajectory .AND. .NOT. MDSiml%HistoryOn) THEN
         WRITE(6,*) "ChOS trajectory written only when history is switched on"
         STOP
      END IF
      
      IF (PRESENT(DeleteHistory)) THEN
         DeleteHistory1=DeleteHistory
      ELSE
         DeleteHistory1=.TRUE.
      END IF
      
      IF (PRESENT(iprint)) THEN
         iprint1=iprint
      ELSE
         iprint1=0
      END IF
      
      IF (PRESENT(IsNeighborListInitialized)) THEN
         IsVLInitialized=IsNeighborListInitialized
      ELSE
         IsVLInitialized=.FALSE.
      END IF
      
      IF (PRESENT(WriteTrajectoryFrequency)) THEN
         WriteTrajectoryFrequency1=WriteTrajectoryFrequency
      ELSE
         WriteTrajectoryFrequency1=nmd1/10
      END IF
      
      WriteNAtoms1=MDSiml%AL%NAtoms
      WritePeriodic1=.TRUE.
      IF (PRESENT(WriteXYZTrajectoryFile)) THEN
         WriteXYZTrajectory=.TRUE.
         OPEN(UNIT=303,FILE=TRIM(WriteXYZTrajectoryFile))
         IF (PRESENT(WriteNAtoms)) WriteNAtoms1=WriteNAtoms
         IF (PRESENT(WritePeriodic)) WritePeriodic1=WritePeriodic
      ELSE
         WriteXYZTrajectory=.FALSE.
      END IF
      
      IF (PRESENT(WriteEnergyTrajectoryFile)) THEN
         WriteEnergyTrajectory=.TRUE.
         OPEN(UNIT=304,FILE=TRIM(WriteEnergyTrajectoryFile))
         WRITE(UNIT=304,FMT='("   MDIter      MDTime    Temperature (set pt)  PotentialE", &
         "     KineticE       TotalE        Temperature(calc)")')
      ELSE
         WriteEnergyTrajectory=.FALSE.
      END IF
      
      
      AL=>MDSiml%AL
      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(*,*) "$Err>> MD object of Atom list has not been setup"
         STOP
      END IF
