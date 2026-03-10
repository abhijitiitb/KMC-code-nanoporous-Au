MODULE NEB_package
   !returns the minimum energy chain of states with energies
   USE VARIABLE_TYPE
   USE db_manipulate
   USE utilities
   USE OptimizationAL
   USE OptimizationNEB
   IMPLICIT NONE
   LOGICAL, PRIVATE :: IsString,IsSplineFit,IsGlobalOptimize,IsClimbing,IsGrowing,TuneSpring
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !suggested changes:
   !instead of using state1 and state2 use chos
   !or can we have both options?
   SUBROUTINE NEB(state1,state2,springconstant,imode,omode,nimages,interpmode,ITMAX, &
      TSEnergy,ntransitionstates,TS,iunitEnergyDiagram,NEBchos,ReuseVL,ChOSFile,errorstatus)
      !imode=1 (NEB)
      !imode=2 (DoubleNEB)
      !imode=3 (SM)
      !add 10 to imode to get spline fits instead of upwind
      !add 100 to imode to get climbing image
      !add 1000 to imode to tune the spring constant -- this can be useful to get the spring to be less floppy
      !omode- optimization mode
      !omode=1 (SD)
      !omode=2 (CG)
      !omode=3 (LBFGS)
      !omode=4 (dynamical-Euler)
      !omode=5 (dynamical-RK method)
      !omode=6 (QuickMin)
      !omode=7 (FIRE)
      !interpmode=1 (Linear interpolation of images)
      !interpmode=2 (Growing mode)
      !add 10 to omode to get global
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: state1,state2,AL,TS
      TYPE(ChOSContainer), POINTER :: chos
      TYPE(ChOSContainer), POINTER, OPTIONAL :: NEBchos
      INTEGER, OPTIONAL :: nimages,imode,omode,interpmode,iunitEnergyDiagram
      INTEGER :: nimages1,omode1,imode1,interpmode1,ITMAX,ntransitionstates,i,errorstatus,NAtoms
      INTEGER :: ntighterspringattempts
      REAL(dp) :: springconstant,springconstant1,TSEnergy
      LOGICAL :: IsMatch,NotConverged,PrintEnergyDiagram,ReuseVL1
      CHARACTER(len=100) :: filename
      CHARACTER(len=100), OPTIONAL :: ChOSFile
      LOGICAL, OPTIONAL :: ReuseVL
      
      errorstatus=0
      
      IF (.NOT. ASSOCIATED(state1)) THEN
         WRITE(6,*) "$Err>> Initial state in NEB calculation is not associated"
         STOP
      END IF

      IF (.NOT. ASSOCIATED(state2)) THEN
         WRITE(6,*) "$Err>> Final state in NEB calculation is not associated"
         STOP
      END IF
      
      IF (PRESENT(iunitEnergyDiagram)) THEN
         PrintEnergyDiagram=.TRUE.
      ELSE
         PrintEnergyDiagram=.FALSE.
      END IF
      
      IF (PRESENT(ReuseVL)) THEN
         ReuseVL1=ReuseVL
      ELSE
         ReuseVL1=.FALSE.
      END IF
      NAtoms=state1%NAtoms
      CALL AddVerletList(state1,NAtoms=NAtoms)
      CALL OptimizeAL(state1,3,NAtoms=state1%NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
         ITMAX=100,MaxDisplacement=0.05_dp)
      !CALL OptimizeAL(state1,3,NAtoms=NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
      !   ITMAX=500,MaxDisplacement=0.05_dp) !required later to find TS Energy
      WRITE(6,*) "State1 optimized energy:",state1%PotentialEnergy
      
      CALL AddVerletList(state2,NAtoms=NAtoms)
      CALL OptimizeAL(state2,3,NAtoms=state2%NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
         ITMAX=100,MaxDisplacement=0.05_dp)
      !CALL OptimizeAL(state2,3,NAtoms=NAtoms,ftol=1.e-8_dp,gtol=1.e-4_dp,xtol=1.e-4_dp, &
      !   ITMAX=500,MaxDisplacement=0.05_dp)
      WRITE(6,*) "State2 optimized energy:",state2%PotentialEnergy
      
      TSEnergy=-9999._dp
      
      IsMatch=CompareALToAL(state1,state2,xtol=0.5_dp)
      
      !makes sure that certain SystemContainer elements arent identical
      IF (IsMatch) THEN
         WRITE(6,*) TRIM(TxtHeader)//"$Err>> images matched in ChOS"
         errorstatus=10
         RETURN
      END IF
      
      IF (PRESENT(imode)) THEN
         imode1=imode
      ELSE
         imode1=1 !NEB method - default
      END IF

      IF (PRESENT(omode)) THEN
         omode1=omode
      ELSE
         omode1=3 !CG method - default
      END IF

      IF (PRESENT(interpmode)) THEN
         interpmode1=interpmode
      ELSE
         interpmode1=1
      END IF
      
      IF (PRESENT(nimages)) THEN
         nimages1=nimages
         IF (nimages1<4) THEN
            WRITE(6,*) "$Err>> Number of images for NEB is too small"
            STOP
         END IF
      ELSE
         IF (interpmode1==2) THEN
            nimages1=4
         ELSE
            nimages1=10
         END IF
      END IF
      
      springconstant1=springconstant
      
      !--------------------------------------
      CALL SetupLogicals(imode1,omode1,interpmode1) !determine what type of calculation is being done
      !---------------------------------------

      NULLIFY(chos) !ensure chos is empty
      CALL AddChOSNEB(state1,state2,chos,interpmode=interpmode1,nimages=nimages1) !NEB/String mode
      chos%SpringConst=springconstant1
      ntransitionstates=0
      
      WRITE(6,*) "Writing initial interpolated MEP"
      filename="chos1.xyz"
      CALL WriteChOS(chos%AL,filename)
      
      NotConverged=.TRUE.
      errorstatus=0

      ntighterspringattempts=0
      IF (IsGlobalOptimize) THEN !global optimization
         IF (IsString) THEN  !string method
            IF (IsGrowing) THEN !growing string method
               IF (IsSplineFit) THEN
               ELSE
               END IF
            ELSE !regular string method
               IF (IsSplineFit) THEN
               ELSE
               END IF
            END IF
         ELSE !NEB
            IF (IsGrowing) THEN !growing NEB
            ELSE
               SELECT CASE (imode1)
               CASE (1)
                  DO WHILE (NotConverged .AND. ntighterspringattempts<5)
                     WRITE(6,*) "NEB>> spring constant used ",springconstant1
                     CALL OptimizeChOSGlobal(chos,omode1,CI_NEB=IsClimbing, &
         ftol=1.0e-9_dp,gtol=2.0e-3_dp,xtol=1.0e-5_dp,ITMAX=ITMAX,errorstatus=errorstatus, &
         ReuseVL=ReuseVL1)
                     IF (errorstatus==0) THEN
                        NotConverged=.FALSE.
                     ELSEIF (TuneSpring .AND. errorstatus/=3) THEN !errorstatus==3 implies initial atoms are overlapping
                        ntighterspringattempts=ntighterspringattempts+1
                        WRITE(6,*) "Using new value of spring constant in NEB ..."; CALL FLUSH(6)
                        IF (ReuseVL1) THEN
                           AL=>chos%AL%NextNeigh
                           DO WHILE (ASSOCIATED(AL))
                              NULLIFY(AL%LL)
                              NULLIFY(AL%VL)
                              AL=>AL%NextNeigh
                           END DO
                        END IF
                        CALL DeleteChos(chos,.TRUE.)
                        NULLIFY(chos)
                        CALL AddChOSNEB(state1,state2,chos,interpmode=interpmode1,nimages=nimages1) !NEB/String mode
                        springconstant1=springconstant1*1.414_dp !prevent images from coming closer
                        WRITE(6,*) "NEB>> Increasing spring constant to ",springconstant1; CALL FLUSH(6)
                        chos%SpringConst=springconstant1  !and con
                     ELSE !there was a problem with force evaluation
                        WRITE(6,*) "Error>> MEP did not converge"
                        RETURN
                     END IF
                  END DO
               CASE (2);  !double NEB
               END SELECT
            END IF
         END IF
      ELSE !image-by-image optimization
         IF (IsString) THEN
            IF (IsGrowing) THEN
               IF (IsSplineFit) THEN
               ELSE
               END IF
            ELSE
               IF (IsSplineFit) THEN
               ELSE
               END IF
            END IF
         ELSE
            IF (IsGrowing) THEN
            ELSE
               SELECT CASE (imode1)
               CASE (1); CALL OptimizeChOS(chos,omode1,CI_NEB=IsClimbing,ITMAX=ITMAX,errorstatus=errorstatus)
               CASE (2);  !double NEB
               END SELECT
            END IF
         END IF
      END IF
      
      WRITE(6,*) "Writing converged MEP"
      filename="chos2.xyz"
      IF (PRESENT(ChOSFile)) filename=ChOSFile
      CALL WriteChOS(chos%AL,filename)
      
      !Find the transition state energy
      AL=>chos%AL%NextNeigh
      ntransitionstates=0
      TSEnergy=-9999._dp
      IF (PrintEnergyDiagram) WRITE(iunitEnergyDiagram,*) 1,chos%AL%PotentialEnergy-chos%AL%PotentialEnergy
      DO i=2,chos%nimages-1
         IF (PrintEnergyDiagram) WRITE(iunitEnergyDiagram,*) i, &
            AL%PotentialEnergy-chos%AL%PotentialEnergy
         IF (AL%PotentialEnergy>=AL%PrevNeigh%PotentialEnergy .AND. &
           AL%PotentialEnergy>=AL%NextNeigh%PotentialEnergy) THEN
            ntransitionstates=ntransitionstates+1
            IF (AL%PotentialEnergy-chos%AL%PotentialEnergy>TSEnergy) THEN
               TSEnergy=MAX(TSEnergy,AL%PotentialEnergy-chos%AL%PotentialEnergy)
               WRITE(UNIT=6,FMT='(" NEB>> Transition state energy: ",F8.3," eV")') TSEnergy
               CALL MakeSize(TS,AL%NAtoms)
               CALL Copy(AL,TS) !need to initialize TS of size Natoms
               !TS=>AL   !using copyal (above) is better than just pointing to AL since soon we shall delete the entire list
            END IF
         END IF
         
         IF (ReuseVL1) THEN
            NULLIFY(AL%LL)
            NULLIFY(AL%VL)
         END IF
         
         AL=>AL%NextNeigh
      END DO
      
      IF (ReuseVL1) THEN
         NULLIFY(AL%LL)
         NULLIFY(AL%VL)
      END IF
      
      IF (PrintEnergyDiagram) WRITE(iunitEnergyDiagram,*) chos%nimages, &
         AL%PotentialEnergy-chos%AL%PotentialEnergy
      
      WRITE(UNIT=6,FMT='("Number of transition states found: ",I2)') ntransitionstates
      
      springconstant=springconstant1 !let the user know what springconstant was used
      
      IF (PRESENT(NEBchos)) THEN
         IF (ASSOCIATED(NEBchos)) THEN
            WRITE(6,*) "$Err>> NEBChoS provided to NEB is associated"
            STOP
         END IF
         NEBchos=>chos
         NULLIFY(chos)
      ELSE
         WRITE(*,*) "Deleting chos ......?"
         CALL DeleteChOS(chos,.TRUE.)
      END IF
   END SUBROUTINE NEB
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetupLogicals(imode,omode,interpmode)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: imode,omode,interpmode
      INTEGER :: ierr
      
      IsString=(MOD(imode,10)==3)
      IsSplineFit=MOD(imode,100)/10==1
      IsClimbing=MOD(imode,1000)/100==1
      IsGrowing=(interpmode==2)
      IsGlobalOptimize=(omode/10==1)
      TuneSpring=MOD(imode,10000)/1000==1
      
      imode=MOD(imode,10)
      omode=MOD(omode,10)
      
      ierr=0
      IF (IsString .AND. IsGlobalOptimize) ierr=1
      IF (.NOT. IsString .AND. IsSplineFit) ierr=1
      IF (IsString .AND. IsClimbing) ierr=1
      IF (.NOT. IsString .AND. IsGrowing) ierr=1
      
         WRITE(UNIT=6,FMT='(" IsString:",L1)') IsString
         WRITE(UNIT=6,FMT='(" IsSplineFit:",L1)') IsSplineFit
         WRITE(UNIT=6,FMT='(" IsClimbing:",L1)') IsClimbing
         WRITE(UNIT=6,FMT='(" IsGrowing:",L1)') IsGrowing
         WRITE(UNIT=6,FMT='(" IsGlobalOptimize:",L1)') IsGlobalOptimize
      
      IF (ierr/=0) THEN
         WRITE(6,*) "$Err>> Options chosen are incompatible"
         STOP
      END IF
   END SUBROUTINE SetupLogicals
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE NEB_package
