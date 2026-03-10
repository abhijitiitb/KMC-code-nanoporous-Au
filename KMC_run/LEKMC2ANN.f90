MODULE LEKMC2PM
!process model (PM)
!a) classify the type of process
!b) find the factors in the process
!c) obtain xdata and ydata and add to appropriate PM model
!PM model usage
!a) find the factors/the local environment
!b) in case the factors are new then
!  i) find the xdata and ydata and add to appropriate model
! ii) use the prediction
!c) if the environment has been seen before then recover the data
!d) if model is available then use it to make prediction
!convert atom coordinate information into factors
!also contains information about PM patterns
   USE KMC_VARIABLE_TYPE
   USE OptimPackage
   IMPLICIT NONE
   INTEGER, PARAMETER :: MaxNFactors=100,MaxNModels=100,MaxNdata=10000,MaxInteractions=21000
   
   !variables used for comparing
   REAL(sp) :: AtomPositionTolerance=0.5_sp !tolerance for detecting one mechanism from the other
   REAL(sp) :: MinInteractionStrength=0.01_sp
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE ProcessModelContainer !model for a particular mechanism and contains info about all factors
      INTEGER :: OptimizationMode=3
      TYPE(ProcessModelContainer), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL()

      !Data
      INTEGER :: ndata=0,nways=2 !ndata points available for training & testing, nways is number of interaction ways to be considered
      INTEGER, DIMENSION(MaxNFactors*MaxNdata/10) :: xdataFactors=0 !this contains the factors present in each data (compact form)
      INTEGER, DIMENSION(MaxNFactors*MaxNdata/10) :: xdataFactorSpecies=0 !this contains the factors species present in each data (compact form)
      INTEGER, DIMENSION(MaxNdata/10) :: xdataFactorRange=0 !tells how many factors are present for each data point
      TYPE(ProcessModel), POINTER :: ActivationBarrier=>NULL() !model for the activation barrier
      TYPE(ProcessModel), POINTER :: Prefactor=>NULL() !model for the prefactor
      
      !Pattern
      INTEGER :: ndisplacedatoms=0 !ndisplacedatoms is number of displaced atoms
      REAL(sp), DIMENSION(:), POINTER :: displacedatomcoord=>NULL() !displacedatomcoord is the coordinate of each displaced atom
      
      !mapping
      INTEGER :: NFactors=0 !number of constant and variable factors for the model
      INTEGER :: NVariableFactors=0 !number of variable factors
      INTEGER, DIMENSION(MaxNFactors) :: NdataFactors=0 !number of data points collected which involve these factors -- this is used for averaging
      INTEGER, DIMENSION(MaxNFactors) :: VariableFactors=0
      REAL(sp), DIMENSION(3*MaxNFactors) :: FactorAverageCoordinates=0. !indexing of factors is based on the order of coordinates given here
      REAL(sp), DIMENSION(3*MaxNFactors) :: FactorAverageMaxDisplacement=0.
      LOGICAL, DIMENSION(MaxNFactors*NSpeciesType) :: FactorSpeciesPresent=.FALSE. !contains information about which species are present
      !here is the indexing for FactorSpeciesPresent (assuming NSpeciesType=4)
      !     FACTOR 1    |    FACTOR 2   | .... factor #
      !  1  2   3   4   |  1  2  3  4   |  .... species #
      !  1  2   3   4      5  6  7  8         factor species index
      LOGICAL, DIMENSION(MaxNFactors) :: FactorIsConstant=.TRUE.
      
      INTEGER :: ModelIndex=0
   END TYPE ProcessModelContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE ProcessModel
      REAL(sp), DIMENSION(21000) :: ydata
      !Here is the way interactions involving factors are implemented for activation barrier
      !  interactions(i) is the interaction strength, interactionfactorrange(i) is the portion which contains the factors
      ! interactionfactor(interactionfactorrange(i-1)+1:interactionfactorrange(i)) are the factors in this interaction
      ! interactionspecies(interactionfactorrange(i-1)+1:interactionfactorrange(i)) are the species in this interaction
      ! Number factors  2-way interactions#  3-way interactions#  2^Numberfactors
      !      50              1275                20875               1e15
      !     100              5050               166750               1e30
      !     200             20100              1333500               1e60
      !Multiply number of species -- note that not all factors will be included
      INTEGER :: ninteractions=0, ninteractionsactive=0 !number of total and active interactions
      REAL(sp), DIMENSION(MaxInteractions+1) :: interactions=0. !contains the interactions - e.g., interaction for [21,28,39]
      !note that the first element of interactions is bias term
      INTEGER, DIMENSION(3*MaxInteractions) :: interactionfactor=0,interactionspecies=0 !corresponding factors and species involved - e.g., factors are [21,28,39] and species are [1,1,2]
      INTEGER, DIMENSION(MaxInteractions) :: interactionfactorrange=0 !gives the how many factors are involved for each interaction type - e.g., count is 3
      LOGICAL, DIMENSION(MaxInteractions) :: interactionisactive=.FALSE. !specifies whether the interaction is active
   END TYPE ProcessModel
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE(ProcessModelContainer), POINTER :: PrMdlList=>NULL()
   TYPE(ProcessModelContainer), POINTER :: ActivePMModel=>NULL() !this is used for training
   INTEGER :: NModels=0 !number of PM models
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !temporary factor variables used by the code
   INTEGER :: TmpNFactors=0
   INTEGER, DIMENSION(MaxNFactors) :: TmpFactorSpecies=0
   LOGICAL, DIMENSION(MaxNFactors) :: TmpFactorIsConstant=.TRUE.
   INTEGER, DIMENSION(MaxNFactors) :: TmpFactors=0 !this is what is needed by PM for training/prediction -- "Factors" contains the index to the actual factor
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !            PM model building routines
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddProcess2PM(ndisplacedatoms,displacedatomcoord, & !process recognition data
      nprcatoms,prcatominitialcoord,prcatommaxdisplacement,prcatomspecies, & !xdata
      activationbarrier,prefactor) !ydata
   !This subroutine is called when the process information needs to be added to later train the ANN model
   !convert the process atom information into PM information that can be later used for training PM
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PMModel
      INTEGER :: ndisplacedatoms,nprcatoms !number of displaced atoms, number of process atoms as identified by LE-KMC
              !process atoms have displacements greater than minimum displacement 
      REAL(sp), DIMENSION(3*ndisplacedatoms) :: displacedatomcoord  !initial coordinate of each displaced atom
      REAL(sp), DIMENSION(3*nprcatoms) :: prcatominitialcoord !initial coordinate
      REAL(sp), DIMENSION(nprcatoms) :: prcatommaxdisplacement !maximum displacement
      INTEGER, DIMENSION(nprcatoms) :: prcatomspecies !species of displaced atoms
      REAL(sp) :: activationbarrier,prefactor
      INTEGER :: i
      LOGICAL :: DoTraining
      
      !step 1:find the matching pattern
      PMModel=>PatternMatchPM(ndisplacedatoms,displacedatomcoord)
      !step 2:create the factors list for the process
      IF (ASSOCIATED(PMModel)) THEN !model found
         CALL ExtractFactors(PMModel,nprcatoms,prcatominitialcoord,prcatomspecies) !create the factor index array
      ELSE !model not found
         CALL AddProcessModelContainer(PMModel,ndisplacedatoms,displacedatomcoord) !add process model and process pattern
      END IF
      
      !step 3: store factors and observables in PM
      IF (DataIsUnique()) CALL AddPMData(PMModel,nprcatoms,prcatominitialcoord, &
         prcatommaxdisplacement,prcatomspecies, & !add data for training
         activationbarrier,prefactor)
      
      !step 4: if conditions are right for training then perform training
      !Conditions:
      !a) 
      DoTraining=.FALSE.
      IF (DoTraining) THEN
         CALL IdentifyVariableFactors(PMModel)
         DoTraining=PMModel%NVariableFactors<2*PMModel%ndata
      END IF
      IF (DoTraining) THEN
         CALL SetUpInteractions(PMModel)
         CALL TrainModel(PMModel)
      END IF
      
   END SUBROUTINE AddProcess2PM
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION PatternMatchPM(ndisplacedatoms,displacedatomcoord)
   !finds whether the atomic process matches a particular PM model mechanism
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PatternMatchPM
      INTEGER :: ndisplacedatoms
      REAL(sp), DIMENSION(3*ndisplacedatoms) :: displacedatomcoord
      LOGICAL :: Found
      
      PatternMatchPM=>PrMdlList
      DO WHILE (ASSOCIATED(PatternMatchPM))
         Found=MatchDisplacedAtoms(ndisplacedatoms,displacedatomcoord,PatternMatchPM)
         IF (Found) RETURN
         PatternMatchPM=>PatternMatchPM%NextNeigh
      END DO
   END FUNCTION PatternMatchPM
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION MatchDisplacedAtoms(ndisplacedatoms,displacedatomcoord,PMModel)
   !finds whether the displacement matches with the information stored in a particular PMModel
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PMModel
      INTEGER :: ndisplacedatoms,iatom,jatom
      REAL(sp), DIMENSION(3*ndisplacedatoms) :: displacedatomcoord
      REAL(sp), DIMENSION(:), POINTER :: displacedatomcoord1
      REAL(sp) :: ri(3),rj(3)
      !INTEGER, DIMENSION(ndisplacedatoms) :: atomspeciesdisplacedatoms
      LOGICAL :: MatchDisplacedAtoms
      
      MatchDisplacedAtoms=.FALSE.
      
      IF (PMModel%ndisplacedatoms/=ndisplacedatoms) RETURN
      
      displacedatomcoord1=>PMModel%displacedatomcoord
      DO iatom=1,ndisplacedatoms
         ri=displacedatomcoord(3*iatom-2:3*iatom)
         DO jatom=1,ndisplacedatoms
            rj=displacedatomcoord1(3*jatom-2:3*jatom) !relative coordinate
            MatchDisplacedAtoms=ALL(ABS(ri-rj)<AtomPositionTolerance)
            IF (MatchDisplacedAtoms) EXIT
         END DO
         IF (.NOT. MatchDisplacedAtoms) THEN !is it possible that the match could not be done because the displacement is close to the tolerance
         END IF
         IF (.NOT. MatchDisplacedAtoms) RETURN !match has failed
      END DO
   END FUNCTION MatchDisplacedAtoms
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddProcessModelContainer(PMModel,ndisplacedatoms,displacedatomcoord)
   !adds a pattern to an PM model and creates the model
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PMModel
      INTEGER :: ndisplacedatoms
      REAL(sp), DIMENSION(3*ndisplacedatoms) :: displacedatomcoord
      
      ALLOCATE(PMModel)
      IF (ASSOCIATED(PrMdlList)) THEN
         PrMdlList%PrevNeigh=>PMModel
         PMModel%NextNeigh=>PrMdlList
      END IF
      PrMdlList=>PMModel
      NModels=NModels+1
      IF (NModels>MaxNModels) THEN
         WRITE(6,*) "Err>> Number of process models exceeds maximum allowable limit"
         STOP
      END IF
      
      !create the pattern
      PMModel%ndisplacedatoms=ndisplacedatoms
      PMModel%displacedatomcoord(1:3*ndisplacedatoms)=displacedatomcoord
      
      !add the ANN models
      ALLOCATE(PMModel%ActivationBarrier)
      ALLOCATE(PMModel%Prefactor)
   END SUBROUTINE AddProcessModelContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteProcessModelContainer(PMModel)
   !delete a model from PMSuite
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PMModel
      
      WRITE(6,*) "Err>> DeleteProcessModelContainer has not been implemented"
      STOP
   END SUBROUTINE DeleteProcessModelContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ExtractFactors(PMModel,nprcatoms,prcatominitialcoord,prcatomspecies)
   !returns the list of atom coordinates and their species
   !displaced atoms are not included as part of the atom information
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PMModel
      INTEGER :: nprcatoms
      REAL(sp), DIMENSION(:), POINTER :: AtomCoord
      REAL(sp), DIMENSION(3*nprcatoms) :: prcatominitialcoord
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      INTEGER, DIMENSION(nprcatoms) :: prcatomspecies
      REAL(sp) :: coord(3)
      INTEGER :: nf,i,species
      
      TmpNFactors=nprcatoms !number of factors in the current process
      nf=PMModel%NFactors !number of factors known to the process model
      AtomCoord=>PMModel%FactorAverageCoordinates(1:3*nf)
      AtomSpecies=>PMModel%xdatafactorSpecies(1:nf)
      DO i=1,nprcatoms
         coord=prcatominitialcoord(3*i-2:3*i)
         species=prcatomspecies(i)
         TmpFactors(i)=LocateCoord()
      END DO
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxx
      FUNCTION LocateCoord()
      !checks whether this coord is unique and has not been added as a process atom
         IMPLICIT NONE
         INTEGER :: LocateCoord,j,nf1
         
         LocateCoord=0
         DO j=1,nf
            IF (ALL(ABS(coord-AtomCoord(3*j-2:3*j))<AtomPositionTolerance .AND. species==AtomSpecies(j))) THEN
               LocateCoord=j
               RETURN
            END IF
         END DO
         !add new factor -- this will be useful in future
         nf1=PMModel%NFactors+1
         PMModel%NFactors=nf1
         PMModel%FactorAverageCoordinates(3*nf1-2:3*nf1)=coord
         PMModel%xdatafactorSpecies(nf1)=species
         LocateCoord=nf1
      END FUNCTION LocateCoord
      !xxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE ExtractFactors
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddPMData(PMModel,nprcatoms,prcatominitialcoord, &
      prcatommaxdisplacement,prcatomspecies,activationbarrier,prefactor)
   !add data for training
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PMModel
      INTEGER :: nprcatoms
      REAL(sp), DIMENSION(3*nprcatoms) :: prcatominitialcoord
      REAL(sp), DIMENSION(nprcatoms) :: prcatommaxdisplacement
      INTEGER, DIMENSION(nprcatoms) :: prcatomspecies
      REAL(sp), DIMENSION(:), POINTER :: FactorAverageCoordinates,FactorAverageMaxDisplacement
      LOGICAL, DIMENSION(:), POINTER :: FactorSpeciesPresent
      LOGICAL, DIMENSION(:), POINTER :: FactorIsConstant
      INTEGER, DIMENSION(:), POINTER :: NdataFactors
      REAL(sp) :: activationbarrier,prefactor
      INTEGER :: istart,istop,ndata,i,ifactor,nfactors,ispecies
      
      ndata=PMModel%ndata
      IF (ndata==0) THEN !factors have not been initialized
         FactorSpeciesPresent=>PMModel%FactorSpeciesPresent
         DO ifactor=1,nprcatoms
            TmpFactors(ifactor)=i
            ispecies=prcatomspecies(ifactor)
            FactorSpeciesPresent((ifactor-1)*NSpeciesType+ispecies)=.TRUE.
         END DO
         TmpFactorSpecies(1:nprcatoms)=prcatomspecies
         
         PMModel%NFactors=nprcatoms
         PMModel%FactorAverageCoordinates(1:3*nprcatoms)=prcatominitialcoord
         PMModel%xdatafactorSpecies(1:nprcatoms)=prcatomspecies
         !create the FactorIndex array so that we can store the data next
         NFactors=nprcatoms
      ELSE !update factor information
         TmpFactorIsConstant=.FALSE.
         nfactors=PMModel%NFactors
         FactorAverageCoordinates=>PMModel%FactorAverageCoordinates
         FactorAverageMaxDisplacement=>PMModel%FactorAverageMaxDisplacement
         NdataFactors=>PMModel%NdataFactors
         DO i=1,TmpNFactors
            ifactor=TmpFactors(i)
            FactorAverageCoordinates(ifactor)=(FactorAverageCoordinates(ifactor)*NdataFactors(ifactor)+ &
               prcatominitialcoord(i))/REAL(NdataFactors(ifactor)+1)
            FactorAverageMaxDisplacement(ifactor)=(FactorAverageMaxDisplacement(ifactor)*NdataFactors(ifactor)+ &
               prcatommaxdisplacement(i))/REAL(NdataFactors(ifactor)+1)
            FactorSpeciesPresent((ifactor-1)*NSpeciesType+ispecies)=.TRUE.
            TmpFactorIsConstant(ifactor)=.TRUE.
         END DO
         FactorIsConstant=>PMModel%FactorIsConstant
         FactorIsConstant(1:nfactors)=FactorIsConstant(1:nfactors) .AND. &
            TmpFactorIsConstant(1:nfactors)
      END IF
      
      !add x-data
      istart=1
      IF (ndata>0) istart=istart+PMModel%xdatafactorRange(ndata)
      istop=istart-1+nprcatoms
      PMModel%xdatafactors(istart:istop)=TmpFactors
      PMModel%xdatafactorSpecies(istart:istop)=TmpFactorSpecies
      PMModel%xdatafactorRange(ndata+1)=istop
      
      ndata=ndata+1
      PMModel%ndata=ndata
      PMModel%ActivationBarrier%ydata(ndata)=activationbarrier
      PMModel%Prefactor%ydata(ndata)=prefactor
      
   END SUBROUTINE AddPMData
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION DataIsUnique()
   !checks whether the data stored in the tmpfactor is unique
      IMPLICIT NONE
      LOGICAL :: DataIsUnique
      
      DataIsUnique=.FALSE.
      WRITE(*,*) "DataIsUnique has not been setup"
      STOP
   END FUNCTION DataIsUnique
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !            PM training routines
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE IdentifyVariableFactors(PMModel)
   !finds number of variable factors
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PMModel
      INTEGER, DIMENSION(:), POINTER :: VariableFactors
      LOGICAL, DIMENSION(:), POINTER :: FactorIsConstant
      INTEGER :: NVariableFactors,i
      
      NVariableFactors=0
      FactorIsConstant=>PMModel%FactorIsConstant
      VariableFactors=>PMModel%VariableFactors
      DO i=1,PMModel%NFactors
         IF (.NOT. FactorIsConstant(i)) THEN
            NVariableFactors=NVariableFactors+1
            VariableFactors(NVariableFactors)=i !ith factor is added as VariableFactor
         END IF
      END DO
      PMModel%NVariableFactors=NVariableFactors
   END SUBROUTINE IdentifyVariableFactors
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetUpInteractions(PMModel)
   !will go through the variable factors and set-up the interactions
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PMModel
      INTEGER, DIMENSION(MaxNFactors) :: VariableFactors,VariableFactorsSpecies
      LOGICAL, DIMENSION(:), POINTER :: FactorSpeciesPresent,FactorIsConstant
      INTEGER, DIMENSION(:), POINTER :: act_interactionfactor,act_interactionspecies,act_interactionrange
      INTEGER, DIMENSION(:), POINTER :: pre_interactionfactor,pre_interactionspecies,pre_interactionrange
      INTEGER :: ifactor,ispecies,interactioncount,nways
      INTEGER :: fac1,fac2,fac3,fac4,ninteractions,pos
      
      !collect the variable factors present
      interactioncount=0
      DO ifactor=1,PMModel%NFactors
         IF (FactorIsConstant(ifactor)) THEN
            DO ispecies=1,NSpeciesType
               IF (FactorSpeciesPresent((ifactor-1)*NSpeciesType+ispecies)) THEN !add
                  interactioncount=interactioncount+1
                  VariableFactors(interactioncount)=VariableFactors(interactioncount)+1
                  VariableFactorsSpecies(interactioncount)= &
                     VariableFactorsSpecies(interactioncount)+1
               END IF
            END DO
         END IF
      END DO
      
      !create the interactions -- this will overwrite all the information in ydata%interactions
      act_interactionfactor=>PMModel%ActivationBarrier%interactionfactor
      act_interactionspecies=>PMModel%ActivationBarrier%interactionspecies
      act_interactionrange=>PMModel%ActivationBarrier%interactionfactorrange
      pre_interactionfactor=>PMModel%Prefactor%interactionfactor
      pre_interactionspecies=>PMModel%Prefactor%interactionspecies
      pre_interactionrange=>PMModel%Prefactor%interactionfactorrange
      
      nways=PMModel%nways
      ninteractions=0
      IF (nways>=1) THEN !add one-way interactions
         DO fac1=1,interactioncount
            ninteractions=ninteractions+1
            pos=pos+1
            act_interactionfactor(ninteractions)=VariableFactors(fac1)
            act_interactionspecies(ninteractions)=VariableFactorsSpecies(fac1)
            pre_interactionfactor(pos)=VariableFactors(fac1)
            pre_interactionspecies(pos)=VariableFactorsSpecies(fac1)
            act_interactionrange(pos)=pos
            pre_interactionrange(pos)=pos
         END DO
      END IF
      
      pos=ninteractions
      IF (nways>=2) THEN !add two-way interactions
         DO fac1=1,interactioncount-1
            DO fac2=2,interactioncount
               ninteractions=ninteractions+1
               pos=pos+1
               act_interactionfactor(pos)=VariableFactors(fac1)
               act_interactionspecies(pos)=VariableFactorsSpecies(fac1)
               pre_interactionfactor(pos)=VariableFactors(fac1)
               pre_interactionspecies(pos)=VariableFactorsSpecies(fac1)
               pos=pos+1
               act_interactionfactor(pos)=VariableFactors(fac2)
               act_interactionspecies(pos)=VariableFactorsSpecies(fac2)
               pre_interactionfactor(pos)=VariableFactors(fac2)
               pre_interactionspecies(pos)=VariableFactorsSpecies(fac2)
               act_interactionrange(pos)=pos
               pre_interactionrange(pos)=pos
            END DO
         END DO
      END IF
      IF (nways>=3) THEN !add three-way interactions
         DO fac1=1,interactioncount-2
            DO fac2=2,interactioncount-1
               DO fac3=3,interactioncount
                  ninteractions=ninteractions+1
                  pos=pos+1
                  act_interactionfactor(pos)=VariableFactors(fac1)
                  act_interactionspecies(pos)=VariableFactorsSpecies(fac1)
                  pre_interactionfactor(pos)=VariableFactors(fac1)
                  pre_interactionspecies(pos)=VariableFactorsSpecies(fac1)
                  pos=pos+1
                  act_interactionfactor(pos)=VariableFactors(fac2)
                  act_interactionspecies(pos)=VariableFactorsSpecies(fac2)
                  pre_interactionfactor(pos)=VariableFactors(fac2)
                  pre_interactionspecies(pos)=VariableFactorsSpecies(fac2)
                  pos=pos+1
                  act_interactionfactor(pos)=VariableFactors(fac3)
                  act_interactionspecies(pos)=VariableFactorsSpecies(fac3)
                  pre_interactionfactor(pos)=VariableFactors(fac3)
                  pre_interactionspecies(pos)=VariableFactorsSpecies(fac3)
                  act_interactionrange(pos)=pos
                  pre_interactionrange(pos)=pos
               END DO
            END DO
         END DO
      END IF
      IF (nways==4) THEN !add four-way interactions
         DO fac1=1,interactioncount-3
            DO fac2=2,interactioncount-2
               DO fac3=3,interactioncount-1
                  DO fac4=4,interactioncount
                     ninteractions=ninteractions+1
                     pos=pos+1
                     act_interactionfactor(pos)=VariableFactors(fac1)
                     act_interactionspecies(pos)=VariableFactorsSpecies(fac1)
                     pre_interactionfactor(pos)=VariableFactors(fac1)
                     pre_interactionspecies(pos)=VariableFactorsSpecies(fac1)
                     pos=pos+1
                     act_interactionfactor(pos)=VariableFactors(fac2)
                     act_interactionspecies(pos)=VariableFactorsSpecies(fac2)
                     pre_interactionfactor(pos)=VariableFactors(fac2)
                     pre_interactionspecies(pos)=VariableFactorsSpecies(fac2)
                     pos=pos+1
                     act_interactionfactor(pos)=VariableFactors(fac3)
                     act_interactionspecies(pos)=VariableFactorsSpecies(fac3)
                     pre_interactionfactor(pos)=VariableFactors(fac3)
                     pre_interactionspecies(pos)=VariableFactorsSpecies(fac3)
                     pos=pos+1
                     act_interactionfactor(pos)=VariableFactors(fac4)
                     act_interactionspecies(pos)=VariableFactorsSpecies(fac4)
                     pre_interactionfactor(pos)=VariableFactors(fac4)
                     pre_interactionspecies(pos)=VariableFactorsSpecies(fac4)
                     act_interactionrange(pos)=pos
                     pre_interactionrange(pos)=pos
                  END DO
               END DO
            END DO
         END DO
      END IF
      
      PMModel%ActivationBarrier%ninteractions=ninteractions
      PMModel%ActivationBarrier%ninteractionsactive=ninteractions
      PMModel%ActivationBarrier%interactionisactive(1:ninteractions)=.TRUE.
      PMModel%Prefactor%ninteractions=ninteractions
      PMModel%Prefactor%ninteractionsactive=ninteractions
      PMModel%Prefactor%interactionisactive(1:ninteractions)=.TRUE.
   END SUBROUTINE SetUpInteractions
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE TrainModel(PMModel)
   !perform model training
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PMModel
      
      !Train the model
      WRITE(6,*) "Training the process model ..."
      
      CALL OptimizeCEModel(PMModel)
      !DO WHILE (EliminateWeakInteractions(PMModel,MinInteractionStrength))
      !   CALL OptimizeCEModel(PMModel)
      !   WRITE(6,*) "++++++++++++++++++++++++++++++"
      !END DO
      !ANN%fdata_selected=ANNOutput(ANN,ANN%xdata_selected,ANN%nselected,ANN%ndimensionsx)
      !RSq= RSquared(ANN%ydata_selected,ANN%fdata_selected,ANN%nselected)
      !WRITE(UNIT=6,FMT='(" ... RSquared value: ",F15.5)') RSq
      
      !find the maximum error
      !maxerror=0._dp
      !WRITE(6,*) "ActualActivationEnergy  PredictedActivationEnergy"
      !DO i=1,ANN%nselected
      !   maxerror=MAX(maxerror,ABS(ANN%ydata_selected(i)-ANN%fdata_selected(i)))
      !   WRITE(6,*) ANN%ydata_selected(i),ANN%fdata_selected(i)
      !END DO
      !WRITE(6,*) "Maximum error is :",maxerror
      !WRITE(6,*) "for total number of data points;",ANN%nselected
      
      !IF (RSq<0.9) THEN
      !   WRITE(6,*) "... Poor fit obtained during training stage"
      !   errorstatus=1
      !   STOP
      !END IF
      
   END SUBROUTINE TrainModel
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE OptimizeCEModel(PMModel)
   !optimize the CE model
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PMModel
      REAL(sp), DIMENSION(:), POINTER :: weights
      REAL(sp) :: fret,ftol,gtol,xtol,SclFac
      INTEGER :: omode,ninteractions,ITMAX,errorstatus,iter
      
      ActivePMModel=>PMModel
      omode=3 !optimization mode
      
      !set-up interactions for activation barrier
      ninteractions=PMModel%ActivationBarrier%ninteractions !we only want to have the possible interactions to make it easy for optimizer
      weights=>PMModel%ActivationBarrier%interactions(1:ninteractions)
      weights=0._sp
      SELECT CASE(omode)
      CASE(1) !steepest descent
         SclFac=0.01_sp !to be used for scaling forces
         ftol=1.0e-9_sp
         gtol=2.0e-6_sp
         xtol=1.0e-5_sp
         CALL SteepestDescent(weights,ftol,gtol,xtol,iter,fret,PMError,PMErrorGradient, &
            dx=0.00001_sp,ITMAX=ITMAX,IsPrint1=.TRUE.,errorstatus=errorstatus)
      CASE(2) !conjugate gradient
         SclFac=0.01_sp
         ftol=1.0e-9_sp !this should come from potential
         gtol=2.0e-6_sp
         xtol=1.0e-5_sp
         CALL ConjugateGradient(weights,ftol,gtol,xtol,iter,fret,PMError,PMErrorGradient, &
            dx=0.0001_sp,ITMAX=ITMAX,IsPrint1=.TRUE.,errorstatus=errorstatus)
      CASE(3) !LBFGS
         ftol=1.0e-9_sp
         gtol=2.0e-3_sp
         xtol=1.0e-5_sp
         ITMAX=200
         CALL Lbfgs(weights,ftol,gtol,xtol,iter,fret,PMError,PMErrorGradient, &
            IsPrint1=.TRUE.,ITMAX=ITMAX,errorstatus=errorstatus)
      !CASE(4) !dynamical Euler
         !CALL 
      !CASE(5) !Runge Kutta
      !CASE(6) !QuickMin
      !CASE(7) !FIRE
      CASE(8) !DFP
         !ftol=1.0e-9_sp
         !gtol=2.0e-6_sp
         !xtol=1.0e-5_sp
         CALL DFP(weights,ftol,gtol,xtol,iter,fret,PMError,PMErrorGradient, &
            IsPrint1=.TRUE.,errorstatus=errorstatus)
      CASE DEFAULT
         WRITE(6,*) "Error>> Optimization type is not known to OptimizeANN"
         STOP
      END SELECT
      WRITE(6,*) "Optimized activation barrier weights:"
      
      !set-up interactions for prefactor
      !weights=>PMModel%prefactor%interaction
      weights=>PMModel%prefactor%interactions
      weights=0._sp
      WRITE(6,*) "Optimized prefactor weights:"
      
   END SUBROUTINE OptimizeCEModel
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION PMError(x,errorstatus)
   !find the error resulting from the current interaction values
   !x(1) contains the bias, x(2:ninteractions+1) contains the interactions
   !many of the interactions could be 0 as they are not active
      IMPLICIT NONE
      REAL(sp) :: PMError
      REAL(sp), DIMENSION(:), INTENT(IN) :: x
      INTEGER :: i,ndata,errorstatus
      
      PMError=0._sp
      DO i=1,ndata
         
      END DO
   END FUNCTION PMError
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION PMErrorGradient(x,errorstatus)
      IMPLICIT NONE
      REAL(sp), DIMENSION(:), INTENT(IN) :: x
      REAL(sp), DIMENSION(SIZE(x)) :: PMErrorGradient
      INTEGER :: errorstatus
      
      PMErrorGradient=0._sp
   END FUNCTION PMErrorGradient
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ArrayElementForInteraction(PMModel,fac,nways)
   !finds the element in the interaction array corresponding to fac
   !fac is a sorted array of size 4 containing the variable factors
   !nways is the number of interactions
      IMPLICIT NONE
      TYPE(ProcessModelContainer), POINTER :: PMModel
      INTEGER :: nways,fac(4),ArrayElementForInteraction
      
      SELECT CASE(nways)
      CASE(1) !1-way interaction
         ArrayElementForInteraction=fac(1)
      CASE(2) !2-way interaction
         ArrayElementForInteraction=PMModel%NVariableFactors
      CASE(3) !3-way interaction
      CASE(4)
      CASE DEFAULT
         WRITE(6,*) "Err>> Number of ways of interaction should lie between 1-4 ..."
         WRITE(6,*) "     ...    in ArrayElementForInteraction"
         STOP
      END SELECT
   END FUNCTION ArrayElementForInteraction
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !            PM model usage routines
END MODULE LEKMC2PM   
