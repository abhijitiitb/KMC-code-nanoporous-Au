MODULE KMCConfigurations
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur

   USE KMC_VARIABLE_TYPE
   USE KMCDbManipulate
   USE KMCUtilities
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateRedundantCfg()
   !redundant configurations at end of the list
   !atoms whose configurations does not matter to the KMC codes
   !for e.g., atoms that are not moving
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: cfg
   
      IF (ASSOCIATED(RedundantCfg)) THEN
         IF (RedundantCfg%ConfigurationIndex<0) RETURN !already exists
      END IF
      
      IF (ASSOCIATED(RecordedCfg)) THEN
         cfg=>RecordedCfg
         DO WHILE (ASSOCIATED(cfg%NextNeigh))
            cfg=>cfg%NextNeigh
         END DO
         ALLOCATE(cfg%NextNeigh); NAllottedConfiguration=NAllottedConfiguration+1
         cfg%NextNeigh%PrevNeigh=>cfg
         cfg=>cfg%NextNeigh
         NULLIFY(cfg%NextNeigh)
      ELSE
         ALLOCATE(RecordedCfg); NAllottedConfiguration=NAllottedConfiguration+1
         NULLIFY(RecordedCfg%PrevNeigh)
         NULLIFY(RecordedCfg%NextNeigh)
         cfg=>RecordedCfg
      END IF
      cfg%ConfigurationIndex=-12321
      cfg%Species=0
      cfg%NAtomsWithCfg=0
      cfg%BondOrder=1000.
      RedundantCfg=>cfg
      write(*,*) "Created redundant cfg"
   END SUBROUTINE CreateRedundantCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CheckAtomForBondOrder(Atom)
   !TRUE if the atom conifguration should be checked further
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: Atom
      LOGICAL :: CheckAtomForBondOrder

      CheckAtomForBondOrder= Atom%IsMoving .AND. Atom%BondOrder<BOCutoffForCfg
   END FUNCTION CheckAtomForBondOrder
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckConfig(Atom)
   !compares atom env to cfglist
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: Config
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCAtomList), POINTER :: AL1
      
      Config=>CompareAtomEnvToCfgList(atom) !finds configuration which matches with the atom cfg
   
      IF (ASSOCIATED(Atom%cfg)) THEN
         IF (Config%ConfigurationIndex/=Atom%cfg%ConfigurationIndex) THEN
            IF (atom%cfg%ConfigurationIndex==0) THEN
               ALLOCATE(AL1); NAllottedKMCAtomList=NAllottedKMCAtomList+1
               AL1%Atom=>Atom
               WRITE(*,*) "Such a situation where CfgIndex==0 cannot occur"
               STOP
            ELSE
               CALL RemoveFromCfg(Atom,AL1) !extract the element AL from Cfg
            END IF
            CALL AddToCfg(Config,Atom,AL1)
         END IF
      ELSE
         ALLOCATE(AL1); NAllottedKMCAtomList=NAllottedKMCAtomList+1
         AL1%Atom=>atom
         CALL AddToCfg(Config,atom,AL1)
      END IF
   END SUBROUTINE CheckConfig
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AssignConfig(Atom)
   !assigns a configuration type to our atom
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: Atom
      TYPE(KMCAtomList), POINTER :: AL1
      LOGICAL :: IsProceed
   
      IsProceed=CheckAtomForBondOrder(Atom)
      IF (IsProceed) THEN
         CALL CheckConfig(Atom)
      ELSE !add to redundant
         IF (ASSOCIATED(Atom%Cfg)) THEN
            IF (Atom%cfg%ConfigurationIndex/=RedundantCfg%ConfigurationIndex) THEN
               CALL RemoveFromCfg(Atom,AL1) !element AL contains atom
               CALL AddToCfg(RedundantCfg,Atom,AL1)
            END IF
         ELSE
            ALLOCATE(AL1); NAllottedKMCAtomList=NAllottedKMCAtomList+1
            AL1%Atom=>Atom
            CALL AddToCfg(RedundantCfg,Atom,AL1)
         END IF
      END IF
   END SUBROUTINE AssignConfig
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE KMCConfigurations
