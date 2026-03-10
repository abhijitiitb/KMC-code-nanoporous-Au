MODULE TersoffInitialize
   USE VARIABLE_TYPE
   USE checkdb
   IMPLICIT NONE
   INTEGER, PARAMETER :: PotentialType=6 !SW is potential type 1
   TYPE(TersoffPotential), POINTER :: Tersoff

   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE TersoffSetup(AL)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: species1,species2,NSpecies,section
      LOGICAL :: Found
      
      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(6,*) "Error>> Atomlist is not associated"
         STOP
      END IF
      IF (.NOT. ASSOCIATED(AL%Potential%Tersoff)) RETURN
      
      NSpecies=NSpecies_global
      
      Tersoff=>AL%Potential%Tersoff
      Found=.FALSE.
      !----------------Add the potentials----------------------
      
      DO species1=1,NSpecies

         DO species2=species1,NSpecies

            section=((2*NSpecies_global-species1)*(species1-1))/2+species2
            AL%Potential%MaxPotentialCutoff(PotentialType)= &
               MAX(AL%Potential%MaxPotentialCutoff(PotentialType),Tersoff%R(section)+Tersoff%RD(section))
            Found=.TRUE.

         END DO

      END DO

      AL%Potential%VLBuffer=MAXVAL(AL%Potential%MaxPotentialCutoff)/30._dp !0.1 !angstrom
      
      IF (Found) WRITE(*,*) "Tersoff>> Tersoff potentials are now available"
      

   END SUBROUTINE TersoffSetup
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE TersoffInitialize
