MODULE BuckinghamInitialize 
   
   USE VARIABLE_TYPE
   USE checkdb
   IMPLICIT NONE
   INTEGER, PARAMETER :: PotentialType=4 !Buckingham is potential type 4
   TYPE(BuckinghamPotential), POINTER :: Buckingham
   
   CONTAINS

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   SUBROUTINE BuckinghamSetUp(AL) 
      !is to be only called by PotentialPackageInitialize
      !... this means all variables have been properly initialized
      
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: species1,species2,NSpecies,section
      LOGICAL :: Found
      
      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(6,*) "Error>> Atomlist is not associated"
         STOP
      END IF
      IF (.NOT. ASSOCIATED(AL%Potential%Buckingham)) RETURN
      
      NSpecies=NSpecies_global

      Buckingham=>AL%Potential%Buckingham
      Found=.FALSE.
      !----------------Add the potentials----------------------
      
      DO species1=1,NSpecies

         DO species2=species1,NSpecies
            section=((2*NSpecies_global-species1)*(species1-1))/2+species2
            !IF (Buckingham%EnabledPP(section)) THEN
               CALL BuckinghamReadPP(species1,species2) !store maxpotentialcutoff for Buckingham
               AL%Potential%MaxPotentialCutoff(PotentialType)= &
                  MAX(AL%Potential%MaxPotentialCutoff(PotentialType),Buckingham%RangePP(2*section))
               Found=.TRUE.
            !END IF

         END DO

      END DO

      AL%Potential%VLBuffer=MAXVAL(AL%Potential%MaxPotentialCutoff)/50._dp !0.1 !angstrom
      
      IF (Found) WRITE(*,*) "Buckingham>> Buckingham potentials are now available"
      
   END SUBROUTINE BuckinghamSetUp

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE BuckinghamReadPP(species1,species2) !note that species1<=species2
      IMPLICIT NONE
      INTEGER :: species1,species2,r1,r2,npt,section,i,atmno1,atmno2,errorstatus
      REAL(dp) :: x(2),dxi,conv
      CHARACTER(len=100) :: FileName
      LOGICAL :: lexist
      
      section=((2*NSpecies_global-species1)*(species1-1))/2+species2
      IF (Buckingham%EnabledPP(section)) THEN
         FileName = Buckingham%TableLocationPP(section)
         
         INQUIRE(FILE=TRIM(FileName),EXIST=lexist)
         IF (.NOT. lexist) THEN
            WRITE(6,*) "$Err>> Unable to open file ",TRIM(FileName)
         END IF
         
         OPEN(UNIT=UnitTmp1,FILE=TRIM(FileName))
         READ(UnitTmp1,*) npt,x(1),x(2),dxi
         IF (npt<1000) THEN
            WRITE(*,*) "$Err>> Too few interp points"
            STOP
         END IF
         IF (npt>MaxBuckTableSize) THEN
            WRITE(*,*) "$Err>> Too many interp points ... increase MaxBuckinghamTableSize"
            STOP
         END IF

         IF (ABS(1._dp/dxi-(x(2)-x(1))/REAL(npt))>1.e-3_dp/dxi) THEN
            WRITE(*,*) "$Err>> Unable to understand dxi",dxi
            STOP
         END IF
      
         Buckingham%SizePP(section)=npt
         Buckingham%RangePP(2*section-1:2*section)=x
         Buckingham%IncrPPinv(section)=dxi

         r1=-3+(section-1)*MaxBuckTableSize
         r2=0+(section-1)*MaxBuckTableSize
         conv=Buckingham%ConvFacPP(section)
         DO i=1,npt/4 !table contains 4 data points in each line
            r1=r1+4
            r2=MIN(r2+4,npt+(section-1)*MaxBuckTableSize)
            READ(UnitTmp1,*) Buckingham%PP(r1:r2)
            Buckingham%PP(r1:r2) = Buckingham%PP(r1:r2) * conv
         END DO
      
         CLOSE(UnitTmp1)
      
         atmno1=SpeciesDirectory_global(species1)
         atmno2=SpeciesDirectory_global(species2)
         WRITE(6,*) "Buckingham PP ",TRIM(SpeciesList%AtomicSymbol(atmno1)),"-", &
         TRIM(SpeciesList%AtomicSymbol(atmno2)),":", x
      ELSE
         Buckingham%SizePP(section)=101
         x(1)=0._dp
         x(2)=3._dp
         Buckingham%RangePP(2*section-1:2*section)=x
         Buckingham%IncrPPinv(section)=100._dp/3._dp
         r1=1+(section-1)*MaxBuckTableSize
         r2=101+(section-1)*MaxBuckTableSize
         Buckingham%PP(r1:r2)=0._dp
      END IF
   END SUBROUTINE BuckinghamReadPP
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
END MODULE BuckinghamInitialize
