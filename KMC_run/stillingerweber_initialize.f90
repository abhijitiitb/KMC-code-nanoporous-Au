MODULE StillingerWeberInitialize
   USE VARIABLE_TYPE
   USE checkdb
   IMPLICIT NONE
   INTEGER, PARAMETER :: PotentialType=5 !SW is potential type 1
   TYPE(StillingerWeberPotential), POINTER :: SW
   
   CONTAINS

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   SUBROUTINE SWSetUp(AL) 
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
      IF (.NOT. ASSOCIATED(AL%Potential%SW)) RETURN
      
      NSpecies=NSpecies_global
      
      SW=>AL%Potential%SW
      Found=.FALSE.
      !----------------Add the potentials----------------------
      
      DO species1=1,NSpecies

         DO species2=species1,NSpecies

            section=((2*NSpecies_global-species1)*(species1-1))/2+species2
            !IF (SW%EnabledPP(section)) THEN
            IF (SW%format==1) THEN
            AL%Potential%MaxPotentialCutoff(PotentialType)= &
                  MAX(AL%Potential%MaxPotentialCutoff(PotentialType),SW%cutoff(section))
            
            
            ELSE
               CALL SWRead(species1,species2) !store maxpotentialcutoff for SW
               AL%Potential%MaxPotentialCutoff(PotentialType)= &
                  MAX(AL%Potential%MaxPotentialCutoff(PotentialType),SW%RangePP(2*section))
               Found=.TRUE.
            !END IF
            END IF

         END DO

      END DO

      AL%Potential%VLBuffer=MAXVAL(AL%Potential%MaxPotentialCutoff)/30._dp !0.1 !angstrom
      
      IF (Found) WRITE(*,*) "SW>> SW potentials are now available"
      
   END SUBROUTINE SWSetUp

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SWRead(species1,species2) !note that species1<=species2
      IMPLICIT NONE
      INTEGER :: species1,species2,r1,r2,npt,section,i,atmno1,atmno2
      REAL(dp) :: x(2),dxi,conv
      CHARACTER(len=100) :: FileName
      CHARACTER :: tmpchar
      
      section=((2*NSpecies_global-species1)*(species1-1))/2+species2
      IF (SW%EnabledPP(section)) THEN
         FileName = SW%TableLocationPP(section)
         WRITE(*,*) TRIM(FileName)
         OPEN(UNIT=UnitTmp1,FILE=TRIM(FileName))
      
        !XXXXXXXXXXXXXXXXXXXXXXXXXXreading pair potential v2 data for pair ij or ik XXXXXXXXXXXXXXXXXXXXXXXXX 
         READ(UnitTmp1,*) npt,x(1),x(2),dxi
         READ(UnitTmp1,*) tmpchar
         
         IF (npt<1000) THEN
            WRITE(*,*) "$Err>> Too few interp points"
            STOP
         END IF
         IF (npt>MaxSWTableSize) THEN
            WRITE(*,*) "$Err>> Too many interp points ... increase MaxLJTableSize"
            STOP
         END IF

         IF (ABS(1._dp/dxi-(x(2)-x(1))/REAL(npt))>1.e-3_dp/dxi) THEN
            WRITE(*,*) "$Err>> Unable to understand dxi",dxi
            STOP
         END IF
         
         SW%SizePP(section)=npt
         SW%RangePP(2*section-1:2*section)=x
         SW%IncrPPinv(section)=dxi

         r1=-3+(section-1)*MaxSWTableSize
         r2=0+(section-1)*MaxSWTableSize
         !conv=SW%ConvFacPP(section)
         DO i=1,npt/4 !table contains 4 data points in each line
            r1=r1+4
            r2=MIN(r2+4,npt+(section-1)*MaxSWTableSize)
            
            READ(UnitTmp1,*) SW%PP(r1:r2)
           ! SW%PP(r1:r2) = SW%PP(r1:r2)* conv
         END DO
         
         IF (MOD(npt,4)/=0) THEN !there must be additional data points left
            READ(UnitTmp1,*) SW%PP(4*(i-1)+1:npt)
            !SW%PP(4*(i-1)+1:npt)=SW%PP(4*(i-1)+1:npt)*conv
         END IF
      
         !XXXXXXXXXXXXXXXXXXXXXXXXXXreading three body potential v3 for pairs ij or ik XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
         READ(UnitTmp1,*) 
         READ(UnitTmp1,*) npt,x(1),x(2),dxi
         READ(UnitTmp1,*) tmpchar
      
      
         IF (npt<1000) THEN
            WRITE(*,*) "$Err>> Too few interp points"
            STOP
         END IF
         IF (npt>MaxSWTableSize) THEN
            WRITE(*,*) "$Err>> Too many interp points ... increase MaxLJTableSize"
            STOP
         END IF

         IF (ABS(1._dp/dxi-(x(2)-x(1))/REAL(npt))>1.e-3_dp/dxi) THEN
            WRITE(*,*) "$Err>> Unable to understand dxi",dxi
            STOP
         END IF
         
         SW%SizeF3(section)=npt
         SW%RangeF3(2*section-1:2*section)=x
         SW%IncrF3inv(section)=dxi

         r1=-3+(section-1)*MaxSWTableSize
         r2=0+(section-1)*MaxSWTableSize
         
         DO i=1,npt/4 !table contains 4 data points in each line
            r1=r1+4
            r2=MIN(r2+4,npt+(section-1)*MaxSWTableSize)
            READ(UnitTmp1,*) SW%F3(r1:r2)
         END DO
         IF (MOD(npt,4)/=0) THEN !there must be additional data points left
            READ(UnitTmp1,*) SW%F3(4*(i-1)+1:npt)
         END IF
         !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX reading values for epssqrt and lam for pairs ij and ik XXXXXXXXXXXXXXXXXXXXXXXX
         READ(UnitTmp1,*) 
         READ(UnitTmp1,*) tmpchar
         READ(UnitTmp1,*) SW%epssqrt(section)
         READ(UnitTmp1,*) SW%lam(section)
        
         CLOSE(UnitTmp1)
      
         atmno1=SpeciesDirectory_global(species1)
         atmno2=SpeciesDirectory_global(species2)
         WRITE(6,*) "SW PP ",TRIM(SpeciesList%AtomicSymbol(atmno1)),"-", &
         TRIM(SpeciesList%AtomicSymbol(atmno2)),":", x
     ! ELSE
      !   LJ%SizePP(section)=101
       !  x(1)=0._dp
       !  x(2)=3._dp
       !  LJ%RangePP(2*section-1:2*section)=x
       !  LJ%IncrPPinv(section)=3._dp/100._dp
       !  r1=1+(section-1)*MaxSWTableSize
       !  r2=101+(section-1)*MaxSWTableSize
       !  LJ%PP(r1:r2)=0._dp
      END IF
   END SUBROUTINE SWRead
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

END MODULE StillingerWeberInitialize
