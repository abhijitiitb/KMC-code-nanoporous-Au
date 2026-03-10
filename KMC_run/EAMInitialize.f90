!For initializing EAM potential - multiple calls possible
!each time a particular files are read for the Atom List

!EAM potential model:  E(i)=EEFn(SUM_ j{DT(r_ij)})+0.5*SUM_j{PP(r_ij)}
! - This code is designed for many-atom species
! - Energy is computed in terms of tables for EEFn, DT and PP
! - Gradient terms are also provided for quickly estimating the forces
! - Note that for more than 1 atomic species it is best to use a full list

MODULE EAMInitialize

   USE VARIABLE_TYPE
   USE checkdb

   IMPLICIT NONE
   INTEGER, PARAMETER :: PotentialType=2 !EAM is potential type 2
   TYPE(EAMPotential), POINTER :: EAM
   LOGICAL :: IsPrint1
   
   CONTAINS

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   SUBROUTINE EAMSetUp(AL,IsPrint) 
      !is to be only called by PotentialPackageInitialize
      !... this means all variables have been properly initialized
      
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      LOGICAL, OPTIONAL :: IsPrint
      LOGICAL :: Found
      
      IF (PRESENT(IsPrint)) THEN
         IsPrint1=IsPrint
      ELSE
         IsPrint1=.TRUE.
      END IF
      
      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(6,*) "Error>> Atomlist is not associated"
         STOP
      END IF
      IF (.NOT. ASSOCIATED(AL%Potential%EAM)) RETURN !it should have been set-up
      
      EAM=>AL%Potential%EAM
      !----------------Add the potentials----------------------
      SELECT CASE (EAM%TableFormat(1:3))
      CASE("VOT")
         Found=ReadVoterFormat(AL)
      CASE("SET")
         Found=ReadSETFLFormat(AL)
      CASE("FUN")
         Found=ReadFUNCFLFormat(AL)
      END SELECT
      AL%Potential%VLBuffer=MAXVAL(AL%Potential%MaxPotentialCutoff)/20._dp       !20. optimized for Ag/Ag(100)
      
      IF (Found .AND. IsPrint1) WRITE(6,*) "EAM>> EAM potentials are now available"
      
   END SUBROUTINE EAMSetUp

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ReadVoterFormat(AL)
      IMPLICIT NONE
      LOGICAL :: ReadVoterFormat
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: species1,species2,NSpecies,section
      
      IF (IsPrint1) WRITE(6,*) "EAM>> Reading Voter potential format ..."
      ReadVoterFormat=.FALSE.
      
      NSpecies=NSpecies_global
      DO species1=1,NSpecies

         IF (EAM%EnabledEE(species1)) THEN
            CALL EAMReadEE(species1)
            ReadVoterFormat=.TRUE.
         END IF

         IF (EAM%EnabledDT(species1)) THEN
            CALL EAMReadDT(species1)
            AL%Potential%MaxPotentialCutoff(PotentialType)= &
               MAX(AL%Potential%MaxPotentialCutoff(PotentialType),EAM%RangeDT(2*species1))
               !EAM is potential type 2
            ReadVoterFormat=.TRUE.
         END IF

         DO species2=species1,NSpecies

            section=((2*NSpecies_global-species1)*(species1-1))/2+species2
            IF (EAM%EnabledPP(section)) THEN
               CALL EAMReadPP(species1,species2) !store maxpotentialcutoff for EAM
               AL%Potential%MaxPotentialCutoff(PotentialType)= &
                  MAX(AL%Potential%MaxPotentialCutoff(PotentialType),EAM%RangePP(2*section))
               ReadVoterFormat=.TRUE.
            END IF

         END DO
      END DO
      
   END FUNCTION ReadVoterFormat
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EAMReadEE(species1)
      IMPLICIT NONE
      INTEGER :: species1,r1,r2,npt,section,i,atmno1,errorstatus
      REAL(dp) :: x(2),dxi,conv
      CHARACTER(len=100) :: FileName
      
      section=species1
      FileName = EAM%TableLocationEE(section)
      OPEN(UNIT=UnitTmp1,FILE=TRIM(FileName),IOSTAT=errorstatus)
      
      READ(UnitTmp1,*) npt,x(1),x(2),dxi
      IF (npt<1000) THEN
         WRITE(*,*) "$Err>> Too few interp points"
         STOP
      END IF
      IF (npt>MaxEAMTableSize) THEN
         WRITE(*,*) "$Err>> Too many interp points ... increase MaxEAMTableSize"
         STOP
      END IF

      IF (ABS(1._dp/dxi-(x(2)-x(1))/REAL(npt))>1.e-3_dp/dxi) THEN
         WRITE(*,*) "$Err>> Unable to understand dxi",dxi
         STOP
      END IF
      
      EAM%SizeEE(section)=npt
      EAM%RangeEE(2*section-1:2*section)=x
      EAM%IncrEEinv(section)=dxi
      

      r1=-3+(section-1)*MaxEAMTableSize
      r2=0+(section-1)*MaxEAMTableSize
      conv=EAM%ConvFacEE(section)
      DO i=1,npt/4 !table contains 4 data points in each line
         r1=r1+4
         r2=MIN(r2+4,npt+(section-1)*MaxEAMTableSize)
         READ(UnitTmp1,*) EAM%EE(r1:r2)
         EAM%EE(r1:r2) = EAM%EE(r1:r2) * conv !Get the correct conversion factor if energy term
      END DO
      
      IF (MOD(npt,4)/=0) THEN !there must be additional data points left
         READ(UnitTmp1,*) EAM%EE(4*(i-1)+1:npt)
         EAM%EE(4*(i-1)+1:npt)=EAM%EE(4*(i-1)+1:npt)*conv
      END IF
      
      CLOSE(UnitTmp1)
      
      atmno1=SpeciesDirectory_global(species1)
      IF (IsPrint1) WRITE(6,*) "Range EAM EE ",TRIM(SpeciesList%AtomicSymbol(atmno1)),":", x
      
   END SUBROUTINE EAMReadEE
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EAMReadDT(species1)
      IMPLICIT NONE
      INTEGER :: species1,r1,r2,npt,section,i,atmno1
      REAL(dp) :: x(2),dxi,conv
      CHARACTER(len=100) :: FileName
      
      section=species1
      FileName = EAM%TableLocationDT(section)
      OPEN(UNIT=UnitTmp1,FILE=TRIM(FileName))
      
      READ(UnitTmp1,*) npt,x(1),x(2),dxi
      IF (npt<1000) THEN
         WRITE(*,*) "$Err>> Too few interp points"
         STOP
      END IF
      IF (npt>MaxEAMTableSize) THEN
         WRITE(*,*) "$Err>> Too many interp points ... increase MaxEAMTableSize"
         STOP
      END IF

      IF (ABS(1._dp/dxi-(x(2)-x(1))/REAL(npt))>1.e-3_dp/dxi) THEN
         WRITE(*,*) "$Err>> Unable to understand dxi",dxi
         STOP
      END IF
      
      EAM%SizeDT(section)=npt
      EAM%RangeDT(2*section-1:2*section)=x
      EAM%IncrDTinv(section)=dxi

      r1=-3+(section-1)*MaxEAMTableSize
      r2=0+(section-1)*MaxEAMTableSize
      conv=EAM%ConvFacDT(section)
      DO i=1,npt/4 !table contains 4 data points in each line
         r1=r1+4
         r2=MIN(r2+4,npt+(section-1)*MaxEAMTableSize)
         READ(UnitTmp1,*) EAM%DT(r1:r2)
         EAM%DT(r1:r2) = EAM%DT(r1:r2) * conv !usually conv is 1._dp
      END DO
      
      IF (MOD(npt,4)/=0) THEN !there must be additional data points left
         READ(UnitTmp1,*) EAM%DT(4*(i-1)+1:npt)
         EAM%DT(4*(i-1)+1:npt)=EAM%DT(4*(i-1)+1:npt)*conv
      END IF
      
      CLOSE(UnitTmp1)
      
      atmno1=SpeciesDirectory_global(species1)
      IF (IsPrint1) WRITE(6,*) "Range EAM DT ",TRIM(SpeciesList%AtomicSymbol(atmno1)),":", x
      
   END SUBROUTINE EAMReadDT
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EAMReadPP(species1,species2) !note that species1<=species2
      IMPLICIT NONE
      INTEGER :: species1,species2,r1,r2,npt,section,i,atmno1,atmno2
      REAL(dp) :: x(2),dxi,conv
      CHARACTER(len=100) :: FileName
      
      section=((2*NSpecies_global-species1)*(species1-1))/2+species2
      FileName = EAM%TableLocationPP(section)
      OPEN(UNIT=UnitTmp1,FILE=TRIM(FileName))
      
      READ(UnitTmp1,*) npt,x(1),x(2),dxi
      IF (npt<1000) THEN
         WRITE(*,*) "$Err>> Too few interp points"
         STOP
      END IF
      IF (npt>MaxEAMTableSize) THEN
         WRITE(*,*) "$Err>> Too many interp points ... increase MaxEAMTableSize"
         STOP
      END IF

      IF (ABS(1._dp/dxi-(x(2)-x(1))/REAL(npt))>1.e-3_dp/dxi) THEN
         WRITE(*,*) "$Err>> Unable to understand dxi",dxi
         STOP
      END IF
      
      EAM%SizePP(section)=npt
      EAM%RangePP(2*section-1:2*section)=x
      EAM%IncrPPinv(section)=dxi

      r1=-3+(section-1)*MaxEAMTableSize
      r2=0+(section-1)*MaxEAMTableSize
      conv=EAM%ConvFacPP(section)
      DO i=1,npt/4 !table contains 4 data points in each line
         r1=r1+4
         r2=MIN(r2+4,npt+(section-1)*MaxEAMTableSize)
         READ(UnitTmp1,*) EAM%PP(r1:r2)
         EAM%PP(r1:r2) = EAM%PP(r1:r2) * conv
      END DO
      
      IF (MOD(npt,4)/=0) THEN !there must be additional data points left
         READ(UnitTmp1,*) EAM%PP(4*(i-1)+1:npt)
         EAM%PP(4*(i-1)+1:npt)=EAM%PP(4*(i-1)+1:npt)*conv
      END IF
      
      CLOSE(UnitTmp1)
      
      atmno1=SpeciesDirectory_global(species1)
      atmno2=SpeciesDirectory_global(species2)
      IF (IsPrint1) WRITE(6,*) "Range EAM PP ",TRIM(SpeciesList%AtomicSymbol(atmno1)),"-", &
        TRIM(SpeciesList%AtomicSymbol(atmno2)),":", x
      !WRITE(6,*) "  ....first few values:",EAM%PP((section-1)*MaxEAMTableSize+1:(section-1)*MaxEAMTableSize+4)
   END SUBROUTINE EAMReadPP
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ReadSETFLFormat(AL)
      IMPLICIT NONE
      LOGICAL :: ReadSETFLFormat
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: species1,species2,r1,r2,section,i,atmno,atmno1
      INTEGER :: nrho,nr,NSpecies
      REAL(dp) :: drho,dr,rcutoff,conv,r
      CHARACTER :: CrapChar
      CHARACTER(len=100) :: FileName
      
      NSpecies=NSpecies_global
      FileName = EAM%TableLocationEE(1) !same file will be used to read other tables
      OPEN(UNIT=UnitTmp1,FILE=TRIM(FileName))
      
      WRITE(6,*) "Reading setfl format ..."
      ReadSETFLFormat=.FALSE.
      READ(UnitTmp1,*) CrapChar
      READ(UnitTmp1,*) CrapChar
      READ(UnitTmp1,*) CrapChar
      READ(UnitTmp1,*) NSpecies
      READ(UnitTmp1,*) nrho,drho,nr,dr,rcutoff
      AL%Potential%MaxPotentialCutoff(PotentialType)= &
         MAX(AL%Potential%MaxPotentialCutoff(PotentialType),rcutoff)
      
      !Read density and embedding functions
      DO species1=1,NSpecies
         section=species1
         READ(UnitTmp1,*) atmno !,atmmass,latticeconst
         atmno1=SpeciesDirectory_global(species1)
         IF (atmno1/=atmno) THEN
            WRITE(6,*) "EAM>> setfl format file has a different numbering for species ..."
            WRITE(6,*) "Expected atomic number: ",atmno1
            WRITE(6,*) "Found atomic number: ",atmno
            STOP
         END IF
         !embedding terms
         EAM%SizeEE(section)=nrho
         EAM%IncrEEinv(section)=1._dp/drho
         r1=-4+(section-1)*MaxEAMTableSize
         r2= 0+(section-1)*MaxEAMTableSize
         conv=EAM%ConvFacEE(section)
         DO i=1,nrho/5
            r1=r1+5
            r2=MIN(r2+5,nrho+(section-1)*MaxEAMTableSize)
            READ(UnitTmp1,*) EAM%EE(r1:r2)
            EAM%EE(r1:r2) = EAM%EE(r1:r2) * conv !Get the conversion factor for energy term
         END DO
         EAM%RangeEE(2*section-1)=drho
         EAM%RangeEE(2*section)=drho*REAL(nrho,dp)
         
         !IF (species1==1) THEN
!         DO i=1,EAM%SizeEE(section)
!            WRITE(*,*) drho*i,EAM%EE(i+(section-1)*MaxEAMTableSize)
!         END DO
         !END IF

         !density terms
         EAM%SizeDT(section)=nr
         EAM%IncrDTinv(section)=1._dp/dr
         r1=-4+(section-1)*MaxEAMTableSize
         r2= 0+(section-1)*MaxEAMTableSize
         conv=EAM%ConvFacDT(section)
         DO i=1,nr/5
            r1=r1+5
            r2=MIN(r2+5,nr+(section-1)*MaxEAMTableSize)
            READ(UnitTmp1,*) EAM%DT(r1:r2)
            EAM%DT(r1:r2) = EAM%DT(r1:r2) * conv !Get the conversion factor for energy term
         END DO
         EAM%RangeDT(2*section-1)=dr
         EAM%RangeDT(2*section)=dr*REAL(nr,dp)
         
      END DO
      
      !Read pair potential terms -- table units are eV*Angstroms
      DO species1=1,NSpecies
         DO species2=species1,NSpecies
            section=((2*NSpecies_global-species1)*(species1-1))/2+species2          
            EAM%SizePP(section)=nr
            EAM%IncrPPinv(section)=1._dp/dr
            r1=-4+(section-1)*MaxEAMTableSize
            r2= 0+(section-1)*MaxEAMTableSize
            conv=EAM%ConvFacPP(section)
            DO i=1,nr/5
               r1=r1+5
               r2=MIN(r2+5,nr+(section-1)*MaxEAMTableSize)
               READ(UnitTmp1,*) EAM%PP(r1:r2)
               !IF (i==1) WRITE(*,*) EAM%PP(r1:r2)
               EAM%PP(r1:r2) = EAM%PP(r1:r2) * conv !Get the conversion factor for energy term
            END DO
            r=0._dp
            DO i=(section-1)*MaxEAMTableSize+1,(section-1)*MaxEAMTableSize+nr
               r=r+dr
               EAM%PP(i)=EAM%PP(i)/r
            !   WRITE(*,*)EAM%PP(i),i
            END DO
         !IF (species1==1 .AND. species2==1) THEN
!         DO i=1,EAM%SizePP(section)
!            WRITE(*,*) dr*i,EAM%PP(i+(section-1)*MaxEAMTableSize)
!         END DO
         !END IF
            EAM%RangePP(2*section-1)=dr
            EAM%RangePP(2*section)=dr*REAL(nr,dp)
         END DO
      END DO
!      stop
      ReadSETFLFormat=.TRUE.
   END FUNCTION ReadSETFLFormat
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ReadFUNCFLFormat(AL)
      IMPLICIT NONE
      LOGICAL :: ReadFUNCFLFormat
      TYPE(SystemContainer), POINTER :: AL
      
      ReadFUNCFLFormat=.FALSE.
      ReadFUNCFLFormat=.TRUE.
   END FUNCTION ReadFUNCFLFormat
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE EAMInitialize
   
