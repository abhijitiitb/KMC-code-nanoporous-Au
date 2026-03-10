MODULE IO

   !reads & prints files
   !standard input in f95: 5,100
   !standard output in f95: 6,101
   !standard error in f95: 0,102

   USE VARIABLE_TYPE
   USE AMD_VARIABLE_TYPE
   USE db_manipulate
   USE utilities
   IMPLICIT NONE
   
   INTERFACE Analyze
      MODULE PROCEDURE AnalyzeAL,AnalyzeChoS
   END INTERFACE Analyze
   INTERFACE PrintInfo
      MODULE PROCEDURE PrintInfo_AtomicNo,PrintInfo_AtomicSymbol
   END INTERFACE PrintInfo
   CONTAINS

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ConvertFile(filein,fileout,form)
      !convert one file format to another

      IMPLICIT NONE
      CHARACTER(len=100) :: filein,fileout1
      CHARACTER(len=100), OPTIONAL :: fileout
      CHARACTER(len=3) :: form
      TYPE(SystemContainer), POINTER :: AL

      IF (PRESENT(fileout)) THEN
         fileout1=fileout
      ELSE
         fileout1=TRIM(filein)//"."//TRIM(form)
      END IF
      
      CALL ReadCoord(AL,filein)

      SELECT CASE(form(1:3))
      CASE('xyz')
         CALL WriteXYZ(AL,fileout1)
      CASE('cls')
         CALL WriteClsman(AL,fileout1)
      CASE('pdb')
      END SELECT

   END SUBROUTINE ConvertFile
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReadCoord(AL,filename,iprint,filetype)
      !reads coordinates from filename into AL
      !if more than one snapshot present in the file then a chain of states is created
      !iprint is print level
      
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      CHARACTER(len=100) :: filename
      INTEGER :: UnitTmp,errorstatus,iprint1
      INTEGER, OPTIONAL :: iprint,filetype
      
      IF (PRESENT(iprint)) THEN
         iprint1=iprint
      ELSE
         iprint1=0
      END IF
      
      UnitTmp=301
      OPEN(UNIT=UnitTmp,FILE=TRIM(filename),IOSTAT=errorstatus)
      IF (errorstatus/=0) THEN
         WRITE(6,*) "$Err: Unable to find file for reading coordinates"
         CALL FLUSH(6)
         STOP
      END IF
      CLOSE(UNIT=UnitTmp)
      
      WRITE(6,*) "IO>> Reading coordinates from file ",TRIM(filename)
      IF (PRESENT(filetype)) THEN
         IF (filetype==1) THEN !xyz
            IF (ReadXYZ(AL=AL,filename=filename,opt=1,iprint=iprint1)) RETURN
            IF (ReadXYZ(AL,filename,opt=2,iprint=iprint1)) RETURN
            IF (ReadXYZ(AL,filename,opt=3,iprint=iprint1)) RETURN
         ELSEIF (filetype==2) THEN !clsman
            IF (ReadClsman(AL,filename,opt=1)) RETURN
            IF (ReadClsman(AL,filename,opt=2)) RETURN
            IF (ReadClsman(AL,filename,opt=3)) RETURN
            IF (ReadClsman(AL,filename,opt=4)) RETURN
            IF (ReadClsman(AL,filename,opt=5)) RETURN
         END IF
      ELSE
         IF (ReadXYZ(AL=AL,filename=filename,opt=1,iprint=iprint1)) RETURN
         IF (ReadXYZ(AL,filename,opt=2,iprint=iprint1)) RETURN
         IF (ReadXYZ(AL,filename,opt=3,iprint=iprint1)) RETURN
         IF (ReadXYZ(AL,filename,opt=4,iprint=iprint1)) RETURN
         IF (ReadXYZ(AL,filename,opt=5,iprint=iprint1)) RETURN
         IF (ReadClsman(AL,filename,opt=1)) RETURN
         IF (ReadClsman(AL,filename,opt=2)) RETURN
         IF (ReadClsman(AL,filename,opt=3)) RETURN
         IF (ReadClsman(AL,filename,opt=4)) RETURN
         IF (ReadClsman(AL,filename,opt=5)) RETURN
      END IF
      WRITE(6,*) "$Err>> Unable to read/recognize input file format: ",TRIM(filename)
      CALL FLUSH(6)
      STOP

   END SUBROUTINE ReadCoord
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ReadClsman(AL,filename,opt)
      !opt=1: Check for Abhijit CLSMAN - when line2 is comment line and no velocities given
      !opt=2: Check for AFV CLSMAN - when line2 is not comment line and velocities given
      !opt=3: Check for AFV CLSMAN - when line2 is not comment line and no velocities given
      !opt=4: Check for AFV CLSMAN with user defn format- when line2 is not comment line and no velocities given

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1
      CHARACTER(len=100) :: filename
      INTEGER :: UnitTmp,errorstatus
      CHARACTER(len=30) :: form,Ctitle
      CHARACTER :: CrapChar
      LOGICAL :: ReadClsman,ReadClsmanPrevImage
      INTEGER :: i,r1,r2,NAtoms,iop,NSpecies,species,NMove(3),nimages
      INTEGER, INTENT(IN) :: opt
      INTEGER, DIMENSION(200) :: Counter
      REAL(dp) :: Box(3)
      
      ReadClsman=.FALSE.
      ReadClsmanPrevImage=.FALSE. 
      
      NMove=0
      UnitTmp=301 !GetUnitNotInUse()
      OPEN(UNIT=UnitTmp,FILE=TRIM(filename))
      SELECT CASE(opt)
      CASE(1); READ(UNIT=UnitTmp,FMT=*,IOSTAT=errorstatus) NAtoms,NMove
      CASE(2,3,4); READ(UNIT=UnitTmp,FMT=*,IOSTAT=errorstatus) NAtoms,NMove,iop
      CASE(5); READ(UNIT=UnitTmp,FMT=*,IOSTAT=errorstatus) NAtoms
      END SELECT
      INCLUDE "io_readclsman_step1.f90"
      
      IF (.NOT. ASSOCIATED(AL)) CALL MakeSize(AL,NAtoms)
      AL%NAtoms=NAtoms
      AL1=>AL
      nimages=0
      
      ReadClsmanPrevImage=.FALSE. !this tracks if previous image could be read
      DO WHILE (ReadClsman)
         
         SELECT CASE(opt)
         CASE(1)
            READ(UNIT=UnitTmp,FMT=*,IOSTAT=errorstatus) CrapChar
         CASE(2,3,4)
            READ(UNIT=UnitTmp,FMT='(a80)',IOSTAT=errorstatus) Ctitle(1:1)
         END SELECT
      
         INCLUDE "io_readclsman_step1.f90"  
      
         READ(UNIT=UnitTmp,FMT='(3f20.8)',IOSTAT=errorstatus) Box
         INCLUDE "io_readclsman_step1.f90"
      
         IF (opt==4) THEN
            READ(UNIT=UnitTmp,FMT='(a)',IOSTAT=errorstatus) form
            INCLUDE "io_readclsman_step1.f90"
         END IF
         
         r1=-2; r2=0
         Counter=0
         
         DO i=1,NAtoms
            r1=r1+3; r2=r2+3
            SELECT CASE(opt)
            CASE(2)
               READ(UNIT=UnitTmp,FMT='(3f20.8,i5,3f20.8)',&
                  IOSTAT=errorstatus) &
                  AL1%AtomCoord(r1:r2),AL1%AtomSpecies(i), &
                  AL1%AtomVelocity(r1:r2)
            CASE(1,3,5)
               READ(UNIT=UnitTmp,FMT='(3f20.8,i5)',&
                  IOSTAT=errorstatus) &
                  AL1%AtomCoord(r1:r2),AL1%AtomSpecies(i)
            CASE(4)
               READ(UNIT=UnitTmp,FMT=TRIM(form),&
                  IOSTAT=errorstatus) &
                  AL1%AtomCoord(r1:r2),AL1%AtomSpecies(i)
            END SELECT
            INCLUDE "io_readclsman_step1.f90"
            Counter(AL1%AtomSpecies(i))=Counter(AL1%AtomSpecies(i))+1
         END DO

         nimages=nimages+1
         IF (ALL(NMove==0)) NMove=NAtoms
         AL1%NMove=NMove
      
         AL1%BoxSize(1:3)=Box
         AL1%BoxSize(4:6)=0._dp
         AL1%AtomIsMoving=.FALSE.
         AL1%AtomIsMoving(1:3*NMove(1)-2:3)=.TRUE.
         AL1%AtomIsMoving(2:3*NMove(2)-1:3)=.TRUE.
         AL1%AtomIsMoving(3:3*NMove(3)  :3)=.TRUE.
         
         IF (nimages==1) THEN
            WRITE(UNIT=6,FMT='(" ...successfully read ",I5," atoms from file")') AL1%NAtoms
            WRITE(UNIT=6,FMT='(" ...number of moving atoms ",3I5," atoms from file")') NMove
            WRITE(UNIT=6,FMT='(" ...identified as clsman format -- type",I2)') opt
            WRITE(UNIT=6,FMT='(" ...number of non-zero momenta, x=[",I4,"],y=[",I4,"],z=[",I4,"]")') NMove
            
            NSpecies=0
            WRITE(UNIT=6,FMT='(" ...species present in system: ")',ADVANCE="NO")
            DO i=1,SIZE(Counter)
               IF (Counter(i)>0) THEN
                  NSpecies=NSpecies+1
                  species=SpeciesDirectory_global(i)
                  WRITE(UNIT=6,FMT='(a3," ")',ADVANCE="NO") &
                     SpeciesList%AtomicSymbol(SpeciesDirectory_global(i))
               END IF
            END DO
            WRITE(6,*) ""
            WRITE(UNIT=6,FMT='(" ...detected ",I2," species")') NSpecies
         END IF
            
         CALL GetAtomListMass(AL1) !add atom mass
         
         !prepare for next read
         !so far things have gone successfully
         ReadClsmanPrevImage=.TRUE.
         SELECT CASE(opt)
         CASE(1)
            READ(UNIT=UnitTmp,FMT=*,IOSTAT=errorstatus) NAtoms,NMove
         CASE(2,3,4)
            READ(UNIT=UnitTmp,FMT=*,IOSTAT=errorstatus) NAtoms,NMove,iop
         CASE(5)
            READ(UNIT=UnitTmp,FMT=*,IOSTAT=errorstatus) NAtoms
         END SELECT
         INCLUDE "io_readclsman_step1.f90"   
      
         IF (.NOT. ASSOCIATED(AL1%NextNeigh)) CALL MakeSize(AL1%NextNeigh,NAtoms)
         AL1%NextNeigh%NAtoms=NAtoms
         AL1=>AL1%NextNeigh
         
      END DO 
      
      WRITE(UNIT=6,FMT='(" ...number of images read ",I4)') nimages
      
      CLOSE(UnitTmp)

   END FUNCTION ReadClsman
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ReadXYZ(AL,filename,opt,iframe_min,iframe_max,nstrides, &
      nframes,fileopen,fileclose,iunit,iprint)
   !Reads an XYZ file given by filename
   !NOTE THAT AT LEAST ONE FRAME WILL BE READ IF THE FILE IS READABLE
   !Optionally the user can open the file outside and pass the file unit number iunit with fileopen and fileclose=.FALSE.
   !iframe_min = read from this frame onwards (if iframe_min it too large, then the code will read only the last frame)
   !iframe_max = read upto this frame
   !nstrides = number of strides after which frame is added
   !nframe = number of frames found in the xyz file (OUTPUT)
   !iprint = print level (0- print all, 1 - print minimum, 2 - print none)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1
      CHARACTER(len=100), OPTIONAL :: filename
      CHARACTER(len=3) :: AtomicSymbol
      INTEGER :: UnitTmp,errorstatus,nimages
      LOGICAL :: ReadXYZ,ReadXYZPrevImage,fileopen1,fileclose1
      INTEGER, DIMENSION(200) :: Counter !in reality we dont need such a big array - but this avoids deallocation 
      ! this is useful when the subroutine returns abruptly
      INTEGER :: i,NSpecies,species,NAtoms,r1,r2,NMove(3),iframe_min1,iframe_max1,iframe
      INTEGER, INTENT(IN) :: opt !provides the reading options
      REAL(dp), DIMENSION(6) :: Box
      INTEGER, OPTIONAL :: iframe_min,iframe_max,iunit,iprint,nframes,nstrides
      INTEGER :: iprint1,nstrides1
      LOGICAL, OPTIONAL :: fileopen,fileclose
      
      ReadXYZ=.FALSE.
      
      IF (PRESENT(iprint)) THEN
         iprint1=iprint
      ELSE
         iprint1=0
      END IF

      IF (PRESENT(fileopen)) THEN
         fileopen1=fileopen
      ELSE
         fileopen1=.TRUE.
      END IF

      IF (fileopen1) THEN 
         UnitTmp=301 !GetUnitNotInUse()
         IF (.NOT. PRESENT(filename)) THEN
            WRITE(6,*) "$Err>> Filename should have been provided in ReadXYZ"
            STOP
         END IF
         OPEN(UNIT=UnitTmp,FILE=TRIM(filename))
      ELSE
         IF (PRESENT(iunit)) THEN
            UnitTmp=iunit
         ELSE
            WRITE(6,*) "$Err>> File unit number should have been specified in ReadXYZ"
            STOP
         END IF
      END IF

      IF (PRESENT(iframe_min)) THEN
         iframe_min1=iframe_min
      ELSE
         iframe_min1=1 !exceeding large number to ensure all frames/images are read
      END IF

      IF (PRESENT(iframe_max)) THEN
         iframe_max1=iframe_max
      ELSE
         iframe_max1=100000000 !exceedingly large number to ensure all frames/images are read
      END IF
      
      IF (PRESENT(nstrides)) THEN
         nstrides1=nstrides !frames that are given by 1 + multiples of 10
      ELSE
         nstrides1=1 !each frame will be added
      END IF
      
      NAtoms=0
      READ(UNIT=UnitTmp,FMT=*,IOSTAT=errorstatus) NAtoms
      INCLUDE "io_readxyz_step1.f90"
      READ(UNIT=UnitTmp,FMT=*,IOSTAT=errorstatus) Box(1:3)
      Box(4:6)=0
      INCLUDE "io_readxyz_step1.f90"
      
      IF (.NOT. ASSOCIATED(AL)) CALL MakeSize(AL,NAtoms)
      AL%NAtoms=NAtoms
      
      AL1=>AL
      iframe=0
      nimages=1 !number of images stored in the file
      
      ReadXYZPrevImage=.FALSE.
      DO WHILE (ReadXYZ)
      
         NMove=0
         Counter=0
      
         r1=-2
         r2=0
         DO i=1,NAtoms
            r1=r1+3
            r2=r2+3
            SELECT CASE(opt)
            CASE(1)
               READ(UNIT=UnitTmp,FMT=*, &
                  IOSTAT=errorstatus) &
                  AtomicSymbol,AL1%AtomCoord(r1:r2),AL1%AtomIsMoving(r1:r2),AL1%AtomCharge(i), &
                  AL1%AtomVelocity(r1:r2)
            CASE(2)
               READ(UNIT=UnitTmp,FMT=*, &
                  IOSTAT=errorstatus) &
                  AtomicSymbol,AL1%AtomCoord(r1:r2),AL1%AtomIsMoving(r1:r2),AL1%AtomVelocity(r1:r2)
            CASE(3)
               READ(UNIT=UnitTmp,FMT=*,&
                  IOSTAT=errorstatus) &
                  AtomicSymbol,AL1%AtomCoord(r1:r2),AL1%AtomIsMoving(r1:r2),AL1%AtomCharge(i)
            CASE(4)
               READ(UNIT=UnitTmp,FMT=*,&
                  IOSTAT=errorstatus) &
                  AtomicSymbol,AL1%AtomCoord(r1:r2),AL1%AtomIsMoving(r1:r2)
            CASE(5)
               READ(UNIT=UnitTmp,FMT=*,&
                  IOSTAT=errorstatus) &
                  AtomicSymbol,AL1%AtomCoord(r1:r2)
               AL1%AtomIsMoving(r1:r2)=.TRUE.
            END SELECT
            INCLUDE "io_readxyz_step1.f90"
         
            species=GetSpeciesSymbol(AtomicSymbol)
            AL1%AtomSpecies(i)=species
            IF (species<1 .OR. species>NSpecies_global) THEN
               WRITE(6,*) "Err>> Species type is not understood ... check input.pot file"
               STOP
            END IF
            Counter(AL1%AtomSpecies(i))=Counter(AL1%AtomSpecies(i))+1
            WHERE (AL1%AtomIsMoving(r1:r2)) NMove=NMove+1
         END DO
         
         AL1%BoxSize=Box
         AL1%NMove=NMove
         
         iframe=iframe+1 !incremented each time
         
         IF (iframe==1) THEN
            IF (iprint1>0) WRITE(UNIT=*,FMT='(" ...succesfully read ",I7," atoms from file")') AL1%NAtoms
            IF (iprint1>0) WRITE(UNIT=6,FMT='(" ...identified as XYZ format -- type",I2)') opt
            IF (iprint1>0) WRITE(UNIT=6,FMT='(" ...number of non-zero momenta, x=[",I7,"],y=[",I7,"],z=[",I7,"]")') NMove
      
            NSpecies=0
            IF (iprint1>0) WRITE(UNIT=6,FMT='(" ...species present in system: ")',ADVANCE="NO")
            DO i=1,SIZE(Counter)
               IF (Counter(i)>0) THEN
                  NSpecies=NSpecies+1
                  species=SpeciesDirectory_global(i)
                  IF (iprint1>0) WRITE(UNIT=6,FMT='(a3," ")',ADVANCE="NO") &
                     SpeciesList%AtomicSymbol(SpeciesDirectory_global(i))
               END IF
            END DO
            IF (iprint1==0) WRITE(6,*) ""
            IF (iprint1>0) WRITE(UNIT=6,FMT='(" ...detected ",I2," species")') NSpecies
         END IF
         
         !add atom mass
         CALL GetAtomListMass(AL1)
       
         !if the number of images read equals the number of frames allowed to read then exit
         IF (iframe==iframe_max1) EXIT
         
         !prepare for next image read
         ReadXYZPrevImage=.TRUE.
         NAtoms=0
         READ(UNIT=UnitTmp,FMT=*,IOSTAT=errorstatus) NAtoms
         INCLUDE "io_readxyz_step1.f90"
      
         READ(UNIT=UnitTmp,FMT=*,IOSTAT=errorstatus) Box
         INCLUDE "io_readxyz_step1.f90"
         
         IF (.NOT. ASSOCIATED(AL1%NextNeigh) .AND. iframe>=iframe_min1 .AND. iframe<iframe_max1 .AND. &
           MOD(iframe-1,nstrides1)==0) THEN
            nimages=nimages+1 !one more image created in AL
            CALL MakeSize(AL1%NextNeigh,NAtoms)
            AL1%NextNeigh%NAtoms=NAtoms
            AL1=>AL1%NextNeigh
         END IF
      END DO
      
      IF (iprint1<2) WRITE(UNIT=6,FMT='("...number of images read ",I4," between frames",I4,"-",I4)') &
         nimages,MIN(iframe_min1,iframe),MIN(iframe_max1,iframe)

      IF (PRESENT(fileclose)) THEN
         fileclose1=fileclose
      ELSE
         fileclose1=.TRUE.
      END IF

      IF (fileclose1) CLOSE(UnitTmp)
      
      IF (PRESENT(nframes)) nframes=iframe
      
   END FUNCTION ReadXYZ
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReadLAMMPSCoord(AL,NAtoms,NMove,BoxSize,filename)
   !reads a LAMMPS input coordinate file with NAtoms present
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      CHARACTER(len=100) :: filename
      INTEGER :: iatom,id,ispecies,q,NAtoms
      REAL(dp) :: x,y,z
      INTEGER, OPTIONAL :: NMove(3)
      REAL(dp), OPTIONAL :: BoxSize(6)
      
      CALL MakeSize(AL,NAtoms)
      AL%NAtoms=NAtoms
      IF (PRESENT(NMove)) THEN
         AL%NMove=NMove
         AL%AtomIsMoving=.FALSE.
         DO iatom=1,NMove(1)
            AL%AtomIsMoving(3*iatom-2)=.TRUE.
         END DO
         DO iatom=1,NMove(2)
            AL%AtomIsMoving(3*iatom-1)=.TRUE.
         END DO
         DO iatom=1,NMove(3)
            AL%AtomIsMoving(3*iatom)=.TRUE.
         END DO
      ELSE
         AL%NMove=NAtoms
         AL%AtomIsMoving(1:3*NAtoms)=.TRUE.
         WRITE(6,*) "Warning>> all atoms can move"
      END IF
      
      IF (PRESENT(BoxSize)) THEN
         AL%BoxSize=BoxSize
      ELSE
         AL%BoxSize=0._dp
         WRITE(6,*) "Warning>> BoxSize is set to zero"
      END IF
      OPEN(UNIT=301,FILE=filename)
      DO iatom=1,NAtoms
         READ(301,*) id,ispecies,q,x,y,z
         AL%AtomCoord(3*id-2)=x
         AL%AtomCoord(3*id-1)=y
         AL%AtomCoord(3*id)=z
         AL%AtomCharge(id)=REAL(q,dp)
         IF (ispecies==0) THEN
            WRITE(6,*) "$Err>> Atom species type found to be zero for id ",id
            STOP
         END IF
         AL%AtomSpecies(id)=ispecies
      END DO
      CLOSE(301)
      WRITE(6,*) ">>> Successfully read the LAMMPS atoms coordinates"
   END SUBROUTINE ReadLAMMPSCoord
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReadPotentialInfo(AL)

   !Reads input.pot -- information related to interaction potentials

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(InteracPotential), POINTER :: Potential
      TYPE(EwaldParameters), POINTER :: Ewald
      INTEGER :: nentries,iunit1,errorstatus,entryno,s1,s2,s3,smin,smax,section,speciesarr(10)
      INTEGER :: atmno1,atmno2,atmno3,itmno1,itmno2,ispecies,nspecies,nspecies1,atmnoarr(10)
      REAL(dp) :: conv,value
      CHARACTER(len=10) :: PotentialName,DataType
      CHARACTER(len=100) :: TableLocation
      CHARACTER(len=2) :: term
      CHARACTER(len=100) :: FileName
      
      Potential=>AL%Potential
      
      FileName="input.pot"
      iunit1=UnitTmp1
      OPEN(UNIT=iunit1,FILE=TRIM(FileName),IOSTAT=errorstatus)
      IF (errorstatus/=0) THEN
         WRITE(*,*) "$Err>> Unable to open ",TRIM(FileName)
         STOP
      END IF

      READ(iunit1,*) nentries,nspecies
      !READ(UNIT=iunit1,FMT='(2I4)',IOSTAT=errorstatus) nentries,nspecies
      IF (errorstatus/=0) THEN
         WRITE(*,*) "$Err>> Unable to read # entries in input.pot"
         STOP
      END IF
      
      IF (NSpecies_global==0) THEN
         WRITE(6,*) "$Err>> Number of species NSpecies_global is zero"
         STOP
      END IF
      IF (NSpecies_global/=nspecies) THEN
         WRITE(6,*) "$Err>> Number of species in input.pot doesnt match NSpecies_global"
         STOP
      END IF
      
      DO entryno=1,nentries
         READ(UNIT=iunit1,FMT='(a10)',IOSTAT=errorstatus) PotentialName
         IF (errorstatus/=0) THEN
            WRITE(*,*) "$Err>> Unable to read PotentialName in input.pot at entry ",entryno
            STOP
         END IF
         SELECT CASE(PotentialName(1:2))
         CASE("LJ") !LJ
            IF (.NOT. ASSOCIATED(Potential%LJ)) ALLOCATE(Potential%LJ)
            READ(iunit1,*) s1,s2,atmno1,atmno2  !convert this to a more general reader which can read nspecies instead of specifically 2
            !itmno1=GetSpecies(atmno1)
            !itmno2=GetSpecies(atmno2)
            smin=MIN(s1,s2)
            smax=MAX(s1,s2)
            section=((2*NSpecies_global-smin)*(smin-1))/2+smax
            READ(UNIT=iunit1,FMT='(a100)') Potential%LJ%TableLocationPP(section)
            READ(iunit1,*) Potential%LJ%ConvFacPP(section)
            Potential%LJ%EnabledPP(section)=.TRUE.
         CASE("EA") !EAM
            IF (.NOT. ASSOCIATED(Potential%EAM)) ALLOCATE(Potential%EAM)
            IF (PotentialName(6:8)/="VOT" .AND. PotentialName(6:8)/="SET" .AND. PotentialName(6:8)/="FUN") THEN
               WRITE(6,*) "Err>> Provide EAM table format type in input.pot, e.g., VOTER, SETFL, FUNCFL"
               STOP
            END IF
            IF (Potential%EAM%TableFormat(1:3)=="xxx") THEN
               Potential%EAM%TableFormat=PotentialName(6:8)
            ELSE
               IF (Potential%EAM%TableFormat(1:3)/=PotentialName(6:8)) THEN
                  WRITE(6,*) "Err>> EAM type already provided in input.pot as ", &
                     Potential%EAM%TableFormat(1:3),"format"
                  WRITE(6,*) "use only one format"
                  STOP
               END IF
            END IF
            IF (PotentialName(1:5)=="EAMPP") THEN
               READ(iunit1,*) s1,s2,atmno1,atmno2
               !itmno1=GetSpecies(atmno1)
               !itmno2=GetSpecies(atmno2)
               smin=MIN(s1,s2)
               smax=MAX(s1,s2)
               section=((2*NSpecies_global-smin)*(smin-1))/2+smax
               READ(UNIT=iunit1,FMT='(a100)') Potential%EAM%TableLocationPP(section)
               READ(iunit1,*) Potential%EAM%ConvFacPP(section)
               Potential%EAM%EnabledPP(section)=.TRUE.
            ELSEIF (PotentialName(1:5)=="EAMDT") THEN !write as EAMDT in input.pot
               READ(iunit1,*) s1,atmno1
               !itmno1=GetSpecies(atmno1)
               section=s1
               READ(UNIT=iunit1,FMT='(a100)') Potential%EAM%TableLocationDT(section)
               READ(iunit1,*) Potential%EAM%ConvFacDT(section)
               Potential%EAM%EnabledDT(section)=.TRUE.
            ELSEIF (PotentialName(1:5)=="EAMEE") THEN
               READ(iunit1,*) s1,atmno1
               !itmno1=GetSpecies(atmno1)
               section=s1
               READ(UNIT=iunit1,FMT='(a100)') Potential%EAM%TableLocationEE(section)
               READ(iunit1,*) Potential%EAM%ConvFacEE(section)
               Potential%EAM%EnabledEE(section)=.TRUE.
            ELSEIF (PotentialName(1:5)=="EAM " .AND. (Potential%EAM%TableFormat(1:3)=="SET" &
               .OR. Potential%EAM%TableFormat(1:3)=="FUN")) THEN !for SETFL and FUNFL there is only one file
               READ(iunit1,*) speciesarr(1:nspecies),atmnoarr(1:nspecies)
               READ(UNIT=iunit1,FMT='(a100)') TableLocation
               READ(iunit1,*) conv
               DO s1=1,nspecies
                  section=s1
                  Potential%EAM%TableLocationEE(section)=TableLocation
                  Potential%EAM%ConvFacEE(section)=conv
                  WRITE(6,*) "EE:",s1,TRIM(Potential%EAM%TableLocationEE(section))
                  Potential%EAM%TableLocationDT(section)=TableLocation
                  Potential%EAM%ConvFacDT(section)=conv
                  WRITE(6,*) "DT:",s1,TRIM(Potential%EAM%TableLocationDT(section))
                  DO s2=s1,nspecies
                     smin=MIN(s1,s2)
                     smax=MAX(s1,s2)
                     section=((2*NSpecies_global-smin)*(smin-1))/2+smax
                     Potential%EAM%TableLocationPP(section)=TableLocation
                     Potential%EAM%ConvFacPP(section)=conv
                     WRITE(6,*) "PP:",s1,s2,TRIM(Potential%EAM%TableLocationPP(section))
                  END DO
               END DO
            ELSE
               WRITE(6,*) "Unable to comprehend EAM potential information"
               WRITE(UNIT=6,FMT='("Potential entry number ",I3," in file ")',ADVANCE="NO") entryno
               WRITE(6,*) TRIM(FileName)
               STOP
            END IF
         CASE("CO") !COULOMB
            IF (.NOT. ASSOCIATED(Potential%Coulomb)) ALLOCATE(Potential%Coulomb)
            READ(iunit1,*) nspecies
            DO ispecies=1,nspecies
               READ(iunit1,*) s1,atmno1
            END DO
            IF (PotentialName(1:9)=="COULOMBEW") THEN
               IF (.NOT. ASSOCIATED(Potential%Coulomb%Ewald)) ALLOCATE(Potential%Coulomb%Ewald)
               Ewald=>Potential%Coulomb%Ewald
               READ(iunit1,*) Ewald%ALPHA,Ewald%TimeRatio,Ewald%ReCutOff,Ewald%ImCutOff,Ewald%Precision
               WRITE(UNIT=6,FMT='("Initialized Ewald parameters: ")')
               IF (Ewald%ALPHA>0._dp) WRITE(6,*) "Alpha:",Ewald%ALPHA
               IF (Ewald%TimeRatio>0._dp)WRITE(6,*) "Time ratio:",Ewald%TimeRatio
               IF (Ewald%ReCutOff>0._dp) WRITE(6,*) "Real space cutoff:",Ewald%ReCutOff
               IF (Ewald%ImCutOff>0) WRITE(6,*) "Imaginary space cutoff:",Ewald%ImCutOff
               IF (Ewald%Precision>0._dp) WRITE(6,*) "Precision:",Ewald%Precision
            END IF
            IF (PotentialName(1:9)=="COULOMBWO") THEN
               IF (.NOT. ASSOCIATED(Potential%Coulomb%Wolf)) ALLOCATE(Potential%Coulomb%Wolf)
               READ(iunit1,*) Potential%Coulomb%Wolf%ReCutOff
               READ(iunit1,*) Potential%Coulomb%Wolf%ALPHA
            END IF
         CASE("BU") !BUCKINGHAM
            IF (.NOT. ASSOCIATED(Potential%Buckingham)) ALLOCATE(Potential%Buckingham)
            READ(iunit1,*) s1,s2,atmno1,atmno2
            !itmno1=GetSpecies(atmno1)
            !itmno2=GetSpecies(atmno2)
            smin=MIN(s1,s2)
            smax=MAX(s1,s2)
            section=((2*NSpecies_global-smin)*(smin-1))/2+smax
            READ(UNIT=iunit1,FMT='(a100)') Potential%Buckingham%TableLocationPP(section)
            WRITE(6,*) TRIM(Potential%Buckingham%TableLocationPP(section))
            READ(iunit1,*) Potential%Buckingham%ConvFacPP(section)
            Potential%Buckingham %EnabledPP(section)=.TRUE.
         CASE("SW") !STILLINGER-WEBER
            IF (.NOT. ASSOCIATED(Potential%SW)) ALLOCATE(Potential%SW)
            READ(iunit1,*) s1,s2,atmno1,atmno2  !convert this to a more general reader which can read nspecies instead of specifically 2
            !itmno1=GetSpecies(atmno1)
            !itmno2=GetSpecies(atmno2)
            smin=MIN(s1,s2)
            smax=MAX(s1,s2)
            section=((2*NSpecies_global-smin)*(smin-1))/2+smax
            IF (PotentialName(5:8)=="ANAL") THEN
               READ(iunit1,FMT=*) Potential%SW%A(section),Potential%SW%B(section),Potential%SW%g(section), &
                  Potential%SW%q(section),Potential%SW%epsln(section),Potential%SW%sigma(section), &
                  Potential%SW%ac(section),Potential%SW%lmd(section),Potential%SW%gamma(section), &
                  Potential%SW%cutoff(section)
               Potential%SW%format=1
            ELSE
               READ(UNIT=iunit1,FMT='(a100)') Potential%SW%TableLocationPP(section)
               Potential%SW%format=2
            END IF
            Potential%SW%EnabledPP(section)=.TRUE.
            !WRITE(UNIT=6,FMT='("File",I2,":")',ADVANCE="NO") section
            !WRITE(UNIT=6,FMT=*) Potential%SW%TableLocationPP(section)
         CASE("TE") !Tersoff potential
            IF (.NOT. ASSOCIATED(Potential%Tersoff)) ALLOCATE(Potential%Tersoff)
            READ(iunit1,*) s1,s2,atmno1,atmno2
            smin=MIN(s1,s2)
            smax=MAX(s1,s2)
            section=((2*NSpecies_global-smin)*(smin-1))/2+smax
            READ(iunit1,*) Potential%Tersoff%A(section),Potential%Tersoff%B(section), &
               Potential%Tersoff%lam1(section),Potential%Tersoff%lam2(section), &
               Potential%Tersoff%lam3(section),Potential%Tersoff%beta(section), &
               Potential%Tersoff%n(section),Potential%Tersoff%c(section), &
               Potential%Tersoff%d(section),Potential%Tersoff%h(section), &
               Potential%Tersoff%R(section),Potential%Tersoff%RD(section)
         CASE("FE") !FENE
         CASE("AI") !Airebo - bond order
         CASE("ME") !MEAM
            IF (.NOT. ASSOCIATED(Potential%Tersoff)) ALLOCATE(Potential%Tersoff)
            
            IF (PotentialName(5:8)=="ANAL") THEN

            READ(UNIT=iunit1,FMT='(a10)',ADVANCE="NO") DataType !PURE or MIX or CMIN or CMAX
            IF (DataType(1:3)=="PUR") THEN
               READ(iunit1,*) s1,atmno1, &
                  Potential%MEAM%A(s1),Potential%MEAM%beta(4*(s1-1):4*s1), &
                  Potential%MEAM%t(3*(s1-1):3*s1)
            ELSEIF (DataType(1:3)=="MIX") THEN
               READ(UNIT=iunit1,FMT='(4I4)',ADVANCE="NO") s1,s2,atmno1,atmno2
               smin=MIN(s1,s2)
               smax=MAX(s1,s2)
               section=((2*NSpecies_global-smin)*(smin-1))/2+smax
               READ(iunit1,*) Potential%MEAM%Ec(section),Potential%MEAM%Re(section), &
                  Potential%MEAM%alfa(section),Potential%MEAM%d(section)
            ELSEIF (DataType(1:4)=="CMIN") THEN
               READ(iunit1,*) s1,s2,s3,atmno1,atmno2,atmno3,value
               !s2 can be any one of the NSpecies_global, so we identify the section belonging to s1-s3
               smin=MIN(s1,s3)
               smax=MAX(s1,s3)
               section=((2*NSpecies_global-smin)*(smin-1))/2+smax
               Potential%MEAM%Cmin(NSpecies_global*(section-1)+s2)=value
            ELSEIF (DataType(1:4)=="CMAX") THEN
               READ(iunit1,*) s1,s2,s3,atmno1,atmno2,atmno3,value
               !s2 can be any one of the NSpecies_global, so we identify the section belonging to s1-s3
               smin=MIN(s1,s3)
               smax=MAX(s1,s3)
               section=((2*NSpecies_global-smin)*(smin-1))/2+smax
               Potential%MEAM%Cmax(NSpecies_global*(section-1)+s2)=value
            ELSE
               WRITE(6,*) "Err>> Unable to understand input.pot"
               STOP
            END IF

            END IF

            IF (PotentialName(5:8)=="TBLE") THEN !parameters to be read
!               READ(iunit1,*) Re
            END IF
         CASE("TB") !Tight-binding
         CASE("OP") !OPLS potential
            IF (.NOT. ASSOCIATED(Potential%OPLS)) ALLOCATE(Potential%OPLS)
         CASE("DF") !DFT calculation
            IF (.NOT. ASSOCIATED(Potential%DFT)) ALLOCATE(Potential%DFT)
            READ(iunit1,*) s1,s2,atmno1,atmno2
         END SELECT
      END DO
      CLOSE(iunit1)
   END SUBROUTINE ReadPotentialInfo
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WriteChOS(History,filename)
   !Print the position of atoms in a chain of states to generate a xyz trajectory file
      TYPE(SystemContainer), POINTER :: History,AL
      CHARACTER(len=100) :: filename
   
      IF (ASSOCIATED(History)) THEN
         AL=>History
         OPEN(UNIT=301,FILE=TRIM(filename)) !do not change unit from 301 -- see SUBROUTINE WriteXYZ
         DO WHILE (ASSOCIATED(AL))
            CALL WriteXYZ(AL,filename,fileopen=.FALSE.,fileclose=.FALSE.,iunit=301,iprint=0)
            AL=>AL%NextNeigh
         END DO
         CLOSE(301)
      ELSE
         WRITE(6,*) "$Err>> No history associated with pointer"
         STOP
      END IF
      
      WRITE(*,*) "IO>> Written ChOS list to ",TRIM(filename)
   END SUBROUTINE WriteChOS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WriteClsman(AL,filename)

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      CHARACTER(len=100) :: filename
      INTEGER :: iunit,i
      
      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(6,*) "$Err>> AL is not associated in WriteClsman"
         STOP
      END IF
      
      iunit=301
      OPEN(UNIT=iunit,FILE=TRIM(filename))
      WRITE(UNIT=iunit,FMT='(4i5)') AL%NAtoms,AL%NMove
      WRITE(iunit,*) "@@@@@"
      WRITE(UNIT=iunit,FMT='(3f20.8)') AL%BoxSize(1:3)
      
      !write moving atoms first
      DO i=1,AL%NAtoms
         IF (.NOT. ALL (AL%AtomIsMoving(3*i-2:3*i).eqv.AL%AtomIsMoving(3*i))) THEN
            WRITE(*,*) "$Err>> Some of the dimn in WriteClsman are moving, others fixed"
            STOP
         END IF
         IF (AL%AtomIsMoving(3*i)) THEN
            WRITE(UNIT=iunit,FMT='(3f20.8,i5)')  &
              AL%AtomCoord(3*i-2:3*i),AL%AtomSpecies(i)
         END IF
         IF (AL%AtomSpecies(i)==0) THEN
            WRITE(6,*) "$Err>> Atom species was found to be zero in WriteClsman"
            STOP
         END IF
      END DO
      !write fixed atoms next
      DO i=1,AL%NAtoms
         IF (.NOT. AL%AtomIsMoving(3*i)) THEN
            WRITE(UNIT=iunit,FMT='(3f20.8,i5)')  &
              AL%AtomCoord(3*i-2:3*i),AL%AtomSpecies(i)
         END IF
      END DO
      CLOSE(iunit)
      
      WRITE(6,*) "IO>> Written atom list to ",TRIM(filename)
   
   END SUBROUTINE WriteClsman
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WriteXYZ(AL,filename,fileopen,fileclose,iunit,iprint,NAtoms,filtermode)
   !fileopen and fileclose has been added for WriteChOS
   !fileopen=.TRUE. asks WriteXYZ to open the file
   !ignore it for most purposes
   !if NAtoms is provided then only the first NAtoms are printed
   
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,ALTmp
      CHARACTER(len=100), OPTIONAL :: filename
      CHARACTER(len=3) :: speciesname
      INTEGER :: iunit1,i,indx,iprint1,NAtoms1,NAtomsRetained
      INTEGER, OPTIONAL :: iprint,iunit,NAtoms,filtermode
      LOGICAL, OPTIONAL :: fileopen,fileclose
      LOGICAL :: fileopen1,fileclose1
      REAL(dp), DIMENSION(3*MaxNFilteredAtoms) :: FilteredAtomCoord
      INTEGER, DIMENSION(MaxNFilteredAtoms) :: FilteredAtomSpecies
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      LOGICAL :: PrintVelocity
      REAL(dp) :: LatticeConst(3)
      CHARACTER(len=10) :: CrystalType
      
      EXTERNAL FilterAtoms
      
      !IF (PRESENT(filtermode)) CALL FilterAtoms(AL,filtermode,FilteredAtomCoord, &
      !   FilteredAtomSpecies,CrystalType,LatticeConst,NAtomsRetained)
      
      INCLUDE "io_writexyz.f90"
      
      IF (PRESENT(NAtoms)) THEN
         NAtoms1=NAtoms
         WRITE(UNIT=iunit1,FMT='(4i10,es20.8," NAtoms Nmove PotentialEnergy")') &
            AL%NAtoms,0,0,0,AL%PotentialEnergy
      ELSE
         NAtoms1=AL%NAtoms
         WRITE(UNIT=iunit1,FMT='(4i10,es20.8," NAtoms Nmove PotentialEnergy")') &
            AL%NAtoms,AL%NMove,AL%PotentialEnergy
      END IF
      WRITE(UNIT=iunit1,FMT='(6f15.8)',ADVANCE="NO") AL%BoxSize
      WRITE(UNIT=iunit1,FMT='("  MD time =",ES12.3," ps")') AL%Time*1.e12_dp
      
      IF (ASSOCIATED(AL%AtomSpecies)) THEN
         AtomSpecies=>AL%AtomSpecies
      ELSE
         ALTmp=>AL
         DO WHILE (ASSOCIATED(ALTmp%PrevNeigh))
            ALTmp=>ALTmp%PrevNeigh
         END DO
         AtomSpecies=>ALTmp%AtomSpecies
      END IF
      
      IF (ASSOCIATED(AL%AtomIsMoving)) THEN
         AtomIsMoving=>AL%AtomIsMoving
      ELSE
         ALTmp=>AL
         DO WHILE (ASSOCIATED(ALTmp%PrevNeigh))
            ALTmp=>ALTmp%PrevNeigh
         END DO
         AtomIsMoving=>ALTmp%AtomIsMoving
      END IF
      
      PrintVelocity=ASSOCIATED(AL%AtomVelocity)
      
      DO i=1,NAtoms1
         indx=AtomSpecies(i)
         IF (indx==0) THEN
            WRITE(6,*) "$Err>> Atom species was found to be zero in WriteXYZ"
            STOP
         END IF
         !speciesname=SpeciesList%AtomicSymbol(AL%SpeciesDirectory(indx))
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(indx))
         IF (PrintVelocity) THEN
            WRITE(UNIT=iunit1,FMT='(a3,3f20.8,"  ",3L2,4f20.8)') &
               !speciesname,AL%AtomCoord(3*i-2:3*i),AtomIsMoving(3*i-2:3*i), &
               !AL%AtomCharge(i),AL%AtomVelocity(3*i-2:3*i)
               speciesname,AL%AtomCoord(3*i-2:3*i)
         ELSE
            WRITE(UNIT=iunit1,FMT='(a3,3f20.8,"  ",3L2,4f20.8)') &
               !speciesname,AL%AtomCoord(3*i-2:3*i),AtomIsMoving(3*i-2:3*i), &
               !AL%AtomCharge(i),0._dp,0._dp,0._dp
               speciesname,AL%AtomCoord(3*i-2:3*i)
         END IF
      END DO
      
      CALL FLUSH(iunit1)
      IF (fileclose1) CLOSE(iunit1)
      
      IF (iprint1>0) THEN
         IF (PRESENT(filename)) THEN
            WRITE(6,*) "IO>> Written atom list to ",TRIM(filename)
         ELSEIF (PRESENT(iunit)) THEN
            WRITE(UNIT=6,FMT='("IO>> Written atom list to file at time",ES12.3)') AL%Time
         END IF
      END IF
      CALL FLUSH(6)
      
   END SUBROUTINE WriteXYZ
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WriteXYZCSP(AL,filename,fileopen,fileclose,iunit,iprint,min_csp,nmax,rcut)
   !fileopen and fileclose has been added for WriteChOS
   !fileopen=.TRUE. asks WriteXYZ to open the file
   !ignore it for most purposes
   !if NAtoms is provided then only the first NAtoms are printed
   !min_csp, nmax and rcut are defined in function CentroSymmetricParameter
   
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      CHARACTER(len=100), OPTIONAL :: filename
      CHARACTER(len=3) :: speciesname
      INTEGER :: iunit1,i,iatom,indx,iprint1,NAtoms1
      INTEGER, OPTIONAL :: iprint,iunit
      LOGICAL, OPTIONAL :: fileopen,fileclose
      LOGICAL :: fileopen1,fileclose1
      REAL(dp) :: LatticeConst(3),csparray(100000)
      CHARACTER(len=10) :: CrystalType
      INTEGER :: ListPrintAtoms(100000),nmax
      REAL :: min_csp,rcut,csp
      
      INCLUDE "io_writexyz.f90"
      
      WRITE(*,*) "The function CentroSymmetricParameter has been modified in amd_utilities.f90"
      STOP
      
      NAtoms1=0
      DO iatom=1,AL%NAtoms
         !csp=CentroSymmetricParameter(AL,iatom,nmax,rcut)
         IF (csp>min_csp) THEN
            NAtoms1=NAtoms1+1
            ListPrintAtoms(NAtoms1)=iatom
            csparray(NAtoms1)=csp
         END IF
      END DO
      
      WRITE(UNIT=iunit1,FMT='(5i8)') NAtoms1,AL%NMove,AL%NAtoms
      WRITE(UNIT=iunit1,FMT='(6f15.8)',ADVANCE="NO") AL%BoxSize
      WRITE(UNIT=iunit1,FMT='("  MD time =",ES12.3," ps")') AL%Time*1.e12_dp
      
      DO i=1,NAtoms1
         iatom=ListPrintAtoms(i)
         csp=csparray(i)
         indx=AL%AtomSpecies(i)
         IF (indx==0) THEN
            WRITE(6,*) "$Err>> Atom species was found to be zero in WriteXYZ"
            STOP
         END IF
         !speciesname=SpeciesList%AtomicSymbol(AL%SpeciesDirectory(indx))
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(indx))
         WRITE(UNIT=iunit1,FMT='(a3,3f20.8,"  ",3L2,4f20.8)') &
           speciesname,AL%AtomCoord(3*iatom-2:3*iatom),AL%AtomIsMoving(3*iatom-2:3*iatom), &
           AL%AtomVelocity(3*iatom-2:3*iatom),csp
      END DO
      
      CALL FLUSH(iunit1)
      IF (fileclose1) CLOSE(iunit1)
      
      IF (iprint1>0) THEN
         IF (PRESENT(filename)) THEN
            WRITE(6,*) "IO>> Written CSP atom list to ",TRIM(filename)
         ELSEIF (PRESENT(iunit)) THEN
            WRITE(UNIT=6,FMT='("IO>> Written CSP atom list to file at time",ES12.3)') AL%Time
         END IF
      END IF
      CALL FLUSH(6)
      
   END SUBROUTINE WriteXYZCSP
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WriteLAMMPSInputCoord(AL,filename,mask)
   !write a LAMMPS input coordinate file with NAtoms present
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      CHARACTER(len=100) :: filename
      INTEGER :: iatom,id,ispecies,q
      LOGICAL, DIMENSION(:), OPTIONAL :: mask
      REAL(dp) :: x,y,z
      LOGICAL :: PrintLine
      
      OPEN(UNIT=301,FILE=filename)
      id=0
      DO iatom=1,AL%NAtoms
         PrintLine=.TRUE.
         IF (PRESENT(mask)) THEN
            IF (mask(iatom)) THEN
               id=id+1
               PrintLine=.TRUE.
            ELSE
               PrintLine=.FALSE.
            END IF
         ELSE
            id=id+1 !always increment id
         END IF
         IF (PrintLine) WRITE(301,*) id,AL%AtomSpecies(iatom),INT(AL%AtomCharge(iatom)), &
            AL%AtomCoord(3*iatom-2:3*iatom)
      END DO
      CLOSE(301)
      WRITE(6,*) ">>> Successfully written the LAMMPS atoms coordinates >> ",TRIM(filename)
   END SUBROUTINE WriteLAMMPSInputCoord
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WriteGlobalAMDStateInfo(state)
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state
      CHARACTER(len=100) :: filename
      CHARACTER(len=NINT2CHAR) :: charstateindx
      INTEGER :: superbasinindex
      
      !Update the state information
      charstateindx=INT2CHAR(state%Index)
      filename="./GlobalAMD/state."//TRIM(charstateindx)//".info"
      OPEN(UNIT=385,FILE=filename)
      superbasinindex=0
      IF (ASSOCIATED(state%superbasin)) superbasinindex=state%superbasin%Index
      WRITE(UNIT=385,FMT='(3I5,ES20.8"  state index, nvisits, superbasin index & potential energy")') &
         state%index,state%NVisits,superbasinindex,state%PotentialEnergy
      WRITE(UNIT=385,FMT='(I5,"  NProcessType")') state%NProcessType
      WRITE(UNIT=385,FMT='(I5,"  NMDTransitions")') state%NMDTransitions
      WRITE(UNIT=385,FMT='(I5,"  NMDFailed")') state%NMDFailed
      WRITE(UNIT=385,FMT='(ES21.8,"  TotalProcessRate")') state%TotalProcessRate
      WRITE(UNIT=385,FMT='(ES15.5,"  delta")') state%delta
      IF (TPMD) THEN
         WRITE(UNIT=385,FMT='(ES21.8,"  basin time elapsed")')state%BasinEscapeTime
      ELSE
         WRITE(UNIT=385,FMT='(ES21.8,"  KMCValidTime")') state%KMCValidTime
      END IF
      WRITE(UNIT=385,FMT='(ES21.8,"  KMCTimeExpired")') state%KMCTimeExpired
      WRITE(UNIT=385,FMT='(ES21.8,"  RateEstimateUnknown")') state%RateEstimateUnknown
      WRITE(UNIT=385,FMT='(ES21.8,"  MDTime")') state%MDTime
      CLOSE(385)
   END SUBROUTINE WriteGlobalAMDStateInfo
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WriteSuperbasinInfo(sb)
      IMPLICIT NONE
      TYPE(SuperbasinContainer), POINTER :: sb
      TYPE(SuperbasinstateContainer), POINTER :: sbstate
      CHARACTER(len=100) :: filename
      CHARACTER(len=NINT2CHAR) :: charisuperbasin
      
      charisuperbasin=INT2CHAR(sb%index)
      filename="./GlobalAMD/superbasin."//TRIM(charisuperbasin)//".info"
      OPEN(UNIT=391,FILE=TRIM(filename))
      WRITE(391,*) sb%Index,sb%NStates," index and # sb states"
      sbstate=>sb%state
      DO WHILE (ASSOCIATED(sbstate))
         WRITE(UNIT=391,FMT='(I6)',ADVANCE="NO") sbstate%AMDstate%Index
         sbstate=>sbstate%NextNeigh
      END DO
      WRITE(391,*) ""
      CLOSE(391)
   END SUBROUTINE WriteSuperbasinInfo
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WriteUbiGraphForAMDStateNetwork(filename,removeorphans)
   !creates a ...
      IMPLICIT NONE
      TYPE(AMDStateContainer), POINTER :: state
      TYPE(AMDProcessContainer), POINTER :: process
     ! TYPE(AMDStateContainer), POINTER :: NextNeigh
     ! TYPE(AMDStateContainer), POINTER :: InitialState, FinalState
      CHARACTER(len=100) :: filename
      CHARACTER(len=NINT2CHAR) :: charstateindx, charstateindx1, charstateindx2
      CHARACTER(len=100) :: statename, statename1
      LOGICAL, OPTIONAL :: removeorphans
      LOGICAL :: removeorphans1
      
      !Step 1. Create the ubigraph file and add initial commands for display
      OPEN(UNIT=381,FILE=TRIM(filename))
      WRITE(UNIT=381,FMT='("import xmlrpclib")')
      WRITE(381,*) "server = xmlrpclib.Server('http://127.0.0.1:20738/RPC2')"
      WRITE(UNIT=381,FMT='("G = server.ubigraph")')

      
      !Step 2. Print vertex and edge information
      
      state=>ListOfGlobalAMDStates
      DO WHILE (ASSOCIATED(state))
         charstateindx=INT2CHAR(state%Index)
         statename="state"//TRIM(charstateindx)
        ! WRITE(6,*) "Adding state ",TRIM(statename)
            
         WRITE(381,*) TRIM(statename)//"=G.new_vertex()"
         WRITE(381,*) "G.set_vertex_attribute("//TRIM(statename)//",'shape','sphere')"
         WRITE(381,*) "G.set_vertex_attribute("//TRIM(statename)//",'color','#ff0000')"
        ! WRITE(381,*) "G.set_vertex_attribute("//TRIM(statename)//",'label','"// &
         !   TRIM(charstateindx)//"')"
         state=>state%NextNeigh
               
      END DO
      
      state=>ListOfGlobalAMDStates
      
      DO WHILE (ASSOCIATED(state))
         charstateindx=INT2CHAR(state%Index)
         process=>state%process
         DO WHILE (ASSOCIATED(process))
            charstateindx2=INT2CHAR(process%InitialState%Index)
            statename="state"//TRIM(charstateindx2)
            charstateindx1=INT2CHAR(process%FinalState%Index)
            statename1="state"//TRIM(charstateindx1)
            WRITE(381,*) "G.set_vertex_attribute("//TRIM(statename)//",'color','#ffff00')"
            !WRITE(381,*) "G.set_vertex_attribute("//TRIM(statename)//",'label','"// &
             !             TRIM(charstateindx2)//"')"
            WRITE(381,*) "G.new_edge("//TRIM(statename),",",TRIM(statename1),")"
            process=>process%NextNeigh
         END DO
         state=>state%NextNeigh
      END DO
         
   END SUBROUTINE WriteUbiGraphForAMDStateNetwork
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ExtractXYZ(filename1,Basename)
   !Purpose: extract XYZ frames from a single file to generate multiple XYZ files 
   !numbered 1 to nframes (nframes is determined by the number of frames in filename1)
      IMPLICIT NONE 

      !CHARACTER (len=100), INTENT(IN):: filename,basename
      CHARACTER (LEN=100):: filename1,Basename,filename2
      !Character(len=100)::Basename
      CHARACTER(len=10)::intchar
      INTEGER :: NAtoms, errorstatus,i,x
      LOGICAL :: NotSatisfied
      REAL, DIMENSION(3) :: AtomCoord	
      CHARACTER(len=3) :: AtomSymbol,tmpchar

      x=0
      NAtoms=0
      errorstatus=0
   
      OPEN(UNIT=351,FILE=TRIM(filename1))
      NotSatisfied = .TRUE.
      DO WHILE(NotSatisfied)
         READ(UNIT=351, FMT=*,IOSTAT=errorstatus) NAtoms
         IF(errorstatus/=0)THEN
            NotSatisfied=.FALSE.
            EXIT
         END IF
      
         x=x+1
         WRITE(intchar,'(I3)')x
         filename2=TRIM(Basename)//TRIM(ADJUSTL(intchar))//".xyz"
         OPEN(unit=352,file=filename2)
         WRITE(352,*)NAtoms
         READ(351,*) tmpchar
         WRITE(352,*) "------"
         DO i=1,NAtoms
            READ(351,*) AtomSymbol,AtomCoord(1:3)
            WRITE(352,*) AtomSymbol, AtomCoord(1:3)
         END DO
         WRITE(6,*) "Written to file:",TRIM(filename2)
         CLOSE(352)
      END DO
      CLOSE(351)
      IF (NotSatisfied) THEN
         WRITE(6,*) "Err>> Error occurred while extracting files in ExtractXYZ"
         STOP
      END IF
   END SUBROUTINE ExtractXYZ
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AnalyzeAL(AL,filename)
      !prints report of the atoms in AL
   
      TYPE(SystemContainer), POINTER :: AL
      CHARACTER(len=20), OPTIONAL :: filename
      INTEGER :: iunit
      
      WRITE(UNIT=*,FMT='(">>Analyzing list of atoms ...")',ADVANCE="NO")
      IF (PRESENT(filename)) THEN
         iunit=301
         WRITE(*,*) " Printing to ",filename(1:len_trim(filename))
         OPEN(UNIT=iunit,FILE=filename(1:len_trim(filename)))
      ELSE
         WRITE(*,*) " Printing to screen"
         iunit=6 !standard output
      END IF
      
      WRITE(iunit,*) "   Details of atomic list"
      WRITE(iunit,*) "   ----------------------"
      WRITE(UNIT=iunit,FMT='("     # atoms  :",i10)') AL%NAtoms
      !WRITE(UNIT=iunit,FMT='("     # species:",i6)') AL%NSpeciesType
      CALL WriteSpecies(AL,iunit)
      WRITE(UNIT=iunit,FMT='("     # moving atoms:",i10)') GetMovingAtoms(AL)
      WRITE(UNIT=iunit,FMT='("     system temperature:",f15.3," K")') GetTemperature(AL)
      WRITE(UNIT=iunit,FMT='("     kinetic energy:",ES12.3," eV")') AL%KineticEnergy
      
      IF (PRESENT(filename)) THEN
         CLOSE(iunit)
      END IF
      
   END SUBROUTINE AnalyzeAL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AnalyzeChoS(chos,iunit)
      TYPE(ChOSContainer), POINTER :: chos
      INTEGER, OPTIONAL :: iunit
      
   END SUBROUTINE AnalyzeChoS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WriteSpecies(AL,iunit)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER :: iunit,i,j
      
      WRITE(UNIT=iunit,FMT='("     species  :")',ADVANCE="NO")
      !DO i=1,AL%NSpeciesType
      DO i=1,NSpecies_global
         !j=AL%SpeciesDirectory(i) !get atomic #
         j=SpeciesDirectory_global(i) !get atomic #
         WRITE(UNIT=iunit,FMT='(a7)',ADVANCE="NO") SpeciesList%AtomicSymbol(j)
      END DO
      WRITE(iunit,*) " "

   END SUBROUTINE WriteSpecies
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WriteVerletList(AL,filename,iunit,iatom)
      !Prints details of Verlet list of atoms in AL

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(VerletListContainer), POINTER :: VL
      INTEGER, DIMENSION(:), OPTIONAL :: iatom
      CHARACTER(len=100) :: filename
      INTEGER, OPTIONAL :: iunit
      INTEGER :: iunit1,iatom1,ListType,MaxAtomPerAtom,NAtoms,i
      
      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(6,*) "$Err>> Atom List is not associated while writing Verlet list"
         RETURN
      ELSE
         IF (.NOT. ASSOCIATED(AL%VL)) THEN
            WRITE(6,*) "$Err>> Verlet List is not associated while writing the list"
            RETURN
         ELSEIF (AL%VL%UnInitialized) THEN
            WRITE(6,*) "$Err>> Verlet List is not initialized while writing the list"
            RETURN
         END IF
      END IF
      
      VL=>AL%VL
      MaxAtomPerAtom=VL%MaxAtomPerAtom
      NAtoms=AL%NAtoms
      IF (VL%ListType/=1 .AND. VL%ListType/=2) THEN
         WRITE(6,*) "Err>> List of unknown type"
         RETURN
      END IF
      
      IF (PRESENT(iunit)) THEN
         iunit1=iunit
      ELSE
         iunit1=302
      END IF
      
      IF (iunit1/=6) OPEN(UNIT=iunit1,FILE=TRIM(filename))
      
      WRITE(UNIT=6,FMT='("Verlet list type: ")',ADVANCE="NO")
      ListType=VL%ListType
      IF (ListType==1) THEN
         WRITE(6,*) "Full List"
      ELSE
         WRITE(6,*) "Half List"
      END IF
      
      IF (PRESENT(iatom)) THEN
         DO i=1,SIZE(iatom)
            iatom1=iatom(i)
            IF (iatom1<1 .OR. iatom1>NAtoms) THEN
               WRITE(6,*) "$Err>> Incorrect list of atoms provided to WriteVerletList"
               STOP
            END IF
            CALL WriteVerletListDriver(iatom1)
         END DO
      ELSE
         DO iatom1=1,NAtoms
            CALL WriteVerletListDriver(iatom1)
         END DO
      END IF
      IF (iunit1/=6) CLOSE(iunit1)
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE WriteVerletListDriver(iatom1)
         IMPLICIT NONE
         INTEGER :: iatom1,r0,r1,r2,ra,j
         
         ra=(iatom1-1)*MaxAtomPerAtom !start index minus 1
         r0=VL%ListRange(iatom1) !range for ith atom
         r1=ra+1
         r2=ra+r0
         
         WRITE(iunit1,*) " "
         WRITE(UNIT=iunit1,FMT='("Position: ",3ES10.3)') AL%AtomCoord(3*iatom1-2:3*iatom1)
         WRITE(UNIT=iunit1,FMT='("Verlet list for atom index: ",i4)') iatom1
         WRITE(UNIT=iunit1,FMT='("Number of neighbors: ",i3)') r0
         WRITE(iunit1,*) " Atom   Distance        Distancex     Distancey      Distancez"
         WRITE(iunit1,*) " ====   ========        =========     =========      ========="
         DO j=r1,r2
            WRITE(UNIT=iunit1,FMT='(i5,4ES15.7)') VL%List(j),VL%drmag(j),VL%dr(3*j-2:3*j)
         END DO
      END SUBROUTINE WriteVerletListDriver
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE WriteVerletList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE WritePotential(AL,iunit,filename,IsDetailed) !,iprint)

      !Prints all details of the potential -- IsDetailed is inactive as of now

      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(InteracPotential), POINTER :: p
      INTEGER, OPTIONAL :: iunit !,iprint
      CHARACTER(len=100), OPTIONAL :: filename
      INTEGER :: iunit1,i,j,species1,species2
      TYPE(EAMPotential), POINTER :: EAM
      LOGICAL, OPTIONAL :: IsDetailed
      
      IF (PRESENT(iunit)) THEN
         iunit1=iunit
         IF (.NOT. PRESENT(filename)) THEN
            WRITE(6,*) "$Err>> Filename expected in WritePotential"
            STOP
         END IF
         OPEN(UNIT=iunit1,FILE=TRIM(filename))
      ELSE
         iunit1=6
      END IF
      
      
      WRITE(iunit1,*) ""
      WRITE(iunit1,*) " Overview of potential for Atom List"
      WRITE(iunit1,*) " +++++++++++++++++++++++++++++++++++"

      p=>AL%Potential
      !NSpecies= p%NSpecies
      WRITE(UNIT=iunit1,FMT='("   NSpeciesType:",i3)') NSpecies_global
      
      !IF (PRESENT(IsDetailed)) THEN
      !   IF (IsDetailed) WRITE(6,*) "EAM potential details are being printed to PotentialType.details.*"
      !END IF
      
      IF (ASSOCIATED(p%EAM)) THEN
         WRITE(6,*) "   EAM potential ................"
         DO species1=1,NSpecies_global
            DO species2=species1,NSpecies_global
            END DO
         END DO
      END IF
      IF (ASSOCIATED(p%LJ)) THEN
         WRITE(6,*) "   LJ potential ................"
         DO species1=1,NSpecies_global
            DO species2=species1,NSpecies_global
            END DO
         END DO
      END IF
      IF (ASSOCIATED(p%Coulomb)) THEN
         WRITE(6,*) "   Coulomb potential ................"
      END IF
      IF (ASSOCIATED(p%Buckingham)) THEN
         WRITE(6,*) "   Buckingham potential ................"
         DO species1=1,NSpecies_global
            DO species2=species1,NSpecies_global
            END DO
         END DO
      END IF
      IF (ASSOCIATED(p%SW)) THEN
         WRITE(6,*) "   SW potential ................"
      END IF
      IF (ASSOCIATED(p%Tersoff)) THEN
         WRITE(6,*) "   Tersoff potential ................"
      END IF

      WRITE(iunit1,*) " +++++++++++++++++++++++++++++++++++"
   END SUBROUTINE WritePotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintInfo_AtomicNo(n)

      !Prints the information for atomic # n

      IMPLICIT NONE
      INTEGER :: n
      
      WRITE(*,*) "++++++++++++PRINTING SPECIES INFORMATION+++++++++++++++++++"
      WRITE(*,*) "   Name:",TRIM(SpeciesList%AtomicName(n))
      WRITE(*,*) "   Symbol:",TRIM(SpeciesList%AtomicSymbol(n))
      WRITE(UNIT=*,FMT='("    AtomicNumber:",I4)') SpeciesList%AtomicNumber(n)
      WRITE(UNIT=*,FMT='("    AtomicMass:",F8.3)') SpeciesList%AtomicMass(n)
      WRITE(UNIT=*,FMT='("    Melting point:",F8.3)') SpeciesList%AtomicMeltPt(n)
      WRITE(UNIT=*,FMT='("    Boiling point:",F8.3)') SpeciesList%AtomicBoilPt(n)
      WRITE(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

   END SUBROUTINE PrintInfo_AtomicNo
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintInfo_AtomicSymbol(name)

      !Prints the information for atomic species

      IMPLICIT NONE
      CHARACTER(len=3) :: name
      INTEGER :: n,l
      LOGICAL :: Found
      
      Found=.FALSE.
      DO n=1,SpeciesList%NSpecies
         l=LEN_TRIM(name)
         IF (name(1:l)==SpeciesList%AtomicSymbol(n)(1:l)) THEN
            Found=.TRUE.
            EXIT
         END IF
      END DO

      IF (.NOT. Found) THEN
         WRITE(*,*) "$Err>> PrintInfo unable to find the atom symbol"
         STOP
      END IF

      WRITE(*,*) "++++++++++++PRINTING SPECIES INFORMATION+++++++++++++++++++"
      WRITE(*,*) "   Name:",TRIM(SpeciesList%AtomicName(n))
      WRITE(*,*) "   Symbol:",TRIM(SpeciesList%AtomicSymbol(n))
      WRITE(UNIT=*,FMT='("    AtomicNumber:",I4)') SpeciesList%AtomicNumber(n)
      WRITE(UNIT=*,FMT='("    AtomicMass:",F8.3)') SpeciesList%AtomicMass(n)
      WRITE(UNIT=*,FMT='("    Melting point:",F8.3)') SpeciesList%AtomicMeltPt(n)
      WRITE(UNIT=*,FMT='("    Boiling point:",F8.3)') SpeciesList%AtomicBoilPt(n)
      WRITE(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

   
   END SUBROUTINE PrintInfo_AtomicSymbol
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

END MODULE IO
