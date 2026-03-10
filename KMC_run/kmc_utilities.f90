MODULE KMCUtilities
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur

   USE KMC_VARIABLE_TYPE
   USE utilities, ONLY : INT2CHAR,NINT2CHAR
   REAL(dp) :: BODimensionlessRMax,BODecayRate,BOSpacing
   
   CONTAINS
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    FUNCTION CheckIfRelaxEachTime(v)
      IMPLICIT NONE
      INTEGER :: v
      LOGICAL :: CheckIfRelaxEachTime
      
      CheckIfRelaxEachTime= (v==1)

    END FUNCTION CheckIfRelaxEachTime
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    FUNCTION CheckIfSnapOn(v)
      IMPLICIT NONE
      INTEGER :: v
      LOGICAL :: CheckIfSnapOn
      
      CheckIfSnapOn= (v==2)

   END FUNCTION CheckIfSnapOn
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE InitializeUtilities()
      IMPLICIT NONE
      !initialization steps before sbrtns in this module can be used
      
      BODimensionlessRMax=1.08
      BODecayRate=2.
      BOSpacing=KMCMaxLatticeSpacing*1.08 !Tolerance
      WRITE(UnitScrap,*) "Ini>> Max lattice spacing:",KMCMaxLatticeSpacing
      WRITE(UnitScrap,*) "Ini>> BO spacing used:",BOSpacing
   END SUBROUTINE InitializeUtilities
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetBondOrder(rad,Atom)
   !Obtains the bond for an atom with another atom based on the spacing rad
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      REAL :: rad
      REAL :: expterm,deltar
      REAL :: GetBondOrder
      TYPE(KMCAtom), POINTER :: Atom
      
      IF (.NOT. ASSOCIATED(Atom)) THEN
         WRITE(UnitScrap,*) "KMCUtilities>>...site passed to fn BondOrder is NULL"
         STOP
      END IF
      
      SELECT CASE (BondOrderFlag)
      CASE(1) !uses step function
         GetBondOrder=0.
         IF (rad<BOSpacing) GetBondOrder=1.
      CASE(2) !use logistic function
         deltar=rad/KMCMaxLatticeSpacing-BODimensionlessRMax
         expterm=exp(-BODecayRate*deltar)
         GetBondOrder=expterm/(1._dp+expterm)
      END SELECT
   END FUNCTION GetBondOrder
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetCellLocation(Coord)
      IMPLICIT NONE
      INTEGER :: GetCellLocation(3)
      REAL, INTENT(IN) :: Coord(3)
      
      !GetCellLocation=MIN(FLOOR(Coord/CellSize)+1,NCells)
      GetCellLocation=CEILING(Coord/KMCCellSize)
      WHERE (GetCellLocation==0) GetCellLocation=1
   END FUNCTION GetCellLocation
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetPTADCoeff()
      IMPLICIT NONE
      INTEGER :: i
      REAL :: kt,Delta
   
      TADAlpha=1.
      Delta=0.25
      DO i=2,100
         TADAlpha=TADAlpha*1.1
         Delta=10.**(TADAlpha*LOG10(TADAlpha)-(TADAlpha+1.)*LOG10(TADAlpha+1.))
         IF (Delta<TADDelta) EXIT
      END DO
 
      kt=LOG((TADAlpha+1)/TADAlpha)
      TADpmi=EXP(-TADAlpha*kt)
      TADpma=1.-EXP(-kt)

      TADNuMin=TADMinPrefac/LOG(1./TADpmi)
      WRITE(UNIT=UnitScrap,FMT='("Ini>>TADAlpha:",f10.3)') TADAlpha
      WRITE(UNIT=UnitScrap,FMT='("Ini>>TADpmi:",f10.3)')  TADpmi
      WRITE(UNIT=UnitScrap,FMT='("Ini>>TADpma:",f10.3)') TADpma
      WRITE(UNIT=UnitScrap,FMT='("Ini>>kt:",f10.3)'),kt
      WRITE(UNIT=UnitScrap,FMT='("Ini>>TADNuMin:" ,es10.3)'),TADNuMin
    END SUBROUTINE GetPTADCoeff
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetOrigin
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      REAL :: MinCoord(3)=0.
      INTEGER :: i
      SAVE MinCoord
  
      WRITE(UNIT=6,FMT='("Ini>> Translating atoms to set origin ...")',ADVANCE="NO")
      atom=>KMC_AL
      MinCoord=atom%coord
      DO i=2,NAtomKMC
         atom=>atom%NextNeigh
         WHERE(MinCoord>atom%Coord) MinCoord=atom%Coord
      END DO
   
      atom=>KMC_AL
      DO i=1,NAtomKMC
         atom%Coord=atom%Coord-MinCoord
         atom=>atom%NextNeigh
      END DO
      WRITE(6,*) " [DONE]"
   END SUBROUTINE SetOrigin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetNeighCell(c,ix,iy,iz)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: c(3),ix,iy,iz
      INTEGER :: dc(3),GetNeighCell(3)
      
      dc=c+(/ix,iy,iz/)
      WHERE (dc<1) dc=dc+NCells
      WHERE (dc>NCells) dc=dc-NCells
      GetNeighCell=dc
   END FUNCTION GetNeighCell
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetkBT(Temperature)
      IMPLICIT NONE
      REAL :: GetkBT
      REAL(dp) :: kB,eV
      REAL, INTENT(IN) :: Temperature

      kB=1.3806503e-23
      eV=1.60217646e-19
      GetkBT=kB*Temperature/eV
   END FUNCTION GetkBT
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetPBCSpacing(coord1,coord2)
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      REAL, DIMENSION(3) :: GetPBCSpacing,vec
      REAL, DIMENSION(3), INTENT(IN) :: coord1,coord2

      vec=coord1-coord2
      WHERE (ABS(vec)>KMCBoxSize-ABS(vec)) vec=SIGN(KMCBoxSize-ABS(vec),-vec) !Min image convention
      GetPBCSpacing=vec
    END FUNCTION GetPBCSpacing
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetKMCSumProcessRates
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessList), POINTER :: ActiveProcess
      
      KMCSumProcessRates=0._dp
      ActiveProcess=>ActiveProcessList
      DO WHILE (ASSOCIATED(ActiveProcess))
         prc=>ActiveProcess%Process
         KMCSumProcessRates=KMCSumProcessRates+prc%Rate*REAL(prc%NumberSubscribers,dp)
         ActiveProcess=>ActiveProcess%NextNeigh
      END DO
   END SUBROUTINE GetKMCSumProcessRates
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetKMCSumProcessRates1
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      
      KMCSumProcessRates=0._dp
      prc=>KMC_PL
      DO WHILE (ASSOCIATED(prc))
         KMCSumProcessRates=KMCSumProcessRates+prc%Rate*REAL(prc%NumberSubscribers,dp)
         prc=>prc%NextNeigh
      END DO
   END SUBROUTINE GetKMCSumProcessRates1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetPBCCoord(coord)
      IMPLICIT NONE
      REAL, DIMENSION(3) :: GetPBCCoord
      REAL, DIMENSION(3), INTENT(IN) :: Coord

      GetPBCCoord=Coord
      WHERE (Coord>KMCBoxSize) GetPBCCoord=Coord-KMCBoxSize
      WHERE (Coord<0.) GetPBCCoord=Coord+KMCBoxSize
      RETURN
    END FUNCTION GetPBCCoord
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetRadiusBin(Radius)
   !Find the bin corresponding to a distance Radius
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      INTEGER :: GetRadiusBin,i
      REAL :: Radius

      i=CEILING(Radius/BinSizeInEnv)
      IF (REAL(i-1)*BinSizeInEnv>=Radius) i=i-1   !Radius must be lying right on bin border
      IF (i<=0 .OR. i>NumberBinsInEnv) THEN
         WRITE(UnitScrap,*) 'Ini>>Radius:',Radius
         STOP
      END IF
      GetRadiusBin=i
   END FUNCTION GetRadiusBin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetCurrentCPUTime(IsNotInitialized)
      IMPLICIT NONE
      REAL, SAVE :: time_init
      REAL :: CPUTime,GetCurrentCPUTime,TARRAY(2)
      LOGICAL, INTENT(IN) :: IsNotInitialized
   
      CALL CPU_TIME(CPUTime)
      !WALL TIME
!      CPUTime=ETIME(TARRAY) !ETIME is present gfortran but not in Intel
      IF (IsNotInitialized) THEN
         time_init=CPUTime
      END IF
      GetCurrentCPUTime=(CPUTime-time_init) ! in seconds
      RETURN
   END FUNCTION GetCurrentCPUTime
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetBinStatus(rad,BinRegion,BinNumber)
   !The bin contains a skin layer to indicate the lack of uncertainty in the position 
   !of an atom. An atom lying in the skin region can belong to both residing bin and in the
   !bin sharing the boundary
   !BinRegion is -1 if neighbor atom lies inside the inner skin of the bin
   !BinRegion is  0 if neighbor atom lies inside the inner skin of the bin
   !BinRegion is  1 if neighbor atom lies inside the inner skin of the bin
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: BinRegion,BinNumber
      REAL, INTENT(IN) :: rad

      BinRegion=-2
      BinNumber=-10
      IF (rad<EnvCutOffRadius) THEN
        BinNumber=GetRadiusBin(rad)
        BinRegion=0
        IF (rad+BinSizeInEnv - BinSizeInEnv*REAL(BinNumber) < RadialTolerance) BinRegion=-1
        IF (REAL(BinNumber)*BinSizeInEnv - rad < RadialTolerance) BinRegion=1
      ELSE
         BinRegion=1
         BinNumber=NumberBinsInEnv !assume this point lies in the outer edge
      END IF
    END SUBROUTINE GetBinStatus
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GetWelcomeMessage
      IMPLICIT NONE
      WRITE(UnitScrap,*) '               ================================='
      WRITE(UnitScrap,*) '                            KMC code            '
      WRITE(UnitScrap,*) '               ================================='
    END SUBROUTINE GetWelcomeMessage
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintSnapshot(snaptype)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: snaptype
      REAL(dp) :: TimeSnap=0._dp
      SAVE

      CALL PrintSystemState()

      SELECT CASE (snaptype)
      CASE (1)
         IF (MOD(KMCStep,iFreqSnapshot)==0 .OR. KMCBlock<0) CALL PrintAtoms(PrintAtomsFormat)
      CASE (2)
         IF (KMCTime>TimeSnap) THEN
            CALL PrintAtoms(PrintAtomsFormat)
            TimeSnap=TimeSnap+TimeFreqSnapshot
         END IF
      CASE (3)
         CALL UserDefnSnapshot  !user defined
      CASE(4)

      END SELECT
   END SUBROUTINE PrintSnapshot
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE UserDefnSnapshot
   END SUBROUTINE UserDefnSnapshot
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintSpeciesNotation()
      IMPLICIT NONE
      INTEGER :: i

      WRITE(UnitScrap,*) " "
      WRITE(UnitScrap,*) "Species  Number"
      DO i=1,NSpeciesType
         WRITE(UnitScrap,*) SpeciesName(i),'  ',i
      END DO
   END SUBROUTINE PrintSpeciesNotation
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintSystemState()
      IMPLICIT NONE
      WRITE(UNIT=UnitTimeSnaps,FMT='(2i14,es20.4,4i14,3f20.4,i14,3es20.4)',ADVANCE="NO") &
           KMCBlock,KMCStep,KMCTime,NumberCfgTypes,NAtomKMC
      CALL FLUSH(UnitTimeSnaps)
    END SUBROUTINE PrintSystemState
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintKMCStatus()
   !prints the current kmc steps, time, number of configurations and processes
      IMPLICIT NONE
      INTEGER(I4B) :: KMCSteps
      
      CPUTime(6)=GetCurrentCPUTime(.FALSE.)
      
      KMCSteps=INT(KMCBlock,I4B)*INT(NKMCStepsPerBlock,I4B)+INT(KMCStep,I4B)
      WRITE(UNIT=UnitTimeSnaps,FMT='(I12,2ES17.8,3I6)') KMCSteps,KMCTime,CPUTime(6),NumberCfgTypes,NumberProcesses,NAtomKMC
      CALL FLUSH(UnitTimeSnaps)
   END SUBROUTINE PrintKMCStatus
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintAtoms(Printformat,fileopen)
      IMPLICIT NONE
      CHARACTER :: Printformat
      TYPE(KMCAtom), POINTER :: atom
      INTEGER :: DecimalPos
      CHARACTER(len=NINT2CHAR) :: CDP
      LOGICAL, OPTIONAL :: fileopen
      LOGICAL :: fileopen1
      
      IF (PRESENT(fileopen)) THEN
         fileopen1=fileopen
      ELSE
         fileopen1=.TRUE.
      END IF

      IF (ForceNoPrint) RETURN
      TotalCalls(7)=TotalCalls(7)+1
      CPUTime(7)=CPUTime(7)-GetCurrentCPUTime(.FALSE.)
      DecimalPos=INT(LOG10(REAL(TotalCalls(7))))

      WRITE(6,*) "Printing lattice"

      IF (fileopen1) THEN
         SnapCounter=SnapCounter+1
         CDP=INT2CHAR(SnapCounter)
      END IF
      !CDP=INT2CHAR(TotalCalls(7))
      FileLatticeSnap='./LatticeSnaps/L'
      FileLatticeSnap='./LatticeSnaps/L.'//TRIM(CDP)//'.dat'
      IF (Printformat=='x') FileLatticeSnap='./LatticeSnaps/L.'//TRIM(CDP)//'.xyz'

      IF (fileopen1) OPEN(Unit=UnitLatticeSnap,FILE=FileLatticeSnap)
      SELECT CASE(Printformat)
      CASE('x')
         !WRITE(UnitLatticeSnap,FMT='(2i14)') NAtomKMC,NMoveKMC
         WRITE(UnitLatticeSnap,FMT='(i14)') NAtomPrint
         WRITE(UnitLatticeSnap,*) 'Snap'
         atom=>KMC_AL
         DO WHILE (ASSOCIATED(atom))
            IF (atom%Index<=NAtomPrint) &
               WRITE(UnitLatticeSnap,*) SpeciesName(atom%Species),atom%Coord
            atom=>atom%NextNeigh
         END DO
      CASE('c') !clsman o/p
         WRITE(UnitLatticeSnap,FMT='(i5)') NAtomKMC !516
         WRITE(UnitLatticeSnap,*) 'Snap>>Current KMC iter=', &
           KMCBlock*NKMCStepsPerBlock+KMCStep,'Current Time=',KMCTime
         WRITE(UnitLatticeSnap,FMT='(3f20.8)') KMCBoxSize
         atom=>KMC_AL
         DO WHILE (ASSOCIATED(atom))
            !IF (atom%Coord(3)>15.8+2.045) THEN
            !   WRITE(UnitLatticeSnap,FMT='(3f20.8,i5)') (atom%Coord(indx),indx=1,3),atom%iSpecies+1
            !ELSE
            WRITE(UnitLatticeSnap,FMT='(3f20.8,i5)') atom%Coord,atom%Species
            !END IF
            atom=>atom%NextNeigh
            !IF (atom%Coord(3)<12.8) EXIT
         END DO
      CASE('o')
         WRITE(UnitLatticeSnap,FMT='(i5)') NAtomKMC !516
         WRITE(UnitLatticeSnap,*) 'Snap>>Current KMC iter=', &
            KMCBlock*NKMCStepsPerBlock+KMCStep,'Current Time=',KMCTime
         WRITE(UnitLatticeSnap,FMT='(3f20.8)') KMCBoxSize
         atom=>KMC_AL
         DO WHILE (ASSOCIATED(atom))
            WRITE(UnitLatticeSnap,FMT='(3f20.8,i5)') atom%Coord,atom%Species
            atom=>atom%NextNeigh
         END DO
      CASE('s')
         WRITE(UnitLatticeSnap,FMT='(3f20.8,i5)') atom%Coord,atom%Species
      END SELECT
      IF (fileopen1) CLOSE (UnitLatticeSnap)
      
      OPEN(UNIT=UnitStoreState,FILE=TRIM(FileStoreState))
      WRITE(UnitStoreState,FMT='(2i14)') NAtomKMC,NMoveKMC
!      WRITE(UnitStoreState,*) 'Last saved snapshot'
      WRITE(UnitStoreState,FMT='(3f20.8)') KMCBoxSize
      atom=>KMC_AL
      DO WHILE (ASSOCIATED(atom))
   WRITE(UnitStoreState,*) TRIM(SpeciesName(atom%Species)),atom%Coord,atom%IsMoving
         atom=>atom%NextNeigh
      END DO
      CLOSE(UNIT=UnitStoreState)
      
      CPUTime(7)=CPUTime(7)+GetCurrentCPUTime(.FALSE.) !update time for printing lattice
   END SUBROUTINE PrintAtoms
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintEnv(atom)
      IMPLICIT NONE
      TYPE (Bin), POINTER :: bin1 !needed if you want to print config
      TYPE (KMCAtom), POINTER :: atom
      INTEGER :: AtomIndex
      INTEGER :: Species
      TYPE (BinAtom), POINTER :: binatom1
 
      AtomIndex=atom%Index
      WRITE(UNIT=UnitScrap,FMT='("Cfg>>------Summary of neighbors for atom # ",I5,"-------")') AtomIndex
      WRITE(UnitScrap,FMT='("Cfg>>Atom Coord:",3(f8.3," "))') atom%Coord
      WRITE(UnitScrap,FMT='("Cfg>>Atom TCell:",3(I5," "))') atom%CellIndex
      WRITE(UnitScrap,*) 'Cfg>>Species',atom%Species
      WRITE(UNIT=UnitScrap,FMT='("Config #:",I4)') atom%Cfg%ConfigurationIndex
      WRITE(UnitScrap,*) 'Cfg>>Local configuration: # neighbors ',atom%Env%ShortRangeConfig%Population
      IF (ASSOCIATED(atom%Env%ShortRangeConfig%BinInfo)) THEN
         bin1=>atom%Env%ShortRangeConfig%BinInfo
         DO WHILE (ASSOCIATED(bin1))
            WRITE(UnitScrap,FMT='("   Bin# [",I5,"]")') bin1%BinNumber
            DO Species=1,NSpeciesType
               WRITE(UnitScrap,FMT='("   Species ",I5," Population(left,center,right) ",3I4)') &
                 Species,bin1%BinEdgeN(Species),bin1%InBin(Species),bin1%BinEdgeP(Species)
            END DO
            binatom1=>bin1%AL
            IF (ASSOCIATED(binatom1)) WRITE(UNIT=UnitScrap,FMT='("   Bin Contents: Atom (distance):")',ADVANCE='NO')
            DO WHILE (ASSOCIATED(binatom1))
               WRITE(UNIT=UnitScrap,FMT='(I5,"(",4f6.3,i4") ")',ADVANCE='YES') &
                 binatom1%Atom%Index,binatom1%distance,binatom1%RelCoord,binatom1%Species
               IF (ASSOCIATED(binatom1%NextNeigh)) THEN
                  IF (binatom1%NextNeigh%PrevNeigh%Atom%Index/=binatom1%Atom%Index) THEN
                     WRITE(6,*) "Cfg>>Mismatch in atom index expected and found"
                     STOP
                  END IF
               END IF
               IF (ASSOCIATED(binatom1%PrevNeigh)) THEN
                  IF (binatom1%PrevNeigh%NextNeigh%Atom%Index/=binatom1%Atom%Index) THEN
                     WRITE(6,*) "Cfg>>Mismatch in atom index expected and found"
                     STOP
                  END IF
               END IF
               binatom1=>binatom1%NextNeigh
            END DO
            WRITE(UnitScrap,*) ' ' !to undo effect of ADVANCE=NO
            bin1=>bin1%NextNeigh
         END DO
      ELSE
         WRITE(UnitScrap,*) 'Cfg>>No neighbor present'
      END IF

      WRITE(UNIT=UnitScrap,FMT='("------End of summary of neighbors for atom # ",I5,"-------")') AtomIndex
   END SUBROUTINE PrintEnv
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintCfgClasses()
      IMPLICIT NONE
      TYPE (Configuration), POINTER :: cfg
      TYPE(KMCAtom), POINTER :: atom1
      INTEGER :: j
      INTEGER, SAVE :: CallNoRecordedConfigsFile=0
      CHARACTER(len=30) :: RCFileName
   
      WRITE(UnitScrap,*) 'Cfg>>---Types of configs present---'
      WRITE(UnitScrap,*) 'Cfg>>No of configs present:',NumberCfgTypes

      cfg=>RecordedCfg
      IF (IsConfigPrint) THEN
         CALL SYSTEM("mkdir -v ./CfgImg/")
         DO WHILE (ASSOCIATED(cfg))
            IF (cfg%ConfigurationIndex/=RedundantCfg%ConfigurationIndex) CALL PrintNewCfgClass(cfg)
            cfg=>cfg%NextNeigh
         END DO
      END IF
   
      !print atoms based on different config #
      atom1=>KMC_AL
      CallNoRecordedConfigsFile=CallNoRecordedConfigsFile+1
      RCFileName='RC.' // & !CHAR(48+(CallNoRecordedConfigsFile/100)) // &
        CHAR(48+(CallNoRecordedConfigsFile/10)) // CHAR(48+MOD(CallNoRecordedConfigsFile,10)) // &
        '.dat'

      OPEN(Unit=130,FILE=RCFileName)
      WRITE(UNIT=130,FMT='(i5)') NAtomKMC
      WRITE(130,*) 'configs'
      WRITE(UNIT=130,FMT='(3f20.8)') KMCBoxSize
      DO WHILE (ASSOCIATED(atom1))
         IF (ASSOCIATED(atom1%Cfg)) THEN
            WRITE(UNIT=130,FMT='(3f20.8,i5)') (atom1%Coord(j), j=1,3), &
              MOD(atom1%Cfg%Configurationindex,10)+1
         ELSE
            WRITE(UNIT=130,FMT='(3f20.8,i5)') (atom1%Coord(j), j=1,3),1
         END IF
         atom1=>atom1%NextNeigh
      END DO
      CLOSE(130)
   END SUBROUTINE PrintCfgClasses
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintNewCfgClass(cfg)
      IMPLICIT NONE
      TYPE (Configuration), POINTER :: cfg
      TYPE (KMCAtomList), POINTER :: AL1
      INTEGER :: j,NoAtoms
      TYPE(Bin), POINTER :: bin1
      TYPE(BinAtom), POINTER :: binatom1
      CHARACTER(len=NINT2CHAR) :: ChCfg
      CHARACTER(len=30) :: FileConfig

  
      ChCfg=INT2CHAR(cfg%ConfigurationIndex)
      FileConfig='./CfgImg/C'//TRIM(ChCfg)
      OPEN(UNIT=130,FILE=TRIM(FileConfig))
      !sitepointer=>cfg%SitesOfAConfigType%SitePointer !print an example of this config
      
      bin1=>cfg%Env%ShortRangeConfig%BinInfo
      NoAtoms=1
      DO WHILE (ASSOCIATED(bin1))
         binatom1=>bin1%AL
         DO WHILE (ASSOCIATED(binatom1))
            NoAtoms=NoAtoms+1
            binatom1=>binatom1%NextNeigh
         END DO
         bin1=>bin1%NextNeigh
      END DO

      bin1=>cfg%Env%ShortRangeConfig%BinInfo
      SELECT CASE (PrintAtomsFormat)
      CASE ('c','s')
         WRITE(UNIT=130,FMT='(i5)') NoAtoms

         WRITE(UNIT=130,FMT='("Cfg>>------Configuration #",i5)',ADVANCE='NO') cfg%ConfigurationIndex
         WRITE(UNIT=130,FMT='(" # atoms:",i5)',ADVANCE='NO') cfg%NAtomsWithCfg
         AL1=>cfg%AtomsWithCfg
         WRITE(UNIT=130,FMT='(" List of atoms#:")',ADVANCE='NO') 
         DO WHILE (ASSOCIATED(AL1))
            WRITE(UNIT=130,FMT='(I5,"(",I5,") ")',ADVANCE='NO') &
            AL1%Atom%Index,AL1%Atom%Cfg%ConfigurationIndex
            AL1=>AL1%NextNeigh
         END DO
         IF (ASSOCIATED(cfg%NextNeigh)) THEN
            IF (cfg%NextNeigh%PrevNeigh%ConfigurationIndex/=cfg%ConfigurationIndex) THEN
               WRITE(UnitScrap,*) 'Cfg>>Mismatch in iConfig expected and found'
               STOP
            END IF
         END IF
         IF (ASSOCIATED(cfg%PrevNeigh)) THEN
            IF (cfg%PrevNeigh%NextNeigh%ConfigurationIndex/=cfg%ConfigurationIndex) THEN
               WRITE(UnitScrap,*) 'Cfg>>Mismatch in iConfig expected and found'
               STOP
            END IF
         END IF
         
         WRITE(130,*) ' '
         WRITE(UNIT=130,FMT='(3f20.8)') (KMCCellSize(j), j=1,3)
         WRITE(UNIT=130,FMT='(3f20.8,i5)') KMCCellSize(1)/2.,KMCCellSize(2)/2.,KMCCellSize(3)/2.,2
         DO WHILE (ASSOCIATED(bin1))
            binatom1=>bin1%AL
            DO WHILE (ASSOCIATED(binatom1))
               WRITE(UNIT=130,FMT='(3f20.8,i5)') &
                 (binatom1%RelCoord(j)+KMCCellSize(j)/2., j=1,3),binatom1%Species
               binatom1=>binatom1%NextNeigh
            END DO
            bin1=>bin1%NextNeigh
         END DO
      CASE ('x')
         WRITE(130,*) ' '
         WRITE(130,*) SpeciesName(2), 0.,0.,0.
         DO WHILE (ASSOCIATED(bin1))
            binatom1=>bin1%AL
            DO WHILE (ASSOCIATED(binatom1))
               WRITE(UnitScrap,*) SpeciesName(1),(binatom1%RelCoord(j), j=1,3)
               binatom1=>binatom1%NextNeigh
            END DO
            bin1=>bin1%NextNeigh
         END DO
      CASE DEFAULT
         WRITE(6,*) 'Cfg>>Incorrect option selected for PrintAtomsFormat'
         STOP
      END SELECT
      CLOSE(130)
   END SUBROUTINE PrintNewCfgClass
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintAffectedCells(iunit)
   !Prints the list of affected cells to a file
      IMPLICIT NONE
      INTEGER, OPTIONAL :: iunit
      INTEGER :: iunit1,i
      
      IF (PRESENT(iunit)) THEN
         iunit1=iunit
      ELSE
         iunit1=6
      END IF
      
      WRITE(UNIT=iunit1,FMT='("Number of affected cells:",I4)') NAffectedCells
      DO i=1,NAffectedCells
         WRITE(UNIT=iunit1,FMT='(3I4)') AffectedCells(3*i-2:3*i)
      END DO
      WRITE(iunit1,*) "End of list"
   END SUBROUTINE PrintAffectedCells
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintAffectedCellAtoms(n,coord)
   !coord is the coordinate which would be origin
      IMPLICIT NONE
      INTEGER :: n,i,indx(3)
      TYPE(KMCAtomList), POINTER :: AL
      REAL :: coord(3),coord1(3)
      
      OPEN(UNIT=391,FILE="./StoreTADKMC/AffectedCellAtoms.xyz")
      WRITE(391,*) "Number of cells to look:",n
      DO i=1,n
         indx=AffectedCells(3*i-2:3*i)
         AL=>Cell(indx(1),indx(2),indx(3))%AL
         DO WHILE (ASSOCIATED(AL))
            coord1=GetPBCSpacing(AL%atom%Coord,coord)
            WRITE(391,*) "Ag  ",coord1
            AL=>AL%NextNeigh
         END DO
      END DO
      CLOSE(391)
   END SUBROUTINE PrintAffectedCellAtoms
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintCfg(cfg)
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: cfg
      TYPE(Histogram), POINTER :: hist1
      TYPE(BinAtom), POINTER :: binatom1
      TYPE(Bin), POINTER :: bin1

      WRITE(UnitScrap,*) "Details of histogram",ASSOCIATED(cfg)
      hist1=>cfg%Env%ShortRangeConfig
      WRITE(UNIT=UnitScrap,FMT='("BO:",f14.4)') cfg%BondOrder
      WRITE(UNIT=UnitScrap,FMT='("Population:",i5)') hist1%Population

      bin1=>hist1%BinInfo
      DO WHILE(ASSOCIATED(bin1))
         WRITE(UNIT=UnitScrap,FMT='("Bin-",i4)') bin1%BinNumber
         binatom1=>bin1%AL
         DO WHILE (ASSOCIATED(binatom1))
            WRITE(UNIT=UnitScrap,FMT='("  ",4f14.3,i4)') &
             binatom1%distance,binatom1%RelCoord,binatom1%Species
            binatom1=>binatom1%NextNeigh
         END DO
         bin1=>bin1%NextNeigh
      END DO
   END SUBROUTINE PrintCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintProcess(process,iunit)
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: process
      TYPE(KMCProcessAtom), POINTER :: atom
      TYPE(KMCAtomList), POINTER :: atom1
      INTEGER, OPTIONAL :: iunit
      INTEGER :: iunit1

      IF (PRESENT(iunit)) THEN
         iunit1=iunit
      ELSE
         iunit1=UnitScrap
      END IF
      
      WRITE(iunit1,*) "Details of process:"
      WRITE(UNIT=iunit1,FMT='("Index=",I8)') process%Index
      WRITE(UNIT=iunit1,FMT='("Configuration Index=",I8)') process%CfgType%ConfigurationIndex
      WRITE(UNIT=iunit1,FMT='("Rate=",1ES14.4," # atoms=",I4)') process%Rate,process%NumberProcessAtoms
      WRITE(UNIT=iunit1,FMT='("Freq=",ES14.4," barriers=",F14.3 )') process%Frequency, &
          process%FwdBarrier
      WRITE(iunit1,*) "Atoms involved:"
      atom=>process%ReactantAL
      DO WHILE (ASSOCIATED(atom))
         WRITE(UNIT=iunit1,FMT='("Initial rcoord:",3f14.4)') atom%RelInitialPos
         WRITE(UNIT=iunit1,FMT='("Final rcoord:",3f14.4)') atom%RelFinalPos
         atom=>atom%NextNeigh
      END DO
      WRITE(UNIT=iunit1,FMT='("Number subscribers=",I8)') process%NumberSubscribers
      WRITE(UNIT=iunit1,FMT='("List of subscribers:")',ADVANCE="NO") 
      atom1=>process%AL
      DO WHILE (ASSOCIATED(atom1))
         WRITE(UNIT=iunit1,FMT='(I5," ")',ADVANCE="NO") atom1%Atom%Index
         atom1=>atom1%NextNeigh
      END DO
      WRITE(iunit1,*) " "
   END SUBROUTINE PrintProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintProcessAssignments()
   !Prints the number of subscribers to each process
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCAtomList), POINTER :: AL
      INTEGER :: nentries
      
      OPEN(UNIT=391,FILE="./StoreTADKMC/ProcessAssignments")
      prc=>KMC_PL
      nentries=0
      DO WHILE (ASSOCIATED(prc))
         WRITE(UNIT=391,FMT='(2I5," Atm:")',ADVANCE="NO") prc%Index,prc%NumberSubscribers
         AL=>prc%AL
         DO WHILE (ASSOCIATED(AL))
            nentries=nentries+prc%NumberProcessAtoms
            WRITE(UNIT=391,FMT='(I6)',ADVANCE="NO") AL%Atom%Index
            AL=>AL%NextNeigh
         END DO
         WRITE(391,*) " "
         prc=>prc%NextNeigh
      END DO
      WRITE(391,*) "Total entries:",nentries
      CLOSE(391)
   END SUBROUTINE PrintProcessAssignments
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintAtomProcessAssignments()
   !Prints the 
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCProcessSubscriptionInfo), POINTER :: prclist
      INTEGER :: iunit,nentries
      
      iunit=UnitScrap
      IF (iunit/=6 .AND. iunit/=UnitScrap) OPEN(UNIT=iunit,FILE="./StoreTADKMC/AtomProcessAssignments")
      atom=>KMC_AL
      nentries=0
      DO WHILE (ASSOCIATED(atom))
         IF (ASSOCIATED(atom%PrcSubscriberInfo)) THEN
            WRITE(UNIT=iunit,FMT='(I5)',ADVANCE="NO") atom%Index
            prclist=>atom%PrcSubscriberInfo
            DO WHILE (ASSOCIATED(prclist))
               nentries=nentries+1
               WRITE(UNIT=iunit,FMT='("[",I5,L1"]")',ADVANCE="NO") prclist%Process%Index,prclist%IsActive
               prclist=>prclist%NextNeigh
            END DO
            WRITE(iunit,*) " "
         END IF
         atom=>atom%NextNeigh
      END DO
      WRITE(iunit,*) "Total entries:",nentries
      IF (iunit/=6 .AND. iunit/=UnitScrap) CLOSE(iunit)
   END SUBROUTINE PrintAtomProcessAssignments
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintConfigurationAssignments()
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: cfg
      TYPE(KMCAtomList), POINTER :: AL
      TYPE(KMCAtom), POINTER :: atom
      INTEGER :: nentries,nentries1,nentries2
      
      cfg=>RecordedCfg
      nentries=0
      nentries1=0
      OPEN(UNIT=391,FILE="./StoreTADKMC/ConfigurationAssignments")
      DO WHILE (ASSOCIATED(cfg))
         IF (cfg%ConfigurationIndex>0 .AND. (ASSOCIATED(cfg%AtomsWithCfg) .OR. cfg%NAtomsWithCfg>0)) THEN
            WRITE(UNIT=391,FMT='(I5,"[",I5,"]")',ADVANCE="NO") cfg%ConfigurationIndex,cfg%NAtomsWithCfg
            nentries1=nentries1+cfg%NAtomsWithCfg
            AL=>cfg%AtomsWithCfg
            nentries2=0
            DO WHILE (ASSOCIATED(AL))
               nentries=nentries+1
               nentries2=nentries2+1
               IF (nentries2<=10) WRITE(UNIT=391,FMT='(I5)',ADVANCE="NO") AL%Atom%Index
               IF (nentries2==10) THEN
                  WRITE(UNIT=391,FMT='("...")',ADVANCE="NO")
               END IF
               AL=>AL%NextNeigh
            END DO
         END IF
         WRITE(391,*) " "
         cfg=>cfg%NextNeigh
      END DO
      WRITE(391,*) "Method 1: Number of entries",nentries
      WRITE(391,*) "Method 2: Number of entries",nentries1
      
      atom=>KMC_AL
      nentries=0
      DO WHILE (ASSOCIATED(atom))
         IF (atom%Cfg%ConfigurationIndex>0) nentries=nentries+1
         atom=>atom%NextNeigh
      END DO
      WRITE(391,*) "Method 3: Number of entries",nentries
      CLOSE(391)
   END SUBROUTINE PrintConfigurationAssignments
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintMemoryAssignments()
   !prints current KMC memory requirements
      IMPLICIT NONE
      LOGICAL, SAVE :: Initialized=.FALSE.
      
      IF (.NOT. Initialized) THEN
         OPEN(UNIT=392,FILE="./StoreTADKMC/MemoryAssignments")
         Initialized=.TRUE.
      END IF
      WRITE(392,*) "---------------------KMCSteps",KMCStep+NKMCStepsPerBlock*KMCBlock,"------"
      WRITE(392,*) "NAllottedCellAtom:",NAllottedCellAtom 
      WRITE(392,*) "NAllottedIsCellAffected:",NAllottedIsCellAffected
      WRITE(392,*) "NAllottedKMCAtom:",NAllottedKMCAtom
      WRITE(392,*) "NAllottedKMCProcessSubscriptionInfo:",NAllottedKMCProcessSubscriptionInfo
      WRITE(392,*) "NAllottedKMCProcess:",NAllottedKMCProcess
      WRITE(392,*) "NAllottedKMCProcessAtom:",NAllottedKMCProcessAtom
      WRITE(392,*) "NAllottedKMCProcessList:",NAllottedKMCProcessList
      WRITE(392,*) "NAllottedEnvironment:",NAllottedEnvironment
      WRITE(392,*) "NAllottedHistogram:",NAllottedHistogram
      WRITE(392,*) "NAllottedBin:",NAllottedBin
      WRITE(392,*) "NAllottedBinAtom:",NAllottedBinAtom
      WRITE(392,*) "NAllottedKMCAtomList:",NAllottedKMCAtomList
      WRITE(392,*) "NAllottedConfiguration:",NAllottedConfiguration
      WRITE(392,*) "NAllottedConfigurationList:",NAllottedConfigurationList
      WRITE(392,*) "NAllottedOthers:",NAllottedOthers
!      CLOSE(391)
      CALL FLUSH(392)
   END SUBROUTINE PrintMemoryAssignments
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintRecordedCfgSequence()
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: ptr
      
      ptr=>RecordedCfg
      DO WHILE (ASSOCIATED(ptr))
         WRITE(UNIT=6,FMT='(I4,"(",F10.3,"), ")',ADVANCE="NO") ptr%ConfigurationIndex,ptr%BondOrder
         ptr=>ptr%NextNeigh
      END DO
      WRITE(*,*) "."
   END SUBROUTINE PrintRecordedCfgSequence
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE KMCUtilities
