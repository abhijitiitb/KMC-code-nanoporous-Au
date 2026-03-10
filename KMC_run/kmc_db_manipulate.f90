MODULE KMCDbManipulate
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
!Perform structural changes to the DB

   USE KMC_VARIABLE_TYPE
   USE KMCUtilities
   IMPLICIT NONE
   INTEGER, PARAMETER :: MaxRadialBins=14
   
   INTERFACE SearchAtom
      MODULE PROCEDURE SearchAtomCoord,SearchAtomIndex
   END INTERFACE SearchAtom
   
   INTERFACE SearchCell
      MODULE PROCEDURE SearchCellAtom,SearchCellCoord
   END INTERFACE SearchCell
   
   CONTAINS
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddAtom(atom)
   !add a new atom after the current element atom, if atom is not associated that the current element is allocated
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
 
      IF (ASSOCIATED(atom)) THEN
         ALLOCATE(atom%NextNeigh); NAllottedKMCAtom=NAllottedKMCAtom+1
         atom%NextNeigh%PrevNeigh=>atom
         atom=>atom%NextNeigh
      ELSE
         ALLOCATE(atom); NAllottedKMCAtom=NAllottedKMCAtom+1
         NULLIFY(atom%PrevNeigh)
      END IF
      NULLIFY(atom%NextNeigh)

      ALLOCATE(atom%Env); NAllottedEnvironment=NAllottedEnvironment+1
      ALLOCATE(atom%Env%ShortRangeConfig); NAllottedHistogram=NAllottedHistogram+1
      NULLIFY(atom%Env%ShortRangeConfig%BinInfo)
      NAtomKMC=NAtomKMC+1
   END SUBROUTINE AddAtom
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddLocalCfg(Atom)
   !adds the local environment to an atom. assumed that environment has not been added before
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE (KMCAtom), POINTER :: Atom
      INTEGER :: RangeNeigh
      INTEGER :: Cellneigh(3),Cell(3)
      INTEGER :: ix,iy,iz
      
      RangeNeigh=NumberNeighborCellsInEnv
      Cell=Atom%CellIndex
      Atom%BondOrder=0.
      Atom%Env%ShortRangeConfig%Population=0
      
      DO ix=-RangeNeigh,RangeNeigh
         DO iy=-RangeNeigh,RangeNeigh
            DO iz=-RangeNeigh,RangeNeigh
               Cellneigh=GetNeighCell(Cell,ix,iy,iz)
               CALL AddHistogramFromCell(CellNeigh,Atom)
            END DO
         END DO
      END DO
      
   END SUBROUTINE AddLocalCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddHistogramFromCell(CellNeigh,Atom)
   !add the atoms in a cell lying the neighborhood of the center atom -- fine even if no neigh found
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      INTEGER :: BinNumber,BinRegion
      INTEGER :: AtomIndex
      INTEGER :: CellNeigh(3)
      REAL :: rad2,coord(3),rad,vec(3),BondOrder
      TYPE(KMCAtom), POINTER :: Atom,Atom1
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(Histogram), POINTER :: hist1
      TYPE(Bin), POINTER :: bin1
      TYPE(BinAtom), POINTER :: bincontent1
     
      AL1=>Cell(CellNeigh(1),CellNeigh(2),CellNeigh(3))%AL
      AtomIndex=Atom%Index
      hist1=>Atom%Env%ShortRangeConfig
      
      coord=Atom%Coord
      BondOrder=0.
      
      DO WHILE (ASSOCIATED(AL1)) !go through the atom list
         Atom1=>AL1%Atom
         IF (Atom1%Index/=AtomIndex) THEN
            vec=GetPBCSpacing(Atom1%Coord,coord)
            rad2=vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
            IF (rad2<Rad2Correction) THEN
               rad=SQRT(rad2)
               IF (rad<RadialTolerance) THEN
                  WRITE(6,*) 'Ini>>Radius is smaller than tolerance:',Rad
                  WRITE(6,*) 'Center atom index:',Atom%Index
                  WRITE(6,*) 'Center atom coordinate:',Atom%Coord
                  WRITE(6,*) 'Neighbor atom index:',Atom1%Index
                  WRITE(6,*) 'Neighbor atom coord:',Atom1%Coord
                  STOP
               END IF
               
               BondOrder=GetBondOrder(rad,Atom)
               Atom%BondOrder=Atom%BondOrder+BondOrder
               hist1%Population(Atom1%Species)=hist1%Population(Atom1%Species)+1
               
               CALL GetBinStatus(rad,BinRegion,BinNumber) !find where Atom1 would be present
               IF (ASSOCIATED(hist1%BinInfo)) THEN
                  bin1=>hist1%BinInfo
                  CALL SearchBinAfter(BinNumber,bin1) !finds the correct bin in best case, or else the neighbor
                  IF (BinNumber>bin1%BinNumber) THEN
                     CALL AddBinAfter(BinNumber,bin1)
                  ELSEIF (BinNumber<bin1%BinNumber) THEN
                     CALL AddBinBefore(BinNumber,bin1,Atom)
                  END IF
               ELSE
                  ALLOCATE(hist1%BinInfo); NAllottedBin=NAllottedBin+1
                  bin1=>hist1%BinInfo
                  bin1%BinNumber=BinNumber
               END IF
               ALLOCATE(bincontent1); NAllottedBinAtom=NAllottedBinAtom+1
               bincontent1%Atom=>Atom1
               bincontent1%Species=Atom1%Species
               bincontent1%RelCoord=vec
               bincontent1%distance=rad
               CALL RefreshBinEdge(bin1,BinRegion,1,INT(Atom1%Species)) !works on InBin,BinEdge
               CALL AddToBin(bincontent1,bin1)
            END IF
         END IF
         AL1=>AL1%NextNeigh
      END DO
      
   END SUBROUTINE AddHistogramFromCell
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddBinAfter(BinNumber,bin1)
   !add a new bin with index BinNumber after bin1
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Bin), POINTER :: bin1,bin2
      INTEGER :: BinNumber

      IF (.NOT. ASSOCIATED(bin1)) THEN
         WRITE(6,*) "Db>>AddBinAfter..Bin pointer is NULL"
         STOP
      END IF
      
      ALLOCATE(bin2); NAllottedBin=NAllottedBin+1
      bin2%BinNumber=BinNumber
      NULLIFY(bin2%AL)
      
      bin2%PrevNeigh=>bin1
      IF (ASSOCIATED(bin1%NextNeigh)) THEN
         bin2%NextNeigh=>bin1%NextNeigh
         bin1%NextNeigh%PrevNeigh=>bin2
      ELSE
         NULLIFY(bin2%NextNeigh)
      END IF
      bin1%NextNeigh=>bin2
      bin1=>bin2
   END SUBROUTINE AddBinAfter
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddBinBefore(BinNumber,bin1,Atom)
   !add a bin before bin1, Atom's bin gets modified
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: Atom
      TYPE(Bin), POINTER :: bin1,bin2
      INTEGER :: BinNumber

      NULLIFY(bin2)
      ALLOCATE(bin2); NAllottedBin=NAllottedBin+1
      bin2%BinNumber=BinNumber
      NULLIFY(bin2%AL)
      
      IF (ASSOCIATED(bin1)) THEN
         bin2%NextNeigh=>bin1
         IF (ASSOCIATED(bin1%PrevNeigh)) THEN
            bin2%PrevNeigh=>bin1%PrevNeigh
            bin1%PrevNeigh%NextNeigh=>bin2
         ELSE
            NULLIFY(bin2%PrevNeigh)
            Atom%Env%ShortRangeConfig%BinInfo=>bin2
         END IF
         bin1%PrevNeigh=>bin2
      ELSE !insert empty element
         IF (.NOT. ASSOCIATED(Atom%Env%ShortRangeConfig%BinInfo)) THEN 
            Atom%Env%ShortRangeConfig%BinInfo=>bin2
         END IF
      END IF
      bin1=>bin2
   END SUBROUTINE AddBinBefore
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddToBin(binatom1,bin1)
   !Adds the binatom1 to bin1
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Bin), POINTER :: bin1
      TYPE(BinAtom), POINTER :: binatom1

      IF (ASSOCIATED(bin1%AL)) THEN
         binatom1%NextNeigh=>bin1%AL
         bin1%AL%PrevNeigh=>binatom1
         bin1%AL=>binatom1
      ELSE
         bin1%AL=>binatom1
         NULLIFY(binatom1%NextNeigh)
      END IF
      NULLIFY(binatom1%PrevNeigh)
   END SUBROUTINE AddToBin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !SUBROUTINE AddFromRedundantBin(atom,binatom1)
   ! !add from site atom's list of redundant bincontents
   !   IMPLICIT NONE
   !   TYPE(KMCAtom), POINTER :: atom
   !   TYPE(BinAtom), POINTER :: binatom1
   !   TYPE(Histogram), POINTER :: hist1

   !   hist1=>atom%Env%ShortRangeConfig
   !   IF (ASSOCIATED(hist1%RedundantNeighbors)) THEN
   !      binatom1=>hist1%RedundantNeighbors
   !      IF (ASSOCIATED(binatom1%NextNeigh)) THEN
   !         hist1%RedundantNeighbors=>binatom1%NextNeigh
   !         NULLIFY(hist1%RedundantNeighbors%PrevNeigh)
   !      ELSE
   !         NULLIFY(hist1%RedundantNeighbors)
   !      END IF
   !   ELSE
   !      ALLOCATE(binatom1); NAllottedBinAtom=NAllottedBinAtom+1
   !   END IF
   !   NULLIFY(binatom1%PrevNeigh)
   !   NULLIFY(binatom1%NextNeigh)
   !END SUBROUTINE AddFromRedundantBin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddFromRedundantCfg(cfg)
    !add from list of redundant configs
    !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: cfg

      ALLOCATE(cfg); NAllottedConfiguration=NAllottedConfiguration+1
   END SUBROUTINE AddFromRedundantCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddCfgToRecordedCfg(cfg1) !cfg1 should already be filled up
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: cfg,cfg1,cfg2,cfg3 !cfg is the proposed place where the new config cfg1 will be placed. cfg2 is a copy of cfg and cfg3 is a copy of neighbor of cfg2
      
      !NumberCfgTypes=NumberCfgTypes+1
      
      !insert configuration cfg1 so that the BondOrder is in ascending order
      cfg2=>RecordedCfg
      
      IF (cfg2%ConfigurationIndex<0) THEN !add cfg1 to the beginning
         CALL AddToRecordedCfgBefore(cfg2,cfg1)
         RETURN
      END IF
      
      IF (cfg2%BondOrder>cfg1%BondOrder .AND. cfg2%ConfigurationIndex>0) THEN
         CALL AddToRecordedCfgBefore(cfg2,cfg1)
      ELSEIF (cfg2%BondOrder<=cfg1%BondOrder .AND. cfg2%ConfigurationIndex>0) THEN !move to right
         cfg3=>cfg2%NextNeigh
         DO WHILE (ASSOCIATED(cfg3))
            IF (cfg3%BondOrder>cfg1%BondOrder .OR. cfg3%ConfigurationIndex<0) EXIT !cfg1 should come before cfg3
            cfg2=>cfg3
            cfg3=>cfg3%NextNeigh
         END DO
         CALL AddToRecordedCfgAfter(cfg2,cfg1)
      END IF
      !CALL PrintRecordedCfgSequence()
   END SUBROUTINE AddCfgToRecordedCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddAtomToRecordedCfg(atom,cfg,EndOfList)
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(Configuration), POINTER :: cfg,cfg1,cfg2,cfg3 !cfg is the proposed place where the new config cfg1 will be placed. cfg2 is a copy of cfg and cfg3 is a copy of neighbor of cfg2
      LOGICAL :: EndOfList
      
      NumberCfgTypes=NumberCfgTypes+1
      !CALL AddFromRedundantCfg(cfg1)
      ALLOCATE(cfg1); NAllottedConfiguration=NAllottedConfiguration+1
      cfg1%ConfigurationIndex=NumberCfgTypes
      CALL CopyAtom2Cfg(atom,cfg1) !cfg1 is the new configuration created
      CALL AddCfgToNewList(cfg1)
      
      IF (EndOfList) THEN
         CALL AddToRecordedCfgBefore(cfg,cfg1) !add cfg1 before cfg
         RETURN
      ELSEIF (.NOT. ASSOCIATED(RecordedCfg)) THEN
         RecordedCfg=>cfg1
         RETURN
      END IF
      
      !insert configuration cfg1 so that the BondOrder is in ascending order
      cfg2=>cfg
      IF (cfg2%BondOrder>cfg1%BondOrder) THEN !move to left
         cfg3=>cfg2%PrevNeigh
         DO WHILE (ASSOCIATED(cfg3))
            IF (cfg3%BondOrder<=cfg1%BondOrder) EXIT !cfg1 should come after cfg3
            cfg2=>cfg3
            cfg3=>cfg3%PrevNeigh
         END DO
         CALL AddToRecordedCfgBefore(cfg2,cfg1)
      ELSEIF (cfg2%BondOrder<cfg1%BondOrder) THEN !move to right
         cfg3=>cfg2%NextNeigh
         DO WHILE (ASSOCIATED(cfg3))
            IF (cfg3%BondOrder>=cfg1%BondOrder) EXIT !cfg1 should come before cfg3
            cfg2=>cfg3
            cfg3=>cfg3%NextNeigh
         END DO
         CALL AddToRecordedCfgAfter(cfg2,cfg1)
      END IF
      cfg=>cfg1
      !CALL PrintRecordedCfgSequence()
   END SUBROUTINE AddAtomToRecordedCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddToRecordedCfgAfter(cfg,cfg1) !add cfg1 after cfg
   !Add a configurationtype after configuration cfg
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: cfg,cfg1

      cfg1%PrevNeigh=>cfg
      IF (ASSOCIATED(cfg%NextNeigh)) THEN
         cfg%NextNeigh%PrevNeigh=>cfg1
         cfg1%NextNeigh=>cfg%NextNeigh
      ELSE
         NULLIFY(cfg1%NextNeigh)
      END IF
      cfg%NextNeigh=>cfg1
   END SUBROUTINE AddToRecordedCfgAfter
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddToRecordedCfgBefore(cfg,cfg1) !add cfg1 before cfg
   !Add new configuration type before configuration type cfg
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: cfg,cfg1

      cfg1%NextNeigh=>cfg
      IF (ASSOCIATED(cfg%PrevNeigh)) THEN
         cfg%PrevNeigh%NextNeigh=>cfg1
         cfg1%PrevNeigh=>cfg%PrevNeigh
      ELSE
         NULLIFY(cfg1%PrevNeigh)
         RecordedCfg=>cfg1
      END IF
      cfg%PrevNeigh=>cfg1
   END SUBROUTINE AddToRecordedCfgBefore
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ResetAffectedCellInfo()
      IMPLICIT NONE
      
      NAffectedCells=0
      IsCellAffected=0
   END SUBROUTINE ResetAffectedCellInfo
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddCellToAffectedList(CellIndex,value)
   !adds CellIndex to the list of affected cells while manually creating a process
   !Last checked: Feb 12, 2011
   
   !value=0: Cell is not affected
   
   !LEKMC+NEB
   !value=1: Cell is affected, some of the atoms in the cell could be moving
   !value=2: Cell is not directly affected, however, the atoms are allowed to move so that the value=1 cells can relax correctly to their min energy state without any residual stresses. Upon optimization, the coordinates of these atoms are not updated
   !value=3: Cell is not directly affected, all atoms in cell are not moving -- this is useful for energy optimization
   
   !LEKMC+MD
   !value=1: Domain region
   !value=2: neighbor domain region
   !value=3: moving buffer region
   !value=4: non-moving buffer region
      IMPLICIT NONE
      INTEGER :: i,ix,iy,iz,CellIndex(3),value
      
      ix=CellIndex(1)
      iy=CellIndex(2)
      iz=CellIndex(3)
      IF (IsCellAffected(ix,iy,iz)==0) THEN
         IsCellAffected(ix,iy,iz)=value
         NAffectedCells=NAffectedCells+1  
         AffectedCells(3*NAffectedCells-2:3*NAffectedCells)=CellIndex 
         IF (NAffectedCells>SizeAffectedCells) THEN
            WRITE(6,*) "Total number of cells:",NCells
            WRITE(6,*) "Number of affected cells:",NAffectedCells
            WRITE(6,*) "Maximum allowed affected cells:",SizeAffectedCells
            STOP
         END IF
      END IF
   END SUBROUTINE AddCellToAffectedList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddCellToAffectedList1(CellIndex)
   !adds CellIndex to the list of affected cells while manually creating a process
   !Last checked: Jan 12, 2011
      IMPLICIT NONE
      INTEGER :: i,CellIndex(3)
      
      DO i=1,NAffectedCells
         IF (ALL(AffectedCells(3*i-2:3*i)==CellIndex)) RETURN
      END DO
      
      NAffectedCells=NAffectedCells+1
      AffectedCells(3*NAffectedCells-2:3*NAffectedCells)=CellIndex
   END SUBROUTINE AddCellToAffectedList1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddToCfg(config,Atom,ptr)
   !adds atom to config%AL
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: config
      TYPE(KMCAtom), POINTER :: Atom
      TYPE(KMCAtomList), POINTER :: ptr

      config%NAtomsWithCfg=config%NAtomsWithCfg+1
      IF (ASSOCIATED(config%AtomsWithCfg)) THEN
         config%AtomsWithCfg%PrevNeigh=>ptr
         ptr%NextNeigh=>config%AtomsWithCfg
         config%AtomsWithCfg=>ptr
         NULLIFY(ptr%PrevNeigh)
      ELSE
         config%AtomsWithCfg=>ptr
         NULLIFY(ptr%PrevNeigh)
         NULLIFY(ptr%NextNeigh)
      END IF
      ptr%Atom=>Atom
      Atom%Cfg=>config
      !CALL AddCfgToNewList(config) !this has been commented to reduce the amont of file writing. As a result, now the latest atoms connected to a configuration are not written
   END SUBROUTINE AddToCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !SUBROUTINE AddToRedundantBin(atom,binatom1)
   ! !add bincontent binatom1 to list of redundants of site atom
   !   IMPLICIT NONE
   !   TYPE(KMCAtom), POINTER :: atom
   !   TYPE(BinAtom), POINTER :: binatom1
   !   TYPE(Histogram), POINTER :: hist1

   !   NULLIFY(binatom1%Atom)
   !   NULLIFY(binatom1%PrevNeigh)
   !   NULLIFY(binatom1%NextNeigh)
   !   binatom1%Species=-1
   !   binatom1%distance=-10.
   !   binatom1%RelCoord=-10.

   !   hist1=>atom%Env%ShortRangeConfig
   !   IF (ASSOCIATED(hist1%RedundantNeighbors)) THEN
   !      hist1%RedundantNeighbors%PrevNeigh=>binatom1
   !      binatom1%NextNeigh=>hist1%RedundantNeighbors
   !   END IF
   !   hist1%RedundantNeighbors=>binatom1
   !END SUBROUTINE AddToRedundantBin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddToCell(AL1,CellIndex)
   !Add element AL containing an atom to CellIndex
   !last checked - Jan 08, 2011
      IMPLICIT NONE
      TYPE (KMCAtomList), POINTER :: AL1
      INTEGER :: CellIndex(3),ix,iy,iz

      IF (.NOT. ASSOCIATED(AL1)) THEN
         WRITE(UnitScrap,*) "Db>>AddToCellList...AL provided is empty"
         STOP
      END IF

      ix=CellIndex(1)
      iy=CellIndex(2)
      iz=CellIndex(3)
      IF (ASSOCIATED(Cell(ix,iy,iz)%AL)) THEN
         Cell(ix,iy,iz)%AL%PrevNeigh=>AL1
         AL1%NextNeigh=>Cell(ix,iy,iz)%AL
         Cell(ix,iy,iz)%AL=>AL1
      ELSE
         Cell(ix,iy,iz)%AL=>AL1
         NULLIFY(AL1%PrevNeigh)
         NULLIFY(AL1%NextNeigh)
      END IF
      Cell(ix,iy,iz)%NAtomsInCell=Cell(ix,iy,iz)%NAtomsInCell+1
      IF (AL1%Atom%IsMoving) THEN
         Cell(ix,iy,iz)%NMoveInCell=Cell(ix,iy,iz)%NMoveInCell+1
      END IF
   END SUBROUTINE AddToCell
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddCfgToNewList(Config)
   !add config to the list which says one of its structure elements has been modified
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: Config
      LOGICAL :: IsConfigExist
   
      IF (Config%ConfigurationIndex<0) RETURN !Redundant config is never printed
      IsConfigExist=SearchCfgFromNewList(Config)
      IF (IsConfigExist) RETURN !Config is already present

      IF (ASSOCIATED(NewCfgList)) THEN
         ALLOCATE(NewCfgList%PrevNeigh); NAllottedConfigurationList=NAllottedConfigurationList+1
         NewCfgList%PrevNeigh%NextNeigh=>NewCfgList
         NewCfgList=>NewCfgList%PrevNeigh
         NULLIFY(NewCfgList%PrevNeigh)
      ELSE
         ALLOCATE(NewCfgList); NAllottedConfigurationList=NAllottedConfigurationList+1
         NULLIFY(NewCfgList%PrevNeigh)
         NULLIFY(NewCfgList%NextNeigh)
      END IF
      NewCfgList%Cfg=>Config
   END SUBROUTINE AddCfgToNewList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddProcessToNewList(prc)
   !add process to NewProcessList only if new process is found
   !Last checked: Jan 12 , 2011
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      LOGICAL :: IsProcessExist
      TYPE(KMCProcessList), POINTER :: prclist
      
      IsProcessExist=SearchProcessFromNewList(prc)
      IF (IsProcessExist) RETURN
     
      IF (prc%Index<=0) THEN
         WRITE(6,*) "Err>> Process index is zero in AddProcessToNewList"
         STOP
      END IF
      IF (ASSOCIATED(NewProcessList)) THEN
         ALLOCATE(NewProcessList%PrevNeigh); NAllottedKMCProcessList=NAllottedKMCProcessList+1
         NewProcessList%PrevNeigh%NextNeigh=>NewProcessList
         NewProcessList=>NewProcessList%PrevNeigh
         NULLIFY(NewProcessList%PrevNeigh)
      ELSE
         ALLOCATE(NewProcessList); NAllottedKMCProcessList=NAllottedKMCProcessList+1
         NULLIFY(NewProcessList%PrevNeigh)
         NULLIFY(NewProcessList%NextNeigh)
      END IF
      NewProcessList%Process=>prc
   END SUBROUTINE AddProcessToNewList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddProcessToActiveList(prc)
   !add process (which was inactive) to ActiveProcessList only if the process become active
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      LOGICAL :: IsProcessExist
      TYPE(KMCProcessList), POINTER :: prclist
      
      NumberActiveProcesses=NumberActiveProcesses+1
      IF (ASSOCIATED(ActiveProcessList)) THEN
         ALLOCATE(ActiveProcessList%PrevNeigh); NAllottedKMCProcessList=NAllottedKMCProcessList+1
         ActiveProcessList%PrevNeigh%NextNeigh=>ActiveProcessList
         ActiveProcessList=>ActiveProcessList%PrevNeigh
         NULLIFY(ActiveProcessList%PrevNeigh)
      ELSE
         ALLOCATE(ActiveProcessList); NAllottedKMCProcessList=NAllottedKMCProcessList+1
         NULLIFY(ActiveProcessList%PrevNeigh)
         NULLIFY(ActiveProcessList%NextNeigh)
      END IF
      ActiveProcessList%Process=>prc
      prc%APLPosition=>ActiveProcessList !this tells us the exact position of the process in the ActiveProcessList -- useful in case we need to remove the process from ActiveProcessList
   END SUBROUTINE AddProcessToActiveList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddProcessToCfgType(prc)
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(Configuration), POINTER :: cfg1
      TYPE(KMCProcessList), POINTER :: PL1
      
      cfg1=>prc%CfgType
      IF (ASSOCIATED(cfg1%PL)) THEN
         PL1=>cfg1%PL
         ALLOCATE(PL1%PrevNeigh); NAllottedKMCProcessList=NAllottedKMCProcessList+1
         PL1=>PL1%PrevNeigh
         PL1%NextNeigh=>cfg1%PL
         cfg1%PL=>PL1
      ELSE
         ALLOCATE(cfg1%PL); NAllottedKMCProcessList=NAllottedKMCProcessList+1
      END IF
      cfg1%PL%Process=>prc
   END SUBROUTINE AddProcessToCfgType
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddProcessToPL(prc)
   !adds a new process to PL
   !Last checked: Jan 12, 2011
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      
      NumberProcesses=NumberProcesses+1
      prc%Index=NumberProcesses
      prc%RecordNumber=NumberProcesses
      prc%FirstAdded=KMCStep+KMCBlock*NKMCStepsPerBlock
      
      IF (ASSOCIATED(KMC_PL)) THEN
         KMC_PL%PrevNeigh=>prc
         prc%NextNeigh=>KMC_PL
      END IF
      KMC_PL=>prc
      CALL AddProcessToNewList(prc)
      CALL AddProcessToCfgType(prc) !adds process to configuration of type CfgType
   END SUBROUTINE AddProcessToPL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddProcessToAtom(prc,atom,IsActive,ActiveAtom)
   !add the information that the process is activated due to this atom
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessSubscriptionInfo), POINTER :: ptr
      TYPE(KMCAtom), POINTER :: atom,ActiveAtom
      LOGICAL :: IsActive
      
      ALLOCATE(ptr); NAllottedKMCProcessSubscriptionInfo=NAllottedKMCProcessSubscriptionInfo+1
      ptr%Process=>prc
      ptr%IsActive=IsActive
      ptr%ActiveAtom=>ActiveAtom
      IF (ASSOCIATED(atom%PrcSubscriberInfo)) THEN
         atom%PrcSubscriberInfo%PrevNeigh=>ptr
         ptr%NextNeigh=>atom%PrcSubscriberInfo
      END IF
      atom%PrcSubscriberInfo=>ptr
   END SUBROUTINE AddProcessToAtom
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddAtomToProcess(prc,relinitialcoord,relfinalcoord,relsaddlecoord,cfgindex,species)
   !adds the details of a reactant atom to a process
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: prcatom
      REAL, DIMENSION(3) :: relinitialcoord,relfinalcoord,relsaddlecoord
      INTEGER :: cfgindex,species
      
      IF (.NOT. ASSOCIATED(prc)) THEN
         WRITE(6,*) "Process has not been allocated any space"
         STOP
      END IF
      
      ALLOCATE(prcatom); NAllottedKMCProcessAtom=NAllottedKMCProcessAtom+1
      prcatom%RelInitialPos=relinitialcoord
      prcatom%RelSaddlePos=relsaddlecoord
      prcatom%RelFinalPos=relfinalcoord
      prcatom%CfgIndex=cfgindex
      prcatom%Species=species
      
      prc%NumberProcessAtoms=prc%NumberProcessAtoms+1
      
      NULLIFY(prcatom%NextNeigh)
      IF (ASSOCIATED(prc%ReactantAL)) THEN
         prc%ReactantAL%PrevNeigh=>prcatom
         prcatom%NextNeigh=>prc%ReactantAL
         prc%ReactantAL=>prcatom
      ELSE
         prc%ReactantAL=>prcatom
      END IF
      NULLIFY(prcatom%PrevNeigh)
   END SUBROUTINE AddAtomToProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AddParticipatingAtomToProcess(atom,prc)
   !adds the atom to prc%AL -- the atom should be active atom of process not the passive one
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCAtomList), POINTER :: ptr
      TYPE(KMCProcess), POINTER :: prc
      
      ALLOCATE(ptr); NAllottedKMCAtomList=NAllottedKMCAtomList+1
      ptr%Atom=>atom
      
      IF(prc%NumberSubscribers==0) CALL AddProcessToActiveList(prc)
      prc%NumberSubscribers=prc%NumberSubscribers+1
      
      IF (ASSOCIATED(prc%AL)) THEN
         prc%AL%PrevNeigh=>ptr
         ptr%NextNeigh=>prc%AL
      END IF
      prc%AL=>ptr
   END SUBROUTINE AddParticipatingAtomToProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AssignProcessFull(prc)
   !Checks all atoms in the system to assign processes to them
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom1,activeatom,atom
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: prcatom
      REAL, DIMENSION(3) :: coord,coord1
      INTEGER :: ActiveAtomIndex,cfgindex1
      REAL :: tol
      LOGICAL :: NoMatch,IsActive
      
      tol=RadialTolerance
      cfgindex1=prc%CfgType%ConfigurationIndex
      atom=>KMC_AL
      DO WHILE (ASSOCIATED(atom))
         NoMatch=.TRUE.
         IF (atom%Cfg%ConfigurationIndex==cfgindex1) THEN
            coord=atom%Coord
            activeatom=>atom
            prcatom=>prc%ReactantAL
            DO WHILE (ASSOCIATED(prcatom))
               coord1=GetPBCCoord(coord+prcatom%RelInitialpos)
               atom1=>SearchAtom(coord1,tol)
               NoMatch=.NOT. ASSOCIATED(atom1)
               IF (NoMatch) EXIT
               NoMatch=atom1%Cfg%ConfigurationIndex/=prcatom%CfgIndex
               IF (NoMatch) EXIT
               prcatom=>prcatom%NextNeigh
            END DO
         END IF
         IF (.NOT. NoMatch .AND. ProcessNotAssigned(activeatom,prc)) THEN !this process should be added
            ActiveAtomIndex=activeatom%Index
            prcatom=>prc%ReactantAL
            DO WHILE (ASSOCIATED(prcatom))
               coord1=GetPBCCoord(coord+prcatom%RelInitialpos)
               atom1=>SearchAtom(coord1,tol)
               IsActive=(atom1%Index==ActiveAtomIndex)
               CALL AddProcessToAtom(prc,atom1,IsActive,activeatom)
               prcatom=>prcatom%NextNeigh
            END DO
            CALL AddParticipatingAtomToProcess(activeatom,prc)
         END IF
         atom=>atom%NextNeigh
      END DO
      
   END SUBROUTINE AssignProcessFull
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AssignProcessInAllAffectedCells()
   !this is a faster version of than the one implemented in kmc_algorithm.f90
      IMPLICIT NONE
      INTEGER :: i,CellIndex(3),ActiveAtomIndex
      TYPE(KMCAtom), POINTER :: atom1,activeatom
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(KMCProcessList), POINTER :: PL1
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: prcatom
      REAL :: tol,coord(3),coord1(3)
      LOGICAL :: NoMatch,IsActive
      
      tol=RadialTolerance
!WRITE(UnitScrap,*) "Assigning processes in affected cells @ KMC step:",KMCStep+NKMCStepsPerBlock*KMCBlock
!WRITE(UnitScrap,*) "NUmber of affected cells while assigning process:",nAffectedCells
      DO i=1,NAffectedCells
         CellIndex=AffectedCells(3*i-2:3*i)
         AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
         DO WHILE (ASSOCIATED(AL1))
            activeatom=>AL1%Atom
            
!IF (activeatom%Index<4)WRITE(UnitScrap,*) "Assigning process to atom:",activeatom%Index

            PL1=>activeatom%Cfg%PL
            coord=activeatom%Coord
            DO WHILE (ASSOCIATED(PL1))
               prc=>PL1%Process
      !WRITE(6,*) "Searching participating atoms for process index:",prc%Index
      !WRITE(6,*) " ...currently # participant atoms:",prc%NumberSubscribers
      !WRITE(6,*) "Attempting to find process atoms for active coord:",coord
               prcatom=>prc%ReactantAL
               NoMatch=.TRUE.
               DO WHILE (ASSOCIATED(prcatom))
                  coord1=GetPBCCoord(coord+prcatom%RelInitialpos)
                  atom1=>SearchAtom(coord1,tol)
                  NoMatch=.NOT. ASSOCIATED(atom1)
      !WRITE(6,*) "Searched coord.",coord1,ASSOCIATED(atom1)
                  IF (NoMatch) EXIT
                  NoMatch=atom1%Cfg%ConfigurationIndex/=prcatom%CfgIndex
                  IF (NoMatch) EXIT
                  prcatom=>prcatom%NextNeigh
               END DO
               IF (.NOT. NoMatch .AND. ProcessNotAssigned(activeatom,prc)) THEN !this process should be added
                  ActiveAtomIndex=activeatom%Index
      !WRITE(6,*) " ...found new participant atoms:",ActiveAtomIndex,activeatom%Coord
                  prcatom=>prc%ReactantAL
                  DO WHILE (ASSOCIATED(prcatom))
                     coord1=GetPBCCoord(coord+prcatom%RelInitialpos)
                     atom1=>SearchAtom(coord1,tol)
                     IsActive=(atom1%Index==ActiveAtomIndex)
                     CALL AddProcessToAtom(prc,atom1,IsActive,activeatom)
                     prcatom=>prcatom%NextNeigh
                  END DO
                  CALL AddParticipatingAtomToProcess(activeatom,prc)
               END IF
               PL1=>PL1%NextNeigh
      !WRITE(6,*) "END OF PARTICIPANT ATOMS"
            END DO
            AL1=>AL1%NextNeigh
         END DO
      END DO
   END SUBROUTINE AssignProcessInAllAffectedCells
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE AssignProcessInAffectedCells(prc)
   !Used in main KMC loop to find new processes in affected cells
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom1,activeatom,atom
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: prcatom
      REAL, DIMENSION(3) :: coord,coord1
      REAL :: tol
      LOGICAL :: NoMatch,IsActive
      INTEGER :: CellIndex(3),i,ActiveAtomIndex,cfgindex1
      
      WRITE(6,*) "Searching participating atoms for process index:",prc%Index
      WRITE(6,*) " ...currently # participant atoms:",prc%NumberSubscribers
      tol=RadialTolerance
      cfgindex1=prc%CfgType%ConfigurationIndex
      DO i=1,NAffectedCells
         CellIndex=AffectedCells(3*i-2:3*i)
         AL1=>Cell(CellIndex(1),CellIndex(2),CellIndex(3))%AL
         DO WHILE (ASSOCIATED(AL1))
            atom=>AL1%Atom
            NoMatch=.TRUE.
            IF (atom%Cfg%ConfigurationIndex==cfgindex1) THEN
               coord=atom%Coord
               activeatom=>atom
               prcatom=>prc%ReactantAL
               DO WHILE (ASSOCIATED(prcatom))
                  coord1=GetPBCCoord(coord+prcatom%RelInitialpos)
                  atom1=>SearchAtom(coord1,tol)
                  NoMatch=.NOT. ASSOCIATED(atom1)
                  IF (NoMatch) EXIT
                  NoMatch=atom1%Cfg%ConfigurationIndex/=prcatom%CfgIndex
                  IF (NoMatch) EXIT
                  prcatom=>prcatom%NextNeigh
               END DO
            END IF
            IF (.NOT. NoMatch .AND. ProcessNotAssigned(activeatom,prc)) THEN !this process should be added
               WRITE(6,*) " ...found new participant atoms:",ActiveAtomIndex,activeatom%Coord
               ActiveAtomIndex=activeatom%Index
               prcatom=>prc%ReactantAL
               DO WHILE (ASSOCIATED(prcatom))
                  coord1=GetPBCCoord(coord+prcatom%RelInitialpos)
                  atom1=>SearchAtom(coord1,tol)
                  IsActive=(atom1%Index==ActiveAtomIndex)
                  CALL AddProcessToAtom(prc,atom1,IsActive,activeatom)
                  prcatom=>prcatom%NextNeigh
               END DO
               CALL AddParticipatingAtomToProcess(activeatom,prc)
            END IF
            AL1=>AL1%NextNeigh
         END DO
      END DO
      WRITE(6,*) "END OF PARTICIPANT ATOMS"
   END SUBROUTINE AssignProcessInAffectedCells
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ConfirmCellIsAffected(CellIndex)
   !confirms that CellIndex is one of the affected cells
      IMPLICIT NONE
      INTEGER :: i,CellIndex(3)
      LOGICAL :: Found
      
      Found=.FALSE.
      DO i=1,NAffectedCells
         Found=ALL(AffectedCells(3*i-2:3*i)==CellIndex)
         IF (Found) EXIT
      END DO
      
      IF (.NOT. Found) THEN
         WRITE(6,*) "CellIndex is not in the list of affected cells"
         WRITE(6,*) "Cell index:",CellIndex
         STOP
      END IF
   END SUBROUTINE ConfirmCellIsAffected
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CheckDisplacement(atom,displacement,tol)
   !check whether an atom if it is displaced will overlap with an existing atom
   !Last checked: Jan 12, 2010
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom,atom1
      LOGICAL :: CheckDisplacement
      REAL :: displacement(3),coord(3),tol
      
      CheckDisplacement=.FALSE.
      coord=GetPBCCoord(atom%Coord+displacement)
      atom1=>SearchAtom(coord,tol)
      CheckDisplacement= .NOT. ASSOCIATED(atom1)
   END FUNCTION CheckDisplacement
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckKMC()
   !performs a number of tests to check if the current KMC database has consistency
   !with each other.
   !Motivation behind this subroutine:
   !The code contains following data structure which share the same piece of information
   !a) If atom participates in process as an active member then process contains it as a subscriber
   !b) If process is active then its number of subscribers must be greater than zero and it should be present in ActiveProcessList
   !c) If process is of type cfg then cfg should contain the process in its list
      IMPLICIT NONE
      
      CALL CheckSumOfRates()
      CALL CheckCellInfo()
      CALL CheckAllLocalEnvironments()
      CALL CheckAllConfigs()
      CALL CheckSubscriptionInfoAllAtoms()
      CALL CheckSubscriberAllProcesses()
      CALL CheckActiveProcessList()
      CALL CheckProcessAPLPointer()
      
      WRITE(6,*) "Checked KMC data structures - no error found"
   END SUBROUTINE CheckKMC
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckSumOfRates()
   END SUBROUTINE CheckSumOfRates
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckCellInfo()
      !checks if each atom has the correct cell location and if the atom is present in the corresponding cell
      IMPLICIT NONE
      INTEGER :: ix,iy,iz,AtomIndex
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(KMCAtom), POINTER :: atom
      LOGICAL :: Found
      
      !check for all cells
      DO iz=1,NCells(3)
         DO iy=1,NCells(2)
            DO ix=1,NCells(1)
               AL1=>Cell(ix,iy,iz)%AL
               DO WHILE (ASSOCIATED(AL1))
                  AtomIndex=AL1%Atom%Index
                  atom=>KMC_AL
                  DO WHILE (ASSOCIATED(atom))
                     IF (atom%Index==AtomIndex) EXIT
                     atom=>atom%NextNeigh
                  END DO
                  IF (ANY(atom%CellIndex/=(/ix,iy,iz/))) THEN
                     WRITE(6,*) "Check>> Cell location of atom is incorrect"
                     WRITE(6,*) "Offending atom:",AtomIndex
                     WRITE(6,*) "Atom coordinate:",Atom%Coord
                     WRITE(6,*) "Cell containing atom:",ix,iy,iz
                     WRITE(6,*) "Cell location stored in atom:",atom%CellIndex
                     STOP
                  END IF
                  AL1=>AL1%NextNeigh
               END DO
            END DO
         END DO
      END DO
      
      !check for all atoms
      atom=>KMC_AL
      DO WHILE (ASSOCIATED(atom))
         AtomIndex=atom%Index
         ix=atom%CellIndex(1)
         iy=atom%CellIndex(2)
         iz=atom%CellIndex(3)
         AL1=>Cell(ix,iy,iz)%AL
         Found=.FALSE.
         DO WHILE (ASSOCIATED(AL1))
            Found=AL1%Atom%Index==AtomIndex
            IF (Found) EXIT
            AL1=>AL1%NextNeigh
         END DO
         IF (.NOT. Found) THEN
            WRITE(6,*) "Check>> Cell location of atom is incorrect"
            WRITE(6,*) "Offending atom:",AtomIndex
            WRITE(6,*) "Atom coordinate:",Atom%Coord
            WRITE(6,*) "Cell containing atom:",ix,iy,iz
            WRITE(6,*) "Cell location stored in atom:",atom%CellIndex
            STOP
         END IF
         atom=>atom%NextNeigh
      END DO
   END SUBROUTINE CheckCellInfo
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckAllLocalEnvironments()
      !makes sure that each atom has the correct local environment and bond order
      IMPLICIT NONE
   END SUBROUTINE CheckAllLocalEnvironments
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckAllConfigs()
      !makes sure that 
      !a) each atom has the correct configuration by comparing atom configuration to local environment
      !b) the data stored in AtomsWithCfg for each configuration is correct
   END SUBROUTINE CheckAllConfigs
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckSubscriptionInfoAllAtoms()
      !checks two types of information -- for all atoms
      !a) if atom is subscriber to a process (as an ActiveAtom) then process should contain that atom as the ActiveAtom in AL
      !b) the configuration of the active atom and the process should match
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCProcessSubscriptionInfo), POINTER :: subscribedprc
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCAtomList), POINTER :: prcatom
      INTEGER :: AtomIndex,NSubscribedAsActive
      LOGICAL :: Found
      
!      IF (kmcstep+nkmcstepsperblock*kmcblock>0) THEN
!      DO AtomIndex=1,3
!         atom=>SearchAtomIndex(AtomIndex)
!         subscribedprc=>atom%PrcSubscriberInfo
!         NSubscribedAsActive=0
!         WRITE(UnitScrap,*) "Processes of atom ",AtomIndex," with cfg.",atom%Cfg%ConfigurationIndex
!         WRITE(UnitScrap,*) "Coord.",atom%Coord
!         DO WHILE (ASSOCIATED(subscribedprc))
!            WRITE(UnitScrap,*) subscribedprc%Process%Index,subscribedprc%IsActive
!            IF (subscribedprc%IsActive) NSubscribedAsActive=NSubscribedAsActive+1
!            subscribedprc=>subscribedprc%NextNeigh
!         END DO
!         IF (NSubscribedAsActive==0) THEN
!            WRITE(UnitScrap,*) "Following atoms does not have any processes:",AtomIndex
!            CALL PrintProcessAssignments()
!            CALL PrintAtomProcessAssignments()
!STOP
!         END IF
!      END DO
!      END IF
      
      atom=>KMC_AL
      DO WHILE (ASSOCIATED(atom))
         !check if an atom with a redundant configuration is subscribing to a process
         IF (ASSOCIATED(atom%PrcSubscriberInfo) .AND. atom%Cfg%ConfigurationIndex<0) THEN
            WRITE(6,*) "Check>> Atom with redundant cfg is subscribing to a process"
            WRITE(6,*) "Offending atom:",atom%Index
            WRITE(6,*) "Atom coordinate:",atom%Coord
            WRITE(6,*) "List of process found:"
            subscribedprc=>atom%PrcSubscriberInfo
            DO WHILE (ASSOCIATED(subscribedprc))
               WRITE(6,*) "Process index:",subscribedprc%Process%Index, &
                 " Active:",subscribedprc%IsActive
               subscribedprc=>subscribedprc%NextNeigh
            END DO
            STOP
         END IF
         !check if the processes atom is subscribing to, contains the atom as a participant
         subscribedprc=>atom%PrcSubscriberInfo
         AtomIndex=atom%Index
         DO WHILE(ASSOCIATED(subscribedprc))
            IF (subscribedprc%IsActive) THEN
               prc=>subscribedprc%Process
               IF (prc%CfgType%ConfigurationIndex/=atom%Cfg%ConfigurationIndex) THEN
                  WRITE(6,*) "Check>> Atom cfg index does not match with the process Cfg index"
                  WRITE(6,*) "Atom index:",AtomIndex
                  WRITE(6,*) "Atom config:",Atom%Cfg%ConfigurationIndex
                  CALL PrintProcess(prc,6)
                  STOP
               END IF
               prcatom=>prc%AL
               Found=.FALSE.
               DO WHILE (ASSOCIATED(prcatom))
                  Found=prcatom%Atom%Index==AtomIndex
                  IF (Found) EXIT
                  prcatom=>prcatom%NextNeigh
               END DO
               IF (.NOT. Found) THEN
                  WRITE(6,*) "Check>> Atom is not found in the subscriber list of process"
                  WRITE(6,*) "Atom index:",AtomIndex
                  WRITE(6,*) "Atom config:",Atom%Cfg%ConfigurationIndex
                  CALL PrintProcess(prc,6)
                  STOP
               END IF
            ELSE
            END IF
            subscribedprc=>subscribedprc%NextNeigh
         END DO
         atom=>atom%NextNeigh
      END DO
   END SUBROUTINE CheckSubscriptionInfoAllAtoms
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckSubscriberAllProcesses()
      !checks two types of information
      !a) if process contains an ActiveAtom then the ActiveAtom should contain the process
      !b) the configuration of the active atom and the process should match
      !c) the number of subscribers to a process should be correct
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCProcessSubscriptionInfo), POINTER :: subscribedprc
      INTEGER :: NSubscribers,prcindex,AtomIndex
      LOGICAL :: Found
      
      prc=>KMC_PL
      DO WHILE (ASSOCIATED(prc))
         NSubscribers=0
         prcindex=prc%Index
         IF (prcindex<=0) THEN
            WRITE(6,*) "$Err>> CheckSubscriberAllProcesses: Process index is zero"
            STOP
         END IF
         AL1=>prc%AL
         DO WHILE (ASSOCIATED(AL1))
            NSubscribers=NSubscribers+1
            AtomIndex=AL1%Atom%Index !atom to search for
            atom=>KMC_AL
            DO WHILE (ASSOCIATED(atom))
               IF (atom%Index==AtomIndex) EXIT
               atom=>atom%NextNeigh
            END DO
            Found=.FALSE.
            subscribedprc=>atom%PrcSubscriberInfo
            DO WHILE (ASSOCIATED(subscribedprc))
               Found=subscribedprc%Process%Index==prcindex
               IF (Found) EXIT
               subscribedprc=>subscribedprc%NextNeigh
            END DO
            IF (.NOT. Found) THEN
               WRITE(6,*) "Check>> Process subscriber list has incorrect atom"
               WRITE(6,*) "Atom index:",AtomIndex
               WRITE(6,*) "Atom config:",Atom%Cfg%ConfigurationIndex
               CALL PrintProcess(prc,6)
               STOP
            END IF
            AL1=>AL1%NextNeigh
         END DO
         IF (NSubscribers/=prc%NumberSubscribers) THEN
            WRITE(6,*) "Check>> Process has incorrect number of subscribers"
            CALL PrintProcess(prc,6)
            STOP
         END IF
         prc=>prc%NextNeigh
      END DO
   END SUBROUTINE CheckSubscriberAllProcesses
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckActiveProcessList()
      !checks if the process in the ActiveProcessList is indeed active by finding if the number of subsrcibers of that process is greater than zero
      IMPLICIT NONE
      TYPE(KMCProcessList), POINTER :: prc
      
      prc=>ActiveProcessList
      DO WHILE (ASSOCIATED(prc))
         IF (prc%Process%NumberSubscribers==0) THEN
            WRITE(6,*) "Process should not have been present in the ActiveProcessList"
            CALL PrintProcess(prc%Process,6)
            STOP
         END IF
         prc=>prc%NextNeigh
      END DO
   END SUBROUTINE CheckActiveProcessList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckProcessAPLPointer()
      !for each process it is checked that when the number of subscribers is greater than zero the process should be present in the ActiveProcessList.
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessList), POINTER :: prclist
      LOGICAL :: Found
      INTEGER :: prcindex
      
      prc=>KMC_PL
      DO WHILE (ASSOCIATED(prc))
         IF (prc%NumberSubscribers>0) THEN
            prcindex=prc%Index
            prclist=>ActiveProcessList
            DO WHILE (ASSOCIATED(prclist))
               Found=prclist%Process%Index==prcindex
               IF (Found) EXIT
               prclist=>prclist%NextNeigh
            END DO
            IF (.NOT. Found) THEN
               WRITE(6,*) "Process should have been present in the ActiveProcessList"
               CALL PrintProcess(prc,6)
               STOP
            END IF
         END IF
         prc=>prc%NextNeigh
      END DO
   END SUBROUTINE CheckProcessAPLPointer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CheckOverlap()
   !assumes that linked list contains the correct location of the atoms and checks for overlap
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom1
      TYPE(KMCAtomList), POINTER :: AL1
      REAL :: Coord(3),tol,v(3)
      INTEGER :: AtomIndex,ix,iy,iz,CellIndex(3),CellIndex1(3)
      
      tol=RadialTolerance
      atom1=>KMC_AL
      DO WHILE (ASSOCIATED(atom1))
         Coord=atom1%Coord
         CellIndex=GetCellLocation(Coord)
         AtomIndex=atom1%Index
         DO ix=-1,1
            DO iy=-1,1
               DO iz=-1,1
                  CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
                  AL1=>Cell(CellIndex1(1),CellIndex1(2),CellIndex1(3))%AL
                  DO WHILE (ASSOCIATED(AL1))
                     IF (AL1%Atom%Index/=AtomIndex) THEN
                        v=GetPBCSpacing(AL1%Atom%Coord,Coord)
                        IF (ALL(ABS(v)<tol)) THEN
                           WRITE(6,*) "Overlap occured:"
                           WRITE(6,*) "AtomIndex and Coord:",AtomIndex,Coord
                           WRITE(6,*) "AtomIndex and Coord:",AL1%Atom%Index,AL1%Atom%Coord
                           STOP
                        END IF
                     END IF
                     AL1=>AL1%NextNeigh
                  END DO
               END DO
            END DO
         END DO
         atom1=>atom1%NextNeigh
      END DO
   END SUBROUTINE CheckOverlap
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE TestCompareAtomEnvToCfg()
   !user can create "fake" configurations and check if the code can compare the two configurations properly
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: Atom,Atom1,Atom2
      TYPE(Configuration), POINTER :: config
      INTEGER :: NAtoms
      
      !read atom1 environment (relative coordinated have to be provided by the user)
      OPEN(UNIT=405,FILE="KMC-TestCompareAtomEnvToCfg-1.xyz")
      ALLOCATE(Atom1); NAllottedKMCAtom=NAllottedKMCAtom+1
      Atom=>Atom1
      ALLOCATE(Atom%Env); NAllottedEnvironment=NAllottedEnvironment+1
      ALLOCATE(Atom%Env%ShortRangeConfig); NAllottedHistogram=NAllottedHistogram+1
      CALL TestCompareAtomEnvToCfgDriver() !read the neighbor coordinates
      CLOSE(405)

      !read atom2 environment (relative coordinated have to be provided by the user)
      OPEN(UNIT=405,FILE="KMC-TestCompareAtomEnvToCfg-2.xyz")
      ALLOCATE(Atom2); NAllottedKMCAtom=NAllottedKMCAtom+1
      Atom=>Atom2
      ALLOCATE(Atom%Env); NAllottedEnvironment=NAllottedEnvironment+1
      ALLOCATE(Atom%Env%ShortRangeConfig); NAllottedHistogram=NAllottedHistogram+1
      CALL TestCompareAtomEnvToCfgDriver() !read the neighbor coordinates
      CLOSE(405)

      !copy atom2 environment to configuration
      ALLOCATE(config); NAllottedConfiguration=NAllottedConfiguration+1
      CALL CopyAtom2Cfg(Atom2,config)
      
      !Do the comparison
      WRITE(6,*) CompareAtomEnvToCfg(Atom1,config)
      !rather than deallocating the arrays we simply stop the code
      STOP
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE TestCompareAtomEnvToCfgDriver()
      !reads an xyz file to get the configuration

         IMPLICIT NONE
         INTEGER :: NAtoms,iatom,BinNumber,BinRegion
         TYPE(Histogram), POINTER :: hist1
         TYPE(Bin), POINTER :: bin1
         TYPE(BinAtom), POINTER :: bincontent1
         REAL :: coord(3),BondOrder,rad,rad2
         INTEGER :: species
         CHARACTER :: CrapChar
         CHARACTER(len=2) :: AtomSymbol

         READ(405,*) NAtoms
         READ(405,*) CrapChar
         BondOrder=0.
         READ(405,*) AtomSymbol,Atom%Coord  !this is info for the main atom
         AtomSymbol="Ag" !imposed
         IF (AtomSymbol=="Ag") Atom%Species=1
         hist1=>Atom%Env%ShortRangeConfig

         DO iatom=2,NAtoms
            READ(405,*) AtomSymbol,coord  !relative coordinate
            SELECT CASE (AtomSymbol(1:2))
            CASE ("Ag"); species=1
            CASE ("Au"); species=2
            END SELECT
            rad2=coord(1)*coord(1)+coord(2)*coord(2)+coord(3)*coord(3)
            IF (rad2<Rad2Correction) THEN
               rad=SQRT(rad2)
               IF (rad<RadialTolerance) THEN
                  WRITE(6,*) 'Ini>>Radius is smaller than tolerance:',Rad
                  WRITE(6,*) 'Center atom coordinate:',Atom%Coord
                  WRITE(6,*) 'Neighbor atom coord:',coord
                  STOP
               END IF
               
               BondOrder=GetBondOrder(rad,Atom1)
               Atom%BondOrder=Atom%BondOrder+BondOrder
               hist1%Population(species)=hist1%Population(species)+1
               
               CALL GetBinStatus(rad,BinRegion,BinNumber) !find where Atom1 would be present
               IF (ASSOCIATED(hist1%BinInfo)) THEN
                  bin1=>hist1%BinInfo
                  CALL SearchBinAfter(BinNumber,bin1) !finds the correct bin in best case, or else the neighbor
                  IF (BinNumber>bin1%BinNumber) THEN
                     CALL AddBinAfter(BinNumber,bin1)
                  ELSEIF (BinNumber<bin1%BinNumber) THEN
                     CALL AddBinBefore(BinNumber,bin1,Atom)
                  END IF
               ELSE
                  ALLOCATE(hist1%BinInfo); NAllottedBin=NAllottedBin+1
                  bin1=>hist1%BinInfo
                  bin1%BinNumber=BinNumber
               END IF
               ALLOCATE(bincontent1); NAllottedBinAtom=NAllottedBinAtom+1
               !bincontent1%Atom=>Atom1 -- Atom pointer is left empty
               bincontent1%Species=species
               bincontent1%RelCoord=coord
               bincontent1%distance=rad
               CALL RefreshBinEdge(bin1,BinRegion,1,INT(Species)) !works on InBin,BinEdge
               CALL AddToBin(bincontent1,bin1)
            ELSE
               WRITE(6,*) "TestCfg:>> Coordinate provided for atom 1 is outside cut off"
               STOP
            END IF
         END DO
         IF (Atom%BondOrder>BOCutOffForCfg) THEN !such bond orders will not be considered by CompareAtomEnvToCfgList
            WRITE(6,*) "Atom bond order:", atom%BondOrder
            WRITE(6,*) "Max  bond order used in code (BOCutoffForCfg):",BOCutoffForCfg
            WRITE(6,*) "Such bond orders are never going to be considered for comparison in code"
         END IF
      END SUBROUTINE TestCompareAtomEnvToCfgDriver
      !xxxxxxxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE TestCompareAtomEnvToCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CompareAtomEnvToCfg(atom,config)
   !compares an atom environment to an existing configuration, and return TRUE if they match
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Environment), POINTER :: env1,env2
      TYPE(KMCAtom), POINTER :: atom
      TYPE(Configuration), POINTER :: config
      LOGICAL :: CompareAtomEnvToCfg
      
      !WRITE(6,*) "Comparing the atom environment to configuration"
      CompareAtomEnvToCfg= atom%Species==config%Species
      !WRITE(6,*) "AtomSpecies:",atom%Species
      !WRITE(6,*) "configSpecies:",config%Species
      IF (.NOT. CompareAtomEnvToCfg) RETURN
      
      env1=>atom%Env
      env2=>config%Env
      CompareAtomEnvToCfg=CompareCfg2Cfg(env1,env2)
   END FUNCTION CompareAtomEnvToCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CompareCfg2Cfg(env1,env2)
      IMPLICIT NONE
      TYPE(Environment), POINTER :: env1,env2
      TYPE(Histogram), POINTER :: hist1,hist2
      TYPE(Bin), POINTER :: bin1,bin2
      TYPE(BinAtom), POINTER :: binatom1,binatom2
      LOGICAL :: CompareCfg2Cfg,IsMatch
      INTEGER :: BinNor1=0,BinNor2=0,Transfer=0,SumPopulation,SpeciesSelected,MinCount
     
      INTEGER, DIMENSION(NSpeciesType,MaxRadialBins) :: Population1,InBin1,BinEdgeP1,BinEdgeN1
      INTEGER, DIMENSION(NSpeciesType,MaxRadialBins) :: Population2,InBin2,BinEdgeP2,BinEdgeN2
      INTEGER, DIMENSION(NSpeciesType,MaxRadialBins) :: Difference
      INTEGER :: PopTolInUse(NSpeciesType),LastBinMaxMismatchCt(NSpeciesType)
      INTEGER :: LastBinMismatchCt(NSpeciesType)
      REAL :: Coord1Selected(3),Coord1Passed(3),Rad1Selected,distance2
      REAL :: TolDistance2,relcoord(3),TargetCoord(3),RadialTolInUse
      INTEGER :: i,ispec
      
      !Obtain the histograms and do preliminary comparison
      hist1=>env1%ShortRangeConfig
      hist2=>env2%ShortRangeConfig
      
      !Check if populations match -- this has been replaced by the overall balance below
      !PopTolInUse=CEILING(REAL(hist1%Population)/10.) !max allowable deviation
      !PopTolInUse=0
      !CompareCfg2Cfg = ALL(ABS(hist1%Population-hist2%Population)<=PopTolInUse)
      !IF (.NOT. CompareCfg2Cfg) RETURN
      
      !Check if radial histograms match
      Population1=0
      Population2=0
      InBin1=0
      InBin2=0
      BinEdgeN1=0
      BinEdgeN2=0
      BinEdgeP1=0
      BinEdgeP2=0
      Difference=0
      
      !Goal: To ensure that all atoms in config within cutoff EnvCutOffRadius are definitely present in atom environment
      !Note: config atoms are only within distance EnvCutOffRadius from the center
      !#1 corresponds to the atom information
      bin1=>hist1%BinInfo
      !WRITE(*,*) "Printing atom bin info ..."
      !IF(NSpeciesType==1) WRITE(6,*) "       BinNo       InBin     BinEdgeP    BinEdgeN"
      !IF(NSpeciesType==2) &
      !   WRITE(6,*) "        BinNor          InBin                     BinEdge                BinEdgeN "
      DO WHILE (ASSOCIATED(bin1))
         BinNor1=bin1%BinNumber
         InBin1(:,BinNor1)=bin1%InBin
         BinEdgeP1(:,BinNor1)=bin1%BinEdgeP
         BinEdgeN1(:,BinNor1)=bin1%BinEdgeN
         !WRITE(6,*) BinNor1,InBin1(:,BinNor1),BinEdgeP1(:,BinNor1),BinEdgeN1(:,BinNor1)
         bin1=>bin1%NextNeigh
      END DO
      Population1=InBin1+BinEdgeP1+BinEdgeN1
      
      !#2 corresponds to the cfg information
      bin2=>hist2%BinInfo
      !WRITE(*,*) "Printing cfg bin info ..."
      !IF(NSpeciesType==1) WRITE(6,*) "       BinNo       InBin     BinEdgeP    BinEdgeN"
      !IF(NSpeciesType==2) &
      !  WRITE(6,*) "        BinNor          InBin                    BinEdge                BinEdgeN "
      DO WHILE (ASSOCIATED(bin2))
         BinNor2=bin2%BinNumber
         InBin2(:,BinNor2)=bin2%InBin
         BinEdgeP2(:,BinNor2)=bin2%BinEdgeP
         BinEdgeN2(:,BinNor2)=bin2%BinEdgeN
         !WRITE(6,*) BinNor2,InBin2(:,BinNor2),BinEdgeP2(:,BinNor2),BinEdgeN2(:,BinNor2)
         bin2=>bin2%NextNeigh
      END DO
      Population2=InBin2+BinEdgeP2+BinEdgeN2
      
      !OVERALL BALANCE
      !making sure that the right number of atoms is present in the cutoff for the cfg
      !Derivation:
      !Denote Population1 is denoted as p1 and Population2 is p2
      !p1 and p2 the number of atoms in the cutoff given by SQRT(Rad2Correction)=EnvCutOffRadius+RadialTolerance
      !Hence p1 is the sum of two numbers x1 and y: x1 is the number of atoms which definitely lie 
      !outside the cutoff EnvCutOffRadius+fuzzy term and y is the number of atoms that could lie inside the 
      !cutoff EnvCutOffRadius+fuzzy term -- here the fuzzy term is the atoms which is outside EnvCutOffRadius in atom env
      !similarly p2=x2+y
      !Note: We want to find the case where y1=y2=y
      !Hence p1=x1+y and p2=x2+y therefore p1-p2=x1-x2
      !Now lets look at the value of BinEdgeP1(Outermost bin,:) -- this contains the outermost atoms. We
      !shall represent BinEdgeP1(Outermost bin,:) as z1. it is z1 and p1 that is know to us, but not x1 and y
      !Nonetheless, z1=x1+w1 here w1 is the number of atoms that are definitely in .
      !Similarly, z2=x2-w2 where w2 can be different from w1. 
      !Hence z1-z2=x1-x2+w1-w2 but 0<=x1<=z1 and 0<=x2<=z2. hence -z2<=x1-x2<z1
      !there -z2<=Population1-Population2<z1
      Difference=Population2-Population1 !note we are using Pop2-Pop1
      !WRITE(6,*) "Difference:"
      !DO i=1,MaxRadialBins
      !   WRITE(6,*) i,Difference(:,i)
      !END DO
      
      !IF(atom%INdex==427) WRITE(6,*) "BINCOM",SUM(Difference,2),NumberBinsInEnv,BinEdgeP2(:,NumberBinsInEnv),&
      !   BinEdgeP1(:,NumberBinsInEnv)
      !CompareCfg2Cfg= ALL(SUM(Difference,2)<=BinEdgeP2(:,BinNor2)) .AND. &
      !   ALL(SUM(Difference,2)>= -BinEdgeP1(:,BinNor1))
      CompareCfg2Cfg= ALL(SUM(Difference,2)<=BinEdgeP2(:,NumberBinsInEnv)) .AND. &
         ALL(SUM(Difference,2)>= -BinEdgeP1(:,NumberBinsInEnv))
      IF (.NOT. CompareCfg2Cfg) RETURN
      
      !WRITE(6,*) "Bin based balance ..."
      !WRITE(6,*) "Number bins in env:",NumberBinsInEnv
      !WRITE(6,*) " "
      !BIN BASED BALANCE
      DO ispec=1,NSpeciesType     
         DO i=1,NumberBinsInEnv-1 !outer bin will be checked separately
            IF (Difference(ispec,i)==0) CYCLE
            !each time there is a difference I will try to bring some atoms/ push some atoms from/to bin i+1
            !WRITE(6,*) "bin:",i,ispec,Difference(ispec,i)
            IF (Difference(ispec,i)<0) THEN
               Transfer=MIN(-Difference(ispec,i),BinEdgeN2(ispec,i+1)+BinEdgeP1(ispec,i)) !we can move some atoms from i+1 bin of config 
                  !to ith bin and similarly move some atoms of ith bin of atom to i+1 bin -- this will be used to lower the difference
               Difference(ispec,i)=Difference(ispec,i)+Transfer
               Difference(ispec,i+1)=Difference(ispec,i+1)-Transfer
            ELSE
               Transfer=MIN(Difference(ispec,i),BinEdgeN1(ispec,i+1)+BinEdgeP2(ispec,i))
               Difference(ispec,i)=Difference(ispec,i)-Transfer
               Difference(ispec,i+1)=Difference(ispec,i+1)+Transfer
               !IF(Transfer/=0) WRITE(*,*) "transferred...",i+1,Difference(ispec,i+1)
            END IF
            IF (Difference(ispec,i)/=0) THEN
               CompareCfg2Cfg=.FALSE.
               RETURN
            END IF
         END DO
      END DO
      !WRITE(6,*) "After bin based balance, CompareCfg2Cfg",CompareCfg2Cfg
      
      !last radial bin
      !Check for last bin is not required -- we know that overall the populations are in agreement
      !DO ispec=1,NSpeciesType
      !   LastBinMaxMismatchCt=0 !# of mismatches allowed at the last bin edge
      !   IF (Difference(NumberBinsInEnv,ispec)==0) CYCLE
      !   IF (Difference(NumberBinsInEnv,ispec)>0) THEN
      !      Transfer=MIN(Difference(NumberBinsInEnv,ispec),BinEdgeP1(NumberBinsInEnv,ispec))
      !      Difference(i,ispec)=Difference(i,ispec)-Transfer !removed excess contribution from BinEdgeP1
      !      LastBinMaxMismatchCt(ispec)=(BinEdgeP1(NumberBinsInEnv,ispec)-Transfer)
      !   ELSE
      !     Transfer=MIN(-Difference(NumberBinsInEnv,ispec),BinEdgeP2(NumberBinsInEnv,ispec))
      !     Difference(i,ispec)=Difference(i,ispec)+Transfer !remove excess contribution from BinEdgeP2
      !      LastBinMaxMismatchCt(ispec)=(BinEdgeP2(NumberBinsInEnv,ispec)-Transfer)
      !   END IF
      !   WHERE (LastBinMaxMismatchCt==0.) LastBinMaxMismatchCt=2.
      !   IF (Difference(NumberBinsInEnv,ispec)/=0) THEN
      !      CompareCfg2Cfg=.FALSE.
      !      RETURN
      !   END IF
      !END DO

      !If the comparison has made it so far that implies that the number of atoms in the
      !bins match nicely. Next we shall do a one-to-one comparison to ensure that the 
      !correct arrangement is present -- note, some tolerance is left in case of a mismatch
      !WRITE(6,*) "Rotational matching ..."
      IsMatch=.FALSE. !incase of only cycling
      RadialTolInUse=RadialTolerance !*REAL(ptrhis1%BinNor)/REAL(BinNor)
      TolDistance2=RadialTolInUse**2
      LastBinMaxMismatchCt=2
      !hist1=>atom%Env%ShortRangeConfig
      !hist2=>config%Env%ShortRangeConfig

      bin1=>hist1%BinInfo !atom environment
      !because we know that the total atoms in each bin match up with atoms in bin of config
      !hence we can simply compare atoms in atom environment to those of config -- most importantly
      !the reverse procedure is not required -- 
      DO WHILE (ASSOCIATED(bin1))
         LastBinMismatchCt=0
         binatom1=>bin1%AL
         DO WHILE (ASSOCIATED(binatom1))
            IsMatch=.FALSE.
            bin2=>hist2%BinInfo
            DO WHILE (ASSOCIATED(bin2))
               IF (ABS(bin1%BinNumber-bin2%BinNumber)<=1) THEN
                  binatom2=>bin2%AL
                  DO WHILE (ASSOCIATED(binatom2))
                     IF (ABS(binatom2%distance-binatom1%distance)<RadialTolInUse  &
                        .AND. binatom2%Species==binatom1%Species) THEN
                        relcoord=binatom1%RelCoord-binatom2%RelCoord
                        distance2=relcoord(1)*relcoord(1)+relcoord(2)*relcoord(2)+relcoord(3)*relcoord(3)
                        IsMatch= (distance2<=TolDistance2)
                     END IF
                     IF (IsMatch) EXIT
                     binatom2=>binatom2%NextNeigh
                  END DO
               END IF
               IF (IsMatch) EXIT
               bin2=>bin2%NextNeigh
            END DO
            !WRITE(6,*) bin1%BinNumber,binatom1%RelCoord,IsMatch
            IF (.NOT. IsMatch) THEN !could not find the atom
               !lets account for error at last bin edges if this atom is at last bin edge
               IsMatch= binatom1%distance>EnvCutOffRadius-RadialTolerance !simpler criterion than the one below
               !IF (binatom1%distance> (EnvCutOffRadius-RadialTolerance)) THEN !last bin edge
               !   LastBinMismatchCt(binatom1%Species)=LastBinMismatchCt(binatom1%Species)+1
               !   IsMatch= (ALL(LastBinMismatchCt<=LastBinMaxMismatchCt))
               !END IF
            END IF
            IF (.NOT. IsMatch) EXIT !match has failed; exit
            binatom1=>binatom1%NextNeigh
         END DO
         IF (.NOT. IsMatch) EXIT !match has failed; exit
         bin1=>bin1%NextNeigh
      END DO
      
      CompareCfg2Cfg=IsMatch
      IF (.NOT. CompareCfg2Cfg) RETURN
      
      !Now compare the cfg environment to the atom environment
      IsMatch=.FALSE. !incase of only cycling
      bin1=>hist2%BinInfo !cfg environment
      DO WHILE (ASSOCIATED(bin1))
         LastBinMismatchCt=0
         binatom1=>bin1%AL
         DO WHILE (ASSOCIATED(binatom1))
            IsMatch=.FALSE.
            bin2=>hist1%BinInfo
            DO WHILE (ASSOCIATED(bin2))
               IF (ABS(bin1%BinNumber-bin2%BinNumber)<=1) THEN
                  binatom2=>bin2%AL
                  DO WHILE (ASSOCIATED(binatom2))
                     IF (ABS(binatom2%distance-binatom1%distance)<RadialTolInUse  &
                        .AND. binatom2%Species==binatom1%Species) THEN
                        relcoord=binatom1%RelCoord-binatom2%RelCoord
                        distance2=relcoord(1)*relcoord(1)+relcoord(2)*relcoord(2)+relcoord(3)*relcoord(3)
                        IsMatch= (distance2<=TolDistance2)
                     END IF
                     IF (IsMatch) EXIT
                     binatom2=>binatom2%NextNeigh
                  END DO
               END IF
               IF (IsMatch) EXIT
               bin2=>bin2%NextNeigh
            END DO
            !WRITE(6,*) bin1%BinNumber,binatom1%RelCoord,IsMatch
            IF (.NOT. IsMatch) THEN !could not find the atom
               !lets account for error at last bin edges if this atom is at last bin edge
               IsMatch= binatom1%distance>EnvCutOffRadius-RadialTolerance !simpler criterion than the one below
               !IF (binatom1%distance> (EnvCutOffRadius-RadialTolerance)) THEN !last bin edge
               !   LastBinMismatchCt(binatom1%Species)=LastBinMismatchCt(binatom1%Species)+1
               !   IsMatch= (ALL(LastBinMismatchCt<=LastBinMaxMismatchCt))
               !END IF
            END IF
            IF (.NOT. IsMatch) EXIT !match has failed; exit
            binatom1=>binatom1%NextNeigh
         END DO
         IF (.NOT. IsMatch) EXIT !match has failed; exit
         bin1=>bin1%NextNeigh
      END DO
      
      CompareCfg2Cfg=IsMatch
      RETURN
   END FUNCTION CompareCfg2Cfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CompareAtomEnvToCfgList(atom)
   !finds which configuration matches with the atom environment
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(Configuration), POINTER :: CompareAtomEnvToCfgList,cfg1
      TYPE(Histogram), POINTER :: hist1
      LOGICAL :: IsConfigExist,EndOfList
      REAL :: ratiobondorder
      INTEGER :: RedundantCfgIndex
      
      IsConfigExist=.FALSE.
      EndOfList=.FALSE.
      IF (ASSOCIATED(atom%Cfg)) THEN
         IF (atom%Cfg%ConfigurationIndex/=RedundantCfg%ConfigurationIndex) &
            IsConfigExist=CompareAtomEnvToCfg(atom,atom%Cfg) !compare with current know atom Cfg
         IF (IsConfigExist) THEN
            CompareAtomEnvToCfgList=>atom%Cfg
            RETURN
         END IF
      END IF
      
      cfg1=>RecordedCfg
      hist1=>atom%Env%ShortRangeConfig
      RedundantCfgIndex=RedundantCfg%ConfigurationIndex

      DO WHILE (ASSOCIATED(cfg1))
         IF (cfg1%ConfigurationIndex==RedundantCfgIndex) EXIT !last configuration
         ratiobondorder=cfg1%BondOrder/atom%BondOrder
         IF (ratiobondorder>1.4) EXIT !Add before current config -- cfgs are stored in increasing bond order
         !IF (ASSOCIATED(atom%cfg)) THEN
         !   IF (cfg1%ConfigurationIndex/=atom%cfg%ConfigurationIndex) THEN
         !      IsConfigExist=CompareAtomEnvToCfg(atom,cfg1)
         !   END IF
         !ELSE
            IsConfigExist=CompareAtomEnvToCfg(atom,cfg1)
         !END IF
         IF (IsConfigExist) EXIT
         IF (ASSOCIATED(cfg1%NextNeigh)) THEN
            cfg1=>cfg1%NextNeigh
         ELSE
            EndOfList=.TRUE.
            EXIT
         END IF
      END DO

      IF (.NOT. IsConfigExist) CALL AddAtomToRecordedCfg(atom,cfg1,EndOfList) !needs to add before cfg1
      CompareAtomEnvToCfgList=>cfg1
      
      !IF (cfg1%ConfigurationIndex>3 .AND. cfg1%ConfigurationIndex<8 .AND. cfg1%NAtomsWithCfg>1) THEN
      !   WRITE(*,*) "PRINTING EXCEPTIONAL CASE"
      !   CALL PrintCfg(cfg1)
      !   CALL PrintEnv(cfg1%AtomsWithCfg%Atom)
      !   CALL PrintEnv(cfg1%AtomsWithCfg%NextNeigh%Atom)
      !   STOP
      !END IF
   END FUNCTION CompareAtomEnvToCfgList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CompareProcess2PL(prc,MatchingProcess,MinDisplacement)
   !FALSE if prc is a new process
   !Last checked: Jan 12, 2011
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc,prc1,MatchingProcess
      LOGICAL :: CompareProcess2PL, Found
      REAL :: tol,MinDisplacement !MinDisplacement is the minimum displacement allowed for the process
      
      Found=.FALSE.
      tol=RadialTolerance*2. !propagation of error of second order
      prc1=>KMC_PL
      NULLIFY(MatchingProcess) !if found stores the process that matched
      DO WHILE (ASSOCIATED(prc1))
         Found=CompareProcess2Process(prc,prc1,tol,MinDisplacement)
         IF (Found) THEN
            MatchingProcess=>prc1
            EXIT
         END IF
         prc1=>prc1%NextNeigh
      END DO
      
      CompareProcess2PL=Found
   END FUNCTION CompareProcess2PL
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE TestCompareProcess2Process()
   !user can create "fake" processes and check if the code can compare the two processes properly
      IMPLICIT NONE
   END SUBROUTINE TestCompareProcess2Process
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION CompareProcess2Process(prc1,prc2,tol,MinDisplacement)
   !TRUE if prc1 and prc2 are identical within tolerance tol
   !Last checked: Jan 12, 2011
      IMPLICIT NONE
      LOGICAL :: CompareProcess2ProcessA,CompareProcess2ProcessB,CompareProcess2Process,NotFound
      TYPE(KMCProcess), POINTER :: prc1,prc2
      TYPE(KMCProcessAtom), POINTER :: AL1,AL2
      REAL :: tol,v1,v2,v3,MinDisplacement,coord(3) !,MaxDisplacement1,MaxDisplacement2
      INTEGER :: cfgindex1
      
      !WRITE(6,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxx"
      CompareProcess2Process=.FALSE.
      IF (prc1%CfgType%ConfigurationIndex/=prc2%CfgType%ConfigurationIndex) RETURN
      !IF (prc1%NumberProcessAtoms/=prc2%NumberProcessAtoms) RETURN !this is a very strict requirement
      
      IF (prc1%Frequency==0. .OR. prc2%Frequency==0.) THEN
         WRITE(6,*) "While comparing processes, at least one process has zero prefactor"
         STOP
      END IF
      IF (MAX(prc1%Frequency,prc2%Frequency)/MIN(prc1%Frequency,prc2%Frequency)>5.) RETURN
      
      IF (prc1%FwdBarrier==0. .OR. prc2%FwdBarrier==0.) THEN
         WRITE(6,*) "While comparing processes, at least one process has zero activation barrier"
         STOP
      END IF
      !WRITE(6,*) "Barriers:",prc1%FwdBarrier,prc2%FwdBarrier
      IF (MAX(prc1%FwdBarrier,prc2%FwdBarrier)-MIN(prc1%FwdBarrier,prc2%FwdBarrier)>0.1) RETURN
      
      !for all process atoms in prc1 find their maximum displacement
      !MaxDisplacement1=prc1%MaxPrcAtomDisplacement
      !MaxDisplacement2=prc2%MaxPrcAtomDisplacement
      !IF (MaxDisplacement1==0. .OR. MaxDisplacement2==0.) THEN
      !   WRITE(6,*) "Maximum displacement of process atoms found to be zero in CompareProcess2Process"
      !   STOP
      !END IF
      !WRITE(6,*) "MaxDisplacements:",MaxDisplacement1,MaxDisplacement2
      
      !atom-by-atom comparison
      AL1=>prc1%ReactantAL
      DO WHILE (ASSOCIATED(AL1))
         cfgindex1=AL1%CfgIndex
         AL2=>prc2%ReactantAL
         NotFound=.TRUE.
         DO WHILE (ASSOCIATED(AL2))
            IF (cfgindex1==AL2%CfgIndex) THEN
               v1=MAXVAL(ABS(AL1%RelInitialPos-AL2%RelInitialPos))
               v2=MAXVAL(ABS(AL1%RelFinalPos-AL2%RelFinalPos))
               v3=MAXVAL(ABS(AL1%RelSaddlePos-AL2%RelSaddlePos))
               NotFound= (MAX(v1,v2,v3)>tol)
            END IF
            IF (.NOT. NotFound) EXIT
            AL2=>AL2%NextNeigh
         END DO
         !WRITE(6,*) "Preliminary search for atom resulted in ",NotFound
         !IF (NotFound) THEN !we need to understand why we did not find the atom
         !   coord=AL1%RelInitialPos-AL1%RelFinalPos
         !   v1=coord(1)*coord(1)+coord(2)*coord(2)+coord(3)*coord(3)
         !   coord=AL1%RelInitialPos-AL1%RelSaddlePos
         !   v1=coord(1)*coord(1)+coord(2)*coord(2)+coord(3)*coord(3)
         !   v1=MAX(SQRT(v1)-MinDisplacement,0.001) !min value of v1 is 0.001
         !   WRITE(6,*) "Displacement of process atom:",v1
         !   WRITE(6,*) "Ratio of displacements:",(MaxDisplacement1-MinDisplacement)/v1
         !   NotFound= (MaxDisplacement1-MinDisplacement)/v1<50. !the amount by which our process atom got displaced is significant
         !END IF
         IF (NotFound) EXIT
         AL1=>AL1%NextNeigh
      END DO
      CompareProcess2ProcessA=.NOT. NotFound
      
      !atom-by-atom comparison by interchanging the processes (in case one of the process has fewer process atoms
      ! although they refer to the same process)
      CompareProcess2ProcessB=CompareProcess2ProcessA
      IF (prc1%NumberProcessAtoms/=prc2%NumberProcessAtoms) THEN
         AL1=>prc2%ReactantAL
         DO WHILE (ASSOCIATED(AL1))
            cfgindex1=AL1%CfgIndex
            AL2=>prc1%ReactantAL
            NotFound=.TRUE.
            DO WHILE (ASSOCIATED(AL2))
               IF (cfgindex1==AL2%CfgIndex) THEN
                  v1=MAXVAL(ABS(AL1%RelInitialPos-AL2%RelInitialPos))
                  v2=MAXVAL(ABS(AL1%RelFinalPos-AL2%RelFinalPos))
                  v3=MAXVAL(ABS(AL1%RelSaddlePos-AL2%RelSaddlePos))
                  NotFound= (MAX(v1,v2,v3)>tol)
               END IF
               IF (.NOT. NotFound) EXIT
               AL2=>AL2%NextNeigh
            END DO
            !IF (NotFound) THEN !we need to understand why we did not find the atom
            !   coord=AL1%RelInitialPos-AL1%RelFinalPos
            !   v1=coord(1)*coord(1)+coord(2)*coord(2)+coord(3)*coord(3)
            !   v1=MAX(SQRT(v1)-MinDisplacement,0.001) !min value of v1 is 0.001
            !   NotFound= MaxDisplacement2/v1<50. !the amount by which our process atom got displaced is significant
            !END IF
            IF (NotFound) EXIT
            AL1=>AL1%NextNeigh
         END DO
         CompareProcess2ProcessB=.NOT. NotFound
!WRITE(6,*) "CompareProcess2ProcessB:",CompareProcess2ProcessB
      END IF
!  CHART OF MIN DISPLACEMENTS AND THE RATIOS
!  dx(actual displacement)  MinDisplacement (for process)  MaxDisplacement (for process) Ratio (MaxDisplacement-MinDisplacement)/(dx-MinDisplacement)
!       0.101                   0.1                             2.1                      2000
!       0.15                    0.1                             2.1                      40
!       0.35                    0.1                             2.1                      8
!       0.35                    0.1                             10.1                     40
!       0.35                    0.1                             100.1                    400
      CompareProcess2Process=CompareProcess2ProcessA .OR. CompareProcess2ProcessB
   END FUNCTION CompareProcess2Process
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CopyCfg2Cfg(cfg1,cfg2)
     TYPE(Configuration), POINTER :: cfg1,cfg2
     TYPE(Bin), POINTER :: bin1,bin2
     TYPE(BinAtom), POINTER :: binatom1,binatom2
     
     ALLOCATE(cfg2); NAllottedConfiguration=NAllottedConfiguration+1
     cfg2%Species=cfg1%Species
     ALLOCATE(cfg2%Env); NAllottedEnvironment=NAllottedEnvironment+1
     ALLOCATE(cfg2%Env%ShortRangeConfig); NAllottedHistogram=NAllottedHistogram+1
     cfg2%Env%ShortRangeConfig%Population=cfg1%Env%ShortRangeConfig%Population
     cfg2%BondOrder=cfg1%BondOrder
     IF (ASSOCIATED(cfg1%Env%ShortRangeConfig%BinInfo)) &
          ALLOCATE(cfg2%Env%ShortRangeConfig%BinInfo); NAllottedBin=NAllottedBin+1
     bin1=>cfg1%Env%ShortRangeConfig%BinInfo
     bin2=>cfg2%Env%ShortRangeConfig%BinInfo
     DO WHILE (ASSOCIATED(bin1))
        bin2%BinNumber=bin1%BinNumber
        bin2%InBin=bin1%InBin
        bin2%BinEdgeP=bin1%BinEdgeP
        bin2%BinEdgeN=bin1%BinEdgeN
        binatom1=>bin1%AL
        IF (ASSOCIATED(binatom1)) THEN
           ALLOCATE(bin2%AL); NAllottedBinAtom=NAllottedBinAtom+1
           binatom2=>bin2%AL
        END IF
        DO WHILE (ASSOCIATED(binatom1))
           binatom2%distance=binatom1%distance
           binatom2%RelCoord=binatom1%RelCoord
           binatom2%Species=binatom1%Species
           binatom1=>binatom1%NextNeigh
           IF (ASSOCIATED(binatom1)) THEN
              ALLOCATE(binatom2%NextNeigh); NAllottedBinAtom=NAllottedBinAtom+1
              binatom2%NextNeigh%PrevNeigh=>binatom2
              binatom2=>binatom2%NextNeigh
              NULLIFY(binatom2%NextNeigh)
           END IF
        END DO
        bin1=>bin1%NextNeigh
        IF (ASSOCIATED(bin1)) THEN
           ALLOCATE(bin2%NextNeigh); NAllottedBin=NAllottedBin+1
           bin2%NextNeigh%PrevNeigh=>bin2
           bin2=>bin2%NextNeigh
           NULLIFY(bin2%NextNeigh)
        END IF
     END DO
   END SUBROUTINE CopyCfg2Cfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CopyAtom2Cfg(atom,cfg)
   !copies the local environment of atom to the new element cfg
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(Configuration), POINTER :: cfg !has already been allocated
      TYPE(Bin), POINTER :: bin1,bin2
      TYPE(BinAtom), POINTER :: binatom1,binatom2
      
      cfg%Species=atom%Species
      cfg%BondOrder=atom%BondOrder
      ALLOCATE(cfg%Env); NAllottedEnvironment=NAllottedEnvironment+1
      ALLOCATE(cfg%Env%ShortRangeConfig); NAllottedHistogram=NAllottedHistogram+1
      cfg%Env%ShortRangeConfig%Population=atom%Env%ShortRangeConfig%Population
      IF (ASSOCIATED(atom%Env%ShortRangeConfig%BinInfo)) &
           ALLOCATE(cfg%Env%ShortRangeConfig%BinInfo); NAllottedBin=NAllottedBin+1
      bin1=>atom%Env%ShortRangeConfig%BinInfo
      bin2=>cfg%Env%ShortRangeConfig%BinInfo
      DO WHILE (ASSOCIATED(bin1))
         bin2%BinNumber=bin1%BinNumber
         bin2%InBin=bin1%InBin
         bin2%BinEdgeP=bin1%BinEdgeP
         bin2%BinEdgeN=bin1%BinEdgeN
         binatom1=>bin1%AL
         IF (ASSOCIATED(binatom1)) THEN
            ALLOCATE(bin2%AL); NAllottedBinAtom=NAllottedBinAtom+1
            binatom2=>bin2%AL
         END IF
         DO WHILE (ASSOCIATED(binatom1))
            binatom2%distance=binatom1%distance
            binatom2%RelCoord=binatom1%RelCoord
            binatom2%Species=binatom1%Species
            binatom1=>binatom1%NextNeigh
            IF (ASSOCIATED(binatom1)) THEN
               ALLOCATE(binatom2%NextNeigh); NAllottedBinAtom=NAllottedBinAtom+1
               binatom2%NextNeigh%PrevNeigh=>binatom2
               binatom2=>binatom2%NextNeigh
               NULLIFY(binatom2%NextNeigh)
            END IF
         END DO
         bin1=>bin1%NextNeigh
         IF (ASSOCIATED(bin1)) THEN
            ALLOCATE(bin2%NextNeigh); NAllottedBin=NAllottedBin
            bin2%NextNeigh%PrevNeigh=>bin2
            bin2=>bin2%NextNeigh
            NULLIFY(bin2%NextNeigh)
         END IF
      END DO
   END SUBROUTINE CopyAtom2Cfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteReactantAtomFromProcess(atom,prc)
   !deletes a reactant atom from the process
   !Last checked: Jan 12, 2011
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: atom
      
      IF (ASSOCIATED(atom%PrevNeigh)) THEN
         IF (ASSOCIATED(atom%NextNeigh)) THEN
            atom%PrevNeigh%NextNeigh=>atom%NextNeigh
            atom%NextNeigh%PrevNeigh=>atom%PrevNeigh
         ELSE
            NULLIFY(atom%PrevNeigh%NextNeigh)
         END IF
      ELSE
         IF (ASSOCIATED(atom%NextNeigh)) THEN
            prc%ReactantAL=>atom%NextNeigh
            NULLIFY(atom%NextNeigh%PrevNeigh)
         ELSE
            NULLIFY(prc%ReactantAL)
         END IF
      END IF
      NULLIFY(atom%PrevNeigh)
      NULLIFY(atom%NextNeigh)
      DEALLOCATE(atom); NAllottedKMCProcessAtom=NAllottedKMCProcessAtom-1
   END SUBROUTINE DeleteReactantAtomFromProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteProcess(prc)
   !Deletes a process element prc that is not needed
   !This process element is not associated with any neighbor
   !Last checked: Jan 12, 2011
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessAtom), POINTER :: atom
      
      IF (ASSOCIATED(prc%AL)) THEN
         WRITE(6,*) "Atoms are associated with the process .. cannot delete"
         STOP
      END IF
      
      IF (ASSOCIATED(prc%PrevNeigh) .OR. ASSOCIATED(prc%NextNeigh)) THEN
         WRITE(6,*) "Process to be deleted has a neighbor associated with it"
         STOP
      END IF
      
      atom=>prc%ReactantAL
      DO WHILE (ASSOCIATED(atom))
         CALL DeleteReactantAtomFromProcess(atom,prc)
         atom=>prc%ReactantAL
      END DO
      
      NULLIFY(prc%CfgType)
      DEALLOCATE(prc); NAllottedKMCProcess=NAllottedKMCProcess-1
   END SUBROUTINE DeleteProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DeleteLocalCfg(atom)
   !deletes the local configuration for an atom
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(Bin), POINTER :: bin1,bin2
      
      bin1=>atom%Env%ShortRangeConfig%BinInfo
      
      DO WHILE (ASSOCIATED(bin1))
         bin2=>bin1%NextNeigh
         CALL DeleteLocalCfgBin() !delete this bin
         DEALLOCATE(bin1); NAllottedBin=NAllottedBin-1
         bin1=>bin2
      END DO
      NULLIFY(atom%Env%ShortRangeConfig%BinInfo)
      atom%Env%ShortRangeConfig%Population=0
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE DeleteLocalCfgBin()
      !To be used only by DeleteLocalCfg
         IMPLICIT NONE
         TYPE(BinAtom), POINTER :: binatom1,binatom2
      
         NULLIFY(bin1%PrevNeigh)
         NULLIFY(bin1%NextNeigh)
      
         binatom1=>bin1%AL
         NULLIFY(bin1%AL)
         DO WHILE (ASSOCIATED(binatom1))
            binatom2=>binatom1%NextNeigh
            NULLIFY(binatom1%Atom)
            DEALLOCATE(binatom1); NAllottedBinAtom=NAllottedBinAtom-1
            binatom1=>binatom2
         END DO
      END SUBROUTINE DeleteLocalCfgBin
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE DeleteLocalCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ExtendAffectedCellLayer(CellIndex1,value,range)
   !adds a layer of affected cells around CellIndex1
      IMPLICIT NONE
      INTEGER :: value,range,ix,iy,iz,rr2
      INTEGER, DIMENSION(3) :: CellIndex1,CellIndex2
      
      DO ix=-range,range
         DO iy=-range,range
            DO iz=-range,range
               rr2=ix*ix+iy*iy+iz*iz
               IF (rr2<=range*range+2) THEN
                  CellIndex2=GetNeighCell(CellIndex1,ix,iy,iz)
                  CALL AddCellToAffectedList(CellIndex2,value=value)
               END IF
            END DO
         END DO
      END DO
      !CellIndex2=GetNeighCell(CellIndex1,-1,-1,-1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 0,-1,-1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 1,-1,-1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1,-1, 0,-1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 0, 0,-1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 1, 0,-1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1,-1, 1,-1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 0, 1,-1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 1, 1,-1); CALL AddCellToAffectedList(CellIndex2,value=value)
      
      !CellIndex2=GetNeighCell(CellIndex1,-1,-1, 0); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 0,-1, 0); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 1,-1, 0); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1,-1, 0, 0); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 1, 0, 0); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1,-1, 1, 0); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 0, 1, 0); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 1, 1, 0); CALL AddCellToAffectedList(CellIndex2,value=value)
                  
      !CellIndex2=GetNeighCell(CellIndex1,-1,-1, 1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 0,-1, 1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 1,-1, 1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1,-1, 0, 1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 0, 0, 1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 1, 0, 1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1,-1, 1, 1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 0, 1, 1); CALL AddCellToAffectedList(CellIndex2,value=value)
      !CellIndex2=GetNeighCell(CellIndex1, 1, 1, 1); CALL AddCellToAffectedList(CellIndex2,value=value)
   END SUBROUTINE ExtendAffectedCellLayer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MarkAffectedCells(atom,displacement,ncutoff)
   !will provide information about cells that are affected by the process
   !cfg of atoms in these cells need to be updated
   !Last checked: Jan 12, 2011
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      REAL :: displacement(3),coord(3)
      INTEGER :: CellIndex(3),CellIndex1(3),CellIndex2(3),ix,iy,iz,i,ncutoff
      
      !NAffectedCells=0
      !Initial position
      CellIndex=GetCellLocation(atom%Coord)
      DO ix=-ncutoff,ncutoff
         DO iy=-ncutoff,ncutoff
            DO iz=-ncutoff,ncutoff
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               CALL AddCellToAffectedList(CellIndex1,value=1)
            END DO
         END DO
      END DO
      !Final position
      coord=GetPBCCoord(atom%Coord+displacement)
      CellIndex2=GetCellLocation(coord)
      IF (ANY(CellIndex/=CellIndex2)) THEN
         DO ix=-ncutoff,ncutoff
            DO iy=-ncutoff,ncutoff
               DO iz=-ncutoff,ncutoff
                  CellIndex1=GetNeighCell(CellIndex2,ix,iy,iz)
                  CALL AddCellToAffectedList(CellIndex1,value=1)
               END DO
            END DO
         END DO
      END IF
      
   END SUBROUTINE MarkAffectedCells
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ProcessNotAssigned(atom,prc)
   !checks whether the process prc is already assigned to atom
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessSubscriptionInfo), POINTER :: prclist
      INTEGER :: prcindex
      LOGICAL :: ProcessNotAssigned
      
      prcindex=prc%Index
      prclist=>atom%PrcSubscriberInfo
      DO WHILE (ASSOCIATED(prclist))
         IF (prclist%Process%Index==prcindex) EXIT
         prclist=>prclist%NextNeigh
      END DO
      ProcessNotAssigned=.NOT. ASSOCIATED(prclist)
   END FUNCTION ProcessNotAssigned
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RefreshCell(atom)
   !Refreshes the cell associated with an atom 
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      INTEGER :: CellIndexOld(3),CellIndexNew(3)
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCAtomList), POINTER :: AL1

      IF (ALL(atom%CellIndex==0)) THEN
         CellIndexNew=GetCellLocation(atom%Coord)
         ALLOCATE(AL1); NAllottedKMCAtomList=NAllottedKMCAtomList+1
         AL1%atom=>atom
         CALL AddToCell(AL1,CellIndexNew)
         AL1%Atom=>atom
         atom%CellIndex=CellIndexNew
      ELSEIF (ALL(atom%CellIndex>0)) THEN !atom should already be attached to a cell
         CellIndexNew=GetCellLocation(atom%Coord)
         CellIndexOld=atom%CellIndex
         IF (ALL(CellIndexOld==CellIndexNew)) RETURN !no changes
         AL1=>SearchCell(atom,CellIndexOld)
         IF (.NOT. ASSOCIATED(AL1)) THEN
            WRITE(6,*) "Db>>SearchCellList..Site not found in subcell list in FindCellLocation"
            STOP
         END IF
         CALL RemoveFromCell(AL1,CellIndexOld)
         CALL AddToCell(AL1,CellIndexNew)
         atom%CellIndex=CellIndexNew
      ELSE
         WRITE(6,*) "Db>>negative cell index encountered"
         STOP
      END IF
   END SUBROUTINE RefreshCell
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RefreshBinEdge(bin1,BinRegion,Sgn,species)
   !works on InBin, BinEdgeP,BinEdgeN depending on the location of the new environment atom
   !Last checked: Jan 08, 2011
      TYPE(Bin), POINTER :: bin1
      INTEGER :: BinRegion,Sgn
      INTEGER :: species

      IF (BinRegion==1) THEN
         bin1%BinEdgeP(species)=bin1%BinEdgeP(species)+Sgn
      ELSEIF (BinRegion==-1) THEN
         bin1%BinEdgeN(species)=bin1%BinEdgeN(species)+Sgn
      ELSE
         bin1%InBin(species)=bin1%InBin(species)+Sgn
      END IF
   END SUBROUTINE RefreshBinEdge
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RemoveProcessFromActiveList(prc)
   !when prc does not have any subscribers then prc is removed from the ActiveProcessList
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessList), POINTER :: APLPosition
      
      NumberActiveProcesses=NumberActiveProcesses-1
      APLPosition=>prc%APLPosition
      IF (ASSOCIATED(APLPosition%PrevNeigh)) THEN
         IF (ASSOCIATED(APLPosition%NextNeigh)) THEN
            APLPosition%PrevNeigh%NextNeigh=>APLPosition%NextNeigh
            APLPosition%NextNeigh%PrevNeigh=>APLPosition%PrevNeigh
         ELSE
            NULLIFY(APLPosition%PrevNeigh%NextNeigh)
         END IF
      ELSE
         IF (ASSOCIATED(APLPosition%NextNeigh)) THEN
            ActiveProcessList=>APLPosition%NextNeigh
            NULLIFY(ActiveProcessList%PrevNeigh)
         ELSE
            NULLIFY(ActiveProcessList)
         END IF
      END IF
      NULLIFY(APLPosition%PrevNeigh)
      NULLIFY(APLPosition%NextNeigh)
      DEALLOCATE(APLPosition); NAllottedKMCProcessList=NAllottedKMCProcessList-1
   END SUBROUTINE RemoveProcessFromActiveList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RemoveLocalCfg(atom)
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: atom
      TYPE(Histogram), POINTER :: hist1
      TYPE(BinAtom), POINTER :: binatom1
      TYPE(Bin), POINTER :: bin1

      hist1=>atom%Env%ShortRangeConfig
      atom%BondOrder=0.
      hist1%Population=0
     
      !Get rid of all redundant neighbors
      !binatom1=>hist1%RedundantNeighbors
      !DO WHILE (ASSOCIATED(binatom1))
      !   hist1%RedundantNeighbors=>binatom1%NextNeigh
      !   NULLIFY(binatom1%PrevNeigh)
      !   NULLIFY(binatom1%NextNeigh)
      !   NULLIFY(binatom1%Atom)
      !   DEALLOCATE(binatom1); NAllottedBinAtom=NAllottedBinAtom-1
      !   binatom1=>hist1%RedundantNeighbors
      !END DO
      !Get rid of all bins and their contents
      bin1=>hist1%BinInfo
      DO WHILE(ASSOCIATED(bin1))
         binatom1=>bin1%AL
         DO WHILE (ASSOCIATED(binatom1))
            bin1%AL=>binatom1%NextNeigh
            NULLIFY(binatom1%PrevNeigh)
            NULLIFY(binatom1%NextNeigh)
            NULLIFY(binatom1%Atom)
            DEALLOCATE(binatom1); NAllottedBinAtom=NAllottedBinAtom-1
            binatom1=>bin1%AL
         END DO
         hist1%BinInfo=>bin1%NextNeigh
         NULLIFY(bin1%NextNeigh)
         NULLIFY(bin1%PrevNeigh)
         DEALLOCATE(bin1); NAllottedBin=NAllottedBin-1
         bin1=>hist1%BinInfo
      END DO
   END SUBROUTINE RemoveLocalCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RemoveLocalCfg1(ptr)
     IMPLICIT NONE
     TYPE(Configuration), POINTER :: ptr
     TYPE(Histogram), POINTER :: hist1
     TYPE(BinAtom), POINTER :: binatom1
     TYPE(Bin), POINTER :: bin1

     hist1=>ptr%Env%ShortRangeConfig
     hist1%Population=0
     
     !Get rid of all bins and their contents
     bin1=>hist1%BinInfo
     DO WHILE(ASSOCIATED(bin1))
        binatom1=>bin1%AL
        DO WHILE (ASSOCIATED(binatom1))
           bin1%AL=>binatom1%NextNeigh
           NULLIFY(binatom1%PrevNeigh)
           NULLIFY(binatom1%NextNeigh)
           NULLIFY(binatom1%Atom)
           DEALLOCATE(binatom1); NAllottedBinAtom=NAllottedBinAtom-1
           binatom1=>bin1%AL
        END DO
        hist1%BinInfo=>bin1%NextNeigh
        NULLIFY(bin1%NextNeigh)
        NULLIFY(bin1%PrevNeigh)
        DEALLOCATE(bin1); NAllottedBin=NAllottedBin-1
        bin1=>hist1%BinInfo
     END DO
     DEALLOCATE(ptr); NAllottedConfiguration=NAllottedConfiguration-1
   END SUBROUTINE RemoveLocalCfg1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RemoveFromCell(AL1,CellIndex)
   !Removes an element from the cell without deleting the element from memory
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE (KMCAtomList), POINTER :: AL1
      INTEGER :: CellIndex(3),ix,iy,iz

      ix=CellIndex(1)
      iy=CellIndex(2)
      iz=CellIndex(3)
      
      IF (ASSOCIATED(AL1%PrevNeigh)) THEN
         IF (ASSOCIATED(AL1%NextNeigh)) THEN
            AL1%PrevNeigh%NextNeigh=>AL1%NextNeigh
            AL1%NextNeigh%PrevNeigh=>AL1%PrevNeigh
         ELSE !last atom
            NULLIFY(AL1%PrevNeigh%NextNeigh)
         END IF
      ELSE
         IF (ASSOCIATED(AL1%NextNeigh)) THEN
            Cell(ix,iy,iz)%AL=>AL1%NextNeigh
            NULLIFY(AL1%NextNeigh%PrevNeigh)
         ELSE
            NULLIFY(Cell(ix,iy,iz)%AL)
         END IF
      END IF
      Cell(ix,iy,iz)%NAtomsInCell=Cell(ix,iy,iz)%NAtomsInCell-1
      IF (AL1%Atom%IsMoving) THEN
         Cell(ix,iy,iz)%NMoveInCell=Cell(ix,iy,iz)%NMoveInCell-1
      END IF
      NULLIFY(AL1%PrevNeigh)
      NULLIFY(AL1%NextNeigh)
   END SUBROUTINE RemoveFromCell
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RemoveFromBin(binatom1,bin1)  !B.7
   !Given a pointer structure Histogram called bin1, the sbrtn adds the LatticeSite binatom1 to it
   !NOTE binatom1 IS NOT DEALLOCATED
      IMPLICIT NONE
      TYPE(Bin), POINTER :: bin1
      TYPE(BinAtom), POINTER :: binatom1

      IF (ASSOCIATED(binatom1%PrevNeigh)) THEN
         IF (ASSOCIATED(binatom1%NextNeigh)) THEN
            binatom1%PrevNeigh%NextNeigh=>binatom1%NextNeigh
            binatom1%NextNeigh%PrevNeigh=>binatom1%PrevNeigh
         ELSE
            NULLIFY(binatom1%PrevNeigh%NextNeigh)
         END IF
      ELSE
         IF (ASSOCIATED(binatom1%NextNeigh)) THEN
            NULLIFY(binatom1%NextNeigh%PrevNeigh)
            bin1%AL=>binatom1%NextNeigh
         ELSE
            NULLIFY(bin1%AL)
         END IF
      END IF
      NULLIFY(binatom1%PrevNeigh)
      NULLIFY(binatom1%NextNeigh)
   END SUBROUTINE RemoveFromBin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RemoveFromCfg(Atom,AL1)
   !Extract AL1 out from the configuration Atom%Cfg
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: cfg
      TYPE(KMCAtom), POINTER :: Atom
      TYPE(KMCAtomList), POINTER :: AL1

      IF (.NOT. ASSOCIATED(Atom%Cfg)) THEN
         WRITE(UnitScrap,*) "Db>>RemoveFromCfg..Atom has no config"
         STOP
      END IF
      
      AL1=>SearchFromCfg(Atom)
      IF (.NOT. ASSOCIATED(AL1)) THEN
         WRITE(UnitScrap,*) "Db>>RemoveFromCfg..Atoms location in config not found"
         STOP
      END IF

      cfg=>Atom%Cfg
      NULLIFY(Atom%Cfg)

      IF (ASSOCIATED(AL1%PrevNeigh)) THEN
         IF (ASSOCIATED(AL1%NextNeigh)) THEN
            AL1%NextNeigh%PrevNeigh=>AL1%PrevNeigh
            AL1%PrevNeigh%NextNeigh=>AL1%NextNeigh
         ELSE
            NULLIFY(AL1%PrevNeigh%NextNeigh)
         END IF
      ELSE
         IF (ASSOCIATED(AL1%NextNeigh)) THEN
            NULLIFY(AL1%NextNeigh%PrevNeigh)
            cfg%AtomsWithCfg=>AL1%NextNeigh
         ELSE
            NULLIFY(cfg%AtomsWithCfg)
         END IF
      END IF
      NULLIFY(AL1%NextNeigh)
      NULLIFY(AL1%PrevNeigh)
      cfg%NAtomsWithCfg=cfg%NAtomsWithCfg-1
      !CALL AddCfgToNewList(cfg) !this has been commented to reduce the amont of file writing. As a result, now the latest atoms connected to a configuration are not written
   END SUBROUTINE RemoveFromCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RemoveAllProcessesFromAtom(atom)
   !Removes the process which the atom has subscribed to
   !Note this is a two-way process of breaking the association
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessSubscriptionInfo), POINTER :: prclist,prclist1
      TYPE(KMCAtom), POINTER :: atom,ActiveAtom
      
      prclist=>atom%PrcSubscriberInfo
      DO WHILE (ASSOCIATED(prclist))
         prclist1=>prclist%NextNeigh
         ActiveAtom=>prclist%ActiveAtom
         prc=>prclist%Process
WRITE(Unitscrap,*) "Original atom for which called for subscriber removal:",atom%Index
         CALL RemoveSubscriberAtomFromProcess(ActiveAtom,prc)
         !IF (prclist%IsActive) THEN
         !   prc=>prclist%Process
         !   CALL RemoveSubscriberAtomFromProcess(atom,prc)
         !END IF
         !NULLIFY(prclist%NextNeigh)
         !NULLIFY(prclist%Process)
         !NULLIFY(prclist%ActiveAtom)
         !DEALLOCATE(prclist); NAllottedKMCProcessSubscriptionInfo=NAllottedKMCProcessSubscriptionInfo-1
         prclist=>prclist1
      END DO
   END SUBROUTINE RemoveAllProcessesFromAtom
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RemoveSubscriberAtomFromProcess(activeatom,prc)
   !Removes the active atom and non-active atoms from a process - the activeatom is a subscriber
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCAtom), POINTER :: atom,activeatom
      TYPE(KMCAtomList), POINTER :: AL1
      TYPE(KMCProcessAtom), POINTER :: ReactantAL
      TYPE(KMCProcessSubscriptionInfo), POINTER :: prclist
      REAL :: tol,coord(3),relcoord(3),coordtosearch(3)
      INTEGER :: prcindex,cellindex(3),activeatomindex
      LOGICAL :: IsActive
      
      tol=RadialTolerance
      coord=activeatom%Coord !this is the reference position
      prcindex=prc%Index
      activeatomindex=activeatom%Index
      ReactantAL=>prc%ReactantAL

!WRITE(UnitScrap,*) "Active atom index:",activeatomindex

      DO WHILE (ASSOCIATED(ReactantAL))
         relcoord=ReactantAL%RelInitialPos
         coordtosearch=GetPBCCoord(coord+relcoord)
         
         atom=>SearchAtom(coordtosearch,tol) !this atom is participating in the process

         IF (.NOT. ASSOCIATED(atom)) THEN
            WRITE(6,*) "$Err>> Process atom not found while deleting process subscribers"
            WRITE(6,*) "Process index:",prcindex
            WRITE(6,*) "Active atom index:",activeatom%Index
            WRITE(6,*) "Searching for coordinate:",coordtosearch
            STOP
         END IF
         IF (atom%Cfg%ConfigurationIndex/=ReactantAL%CfgIndex) THEN
            WRITE(6,*) "$Err>> Atom cfg differs from target cfg while deleting process subscribers"
            STOP
         END IF
!  WRITE(UnitScrap,*) "Deleting atom ",atom%Index," from process index ",prcindex
!  WRITE(UnitScrap,*) "Rel coordinate:",relcoord
!  CALL PrintAtomProcessAssignments()
         !now search position in list of processes currently for atom
         prclist=>atom%PrcSubscriberInfo
         
         IsActive=ALL(ABS(relcoord)<RadialTolerance)

!WRITE(UnitScrap,*) "Atom is active:",IsActive

         DO WHILE (ASSOCIATED(prclist)) !one process can show up multiple times; so we want to delete the correct Active type
            IF (prclist%Process%Index==prcindex .AND. (prclist%IsActive.EQV.IsActive) &
              .AND. prclist%ActiveAtom%Index==activeatomindex) EXIT
            prclist=>prclist%NextNeigh
         END DO
         IF (.NOT. ASSOCIATED(prclist)) THEN
            WRITE(6,*) "$Err>> While deleting process subscribers the process was not found in atom"
            WRITE(6,*) "Atom thought to be subscribed:",atom%Index
            WRITE(6,*) "Process index:",prcindex
            STOP
         END IF
         !now delete the entry for the process
         IF (ASSOCIATED(prclist%PrevNeigh)) THEN
            IF (ASSOCIATED(prclist%NextNeigh)) THEN
               prclist%NextNeigh%PrevNeigh=>prclist%PrevNeigh
               prclist%PrevNeigh%NextNeigh=>prclist%NextNeigh
            ELSE
               NULLIFY(prclist%PrevNeigh%NextNeigh)
            END IF
         ELSE
            IF (ASSOCIATED(prclist%NextNeigh)) THEN
               atom%PrcSubscriberInfo=>prclist%NextNeigh
               NULLIFY(prclist%NextNeigh%PrevNeigh)
            ELSE
               NULLIFY(atom%PrcSubscriberInfo)
            END IF
         END IF
         NULLIFY(prclist%PrevNeigh)
         NULLIFY(prclist%NextNeigh)
!WRITE(UnitScrap,*) "Removed process:",prclist%Process%Index
         DEALLOCATE(prclist); NAllottedKMCProcessSubscriptionInfo=NAllottedKMCProcessSubscriptionInfo-1
         ReactantAL=>ReactantAL%NextNeigh
      END DO
      
      AL1=>SearchAtomInProcess(activeatom,prc) !remove active atom from prc%AL structure containing subscribers
      IF (ASSOCIATED(AL1%PrevNeigh)) THEN
         IF (ASSOCIATED(AL1%NextNeigh)) THEN
            AL1%PrevNeigh%NextNeigh=>AL1%NextNeigh
            AL1%NextNeigh%PrevNeigh=>AL1%PrevNeigh
         ELSE
            NULLIFY(AL1%PrevNeigh%NextNeigh)
         END IF
      ELSE
         IF (ASSOCIATED(AL1%NextNeigh)) THEN
            prc%AL=>AL1%NextNeigh
            NULLIFY(AL1%NextNeigh%PrevNeigh)
         ELSE
            NULLIFY(prc%AL)
         END IF
      END IF
      DEALLOCATE(AL1); NAllottedKMCAtomList=NAllottedKMCAtomList-1
      prc%NumberSubscribers=prc%NumberSubscribers-1
      IF (prc%NumberSubscribers==0) CALL RemoveProcessFromActiveList(prc)
   END SUBROUTINE RemoveSubscriberAtomFromProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReplaceProcess(prc1,prc2)
   !exchanges the process information so that the process info of prc2 now belong to prc1 
   !while that of prc1 belongs to prc2
   !in addition if prc1 is stored in NewProcessList then its place is taken by prc2
   
   !NOTE: prc1 is part of KMC_PL and it will stay as a part of the list, 
   !only the contents of prc2 are copied. similary if prc1 is part of NewProcessList
   !then it will continue to be a part
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc1,prc2
      TYPE(KMCProcessAtom), POINTER :: ReactantAL1,ReactantAL2
      TYPE(Configuration), POINTER :: CfgType1,CfgType2
      !TYPE(KMCProcessList), POINTER :: APLPosition1,APLPosition2
      TYPE(KMCProcessList), POINTER :: prclist
      TYPE(KMCAtomList), POINTER :: AL1
      INTEGER :: NumberProcessAtoms1,NumberProcessAtoms2
      !INTEGER :: NumberSubscribers1,NumberSubscribers2
      !INTEGER :: FiringCount1,FiringCount2,MaxFiringAllowed1,MaxFiringAllowed2
      REAL(dp) :: Rate1,Rate2,Frequency1,Frequency2,FwdBarrier1,FwdBarrier2
      !INTEGER :: Index1,Index2,RecordNumber1,RecordNumber2
      !INTEGER :: LastUpdate1,LastUpdate2,FirstAdded1,FirstAdded2

      INTEGER :: ndifferent,CfgDifferent(100),prcindex1,i
      REAL :: RelCoordIniDifferent(300),coord(3),coordtosearch(3),dr(3)
      TYPE(KMCAtom), POINTER :: atom1
      TYPE(KMCProcessSubscriptionInfo), POINTER :: prcsubscribed
      LOGICAL :: EntryDeletion(100),NotFound
      
WRITE(UnitScrap,*) "Replacing process index ",prc1%Index
            
      !NOTE: keep in mind that the number of process atoms can change soon
      !Three possibilities exist for number of process atoms: a) decrease, b) increase, or c) no change
      !if the number decreases/increases we need to inform the subscribed atoms about their change in subscription
      !Step 1. Identify how many atoms are different
      ndifferent=0 !number of atoms that are different
      ReactantAL1=>prc1%ReactantAL !old process
      DO WHILE (ASSOCIATED(ReactantAL1))
         NotFound=.TRUE.
         ReactantAL2=>prc2%ReactantAL !new process to replace the old one
         DO WHILE (ASSOCIATED(ReactantAL2))
            IF (ReactantAL1%CfgIndex==ReactantAL2%CfgIndex) THEN
               dr=ReactantAL1%RelInitialPos-ReactantAL2%RelInitialPos
               IF (ALL(ABS(dr)<RadialTolerance)) THEN !we have found the entry
                  NotFound=.FALSE.
                  EXIT
               END IF
            END IF 
            ReactantAL2=>ReactantAL2%NextNeigh
         END DO
         !now we know if we were able to find the entry in ReactantAL1 or not
         IF (NotFound) THEN
            ndifferent=ndifferent+1
            IF (ndifferent>100) THEN
               WRITE(6,*) "Increase the size of CfgDifferent in ReplaceProcess"
               STOP
            END IF
            RelCoordIniDifferent(3*ndifferent-2:3*ndifferent)=ReactantAL1%RelInitialPos
            CfgDifferent(ndifferent)=ReactantAL1%CfgIndex
            EntryDeletion(ndifferent)=.TRUE. !this entry gets deleted from the process atoms
         END IF
         ReactantAL1=>ReactantAL1%NextNeigh
      END DO
      
      ReactantAL2=>prc2%ReactantAL !new process
      DO WHILE (ASSOCIATED(ReactantAL2))
         NotFound=.TRUE.
         ReactantAL1=>prc1%ReactantAL !old process to be replaced by the new one
         DO WHILE (ASSOCIATED(ReactantAL1))
            IF (ReactantAL1%CfgIndex==ReactantAL2%CfgIndex) THEN
               dr=ReactantAL1%RelInitialPos-ReactantAL2%RelInitialPos
               IF (ALL(ABS(dr)<RadialTolerance)) THEN !we have found the entry
                  NotFound=.FALSE.
                  EXIT
               END IF
            END IF 
            ReactantAL1=>ReactantAL1%NextNeigh
         END DO
         !now we know if we were able to find the entry in ReactantAL2 or not
         IF (NotFound) THEN
            ndifferent=ndifferent+1
            IF (ndifferent>100) THEN
               WRITE(6,*) "Increase the size of CfgDifferent in ReplaceProcess"
               STOP
            END IF
            RelCoordIniDifferent(3*ndifferent-2:3*ndifferent)=ReactantAL2%RelInitialPos
            CfgDifferent(ndifferent)=ReactantAL2%CfgIndex
            EntryDeletion(ndifferent)=.FALSE.
         END IF
         ReactantAL2=>ReactantAL2%NextNeigh
      END DO
      
      !now find all the atoms with subscribscription change to the process and make the changes
      prcindex1=prc1%Index
      IF (ndifferent>0) THEN
         AL1=>prc1%AL !subscribers of AL1
         DO WHILE (ASSOCIATED(AL1))
            coord=AL1%Atom%Coord !coordinate of active atom
            DO i=1,ndifferent

               dr=RelCoordIniDifferent(3*i-2:3*i)
               IF (ALL(ABS(dr)<RadialTolerance)) THEN !this is the active atom
                  WRITE(6,*) "Err> ReplaceProces:: Subscription change for active atom"
                  STOP
               END IF

               coordtosearch=GetPBCCoord(coord+dr)
               atom1=>SearchAtom(coordtosearch,RadialTolerance) !atom with change in subscription

               IF (ASSOCIATED(atom1)) THEN
                  IF (atom1%Cfg%ConfigurationIndex/=CfgDifferent(i)) THEN
                     WRITE(6,*) "Err>> ReplaceProcess: Configurations are different"
                     STOP
                  END IF
               ELSE
                  WRITE(6,*) "Err>> ReplaceProcess: Unable to locate atom"
                  STOP
               END IF


!WRITE(UnitScrap,*) "Removing atom index:",atom1%Index," from subscrip of prc ",prcindex1," with rel init coord",dr
               IF (EntryDeletion(i)) THEN !delete subscription
                  prcsubscribed=>atom1%PrcSubscriberInfo !search for the process
                  DO WHILE (ASSOCIATED(prcsubscribed))
                     IF (prcsubscribed%Process%Index==prcindex1 .AND. &
                        prcsubscribed%ActiveAtom%Index==AL1%atom%Index) EXIT
                     prcsubscribed=>prcsubscribed%NextNeigh
                  END DO
                  IF (prcsubscribed%IsActive) THEN
                     WRITE(6,*) "Err>> Subscription is being removed for an active atom"
                     STOP
                  END IF
                  IF (ASSOCIATED(prcsubscribed)) THEN !delete
                     IF (ASSOCIATED(prcsubscribed%PrevNeigh)) THEN
                        IF (ASSOCIATED(prcsubscribed%NextNeigh)) THEN
                           prcsubscribed%PrevNeigh%NextNeigh=>prcsubscribed%NextNeigh
                           prcsubscribed%NextNeigh%PrevNeigh=>prcsubscribed%PrevNeigh
                        ELSE
                           NULLIFY(prcsubscribed%PrevNeigh%NextNeigh)
                        END IF
                     ELSE
                        IF (ASSOCIATED(prcsubscribed%NextNeigh)) THEN
                           NULLIFY(prcsubscribed%NextNeigh%PrevNeigh)
                           atom1%PrcSubscriberInfo=>prcsubscribed%NextNeigh
                        ELSE
                           NULLIFY(atom1%PrcSubscriberInfo)
                        END IF
                     END IF
                     NULLIFY(prcsubscribed%PrevNeigh)
                     NULLIFY(prcsubscribed%NextNeigh)
                     DEALLOCATE(prcsubscribed); NAllottedKMCProcessSubscriptionInfo=NAllottedKMCProcessSubscriptionInfo-1
                  ELSE !error has occured
                     WRITE(6,*) "Err>> Unable to find process in non-active atom while replacing process"
                     STOP
                  END IF
               ELSE !add subscription to the atom
                  WRITE(6,*) "Adding subscription is not available to the code right now."
                  STOP
               END IF
            END DO
            AL1=>AL1%NextNeigh
         END DO
      END IF
         
      ReactantAL1=>prc1%ReactantAL
      CfgType1=>prc1%CfgType
      !APLPosition1=>prc1%APLPosition
      !AL1=>prc1%AL
      NumberProcessAtoms1=prc1%NumberProcessAtoms
      !NumberSubscribers1=prc1%NumberSubscribers
      !Index1=prc1%Index
      !FiringCount1=prc1%FiringCount
      !MaxFiringAllowed1=prc1%MaxFiringAllowed
      Rate1=prc1%Rate
      Frequency1=prc1%Frequency
      FwdBarrier1=prc1%FwdBarrier
      !RecordNumber1=prc1%RecordNumber
      !LastUpdate1=prc1%LastUpdate
      !FirstAdded1=prc1%FirstAdded
      
      
      ReactantAL2=>prc2%ReactantAL
      CfgType2=>prc2%CfgType
      !APLPosition2=>prc2%APLPosition
      !AL2=>prc2%AL
      NumberProcessAtoms2=prc2%NumberProcessAtoms
      !NumberSubscribers2=prc2%NumberSubscribers
      !Index2=prc2%Index
      !FiringCount2=prc2%FiringCount
      !MaxFiringAllowed2=prc2%MaxFiringAllowed
      Rate2=prc2%Rate
      Frequency2=prc2%Frequency
      FwdBarrier2=prc2%FwdBarrier
      !RecordNumber2=prc2%RecordNumber
      !LastUpdate2=prc2%LastUpdate
      !FirstAdded2=prc2%FirstAdded
      
      !search for the process prc in process list
      !prclist=>NewProcessList
      !DO WHILE (ASSOCIATED(prclist))
      !   IF (prclist%Process%Index==prc1%Index) THEN
      !      WRITE(UnitScrap,*) "Replacing process in NewProcessList ..."
      !      WRITE(UnitScrap,*) "Original process ..."
      !      CALL PrintProcess(prc1,UnitScrap)
      !      WRITE(UnitScrap,*) "New process ..."
      !      CALL PrintProcess(prc2,UnitScrap)
      !      prclist%Process=>prc2
      !      EXIT
      !   END IF
      !   prclist=>prclist%NextNeigh
      !END DO
      
      !now exchange the information...
      prc1%ReactantAL=>ReactantAL2
      IF (CfgType1%ConfigurationIndex/=CfgType2%ConfigurationIndex) THEN
         WRITE(6,*) "Err>> While replacing the process it is found that configuration indices do not match"
         STOP
      END IF
      !prc1%CfgType=>CfgType2
      !prc1%APLPosition=>APLPosition2
      !prc1%AL=>AL2
      prc1%NumberProcessAtoms=NumberProcessAtoms2
      !prc1%NumberSubscribers=NumberSubscribers2
      !prc1%Index=Index2
      !prc1%FiringCount=FiringCount2
      !prc1%MaxFiringAllowed=MaxFiringAllowed2
      prc1%Rate=Rate2
      prc1%Frequency=Frequency2
      prc1%FwdBarrier=FwdBarrier2
      !prc1%RecordNumber=RecordNumber2
      !prc1%LastUpdate=LastUpdate2
      !prc1%FirstAdded=FirstAdded2
      
      prc2%ReactantAL=>ReactantAL1
      !prc2%CfgType=>CfgType1
      !prc2%APLPosition=>APLPosition1
      !prc2%AL=>AL1
      prc2%NumberProcessAtoms=NumberProcessAtoms1
      !prc2%NumberSubscribers=NumberSubscribers1
      !prc2%Index=Index1
      !prc2%FiringCount=FiringCount1
      !prc2%MaxFiringAllowed=MaxFiringAllowed1
      prc2%Rate=Rate1
      prc2%Frequency=Frequency1
      prc2%FwdBarrier=FwdBarrier1
      !prc2%RecordNumber=RecordNumber1
      !prc2%LastUpdate=LastUpdate1
      !prc2%FirstAdded=FirstAdded1
      
      !now get the neighbors right
      !IF (ASSOCIATED(prc1%PrevNeigh)) THEN
      !   IF (ASSOCIATED(prc1%NextNeigh)) THEN
      !      prc2%PrevNeigh=>prc1%PrevNeigh
      !      prc2%PrevNeigh%NextNeigh=>prc2
      !      prc2%NextNeigh=>prc1%NextNeigh
      !      prc2%NextNeigh%PrevNeigh=>prc2
      !   ELSE
      !      prc2%PrevNeigh=>prc1%PrevNeigh
      !      prc2%PrevNeigh%NextNeigh=>prc2
      !      NULLIFY(prc2%NextNeigh)
      !   END IF
      !ELSE
      !   IF (ASSOCIATED(prc1%NextNeigh)) THEN
      !      prc2%NextNeigh=>prc1%NextNeigh
      !      prc2%NextNeigh%PrevNeigh=>prc2
      !      KMC_PL=>prc2
      !      NULLIFY(prc2%PrevNeigh)
      !   ELSE
      !      KMC_PL=>prc2
      !      NULLIFY(prc2%PrevNeigh)
      !      NULLIFY(prc2%NextNeigh)
      !   END IF
      !END IF
      !NULLIFY(prc1%PrevNeigh)
      !NULLIFY(prc1%NextNeigh)
      
      WRITE(UnitScrap,*) "Concluding replacing of process" ! ... here are the details of new process"
      !CALL PrintProcess(prc1,UnitScrap)
      !WRITE(UnitScrap,*) "...and here are details of the old process"
      !CALL PrintProcess(prc2,UnitScrap)
      
   END SUBROUTINE ReplaceProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchAtomInProcess(atom,prc)
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCAtom), POINTER :: atom
      TYPE(KMCAtomList), POINTER :: AL1,SearchAtomInProcess
      INTEGER :: AtomIndex
      LOGICAL :: Found
      
      AtomIndex=atom%Index
      AL1=>prc%AL
      DO WHILE (ASSOCIATED(AL1))
         Found= AL1%Atom%Index==AtomIndex
         IF (Found) EXIT
         AL1=>AL1%NextNeigh
      END DO
      SearchAtomInProcess=>AL1
   END FUNCTION SearchAtomInProcess
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchAtomIndex(Index)
      IMPLICIT NONE
      INTEGER :: Index
      TYPE (KMCAtom), POINTER :: atom,SearchAtomIndex

      atom=>KMC_AL
      DO WHILE (ASSOCIATED(atom))
         IF (atom%Index==Index) EXIT
         atom=>atom%NextNeigh
      END DO
      SearchAtomIndex=>atom
   END FUNCTION SearchAtomIndex
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchAtomCoord(Coord,tol)
   !Searches an atom with a particular coordinate within a tolerance tol
      IMPLICIT NONE
      TYPE(KMCAtom), POINTER :: SearchAtomCoord
      TYPE(KMCAtomList), POINTER :: AL1
      INTEGER :: ix,iy,iz,nx,ny,nz
      INTEGER, DIMENSION(3) :: CellIndex,CellIndex1
      REAL :: Coord(3),tol
      LOGICAL :: Found
      
      CellIndex=GetCellLocation(Coord)
      
      nx=CEILING(tol/KMCCellSize(1))
      ny=CEILING(tol/KMCCellSize(2))
      nz=CEILING(tol/KMCCellSize(3))
      NULLIFY(AL1)
      Found=.FALSE.
      DO ix=-nx,nx
         DO iy=-ny,ny
            DO iz=-nz,nz
               CellIndex1=GetNeighCell(CellIndex,ix,iy,iz)
               AL1=>SearchCell(coord,CellIndex1,tol)
               Found=ASSOCIATED(AL1)
               IF (Found) EXIT
            END DO
            IF (Found) EXIT
         END DO
         IF (Found) EXIT
      END DO
      NULLIFY(SearchAtomCoord)
      IF (ASSOCIATED(AL1)) THEN
         SearchAtomCoord=>AL1%Atom
      END IF
   END FUNCTION SearchAtomCoord
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchCfg(Index)
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: SearchCfg
      INTEGER :: Index

      SearchCfg=>RecordedCfg
      DO WHILE (ASSOCIATED(SearchCfg))
         IF (SearchCfg%ConfigurationIndex==Index) EXIT
         SearchCfg=>SearchCfg%NextNeigh
      END DO
   END FUNCTION SearchCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SearchFromBin(atom,bin1,binatom1)
   !Given a pointer structure Histogram called bin1, the sbrtn searches for the BinContent binatom1
   !with atom atom
      IMPLICIT NONE
      TYPE(Bin), POINTER :: bin1
      TYPE(BinAtom), POINTER :: binatom1
      TYPE(KMCAtom), POINTER :: atom
      INTEGER :: AtomIndex
   
      binatom1=>bin1%AL
      AtomIndex=atom%Index
      DO WHILE (ASSOCIATED(binatom1))
         IF (AtomIndex==binatom1%Atom%Index) RETURN
         binatom1=>binatom1%NextNeigh
      END DO
   END SUBROUTINE SearchFromBin
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchFromCfg(Atom)
   !Locates the atom from the configuration Atom%Cfg
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(KMCAtomList), POINTER :: SearchFromCfg,AL1
      TYPE(KMCAtom), POINTER :: Atom
      INTEGER :: AtomIndex

      AtomIndex=Atom%Index

      AL1=>Atom%Cfg%AtomsWithCfg
      DO WHILE (ASSOCIATED(AL1))
         IF (AL1%Atom%Index==AtomIndex) EXIT
         AL1=>AL1%NextNeigh
      END DO
      SearchFromCfg=>AL1
   END FUNCTION SearchFromCfg
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SearchFromHistogram(atom1,atom,bin1,binatom1)
   !Searches for element binatom1 in bin- bin1 of atom's histogram such that
   !binatom1 points to atom- atom1
      IMPLICIT NONE
      TYPE(BinAtom), POINTER :: binatom1
      TYPE(Bin), POINTER :: bin1
      TYPE(KMCAtom), POINTER :: atom,atom1
      INTEGER :: AtomIndex

      bin1=>atom%Env%ShortRangeConfig%BinInfo
      AtomIndex=atom1%Index
      DO WHILE (ASSOCIATED(bin1))
         binatom1=>bin1%AL
         DO WHILE (ASSOCIATED(binatom1))
            IF (AtomIndex==binatom1%atom%Index) RETURN
            binatom1=>binatom1%NextNeigh
         END DO
         bin1=>bin1%NextNeigh
      END DO
   END SUBROUTINE SearchFromHistogram
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchCfgFromNewList(config)
   !TRUE if config is already present in NewCfgList
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE(Configuration), POINTER :: config
      TYPE(ConfigurationList), POINTER :: ptr
      LOGICAL :: SearchCfgFromNewList
      
      SearchCfgFromNewList=.FALSE.
      ptr=>NewCfgList
      DO WHILE (ASSOCIATED(ptr))
         IF (ptr%Cfg%ConfigurationIndex==config%ConfigurationIndex) THEN
             SearchCfgFromNewList=.TRUE.
             RETURN
         END IF
         ptr=>ptr%NextNeigh
      END DO
    END FUNCTION SearchCfgFromNewList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchProcessFromNewList(prc)
   !TRUE if prc is already present in NewProcessList
   !Last checked: Jan 12, 2011
      IMPLICIT NONE
      TYPE(KMCProcess), POINTER :: prc
      TYPE(KMCProcessList), POINTER :: ptr
      LOGICAL :: SearchProcessFromNewList
      REAL :: tol
      
      SearchProcessFromNewList=.FALSE.
      ptr=>NewProcessList
      tol=RadialTolerance
      DO WHILE (ASSOCIATED(ptr))
         !IF (CompareProcess2Process(prc,ptr%Process,tol)) THEN
         IF (prc%Index==ptr%Process%Index) THEN
            SearchProcessFromNewList=.TRUE.
            RETURN
         END IF
         ptr=>ptr%NextNeigh
      END DO
   END FUNCTION SearchProcessFromNewList
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchCellAtom(atom,cellindex)
   !search atom from cell with index cellindex
   !Last checked: Jan 08, 2011
      IMPLICIT NONE
      TYPE (KMCAtom), POINTER :: atom
      TYPE (KMCAtomList), POINTER :: SearchCellAtom
      INTEGER :: AtomIndex
      INTEGER :: cellindex(3)

      IF (.NOT. ASSOCIATED(atom)) THEN
         WRITE(6,*) "Db>>SearchCell..atom pointer is empty"
         STOP
      END IF

      AtomIndex=atom%Index
      SearchCellAtom=>Cell(cellindex(1),cellindex(2),cellindex(3))%AL
      DO WHILE (ASSOCIATED(SearchCellAtom))
         IF (SearchCellAtom%Atom%Index==AtomIndex) EXIT
         SearchCellAtom=>SearchCellAtom%NextNeigh
      END DO
   END FUNCTION SearchCellAtom
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION SearchCellCoord(coord,cellindex,tol)
   !Searches in cell if any atom has coord within tol
      IMPLICIT NONE
      REAL :: coord(3)
      INTEGER :: cellindex(3)
      TYPE (KMCAtomList), POINTER :: SearchCellCoord
      REAL :: v(3),distance,tol
      
      SearchCellCoord=>Cell(cellindex(1),cellindex(2),cellindex(3))%AL
      DO WHILE (ASSOCIATED(SearchCellCoord))
         v=GetPBCSpacing(SearchCellCoord%Atom%Coord,coord)
         !distance=SQRT(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
         distance=MAXVAL(ABS(v))
         IF (distance<=tol) EXIT
         SearchCellCoord=>SearchCellCoord%NextNeigh
      END DO
   END FUNCTION SearchCellCoord
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SearchBinAfter(binnumber,bin1)
   !tries to find the bin (given by bin1) that has a BinNumber > the input binnumber
   !Note: in some condition such a bin cannot be found
   !Last checked: Jan 08, 2011 -- code is fishy
      IMPLICIT NONE
      TYPE (Bin), POINTER :: bin1
      INTEGER :: binnumber

      IF (.NOT. ASSOCIATED(bin1)) THEN
         WRITE(6,*) "Db>>SearchHistogram.. bin is NULL"
         STOP
      END IF

      IF (bin1%BinNumber>binnumber) RETURN
      DO WHILE (ASSOCIATED(bin1))
         IF (bin1%BinNumber==binnumber) THEN
            EXIT
         ELSE
            IF (ASSOCIATED(bin1%NextNeigh)) THEN
               IF (bin1%NextNeigh%BinNumber<=binnumber) THEN
                  bin1=>bin1%NextNeigh
               ELSE
                  EXIT
               END IF
            ELSE
               EXIT
            END IF
         END IF
      END DO
   END SUBROUTINE SearchBinAfter
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SearchBinBefore(binnumber,bin1)
   !tries to find the bin (given by bin1) that has a BinNumber > the input binnumber
   !Note: in some condition such a bin cannot be found
   !Last checked: Jan 08, 2011 -- code is fishy
      IMPLICIT NONE
      TYPE (Bin), POINTER :: bin1
      INTEGER :: binnumber

      IF (.NOT. ASSOCIATED(bin1)) THEN
         WRITE(6,*) "Db>>SearchHistogram.. bin is NULL"
         STOP
      END IF
    
      IF (bin1%BinNumber<binnumber) RETURN
      DO WHILE (ASSOCIATED(bin1))
         IF (bin1%BinNumber==binnumber) THEN
            EXIT
         ELSE
            IF (ASSOCIATED(bin1%PrevNeigh)) THEN
               IF (bin1%PrevNeigh%BinNumber>=binnumber) THEN
                  bin1=>bin1%PrevNeigh
               ELSE
                  EXIT !IHist is smaller than IHistogram before this
               END IF
            ELSE
               EXIT !end of list before this
            END IF
         END IF
      END DO
   END SUBROUTINE SearchBinBefore
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 END MODULE KMCDbManipulate
