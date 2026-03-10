MODULE mdmpi_domaindecomposition

!Version date: Sep 20, 2011

!Usage:
!mpiexec -n (PRODUCT(DOMAINS)+1 ./amd.par.x
!Nomenclature:
!-------------
!System -- all atoms including size of the system given by SystemSize
!Domain -- atoms contained inside the system minus the buffer
!Buffer -- atoms forming the buffer region around the domain
!Common variables in use:
!------------------------
!Domains - number of domains --Number of domains is equal to the number of slaves
!SystemSize(6) - enforced as 3 times the DomainSize
!DomainSize(6) - All domains are assumed to be of equal size
!CellSize(3) - size of the cell
!NCell(3) - number of cells in the system 
!NCellDomain(3) - number of cells in domain
!NCellThickBuffer(3) - thickness of buffer in terms of number of cells
!CornerCellDomain(3) - cell index 
!DomainNeighbors(6) - domain/slave number for neighbors around the central domain
!AL - Atom list -- domain atoms are listed first followed by the buffer atoms
!MaxBufferCellsPerDomain - maximum buffer cells that can be sent from the neighbor to the center domain
!CellExchangeList(12*MaxBufferCellsPerDomain) - consists of 6 sections of length 2*MaxBufferCellsPerDomain each, 
!  which provides the cell index of the neighbor domain that will be sent (in the odd 
!  numbered elements of CellExchangeList) and the corresponding cell index in the current domain 
!  where this information shall be filled (in the even numbered elements of CellExchangeList). 
!  This information is known to both domains prior to MD simulation.
!CellExchangeListRange(6) - the number of cells to be sent by the neighbor
!BufferCellFilledList(6*MaxBufferCellsPerDomain) - the neighbor domain knows which cell will get filled
! up in the center domain. But the center domain does not know this sequence. Hence, the correct order
! in which the cells are to be filled by the center is intimated by the neighbor
!DomainCellList() -- list of cells that form the domain in AL
!BufferCellList() - of size MaxBufferCellsPerDomain - list of cells that form the buffer in AL
!BufferCellAtomCount(1) - of size MaxBufferCellsPerDomain, which provides the number of atoms in each 
! buffer cell that is being sent by the neighbor there are 6 domains surrounding the center 
! domain in 3D. Each neighbor domain will send the buffer information in a prespecified order. Each 
! time a particular neighbor sends the buffer information, the center domain knows, 
! which cells the information needs to be filled into by looking up even numbered elements of CellList. 
! IMPORTANT: BufferCellAtomCount is used for both sending and receiving the buffer information
! by toggling between sending and receiving.
!BufferAtomSequence() - of size MaxBufferAtomsPerDomain -- sequence of atoms 
! prepared by the neighbor while creating the buffer of center domain
! which shall be used by the center to only cross-check. The size of the 
! message sent is compact. However, which atom index comes from which cell of the 
! neighbor domain can be deciphered using BufferCellAtomCount
!AtomCoordStored(1) [implemented as AtomCoord]- of size 3*MaxBufferAtomsPerDomain -- 
! coordinates of atoms that needs to be sent/received. 
! The cell from/to which the atoms comes from/go to can be decipehered from 
! BufferCellAtomCount. Same is true for BufferAtomSpecies, 
! BufferAtomCharge. Depending on the prespecified order in which the atoms are received they are added into AL
!BufferAllAtomSequence(1) - of size 6*MaxBufferAtomsPerDomain -- sequence of 
! atom indices prepared by neighbors while creating the buffer for the center domain
! which is used to crosscheck if any atom has moved out of the buffer. Each time an 
! atom moves out of a domain or moves out of the buffer. The sequence followed in 
! AL%AtomCoord is disrupted and the VL is consequently affected. Hence it is important 
! to recognize when an atom has moved out from its domain/buffer region
!NeighborDomainSend of a domain I is of size 3^3-1 and it contains the domains neighboring I
   USE VARIABLE_TYPE
   USE neighborlist
   USE db_manipulate
   USE utilities
   USE Species
   USE PotentialPackage
   USE mpi
   IMPLICIT NONE
   !INCLUDE "mpif.h"
   INTEGER, PARAMETER, PRIVATE :: MaxBufferAtomsPerDomain=2000,MaxBufferCellsPerDomain=1000
   INTEGER, PARAMETER, PRIVATE :: MaxAtomsMigratedPerDomain=100,MaxCellsPerNeighbor=10000 
   !128 is set as this is set from the LLFill.. there has to be a better way.
   INTEGER, PRIVATE :: ierr,irank,num_procs,actualrank,TAG
   INTEGER, DIMENSION(MPI_STATUS_SIZE), PRIVATE :: status
   CHARACTER(len=NINT2CHAR) :: intchar
   
   TYPE(SystemContainer), POINTER, PRIVATE :: AL
   !used by to exchange info between master and slaves
   REAL(dp), DIMENSION(3*MaxSlaveNAtoms), PRIVATE :: AtomCoordStored,AtomVelocityStored,AtomMassStored
   INTEGER, DIMENSION(MaxSlaveNAtoms), PRIVATE :: AtomSpeciesStored,AtomSequenceStored,AtomSequence !atomsequence for master and stored for slaves
   REAL(dp), DIMENSION(MaxSlaveNAtoms), PRIVATE :: AtomChargeStored
   REAL(dp), DIMENSION(26*3*MaxAtomsMigratedPerDomain), PRIVATE :: &
      MigratedAtomCoord,MigratedAtomVelocity,MigratedAtomMass!26 neighbours exist &
      !where the atom may migrate to. There would be 26 sections in these arrays &
      !corresponding to different domains
   REAL(dp), DIMENSION(26*MaxAtomsMigratedPerDomain), PRIVATE :: MigratedAtomCharge
   REAL(dp), DIMENSION(3*MaxSlaveNAtoms), PRIVATE :: xp,xc,xf,m
   REAL(dp), PRIVATE :: dtinv,Time,fv
   REAL(dp), PRIVATE :: BoxSize(3)
   INTEGER :: WriteTrajectoryFrequency1
   INTEGER, DIMENSION(6*MaxAtomsMigratedPerDomain), PRIVATE :: MigratedAtomSpecies
   INTEGER, DIMENSION(26), PRIVATE :: MigratedAtomCount,MigratedAtomCountRecv
   INTEGER, DIMENSION(12*MaxCellsPerNeighbor), PRIVATE :: CellExchangeList
   INTEGER, DIMENSION(6*MaxCellsPerNeighbor), PRIVATE :: CellExchangeListRange
   LOGICAL, DIMENSION(3*MaxSlaveNAtoms), PRIVATE :: AtomIsMovingStored
   INTEGER, PRIVATE :: NAtomsStored,NDomainAtoms,NDomainMove(3),NBufferAtoms, &
      NBufferAtoms1,NAtomsToSend,NAtomsToAdd,NAtmSoFar !nbufferatoms1 is the number of atoms in buffer after md
   INTEGER, PRIVATE :: NMigratedAtoms,NAtoms,NAtomsInDomain,TotalCellsInDomain,&
   size_buffersequence
   INTEGER, DIMENSION(3), PRIVATE :: NCellsInDomain,DomainIndex,NCellsBufferThickness,NCell !number of cells in Domain, buffer wall and system
   REAL(dp), PRIVATE :: DomainSize(3),CellSize(3),SystemSize(6)
   INTEGER, DIMENSION(26), PRIVATE :: NeighborDomainSend
   INTEGER, DIMENSION(MaxSlaveNAtoms), PRIVATE :: BufferAllAtomSequence !stores the sequence of atom indices
   INTEGER, DIMENSION(100), PRIVATE :: DomainAssignments,IRankAssingments
   INTEGER, DIMENSION(3), PRIVATE :: SlaveCoord
   LOGICAL, PRIVATE :: MigratedGlobal,MigratedGlobalf
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE md_mpi_master(MDSiml,nmd,Domains,num_procs,irank1,dt,WriteXYZTrajectoryFile,WriteTrajectoryFrequency)
   !Domains is the number of domains along x, y and z direction -- product(Domains)==Nprocessor-1
   !Domain contains the active atoms for a slave
   !Buffer constitutes shell of atoms around the domain
   !System=Buffer+Domain -- note that System for a slave is different from AL
   !num_procs is the number of slave processors doing MD
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      TYPE(LinkedListContainer), POINTER :: LL
      INTEGER, DIMENSION(:), POINTER :: LLListRange,LLList
      INTEGER :: idomain,jdomain,kdomain,i,slave,master,NAtoms,imd,nmd,num_procs,irank1
      INTEGER :: WriteTrajectoryFrequency,errorstatus,cellindx,NextInit
      INTEGER, DIMENSION(3) :: Domains,CornerCell,CornerCellForDomain,slave3d
      REAL(dp) :: dt
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomVelocity,AtomForce
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      LOGICAL :: IsVLInitialized
      CHARACTER(len=100) :: WriteXYZTrajectoryFile
      
      irank=irank1
      write(6,*) "my rank is ", irank!,num_procs,Domains
      WriteTrajectoryFrequency1=WriteTrajectoryFrequency
      master=0
      !IF (irank==master) THEN !master initializes the job
         AL=>MDSiml%AL !does master really receive AL?? can slave also receive the same AL??
         NAtoms=AL%NAtoms
         BoxSize=AL%BoxSize(1:3)-AL%BoxSize(4:6)
         DomainSize=BoxSize/REAL(Domains,dp)
         CellSize=MAXVAL(AL%Potential%MaxPotentialCutOff)
         NCellsInDomain=FLOOR(REAL(FLOOR(BoxSize/CellSize))/Domains)
         CellSize=BoxSize/REAL(NCellsInDomain*Domains,dp)
         
         IF (ASSOCIATED(AL%LL)) CALL DeleteLL(AL%LL)
         NULLIFY(LL)
         CALL LLInitializeCustomized(AL,LL,CellSize,NAtoms) !LL is the linked list for AL
         CALL LLFill(AL,LL,NAtoms)
          
         NCellsBufferThickness=1 !this is the thickness of the buffer wall (actual thickness in each direction will be NCellBufferThickness)
         
         CALL CreateDomainAssignments(Domains) !done more info in the subroutine.
         CornerCell=(/1,1,1/) !corner for the large system
         
         
         DO kdomain=1,Domains(3)
            DO jdomain=1,Domains(2)
               DO idomain=1,Domains(1)
                  
                  !Setup the corner cell for each domain
                  CornerCellForDomain=CornerCell+(/idomain-1,jdomain-1,kdomain-1/)*NCellsInDomain
                  
                  slave=(jdomain-1)*Domains(3)+kdomain+(idomain-1)*Domains(2)*Domains(3) !slave processor for the domain
                  CALL CreateNeighborDomainSend(idomain,jdomain,kdomain,Domains)
                  !creates the array which stores the neighbouring domains
                  !Master creates the array of the iranks which have the neighbouring domains with repesct to each slave
                  CALL SubmitJobToSlave(AL,DomainSize,slave,CornerCellForDomain,& 
                    NeighborDomainSend,Domains,CellSize,idomain,jdomain,kdomain) 
               END DO
            END DO
         END DO
         STOP
         write(*,*) "hi I was master and I have submittied the job properly."
         CALL ReceiveFromSlaves(Domains,AL) !done
         !now all slaves have the atom information that they need and they can perform MD
         !at the end delete LL
         !we need to have the receiving part after this
   END SUBROUTINE md_mpi_master
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE md_mpi(irank1)
      !ELSE !slave takes care of all of the MD steps and communication needed. 
         !Perform the following tasks:
         ! 1. Recieve atom coordinates and velocity from the master
         ! 2. Compute the VL, LL, forces for atoms
         ! 3. Informs a neighbor slave about the updated atom coordinates 
         !   in its domain, for atoms which also happen to be present  -- sending linked list
         !   in buffer of the neighbor slave
         ! 4. Reports any changes in buffer of neighbor slaves -- same as 3.
         ! 5. Reports any changes in domain of neighbor slaves 
         
         ! 1. Recieve atom coordinates, velocity, species ... from the master
         
      TYPE(LinkedListContainer), POINTER :: LL
      INTEGER, DIMENSION(:), POINTER :: LLListRange,LLList
      INTEGER :: i,slave,master,NAtoms,imd,nmd,num_procs,irank1
      INTEGER :: WriteTrajectoryFrequency,errorstatus,cellindx,NextInit
      INTEGER, DIMENSION(3) :: Domains,CornerCell,CornerCellForDomain,slave3d
      REAL(dp) :: dt
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomVelocity,AtomForce
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      LOGICAL :: IsVLInitialized
      CHARACTER(len=100) :: WriteXYZTrajectoryFile
      
      !PRESENTLY THE SLAVES SIMPLY STOP AFTER THEIR JOB -- BUT WHEN FILES NEED TO BE WRITTEN
      ! WITH THE COORDINATES, THE SLAVES NEED TO PAUSE WHAT SHOULD WE DO WITH THIS???
      
         !return
         !CALL Welcome(IsPrint=.FALSE.)
         !write(*,*) "hello there."
         !CALL SetUpSpecies()
         !return
         size_buffersequence=MaxSlaveNAtoms/PRODUCT(Domains)
         CALL MakeSize(AL,MaxSlaveNAtoms)
         CALL PotentialInitialize(AL)
         NCellsBufferThickness=2
         CALL RecieveJobFromMaster(Domains)
      
!AL=>MDSiml%AL 
        !return
        
         NAtoms=AL%NAtoms
         AtomCoord=>AL%AtomCoord
         AtomVelocity=>AL%AtomVelocity
         AtomIsMoving=>AL%AtomIsMoving
         AtomForce=>AL%AtomForce
         NULLIFY(LL)
         IF(.NOT. ASSOCIATED(LL)) CALL LLInitializeCustomized(AL,LL,CellSize,MaxSlaveNAtoms) !LL is the linked list for AL
         !to begin with master has already sent the buffer atoms as well to the slaves
         CALL LLFill(AL,LL,NAtoms=(NDomainAtoms+NBufferAtoms))
          
         NCellsBufferThickness=1
         CALL VLInitialize(AL,NAtoms=MaxSlaveNAtoms)  
         CALL AddVerletList(AL,IsUninitialized=.FALSE.,NAtoms=(NDomainAtoms+NBufferAtoms))
         NCell=LL%NumberCell3D
         slave=IRankAssingments(irank) !slave is the domain number given to the irank
         CornerCell=(/1,1,1/)
         CornerCellForDomain=CornerCell+NCellsInDomain
         CALL SetupArrays(CornerCellForDomain,NCellsInDomain,NCellsBufferThickness,NCell)
         !CALL GetForcesVL(AL,NAtoms,errorstatus)
         IF (errorstatus/=0) STOP
         intchar=INT2CHAR(irank)
         OPEN(UNIT=200,FILE='Output.amd_mpi'//TRIM(intchar))
         fv=SQRT(fev/famu)/fang !conversion factor to get 
         !Angstrom from time*velocity (time=sec, velocity=SQRT(eV/amu))
         ! as well as Velocity from time*accl
         LLList=>LL%List
         LLListRange=>LL%ListRange
         xc=0._dp
         xc(1:3*NAtoms)=AtomCoord(1:3*NAtoms)
         m=AL%AtomMass/dt/dt/fv/fv !now force/m is Ang units
         xp=0._dp
         WHERE (AtomIsMoving(1:3*NDomainAtoms)) !Major problem will occur when atoms are added to the domain
            xp=xc-AtomVelocity*dt*fv
         ELSEWHERE
            xp=xc
            AtomVelocity=0._dp
         END WHERE
         
         xf=xc
         dtinv=0.5_dp/dt/fv !from ang/s to SQRT(eV/amu)
         !nmd=7000 !total number of md steps. I don't know what should be here but for testing purpose
         NAtomsInDomain=NDomainAtoms+NBufferAtoms
         NextInit=0  
         !WRITE(*,*) "BoxSize",AL%BoxSize(1:3)-AL%BoxSize(4:6)
         DO imd=1,nmd  
            NextInit=NextInit+1
            IF(MigratedGlobalf .OR. NextInit==10) THEN !VLinitialize.
               !WRITE(*,*) "migratedor not",MigratedGlobalf,NextInit,irank
               !CALL AddVerletList(AL,IsUninitialized=.FALSE.,NAtoms=(NAtomsInDomain))
               CALL AddVerletListQuick(AL,NAtomsInDomain) !-- is this needed since we are separately updating VL
               NextInit=0
            END IF
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            !MDstep
            CALL GetForcesVL(AL,NAtomsInDomain,errorstatus)
            !DO i=1,SIZE(AL%AtomForce)
            IF(imd==nmd) THEN
               WRITE(*,*) "error",AL%AtomForce(1:6)  
              
            END IF
           ! END DO
            xc(1:3*NAtoms)=AtomCoord(1:3*NAtoms)
            IF (errorstatus/=0) STOP
          !! IF(irank==1) THEN 
           ! WRITE(*,*) "domain",AL%AtomForce(1:15)
           ! stop
         !END IF
            WHERE (AtomIsMoving(1:3*NDomainAtoms)) !NAtoms is defined as a global
               xf=xc+xc-xp+AtomForce/m
            ELSEWHERE
               xf=xc
            END WHERE
       
            Time=Time+dt
         
            WHERE (AtomIsMoving(1:3*NAtoms))
               AtomVelocity=dtinv*(xf-xp)
            ELSEWHERE
               AtomVelocity=0._dp
            END WHERE
      
            xp(1:3*NDomainAtoms)=xc(1:3*NDomainAtoms)
            xc(1:3*NAtoms)=xf(1:3*NAtoms)
            AtomCoord(1:3*NAtoms)=xc(1:3*NAtoms)  
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            
            MigratedGlobalf=.FALSE.
            MigratedGlobal=.FALSE.
            CALL ExchangeAtomInfo(AL,Domains,NCell,SUM(LLListRange), &
            DomainIndex,DomainSize,NBufferAtoms,CornerCellForDomain,&
            NeighborDomainSend,BufferAllAtomSequence,AtomSequenceStored,imd)
            NAtomsInDomain=NDomainAtoms+NBufferAtoms
           ! WRITE(*,*) "NAtomsindomain",NDomainAtoms,NBufferAtoms,irank
           

            !CALL LLFill(AL,LL,NAtoms=(NDomainAtoms+NBufferAtoms))  
           
         END DO
        !write(*,*) "natoms",NDomainAtoms,NBufferAtoms,irank
        !WRITE(*,*) "they are:",AtomCoord(1:3*NAtomsInDomain),irank
         CALL SendToMaster(AtomSequenceStored,Domains)
         !at the end delete AL        
      !END IF
      !CALL MPI_FINALIZE(ierr)
      IF (ierr/=0) THEN
         WRITE(6,*) "$Err>> MPI_MD_DOMAINDECOMPOSITION failed to finalize MPI"
      END IF
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE MDStep(AtomForce,imd,dtinv,Time,dt)
         IMPLICIT NONE
         REAL(dp), DIMENSION(:), POINTER :: AtomForce
         INTEGER :: i,j,iblock,errorstatus,imd
         REAL(dp) :: invtime,Time,dtinv
         CHARACTER(len=30) :: filename
         LOGICAL :: WriteXYZTrajectory
         INTEGER :: atmaxdisp,NextInit=0
         REAL(dp) :: maxdisp,dt
         
         WriteXYZTrajectory=.FALSE.
         NAtoms=NAtomsInDomain
         invtime=dtinv
         NextInit=NextInit+1
      
         iblock=0  
         
         IF(MigratedGlobalf .OR. NextInit==19) THEN !VLinitialize.
            WRITE(*,*) "migratedor not",MigratedGlobalf,NextInit,irank,NAtoms
            CALL AddVerletList(AL,IsUninitialized=.FALSE.,NAtoms=(NAtomsInDomain))
            !CALL AddVerletListQuick(AL,NAtoms) !-- is this needed since we are separately updating VL
            NextInit=0
         END IF
         CALL GetForcesVL(AL,NAtoms,errorstatus)
         xc(1:3*NAtoms)=AtomCoord(1:3*NAtoms)
         IF (errorstatus/=0) STOP
         
         WHERE (AtomIsMoving(1:3*NDomainAtoms)) !NAtoms is defined as a global
            xf=xc+xc-xp+AtomForce/m
         ELSEWHERE
            xf=xc
         END WHERE
       
         iblock=iblock+1

         Time=Time+dt
         IF (iblock==WriteTrajectoryFrequency1) THEN
            iblock=0
            WriteXYZTrajectory=.TRUE.
            !AL%Time=Time
            WHERE (AtomIsMoving(1:3*NAtoms))
               AtomVelocity=invtime*(xf-xp)
            ELSEWHERE
               AtomVelocity=0._dp
            END WHERE
      !   CALL PBC(AL)
!FILE HAS TO BE OPENED BY THE SLAVE IN MD_MPI SUBROUTINE
        ! IF (WriteXYZTrajectory) CALL !WriteXYZ(AL,fileopen=.FALSE.,fileclose=.FALSE.,iunit=200,iprint=1)
       !  IF (WriteChOSTrajectory) CALL AddALToMDHistory(MDSiml,iprint=0)
       !  IF (WriteEnergyTrajectory) THEN
       !     CALL GetKineticEnergy(AL)
       !     WRITE(UNIT=304,FMT='(I8,6ES16.6)') i,AL%Time,MDSiml%Temperature, &
       !          !AL%PotentialEnergy,AL%KineticEnergy,AL%KineticEnergy+AL%PotentialEnergy,GetTemperatur!e(AL)
       !  END IF
         END IF  
         xp(1:3*NDomainAtoms)=xc(1:3*NDomainAtoms)
         xc(1:3*NAtoms)=xf(1:3*NAtoms)
         AtomCoord(1:3*NAtoms)=xc(1:3*NAtoms)
      !some criteria to make sure that the atom hasnot gone out of the huge system
     ! IF(NAtomsInDomain < SUM(LLListRange)) THEN !error exists.. there is some problem
     !    WRITE (*,*) "Error: There seem to be wrong calculation of forces. "
     !    WRITE (*,*) "Atom has covered a large displacement."
     !    RETURN
     ! END IF
      !write(*,*) "difference in coordinates",AtomVelocity(1:3). 

     ! write(*,*) "difference in coordinates",(xc(1:3*NAtoms)-xp(1:3*NAtoms))
      END SUBROUTINE MDStep
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      FUNCTION INT2CHAR(Natural)
         IMPLICIT NONE
         INTEGER :: Natural,DecimalInFocus(NINT2CHAR),DecimalPos
         CHARACTER(len=NINT2CHAR)::INT2CHAR
         CHARACTER :: Ch(NINT2CHAR)
         INTEGER :: i

         IF (Natural<=0) THEN
            WRITE(*,*) "Err>> Only natural numbers allowed in INT2CHAR not ",Natural
            STOP
         END IF
         Ch=' '
         DecimalPos=INT(LOG10(REAL(Natural))) 
         IF (DecimalPos+1>NINT2CHAR) THEN
            WRITE(*,*) 'Exceeded maximum limit handled by INT2CHAR'
            STOP
         END IF
         DO i=1,NINT2CHAR
            DecimalInFocus(i)=MOD(Natural,10**(i+1))/10**i
         END DO
         DO i=1,DecimalPos
            Ch(i)=CHAR(48+DecimalInFocus(DecimalPos+1-i))
         END DO
         Ch(DecimalPos+1)=CHAR(48+MOD(Natural,10))
         INT2CHAR=Ch(1)
         DO i=1,DecimalPos
            INT2CHAR=INT2CHAR(1:len_trim(INT2CHAR))//Ch(i+1)
         END DO
        ! RETURN
      END FUNCTION INT2CHAR
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx	
   END SUBROUTINE md_mpi
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   !think about this carefully
   SUBROUTINE CreateDomainAssignments(Domains)
   !this subroutine assigns the domains to the slaves and keeps track of the domains
   !corresponding to the slaves
   !domainassignments: given a domain, i can find the corresponding irank
   !IRankAssingments: given a irank, i can find the domain
      IMPLICIT NONE 
      INTEGER :: i   
      INTEGER, DIMENSION(3) :: Domains
      
      DomainAssignments=0
      DO i=1,PRODUCT(Domains)
         DomainAssignments(i)=i !right now the most basic arrangement. The i rank is the domain
         IRankAssingments(i)=i
      END DO
   END SUBROUTINE CreateDomainAssignments
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SubmitJobToSlave(AL,DomainSize,slave,CornerCell,NeighborDomainSend,&
      Domains,CellSize,idomain,jdomain,kdomain)
   !master creating the list of atoms to be sent to the slaves
   !accessed only by master
      !IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      TYPE(LinkedListContainer), POINTER :: LL
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomVelocity,AtomMass,AtomCharge
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies,LLList,LLListRange
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      INTEGER, DIMENSION(3) :: CornerCell,Domains,DomainIndex !not declared anywhere
      INTEGER :: MaxAtomPerCell,indx,ixmin,ixmax,iymin,iymax,izmin,izmax,slave,K,i,idomain,jdomain,kdomain
      LOGICAL :: IsBuffer
      REAL(dp), DIMENSION(6) :: SystemSize
      REAL(dp), DIMENSION(3) :: CellSize,DomainSize
      INTEGER, DIMENSION(26) :: NeighborDomainSend
      INTEGER :: cellindx,count1,bufcount
      
      LL=>AL%LL
      NDomainAtoms=0
      NDomainMove=0
      NBufferAtoms=0
      DomainIndex=(/idomain,jdomain,kdomain/)
      NCell=LL%NumberCell3D !number of cells in original (full) system
      MaxAtomPerCell=LL%MaxAtomPerCell
      LLList=>LL%List
      LLListRange=>LL%ListRange
      
      AtomCoord=>AL%AtomCoord
      AtomVelocity=>AL%AtomVelocity
      AtomMass=>AL%AtomMass
      AtomIsMoving=>AL%AtomIsMoving
      AtomSpecies=>AL%AtomSpecies
      AtomCharge=>AL%AtomCharge
      !AtomForce=>AL%AtomForce
      !WRITE(*,*) "master",AtomForce(1:3)
      NAtomsStored=0
      
           ! !CALL GetForcesVL(AL,NAtomsInDomain,errorstatus)
            !WRITE(*,*) "hello",AL%AtomForce(1:18)
      count1=1
      DO i=1,SIZE(LLList)
         IF(LLList(i)>0) THEN
            !cellindx=i/MaxAtomPerCell+1
            AtomSequence(count1)=LLList(i) !stores the global atomindex (master)
            count1=count1+1
         END IF
      END DO
      
      bufcount=0
      izmin=1
      izmax=NCellsInDomain(3)
      iymin=1
      iymax=NCellsInDomain(2)
      ixmin=1
      ixmax=NCellsInDomain(1)
      IsBuffer=.FALSE.
      CALL StoreAtoms()   !store atoms present in the domain to be sent
     
      NDomainAtoms=NAtomsStored
      
      !include buffer region -- FYI - here we are describing the different buffer regions
      !Buffer region 1 - left region
      izmin=1
      izmax=NCellsInDomain(3)
      iymin=1
      iymax=NCellsInDomain(2)
      ixmin=-NCellsBufferThickness(1)+1
      ixmax=0
      IsBuffer=.TRUE.
      CALL StoreAtoms()
      
      !Buffer region 2 - right region
      izmin=1
      izmax=NCellsInDomain(3)
      iymin=1
      iymax=NCellsInDomain(2)
      ixmin=NCellsInDomain(1)+1
      ixmax=NCellsInDomain(1)+NCellsBufferThickness(1) 
      IsBuffer=.TRUE.
      CALL StoreAtoms()
      
      !Buffer region 3 - top region
      izmin=1-NCellsBufferThickness(3)
      izmax=0
      iymin=1
      iymax=NCellsInDomain(2)
      ixmin=1-NCellsBufferThickness(1)
      ixmax=NCellsInDomain(1)+NCellsBufferThickness(1)  
      IsBuffer=.TRUE. 
      CALL StoreAtoms()
      
      !Buffer region 4 - bottom region
      izmin=NCellsInDomain(3)+1
      izmax=NCellsInDomain(3)+NCellsBufferThickness(3)
      iymin=1
      iymax=NCellsInDomain(2)
      ixmin=1-NCellsBufferThickness(1)
      ixmax=NCellsInDomain(1)+NCellsBufferThickness(1)  
      IsBuffer=.TRUE.
      CALL StoreAtoms()
       
      !Buffer region 5 - front region
      izmin=1-NCellsBufferThickness(3)
      izmax=NCellsInDomain(3)+NCellsBufferThickness(3)
      iymin=1-NCellsBufferThickness(2)
      iymax=0
      ixmin=1-NCellsBufferThickness(1)
      ixmax=NCellsInDomain(1)+NCellsBufferThickness(1)   
      IsBuffer=.TRUE.
      CALL StoreAtoms()
      
      !Buffer region 6 - back region
      izmin=1-NCellsBufferThickness(3)
      izmax=NCellsInDomain(3)+NCellsBufferThickness(3)
      iymin=NCellsInDomain(2)+1
      iymax=NCellsInDomain(2)+NCellsBufferThickness(2)
      ixmin=1-NCellsBufferThickness(1)
      ixmax=NCellsInDomain(1)+NCellsBufferThickness(1)  
      IsBuffer=.TRUE.
      CALL StoreAtoms()
       
      NBufferAtoms=NAtomsStored-NDomainAtoms
      
      !write(*,*) "sending to slave buffer sequence",slave,BufferAllAtomSequence(1:bufcount)
     
      CALL MPI_SEND(DomainSize,3,MPI_DOUBLE_PRECISION,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(CellSize,3,MPI_DOUBLE_PRECISION,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(NeighborDomainSend,26,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(NCellsInDomain,3,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(NDomainAtoms,1,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(NBufferAtoms,1,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(NDomainMove,3,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(AtomCoordStored,3*NAtomsStored,MPI_DOUBLE_PRECISION,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(AtomVelocityStored(1:3*NAtomsStored),3*NAtomsStored,MPI_DOUBLE_PRECISION,slave, &
         TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(AtomMassStored,3*NAtomsStored,MPI_DOUBLE_PRECISION,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(AtomIsMovingStored,3*NAtomsStored,MPI_LOGICAL,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(AtomSpeciesStored,NAtomsStored,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(AtomChargeStored,NAtomsStored,MPI_DOUBLE_PRECISION,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(BufferAllAtomSequence,MaxSlaveNAtoms,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(AtomSequenceStored,MaxSlaveNAtoms,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
      SystemSize(4:6)=(2*AL%BoxSize(4:6)-AL%BoxSize(1:3))/Domains+(DomainIndex-1)*BoxSize/Domains
      SystemSize(1:3)=(2*AL%BoxSize(1:3)-AL%BoxSize(4:6))/Domains+(DomainIndex-1)*BoxSize/Domains
      !This system size includes the vacuum region and also takes care of the fact that the domain coordinates are the reference
    
      CALL MPI_SEND(SystemSize,6,MPI_DOUBLE_PRECISION,slave,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(DomainIndex,3,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
      
      !i have sent and received the buffer sequence and the atom sequence. Buffer all atom sequence contin the 
      !buffer atom indices for all the domains. there are equal sections and they can be each accessed.
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE StoreAtoms()
       !  IMPLICIT NONE
         INTEGER :: ix,iy,iz,jstart,jend,j
         INTEGER :: c(3),cellindx,atomindx, sourcedomain(3),i,expectedcellindx(3)
         LOGICAL :: IsMoving(3)
         REAL(dp), DIMENSION(3) :: Transform
         i=1
         
         DO iz=izmin,izmax
            DO iy=iymin,iymax
               DO ix=ixmin,ixmax
                  c=CornerCell-1+(/ix,iy,iz/)
                  expectedcellindx=c
                  c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
                  
                  cellindx=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(2)*NCell(3) !conversion
                  !cellindx=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !conversion
                  !I have changed the LLFactor. I can change it back to what it was earlier and make the necessary 
                  !changes in my code. Right now I have changed the LLFactor from neighborlist.f90
                  IF(LLListRange(cellindx)>0) THEN
                     jstart=(cellindx-1)*MaxAtomPerCell+1
                     jend=jstart+LLListRange(cellindx)-1
                     
                     DO j=jstart,jend
                        atomindx=LLList(j)
                        Transform=0._dp !this is to shift the buffer coordinate to the new domain
                        NAtomsStored=NAtomsStored+1
                        IF(IsBuffer) THEN
                           bufcount=bufcount+1
                           Transform=(expectedcellindx-c)*CellSize
                           BufferAllAtomSequence(bufcount)=atomindx
                        ELSE
                           AtomSequenceStored(NAtomsStored)=atomindx
                        END IF
                        AtomCoordStored(3*NAtomsStored-2:3*NAtomsStored)=AtomCoord(3*atomindx-2:3*atomindx)+Transform
                        AtomVelocityStored(3*NAtomsStored-2:3*NAtomsStored)= &
                           AtomVelocity(3*atomindx-2:3*atomindx)
                        AtomMassStored(3*NAtomsStored-2:3*NAtomsStored)=AtomMass(3*atomindx-2:3*atomindx)
                                               
                        IsMoving=AtomIsMoving(3*atomindx-2:3*atomindx) .AND. .NOT. IsBuffer
                        WHERE (IsMoving) NDomainMove=NDomainMove+1
                        AtomIsMovingStored(3*NAtomsStored-2:3*NAtomsStored)=IsMoving
                        AtomSpeciesStored(NAtomsStored)=AtomSpecies(atomindx)
                        AtomChargeStored(NAtomsStored)=AtomCharge(atomindx)
                     END DO
                  END IF
               END DO
            END DO
         END DO
      END SUBROUTINE StoreAtoms
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE SubmitJobToSlave
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReceiveFromSlaves(Domains,AL) !this is done finally after completing the md steps
      INTEGER :: i,rank,DomainReceived,j
      REAL(dp), DIMENSION(300) :: AtomCoordStored,AtomVelocityStored,AtomMassStored,AtomChargeStored
      INTEGER, DIMENSION(100) :: AtomSpeciesStored
      LOGICAL, DIMENSION(300) :: AtomIsMovingStored
      INTEGER, DIMENSION(3) :: Domains
      TYPE(SystemContainer), POINTER :: AL
      TYPE(LinkedListContainer), POINTER :: LL
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomVelocity,AtomMass,AtomCharge
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies,LLList,LLListRange
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      LL=>AL%LL      
      LLList=>LL%List      
      LLListRange=>LL%ListRange
      AtomCoord=>AL%AtomCoord
      AtomVelocity=>AL%AtomVelocity
      AtomMass=>AL%AtomMass
      AtomIsMoving=>AL%AtomIsMoving
      AtomSpecies=>AL%AtomSpecies
      AtomCharge=>AL%AtomCharge
      NAtoms=0
      DO i=1,Domains(1)*Domains(2)*Domains(3)
         CALL MPI_RECV(DomainReceived,1,MPI_INTEGER,MPI_ANY_SOURCE,TAG,MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(NDomainAtoms,1,MPI_INTEGER,DomainReceived,TAG,MPI_COMM_WORLD,status,ierr)
         NAtoms=NAtoms+NDomainAtoms
         CALL MPI_RECV(AtomSequenceStored(1:NDomainAtoms),NDomainAtoms,MPI_INTEGER,DomainReceived, &
            TAG,MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(AtomCoordStored(1:3*NdomainAtoms),3*NDomainAtoms,MPI_DOUBLE_PRECISION,DomainReceived, & 
            TAG,MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(AtomVelocityStored(1:3*NdomainAtoms),3*NDomainAtoms,MPI_DOUBLE_PRECISION,DomainReceived, &
            TAG,MPI_COMM_WORLD,status,ierr)
         ! I dont think it is required to send or receive this following infromation 
         ! assuming the the atom species or mass or charge is not changed. For a system where
         !these quantities change, they can be received or sent
         !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         !CALL MPI_RECV(AtomMassStored(1:3*NdomainAtoms),3*NDomainAtoms,MPI_DOUBLE_PRECISION,DomainReceived, &
         !   TAG,MPI_COMM_WORLD,status,ierr)
         !CALL MPI_RECV(AtomIsMovingStored(1:3*NDomainAtoms),3*NDomainAtoms,MPI_LOGICAL,DomainReceived, &
         !   TAG,MPI_COMM_WORLD,status,ierr)
         !CALL MPI_RECV(AtomSpeciesStored(1:NDomainAtoms),NDomainAtoms,MPI_INTEGER,DomainReceived,TAG, & 
         !   MPI_COMM_WORLD,status,ierr)
         !CALL MPI_RECV(AtomChargeStored(1:NDomainAtoms),NDomainAtoms,MPI_DOUBLE_PRECISION,DomainReceived,&
         !   TAG,MPI_COMM_WORLD,status,ierr)
         !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         !append to the existing linked list of the initial setting
         
         IF(NDomainAtoms>0) THEN
            DO j=1,NDomainAtoms
               WRITE(*,*) "Atomsequence:",AtomSequenceStored(j),"i am",DomainReceived
               AtomCoord(3*AtomSequenceStored(j)-2:3*AtomSequenceStored(j))= & 
               AtomCoordStored(3*j-2:3*j)
               AtomVelocity(3*AtomSequenceStored(j)-2:3*AtomSequenceStored(j))=& 
               AtomVelocityStored(3*j-2:3*j)
            END Do
            
      !write(*,*) "Finally, the coordinates are:",AtomCoord(AtomSequenceStored(1:NDomainAtoms)
         END IF
         
      END DO
      !write(*,*) "Finally, the velocities are:",AtomVelocity(1:18),NAtoms
      write(*,*) "Finally, the coordinates are:",AtomCoord(1:3*NAtoms),NAtoms
      !finally, the information is retrieved   
   END SUBROUTINE ReceiveFromSlaves
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SendToMaster(AtomSequenceStored,Domains) !this is done finally after completing the md steps
      INTEGER :: status,AtomSequenceStored(MaxSlaveNAtoms),Domains(3)
      REAL(dp), DIMENSION(:), POINTER :: AtomCoordStored,AtomVelocityStored,AtomMassStored,AtomChargeStored
      INTEGER, DIMENSION(:), POINTER :: AtomSpeciesStored
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMovingStored       
      !AL will contain info about buffer also which need not be sent. Select only info of the domain and create a compact list
      AtomCoordStored=>AL%AtomCoord(1:3*NdomainAtoms)
      AtomVelocityStored=>AL%AtomVelocity(1:3*NdomainAtoms)
      AtomMassStored=>AL%AtomMass(1:3*NdomainAtoms)
      AtomIsMovingStored=>AL%AtomIsMoving(1:3*NdomainAtoms)
      AtomSpeciesStored=>AL%AtomSpecies(1:NDomainAtoms)
      AtomChargeStored=>AL%AtomCharge(1:NDomainAtoms) 
      CALL MPI_SEND(irank,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,ierr)  ! this is sent so that master would know from 
      !whom to receive first as it is receiving from any sender
      !WRITE(*,*)"sending to master", NDomainAtoms,NBufferAtoms,NAtoms
      CALL MPI_SEND(NDomainAtoms,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,ierr)      
      CALL MPI_SEND(AtomSequenceStored(1:NDomainAtoms),NDomainAtoms,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,ierr)      
      CALL MPI_SEND(AtomCoordStored(1:3*NDomainAtoms),3*NDomainAtoms,MPI_DOUBLE_PRECISION,0,TAG,MPI_COMM_WORLD,ierr)
      CALL MPI_SEND(AtomVelocityStored(1:3*NDomainAtoms),3*NDomainAtoms,MPI_DOUBLE_PRECISION,0,TAG,MPI_COMM_WORLD,ierr)
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !CALL MPI_SEND(AtomMassStored(1:3*NDomainAtoms),3*NDomainAtoms,MPI_DOUBLE_PRECISION,0,TAG,MPI_COMM_WORLD,ierr)
      !CALL MPI_SEND(AtomIsMovingStored(1:3*NDomainAtoms),3*NDomainAtoms,MPI_LOGICAL,0,TAG,MPI_COMM_WORLD,ierr)
      !CALL MPI_SEND(AtomSpeciesStored(1:NDomainAtoms),NDomainAtoms,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,ierr)
      !CALL MPI_SEND(AtomChargeStored(1:NDomainAtoms),NDomainAtoms,MPI_DOUBLE_PRECISION,0,TAG,MPI_COMM_WORLD,ierr)
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !write(*,*) "sent back to master, I am:",irank
   END SUBROUTINE SendToMaster
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE RecieveJobFromMaster(Domains)
      INTEGER :: master,BufferSize,size_buffersequence,Domains(3)
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomVelocity,AtomMass,AtomCharge
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      
      master=0
      AtomCoord=>AL%AtomCoord
      AtomVelocity=>AL%AtomVelocity
      AtomMass=>AL%AtomMass
      AtomIsMoving=>AL%AtomIsMoving
      AtomSpecies=>AL%AtomSpecies
      AtomCharge=>AL%AtomCharge
       
      !receive atom and system information, to be received just once
      CALL MPI_RECV(DomainSize,3,MPI_DOUBLE_PRECISION,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(CellSize,3,MPI_DOUBLE_PRECISION,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(NeighborDomainSend,26,MPI_INTEGER,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(NCellsInDomain,3,MPI_INTEGER,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(NDomainAtoms,1,MPI_INTEGER,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(NBufferAtoms,1,MPI_INTEGER,master,TAG,MPI_COMM_WORLD,status,ierr)
      NAtoms=NDomainAtoms+NBufferAtoms
      AL%NAtoms=NAtoms
      CALL MPI_RECV(AL%NMove(1:3),3,MPI_INTEGER,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(AtomCoord(1:3*NAtoms),3*NAtoms,MPI_DOUBLE_PRECISION,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(AtomVelocity(1:3*NAtoms),3*NAtoms,MPI_DOUBLE_PRECISION,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(AtomMass(1:3*NAtoms),3*NAtoms,MPI_DOUBLE_PRECISION,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(AtomIsMoving(1:3*NAtoms),3*NAtoms,MPI_LOGICAL,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(AtomSpecies(1:NAtoms),NAtoms,MPI_INTEGER,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(AtomCharge(1:NAtoms),NAtoms,MPI_DOUBLE_PRECISION,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(BufferAllAtomSequence,MaxSlaveNAtoms,MPI_INTEGER,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(AtomSequenceStored,MaxSlaveNAtoms,MPI_INTEGER,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(SystemSize,6,MPI_DOUBLE_PRECISION,master,TAG,MPI_COMM_WORLD,status,ierr)
      CALL MPI_RECV(DomainIndex,3,MPI_INTEGER,master,TAG,MPI_COMM_WORLD,status,ierr)
      AL%BoxSize=SystemSize
     
   END SUBROUTINE RecieveJobFromMaster
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ExchangeAtomInfo(AL,Domains,NCell,LLListRangeSum,DomainIndex,DomainSize,NBufferAtoms,CornerCell,NeighborDomainSend,& 
     BufferAllAtomSequence,AtomSequenceStored,imd)
   !exchanges buffer information and migrated atom info with neighbors and hence, makes the necessary changes
     ! IMPLICIT NONE
      TYPE(LinkedListContainer), POINTER :: LL
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomVelocity,AtomMass,AtomCharge
      REAL(dp), DIMENSION(:), POINTER :: BufferAtomCoord,BufferAtomMass
      REAL(dp), DIMENSION(26*3*MaxAtomsMigratedPerDomain) :: MigratedAtomCoord,&
         MigratedAtomVelocity,MigratedAtomMass,MigratedAtomCoord_p
      REAL(dp), DIMENSION(26*MaxAtomsMigratedPerDomain) ::  MigratedAtomCharge,BufferAtomCharge
      INTEGER, DIMENSION(26*MaxAtomsMigratedPerDomain) :: MigratedAtomSpecies
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies,LLList,LLListRange,BufferAtomSpecies,BufferAtomSequence
      INTEGER, DIMENSION(26*MaxAtomsMigratedPerDomain) :: ListOfMigratedAtom
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving
      INTEGER, DIMENSION(3) :: CornerCell,NCell,DomainIndex,Domains
      REAL(dp), DIMENSION(3) :: coord,DomainSize,xp_migrated
      INTEGER :: i,ra,rb,j,cellindx,count,atomindx,count1,sa,sb,slave,slaveno1(3), &
         atomindx1,NAtoms,ix,iy,iz,actualslaveno1(3)
      INTEGER :: LLListRangeSum,NBufferAtoms,NDomainAtoms1,TEST,NextDomain,imd,MaxAtomPerCell
      REAL(dp) :: Transform(3)
      INTEGER, DIMENSION(26) :: MigratedAtomCount,NeighborDomainSend
      LOGICAL :: AtomMigrated,ChangeInBuffer!,MigratedGlobalf,MigratedGlobal !think about not using change in buffer..
      LOGICAL, DIMENSION(26*MaxAtomsMigratedPerDomain) :: MigratedAtomIsMoving
      INTEGER, DIMENSION(3) :: c,cc,bct !cellindex and buffercellthickness
      INTEGER, DIMENSION(MaxSlaveNAtoms) :: BufferAllAtomSequence,AtomSequenceStored,MigratedAtomIndex
      
      AtomCoord=>AL%AtomCoord !have to initialize a few things
      AtomVelocity=>AL%AtomVelocity
      AtomMass=>AL%AtomMass
      AtomSpecies=>AL%AtomSpecies
      AtomCharge=>AL%AtomCharge
      AtomIsMoving=>AL%AtomIsMoving
      LL=>AL%LL
      NAtoms=AL%NAtoms
      LLList=>LL%List
      LLListRange=>LL%ListRange
      MaxAtomPerCell=LL%MaxAtomPerCell
      slave=IRankAssingments(irank) !given a irank, IRankAssingments sgives the domain
      !report any atom moving out of domain
      !Note: if any has moved out of domain into what was supposed to be buffer then the
      
      ! buffer sequence will be of no use
      CALL LLFill(AL,LL,NAtoms=(NDomainAtoms+NBufferAtoms)) !only updates cell location of Atoms
      MigratedAtomCount=0
       !To find the migrated atoms, check whether the number of atoms in domain has decreased  
      NDomainAtoms1=0
      !Search for migrated atoms
      
      AtomMigrated=.FALSE.
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !=>=> error fro avoiding large displacement. condition to be specified.
     ! WRITE(*,*) 'Initially, Domain,buffer,irank',NDomainAtoms,NBufferAtoms,irank
      !WRITE(*,*) 'Repeat rank:+sum lllist range',NDomainAtoms,NBufferAtoms,irank
      ! IF(NAtomsInDomain < SUM(LLListRange)) THEN !error exists.. there is some problem
      !    WRITE (*,*) "Error: There seem to be wrong calculation of forces. "
      !    WRITE (*,*) "Atom has covered a large displacement.",NAtomsInDomain,SUM(LLListRange)
      !    STOP
      ! END IF   
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    
      ! To know which atoms are buffer and which are domain. No need in actual code, for testing purpose
      !xxxxxxxxxxxxxx
      TEST=0
      IF(TEST>0) THEN
      DO i=1,SIZE(LLList)
         IF(LLList(i)>0) THEN
            cellindx=(i-1)/MaxAtomPerCell+1 !I dont think this is used..
            IF(LLList(i)<=NDomainAtoms) THEN
              ! write(*,*) " atomcount",AtomSequenceStored(LLList(i)),irank,cellindx
               
            ELSE
              ! write(*,*) "bufcount",BufferAllAtomSequence(LLList(i)-NDomainAtoms),irank
               
            END IF
         END IF
      END DO
      END IF
      !return
      !xxxxxxxxxxxxx
      !search for migrated atom. Look for the atoms in the domain. If the number of atoms 
      !decreases then, it has migrated. Better would be to search in the buffer region
      NDomainAtoms1=0 !just for checking purpose
      IF(NDomainAtoms1==0) THEN
      
      bct=NCellsBufferThickness      
      DO iz=1,NCellsInDomain(3)
         DO iy=1,NCellsInDomain(2)
            DO ix=1,NCellsInDomain(1)
               c=CornerCell-1+(/ix,iy,iz/)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               cellindx=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(2)*NCell(3) !conversion
               !cellindx=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(2)*NCell(1) 
               NDomainAtoms1=NDomainAtoms1+LLListRange(cellindx)
            END DO
         END DO
      END DO
      NMigratedAtoms=NDomainAtoms-NDomainAtoms1 ! this is working
            
      END IF
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !searching only in buffer region. What if atom crossed buffer region and is not in 
      !buffer?? then its better to do the above search in the domain. But it may be expensive 
      !if my domain is huge. To ensure that the atom is not moving huge distance
      !a stoping can be done if it moves more than a certain limit
      !To be checked
      NBufferAtoms1=1
      IF(NBufferAtoms1==0) THEN
      
      DO iz=1,6
         ra=2*MaxCellsPerNeighbor*(iz-1)+1
         rb=ra+2*CellExchangeListRange(iz)-1
         DO  iy=ra+1,rb,2 !this gives the cells which belong to buffer region
            cellindx=CellExchangeList(iy)
            NBufferAtoms1=NBufferAtoms1+LLListRange(cellindx)
         END DO
      END DO
      NMigratedAtoms=NBufferAtoms1-NBufferAtoms
      
      END IF
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     ! write(*,*) "migrated",NMigratedAtoms,irank
      count1=0
      IF(NMigratedAtoms>0) THEN
      !delete buffer, unnecessary. as i will have to build it again
         
         AtomMigrated=.TRUE. !this to be sent only to concerned slave
         MigratedGlobal=.TRUE. !to be sent to everyone  
         DO i=1,NDomainAtoms
            coord=AtomCoord(3*i-2:3*i)
            slaveno1=slaveno(coord,DomainIndex,DomainSize,Domains) !gives me location of atom
            !cellindx=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(2)*NCell(3) !conversion
            NextDomain=(slaveno1(2)-1)*Domains(3)+slaveno1(3)+(slaveno1(1)-1)*Domains(3)*Domains(2)
            IF(NextDomain/=irank) THEN ! this means migrated
               actualslaveno1=actualslaveno(coord,DomainIndex,DomainSize,Domains) 
               !this gives me location without periodic boudary condition
               Transform=(slaveno1-actualslaveno1)*DomainSize
               xp_migrated=xp(3*i-2:3*i)
               
               MigratedAtomCount(NextDomain)=MigratedAtomCount(NextDomain)+1
               count=MigratedAtomCount(NextDomain)
               count1=count1+1
               sa=MaxAtomsMigratedPerDomain*(NextDomain-1) !slave is my domain or irank 
               count=count+sa
               ListOfMigratedAtom(count1)=i 
               MigratedAtomCoord(3*count-2:3*count)=coord+Transform
               MigratedAtomCoord_p(3*count-2:3*count)=xp_migrated+Transform
               MigratedAtomVelocity(3*count-2:3*count)=AtomVelocity(3*i-2:3*i)
               MigratedAtomMass(3*count-2:3*count)=AtomMass(3*i-2:3*i)
               MigratedAtomSpecies(count)=AtomSpecies(i)
               MigratedAtomCharge(count)=AtomCharge(i)
               MigratedAtomIndex(count)=AtomSequenceStored(i) !for slave
               
            END IF
         END DO
         !now compact AL to account for migrated atoms
         CALL Sort(ListOfMigratedAtom(1:NMigratedAtoms),NMigratedAtoms) !to go through the AL in the correct order 
         !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         !i havent moved 2 atoms so i have to be sure this work
         IF(NMigratedAtoms>=2) THEN
            DO i=1,NMigratedAtoms-1
               atomindx=ListOfMigratedAtom(i)! for the domain
               atomindx1=ListOfMigratedAtom(i+1)
               AtomCoord(3*(atomindx-i+1)-2:3*(atomindx1-i-1))=AtomCoord(3*(atomindx+1)-2:3*(atomindx1-1))
               xp(3*(atomindx-i+1)-2:3*(atomindx1-i-1))=xp(3*(atomindx+1)-2:3*(atomindx1-1))
               AtomVelocity(3*(atomindx-i+1):3*(atomindx1-i-1))=AtomVelocity(3*(atomindx+1):3*(atomindx1-1))
               AtomMass(3*(atomindx-i+1):3*(atomindx1-i-1))=AtomMass(3*(atomindx+1):3*(atomindx1-1))
               AtomIsMoving(3*(atomindx-i+1):3*(atomindx1-i-1))=AtomIsMoving(3*(atomindx+1):3*(atomindx1-1))
               AtomSpecies(atomindx-i+1:atomindx1-i-1)=AtomSpecies(atomindx+1:atomindx1-1)
               AtomCharge(atomindx-i+1:atomindx1-i-1)=AtomCharge(atomindx+1:atomindx1-1)
               AtomSequenceStored(atomindx-i+1:atomindx1-i-1)= & 
                  AtomSequenceStored(atomindx+1:atomindx1-1) 
               !this part is not checked as only one atom has made to migrated. Migrate 2
               !atoms if next level which I would do 
            END DO
         END IF
         !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
         IF (ListOfMigratedAtom(NMigratedAtoms)/=NAtoms) THEN
            atomindx1=ListOfMigratedAtom(NMigratedAtoms)  
            AtomCoord(3*(atomindx1-NMigratedAtoms+1)-2:3*(NAtoms-NMigratedAtoms))= &
               AtomCoord(3*(atomindx1+1)-2:3*NAtoms)
            xp(3*(atomindx1-NMigratedAtoms+1)-2:3*(NAtoms-NMigratedAtoms))= &
               xp(3*(atomindx1+1)-2:3*NAtoms)
            AtomVelocity(3*(atomindx1-NMigratedAtoms+1)-2:3*(NAtoms-NMigratedAtoms))= &
               AtomVelocity(3*(atomindx1+1)-2:3*NAtoms)
            AtomMass(3*(atomindx1-NMigratedAtoms+1)-2:3*(NAtoms-NMigratedAtoms))= &
               AtomMass(3*(atomindx1+1)-2:3*NAtoms)
            AtomIsMoving(3*(atomindx1-NMigratedAtoms+1)-2:3*(NAtoms-NMigratedAtoms))= &
               AtomIsMoving(3*(atomindx1+1)-2:3*NAtoms)
            AtomSpecies(atomindx1-NMigratedAtoms+1:NAtoms-NMigratedAtoms)=AtomSpecies(atomindx1+1:NAtoms)
            AtomCharge(atomindx1-NMigratedAtoms+1:NAtoms-NMigratedAtoms)=AtomCharge(atomindx1+1:NAtoms)
            AtomSequenceStored(atomindx1-NMigratedAtoms+1:NAtoms-NMigratedAtoms)= &
               AtomSequenceStored(atomindx1+1:NAtoms)
            !checked and working
         END IF
            
      END IF
      
      NDomainAtoms1=NDomainAtoms
      IF(AtomMigrated) THEN
         NDomainAtoms=NDomainAtoms-NMigratedAtoms              
      END IF
      CALL SendAndReceiveMigratedAtoms()
      !Initialization 
      !exclude NBufferAtoms
      IF(AtomMigrated) THEN
         MigratedGlobalf=.TRUE.
      END IF
      NBufferAtoms=0
      IF (AtomMigrated) THEN
         NAtoms=NDomainAtoms+NBufferAtoms
         LL=>AL%LL
         LLList=>LL%List
         LLListRange=>LL%ListRange
         CALL LLFill(AL,LL,NAtoms=(NDomainAtoms+NBufferAtoms)) 
        ! write(*,*) "initializing after migration"
         !to check whether initialization is done properly or not
         iz=0
         DO i=1,SIZE(LLList)
            IF(LLList(i)>0) THEN 
               cellindx=(i-1)/MaxAtomPerCell+1
            !   write(*,*) "atoms present are:",LLList(i),cellindx,irank
               iz=iz+1
            END IF
         END DO
      END IF
      
      IF(.NOT. AtomMigrated) THEN
         NAtoms=NDomainAtoms+NBufferAtoms
         LL=>AL%LL
         LLList=>LL%List
         LLListRange=>LL%ListRange
         CALL LLFill(AL,LL,NAtoms=(NDomainAtoms+NBufferAtoms))
        ! write(*,*) "initializing after ",SUM(LLListRange)
      END IF
            
      BufferAllAtomSequence=0 
      CALL SendAndReceiveBufferAtoms() 
      !I am updating completly
      !Reasons are: 1) Once the atom is migrated, I have to make sure that it is deleted
      !from all the buffersequences. I have to in short change the buffer sequences(compress)
      ! and all before actually reaching to send and receive buffer where I will receive new 
      !sequences and this may lead to confusions. Moreover, the infromation about
      !migrated has to be informed to everyone and this may be an issue.
      !Either I can involve master in between or just do it between the slaves which would 
      !be complex
      
      !return
      !Initialization after buffer change
      !NAtoms=NDomainAtoms+NBufferAtoms
      !CALL AddVerletList(AL,NAtoms=(NDomainAtoms+NBufferAtoms))
      !LL=>AL%LL
      !LLList=>LL%List
      !LLListRange=>LL%ListRange
      !CALL LLFill(AL,LL,NAtoms=(NDomainAtoms+NBufferAtoms))
      !write(*,*) "initializing after buffer"
      xc=0
      xc(1:3*NAtoms)=AtomCoord(1:3*NAtoms)
     
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      FUNCTION slaveno(Coord,DomainIndex,DomainSize,Domains) !returns me the new domain of the migrated atom
         INTEGER :: slaveno(3) 
         REAL(dp), DIMENSION(3) :: Coord,DimlessCoord,DomainSize,MinCoord
         REAL(dp), DIMENSION(6) :: BoxSize
         INTEGER, DIMENSION(3) :: c,DomainIndex,CoordIndex,Domains
         MinCoord=(DomainIndex-1)*DomainSize
         CoordIndex=DomainIndex+FLOOR((Coord-MinCoord)/DomainSize)
        !peridic boundary condition
         CoordIndex=CoordIndex-FLOOR(REAL(CoordIndex-1)/REAL(Domains))*Domains
         slaveno=CoordIndex !slaveno=(CoordIndex(2)-1)*Domains(3)+CoordIndex(3)+(CoordIndex(1)-1)*Domains(3)*Domains(2)
       END FUNCTION slaveno
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      FUNCTION actualslaveno(Coord,DomainIndex,DomainSize,Domains) !actual domainw/o periodic bdry condition
      !for migrated atoms
         INTEGER :: actualslaveno(3)
         REAL(dp), DIMENSION(3) :: Coord,DimlessCoord,DomainSize,MinCoord
         REAL(dp), DIMENSION(6) :: BoxSize
         INTEGER, DIMENSION(3) :: c,DomainIndex,CoordIndex,Domains
         MinCoord=(DomainIndex-1)*DomainSize
         CoordIndex=DomainIndex+FLOOR((Coord-MinCoord)/DomainSize)
        actualslaveno=CoordIndex !actualslaveno=(CoordIndex(2)-1)*Domains(3)+CoordIndex(3)+(CoordIndex(1)-1)*Domains(3)*Domains(2)
       END FUNCTION actualslaveno
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE SendAndReceiveMigratedAtoms()
         INTEGER :: i
         i=1
         DO WHILE(i<26)
            CALL SendMigrated(i)
            CALL RecvMigrated(i+1) 
            CALL SendMigrated(i+1)
            CALL RecvMigrated(i)
            i=i+2
         END DO 
               
      END SUBROUTINE SendAndReceiveMigratedAtoms
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE SendAndReceiveBufferAtoms()      
         INTEGER :: i
         !Steps:
         !1) send to left and recv from right
         !2) send to right and receive from left
         !3) send to back and receive from front
         !4) send to front and receive from back
         !5) send to bottom and receive from top
         !6) send to top and receive from bottom
         !1=>left, 2=>right, 3=>back, 4=> front, 5=> bottom, 6=> top 
         i=1
         DO WHILE(i<6)
            CALL SendBuffer(i)
            CALL ReceiveBuffer(i+1)
            CALL SendBuffer(i+1)
            CALL ReceiveBuffer(i)
            i=i+2
         END DO
      END SUBROUTINE SendAndReceiveBufferAtoms
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE SendMigrated(slaveno)
         INTEGER :: slaveno
         slave=NeighborDomainSend(slaveno) !gives the irank which has the neighbouring buffer
         NAtomsToSend=MigratedAtomCount(slave)
         MigratedAtomCount(slave)=0
         
         CALL MPI_SEND(NAtomsToSend,1,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr) 
         CALL MPI_SEND(MigratedGlobal,1,MPI_LOGICAL,slave,TAG,MPI_COMM_WORLD,ierr)  
         IF (NAtomsToSend>0) THEN
            sa=MaxAtomsMigratedPerDomain*(slave-1)+1
            sb=sa-1+NAtomsToSend
            CALL MPI_SEND(MigratedAtomCoord(3*sa-2:3*sb),3*NAtomsToSend,MPI_DOUBLE_PRECISION,slave, &
               TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(MigratedAtomCoord_p(3*sa-2:3*sb),3*NAtomsToSend,MPI_DOUBLE_PRECISION,slave, &
               TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(MigratedAtomVelocity(3*sa-2:3*sb),3*NAtomsToSend,MPI_DOUBLE_PRECISION,slave, &
               TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(MigratedAtomMass(3*sa-2:3*sb),3*NAtomsToSend,MPI_DOUBLE_PRECISION,slave, &
               TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(MigratedAtomSpecies(sa:sb),NAtomsToSend,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(MigratedAtomCharge(sa:sb),NAtomsToSend,MPI_DOUBLE_PRECISION,slave,TAG,MPI_COMM_WORLD,ierr)
            MigratedAtomIsMoving(3*sa-2:3*sb)=.TRUE.
            CALL MPI_SEND(MigratedAtomIsMoving(3*sa-2:3*sb),3*NAtomsToSend,MPI_LOGICAL,slave, &
               TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(MigratedAtomIndex(sa:sb),NAtomsToSend,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
            AtomMigrated=.TRUE.
            !WRITE(*,*) "atom migrated to:",slave,"from",irank,MigratedAtomCoord(3*sa-2:3*sb)
         END IF 
      END SUBROUTINE SendMigrated
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE RecvMigrated(slaveno)
         INTEGER :: NAtomsRecv,NDomainAtoms1,slaveno
         
         slave=NeighborDomainSend(slaveno)
         CALL MPI_RECV(NAtomsRecv,1,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(MigratedGlobal,1,MPI_LOGICAL,slave,TAG,MPI_COMM_WORLD,status,ierr)
        
         IF(MigratedGlobal) THEN
            MigratedGlobalf=MigratedGlobal
         END IF
         MigratedAtomCountRecv(slave)=NAtomsRecv
         IF (NAtomsRecv>0) THEN
           
            NDomainAtoms1=NDomainAtoms+NAtomsRecv
            !make way for new atoms
            AtomCoord(3*(NDomainAtoms1+1)-2:3*(NDomainAtoms1+NBufferAtoms))= & 
               AtomCoord(3*(NDomainAtoms+1)-2:3*(NDomainAtoms+NBufferAtoms))
            AtomMass(3*(NDomainAtoms1+1)-2:3*(NDomainAtoms1+NBufferAtoms))= & 
               AtomMass(3*(NDomainAtoms+1)-2:3*(NDomainAtoms+NBufferAtoms))
            AtomSpecies((NDomainAtoms1+1):(NDomainAtoms1+NBufferAtoms))= & 
               AtomSpecies((NDomainAtoms+1):(NDomainAtoms+NBufferAtoms))
            AtomCharge((NDomainAtoms1+1):(NDomainAtoms1+NBufferAtoms))= & 
               AtomCharge((NDomainAtoms+1):(NDomainAtoms+NBufferAtoms))
            AtomIsMoving(3*(NDomainAtoms1+1)-2:3*(NDomainAtoms1+NBufferAtoms))= & 
               AtomIsMoving(3*(NDomainAtoms+1)-2:3*(NDomainAtoms+NBufferAtoms)) 
            
            CALL MPI_RECV(AtomCoord(3*NDomainAtoms+1:3*NDomainAtoms1),3*NAtomsRecv,MPI_DOUBLE_PRECISION, &
               slave,TAG,MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(xp(3*NDomainAtoms+1:3*NDomainAtoms1),3*NAtomsRecv,MPI_DOUBLE_PRECISION, &
               slave,TAG,MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(AtomVelocity(3*NDomainAtoms+1:3*NDomainAtoms1),3*NAtomsRecv,MPI_DOUBLE_PRECISION, &
               slave,TAG,MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(AtomMass(3*NDomainAtoms+1:3*NDomainAtoms1),3*NAtomsRecv,MPI_DOUBLE_PRECISION,slave, &
               TAG,MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(AtomSpecies(NDomainAtoms+1:NDomainAtoms1),NAtomsRecv,MPI_INTEGER,slave, &
               TAG,MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(AtomCharge(NDomainAtoms+1:NDomainAtoms1),NAtomsRecv,MPI_DOUBLE_PRECISION,slave, &
               TAG,MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(AtomIsMoving(3*NDomainAtoms+1:3*NDomainAtoms1),3*NAtomsRecv,MPI_LOGICAL,slave, &
               TAG,MPI_COMM_WORLD,status,ierr)            
            CALL MPI_RECV(AtomSequenceStored(NDomainAtoms+1:NDomainAtoms1),NAtomsRecv,MPI_INTEGER,slave, &
               TAG,MPI_COMM_WORLD,status,ierr)
            AtomMigrated=.TRUE.  
           ! write(*,*) "received migrated atoms from:",AtomSequenceStored(1:NDomainAtoms1),irank,slave
            NDomainAtoms=NDomainAtoms1
         END IF
         
      END SUBROUTINE RecvMigrated
      !send and receive migrated is done
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE SendBuffer(slaveno) 
         !send buffer info (info about the domain atoms which may lie in the buffer of other neighbours
         INTEGER, DIMENSION(MaxCellsPerNeighbor) :: BufferCellindx,BufferCellRange
         REAL(dp), DIMENSION(3*MaxBufferAtomsPerDomain) :: BufferAtomCoord,BufferAtomMass,BufferAtomCharge
         INTEGER, DIMENSION(MaxBufferAtomsPerDomain) :: BufferAtomSpecies
         INTEGER, DIMENSION(MaxSlaveNAtoms) :: BufferAtomSequence
        ! REAL(dp), DIMENSION(100) :: BufferAtomCharge
         INTEGER :: NAtomsToAdd,NAtomsFilled,atomindx,slaveno,cellno,sumtemp,NAtomsToStore,NAtomsFilled1
         INTEGER :: cellindx_neighbour,c(3),c_neighbour(3),slave3d(3),irank3d(3), &
         ExpectedSlave3d(3),ExpectedSlave,slaveno3d(3)
         REAL(dp), DIMENSION(3) :: Transform1
         
         slave=NeighborDomainSend(slaveno)!this is the corresponding domain to shich it is to be sent
         ra=2*MaxCellsPerNeighbor*(slaveno-1) !setting range for the cells to be taken
         !cellexchangelist stores the cellindex of 
         ! both cells belonging to my domain and to the neighbor domain.So 2 factor comes
         rb=ra+2*CellExchangeListRange(slaveno)
         atomindx=0
         NAtomsFilled=0
         cellno=0
         irank3d=DomainIndex
         !think about this
         slave3d(1)=(slave-1)/(Domains(2)*Domains(3))+1
         slave3d(3)=MOD((slave-1),Domains(3))+1
         slave3d(2)=(slave-1)/Domains(3)+1
         
         !slave3d(3)=(slave-1)/(Domains(2)*Domains(1))+1
         !slave3d(1)=MOD((slave-1),Domains(1))+1
         !slave3d(2)=(slave-1)/Domains(1)+1
         slave3d=slave3d-Domains*FLOOR(REAL(slave3d-1)/REAL(Domains))!periodic boundary condition
         !slave3d=MOD(slave3d-1,Domains)+1
         ExpectedSlave3d=ExpectedDomain(slaveno)
         !ExpectedSlave=(ExpectedSlave3d(2)-1)*Domains(1)+ExpectedSlave3d(1)+(ExpectedSlave3d(3)-1)*Domains(1)*Domains(2)
         ExpectedSlave=(ExpectedSlave3d(2)-1)*Domains(3)+ExpectedSlave3d(3)+(ExpectedSlave3d(1)-1)*Domains(3)*Domains(2)
         
         Transform1=(slave3d-ExpectedSlave3d)*DomainSize
         cellno=0
         DO i=ra+1,rb,2 !every odd gives the cellindex which i have to send
            cellindx=CellExchangeList(i) !cell from which atoms need to be sent
            IF(LLListRange(cellindx)>0) THEN !this implies there is an atom in the cell
               sa=(cellindx-1)*LL%MaxAtomPerCell+1 !for storing the atoms information, the range in the array showing the atom location in LLList
               NAtomsToAdd=LLListRange(cellindx) !for this cell only   
               cellno=cellno+1
               BufferCellindx(cellno)= CellExchangeList(i+1)!This is the cell in the neighbor where I have to paste
               
               BufferCellRange(cellno)=NAtomsToAdd
               sb=sa-1+NAtomsToAdd
               atomindx=0
               DO j=sa,sb
                  atomindx=LLList(j)
                  NAtomsFilled=NAtomsFilled+1
                  IF(atomindx<=NDomainAtoms) THEN
                     BufferAtomSequence(NAtomsFilled)=AtomSequenceStored(atomindx) !there seems to be some initialization error
                    ! write(*,*) "hi, i have entered this loop",irank,slave,BufferAtomSequence(NAtomsFilled)
                  ELSEIF(atomindx>NDomainAtoms) THEN                      
                     BufferAtomSequence(NAtomsFilled)=BufferAllAtomSequence(atomindx-NDomainAtoms)                        
                    ! write(*,*) "hi, i have entered this loopbuffer",irank,slave,BufferAtomSequence(NAtomsFilled)
                  END IF
                  BufferAtomCoord(3*(NAtomsFilled)-2:3*(NAtomsFilled))=AtomCoord(3*atomindx-2:3*atomindx)+Transform1
                  BufferAtomMass(3*(NAtomsFilled+1)-2:3*(NAtomsFilled+1))=AtomMass(3*atomindx-2:3*atomindx)
                  BufferAtomSpecies(NAtomsFilled)=AtomSpecies(atomindx)
                  BufferAtomCharge(NAtomsFilled)=AtomCharge(atomindx)
                  
               END DO
            END IF
         END DO 
       
         IF(NAtomsFilled>MaxBufferAtomsPerDomain) THEN
            WRITE(*,*) "Increase the size of the buffer arrays"
            STOP
         END IF
         CALL MPI_SEND(cellno,1,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(BufferCellRange(1:cellno),cellno,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(BufferCellindx(1:cellno),cellno,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(NAtomsFilled,1,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(BufferAtomSequence(1:NAtomsFilled),NAtomsFilled,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(BufferAtomCoord(1:3*NAtomsFilled),3*NAtomsFilled,MPI_DOUBLE_PRECISION,slave,TAG, &
           MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(AtomMassStored(1:3*NAtomsFilled),3*NAtomsFilled,MPI_DOUBLE_PRECISION,slave, &
            TAG,MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(BufferAtomSpecies(1:NAtomsFilled),NAtomsFilled,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,ierr)
         CALL MPI_SEND(BufferAtomCharge(1:NAtomsFilled),NAtomsFilled,MPI_DOUBLE_PRECISION,slave,TAG,MPI_COMM_WORLD,ierr)
        ! WRITE(*,*) "cellno",SUM(LLListRange),irank,slave
         
      END SUBROUTINE SendBuffer
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE ReceiveBuffer(slaveno)
      !receives the buffer info send by neighbor slave/domain given by slave
         INTEGER, DIMENSION(MaxCellsPerNeighbor) :: BufferCellindx,BufferCellRange
         REAL(dp), DIMENSION(3*MaxBufferAtomsPerDomain) :: BufferAtomCoord,BufferAtomMass,BufferAtomCharge
         INTEGER, DIMENSION(MaxBufferAtomsPerDomain) :: BufferAtomSpecies
         INTEGER, DIMENSION(MaxSlaveNAtoms) :: BufferAtomSequence
         INTEGER :: NAtomsToAdd,NAtomsFilled,atomindx,slaveno,cellno,sumtemp,NAtomsToStore,NAtomsFilled1
         INTEGER :: cellindx_neighbour,c(3),c_neighbour(3),AtomLocation,AtomUpdated
         LOGICAL :: ChangeInBuffer=.FALSE.
         slave=NeighborDomainSend(slaveno)
         
         CALL MPI_RECV(cellno,1,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(BufferCellRange(1:cellno),cellno,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(BufferCellindx(1:cellno),cellno,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(NAtomsFilled,1,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,status,ierr)
         
         CALL MPI_RECV(BufferAtomSequence(1:NAtomsFilled),NAtomsFilled,MPI_INTEGER, & 
         slave,TAG,MPI_COMM_WORLD,status,ierr)

         CALL MPI_RECV(BufferAtomCoord(1:3*NAtomsFilled),3*NAtomsFilled,MPI_DOUBLE_PRECISION,slave,TAG, &
         MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(BufferAtomMass(1:3*NAtomsFilled),3*NAtomsFilled,MPI_DOUBLE_PRECISION,slave,TAG, &
            MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(BufferAtomSpecies(1:NAtomsFilled),NAtomsFilled,MPI_INTEGER,slave,TAG,MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(BufferAtomCharge(1:NAtomsFilled),NAtomsFilled,MPI_DOUBLE_PRECISION,slave,TAG,MPI_COMM_WORLD,status,ierr)
        ! WRITE(*,*) "slave",slave,"receiver", irank,"cellno",BufferAtomCoord(1:3*NAtomsFilled),"NAtomsFilled",NAtomsFilled
   
         
         IF(NAtomsFilled>0) THEN
            NBufferAtoms=NBufferAtoms+NAtomsFilled  
            BufferAllAtomSequence(NBufferatoms-NAtomsFilled+1:NBufferAtoms)= &
               BufferAtomSequence(1:NAtomsFilled)
            NAtoms=NDomainAtoms+NBufferAtoms
            AtomCoord(3*(NAtoms-NAtomsFilled+1)-2:3*(NAtoms))=BufferAtomCoord(1:3*NAtomsFilled)
            AtomMass(3*(NAtoms-NAtomsFilled+1)-2:3*(NAtoms))=BufferAtomMass(1:3*NAtomsFilled)
            AtomSpecies(NAtoms-NAtomsFilled+1:NAtoms)=BufferAtomSpecies(1:NAtomsFilled)
            AtomCharge(NAtoms-NAtomsFilled+1:NAtoms)=BufferAtomCharge(1:NAtomsFilled)
            !buffercellindx can also be found from cellexchangelist
            !Globally, I need some count of how mny atoms arefilled 
            AtomUpdated=0
            
            CALL LLFill(AL,LL,NAtoms=(NDomainAtoms+NBufferAtoms)) 
            
            
            !updating LL without doing LLFill there is some problem. so updating ll 
            !using llfill
            !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            IF(AtomUpdated>0) THEN
            DO i=1,cellno
            
               cellindx=BufferCellindx(i)
               AtomLocation=LLListRange(cellindx)+(cellindx-1)*MaxAtomPerCell+1 
               DO j=AtomLocation,AtomLocation+BufferCellRange(i)-1
                  AtomUpdated=AtomUpdated+1
                  LLList(j)=NAtoms-NAtomsFilled+AtomUpdated
                  !WRITE(*,*) "buffer index received",LLList(j),irank,slave
               END DO
               LLListRange(cellindx)=LLListRange(cellindx)+BufferCellRange(i)
               IF(AtomUpdated>NAtomsFilled) THEN
                  WRITE(*,*) "There is some error:AtomUpdated went out",AtomUpdated
                  STOP
               END IF
            END DO
            END IF
         END IF
              
      END SUBROUTINE ReceiveBuffer
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      FUNCTION ExpectedDomain(slaveno) !actual domainw/o periodic bdry condition
      !for buffer atoms
         INTEGER :: slaveno ! suppose I manage a way to find the coordinates, I have the cornercell coordinates
         
         INTEGER, DIMENSION(3) :: ExpectedDomain,irank3d
       !  write(*,*)"BoxSize:",DomainIndex,Domains
         irank3d=DomainIndex
         IF(slaveno==1) THEN
            ExpectedDomain=irank3d+(/-1,0,0/)
         ELSEIF(slaveno==2) THEN
            ExpectedDomain=irank3d+(/1,0,0/)
         ELSEIF(slaveno==3) THEN
            ExpectedDomain=irank3d+(/0,1,0/)
         ELSEIF(slaveno==4) THEN
            ExpectedDomain=irank3d+(/0,-1,0/)
         ELSEIF(slaveno==5) THEN
            ExpectedDomain=irank3d+(/0,0,1/)
         ELSEIF(slaveno==6) THEN
            ExpectedDomain=irank3d+(/0,0,-1/)
         END IF      
         
       END FUNCTION ExpectedDomain
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     
   END SUBROUTINE ExchangeAtomInfo 
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Sort(ListOfMigratedAtom,count) !sorts the list of migrated atoms
      IMPLICIT NONE
      INTEGER :: count
      INTEGER, DIMENSION(count) :: ListOfMigratedAtom,indx
      
      !check if the ListOfMigratedAtom is automatically updated when the sorting is completed, i.e,
      !after sorting the ListOfMigratedAtom should have indices in the ascending order and we can ignore
      !the array indx
      IF (count<20) THEN
         WRITE(UNIT=*,FMT='(" Regular sort")')
        !WRITE(*,*) "soirts:",irank
         CALL sort_pick_int(ListOfMigratedAtom,indx)
      ELSEIF (count<1000) THEN
         WRITE(UNIT=*,FMT='(" Heap sort")') 
         CALL sort_heap_int(ListOfMigratedAtom,indx)
      ELSE
         WRITE(UNIT=*,FMT='(" Quick sort")') 
         CALL sort_quick_int(ListOfMigratedAtom,indx)
      END IF
   END SUBROUTINE Sort
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CreateNeighborDomainSend(idomain,jdomain,kdomain,Domains)
      IMPLICIT NONE
      INTEGER :: slave,nx,ny,nz,nx1,ny1,nz1,idomain,jdomain,kdomain,pos
      !INTEGER, DIMENSION(26) :: NeighborDomainSend
      INTEGER, DIMENSION(3) :: Domains      
      !creates the array which contains the list of neighbouring domains.(26 neighbours)
      
      nx=idomain
      ny=jdomain
      nz=kdomain
      nx1=nx-1
      ny1=ny
      nz1=nz
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      NeighborDomainSend(1)=DomainAssignments(pos)
          
      nx1=nx+1
      ny1=ny
      nz1=nz
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      NeighborDomainSend(2)=DomainAssignments(pos)
     
      nx1=nx
      nz1=nz
      ny1=ny+1
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      NeighborDomainSend(3)=DomainAssignments(pos)
      
      nx1=nx
      nz1=nz
      ny1=ny-1
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      NeighborDomainSend(4)=DomainAssignments(pos)
      
      nx1=nx
      nz1=nz+1
      ny1=ny
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      NeighborDomainSend(5)=DomainAssignments(pos)

      nx1=nx
      nz1=nz-1
      ny1=ny
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      NeighborDomainSend(6)=DomainAssignments(pos)
      !back layer => ny+1, upper line=>nz-1, top layer => ny-1, lower line=>nz+1 middle layer=>ny, middle line =>nz
      nx1=nx-1
      ny1=ny-1
      nz1=nz-1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      NeighborDomainSend(7)=DomainAssignments(pos)
       !(in top layer)
      
      nx1=nx+1
      ny1=ny+1
      nz1=nz+1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      NeighborDomainSend(8)=DomainAssignments(pos)
      
      
      nx1=nx+1
      nz1=nz-1
      ny1=ny-1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      NeighborDomainSend(9)=DomainAssignments(pos)
       !(in top layer)
      
      nx1=nx-1
      nz1=nz+1
      ny1=ny+1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      NeighborDomainSend(10)=DomainAssignments(pos)
      
      nx1=nx-1
      nz1=nz+1
      ny1=ny-1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
    !  pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(11)=DomainAssignments(pos)
            
      nx1=nx+1
      nz1=nz-1
      ny1=ny+1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
    !  pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(12)=DomainAssignments(pos)
            
      nx1=nx+1
      ny1=ny-1
      nz1=nz+1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
      !pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(13)=DomainAssignments(pos)
      
      nx1=nx-1
      ny1=ny+1
      nz1=nz-1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
    !  pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(14)=DomainAssignments(pos)
      
      nx1=nx
      nz1=nz-1
      ny1=ny-1
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
     ! pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(15)=DomainAssignments(pos)
      
      nx1=nx
      nz1=nz+1
      ny1=ny+1
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
      !pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(16)=DomainAssignments(pos)
      
      nx1=nx
      nz1=nz+1
      ny1=ny-1
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
      !pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(17)=DomainAssignments(pos)
      
      nx1=nx
      nz1=nz-1
      ny1=ny+1
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
      !pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(18)=DomainAssignments(pos)
      
      nx1=nx-1
      ny1=ny-1
      nz1=nz
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
     ! pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(19)=DomainAssignments(pos)
      
      nx1=nx+1
      ny1=ny+1
      nz1=nz
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      !pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(20)=DomainAssignments(pos)
      
      nx1=nx+1
      nz1=nz
      ny1=ny-1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
     ! pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(21)=DomainAssignments(pos)
      
      nx1=nx-1
      nz1=nz
      ny1=ny+1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      ny1=ny1-FLOOR(REAL(ny1-1)/REAL(Domains(2)))*Domains(2) !periodic BC
      !pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(22)=DomainAssignments(pos)
      
      nx1=nx-1
      nz1=nz-1
      ny1=ny
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
     ! pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(23)=DomainAssignments(pos)
      
      nx1=nx+1
      nz1=nz+1
      ny1=ny
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
      !pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(24)=DomainAssignments(pos)
      
      nx1=nx+1
      ny1=ny
      nz1=nz-1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
     ! pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(25)=DomainAssignments(pos)
      
      nx1=nx-1
      ny1=ny
      nz1=nz+1
      nx1=nx1-FLOOR(REAL(nx1-1)/REAL(Domains(1)))*Domains(1) !periodic BC
      nz1=nz1-FLOOR(REAL(nz1-1)/REAL(Domains(3)))*Domains(3) !periodic BC
     ! pos= ((nx1-1)*Domains(2)+ny1-1)*Domains(3)+nz1
      !pos=(ny1-1)*Domains(1)+nx1+(nz1-1)*Domains(1)*Domains(2)
      pos=(ny1-1)*Domains(3)+nz1+(nx1-1)*Domains(3)*Domains(2)
      !pos=(nx1-1)*Domains(2)+ny1+(nz1-1)*Domains(1)*Domains(2)
      NeighborDomainSend(26)=DomainAssignments(pos) 
   END SUBROUTINE CreateNeighborDomainSend
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !1->left, 2->right, 3-> top, 4-> bottom, 5-> front, 6-> back
   SUBROUTINE SetupArrays(CornerCellForDomain,DomainNCell,BufferCellThickness,NCell)
   !CornerCellIndex gives location of the corner cell in the domain
   !CellSize is the size of each cell in the linked list
   !DomainNCell is the number of cells in the domain
   !BufferCellThickness is the number of cells forming the buffer
   !NCell is the number of cells in linked list
   
   ! PURPOSE: fills in array ListOfCells and ListRange with the following information
   ! 1) there are six sections for each kind of transition
   ! 2) 1st section is for sending to left and receiving from right
   ! 3) 2nd section is for sending to right and receiving from left
   ! 4) 3rd section is for sending bottom and receiving from top
   ! 5) 4rth section is for sending to top and receiving from bottom
   ! 6) 5th section is for sending to bottom and receiving from top
   ! Each section is divided into two parts. Odd indices show the sending part
   ! Even indices show the receiving cell indx. 
   ! when I do CellExchangeList(i) and i is odd, then I am sending info from my domain
   ! about the atoms which may be in the buffer for my neighbor
   ! when I do CellExchangeList(i+1) (where i is from previous example), I am referring to 
   ! the new cellindx in the neighbor where the information is to be sent
      
      IMPLICIT NONE
      INTEGER :: pos,ix,iy,iz
      TYPE(LinkedListContainer), POINTER :: LL
      INTEGER, DIMENSION(:), POINTER :: LLListRange
      INTEGER, DIMENSION(3) :: dnc,bct,c,DomainNCell,NCell,cc,BufferCellThickness,CornerCellForDomain
      LL=>AL%LL
      !LLList=>LL%List
      LLListRange=>LL%ListRange
      !DomainCell is the number of cells in the x, y and z direction
      !from 1 to center -- from left
      !cc(1:3)=CornerCellIndex(1:3) !location of corner cell for ANY domain
      dnc=DomainNCell
      bct=BufferCellThickness
      cc=CornerCellForDomain
      !write(*,*) "hello", MaxCellsPerNeighbor
      pos=0 !slave 1 (left slave starts from position 0)
      !if the command is Send Buffer (1)=> send to your left
      DO iz=1,dnc(3)
         DO iy=1,dnc(2)
            DO ix=1,bct(1) !range for the neighbor
               c=cc-1+(/ix,iy,iz/)
               !c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !WHY??this is the cell for neighbor domain
               pos=pos+1
               IF (pos>2*MaxCellsPerNeighbor-1) THEN
                  WRITE(6,*) "Err>> Increase the size of CLL in MD domain decomposition"
                  STOP
               END IF
               
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !c=MOD(c-1,NCell)+1
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version
               !center domain
               !Recv buffer(2)=> receive from your right cellindx corresponding to location of buffer cell in my domain
               !the region which is sent to my left is being attached on the right buffer region
               pos=pos+1
               c=c+(/dnc(1),0,0/) !shift the cell location -- this corresponds to position in center domain
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version          
            END DO
         END DO
      END DO
      CellExchangeListRange(1)=pos/2 !these many cells are sent to left and received from left
      ! DO iz=1,pos,2 
      !       write(*,*) "buffer atom is in cellindex:",CellExchangeList(iz)
      !end do
      !return
      pos=2*MaxCellsPerNeighbor
      DO iz=1,dnc(3) ! when send buffer(2) sending to right
         DO iy=1,dnc(2)
            DO ix=dnc(1)-bct(1)+1,dnc(1) !range for the neighbor
               c=cc-1+(/ix,iy,iz/) !this is the cell for neighbor domain
               pos=pos+1
               IF (pos>4*MaxCellsPerNeighbor-1) THEN
                  WRITE(6,*) "Err>> Increase the size of CLL in MD domain decomposition"
                  STOP
               END IF
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version
               !center domain
               !when recv buffer(1) receive from left,
               pos=pos+1
               c=c+(/-dnc(1),0,0/) !shift the cell location -- this corresponds to position in center domain
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !CellExchangeList(pos)=((c(1)-1)*NCell(2)+c(2)-1)*NCell(3)+c(3)
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version        
            END DO
         END DO
      END DO
      !return
      CellExchangeListRange(2)=(pos-2*MaxCellsPerNeighbor)/2 
      pos=4*MaxCellsPerNeighbor
      DO iz=1,dnc(3) ! send buffer(3) sending behind
         DO iy=dnc(2)-bct(2)+1,dnc(2)
            DO ix=1-bct(1),bct(1)+dnc(1) !range for the neighbor
               c=cc-1+(/ix,iy,iz/) !this is the cell for neighbor domain
               pos=pos+1
               IF (pos>6*MaxCellsPerNeighbor-1) THEN
                  WRITE(6,*) "Err>> Increase the size of CLL in MD domain decomposition"
                  STOP
               END IF
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !CellExchangeList(pos)=((c(1)-1)*NCell(2)+c(2)-1)*NCell(3)+c(3) !convert c into 1D version
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version 
               !center domain
               !recvbuffer(4) front buffer layer
               pos=pos+1
               c=c+(/0,-dnc(2),0/) !shift the cell location -- this corresponds to position in center domain
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !recv buffer(4) means receive from front these cellindices are the cells where i will paste the received info
               !CellExchangeList(pos)=((c(1)-1)*NCell(2)+c(2)-1)*NCell(3)+c(3)
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version         
            END DO
         END DO
      END DO
      CellExchangeListRange(3)=(pos-4*MaxCellsPerNeighbor)/2
      pos=6*MaxCellsPerNeighbor
      !from 4 to center -- 4 sending front and 3 receiving from back
      !pos=4*MaxCellsPerNeighbor
      DO iz=1,dnc(3) ! send buffer(3) sending front
         DO iy=1,bct(2)
            DO ix=1-bct(1),bct(1)+dnc(1) !range for the neighbor
               c=cc-1+(/ix,iy,iz/) !this is the cell for neighbor domain
               pos=pos+1
               IF (pos>8*MaxCellsPerNeighbor-1) THEN
                  WRITE(6,*) "Err>> Increase the size of CLL in MD domain decomposition"
                  STOP
               END IF
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !CellExchangeList(pos)=((c(1)-1)*NCell(2)+c(2)-1)*NCell(3)+c(3) !convert c into 1D version
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version 
               !center domain
               pos=pos+1 !recvbuffer(3) means receive from back, to form the back layer
               c=c+(/0,dnc(2),0/) !shift the cell location -- this corresponds to position in center domain
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !recv buffer(3) means receive from behind these cellindices are the cells where i will paste the received info
               !CellExchangeList(pos)=((c(1)-1)*NCell(2)+c(2)-1)*NCell(3)+c(3)
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version             
            END DO
         END DO
      END DO
      CellExchangeListRange(4)=(pos-6*MaxCellsPerNeighbor)/2 
     
      pos=8*MaxCellsPerNeighbor
      !for 5 to center -- 5 sending to bottom and 6 receiving from top
      DO iz=dnc(3)-bct(3)+1,dnc(3)
         DO iy=1-bct(2),bct(2)+dnc(2) !sendbuffer(5) send info from foll cells to bottom
            DO ix=1-bct(1),bct(1)+dnc(1) !range for the neighbor
               c=cc-1+(/ix,iy,iz/) !this is the cell for neighbor domain
               pos=pos+1
               IF (pos>10*MaxCellsPerNeighbor-1) THEN
                  WRITE(6,*) "Err>> Increase the size of CLL in MD domain decomposition"
                  STOP
               END IF
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !CellExchangeList(pos)=((c(1)-1)*NCell(2)+c(2)-1)*NCell(3)+c(3) !convert c into 1D version
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version
               !center domain
               pos=pos+1
               c=c+(/0,0,-dnc(3)/) !shift the cell location -- this corresponds to top layer
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !CellExchangeList(pos)=((c(1)-1)*NCell(2)+c(2)-1)*NCell(3)+c(3)
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version             
            END DO
         END DO
      END DO
      CellExchangeListRange(5)=(pos-8*MaxCellsPerNeighbor)/2 
     
      pos=10*MaxCellsPerNeighbor
      !6 send to top and 5 receive from bottom
      DO iz=1,bct(3)
         DO iy=1-bct(2),bct(2)+dnc(2) !sendbuffer(6) send info from foll cells to top
            DO ix=1-bct(1),bct(1)+dnc(1) !range for the neighbor
               c=cc-1+(/ix,iy,iz/) !this is the cell for neighbor domain
               pos=pos+1
               IF (pos>12*MaxCellsPerNeighbor-1) THEN
                  WRITE(6,*) "Err>> Increase the size of CLL in MD domain decomposition"
                  STOP
               END IF
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !CellExchangeList(pos)=((c(1)-1)*NCell(2)+c(2)-1)*NCell(3)+c(3) !convert c into 1D version
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version 
               !center domain
               pos=pos+1 !recv buffer(5) means receive from bottom
               c=c+(/0,0,dnc(3)/) 
               !c(2)=c(2)+(c(3)-1)/NCell(3) !periodic BC
               !c(1)=c(1)+(c(2)-1)/NCell(2)
               c=c-NCell*FLOOR(REAL(c-1)/REAL(NCell)) !wrap the box using periodic BC
               !CellExchangeList(pos)=((c(1)-1)*NCell(2)+c(2)-1)*NCell(3)+c(3)
               !CellExchangeList(pos)=(c(2)-1)*NCell(1)+c(1)+(c(3)-1)*NCell(1)*NCell(2) !convert c into 1D version
               CellExchangeList(pos)=(c(2)-1)*NCell(3)+c(3)+(c(1)-1)*NCell(3)*NCell(2) !convert c into 1D version           
            END DO
         END DO
      END DO
      
      CellExchangeListRange(6)=(pos-10*MaxCellsPerNeighbor)/2 !these many cells are present from neighbor 6 
     
     
     !i have tested cellexchangelist.. work for send buffer and recev buffer and migrated atom and receive from slaves  
      
   END SUBROUTINE SetupArrays
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE sort_pick_int(arr,brr)
      !Trivial sort - use when SIZE(arr)<20
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(INOUT) :: arr
      INTEGER, DIMENSION(:), INTENT(INOUT) :: brr
      INTEGER :: i,j,n,b
      REAL(DP) :: a
      
      n=SIZE(arr)
      
      DO j=2,n
         a=arr(j)
         b=brr(j)
         DO i=j-1,1,-1
            IF (arr(i) <= a) EXIT
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
         END DO
         arr(i+1)=a
         brr(i+1)=b
      END DO
   END SUBROUTINE sort_pick_int
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
   SUBROUTINE sort_heap_int(arr,brr) !sorting..
      !Heap sort - use when SIZE(arr)<1000
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(INOUT) :: arr
      INTEGER, DIMENSION(:), INTENT(INOUT) :: brr
      INTEGER(I4B) :: i,n
      
      n=SIZE(arr)
      DO i=n/2,1,-1
         CALL sift_down(i,n)
      END DO
      DO i=n,2,-1
         CALL swap(arr(1),arr(i))
         CALL swap(brr(1),brr(i))
         CALL sift_down(1,i-1)
      END DO
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE sift_down(l,r)
         INTEGER(I4B), INTENT(IN) :: l,r
         INTEGER(I4B) :: j,jold
         REAL(DP) :: a,b
         
         a=arr(l)
         b=brr(l)
         jold=l
         j=l+l
         DO
            IF (j > r) EXIT
            IF (j < r) THEN
               IF (arr(j) < arr(j+1)) j=j+1
            END IF
            IF (a >= arr(j)) EXIT
            arr(jold)=arr(j)
            brr(jold)=brr(j)
            jold=j
            j=j+j
         END DO
         arr(jold)=a
         brr(jold)=b
      END SUBROUTINE sift_down
   END SUBROUTINE sort_heap_int
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE sort_quick_int(arr,brr)
      !Quick sort - use when SIZE(arr)>1000
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(INOUT) :: arr
      INTEGER, DIMENSION(:), INTENT(INOUT) :: brr
      INTEGER, PARAMETER :: NN=15, NSTACK=50
      REAL(DP) :: a
      INTEGER :: b,n,k,i,j,jstack,l,r
      INTEGER, DIMENSION(NSTACK) :: istack
      
      n=SIZE(arr)
      jstack=0
      l=1
      r=n
      DO
         IF (r-l < NN) THEN
            DO j=l+1,r
               a=arr(j)
               b=brr(j)
               DO i=j-1,l,-1
                  IF (arr(i) <= a) EXIT
                  arr(i+1)=arr(i)
                  brr(i+1)=brr(i)
               END DO
               arr(i+1)=a
               brr(i+1)=b
            END DO
            IF (jstack == 0) RETURN
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
         ELSE
            k=(l+r)/2
            CALL swap(arr(k),arr(l+1))
            CALL swap(brr(k),brr(l+1))
            IF (arr(l)>arr(r)) THEN
               CALL swap(arr(l),arr(r))
               CALL swap(brr(l),brr(r))
            END IF
            IF (arr(l+1)>arr(r)) THEN
               CALL swap(arr(l+1),arr(r))
               CALL swap(brr(l+1),brr(r))
            END IF
            IF (arr(l)>arr(l+1)) THEN
               CALL swap(arr(l),arr(l+1))
               CALL swap(brr(l),brr(l+1))
            END IF
            i=l+1
            j=r
            a=arr(l+1)
            b=brr(l+1)
            DO
               DO
                  i=i+1
                  IF (arr(i) >= a) EXIT
               END DO
               DO
                  j=j-1
                  IF (arr(j) <= a) EXIT
               END DO
               IF (j < i) EXIT
               CALL swap(arr(i),arr(j))
               CALL swap(brr(i),brr(j))
            END DO
            arr(l+1)=arr(j)
            brr(l+1)=brr(j)
            arr(j)=a
            brr(j)=b
            jstack=jstack+2
            IF (jstack > NSTACK) THEN
               WRITE(*,*) "$Err>> Stack size small in sort_quick_int"
               STOP
            END IF
            IF (r-i+1 >= j-l) THEN
               istack(jstack)=r
               istack(jstack-1)=i
               r=j-1
            ELSE
               istack(jstack)=j-1
               istack(jstack-1)=l
               l=i
            END IF
         END IF
      END DO
   END SUBROUTINE sort_quick_int
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
END MODULE mdmpi_domaindecomposition
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
