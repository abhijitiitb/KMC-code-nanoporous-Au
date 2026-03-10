MODULE VARIABLE_TYPE
   IMPLICIT NONE
   INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
   INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
   INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
   INTEGER, PARAMETER :: SP = KIND(1.0)
   INTEGER, PARAMETER :: DP = KIND(1.0D0)
   INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
   INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
   INTEGER, PARAMETER :: LGT = KIND(.true.)
   REAL(DP), PARAMETER :: PI=3.141592653589793238462643383279502884197_dp
   REAL(DP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_dp
   REAL(DP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_dp
   REAL(DP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_dp
   REAL(DP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_dp
   
   !Units used in this work:
   !Mass : amu (atomic mass unit)
   !Energy : eV (electron volt); Velocity=sqrt(eV/amu); Momentum=sqrt(amu.eV)
   !Distance : Ang; Force=eV/Ang
   !Time : second; prefactor=1/second
   !Conversion factors -
   REAL(DP), PARAMETER :: famu=1.660538772e-27_dp !amu to kg
   REAL(DP), PARAMETER :: femu=9.1093826e-31_dp !electron rest mass to kg
   REAL(DP), PARAMETER :: fev=1.60217653e-19_dp !electron charge in C
   REAL(DP), PARAMETER :: fang=1.0e-10_dp !angstroms to m
   REAL(DP), PARAMETER :: fbohr=5.291772108e-11_dp !Bohr radius to m
   REAL(DP), PARAMETER :: fhb=1.05457168e-34 !Reduced Planck's const (h/2pi) in Js
   REAL(DP), PARAMETER :: fhart=4.35974417e-18 !Hartree in J
   REAL(DP), PARAMETER :: Rgas=8.314472_dp !in J/K.mol
   REAL(DP), PARAMETER :: ftime=2.418884326505e-17_dp !in sec = fhb/fhart
   REAL(DP), PARAMETER :: fvel=2.1876912633e6_dp !in m/s = fbohr*fhart/fhb
   REAL(DP), PARAMETER :: ftemp=3.1577464e5 !in K = fhart/kb
   REAL(DP), PARAMETER :: fforce=8.2387225e-8 !in N = fhart/fbohr
   REAL(DP), PARAMETER :: fpressure=2.9421912e13 !in N/m2 = fhart/fbohr^3
   REAL(DP), PARAMETER :: NAvagadro=6.02214179e23_dp !Avagadro number
   REAL(DP), PARAMETER :: kboltzmann=8.617343e-5_dp !boltzmann const in eV/K
   
   INTEGER, PARAMETER :: UnitTmp1=301,UnitTmp2=302,UnitTmp3=303
   INTEGER, PARAMETER :: NPotentialType=10 !see Potential
   INTEGER, PARAMETER :: MaxNAtoms=150000 !used by Optimization codes as size for arrays
   INTEGER, PARAMETER :: MaxNEBAtoms=15000 !used by NEB codes as size for arrays
   INTEGER, PARAMETER :: MaxNFilteredAtoms=20000 !used by FilterAtoms
   INTEGER, PARAMETER :: MaxSlaveNAtoms=5000 !use for parallel codes
   INTEGER, PARAMETER :: MaxEAMSpecies=5,MaxEAMTableSize=4000
   INTEGER, PARAMETER :: MaxLJSpecies=5,MaxLJTableSize=4000
   INTEGER, PARAMETER :: MaxBuckSpecies=5,MaxBuckTableSize=4000
   INTEGER, PARAMETER :: MaxSWSpecies=5,MaxSWTableSize=4000
   INTEGER, PARAMETER :: MaxMEAMSpecies=2,MaxMEAMTableSize=4000
   INTEGER, PARAMETER :: MaxOPLSSpecies=10
   INTEGER, PARAMETER :: MaxTersoffSpecies=5
   INTEGER, PARAMETER :: MaxNTicks=200 !IMPORTANT : keep this more than 10 because tvalidBandScan needs it to be atleast 11
   INTEGER, PARAMETER :: MaxPolymerAtomConnectivity=6

   INTEGER, PARAMETER :: NINT2CHAR=5

   INTEGER, DIMENSION(:), POINTER :: SpeciesDirectory_global
   INTEGER :: NSpecies_global=0 !total number of species
   CHARACTER(len=6) :: TxtHeader 
   REAL(dp) :: Atom_Density(200000)=0._dp
   INTEGER :: LatticeSpeciesArray(200000)=0,LatticeEnvArray(200000)=0
   
   LOGICAL :: IntelligentPrint=.FALSE. !used to print special circumstances
 
   !master controls
   INTEGER :: VLListTypeDefault=1 !1= full list, 2= half list
   REAL(dp), DIMENSION(3) :: ElectricField=0._dp
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE MDContainer   
   !Last checked: July 03, 2009
   !Last modified: 
   !Last modified by: AC

      CHARACTER(len=10) :: Integrator
      
      INTEGER :: HistorySize=0,HistorySizeMax=0
         !# of AL chain of states in History - 
         ! e.g., if HistorySize is 10 then a 
         ! stack of 10 is created, oldest snaps are lost
         ! the current snapshot will correspond to 
         ! (HistorySize+1)th snapshot
         !Note that current snapshot is not connected
         ! to the History
      !INTEGER :: HistoryFreq=0 !how often snaps should be stored in history
      
      LOGICAL :: HistoryOn=.FALSE.
      
      REAL(dp) :: dt=0._dp !MD time increment
      REAL(dp) :: Temperature=0._dp !system temperature
      REAL(dp) :: LangevinCoeff=0._dp !Langevin dynamics
      REAL(dp) :: TemperatureRamp=0 !in K/md step with which the temperature is modified
        !for instance if TemperatureRamp=1.e-3 K/MD step then the temperature will increase by 1 K after 1000 MD steps
        !if TemperatureRamp=-1.e-3 K/MD step then the temperature will decrease by 1 K
      
      TYPE(SystemContainer), POINTER :: AL=>NULL() !more detailed system
      TYPE(SystemContainer), POINTER :: History=>NULL() !chain of states
         !History is begining of the ChOS, while AL points to end
         
   END TYPE MDContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE MDProcessContainer !contains the process found using MD
      REAL(dp) :: prefactor=0._dp,Eact=0._dp
      INTEGER :: ntransitionstates=0
   END TYPE MDProcessContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE SystemContainer   
   !Last checked: July 03, 2009
   !Last modified: 
   !Last modified by: AC
      
      TYPE(SystemContainer), POINTER :: PrevNeigh=>NULL() !used for creating a chain of AL
      TYPE(SystemContainer), POINTER :: NextNeigh=>NULL() !used for creating a chain of AL
      TYPE(LinkedListContainer), POINTER :: LL=>NULL() !list of linked list neighbours
      TYPE(VerletListContainer), POINTER :: VL=>NULL() !list of verlet list neighbours
      
      !AtomPolarize is the displacement of the shell charge from the core
      
      !Scalar
      REAL(dp) :: Time=0._dp !elapsed time
      REAL(dp) :: PotentialEnergy=0._dp
      REAL(dp) :: KineticEnergy=0._dp
      INTEGER :: NAtoms=0
      !INTEGER :: NSpecies=0
      INTEGER, DIMENSION(3) :: NMove=0 !# dimensions moving along x,y,z
      CHARACTER(len=30) :: BoundaryCondX
      CHARACTER(len=30) :: BoundaryCondY
      CHARACTER(len=30) :: BoundaryCondZ
      
      !Array-other
      !INTEGER, DIMENSION(:), POINTER :: SpeciesDirectory=>NULL() 
         !contains ref. to atomic number
         !  for e.g., SpeciesDirectory(1)=12 indicates that
         !  species type 1 is C atom
      REAL(dp), DIMENSION(6) :: BoxSize=0._dp
         !includes [upper_limit lower_limit]
      
      !Arrays of Size(NAtoms)
      REAL(dp), DIMENSION(:), POINTER :: AtomCharge=>NULL()
      INTEGER, DIMENSION(:), POINTER :: AtomSpecies=>NULL()
      
      !Arrays of Size(3*NAtoms)
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: AtomVelocity=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: AtomDriftVelocity=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: AtomMass=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: AtomForce=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: AtomShiftCoord=>NULL()
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: ChOSTangent=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: ChOSForce=>NULL()
      
      !Arrays of Size(3*NAtoms*(1+3*NAtoms)/2)
      !REAL(dp), DIMENSION(:), POINTER :: AtomEnergyHessian=>NULL()
      
      !InteracPotentials for this system
      TYPE(InteracPotential), POINTER :: Potential=>NULL()
      TYPE(PolymerAtomInfo), POINTER :: PolymerConnectivity=>NULL()
      
   END TYPE SystemContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE ChOSContainer

      TYPE(SystemContainer), POINTER :: AL=>NULL()
      REAL(dp) :: SpringConst=0._dp,StringWeight=0._dp
      !REAL(dp), DIMENSION(:), POINTER :: g,gold,fold,spline !g is search direction
      INTEGER :: nimages=0

   END TYPE ChOSContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE ListOfScalars
      
      TYPE(ListOfScalars), POINTER :: PrevNeigh=>NULL()
      TYPE(ListOfScalars), POINTER :: NextNeigh=>NULL()
      REAL(dp) :: value=0._dp
      
   END TYPE ListOfScalars
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE LinkedListContainer

      !Each system is divided into NumberCell cells
      !alternatively in xyz NumberCell3D(3) cells
      !each cell can have MaxAtom atoms in them
      !these atom indicess are stored in List(r1:r2) for cell i
      !here r1=MaxAtom*(i-1)+1, r2=MaxAtom*(i-1)+ListRange(i)
      !in other words ListRange(i) gives # of atoms in cell i

      LOGICAL :: UnInitialized=.TRUE.
      
      TYPE(LinkedListContainer), POINTER :: PrevNeigh=>NULL()
      TYPE(LinkedListContainer), POINTER :: NextNeigh=>NULL()
      
      REAL(dp), DIMENSION(3) :: CellSize !gives an idea of the potential cutoff
      
      INTEGER :: MaxAtom=0
      INTEGER :: NAtoms
      INTEGER :: MaxAtomPerCell=0
      INTEGER :: NumberCell3D(3)=0
      INTEGER :: NumberCell=0
      
      INTEGER, DIMENSION(:), POINTER :: ListRange=>NULL()
      INTEGER, DIMENSION(:), POINTER :: List=>NULL()
      INTEGER, DIMENSION(:), POINTER :: AtomCellIndx=>NULL() !this is not serving much purpose
      
   END TYPE LinkedListContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE VerletListContainer

      !Each atom can have a record of size MaxAtom of neighbors
      !The indices of these neighbors are stored in List
      !For atom i the neighbors are stored in List(r1:r2)
      !here r1=MaxAtom*(i-1)+1, and r2=MaxAtom*(i-1)+ListRange(i)
      !in other words ListRange(i) gives # of neighbors of atom i

      LOGICAL :: UnInitialized=.FALSE.
   
      REAL(dp), DIMENSION(:), POINTER :: drmag=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: OrigCoord=>NULL()
      !REAL(dp), DIMENSION(3*MaxNAtoms) :: OrigCoord
      REAL(dp), DIMENSION(:), POINTER :: dr=>NULL()
      
      REAL(dp) :: CutOff=0._dp
      REAL(dp) :: Buffer=0._dp

      INTEGER :: MaxAtom=0
      INTEGER :: NAtoms=0
      INTEGER :: MaxAtomPerAtom=0
      INTEGER :: ListType=0 !0=not known, 1=full list, 2=half list
      INTEGER, DIMENSION(:), POINTER :: ListRange=>NULL()
      INTEGER, DIMENSION(:), POINTER :: List=>NULL()
      LOGICAL, DIMENSION(:), POINTER :: ListDomainAtom=>NULL()

   END TYPE VerletListContainer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE PolymerAtomInfo
      INTEGER, DIMENSION(:), POINTER ::Connectivity=>NULL(),NumberBonds=>NULL() !provides
      !connectivity information between atoms in the polymer chain
      INTEGER, DIMENSION(:), POINTER :: Hybridization=>NULL() !following
      !notation is used: 1-sp3, 2-sp2, 3-sp, 0-unknown
   END TYPE PolymerAtomInfo
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !                           Interaction Potential description
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !EAM-all types of SMA-TB
   TYPE EAMPotential
      !EE=EmbeddingEnergy PP=PairPotential DT=DensityTerm
      !for upto n=5 species there will be n(n+1)/2=15 required for PP and n=5 for EE and DT
      CHARACTER(len=10) :: TableFormat="xxxxxxxxxx" !e.g., VOTER, SETFL, FUNCFL -- only one of these type is allowed
      CHARACTER(len=100), DIMENSION((MaxEAMSpecies*(MaxEAMSpecies+1))/2) :: TableLocationPP
      CHARACTER(len=100), DIMENSION(MaxEAMSpecies) :: TableLocationDT,TableLocationEE
      INTEGER, DIMENSION((MaxEAMSpecies*(MaxEAMSpecies+1))/2) :: SizePP=0
      INTEGER, DIMENSION(MaxEAMSpecies) :: SizeDT=0,SizeEE=0  !size of table, i.e., number of points
      REAL(dp), DIMENSION(MaxEAMSpecies*(MaxEAMSpecies+1)) :: RangePP=0._dp
      REAL(dp), DIMENSION(2*MaxEAMSpecies) :: RangeDT=0._dp,RangeEE=0._dp !start and end cutoff
      REAL(dp), DIMENSION((MaxEAMSpecies*(MaxEAMSpecies+1))/2) :: IncrPPinv=0._dp
      REAL(dp), DIMENSION(MaxEAMSpecies) :: IncrEEinv(5)=0._dp,IncrDTinv(5)=0._dp
      REAL(dp), DIMENSION(MaxEAMTableSize*(MaxEAMSpecies*(MaxEAMSpecies+1))/2) :: PP=0._dp
      REAL(dp), DIMENSION(MaxEAMSpecies*MaxEAMTableSize) :: EE=0._dp,DT=0._dp !each pair can store upto MaxEAMTableSize pts
      REAL(dp), DIMENSION((MaxEAMSpecies*(MaxEAMSpecies+1))/2) :: ConvFacPP=1._dp
      REAL(dp), DIMENSION(MaxEAMSpecies) :: ConvFacDT=1._dp,ConvFacEE=1._dp
      LOGICAL, DIMENSION((MaxEAMSpecies*(MaxEAMSpecies+1))/2) :: EnabledPP=.FALSE.
      LOGICAL, DIMENSION(MaxEAMSpecies) :: EnabledDT=.FALSE.,EnabledEE=.FALSE.
      
   END TYPE EAMPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !Lennard Jones
   TYPE LJPotential
      !PP=PairPotential
      !for upto n=5 species there will be n(n+1)/2=15 required for PP
      CHARACTER(len=100), DIMENSION((MaxLJSpecies*(MaxLJSpecies+1))/2):: TableLocationPP
      INTEGER, DIMENSION((MaxLJSpecies*(MaxLJSpecies+1))/2) :: SizePP=0 !size of table, i.e., number of points
      REAL(dp), DIMENSION(MaxLJSpecies*(MaxLJSpecies+1)) :: RangePP=0._dp !start and end cutoff
      REAL(dp), DIMENSION((MaxLJSpecies*(MaxLJSpecies+1))/2) :: IncrPPinv=0._dp
      REAL(dp), DIMENSION(MaxLJTableSize*(MaxLJSpecies*(MaxLJSpecies+1))/2) :: PP=0._dp !each pair can store upto MaxLJTableSize pts
      REAL(dp), DIMENSION((MaxLJSpecies*(MaxLJSpecies+1))/2) :: ConvFacPP=1._dp
      LOGICAL, DIMENSION((MaxLJSpecies*(MaxLJSpecies+1))/2) :: EnabledPP=.FALSE.
   END TYPE LJPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !Coulomb
   TYPE EwaldParameters
      !LOGICAL :: IsUnInitialized=.TRUE. !Ewald parameters need to be read only once
      INTEGER :: NAtoms=0 !number of atoms for which the Ewald is initialized
      INTEGER :: ImCutOff=0 !cutoff for imaginary part
      REAL(dp) :: ReCutOff=0._dp !cutoff for real part -- this is akin to PotentialCutOff in LJ and EAM
      REAL(dp) :: ALPHA=0._dp
      REAL(dp) :: TimeRatio=0._dp !ratio CPU time(real time)/CPU time (reciprocal space)
      REAL(dp) :: Precision=0._dp !accuracy we desire from the calculation
      INTEGER :: TOTk=0 !TOTr=0, !this is the information size used for arrays KVEC and RVEC, respectively
      INTEGER, DIMENSION(:), POINTER :: KVector=>NULL() !RVEC=>NULL(),
      !COMPLEX(dp), DIMENSION(:,:), POINTER :: EIKX=>NULL()
      !COMPLEX(dp), DIMENSION(:,:), POINTER :: EIKY=>NULL()
      !COMPLEX(dp), DIMENSION(:,:), POINTER :: EIKZ=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: KVEC=>NULL()
   END TYPE EwaldParameters
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE WolfParameters
      REAL(dp) :: ALPHA=0._dp
      REAL(dp) :: ReCutOff=0._dp
   END TYPE WolfParameters
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE PPPMParameters
      INTEGER :: PPM
   END TYPE PPPMParameters
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE FMMParameters
      INTEGER :: FMM
   END TYPE FMMParameters
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE CoulombPotential
      TYPE(EwaldParameters), POINTER :: Ewald=>NULL()
      TYPE(WolfParameters), POINTER :: Wolf=>NULL()
      TYPE(PPPMParameters), POINTER :: PPPM=>NULL()
      TYPE(FMMParameters), POINTER :: FMM=>NULL() !Fast multipole method
   END TYPE CoulombPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !Buckingham - potential type 4
   TYPE BuckinghamPotential
      !PP=PairPotential
      !for upto n=5 species there will be n(n+1)/2=15 required for PP
      CHARACTER(len=100), DIMENSION((MaxBuckSpecies*(MaxBuckSpecies+1))/2):: TableLocationPP
      INTEGER, DIMENSION((MaxBuckSpecies*(MaxBuckSpecies+1))/2) :: SizePP=0 !size of table, i.e., number of points
      REAL(dp), DIMENSION(MaxBuckSpecies*(MaxBuckSpecies+1)) :: RangePP=0._dp !start and end cutoff
      REAL(dp), DIMENSION((MaxBuckSpecies*(MaxBuckSpecies+1))/2) :: IncrPPinv=0._dp
      REAL(dp), DIMENSION(MaxBuckTableSize*(MaxBuckSpecies*(MaxBuckSpecies+1))/2) :: PP=0._dp !each pair can store upto MaxBuckTableSize pts
      REAL(dp), DIMENSION((MaxBuckSpecies*(MaxBuckSpecies+1))/2) :: ConvFacPP=1._dp
      LOGICAL, DIMENSION((MaxBuckSpecies*(MaxBuckSpecies+1))/2) :: EnabledPP=.FALSE.
   END TYPE BuckinghamPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !Stillinger-Weber
   TYPE StillingerWeberPotential
      INTEGER :: format=1 !(1)analytical or (2)table
      !analytical evaluation
      REAL(dp), DIMENSION(MaxSWSpecies*(MaxSWSpecies+1)/2) :: A,B,g,q,epsln,sigma,ac,lmd,gamma,cutoff
      
      !table-based evaluation
      !PP=PairPotential
      !for upto n=5 species there will be n(n+1)/2=15 required for PP
      CHARACTER(len=100), DIMENSION((MaxSWSpecies*(MaxSWSpecies+1))/2):: TableLocationPP
      INTEGER, DIMENSION((MaxSWSpecies*(MaxSWSpecies+1))/2) :: SizePP=0,SizeF3=0 !size of table, i.e., number of points
      REAL(dp), DIMENSION(MaxSWSpecies*(MaxSWSpecies+1)) :: RangePP=0._dp,RangeF3=0._dp !start and end cutoff
      REAL(dp), DIMENSION((MaxSWSpecies*(MaxSWSpecies+1))/2) :: IncrPPinv=0._dp,IncrF3inv=0
      REAL(dp), DIMENSION(MaxSWTableSize*(MaxSWSpecies*(MaxSWSpecies+1))/2) :: PP=0._dp,F3=0._dp !each pair can store upto MaxLJTableSize pts
      REAL(dp), DIMENSION((MaxSWSpecies*(MaxSWSpecies+1))/2) :: epssqrt=0._dp,lam=0._dp
      REAL(dp), DIMENSION((MaxSWSpecies*(MaxSWSpecies+1))/2) :: ConvFacPP=1._dp
      LOGICAL, DIMENSION((MaxSWSpecies*(MaxSWSpecies+1))/2) :: EnabledPP=.FALSE.
   END TYPE StillingerWeberPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !Tersoff - no table used
   TYPE TersoffPotential !based on Tersoff, PRB, 37, 12, 6991 (1998) 
      REAL(dp), DIMENSION(MaxTersoffSpecies*(MaxTersoffSpecies+1)/2) :: &
         A,B,lam1,lam2,lam3,beta,R,RD,c,d,h,n
      REAL(dp), DIMENSION(MaxTersoffSpecies*(MaxTersoffSpecies+1)/2) :: Rangefc=0._dp !maximum cutoff for potential
   END TYPE TersoffPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !FENE
   TYPE FENEPotential
      INTEGER :: Fill
   END TYPE FENEPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !AIREBO
   TYPE AIREBOPotential
      INTEGER :: Fill
   END TYPE AIREBOPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !MEAM
   TYPE MEAMPotential
      !EE=embedding energy, DT=density term, PP=pair potential
      CHARACTER(len=10) :: TableFormat="xxxxxxxxxx"
      !CHARACTER(len=100), DIMENSION((MaxMEAMSpecies*(MaxMEAMSpecies+1))/2) :: TableLocationPP
      !CHARACTER(len=100), DIMENSION(MaxMEAMSpecies) :: TableLocationDT,TableLocationEE
      INTEGER, DIMENSION((MaxMEAMSpecies*(MaxMEAMSpecies+1))/2) :: SizePP=0
      INTEGER, DIMENSION(MaxMEAMSpecies) :: SizeDT=0,SizeEE=0  !size of table, i.e., number of points
      REAL(dp), DIMENSION(MaxMEAMSpecies*(MaxMEAMSpecies+1)) :: RangePP=0._dp
      !REAL(dp), DIMENSION(2*MaxMEAMSpecies) :: RangeDT=0._dp,RangeEE=0._dp !start and end cutoff
      REAL(dp), DIMENSION((MaxMEAMSpecies*(MaxMEAMSpecies+1))/2) :: IncrPPinv=0._dp
      !REAL(dp), DIMENSION(MaxMEAMSpecies) :: IncrEEinv(5)=0._dp,IncrDTinv(5)=0._dp
      REAL(dp), DIMENSION(MaxMEAMTableSize*(MaxMEAMSpecies*(MaxMEAMSpecies+1))/2) :: PP=0._dp
      !REAL(dp), DIMENSION(MaxMEAMSpecies*MaxMEAMTableSize) :: EE=0._dp,DT=0._dp !each pair can store upto MaxEAMTableSize pts
      REAL(dp), DIMENSION((MaxMEAMSpecies*(MaxMEAMSpecies+1))/2) :: ConvFacPP=1._dp
      !REAL(dp), DIMENSION(MaxMEAMSpecies) :: ConvFacDT=1._dp,ConvFacEE=1._dp
      LOGICAL, DIMENSION((MaxMEAMSpecies*(MaxMEAMSpecies+1))/2) :: EnabledPP=.FALSE.
      !LOGICAL, DIMENSION(MaxMEAMSpecies) :: EnabledDT=.FALSE.,EnabledEE=.FALSE.
      REAL(dp), DIMENSION(MaxMEAMSpecies) :: A,rho0,beta,t
      REAL(dp), DIMENSION((MaxMEAMSpecies*(MaxMEAMSpecies+1))/2) :: Ec,Re,d,alfa !d,alfa=sqrt(9*B*omega/Ec) used for pair potential
      REAL(dp), DIMENSION(MaxMEAMSpecies*(MaxMEAMSpecies*(MaxMEAMSpecies+1))/2) :: Cmin,Cmax !see ReadPotentialInfo for more info

   END TYPE MEAMPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !Tight binding
   TYPE TBPotential
      INTEGER :: Fill
   END TYPE TBPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE QMPotential
      INTEGER :: Fill
   END TYPE QMPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE OPLSPotential 
      REAL(dp), DIMENSION((MaxOPLSSpecies*(MaxOPLSSpecies+1))/2) :: kSpring,Req
      !REAL(dp), DIMENSION( :: angle,dihedral term
   END TYPE OPLSPotential
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE InteracPotential

      !TYPE (InteracPotential), POINTER :: NextNeigh=>NULL(), PrevNeigh=>NULL()
         !when two species interact with more than one potential form -- this is disabled as of now so that within this InteracPotential the user can find all possible potential options. For e.g., if the atom has LJ, Buckingham and Coulomb terms then all three can be enabled below
      
      LOGICAL :: UnInitialized=.TRUE.

      REAL(dp) :: VLBuffer=0._dp

      REAL(dp), DIMENSION(NPotentialType) :: MaxPotentialCutoff=0._dp !max cutoff for pot # i

      TYPE(LJPotential), POINTER :: LJ=>NULL() !pot # 1
      TYPE(EAMPotential), POINTER :: EAM=>NULL() !pot # 2
      TYPE(CoulombPotential), POINTER :: Coulomb=>NULL() !pot # 3, e.g. CouEw,CouWo,CouPP
      TYPE(BuckinghamPotential), POINTER :: Buckingham=>NULL() !pot # 4
      TYPE(StillingerWeberPotential), POINTER :: SW=>NULL() !pot # 5
      TYPE(TersoffPotential), POINTER :: Tersoff=>NULL() !pot # 6
      TYPE(FENEPotential), POINTER :: FENE=>NULL() !pot # 7
      TYPE(AIREBOPotential), POINTER :: AIREBO=>NULL() !pot # 8
      TYPE(MEAMPotential), POINTER :: MEAM=>NULL() !pot # 9
      TYPE(TBPotential), POINTER :: TB=>NULL() !pot # 10
      !TYPE(LenoskyPotential), POINTER :: Lenosky=>NULL() !pot # 11
      !TYPE(ReaxFF), POINTER :: ReaxFF=>NULL() !pot # 12
      !TYPE(EDIPPotential), POINTER :: EDIP=>NULL() !pot # 13
      TYPE(QMPotential), POINTER :: DFT=>NULL() !pot # 14
      TYPE(OPLSPotential), POINTER :: OPLS=>NULL() !pot #15   

   END TYPE InteracPotential

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !                          Species information
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   TYPE SpeciesDetails

      INTEGER :: NSpecies
      
      CHARACTER(len=3), DIMENSION(:), POINTER :: AtomicSymbol=>NULL()
      CHARACTER(len=12), DIMENSION(:), POINTER :: AtomicGroup=>NULL()
      CHARACTER(len=12), DIMENSION(:), POINTER :: AtomicOxidState=>NULL()
      CHARACTER(len=20), DIMENSION(:), POINTER :: AtomicName=>NULL()
      
      INTEGER, DIMENSION(:), POINTER :: AtomicNumber=>NULL()
      
      REAL(dp), DIMENSION(:), POINTER :: AtomicMass=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: AtomicMeltPt=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: AtomicBoilPt=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: AtomElectroneg_Pauling=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: Atomic1stIonization=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: AtomicDensity=>NULL()
      REAL(dp), DIMENSION(:), POINTER :: AtomicRadius=>NULL()
      
   END TYPE SpeciesDetails
   
   !xxxxxxxxxxxxxxxxxxxxxxxxx
   
   TYPE(SpeciesDetails), POINTER :: SpeciesList=>NULL()
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   TYPE RDFHistogram !for real valued variables
      REAL(dp), DIMENSION(:), POINTER :: Hist=>NULL()
      REAL(dp) :: xstart,xend
      INTEGER :: nbin
   END TYPE RDFHistogram
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE VARIABLE_TYPE

