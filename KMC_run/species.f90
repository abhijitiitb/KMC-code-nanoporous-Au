MODULE Species

   !contains db of chemical species
   
   USE VARIABLE_TYPE
   USE db_manipulate
   IMPLICIT NONE
   
   CONTAINS
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE SetUpSpecies()

      IMPLICIT NONE
      INTEGER :: NSpecies
      
      WRITE(*,*) "SpeciesData ... setting up species list"
      NSpecies=117
      ALLOCATE(SpeciesList)
      SpeciesList%NSpecies=NSpecies
      CALL MakeSize(SpeciesList%AtomicSymbol,NSpecies)
      CALL MakeSize(SpeciesList%AtomicGroup,NSpecies)
      CALL MakeSize(SpeciesList%AtomicOxidState,NSpecies)
      CALL MakeSize(SpeciesList%AtomicName,NSpecies)
      CALL MakeSize(SpeciesList%AtomicNumber,NSpecies)
      CALL MakeSize(SpeciesList%AtomicMass,NSpecies)
      CALL MakeSize(SpeciesList%AtomicMeltPt,NSpecies)
      CALL MakeSize(SpeciesList%AtomicBoilPt,NSpecies)
      CALL MakeSize(SpeciesList%AtomElectroneg_Pauling,NSpecies)
      CALL MakeSize(SpeciesList%Atomic1stIonization,NSpecies)
      CALL MakeSize(SpeciesList%AtomicDensity,NSpecies)
      CALL MakeSize(SpeciesList%AtomicRadius,NSpecies)

      SpeciesList%AtomicNumber=(/1,2,3,4,5,6,7,8,9,10,11,&
         12,13,14,15,16,17,18,19,20,21,22,23, &
         24,25,26,27,28,29,30,31,32,33,34,35,36, &
         37,38,39,40,41,42,43,44,45,46,47,48, &
         49,50,51,52,53,54,55,56,57,58,59,60,61, &
         62,63,64,65,66,67,68,69,70,71,72,73,74, &
         75,76,77,78,79,80,81,82,83,84,85,86, &
         87,88,89,90,91,92,93,94,95,96,97,98,99, &
         100,101,102,103,104,105,106,107,108,109,&
         110,111,112,113,114,115,116,118/);

      SpeciesList%AtomicSymbol=(/"H  ","He ","Li ","Be ","B  ","C  ","N  ","O  ", &
         "F  ","Ne ","Na ","Mg ","Al ","Si ","P  ","S  ","Cl ","Ar ", &
         "K  ","Ca ","Sc ","Ti ","V  ","Cr ","Mn ","Fe ","Co ","Ni ", &
         "Cu ","Zn ","Ga ","Ge ","As ","Se ","Br ","Kr ","Rb ","Sr ", &
         "Y  ","Zr ","Nb ","Mo ","Tc ","Ru ","Rh ","Pd ","Ag ","Cd ", &
         "In ","Sn ","Sb ","Te ","I  ","Xe ","Cs ","Ba ","La ","Ce ", &
         "Pr ","Nd ","Pm ","Sm ","Eu ","Gd ","Tb ","Dy ","Ho ","Er ", &
         "Tm ","Yb ","Lu ","Hf ","Ta ","W  ","Re ","Os ","Ir ","Pt ", &
         "Au ","Hg ","Tl ","Pb ","Bi ","Po ","At ","Rn ","Fr ","Ra ", &
         "Ac ","Th ","Pa ","U  ","Np ","Pu ","Am ","Cm ","Bk ","Cf ", &
         "Es ","Fm ","Md ","No ","Lr ","Rf ","Db ","Sg ","Bh ","Hs ", &
         "Mt ","Ds ","Rg ","Uub","Uut","Uuq","Uup","Uuh","Uuo"/)

      SpeciesList%AtomicMass=(/1.008_dp,4.003_dp,6.941_dp,9.012_dp,&
         10.81_dp,12.01_dp,14.01_dp,16.00_dp,  &
         19.00_dp,20.18_dp,22.99_dp,24.31_dp, & 
         26.98_dp,28.09_dp,30.97_dp,32.07_dp,35.45_dp,39.95_dp, &
         39.10_dp,40.08_dp,44.96_dp,47.87_dp,50.94_dp, &
         52.00_dp,54.94_dp,55.85_dp,58.93_dp,58.69_dp,&
         63.55_dp,65.41_dp,69.72_dp,72.64_dp,74.92_dp, &
         78.96_dp,79.90_dp,83.80_dp,85.47_dp,87.62_dp,&
         88.91_dp,91.22_dp,92.91_dp,95.94_dp,98._dp, &
         101.1_dp,102.9_dp,106.4_dp,107.9_dp,112.4_dp,   &
         114.8_dp,118.7_dp,121.8_dp,127.6_dp,126.9_dp, &
         131.3_dp,132.9_dp,137.3_dp,138.9_dp,140.1_dp,&
         140.9_dp,144.2_dp,145._dp,150.4_dp,152.0_dp, &
         157.3_dp,158.9_dp,162.5_dp,164.9_dp,167.3_dp,  &
         168.9_dp,173.0_dp,175.0_dp,178.5_dp,180.9_dp, &
         183.8_dp,186.2_dp,190.2_dp,192.2_dp,195.1_dp,&
         197.0_dp,200.6_dp,204.4_dp,207.2_dp,209.0_dp, &
         209._dp,210._dp,222._dp,223._dp,226._dp,227._dp,232.0_dp,&
         231.0_dp,238.0_dp,237._dp,239._dp,243._dp,247._dp,&
         247._dp,251._dp,252._dp,257._dp,258._dp,259._dp,262._dp,261._dp,&
         262._dp,266._dp,264._dp,277._dp,268._dp,281._dp,&
         272._dp,285._dp,284._dp,289._dp,288._dp,292._dp,294._dp/) !in amu

      SpeciesList%AtomicMeltPt=(/-259.1_dp,-272.2_dp,180.5_dp,1278._dp,2079._dp, &
        3367._dp,-209.9_dp,-218.4_dp,-219.8_dp,-248._dp,97.8_dp,649._dp,660._dp, &
        1410._dp,44.1_dp,112.8_dp,-101._dp,-189.2_dp,63.25_dp,839._dp,1541._dp,  &
        1660._dp,1890._dp,1857._dp,1244._dp,1535._dp,1495._dp,1453._dp,1083._dp,     &
        419.6_dp,29.8_dp,947.4_dp,817._dp,217._dp,-7.2_dp,-157._dp,38.9_dp,     &
        769._dp,1523._dp,1852._dp,2468._dp,2617._dp,2172._dp,2310._dp,1966._dp,      &
        1554._dp,962._dp,320.9_dp,156.6_dp,232._dp,631._dp,449.5_dp,113.5_dp,    &
        -111.8_dp,28.4_dp,725._dp,920._dp,798._dp,931._dp,1016._dp,1042._dp,1074._dp,&
        822._dp,1313._dp,1365._dp,1412._dp,1474._dp,1529._dp,1545._dp,819._dp,1663._dp,&
        2227._dp,2996._dp,3410._dp,3180._dp,3045._dp,2410._dp,1772._dp,1064._dp,     &
        -38.9_dp,303._dp,327.5_dp,271._dp,254._dp,302._dp,-71._dp,27._dp,700._dp,    &
        1050._dp,1750._dp,1570._dp,1132._dp,640._dp,641._dp,994._dp,1340._dp,986._dp,  &
        -999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,     &
        -999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,     &
        -999._dp,-999._dp,-999._dp,-999._dp/) !CHECK THIS DATA MAY BE INCRRCT 

      SpeciesList%AtomicBoilPt=(/-252.9_dp,-268.9_dp,1342._dp,2970._dp,2550._dp,4827._dp, &
        -195.8_dp,-183._dp,-188.1_dp,-248.7_dp,883._dp,1090._dp, &
        2467._dp,2355._dp,280._dp,  &
        444.7_dp,-34.6_dp,-185.7_dp,760._dp,1484._dp,&
        2832._dp,3287._dp,3380._dp,2672._dp, &
        1962._dp,2750._dp,2870._dp,2730._dp,2567._dp,906._dp,&
        2403._dp,2830._dp,617._dp,685._dp,&
        58.8_dp,-152._dp,686._dp,1384._dp,3337._dp,4377._dp,&
        4742._dp,4612._dp,4877._dp,3900._dp, &
        3727._dp,3140._dp,2212._dp,765._dp,2080._dp,2270._dp,&
        1950._dp,989.8_dp,184._dp,-107.1_dp,&
        669._dp,1640._dp,3454._dp,3257._dp,3017._dp,3127._dp,&
        3000._dp,1794._dp,1529._dp,3273._dp,&
        3230._dp,2567._dp,2700._dp,2868._dp,1950._dp,1196._dp,&
        3402._dp,4600._dp,5425._dp,5660._dp,&
        5600._dp,5030._dp,4130._dp,3827._dp,3080._dp,357._dp,&
        1457._dp,1740._dp,1560._dp,962._dp, &
        337._dp,-61.8_dp,677._dp,1140._dp,3200._dp,&
        4790._dp,4000._dp,3818._dp,3900._dp,3232._dp,&
        2607._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,&
        -999._dp,-999._dp,-999._dp,-999._dp,&
        -999._dp,-999._dp,-999._dp,-999._dp,-999._dp,&
        -999._dp,-999._dp,-999._dp,-999._dp,-999._dp,&
        -999._dp,-999._dp,-999._dp/) !in deg celcius !CHECK THIS DATA MAY BE INCRRCT 


!one extra species is present
!       SpeciesList%AtomicGroup=(/"1","18","1","2","13","14","15","16","17","18", &
!          "1","2","13","14","15","16","17","18","1","2","3","4","5",  &
!          "6","7","8","9","10","11","12","13","14","15","16","17",    &
!          "18","1","2","3","4","5","6","7","8","9","10","11","12",    &
!          "13","14","15","16","17","18","1","2","lanthanides","lanthanides",&
!          "lanthanides","lanthanides","lanthanides","lanthanides","lanthanides",&
!          "lanthanides","lanthanides","lanthanides","lanthanides","lanthanides",&
!          "lanthanides","lanthanides","3","4","5","6","7","8","9","10","11","12",&
!          "13","14","15","16","17","18","1","2","actinides","actinides","actinides",&
!          "actinides","actinides","actinides","actinides","actinides","actinides",&
!          "actinides","actinides","actinides","actinides","actinides","3","4","5",&
!          "6","7","8","9","10","11","12","13","14","15","16","17","18"/)

      SpeciesList%AtomicName=(/ &
         "hydrogen            ", &
         "helium              ", &
         "lithium             ", &
         "beryllium           ", &
         "boron               ", &
         "carbon              ", &
         "nitrogen            ", &
         "oxygen              ", &
         "fluorine            ", &
         "neon                ", &
         "sodium              ", &
         "magnesium           ", &
         "aluminum            ", &
         "silicon             ", &
         "phosphorus          ", &
         "sulfur              ", &
         "chlorine            ", &
         "argon               ", &
         "potassium           ", &
         "calcium             ", &
         "scandium            ", &
         "titanium            ", &
         "vanadium            ", &
         "chromium            ", &
         "manganese           ", &
         "iron                ", &
         "cobalt              ", &
         "nickel              ", &
         "copper              ", &
         "zinc                ", &
         "gallium             ", &
         "germanium           ", &
         "arsenic             ", &
         "selenium            ", &
         "bromine             ", &
         "krypton             ", &
         "rubidium            ", &
         "strontium           ", &
         "yttrium             ", &
         "zirconium           ", &
         "niobium             ", &
         "molybdenum          ", &
         "technetium          ", &
         "ruthenium           ", &
         "rhodium             ", &
         "palladium           ", &
         "silver              ", &
         "cadmium             ", &
         "indium              ", &
         "tin                 ", &
         "antimony            ", &
         "tellurium           ", &
         "iodine              ", &
         "xenon               ", &
         "cesium              ", &
         "barium              ", &
         "lanthanum           ", &
         "cerium              ", &
         "praseodymium        ", &
         "neodymium           ", &
         "promethium          ", &
         "samarium            ", &
         "europium            ", &
         "gadolinium          ", &
         "terbium             ", &
         "dysprosium          ", &
         "holmium             ", &
         "erbium              ", &
         "thulium             ", &
         "ytterbium           ", &
         "lutetium            ", &
         "hafnium             ", &
         "tantalum            ", &
         "tungsten            ", &
         "rhenium             ", &
         "osmium              ", &
         "iridium             ", &
         "platinum            ", &
         "gold                ", &
         "mercury             ", &
         "thallium            ", &
         "lead                ", &
         "bismuth             ", &
         "polonium            ", &
         "astatine            ", &
         "radon               ", &
         "francium            ", &
         "radium              ", &
         "actinium            ", &
         "thorium             ", &
         "protactinium        ", &
         "uranium             ", &
         "neptunium           ", &
         "plutonium           ", &
         "americium           ", &
         "curium              ", &
         "berkelium           ", &
         "californium         ", &
         "einsteinium         ", &
         "fermium             ", &
         "mendelevium         ", &
         "nobelium            ", &
         "lawrencium          ", &
         "rutherfordium       ", &
         "dubnium             ", &
         "seaborgium          ", &
         "bohrium             ", &
         "hassium             ", &
         "meitnerium          ", &
         "darmstadtium        ", &
         "roentgentium        ", &
         "ununbium            ", &
         "ununtrium           ", &
         "ununquadium         ", &
         "ununpentium         ", &
         "ununhexium          ", &
         "ununseptium         "/) !MAY BE INCRRCT

!       SpeciesList%AtomElectroneg_Pauling=(/2.2_dp,-999._dp,0.98_dp,1.57_dp,2.04_dp,2.55_dp,&
!          3.04_dp,3.44_dp,3.98_dp,-999._dp,0.93_dp,1.31_dp,1.61_dp,1.9_dp,2.19_dp,&
!          2.58_dp,3.16_dp,-999._dp,0.82_dp,1._dp,1.36_dp,1.54_dp,1.63_dp,1.66_dp, &
!          1.55_dp,1.83_dp,1.88_dp,1.91_dp,1.9_dp,1.65_dp,1.81_dp,2.01_dp,2.18_dp,&
!          2.55_dp,2.96_dp,3._dp,0.82_dp,0.95_dp,1.22_dp,1.33_dp,1.6_dp,2.16_dp,1.9_dp,&
!          2.2_dp,2.28_dp,2.2_dp,1.93_dp,1.69_dp,1.78_dp,1.96_dp,2.05_dp,2.1_dp,2.66_dp,&
!          2.6_dp,0.79_dp,0.89_dp,1.1_dp,1.12_dp,1.13_dp,1.14_dp,-999._dp,1.17_dp,&
!          -999._dp,1.2_dp,-999._dp,1.22_dp,1.23_dp,1.24_dp,1.25_dp,-999._dp,1.27_dp,&
!          1.3_dp,1.5_dp,2.36_dp,1.9_dp,2.2_dp,2.2_dp,2.28_dp,2.54_dp,2._dp,1.62_dp,&
!          2.33_dp,2.02_dp,2._dp,2.2_dp,-999._dp,0.7_dp,0.9_dp,1.1_dp,1.3_dp,1.5_dp,&
!          1.38_dp,1.36_dp,1.28_dp,1.3_dp,1.3_dp,1.3_dp,1.3_dp,1.3_dp,1.3_dp,1.3_dp,&
!          1.3_dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,&
!          -999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp/)

!       SpeciesList%Atomic1stIonization=(/1312._dp,2372._dp,520._dp,899._dp,&
!          801._dp,1086._dp,1402._dp,&
!          1314._dp,1681._dp,2081._dp,496._dp,738._dp,578._dp,&
!          787._dp,1012._dp,1000._dp,1251._dp,&
!          1521._dp,419._dp,590._dp,633._dp,659._dp,651._dp,&
!          653._dp,717._dp,762._dp,760._dp,737._dp,&
!          745._dp,906._dp,579._dp,762._dp,947._dp,941._dp,&
!          1140._dp,1351._dp,403._dp,549._dp,600._dp,&
!          640._dp,652._dp,684._dp,702._dp,710._dp,720._dp,&
!          804._dp,731._dp,868._dp,558._dp,709._dp,&
!          834._dp,869._dp,1008._dp,1170._dp,376._dp,503._dp,&
!          538._dp,534._dp,527._dp,533._dp,535._dp,&
!          545._dp,547._dp,593._dp,569._dp,573._dp,581._dp,&
!          589._dp,597._dp,603._dp,524._dp,659._dp,&
!          761._dp,770._dp,760._dp,839._dp,878._dp,868._dp,&
!          890._dp,1007._dp,589._dp,716._dp,703._dp,&
!          812._dp,917._dp,1037._dp,380._dp,509._dp,499._dp,&
!          587._dp,568._dp,598._dp,605._dp,585._dp,&
!          578._dp,581._dp,601._dp,608._dp,619._dp,627._dp,&
!          635._dp,642._dp,-999._dp,-999._dp,-999._dp,&
!          -999._dp,-999._dp,-999._dp,-999._dp,-999._dp,&
!          -999._dp,-999._dp,-999._dp,-999._dp,-999._dp,&
!          -999._dp,-999._dp,-999._dp/) !in kJ/mol

!       SpeciesList%AtomicDensity=(/0.0000699_dp,0.000179_dp,0.543_dp, &
!          1.85_dp,2.34_dp,2.25_dp,0.00125_dp,0.00143_dp,0.0017_dp, &
!          0.0009_dp,0.971_dp,1.74_dp,2.7_dp,2.33_dp,1.82_dp,2.07_dp, &
!          0.00321_dp,0.00178_dp,0.86_dp,1.55_dp,2.99_dp,4.54_dp, &
!          6.11_dp,7.19_dp,7.43_dp,7.86_dp,8.9_dp,8.9_dp,8.96_dp, &
!          7.13_dp,5.9_dp,5.32_dp,5.73_dp,4.79_dp,3.12_dp,0.00374_dp,&
!          1.53_dp,2.54_dp,4.47_dp,6.51_dp,8.57_dp,10.2_dp,11.5_dp,&
!          12.4_dp,12.4_dp,12._dp,10.5_dp,8.65_dp,7.31_dp,7.31_dp,&
!          6.69_dp,6.24_dp,4.93_dp,0.00589_dp,1.87_dp,3.5_dp,6.15_dp,&
!          6.66_dp,6.77_dp,7._dp,7.26_dp,7.52_dp,5.24_dp,7.9_dp,8.23_dp,&
!          8.55_dp,8.8_dp,9.07_dp,9.32_dp,6.97_dp,9.84_dp,13.3_dp,&
!          16.6_dp,19.3_dp,21._dp,22.6_dp,22.4_dp,21.4_dp,19.3_dp,13.5_dp,&
!          11.9_dp,11.4_dp,9.75_dp,9.32_dp,-999._dp,0.00973_dp,-999._dp,&
!          5._dp,10.1_dp,11.7_dp,15.4_dp,19._dp,20.2_dp,19.8_dp,13.7_dp,&
!          13.5_dp,14._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,&
!          -999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,&
!          -999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp/) !in g/ml

!       SpeciesList%AtomicRadius=(/37.1_dp,31._dp,152._dp,112._dp,85._dp,77.2_dp, &
!          70._dp,73._dp,72._dp,71._dp,186._dp,160._dp,143._dp,117.6_dp, &
!          110._dp,103._dp,100._dp,98._dp,227._dp,197._dp,162._dp,147._dp,&
!          134._dp,128._dp,127._dp,126._dp,125._dp,124._dp,128._dp,134._dp,&
!          135._dp,122.3_dp,120._dp,119._dp,114._dp,112._dp,248._dp,&
!          215._dp,180._dp,160._dp,146._dp,139._dp,136._dp,134._dp,&
!          134._dp,137._dp,144._dp,151._dp,167._dp,140.5_dp,140._dp,&
!          142._dp,133._dp,131._dp,265._dp,222._dp,187._dp,182._dp,&
!          182._dp,181._dp,183._dp,180._dp,208._dp,180._dp,177._dp,178._dp,&
!          176._dp,176._dp,176._dp,193._dp,174._dp,159._dp,146._dp,&
!          139._dp,137._dp,135._dp,136._dp,139._dp,144._dp,151._dp,&
!          170._dp,146._dp,150._dp,168._dp,-999._dp,-999._dp,-999._dp,&
!          -999._dp,-999._dp,179._dp,163._dp,156._dp,155._dp,159._dp,&
!          173._dp,174._dp,170._dp,186._dp,186._dp,-999._dp,-999._dp,&
!          -999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,&
!          -999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,-999._dp,&
!          -999._dp,-999._dp,-999._dp/) !in picometer

      WHERE(SpeciesList%AtomicRadius>0._dp) SpeciesList%AtomicRadius=SpeciesList%AtomicRadius*1.e-2_dp !in Ang

      SpeciesList%AtomicOxidState=(/ &
         "x1          ", &
         "0           ", &
         "+1          ", &
         "+2          ", &
         "+3          ", &
         "x4          ", &
         "-3          ", &
         "-2          ", &
         "-1          ", &
         "0           ", &
         "+1          ", &
         "+2          ", &
         "+3          ", &
         "x4          ", &
         "-3          ", &
         "-2          ", &
         "-1          ", &
         "0           ", &
         "+1          ", &
         "+2          ", &
         "+3          ", &
         "+4,3,2      ", &
         "+5,2,3,4    ", &
         "+3,2,6      ", &
         "+2,3,4,6,7  ", &
         "+3,2        ", &
         "+2,3        ", &
         "+2,3        ", &
         "+2,1        ", &
         "+2          ", &
         "+3          ", &
         "+4,2        ", &
         "x3,+5       ", &
         "+4,-2,+6    ", &
         "x1,+5       ", &
         "0           ", &
         "+1          ", &
         "+2          ", &
         "+3          ", &
         "+4          ", &
         "+5,3        ", &
         "+6,3,5      ", &
         "+7,4,6      ", &
         "+4,3,6,8    ", &
         "+3,4,6      ", &
         "+2,4        ", &
         "+1          ", &
         "+2          ", &
         "+3          ", &
         "+4,2        ", &
         "+3,5        ", &
         "+4,6,-2     ", &
         "-1,+5,7     ", &
         "0           ", &
         "+1          ", &
         "+2          ", &
         "+3          ", &
         "+3,4        ", &
         "+3,4        ", &
         "+3          ", &
         "+3          ", &
         "+3,2        ", &
         "+3,2        ", &
         "+3          ", &
         "+3,4        ", &
         "+3          ", &
         "+3          ", &
         "+3          ", &
         "+3,2        ", &
         "+3,2        ", &
         "+3          ", &
         "+4          ", &
         "+5          ", &
         "+6,4        ", &
         "+7,4,6      ", &
         "+4,6,8      ", &
         "+4,3,6      ", &
         "+4,2        ", &
         "+3,1        ", &
         "+2,1        ", &
         "+1,3        ", &
         "+2,4        ", &
         "+3,5        ", &
         "+4,2        ", &
         "-999        ", &
         "0           ", &
         "+1          ", &
         "+2          ", &
         "+3          ", &
         "+4          ", &
         "+5,4        ", &
         "+6,3,4,5    ", &
         "+5,3,4,6    ", &
         "+4,3,5,6    ", &
         "+3,4,5,6    ", &
         "+3          ", &
         "+3,4        ", &
         "+3          ", &
         "+3          ", &
         "+3          ", &
         "+3,2        ", &
         "+2,3        ", &
         "+3          ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        ", &
         "-999        "/)

   END SUBROUTINE SetUpSpecies
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE PrintSpeciesDb()
   
      IMPLICIT NONE
      INTEGER :: i
      
      IF (.NOT. ASSOCIATED(SpeciesList)) THEN
         WRITE(*,*) "$Err>> Species list is not initialized"
         STOP
      END IF
      
      WRITE(*,*) ">> Entering species database ..."
      WRITE(*,*) ">> number of elements in db:",SpeciesList%NSpecies
      
      DO i=1,SpeciesList%NSpecies
         WRITE(UNIT=*,FMT='(i5,"    |",a5)') SpeciesList%AtomicNumber(i),SpeciesList%AtomicSymbol(i)
      END DO
      
   END SUBROUTINE PrintSpeciesDb
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
END MODULE Species
