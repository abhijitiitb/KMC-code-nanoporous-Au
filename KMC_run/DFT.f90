MODULE DFTSubroutines
!interface for performing DFT calculations
   USE VARIABLE_TYPE
   IMPLICIT NONE
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DFTForces(AL) !,filename_in,filename_out)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord,AtomForce
      INTEGER :: dir,iatom,NAtoms
      CHARACTER(len=100) :: filename_in,filename_out
      CHARACTER(len=3) :: speciesname
      LOGICAL :: Satisfied

      CALL SYSTEM("rm dft_status")
      filename_in="Si.scf.inp"
      filename_out="Si.scf.out"

      !create the input file
      OPEN(UNIT=423,FILE=filename_in)
      WRITE(423,*) "&CONTROL"
      WRITE(423,*) "calculation = 'scf' ,"
      WRITE(423,*) "restart_mode = 'from_scratch' ,"
      WRITE(423,*) "prefix = 'silicon' ,"
      WRITE(423,*) "tstress = .false. ,"
      WRITE(423,*) "tprnfor = .true. ,"
      WRITE(423,*) "pseudo_dir = '/opt/apps/espresso-5.0.1/pseudo',"
      WRITE(423,*) "outdir= '/home/abhijit/tmp'"
      WRITE(423,*) "/"
      
      WRITE(423,*) "&SYSTEM"
      WRITE(423,*) "ibrav = 1,"
      WRITE(UNIT=423,FMT='(" celldm(1) = ",3F6.2)') AL%BoxSize(1)/0.529177249_dp
      WRITE(UNIT=423,FMT='(" nat = ",i4,",")') AL%NAtoms
      WRITE(423,*) "ntyp = 1,"
      WRITE(423,*) "ecutwfc = 20,"
      WRITE(423,*) "/"
      
      WRITE(423,*) "&ELECTRONS"
      WRITE(423,*) "diagonalization='cg',"
      WRITE(423,*) "mixing_mode = 'plain' ,"
      WRITE(423,*) "conv_thr = 1.0d-8"
      WRITE(423,*) "/"
      WRITE(423,*) "ATOMIC_SPECIES"
      WRITE(423,*) "Si 28.086 Si.rel-pbe-rrkj.UPF"
      WRITE(423,*) "ATOMIC_POSITIONS crystal"
      DO iatom = 1,AL%NAtoms
         speciesname=SpeciesList%AtomicSymbol(SpeciesDirectory_global(AL%AtomSpecies(iatom)))
         WRITE(UNIT=423,FMT='(a3,3f8.3)') speciesname, &
         (AL%AtomCoord(3*iatom-2:3*iatom)-AL%BoxSize(4:6))/ &
         (AL%BoxSize(1:3)-AL%BoxSize(4:6))
      END DO
      WRITE(423,*) "K_POINTS automatic"
      WRITE(423,*) "2 2 2  0 0 0"
      CLOSE(423)
       
      !run Q-espresso
      CALL SYSTEM("echo 'mpiexec -n 1 /opt/apps/espresso-5.0.1/bin/pw.x < "// &
        TRIM(filename_in)//" > "//TRIM(filename_out)//"' > qe_script.x")
      CALL SYSTEM("echo 'echo completed > dft_status' | cat >> qe_script.x")
      CALL SYSTEM("chmod +x qe_script.x; ./qe_script.x")
      
      Satisfied=.FALSE.
      DO WHILE (.NOT. Satisfied)
         INQUIRE(FILE="dft_status",EXIST=Satisfied)
      END DO

      !read the output file
      CALL SYSTEM("grep 'force =' "//TRIM(filename_out)//" > "//TRIM(filename_out)//".tmp")
      CALL SYSTEM("awk '{print $7}' "//TRIM(filename_out)//".tmp > "//TRIM(filename_out)//".tmp1") !x-forces
      CALL SYSTEM("awk '{print $8}' "//TRIM(filename_out)//".tmp | cat >> "//TRIM(filename_out)//".tmp1") !y-forces
      CALL SYSTEM("awk '{print $9}' "//TRIM(filename_out)//".tmp | cat >> "//TRIM(filename_out)//".tmp1") !z-forces
      
      !obtain energy
      CALL SYSTEM("grep '!    total energy' "//TRIM(filename_out)//" > "//TRIM(filename_out)//".tmp")
      CALL SYSTEM("awk '{print $5}' "//TRIM(filename_out)//".tmp > "//TRIM(filename_out)//".tmp1")
      OPEN(UNIT=423,FILE=TRIM(filename_out)//".tmp1")
      READ(423,*) AL%PotentialEnergy
      AL%PotentialEnergy=AL%PotentialEnergy*13.6056981
      CLOSE(423)
      
      WRITE(*,*) AL%PotentialEnergy
      stop
      !obtain forces
      NAtoms=AL%NAtoms
      AtomForce=>AL%AtomForce
      OPEN(UNIT=423,FILE=TRIM(filename_out)//".tmp1")
      DO dir=1,3
         DO iatom=1,NAtoms
            READ(423,*) AtomForce(3*(iatom-1)+dir)
            AtomForce(3*(iatom-1)+dir)=AtomForce(3*(iatom-1)+dir)*13.6056981/0.529177249 !from Ry/bohr to eV/A
         END DO
      END DO
      CLOSE(423)
   END SUBROUTINE DFTForces
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE DFTSubroutines

