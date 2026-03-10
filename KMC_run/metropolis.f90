MODULE metropolis
   USE VARIABLE_TYPE
   USE PotentialPackage
   USE Ecuyer_random
   
   IMPLICIT NONE
   
   CONTAINS
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE MetropolisNVT(AL,ntherm,niter,drmax,Temperature,iseed,iprint,iprintfreq)
   
   
   !To do: Make only moving atoms get selected
   
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL
      INTEGER, OPTIONAL :: iprint,iprintfreq
      INTEGER :: niter,ntherm,iter,iatom,NAtoms,iseed,iprint1,iprintfreq1,nreject
      REAL(dp) :: drmax,dr(3),random_no,Temperature,coord(3),enew,eprev,kbT,probability
      REAL(dp), DIMENSION(:), POINTER :: AtomCoord
      REAL(dp) :: en,en2
      
      IF (PRESENT(iprint)) THEN
         iprint1=iprint
      ELSE
         iprint1=0
      END IF
      
      IF (PRESENT(iprintfreq)) THEN
         iprintfreq1=iprintfreq
      ELSE
         iprintfreq1=0
      END IF
      
      NAtoms=AL%NAtoms
      kbT=kboltzmann*Temperature
      
      en=0._dp
      en2=0._dp
      
      WRITE(6,*) "Performing Metropolis NVT calculation ..."
      AtomCoord=>AL%AtomCoord
      
      CALL InitializeTaus88(iseed)
      
      CALL AddVerletList(AL,NAtoms=NAtoms)
      CALL GetEnergyVL(AL,NAtoms=NAtoms)
      
      nreject=0
      WRITE(6,*) " ... thermalizing system (printing energy (eV))"
      DO iter=1,ntherm
         random_no=taus88()
         iatom=CEILING(random_no*REAL(NAtoms))
         random_no=taus88()
         dr(1)=drmax*(2._dp*random_no-1)
         random_no=taus88()
         dr(2)=drmax*(2._dp*random_no-1)
         random_no=taus88()
         dr(3)=drmax*(2._dp*random_no-1)
         !write(*,*) dr
         !stop
         coord=AtomCoord(3*iatom-2:3*iatom)
         eprev=AL%PotentialEnergy
         
         AtomCoord(3*iatom-2:3*iatom)=coord+dr
         CALL AddVerletList(AL,NAtoms=NAtoms)
         CALL GetEnergyVL(AL,NAtoms=NAtoms)
         enew=AL%PotentialEnergy
         
         probability=EXP(-(enew-eprev)/kBT)
         random_no=taus88()
         IF (random_no>probability) THEN
            AtomCoord(3*iatom-2:3*iatom)=coord !reject the move
            nreject=nreject+1
         END IF
         
         IF (MOD(iter,iprintfreq1)==1) WRITE(6,*) enew
      END DO
      WRITE(6,*) " ... fraction moves rejected (%)",REAL(nreject)*100./REAL(ntherm)
      
      nreject=0
      WRITE(6,*) " "
      WRITE(6,*) " ... sampling configuration space "
      DO iter=1,niter
         random_no=taus88()
         iatom=CEILING(random_no*REAL(NAtoms))
         random_no=taus88()
         dr(1)=drmax*random_no
         random_no=taus88()
         dr(2)=drmax*random_no
         random_no=taus88()
         dr(3)=drmax*random_no
         coord=AtomCoord(3*iatom-2:3*iatom)
         eprev=AL%PotentialEnergy
         
         AtomCoord(3*iatom-2:3*iatom)=coord+dr
         CALL AddVerletList(AL,NAtoms=NAtoms)
         CALL GetEnergyVL(AL,NAtoms=NAtoms)
         enew=AL%PotentialEnergy
         
         probability=EXP(-(enew-eprev)/kBT)
         random_no=taus88()
         IF (random_no>probability) THEN
            AtomCoord(3*iatom-2:3*iatom)=coord
            nreject=nreject+1
         END IF
         IF (MOD(iter,iprintfreq1)==1 .AND. iprint1>0) WRITE(6,*) enew
         en=en+enew
         en2=en2+enew*enew
      END DO
      WRITE(6,*) " ... fraction moves rejected in sampling stage(%)",REAL(nreject)*100./REAL(ntherm)
      WRITE(*,*) "Sum energy:",en
      WRITE(*,*) "Sum energy2:",en2
   END SUBROUTINE MetropolisNVT
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
END MODULE metropolis
