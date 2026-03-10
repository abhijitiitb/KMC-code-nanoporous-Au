!uses the improved tangent generation procedure

ChOSTangent(1:3*NAtoms)=0._dp
Vprev=AL%PrevNeigh%PotentialEnergy
Vcurr=AL%PotentialEnergy
Vnext=AL%NextNeigh%PotentialEnergy
Vnext=AL%NextNeigh%PotentialEnergy
Vcurr=AL%PotentialEnergy
Vprev=AL%PrevNeigh%PotentialEnergy
IF (Vnext>Vcurr .AND. Vcurr>Vprev) THEN !reactant side
   ChOSTangent(1:3*NAtoms)=PBCdistance(AL%NextNeigh,AL)
   !r1=-2
   !r2=0
   !DO i=1,NAtoms
   !   r1=r1+3
   !   r2=r2+3
   !   ChOSTangent(r1:r2)= &
   !      PBCdistance(AL,AtomCoordNext(r1:r2),AL%AtomCoord(r1:r2))
   !END DO
ELSEIF (Vnext<Vcurr .AND. Vcurr<Vprev) THEN !product side
   ChOSTangent(1:3*NAtoms)=PBCdistance(AL,AL%PrevNeigh)
   !r1=-2
   !r2=0
   !DO i=1,NAtoms
   !   r1=r1+3
   !   r2=r2+3
   !   WHERE (AL%AtomIsMoving(r1:r2)) ChOSTangent(r1:r2)= &
   !      PBCdistance(AL,AL%AtomCoord(r1:r2),AL%PrevNeigh%AtomCoord(r1:r2))
   !END DO
ELSEIF (Vnext>Vprev) THEN
   dVmax=MAX(ABS(Vnext-Vcurr),ABS(Vprev-Vcurr))
   dVmin=MIN(ABS(Vnext-Vcurr),ABS(Vprev-Vcurr))
   ChOSTangent(1:3*NAtoms)=PBCdistance(AL%NextNeigh,AL)*dVmax+PBCdistance(AL,AL%PrevNeigh)*dVmin
   !r1=-2
   !r2=0
   !DO i=1,NAtoms
   !   r1=r1+3
   !   r2=r2+3
   !   WHERE (AL%AtomIsMoving(r1:r2)) ChOSTangent(r1:r2)= &
   !      PBCdistance(AL,AL%NextNeigh%AtomCoord(r1:r2),AL%AtomCoord(r1:r2))*dVmax+ &
   !      PBCdistance(AL,AL%AtomCoord(r1:r2),AL%PrevNeigh%AtomCoord(r1:r2))*dVmin
   !END DO
ELSEIF (Vnext<Vprev) THEN !maximum
   dVmax=MAX(ABS(Vnext-Vcurr),ABS(Vprev-Vcurr))
   dVmin=MIN(ABS(Vnext-Vcurr),ABS(Vprev-Vcurr))
   ChOSTangent(1:3*NAtoms)=PBCdistance(AL%NextNeigh,AL)*dVmin+PBCdistance(AL,AL%PrevNeigh)*dVmax
   !r1=-2
   !r2=0
   !DO i=1,NAtoms
   !   r1=r1+3
   !   r2=r2+3
   !   WHERE (AL%AtomIsMoving(r1:r2)) ChOSTangent(r1:r2)= &
   !      PBCdistance(AL,AL%NextNeigh%AtomCoord(r1:r2),AL%AtomCoord(r1:r2))*dVmin+ &
   !      PBCdistance(AL,AL%AtomCoord(r1:r2),AL%PrevNeigh%AtomCoord(r1:r2))*dVmax
   !END DO
ELSE
   WRITE(6,*) "$Err>> Exception occured while estimating tangent using the NEB improved tangent method"
   WRITE(UNIT=6,FMT='(" ...Energy of left   image:",ES15.6)') Vprev
   WRITE(UNIT=6,FMT='(" ...Energy of center image:",ES15.6)') Vcurr
   WRITE(UNIT=6,FMT='(" ...Energy of right  image:",ES15.6)') Vnext
   STOP
END IF

normt=NORM(ChOSTangent,3*NAtoms)
ChOSTangent(1:3*NAtoms)=ChOSTangent(1:3*NAtoms)/normt
