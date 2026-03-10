MODULE InteratomicPotentialsAssess
!contains subroutines to assess various aspects of an interatomic potential
   USE VARIABLE_TYPE
   USE PotentialPackage
   USE NeighborList
   USE Crystal
   IMPLICIT NONE
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EnergyDensityPlot(AL,volume_min,volume_max,npoints,iunit)
   !changes the existing density and finds the energy of the resulting system
   !volume_max and volume_min are the relative volumes with respect to the current volume
   !npoints is the number of points for which the energy has to be computed
      IMPLICIT NONE
      INTEGER, PARAMETER :: maxnpoints=100
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp) :: volume_max,volume_min,dv,sclfac1(3)
      REAL(dp) :: volume(maxnpoints)
      INTEGER, OPTIONAL :: iunit
      INTEGER :: iunit1,npoints,iexpansion,ipoint,errorstatus
      
      IF (npoints>maxnpoints) THEN
         WRITE(6,*) "Err>> Increase maxnpoints in EnergyDensityPlot"
         STOP
      END IF
      
      IF (volume_min>=volume_max) THEN
         WRITE(6,*) "Err>> volume_min>=volume_max in EnergyDensityPlot"
         STOP
      END IF
      
      IF (PRESENT(iunit)) THEN
         iunit1=iunit
      ELSE
         iunit1=6
      END IF
      
      !set up the volumes
      dv=(volume_max-volume_min)/REAL(npoints-2,dp)
      volume(1)=volume_min
      iexpansion=0
      IF (iexpansion==0 .AND. volume(1)>1._dp) iexpansion=1 !expansion occurs here
      DO ipoint=2,npoints-1
         volume(ipoint)=volume(ipoint-1)+dv
         IF (iexpansion==0 .AND. volume(ipoint)>1._dp) iexpansion=ipoint !expansion occurs here
      END DO
      IF (iexpansion>0) THEN
         volume(iexpansion+1:npoints)=volume(iexpansion:npoints-1)
         volume(iexpansion)=1._dp
      ELSE
         volume(npoints)=1._dp
      END IF
      
      DO ipoint=1,npoints
         IF (ipoint==1) THEN
            sclfac1=volume(1)
         ELSE
            sclfac1=volume(ipoint)/volume(ipoint-1)
         END IF
         CALL Scl(AL,sclfac1)
         CALL AddVerletList(AL,ListType=AL%VL%ListType,NAtoms=AL%NAtoms)
         CALL GetForcesVL(AL,AL%NAtoms,errorstatus)
         IF (errorstatus/=0) THEN
            WRITE(6,*) "Err>> Error encountered while finding energy"
            STOP
         END IF
         WRITE(iunit1,*) volume(ipoint),AL%PotentialEnergy
      END DO
   END SUBROUTINE EnergyDensityPlot
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE InteratomicPotentialsAssess
