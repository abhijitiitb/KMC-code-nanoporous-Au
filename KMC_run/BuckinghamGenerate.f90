MODULE BuckinghamGenerate

   IMPLICIT NONE
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateBuckinghamTables()
      IMPLICIT NONE
      INTEGER, PARAMETER :: dp = KIND(1.0d0)
      REAL(dp) :: xstart,xstop,dx,E1,E2,E3,E4,r1,r2,r3,r4,delr
      INTEGER :: npts,j
      
      REAL(dp), PARAMETER :: AZR=1024.6_dp
      REAL(dp), PARAMETER :: rhoZR=0.376_dp
      REAL(dp), PARAMETER :: cZR=0.0_dp
      REAL(dp), PARAMETER :: AY=1325.6_dp
      REAL(dp), PARAMETER :: rhoY=0.3461_dp
      REAL(dp), PARAMETER :: cY=0._dp
      REAL(dp), PARAMETER :: AO=22764.3_dp
      REAL(dp), PARAMETER :: rhoO=0.149_dp
      REAL(dp), PARAMETER :: cO=27.89_dp
   
   
      npts=3000
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Create phi file - pair potential
   
      OPEN(UNIT=200,FILE="Y_O.TXT")
      WRITE(200,*) "Y_O..."
      xstart=0.001_dp !Angstrom
      xstop=6.4729_dp !cut off radius for Yttrium-Oxygen at E=10^(-5)ev 
      r1=xstart
      dx=(xstop-xstart)/REAL(npts,dp)
      delr=(1/dx)
      r2=xstart+dx
      r3=xstart+2*dx
      r4=xstart+3*dx
      WRITE(200,*) npts, xstart,xstop,delr 
      
      
      
      DO j=1,3000,4
         E1=AY*exp(-r1/rhoY)-(cY/r1**6) !eV
         E2=AY*exp(-r2/rhoY)-(cY/r2**6) !ev 
         E3=AY*exp(-r3/rhoY)-(cY/r3**6) !ev 
         E4=AY*exp(-r4/rhoY)-(cY/r4**6) !ev 
         WRITE(200,75) E1,E2,E3,E4
         75 FORMAT(1X, E20.12, E20.12, E20.12, E20.12)
         r1=r1+(4)*dx
         r2=r2+(4)*dx
         r3=r3+(4)*dx
         r4=r4+(4)*dx
     
      END DO
      close(200)
   
      OPEN(UNIT=200,FILE="O_O.TXT")
      WRITE(200,*) "O_O..."
      xstart=0.001_dp !Angstrom
      xstop=11.9_dp !cut off radius for Oxygen-Oxygen at E=10^(-5)ev 
      r1=xstart
      dx=(xstop-xstart)/REAL(npts,dp)
      delr=(1/dx)
      r2=xstart+dx
      r3=xstart+2*dx
      r4=xstart+3*dx
      WRITE(200,*) npts, xstart,xstop,delr
   
  
      DO j=1,3000,4
         E1=AO*exp(-r1/rhoO)-(cO/r1**6) !eV
         E2=AO*exp(-r2/rhoO)-(cO/r2**6) !eV
         E3=AO*exp(-r3/rhoO)-(cO/r3**6) !eV
         E4=AO*exp(-r4/rhoO)-(cO/r4**6) !eV
         WRITE(200,70) E1,E2,E3,E4
         70 FORMAT(1X, E20.12, E20.12, E20.12, E20.12)
         r1=r1+(4)*dx
         r2=r2+(4)*dx
         r3=r3+(4)*dx
         r4=r4+(4)*dx
     
      END DO
      close(200)
   
   
      OPEN(UNIT=200,FILE="Zr_O.TXT") 
      WRITE(200,*) "Zr_O..."
      xstart=0.001_dp !Angstrom
      xstop=6.9353_dp !cut off radius for zirconium-Oxygen at E=10^(-5)ev 
      r1=xstart
      dx=(xstop-xstart)/REAL(npts,dp)
      delr=(1/dx)
      r2=xstart+dx
      r3=xstart+2*dx
      r4=xstart+3*dx
      WRITE(200,*) npts, xstart,xstop,delr 
   

   
      DO j=1,3000,4 
         E1=AZR*exp(-r1/rhoZR)-(cZR/r1**6) !eV
         E2=AZR*exp(-r2/rhoZR)-(cZR/r2**6) !eV
         E3=AZR*exp(-r3/rhoZR)-(cZR/r3**6) !eV
         E4=AZR*exp(-r4/rhoZR)-(cZR/r4**6) !eV
         WRITE(200,80) E1,E2,E3,E4
         80 FORMAT(1X, E20.12, E20.12, E20.12, E20.12)
         r1=r1+(4)*dx
         r2=r2+(4)*dx
         r3=r3+(4)*dx
         r4=r4+(4)*dx
     
      END DO
      close(200)
   END SUBROUTINE GenerateBuckinghamTables
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE BuckinghamGenerate
