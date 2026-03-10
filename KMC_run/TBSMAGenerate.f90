MODULE TBSMAGenerate
!pgm to calculate out the Tight Binding second moment approximation of Ag-Ag, Pt-Pt Ag-Pt

   IMPLICIT NONE
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateTBSMATables()
      IMPLICIT NONE
      INTEGER, PARAMETER :: dp = KIND(1.0d0)
      REAL(dp) :: xstart,xstop,dx,E1,E2,E3,E4,r1,r2,r3,r4,delr
      INTEGER :: npts,j
      
      REAL(dp), PARAMETER :: A_AgAg=1.9994_dp !lambda
      REAL(dp), PARAMETER :: A_AgPt=2.4136_dp
      REAL(dp), PARAMETER :: p_AgAg=6.2177_dp
      REAL(dp), PARAMETER :: p_AgPt=6.5919_dp
      REAL(dp), PARAMETER :: r0_AgAg=2.89_dp
      REAL(dp), PARAMETER :: r0_AgPt=2.77_dp
      REAL(dp), PARAMETER :: q_AgAg=2.9138_dp
      REAL(dp), PARAMETER :: q_AgPt=3.0612_dp
      REAL(dp), PARAMETER :: epsilon_AgAg=2.5096_dp
      REAL(dp), PARAMETER :: epsilon_AgPt=3.1546_dp
      REAL(dp), PARAMETER :: alpha_AgAg=0.9518_dp
      REAL(dp), PARAMETER :: alpha_AgPt=0.9427_dp
      
      npts=2000
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Create phi file - pair potential
   
      !Ag-Ag
      OPEN(UNIT=200,FILE="AgAg_phi.TXT")
      WRITE(6,*) "AgAg_phi..."
      xstart=0.001_dp !Angstrom
      xstop=4.08_dp !cut off radius for Ag-Ag
      r1=xstart
      dx=(xstop-xstart)/REAL(npts,dp)
      delr=(1/dx)
      r2=xstart+dx
      r3=xstart+2*dx
      r4=xstart+3*dx
      WRITE(200,*) npts, xstart,xstop,delr 
   
      DO j=1,npts,4
         E1=A_AgAg*exp(-p_AgAg*(r1/r0_AgAg-1)) !eV 
         E2=A_AgAg*exp(-p_AgAg*(r2/r0_AgAg-1)) !ev 
         E3=A_AgAg*exp(-p_AgAg*(r3/r0_AgAg-1)) !ev 
         E4=A_AgAg*exp(-p_AgAg*(r4/r0_AgAg-1)) !ev 
         WRITE(200,75) E1,E2,E3,E4
         75 FORMAT(1X, E20.12, E20.12, E20.12, E20.12)
         r1=r1+(4)*dx
         r2=r2+(4)*dx
         r3=r3+(4)*dx
         r4=r4+(4)*dx
     
      END DO
      close(200)
   
      !Ag-Pt 
      OPEN(UNIT=300,FILE="AgPt_phi.TXT")
      WRITE(6,*) "AgPt_phi..."
      xstart=0.001_dp !Angstrom
      xstop=4.08_dp !cut off radius for Ag-Pt
      r1=xstart
      dx=(xstop-xstart)/REAL(npts,dp)
      delr=(1/dx)
      r2=xstart+dx
      r3=xstart+2*dx
      r4=xstart+3*dx
      WRITE(300,*) npts, xstart,xstop,delr 
      
      DO j=1,npts,4
         E1=A_AgPt*exp(-p_AgPt*(r1/r0_AgPt-1)) !eV 
         E2=A_AgPt*exp(-p_AgPt*(r2/r0_AgPt-1)) !ev 
         E3=A_AgPt*exp(-p_AgPt*(r3/r0_AgPt-1)) !ev 
         E4=A_AgPt*exp(-p_AgPt*(r4/r0_AgPt-1)) !ev 
         WRITE(300,80) E1,E2,E3,E4
         80 FORMAT(1X, E20.12, E20.12, E20.12, E20.12)
         r1=r1+(4)*dx
         r2=r2+(4)*dx
         r3=r3+(4)*dx
         r4=r4+(4)*dx
     
      END DO
      close(300)
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Create rho file 
  
      !Ag-Ag
      OPEN(UNIT=200,FILE="AgAg_rho.TXT")
      WRITE(6,*) "AgAg_rho..."
      xstart=0.001_dp !Angstrom
      xstop=4.08_dp !cut off radius for Ag-Ag
      r1=xstart
      dx=(xstop-xstart)/REAL(npts,dp)
      delr=(1/dx)
      r2=xstart+dx
      r3=xstart+2*dx
      r4=xstart+3*dx
      WRITE(200,*) npts, xstart,xstop,delr
   
      DO j=1,npts,4
         E1=exp(-2*q_AgAg*(r1/r0_AgAg-1)) !eV
         E2=exp(-2*q_AgAg*(r2/r0_AgAg-1)) !eV
         E3=exp(-2*q_AgAg*(r3/r0_AgAg-1)) !eV
         E4=exp(-2*q_AgAg*(r4/r0_AgAg-1)) !eV
         WRITE(200,85) E1,E2,E3,E4
         85 FORMAT(1X, E20.12, E20.12, E20.12, E20.12)
         r1=r1+(4)*dx
         r2=r2+(4)*dx
         r3=r3+(4)*dx
         r4=r4+(4)*dx
     
      END DO
      close(200)
 
      !Ag-Pt
      OPEN(UNIT=300,FILE="AgPt_rho.TXT")
      WRITE(6,*) "AgPt_rho..."
      xstart=0.001_dp !Angstrom
      xstop=4.08_dp !cut off radius for Ag-Pt
      r1=xstart
      dx=(xstop-xstart)/REAL(npts,dp)
      delr=(1/dx)
      r2=xstart+dx
      r3=xstart+2*dx
      r4=xstart+3*dx
      WRITE(300,*) npts, xstart,xstop,delr
   
      DO j=1,npts,4
         E1=exp(-2*q_AgPt*(r1/r0_AgPt-1)) !eV   
         E2=exp(-2*q_AgPt*(r2/r0_AgPt-1)) !eV
         E3=exp(-2*q_AgPt*(r3/r0_AgPt-1)) !eV
         E4=exp(-2*q_AgPt*(r4/r0_AgPt-1)) !eV
         WRITE(300,90) E1,E2,E3,E4
         90 FORMAT(1X, E20.12, E20.12, E20.12, E20.12)
         r1=r1+(4)*dx
         r2=r2+(4)*dx
         r3=r3+(4)*dx
         r4=r4+(4)*dx
     
      END DO
      close(200)
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Create f file  
  
      !Ag-Ag
   
      OPEN(UNIT=200,FILE="AgAg_f.TXT") 
      WRITE(6,*) "AgAg_f..."
      xstart=338.859_dp ! corresponding to r=0.001 Angstrong
      xstop=0.0907_dp !corresponding to cut off radius for Ag-Ag r_c=4.08 Angstrong
      r1=xstart
      dx=(xstop-xstart)/REAL(npts,dp)
      delr=(1/dx)
      r2=xstart+dx
      r3=xstart+2*dx
      r4=xstart+3*dx
      WRITE(200,*) npts, xstart,xstop,delr 
   
      DO j=1,npts,4 
         E1=epsilon_AgAg*((exp(-2*q_AgAg*(r1/r0_AgAg-1))*12)**alpha_AgAg) !eV
         E2=epsilon_AgAg*((exp(-2*q_AgAg*(r2/r0_AgAg-1))*12)**alpha_AgAg) !eV
         E3=epsilon_AgAg*((exp(-2*q_AgAg*(r3/r0_AgAg-1))*12)**alpha_AgAg) !eV
         E4=epsilon_AgAg*((exp(-2*q_AgAg*(r4/r0_AgAg-1))*12)**alpha_AgAg) !eV
         WRITE(200,95) E1,E2,E3,E4
         95 FORMAT(1X, E20.12, E20.12, E20.12, E20.12)
         r1=r1+(4)*dx
         r2=r2+(4)*dx
         r3=r3+(4)*dx
         r4=r4+(4)*dx
      
      END DO
      close(200)
   
   
      !Ag-Pt
   
      OPEN(UNIT=300,FILE="AgPt_f.TXT") 
      WRITE(6,*) "AgPt_f..."
      xstart=338.859_dp ! corresponding to r=0.001 Angstrong
      xstop=0.0907_dp !corresponding to cut off radius for Ag-Ag r_c=4.08 Angstrong
      r1=xstart
      dx=(xstop-xstart)/REAL(npts,dp)
      delr=(1/dx)
      r2=xstart+dx
      r3=xstart+2*dx
      r4=xstart+3*dx
      WRITE(300,*) npts, xstart,xstop,delr 
      
      DO j=1,npts,4 
         E1=epsilon_AgAg*((exp(-2*q_AgAg*(r1/r0_AgPt-1))*12)**alpha_AgPt) !eV
         E2=epsilon_AgAg*((exp(-2*q_AgAg*(r2/r0_AgPt-1))*12)**alpha_AgPt) !eV
         E3=epsilon_AgAg*((exp(-2*q_AgAg*(r3/r0_AgPt-1))*12)**alpha_AgPt) !eV
         E4=epsilon_AgAg*((exp(-2*q_AgAg*(r4/r0_AgPt-1))*12)**alpha_AgPt) !eV
         WRITE(300,100) E1,E2,E3,E4
         100 FORMAT(1X, E20.12, E20.12, E20.12, E20.12)
         r1=r1+(4)*dx
         r2=r2+(4)*dx
         r3=r3+(4)*dx
         r4=r4+(4)*dx
     
      END DO
      close(300)
   END SUBROUTINE GenerateTBSMATables
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
END MODULE TBSMAGenerate
