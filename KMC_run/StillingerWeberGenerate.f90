MODULE StillingerWeberGenerate
   USE VARIABLE_TYPE
   IMPLICIT NONE
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GenerateStillingerWeberSubrtTables()
      IMPLICIT NONE
      CHARACTER(len=100):: Filename
      
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      !Create potential file for Si_Si interaction
      WRITE(UNIT=6,FMT='("SWInitialize>> Creating SW table")',ADVANCE="NO")
      Filename="SiSi.SW.txt"
      CALL SWSetup(npts=3000,A=7.049556277_dp,B=0.6022245584_dp,p=4._dp,q=0._dp, &
         eps=2.1678_dp,lam1=21.0_dp,lam2=21.0_dp,gamma=1.2_dp,ac=1.8_dp,theta0=0.333333333333_dp, &
         sigma=2.0951_dp,xstart=1.2_dp,xstop=3.74_dp,Filename=Filename)  !xstop=3.77118_dp
         !parameter values have taken from reference 
         !F.H.Stillinger Phy.Rev.B,31, 5262
         !Mohamed Laradji, D.P.Landau Phy.Rev.B,51, 4894
         !Zi Jian Phy.Rev.B,vol 41
      
      ! create potential file for Ge_Ge interaction
      WRITE(UNIT=6,FMT='("SWInitialize>> Creating SW table")',ADVANCE="NO")
      Filename="GeGe.SW.txt"
      CALL SWSetup(npts=3000,A=7.049556277_dp,B=0.6022245584_dp,p=4.0_dp,q=0.0_dp, &
         eps=1.93_dp,lam1=31.0_dp,lam2=31.0_dp,gamma=1.2_dp,ac=1.8_dp,theta0=0.333333333333_dp, &
         sigma=2.181_dp,xstart=1.2_dp,xstop=3.90_dp,Filename=Filename)  !xstop=3.9258_dp
         !parameter values have taken from reference 
         !F.H.Stillinger Phy.Rev.B,31, 5262
         !Mohamed Laradji, D.P.Landau Phy.Rev.B,51, 4894
         !Zi Jian Phy.Rev.B,vol 41
      
      ! create potential file for Si_Ge interaction
      WRITE(UNIT=6,FMT='("SWInitialize>> Creating SW table")',ADVANCE="NO")
      Filename="SiGe.SW.txt"
      CALL SWSetup(npts=3000,A=7.049556277_dp,B=0.6022245584_dp,p=4.0_dp,q=0.0_dp, &
         eps=2.0427_dp,lam1=21.0_dp,lam2=31.0_dp,gamma=1.2_dp,ac=1.8_dp,theta0=0.333333333333_dp, &
         sigma=2.1353_dp,xstart=1.2_dp,xstop=3.82_dp,Filename=Filename) !xstop=3.84354_dp
         !parameter values have taken from reference 
         !F.H.Stillinger Phy.Rev.B,31, 5262
         !Mohamed Laradji, D.P.Landau Phy.Rev.B,51, 4894
         !Zi Jian Phy.Rev.B,vol 41
   END SUBROUTINE GenerateStillingerWeberSubrtTables
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SWSetup(npts,A,B,p,q,eps,lam1,lam2,gamma,ac,theta0,sigma,xstart,xstop,Filename)
      IMPLICIT NONE
      REAL(dp) :: xstart,xstop,dx,E1,E2,E3,E4,r1,r2,r3,r4,delr,A,B,p,q,eps,lam,gamma, &
         ac,theta0,sigma,lam1,lam2,epssqrt,F3a,F3b,F3c,F3d
      INTEGER :: npts,j
      CHARACTER(len=100):: Filename
   
      OPEN(UNIT=200,FILE=Filename)
      r1=xstart  !angstrom
      dx=(xstop-xstart)/REAL(npts,dp)
      delr=1./dx
      r2=xstart+dx
      r3=xstart+2.*dx
      r4=xstart+3.*dx
      !XXXXXXXXXX this is for pair potential v2 XXXXXXXXX
      WRITE(200,*) npts,xstart,xstop,delr 
      WRITE(UNIT=200,fmt=*) "comment line : Instructions for Pair Potential using formula v2=", &
      "E=eps*A*(B*(r/sigma)**(-p)-(r/sigma)**(-q))*exp(1/(r/sigma-ac) Ref:F.H.Stillinger Phy.Rev.B,31, 5262"
       
      WRITE(6,*) "Creating file ",TRIM(Filename)
      WRITE(6,*) "Parameter values"
      WRITE(6,*) "----------------"
      WRITE(UNIT=6,FMT='("A: ",F15.5)') A
      WRITE(UNIT=6,FMT='("B: ",F15.5)') B
      WRITE(UNIT=6,FMT='("p: ",F15.5)') p
      WRITE(UNIT=6,FMT='("q: ",F15.5)') q
      WRITE(UNIT=6,FMT='("eps: ",F15.5)') eps
      WRITE(UNIT=6,FMT='("lam1,lam2 ",2F15.5)') lam1,lam2
      WRITE(UNIT=6,FMT='("gamma: ",F15.5)') gamma
      WRITE(UNIT=6,FMT='("ac: ",F15.5)') ac
      WRITE(UNIT=6,FMT='("sigma: ",F15.5)') sigma
      WRITE(UNIT=6,FMT='("xstart,xstop: ",F15.5)') xstart,xstop
      DO j=1,3000,4
      ! equation taken from reference F.H.Stillinger Phy.Rev.B,31, 5262
         E1=v2(eps,A,B,r1/sigma,p,q,ac) !eV
         !eps*A*(B*(r1/sigma)**(-p)-(r1/sigma)**(-q))*exp(1._dp/((r1/sigma)-ac))
         E2=v2(eps,A,B,r2/sigma,p,q,ac) !ev 
         E3=v2(eps,A,B,r3/sigma,p,q,ac) !ev 
         E4=v2(eps,A,B,r4/sigma,p,q,ac) !ev
         WRITE(UNIT=200,FMT='(1X,4E20.12)') E1,E2,E3,E4
         
         r1=r1+4.*dx
         r2=r2+4.*dx
         r3=r3+4.*dx
         r4=r4+4.*dx
      END DO
      !XXXXXXXXXXXX this is for three body potential v3 XXXXXXXXXXXXXXXX
       
      r1=xstart  !angstrom
      dx=(xstop-xstart)/REAL(npts,dp)
      delr=1./dx
      r2=xstart+dx
      r3=xstart+2.*dx
      r4=xstart+3.*dx
      WRITE(200,*) "  "
      WRITE(200,*) npts,xstart,xstop,delr 
      WRITE(UNIT=200,fmt=*) "Comment line: Instructions for Three body potential", &
         "v3=eps*(sum of h)and h=epssqrtij*epssqrtik*lamij*lamik*F3ij*F3ik*costerm Ref:D.P.Landau Phy.Rev.B,51, 4894"
      DO j=1,3000,4
         F3a=exp(gamma/((r1/sigma)-ac))       !gives values of F3 for pairs ij or ik Ref:D.P.Landau Phy.Rev.B,51
         F3b=exp(gamma/((r2/sigma)-ac))
         F3c=exp(gamma/((r3/sigma)-ac))
         F3d=exp(gamma/((r4/sigma)-ac))
         WRITE(UNIT=200,FMT='(1X,4E20.12)') F3a,F3b,F3c,F3d
         r1=r1+4.*dx
         r2=r2+4.*dx
         r3=r3+4.*dx
         r4=r4+4.*dx
      END DO
       
      WRITE(200,*) " "
      WRITE(200,*) "Comment line : Instructions for epssqrt and lam for pairs ij or ik, Ref:D.P.Landau Phy.rev.B,51, 4894"
      lam=(lam1*lam2)**0.25
      epssqrt=eps**0.5
      WRITE(UNIT=200,FMT='(1X,E20.12)') epssqrt
      WRITE(UNIT=200,FMT='(1X,E20.12)') lam
    
      CLOSE(200)   
      WRITE(6,*) " ... ",TRIM(Filename)," created."
   END SUBROUTINE SWSetup
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION v2(eps,A,B,r,p,q,ac) !v2 term of SW
      IMPLICIT NONE
      REAL(dp) :: eps,A,B,r,p,q,ac,v2
      
      IF (r<ac) THEN
         v2=eps*A*(B*r**(-p)-r**(-q))*exp(1._dp/(r-ac))
      ELSEIF (r>=ac) THEN
         v2=0._dp
      ELSE
         WRITE(6,*) "Err>> Interatomic spacing is negative in v2 function"
         STOP
      END IF
   END FUNCTION v2
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE StillingerWeberGenerate
