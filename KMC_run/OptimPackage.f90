MODULE OptimPackage
   !contains codes for steepest-descent,conjugate-gradient,dfp,lbfgs,simplex,sa
   !Usage
   USE VARIABLE_TYPE
   USE nrprocedures
   IMPLICIT NONE
   
   INTERFACE ConjugateGradient
      MODULE PROCEDURE  ConjugateGradientdp,ConjugateGradientsp
   END INTERFACE ConjugateGradient
   INTERFACE DFP
      MODULE PROCEDURE DFPdp,DFPsp
   END INTERFACE DFP
   INTERFACE Lbfgs
      MODULE PROCEDURE Lbfgsdp,Lbfgssp
   END INTERFACE Lbfgs
   INTERFACE SteepestDescent
      MODULE PROCEDURE SteepestDescentdp,SteepestDescentsp
   END INTERFACE SteepestDescent
   INTERFACE brent
      MODULE PROCEDURE brentdp,brentsp
   END INTERFACE brent
   INTERFACE dbrent
      MODULE PROCEDURE dbrentdp,dbrentsp
   END INTERFACE dbrent
   INTERFACE linmin
      MODULE PROCEDURE linmindp,linminsp
   END INTERFACE linmin
   INTERFACE lnsrch
      MODULE PROCEDURE lnsrchdp,lnsrchsp
   END INTERFACE lnsrch
   INTERFACE mnbrak
      MODULE PROCEDURE mnbrakdp,mnbraksp
   END INTERFACE mnbrak
   INTERFACE mnbrak1d
      MODULE PROCEDURE mnbrak1d_dp,mnbrak1d_sp
   END INTERFACE mnbrak1d
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Newton1D(x,tol,iter,func,dfunc,errorstatus)
      !my implementation of the Newton method for 1D
      !finds x where func=0, i.e., this is nonlinear eq solver
      IMPLICIT NONE
      REAL(dp), INTENT(INOUT) :: x
      REAL(dp), INTENT(IN) :: tol
      INTEGER, INTENT(IN) :: iter
      INTEGER :: it,errorstatus
      REAL(dp) :: f,g,fold
      INTERFACE
         FUNCTION func(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), INTENT(IN) :: p
            REAL(dp) :: func
         END FUNCTION func
         FUNCTION dfunc(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), INTENT(IN) :: p
            REAL(dp) :: dfunc
         END FUNCTION dfunc
      END INTERFACE
      
      DO it=1,iter
         f=func(x,errorstatus)
         IF (errorstatus/=0) RETURN
         IF (it>4 .AND. ABS(f)<tol) RETURN
         fold=f
         g=dfunc(x,errorstatus)
         IF (errorstatus/=0) RETURN
         write(*,*) "Newton x,f,g:",x,f,g
         x=x-0.5_dp*f/g
      END DO
   END SUBROUTINE Newton1D
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SteepestDescentdp(p,ftol,gtol,xtol,iter,fret,func,dfunc,dx,ITMAX,IsPrint1,errorstatus)
      !Written by Abhijit based on ConjugateGradient sbrtn
      !p(1:n)=Initial guess
      !gtol=Convergence requirement on zeroing the gradient
      !iter=Return # iterations performed
      !fret=Min value of fn found
      !func=User function
      !dfunc=User gradient function
      IMPLICIT NONE
      INTEGER :: errorstatus
      INTEGER(I4B), INTENT(OUT) :: iter
      REAL(dp), INTENT(IN) :: ftol,gtol,xtol
      REAL(dp), INTENT(OUT) :: fret
      REAL(dp), DIMENSION(:), INTENT(INOUT) :: p
      INTERFACE
         FUNCTION func(p,errorstatus) !energy
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), DIMENSION(:), INTENT(IN) :: p
            REAL(dp) :: func
         END FUNCTION func

         FUNCTION dfunc(p,errorstatus) !
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), DIMENSION(:), INTENT(IN) :: p
            REAL(dp), DIMENSION(size(p)) :: dfunc
         END FUNCTION dfunc
      END INTERFACE
      !INTEGER(I4B), PARAMETER :: ITMAX=100
      INTEGER :: ITMAX
      REAL(dp), PARAMETER :: STPMX=100.0_dp,EPS=epsilon(p),TOLX=4.0_dp*EPS
      INTEGER(I4B) :: its
      REAL(dp) :: gg,fp
      REAL(dp) :: ggmax,dpmax,dV,fpold
      REAL(dp), DIMENSION(size(p)) :: h,xi,pold
      REAL(dp), OPTIONAL :: dx
      LOGICAL, OPTIONAL :: IsPrint1
      LOGICAL :: IsPrint
      
      IF (PRESENT(IsPrint1)) THEN
         IsPrint=IsPrint1
      ELSE
         IsPrint=.FALSE.
      END IF
      
      fp=func(p,errorstatus)
      IF (errorstatus/=0) RETURN
      xi=-dfunc(p,errorstatus) !get the energy min direction
      IF (errorstatus/=0) RETURN
      DO its=1,ITMAX
         iter=its
         pold=p
         fpold=fp
         IF (PRESENT(dx)) THEN
            CALL linmin(p,xi,fret,func,dx,errorstatus) !dx gives the distance for line search
            IF (errorstatus/=0) RETURN
         ELSE
            CALL linmin(p,xi,fret,func,errorstatus=errorstatus)
            IF (errorstatus/=0) RETURN
         END IF
         IF (2.0_dp*ABS(fret-fp) <= ftol*(ABS(fret)+ABS(fp)+EPS)) THEN
            !WRITE(*,*) "Fn value change too small"
            RETURN
         END IF
         IF (errorstatus/=0) RETURN
         dV=fp-fpold
         xi=-dfunc(p,errorstatus)
         IF (errorstatus/=0) RETURN
         gg=DOT_PRODUCT(xi,xi)
         IF (gg==0.0) RETURN
         dpmax=MAXVAL(ABS(p-pold))
         IF (IsPrint) WRITE(UNIT=*,FMT='("...",I3," fp: ",ES14.6, &
            " dpmax:",ES14.6, &
            " gg:",ES14.6," dV:",ES14.6)') iter,fp,dpmax,gg,dV
         IF (dpmax<xtol .AND. gg<gtol .AND. ABS(dV)<ftol) RETURN
      END DO
      IF (IsPrint) WRITE(*,*) "$Err>>Steepest descent: maximum iterations exceeded"
   END SUBROUTINE SteepestDescentdp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SteepestDescentsp(p,ftol,gtol,xtol,iter,fret,func,dfunc,dx,ITMAX,IsPrint1,errorstatus)
      !Written by Abhijit based on ConjugateGradient sbrtn
      !p(1:n)=Initial guess
      !gtol=Convergence requirement on zeroing the gradient
      !iter=Return # iterations performed
      !fret=Min value of fn found
      !func=User function
      !dfunc=User gradient function
      IMPLICIT NONE
      INTEGER :: errorstatus
      INTEGER(I4B), INTENT(OUT) :: iter
      REAL(sp), INTENT(IN) :: ftol,gtol,xtol
      REAL(sp), INTENT(OUT) :: fret
      REAL(sp), DIMENSION(:), INTENT(INOUT) :: p
      INTERFACE
         FUNCTION func(p,errorstatus) !energy
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), DIMENSION(:), INTENT(IN) :: p
            REAL(sp) :: func
         END FUNCTION func

         FUNCTION dfunc(p,errorstatus) !
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), DIMENSION(:), INTENT(IN) :: p
            REAL(sp), DIMENSION(size(p)) :: dfunc
         END FUNCTION dfunc
      END INTERFACE
      !INTEGER(I4B), PARAMETER :: ITMAX=100
      INTEGER :: ITMAX
      REAL(sp), PARAMETER :: STPMX=100.0_sp,EPS=epsilon(p),TOLX=4.0_sp*EPS
      INTEGER(I4B) :: its
      REAL(sp) :: gg,fp
      REAL(sp) :: ggmax,dpmax,dV,fpold
      REAL(sp), DIMENSION(size(p)) :: h,xi,pold
      REAL(sp), OPTIONAL :: dx
      LOGICAL, OPTIONAL :: IsPrint1
      LOGICAL :: IsPrint
      
      IF (PRESENT(IsPrint1)) THEN
         IsPrint=IsPrint1
      ELSE
         IsPrint=.FALSE.
      END IF
      
      fp=func(p,errorstatus)
      IF (errorstatus/=0) RETURN
      xi=-dfunc(p,errorstatus) !get the energy min direction
      IF (errorstatus/=0) RETURN
      DO its=1,ITMAX
         iter=its
         pold=p
         fpold=fp
         IF (PRESENT(dx)) THEN
            CALL linmin(p,xi,fret,func,dx,errorstatus) !dx gives the distance for line search
            IF (errorstatus/=0) RETURN
         ELSE
            CALL linmin(p,xi,fret,func,errorstatus=errorstatus)
            IF (errorstatus/=0) RETURN
         END IF
         IF (2.0_sp*ABS(fret-fp) <= ftol*(ABS(fret)+ABS(fp)+EPS)) THEN
            !WRITE(*,*) "Fn value change too small"
            RETURN
         END IF
         IF (errorstatus/=0) RETURN
         dV=fp-fpold
         xi=-dfunc(p,errorstatus)
         IF (errorstatus/=0) RETURN
         gg=DOT_PRODUCT(xi,xi)
         IF (gg==0.0) RETURN
         dpmax=MAXVAL(ABS(p-pold))
         IF (IsPrint) WRITE(UNIT=*,FMT='("...",I3," fp: ",ES14.6, &
            " dpmax:",ES14.6, &
            " gg:",ES14.6," dV:",ES14.6)') iter,fp,dpmax,gg,dV
         IF (dpmax<xtol .AND. gg<gtol .AND. ABS(dV)<ftol) RETURN
      END DO
      WRITE(*,*) "$Err>>Steepest descent: maximum iterations exceeded"
   END SUBROUTINE SteepestDescentsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DFPdp(p,ftol,gtol,xtol,iter,fret,func,dfunc,IsPrint1,errorstatus)
      !p(1:n)=Initial guess
      !gtol=Convergence requirement on zeroing the gradient
      !iter=Return # iterations performed
      !fret=Min value of fn found
      !func=User function
      !dfunc=User gradient function
      IMPLICIT NONE
      INTEGER(I4B), INTENT(OUT) :: iter
      INTEGER :: errorstatus
      REAL(dp), INTENT(IN) :: ftol,gtol,xtol
      REAL(dp), INTENT(OUT) :: fret
      REAL(dp), DIMENSION(:), INTENT(INOUT) :: p
      INTERFACE
         FUNCTION func(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), DIMENSION(:), INTENT(IN) :: p
            REAL(dp) :: func
         END FUNCTION func

         FUNCTION dfunc(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), DIMENSION(:), INTENT(IN) :: p
            REAL(dp), DIMENSION(size(p)) :: dfunc
         END FUNCTION dfunc
      END INTERFACE
      INTEGER(I4B), PARAMETER :: ITMAX=200
      REAL(dp), PARAMETER :: STPMX=100.0_dp,EPS=epsilon(p),TOLX=4.0_dp*EPS
      INTEGER(I4B) :: its
      LOGICAL :: check
      REAL(dp) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
      REAL(dp), DIMENSION(size(p)) :: dg,g,hdg,pnew,xi
      REAL(dp), DIMENSION(size(p),size(p)) :: hessin
      REAL(dp) :: dpmax,gg,dV,fpold
      LOGICAL, OPTIONAL :: IsPrint1
      LOGICAL :: IsPrint
      
      IF (PRESENT(IsPrint1)) THEN
         IsPrint=IsPrint1
      ELSE
         IsPrint=.FALSE.
      END IF
      
      fp=func(p,errorstatus)
      IF (errorstatus/=0) RETURN
      g=dfunc(p,errorstatus)
      IF (errorstatus/=0) RETURN

      CALL unit_matrix(hessin)
      xi=-g
      stpmax=STPMX*MAX(vABS(p),real(size(p),sp))
      DO its=1,ITMAX
         iter=its
         fpold=fp
         CALL lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,func,errorstatus)
         IF (errorstatus/=0) RETURN
         fp=fret
         dV=fp-fpold
         xi=pnew-p
         dpmax=MAXVAL(ABS(xi))
         p=pnew
         IF (MAXVAL(ABS(xi)/MAX(ABS(p),1.0_dp)) < TOLX) THEN
            WRITE(*,*) "Insignificant change in x."
            RETURN
         END IF
         dg=g
         g=dfunc(p,errorstatus)
         IF (errorstatus/=0) RETURN
         den=MAX(fret,1.0_dp)
         IF (MAXVAL(ABS(g)*MAX(ABS(p),1.0_dp)/den) < gtol) RETURN
         dg=g-dg
         hdg=MATMUL(hessin,dg)
         fac=DOT_PRODUCT(dg,xi)
         fae=DOT_PRODUCT(dg,hdg)
         sumdg=DOT_PRODUCT(dg,dg)
         sumxi=DOT_PRODUCT(xi,xi)
         IF (fac**2 > EPS*sumdg*sumxi) THEN
            fac=1.0_dp/fac
            fad=1.0_dp/fae
            dg=fac*xi-fad*hdg
            hessin=hessin+fac*outerprod(xi,xi)-&
               fad*outerprod(hdg,hdg)+fae*outerprod(dg,dg)
         END IF
         xi=-MATMUL(hessin,g)
         IF (IsPrint) WRITE(UNIT=*,FMT='("...",I3," fp: ",ES14.6, &
            " dpmax:",ES14.6, &
            " gg:",ES14.6," dV:",ES14.6)') iter,fp,dpmax,gg,dV
         IF (dpmax<xtol .AND. gg<gtol .AND. ABS(dV)<ftol) RETURN
      END DO
      WRITE(*,*) "Err>>DFP: too many iterations"
   END SUBROUTINE DFPdp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE DFPsp(p,ftol,gtol,xtol,iter,fret,func,dfunc,IsPrint1,errorstatus)
      !p(1:n)=Initial guess
      !gtol=Convergence requirement on zeroing the gradient
      !iter=Return # iterations performed
      !fret=Min value of fn found
      !func=User function
      !dfunc=User gradient function
      IMPLICIT NONE
      INTEGER(I4B), INTENT(OUT) :: iter
      INTEGER :: errorstatus
      REAL(sp), INTENT(IN) :: ftol,gtol,xtol
      REAL(sp), INTENT(OUT) :: fret
      REAL(sp), DIMENSION(:), INTENT(INOUT) :: p
      INTERFACE
         FUNCTION func(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), DIMENSION(:), INTENT(IN) :: p
            REAL(sp) :: func
         END FUNCTION func

         FUNCTION dfunc(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), DIMENSION(:), INTENT(IN) :: p
            REAL(sp), DIMENSION(size(p)) :: dfunc
         END FUNCTION dfunc
      END INTERFACE
      INTEGER(I4B), PARAMETER :: ITMAX=200
      REAL(sp), PARAMETER :: STPMX=100.0_sp,EPS=epsilon(p),TOLX=4.0_sp*EPS
      INTEGER(I4B) :: its
      LOGICAL :: check
      REAL(sp) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
      REAL(sp), DIMENSION(size(p)) :: dg,g,hdg,pnew,xi
      REAL(sp), DIMENSION(size(p),size(p)) :: hessin
      REAL(sp) :: dpmax,gg,dV,fpold
      LOGICAL, OPTIONAL :: IsPrint1
      LOGICAL :: IsPrint
      
      IF (PRESENT(IsPrint1)) THEN
         IsPrint=IsPrint1
      ELSE
         IsPrint=.FALSE.
      END IF
      
      fp=func(p,errorstatus)
      IF (errorstatus/=0) RETURN
      g=dfunc(p,errorstatus)
      IF (errorstatus/=0) RETURN

      CALL unit_matrixsp(hessin)
      xi=-g
      stpmax=STPMX*MAX(vABSsp(p),real(size(p),sp))
      DO its=1,ITMAX
         iter=its
         fpold=fp
         CALL lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,func,errorstatus)
         IF (errorstatus/=0) RETURN
         fp=fret
         dV=fp-fpold
         xi=pnew-p
         dpmax=MAXVAL(ABS(xi))
         p=pnew
         IF (MAXVAL(ABS(xi)/MAX(ABS(p),1.0_sp)) < TOLX) THEN
            WRITE(*,*) "Insignificant change in x."
            RETURN
         END IF
         dg=g
         g=dfunc(p,errorstatus)
         IF (errorstatus/=0) RETURN
         den=MAX(fret,1.0_sp)
         IF (MAXVAL(ABS(g)*MAX(ABS(p),1.0_sp)/den) < gtol) RETURN
         dg=g-dg
         hdg=MATMUL(hessin,dg)
         fac=DOT_PRODUCT(dg,xi)
         fae=DOT_PRODUCT(dg,hdg)
         sumdg=DOT_PRODUCT(dg,dg)
         sumxi=DOT_PRODUCT(xi,xi)
         IF (fac**2 > EPS*sumdg*sumxi) THEN
            fac=1.0_sp/fac
            fad=1.0_sp/fae
            dg=fac*xi-fad*hdg
            hessin=hessin+fac*outerprodsp(xi,xi)-&
               fad*outerprodsp(hdg,hdg)+fae*outerprodsp(dg,dg)
         END IF
         xi=-MATMUL(hessin,g)
         IF (IsPrint) WRITE(UNIT=*,FMT='("...",I3," fp: ",ES14.6, &
            " dpmax:",ES14.6, &
            " gg:",ES14.6," dV:",ES14.6)') iter,fp,dpmax,gg,dV
         IF (dpmax<xtol .AND. gg<gtol .AND. ABS(dV)<ftol) RETURN
      END DO
      WRITE(*,*) "Err>>DFP: too many iterations"
   END SUBROUTINE DFPsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ConjugateGradientdp(p,ftol,gtol,xtol,iter,fret,func,dfunc,dx,ITMAX,IsPrint1,errorstatus)
      !p(1:n)=Initial guess
      !ftol=Convergence requirement on zeroing the gradient
      !iter=Return # iterations performed
      !fret=Min value of fn found
      !func=User function
      !dfunc=User gradient function
      IMPLICIT NONE
      INTEGER :: errorstatus
      INTEGER(I4B), INTENT(OUT) :: iter
      REAL(dp), INTENT(IN) :: ftol,gtol,xtol
      REAL(dp), INTENT(OUT) :: fret
      REAL(dp), DIMENSION(:), INTENT(INOUT) :: p
      INTERFACE
         FUNCTION func(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), DIMENSION(:), INTENT(IN) :: p
            REAL(dp) :: func
         END FUNCTION func
   
         FUNCTION dfunc(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), DIMENSION(:), INTENT(IN) :: p
            REAL(dp), DIMENSION(size(p)) :: dfunc
         END FUNCTION dfunc
      END INTERFACE
      !INTEGER(I4B), PARAMETER :: ITMAX=200
      INTEGER :: ITMAX 
      REAL(dp), PARAMETER :: EPS=1.0e-10_dp
      INTEGER(I4B) :: its
      REAL(dp) :: dgg,fp,gam,gg
      REAL(dp) :: ggmax,dpmax,dV,fpold
      REAL(dp), DIMENSION(size(p)) :: g,h,xi,pold
      REAL(dp), OPTIONAL :: dx
      LOGICAL, OPTIONAL :: IsPrint1
      LOGICAL :: IsPrint
      
      IF (PRESENT(IsPrint1)) THEN
         IsPrint=IsPrint1
      ELSE
         IsPrint=.FALSE.
      END IF
      
      fp=func(p,errorstatus)
      IF (errorstatus/=0) RETURN
      xi=dfunc(p,errorstatus)
      IF (errorstatus/=0) RETURN
      g=-xi
      h=g
      xi=h
      DO its=1,ITMAX
         iter=its
         pold=p
         fpold=fp
         IF (PRESENT(dx)) THEN
            CALL linmin(p,xi,fret,func,dx,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE
            CALL linmin(p,xi,fret,func,errorstatus=errorstatus)
            IF (errorstatus/=0) RETURN
         END IF
         IF (2.0_dp*ABS(fret-fp) <= ftol*(ABS(fret)+ABS(fp)+EPS)) THEN
            !WRITE(*,*) "Fn value change too small"
            RETURN
         END IF
         fp=func(p,errorstatus)
         IF (errorstatus/=0) RETURN
         dV=fp-fpold
         xi=dfunc(p,errorstatus)
         IF (errorstatus/=0) RETURN
         gg=DOT_PRODUCT(g,g)
   !     dgg=DOT_PRODUCT(xi,xi)
         dgg=DOT_PRODUCT(xi+g,xi)
         IF (gg == 0.0) RETURN
         gam=dgg/gg
         g=-xi
         h=g+gam*h
         xi=h
         dpmax=MAXVAL(ABS(p-pold))
         IF (IsPrint .AND. MOD(its,10)==1) WRITE(UNIT=*,FMT='("...",I3," fp: ",ES14.6, &
            " dpmax:",ES14.6, &
            " gg:",ES14.6," dV:",ES14.6)') iter,fp,dpmax,gg,dV
         IF (dpmax<xtol .AND. gg<gtol .AND. ABS(dV)<ftol) RETURN
      END DO
      WRITE(*,*) "Err>>ConjugateGradient: maximum iterations exceeded"
   END SUBROUTINE ConjugateGradientdp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ConjugateGradientsp(p,ftol,gtol,xtol,iter,fret,func,dfunc,dx,ITMAX,IsPrint1,errorstatus)
      !p(1:n)=Initial guess
      !ftol=Convergence requirement on zeroing the gradient
      !iter=Return # iterations performed
      !fret=Min value of fn found
      !func=User function
      !dfunc=User gradient function
      IMPLICIT NONE
      INTEGER :: errorstatus
      INTEGER(I4B), INTENT(OUT) :: iter
      REAL(sp), INTENT(IN) :: ftol,gtol,xtol
      REAL(sp), INTENT(OUT) :: fret
      REAL(sp), DIMENSION(:), INTENT(INOUT) :: p
      INTERFACE
         FUNCTION func(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), DIMENSION(:), INTENT(IN) :: p
            REAL(sp) :: func
         END FUNCTION func
   
         FUNCTION dfunc(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), DIMENSION(:), INTENT(IN) :: p
            REAL(sp), DIMENSION(size(p)) :: dfunc
         END FUNCTION dfunc
      END INTERFACE
      !INTEGER(I4B), PARAMETER :: ITMAX=200
      INTEGER :: ITMAX
      REAL(sp), PARAMETER :: EPS=1.0e-10_sp
      INTEGER(I4B) :: its
      REAL(sp) :: dgg,fp,gam,gg
      REAL(sp) :: ggmax,dpmax,dV,fpold
      REAL(sp), DIMENSION(size(p)) :: g,h,xi,pold
      REAL(sp), OPTIONAL :: dx
      LOGICAL, OPTIONAL :: IsPrint1
      LOGICAL :: IsPrint
      
      IF (PRESENT(IsPrint1)) THEN
         IsPrint=IsPrint1
      ELSE
         IsPrint=.FALSE.
      END IF
      
      fp=func(p,errorstatus)
      IF (errorstatus/=0) RETURN
      xi=dfunc(p,errorstatus)
      IF (errorstatus/=0) RETURN
      g=-xi
      h=g
      xi=h
      DO its=1,ITMAX
         iter=its
         pold=p
         fpold=fp
         IF (PRESENT(dx)) THEN
            CALL linmin(p,xi,fret,func,dx,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE
            CALL linmin(p,xi,fret,func,errorstatus=errorstatus)
            IF (errorstatus/=0) RETURN
         END IF
         IF (2.0_sp*ABS(fret-fp) <= ftol*(ABS(fret)+ABS(fp)+EPS)) THEN
            !WRITE(*,*) "Fn value change too small"
            RETURN
         END IF
         fp=func(p,errorstatus)
         IF (errorstatus/=0) RETURN
         dV=fp-fpold
         xi=dfunc(p,errorstatus)
         IF (errorstatus/=0) RETURN
         gg=DOT_PRODUCT(g,g)
   !     dgg=DOT_PRODUCT(xi,xi)
         dgg=DOT_PRODUCT(xi+g,xi)
         IF (gg == 0.0) RETURN
         gam=dgg/gg
         g=-xi
         h=g+gam*h
         xi=h
         dpmax=MAXVAL(ABS(p-pold))
         IF (IsPrint .AND. MOD(its,10)==1) WRITE(UNIT=*,FMT='("...",I3," fp: ",ES14.6, &
            " dpmax:",ES14.6, &
            " gg:",ES14.6," dV:",ES14.6)') iter,fp,dpmax,gg,dV
         IF (dpmax<xtol .AND. gg<gtol .AND. ABS(dV)<ftol) RETURN
      END DO
      WRITE(*,*) "Err>>ConjugateGradient: maximum iterations exceeded"
   END SUBROUTINE ConjugateGradientsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE linmindp(p,xi,fret,func,dxri,errorstatus)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(dp), INTENT(OUT) :: fret
      REAL(dp), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
      REAL(dp), PARAMETER :: TOL=1.0e-4_dp
      REAL(dp) :: ax,bx,fa,fb,fx,xmin,xx,dxr
      REAL(dp), OPTIONAL :: dxri
      INTERFACE
         FUNCTION func(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), DIMENSION(:), INTENT(IN) :: p
            REAL(dp) :: func
         END FUNCTION func
      END INTERFACE

      IF (PRESENT(dxri)) THEN
         dxr=dxri
      ELSE
         dxr=1.0_dp
      END IF
      
      ax=0.0
      xx=1.0*dxr
      CALL mnbrak(ax,xx,bx,fa,fx,fb,func,p,xi,errorstatus)
      IF (errorstatus/=0) RETURN
      fret=brent(ax,xx,bx,func,p,xi,TOL,xmin,errorstatus)
      IF (errorstatus/=0) RETURN
      !fret=dbrent(ax,xx,bx,f1dim,TOL,xmin) !NR recommends using dbrent
      xi=xmin*xi
      p=p+xi
   END SUBROUTINE linmindp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE linminsp(p,xi,fret,func,dxri,errorstatus)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(sp), INTENT(OUT) :: fret
      REAL(sp), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
      REAL(sp), PARAMETER :: TOL=1.0e-4_sp
      REAL(sp) :: ax,bx,fa,fb,fx,xmin,xx,dxr
      REAL(sp), OPTIONAL :: dxri
      INTERFACE
         FUNCTION func(p,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), DIMENSION(:), INTENT(IN) :: p
            REAL(sp) :: func
         END FUNCTION func
      END INTERFACE

      IF (PRESENT(dxri)) THEN
         dxr=dxri
      ELSE
         dxr=1.0_sp
      END IF
      
      ax=0.0
      xx=1.0*dxr
      CALL mnbrak(ax,xx,bx,fa,fx,fb,func,p,xi,errorstatus)
      IF (errorstatus/=0) RETURN
      fret=brent(ax,xx,bx,func,p,xi,TOL,xmin,errorstatus)
      IF (errorstatus/=0) RETURN
      !fret=dbrent(ax,xx,bx,f1dim,TOL,xmin) !NR recommends using dbrent
      xi=xmin*xi
      p=p+xi
   END SUBROUTINE linminsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE mnbrak1d_dp(ax,bx,cx,fa,fb,fc,func,errorstatus) !to bracket with linmin
      !original version of NR
      !p is the initial vector position of the state (pls. check)
      !xi is the search direction (pls. check)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(dp), INTENT(INOUT) :: ax,bx
      REAL(dp), INTENT(OUT) :: cx,fa,fb,fc
      REAL(dp) :: fd
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), INTENT(IN) :: x
            REAL(dp) :: func
         END FUNCTION func
      END INTERFACE
      REAL(dp), PARAMETER :: GOLD=1.618034_dp,GLIMIT=100.0_dp,TINY=1.0e-20_dp !original one
      !REAL(dp), PARAMETER :: GOLD=1.2_dp,GLIMIT=100.0_dp,TINY=1.0e-20_dp !abhijit
      REAL(dp) :: fu,q,r,u,ulim

      fa=func(ax,errorstatus)
      IF (errorstatus/=0) RETURN
      fb=func(bx,errorstatus)
      IF (errorstatus/=0) RETURN

      IF (fb > fa) THEN
         CALL swap(ax,bx)
         CALL swap(fa,fb)
      END IF
      cx=bx+GOLD*(bx-ax)
      fc=func(cx,errorstatus)
      IF (errorstatus/=0) RETURN

      DO
         IF (fb < fc) RETURN
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*SIGN(MAX(ABS(q-r),TINY),q-r))
         ulim=bx+GLIMIT*(cx-bx)
         IF ((bx-u)*(u-cx) > 0.0) THEN
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
            IF (fu < fc) THEN
               ax=bx
               fa=fb
               bx=u
               fb=fu
               RETURN
            ELSE IF (fu > fb) THEN
               cx=u
               fc=fu
               RETURN
            END IF
            u=cx+GOLD*(cx-bx)
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE IF ((cx-u)*(u-ulim) > 0.0) THEN
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
            IF (fu < fc) THEN
               bx=cx
               cx=u
               u=cx+GOLD*(cx-bx)
               fd=func(u,errorstatus)
               IF (errorstatus/=0) RETURN
               CALL shft(fb,fc,fu,fd)
            END IF
         ELSE IF ((u-ulim)*(ulim-cx) >= 0.0) THEN
            u=ulim
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE
            u=cx+GOLD*(cx-bx)
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
         END IF
         CALL shft(ax,bx,cx,u)
         CALL shft(fa,fb,fc,fu)
      END DO
   END SUBROUTINE mnbrak1d_dp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE mnbrak1d_sp(ax,bx,cx,fa,fb,fc,func,errorstatus) !to bracket with linmin
      !original version of NR
      !p is the initial vector position of the state (pls. check)
      !xi is the search direction (pls. check)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(sp), INTENT(INOUT) :: ax,bx
      REAL(sp), INTENT(OUT) :: cx,fa,fb,fc
      REAL(sp) :: fd
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), INTENT(IN) :: x
            REAL(sp) :: func
         END FUNCTION func
      END INTERFACE
      REAL(sp), PARAMETER :: GOLD=1.618034_sp,GLIMIT=100.0_sp,TINY=1.0e-20_sp !original one
      !REAL(sp), PARAMETER :: GOLD=1.2_sp,GLIMIT=100.0_sp,TINY=1.0e-20_sp !abhijit
      REAL(sp) :: fu,q,r,u,ulim

      fa=func(ax,errorstatus)
      IF (errorstatus/=0) RETURN
      fb=func(bx,errorstatus)
      IF (errorstatus/=0) RETURN

      IF (fb > fa) THEN
         CALL swap(ax,bx)
         CALL swap(fa,fb)
      END IF
      cx=bx+GOLD*(bx-ax)
      fc=func(cx,errorstatus)
      IF (errorstatus/=0) RETURN

      DO
         IF (fb < fc) RETURN
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_sp*SIGN(MAX(ABS(q-r),TINY),q-r))
         ulim=bx+GLIMIT*(cx-bx)
         IF ((bx-u)*(u-cx) > 0.0) THEN
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
            IF (fu < fc) THEN
               ax=bx
               fa=fb
               bx=u
               fb=fu
               RETURN
            ELSE IF (fu > fb) THEN
               cx=u
               fc=fu
               RETURN
            END IF
            u=cx+GOLD*(cx-bx)
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE IF ((cx-u)*(u-ulim) > 0.0) THEN
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
            IF (fu < fc) THEN
               bx=cx
               cx=u
               u=cx+GOLD*(cx-bx)
               fd=func(u,errorstatus)
               IF (errorstatus/=0) RETURN
               CALL shftsp(fb,fc,fu,fd)
            END IF
         ELSE IF ((u-ulim)*(ulim-cx) >= 0.0) THEN
            u=ulim
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE
            u=cx+GOLD*(cx-bx)
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
         END IF
         CALL shftsp(ax,bx,cx,u)
         CALL shftsp(fa,fb,fc,fu)
      END DO
   END SUBROUTINE mnbrak1d_sp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE mnbrakdp(ax,bx,cx,fa,fb,fc,func,p,xi,errorstatus) !to bracket with linmin
      !modified version of NR
      !p is the initial vector position of the state (pls. check)
      !xi is the search direction (pls. check)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(dp), INTENT(INOUT) :: ax,bx
      REAL(dp), DIMENSION(:), INTENT(IN) :: p,xi
      REAL(dp), INTENT(OUT) :: cx,fa,fb,fc
      REAL(dp) :: fd
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), DIMENSION(:), INTENT(IN) :: x
            REAL(dp) :: func
         END FUNCTION func
      END INTERFACE
      REAL(dp), PARAMETER :: GOLD=1.618034_dp,GLIMIT=100.0_dp,TINY=1.0e-20_dp !original one
      !REAL(dp), PARAMETER :: GOLD=1.2_dp,GLIMIT=100.0_dp,TINY=1.0e-20_dp !abhijit
      REAL(dp) :: fu,q,r,u,ulim

      fa=func(p+ax*xi,errorstatus)
      IF (errorstatus/=0) RETURN
      fb=func(p+bx*xi,errorstatus)
      IF (errorstatus/=0) RETURN

      IF (fb > fa) THEN
         CALL swap(ax,bx)
         CALL swap(fa,fb)
      END IF
      cx=bx+GOLD*(bx-ax)
      fc=func(p+cx*xi,errorstatus)
      IF (errorstatus/=0) RETURN

      DO
         IF (fb < fc) RETURN
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*SIGN(MAX(ABS(q-r),TINY),q-r))
         ulim=bx+GLIMIT*(cx-bx)
         IF ((bx-u)*(u-cx) > 0.0) THEN
            fu=func(p+u*xi,errorstatus)
            IF (errorstatus/=0) RETURN
            IF (fu < fc) THEN
               ax=bx
               fa=fb
               bx=u
               fb=fu
               RETURN
            ELSE IF (fu > fb) THEN
               cx=u
               fc=fu
               RETURN
            END IF
            u=cx+GOLD*(cx-bx)
            fu=func(p+u*xi,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE IF ((cx-u)*(u-ulim) > 0.0) THEN
            fu=func(p+u*xi,errorstatus)
            IF (errorstatus/=0) RETURN
            IF (fu < fc) THEN
               bx=cx
               cx=u
               u=cx+GOLD*(cx-bx)
               fd=func(p+u*xi,errorstatus)
               IF (errorstatus/=0) RETURN
               CALL shft(fb,fc,fu,fd)
            END IF
         ELSE IF ((u-ulim)*(ulim-cx) >= 0.0) THEN
            u=ulim
            fu=func(p+u*xi,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE
            u=cx+GOLD*(cx-bx)
            fu=func(p+u*xi,errorstatus)
            IF (errorstatus/=0) RETURN
         END IF
         CALL shft(ax,bx,cx,u)
         CALL shft(fa,fb,fc,fu)
      END DO
   END SUBROUTINE mnbrakdp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE mnbraksp(ax,bx,cx,fa,fb,fc,func,p,xi,errorstatus) !to bracket with linmin
      !modified version of NR
      !p is the initial vector position of the state (pls. check)
      !xi is the search direction (pls. check)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(sp), INTENT(INOUT) :: ax,bx
      REAL(sp), DIMENSION(:), INTENT(IN) :: p,xi
      REAL(sp), INTENT(OUT) :: cx,fa,fb,fc
      REAL(sp) :: fd
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), DIMENSION(:), INTENT(IN) :: x
            REAL(sp) :: func
         END FUNCTION func
      END INTERFACE
      REAL(sp), PARAMETER :: GOLD=1.618034_sp,GLIMIT=100.0_sp,TINY=1.0e-20_sp !original one
      !REAL(sp), PARAMETER :: GOLD=1.2_sp,GLIMIT=100.0_sp,TINY=1.0e-20_sp !abhijit
      REAL(sp) :: fu,q,r,u,ulim

      fa=func(p+ax*xi,errorstatus)
      IF (errorstatus/=0) RETURN
      fb=func(p+bx*xi,errorstatus)
      IF (errorstatus/=0) RETURN

      IF (fb > fa) THEN
         CALL swap(ax,bx)
         CALL swap(fa,fb)
      END IF
      cx=bx+GOLD*(bx-ax)
      fc=func(p+cx*xi,errorstatus)
      IF (errorstatus/=0) RETURN

      DO
         IF (fb < fc) RETURN
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_sp*SIGN(MAX(ABS(q-r),TINY),q-r))
         ulim=bx+GLIMIT*(cx-bx)
         IF ((bx-u)*(u-cx) > 0.0) THEN
            fu=func(p+u*xi,errorstatus)
            IF (errorstatus/=0) RETURN
            IF (fu < fc) THEN
               ax=bx
               fa=fb
               bx=u
               fb=fu
               RETURN
            ELSE IF (fu > fb) THEN
               cx=u
               fc=fu
               RETURN
            END IF
            u=cx+GOLD*(cx-bx)
            fu=func(p+u*xi,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE IF ((cx-u)*(u-ulim) > 0.0) THEN
            fu=func(p+u*xi,errorstatus)
            IF (errorstatus/=0) RETURN
            IF (fu < fc) THEN
               bx=cx
               cx=u
               u=cx+GOLD*(cx-bx)
               fd=func(p+u*xi,errorstatus)
               IF (errorstatus/=0) RETURN
               CALL shftsp(fb,fc,fu,fd)
            END IF
         ELSE IF ((u-ulim)*(ulim-cx) >= 0.0) THEN
            u=ulim
            fu=func(p+u*xi,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE
            u=cx+GOLD*(cx-bx)
            fu=func(p+u*xi,errorstatus)
            IF (errorstatus/=0) RETURN
         END IF
         CALL shftsp(ax,bx,cx,u)
         CALL shftsp(fa,fb,fc,fu)
      END DO
   END SUBROUTINE mnbraksp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION brentdp(ax,bx,cx,func,pxi,xi,tol,xmin,errorstatus)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(dp), INTENT(IN) :: ax,bx,cx,tol
      REAL(dp), DIMENSION(:), INTENT(IN) :: pxi,xi
      REAL(dp), INTENT(OUT) :: xmin
      REAL(dp) :: brentdp
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), DIMENSION(:), INTENT(IN) :: x
            REAL(dp) :: func
         END FUNCTION func
      END INTERFACE
      INTEGER(I4B), PARAMETER :: ITMAX=100
      REAL(dp), PARAMETER :: CGOLD=0.3819660_dp,ZEPS=1.0e-3_dp*epsilon(ax)
      INTEGER(I4B) :: iter
      REAL(dp) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      
      a=min(ax,cx)
      b=MAX(ax,cx)
      v=bx
      w=v
      x=v
      e=0.0
      fx=func(pxi+x*xi,errorstatus)
      IF (errorstatus/=0) RETURN
      fv=fx
      fw=fx
      DO iter=1,ITMAX
         xm=0.5_dp*(a+b)
         tol1=tol*ABS(x)+ZEPS
         tol2=2.0_dp*tol1
         IF (ABS(x-xm) <= (tol2-0.5_dp*(b-a))) THEN
            xmin=x
            brentdp=fx
            RETURN
         END IF
         IF (ABS(e) > tol1) THEN
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.0_dp*(q-r)
            IF (q > 0.0) p=-p
            q=ABS(q)
            etemp=e
            e=d
            IF (ABS(p) >= ABS(0.5_dp*q*etemp) .or. &
               p <= q*(a-x) .or. p >= q*(b-x)) THEN
               e=MERGE(a-x,b-x, x >= xm )
               d=CGOLD*e
            ELSE
               d=p/q
               u=x+d
               IF (u-a < tol2 .or. b-u < tol2) d=SIGN(tol1,xm-x)
            END IF
         ELSE
            e=MERGE(a-x,b-x, x >= xm )
            d=CGOLD*e
         END IF
         u=MERGE(x+d,x+SIGN(tol1,d), ABS(d) >= tol1 )
         fu=func(pxi+u*xi,errorstatus)
         IF (errorstatus/=0) RETURN
         IF (fu <= fx) THEN
            IF (u >= x) THEN
               a=x
            ELSE
               b=x
            END IF
            CALL shft(v,w,x,u)
            CALL shft(fv,fw,fx,fu)
         ELSE
            IF (u < x) THEN
               a=u
            ELSE
               b=u
            END IF
            IF (fu <= fw .or. w == x) THEN
               v=w
               fv=fw
               w=u
               fw=fu
            ELSE IF (fu <= fv .or. v == x .or. v == w) THEN
               v=u
               fv=fu
            END IF
         END IF
      END DO
      WRITE(*,*) "Err>> brentdp: exceed maximum iterations"
   END FUNCTION brentdp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION brentsp(ax,bx,cx,func,pxi,xi,tol,xmin,errorstatus)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(sp), INTENT(IN) :: ax,bx,cx,tol
      REAL(sp), DIMENSION(:), INTENT(IN) :: pxi,xi
      REAL(sp), INTENT(OUT) :: xmin
      REAL(sp) :: brentsp
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), DIMENSION(:), INTENT(IN) :: x
            REAL(sp) :: func
         END FUNCTION func
      END INTERFACE
      INTEGER(I4B), PARAMETER :: ITMAX=100
      REAL(sp), PARAMETER :: CGOLD=0.3819660_sp,ZEPS=1.0e-3_sp*epsilon(ax)
      INTEGER(I4B) :: iter
      REAL(sp) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      
      a=min(ax,cx)
      b=MAX(ax,cx)
      v=bx
      w=v
      x=v
      e=0.0
      fx=func(pxi+x*xi,errorstatus)
      IF (errorstatus/=0) RETURN
      fv=fx
      fw=fx
      DO iter=1,ITMAX
         xm=0.5_sp*(a+b)
         tol1=tol*ABS(x)+ZEPS
         tol2=2.0_sp*tol1
         IF (ABS(x-xm) <= (tol2-0.5_sp*(b-a))) THEN
            xmin=x
            brentsp=fx
            RETURN
         END IF
         IF (ABS(e) > tol1) THEN
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.0_sp*(q-r)
            IF (q > 0.0) p=-p
            q=ABS(q)
            etemp=e
            e=d
            IF (ABS(p) >= ABS(0.5_sp*q*etemp) .or. &
               p <= q*(a-x) .or. p >= q*(b-x)) THEN
               e=MERGE(a-x,b-x, x >= xm )
               d=CGOLD*e
            ELSE
               d=p/q
               u=x+d
               IF (u-a < tol2 .or. b-u < tol2) d=SIGN(tol1,xm-x)
            END IF
         ELSE
            e=MERGE(a-x,b-x, x >= xm )
            d=CGOLD*e
         END IF
         u=MERGE(x+d,x+SIGN(tol1,d), ABS(d) >= tol1 )
         fu=func(pxi+u*xi,errorstatus)
         IF (errorstatus/=0) RETURN
         IF (fu <= fx) THEN
            IF (u >= x) THEN
               a=x
            ELSE
               b=x
            END IF
            CALL shftsp(v,w,x,u)
            CALL shftsp(fv,fw,fx,fu)
         ELSE
            IF (u < x) THEN
               a=u
            ELSE
               b=u
            END IF
            IF (fu <= fw .or. w == x) THEN
               v=w
               fv=fw
               w=u
               fw=fu
            ELSE IF (fu <= fv .or. v == x .or. v == w) THEN
               v=u
               fv=fu
            END IF
         END IF
      END DO
      WRITE(*,*) "Err>> brentsp: exceed maximum iterations"
   END FUNCTION brentsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION dbrentdp(ax,bx,cx,func,dfunc,tol,xmin,errorstatus)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(dp), INTENT(IN) :: ax,bx,cx,tol
      REAL(dp), INTENT(OUT) :: xmin
      REAL(dp) :: dbrentdp
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), INTENT(IN) :: x
            REAL(dp) :: func
         END FUNCTION func
         
         FUNCTION dfunc(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), INTENT(IN) :: x
            REAL(dp) :: dfunc
         END FUNCTION dfunc
      END INTERFACE
      INTEGER(I4B), PARAMETER :: ITMAX=100
      REAL(dp), PARAMETER :: ZEPS=1.0e-3_dp*epsilon(ax)
      INTEGER(I4B) :: iter
      REAL(dp) :: a,b,d,d1,d2,du,dv,dw,dx,e
      REAL(dp) :: fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm
      LOGICAL :: ok1,ok2
      
      a=MIN(ax,cx)
      b=MAX(ax,cx)
      v=bx
      w=v
      x=v
      e=0.0
      fx=func(x,errorstatus)
      IF (errorstatus/=0) RETURN
      fv=fx
      fw=fx
      dx=dfunc(x,errorstatus)
      IF (errorstatus/=0) RETURN
      dv=dx
      dw=dx
      DO iter=1,ITMAX
         xm=0.5_dp*(a+b)
         tol1=tol*ABS(x)+ZEPS
         tol2=2.0_dp*tol1
         IF (ABS(x-xm) <= (tol2-0.5_dp*(b-a))) EXIT
         IF (ABS(e) > tol1) THEN
            d1=2.0_dp*(b-a)
            d2=d1
            IF (dw /= dx) d1=(w-x)*dx/(dx-dw)
            IF (dv /= dx) d2=(v-x)*dx/(dx-dv)
            u1=x+d1
            u2=x+d2
            ok1=((a-u1)*(u1-b) > 0.0) .and. (dx*d1 <= 0.0)
            ok2=((a-u2)*(u2-b) > 0.0) .and. (dx*d2 <= 0.0)
            olde=e
            e=d
            IF (ok1 .or. ok2) THEN
               IF (ok1 .and. ok2) THEN
                  d=MERGE(d1,d2, ABS(d1) < ABS(d2))
               ELSE
                  d=MERGE(d1,d2,ok1)
               END IF
               IF (ABS(d) <= ABS(0.5_dp*olde)) THEN
                  u=x+d
                  IF (u-a < tol2 .or. b-u < tol2) &
                     d=SIGN(tol1,xm-x)
               ELSE
                  e=MERGE(a,b, dx >= 0.0)-x
                  d=0.5_dp*e
               END IF
            ELSE
               e=MERGE(a,b, dx >= 0.0)-x
               d=0.5_dp*e
            END IF
         ELSE
            e=MERGE(a,b, dx >= 0.0)-x
            d=0.5_dp*e
         END IF
         IF (ABS(d) >= tol1) THEN
            u=x+d
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE
            u=x+SIGN(tol1,d)
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
            IF (fu > fx) exit
         END IF
         du=dfunc(u,errorstatus)
         IF (errorstatus/=0) RETURN
         IF (fu <= fx) THEN
            IF (u >= x) THEN
               a=x
            ELSE
               b=x
            END IF
            CALL mov3(v,fv,dv,w,fw,dw)
            CALL mov3(w,fw,dw,x,fx,dx)
            CALL mov3(x,fx,dx,u,fu,du)
         ELSE
            IF (u < x) THEN
               a=u
            ELSE
               b=u
            END IF
            IF (fu <= fw .or. w == x) THEN
               CALL mov3(v,fv,dv,w,fw,dw)
               CALL mov3(w,fw,dw,u,fu,du)
            ELSE IF (fu <= fv .or. v == x .or. v == w) THEN
               CALL mov3(v,fv,dv,u,fu,du)
            END IF
         END IF
      END DO
      IF (iter > ITMAX) THEN
         WRITE(*,*) "Err>> dbrentdp: exceeded maximum iterations"
         STOP
      END IF
      xmin=x
      dbrentdp=fx
   END FUNCTION dbrentdp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION dbrentsp(ax,bx,cx,func,dfunc,tol,xmin,errorstatus)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(sp), INTENT(IN) :: ax,bx,cx,tol
      REAL(sp), INTENT(OUT) :: xmin
      REAL(sp) :: dbrentsp
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), INTENT(IN) :: x
            REAL(sp) :: func
         END FUNCTION func
         
         FUNCTION dfunc(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), INTENT(IN) :: x
            REAL(sp) :: dfunc
         END FUNCTION dfunc
      END INTERFACE
      INTEGER(I4B), PARAMETER :: ITMAX=100
      REAL(sp), PARAMETER :: ZEPS=1.0e-3_sp*epsilon(ax)
      INTEGER(I4B) :: iter
      REAL(sp) :: a,b,d,d1,d2,du,dv,dw,dx,e
      REAL(sp) :: fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm
      LOGICAL :: ok1,ok2
      
      a=MIN(ax,cx)
      b=MAX(ax,cx)
      v=bx
      w=v
      x=v
      e=0.0
      fx=func(x,errorstatus)
      IF (errorstatus/=0) RETURN
      fv=fx
      fw=fx
      dx=dfunc(x,errorstatus)
      IF (errorstatus/=0) RETURN
      dv=dx
      dw=dx
      DO iter=1,ITMAX
         xm=0.5_sp*(a+b)
         tol1=tol*ABS(x)+ZEPS
         tol2=2.0_sp*tol1
         IF (ABS(x-xm) <= (tol2-0.5_sp*(b-a))) EXIT
         IF (ABS(e) > tol1) THEN
            d1=2.0_sp*(b-a)
            d2=d1
            IF (dw /= dx) d1=(w-x)*dx/(dx-dw)
            IF (dv /= dx) d2=(v-x)*dx/(dx-dv)
            u1=x+d1
            u2=x+d2
            ok1=((a-u1)*(u1-b) > 0.0) .and. (dx*d1 <= 0.0)
            ok2=((a-u2)*(u2-b) > 0.0) .and. (dx*d2 <= 0.0)
            olde=e
            e=d
            IF (ok1 .or. ok2) THEN
               IF (ok1 .and. ok2) THEN
                  d=MERGE(d1,d2, ABS(d1) < ABS(d2))
               ELSE
                  d=MERGE(d1,d2,ok1)
               END IF
               IF (ABS(d) <= ABS(0.5_sp*olde)) THEN
                  u=x+d
                  IF (u-a < tol2 .or. b-u < tol2) &
                     d=SIGN(tol1,xm-x)
               ELSE
                  e=MERGE(a,b, dx >= 0.0)-x
                  d=0.5_sp*e
               END IF
            ELSE
               e=MERGE(a,b, dx >= 0.0)-x
               d=0.5_sp*e
            END IF
         ELSE
            e=MERGE(a,b, dx >= 0.0)-x
            d=0.5_sp*e
         END IF
         IF (ABS(d) >= tol1) THEN
            u=x+d
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
         ELSE
            u=x+SIGN(tol1,d)
            fu=func(u,errorstatus)
            IF (errorstatus/=0) RETURN
            IF (fu > fx) exit
         END IF
         du=dfunc(u,errorstatus)
         IF (errorstatus/=0) RETURN
         IF (fu <= fx) THEN
            IF (u >= x) THEN
               a=x
            ELSE
               b=x
            END IF
            CALL mov3sp(v,fv,dv,w,fw,dw)
            CALL mov3sp(w,fw,dw,x,fx,dx)
            CALL mov3sp(x,fx,dx,u,fu,du)
         ELSE
            IF (u < x) THEN
               a=u
            ELSE
               b=u
            END IF
            IF (fu <= fw .or. w == x) THEN
               CALL mov3sp(v,fv,dv,w,fw,dw)
               CALL mov3sp(w,fw,dw,u,fu,du)
            ELSE IF (fu <= fv .or. v == x .or. v == w) THEN
               CALL mov3sp(v,fv,dv,u,fu,du)
            END IF
         END IF
      END DO
      IF (iter > ITMAX) THEN
         WRITE(*,*) "Err>> dbrentsp: exceeded maximum iterations"
         STOP
      END IF
      xmin=x
      dbrentsp=fx
   END FUNCTION dbrentsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE lnsrchdp(xold,fold,g,p,x,f,stpmax,check,func,errorstatus)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(dp), DIMENSION(:), INTENT(IN) :: xold,g
      REAL(dp), DIMENSION(:), INTENT(INOUT) :: p
      REAL(dp), INTENT(IN) :: fold,stpmax
      REAL(dp), DIMENSION(:), INTENT(OUT) :: x
      REAL(dp), INTENT(OUT) :: f
      LOGICAL(LGT), INTENT(OUT) :: check
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp) :: func
            REAL(dp), DIMENSION(:), INTENT(IN) :: x
         END FUNCTION func
      END INTERFACE
      REAL(dp), PARAMETER :: ALF=1.0e-4_dp,TOLX=epsilon(x)
      INTEGER(I4B) :: ndum
      REAL(dp) :: a,alam,alam2,alamin,b,disc,f2,fold2,pABS,rhs1,rhs2,slope,&
         tmplam
         
      ndum=assert_eq4(size(g),size(p),size(x),size(xold),'lnsrch')
      check=.false.
      pABS=vABS(p(:))
      IF (pABS > stpmax) p(:)=p(:)*stpmax/pABS
      slope=DOT_PRODUCT(g,p)
      alamin=TOLX/MAXVAL(ABS(p(:))/MAX(ABS(xold(:)),1.0_dp))
      alam=1.0
      DO
         x(:)=xold(:)+alam*p(:)
         f=func(x,errorstatus)
         IF (errorstatus/=0) RETURN
         IF (alam < alamin) THEN
            x(:)=xold(:)
            check=.true.
            RETURN
         ELSE IF (f <= fold+ALF*alam*slope) THEN
            RETURN
         ELSE
            IF (alam == 1.0) THEN
               tmplam=-slope/(2.0_dp*(f-fold-slope))
            ELSE
               rhs1=f-fold-alam*slope
               rhs2=f2-fold2-alam2*slope
               a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
               b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                  (alam-alam2)
               IF (a == 0.0) THEN
                  tmplam=-slope/(2.0_dp*b)
               ELSE
                  disc=b*b-3.0_dp*a*slope
                  IF (disc < 0.0) THEN
                     WRITE(*,*) "Err>>lnsrch: roundoff problem"
                     STOP
                  END IF
                  tmplam=(-b+sqrt(disc))/(3.0_dp*a)
               END IF
               IF (tmplam > 0.5_dp*alam) tmplam=0.5_dp*alam
            END IF
         END IF
         alam2=alam
         f2=f
         fold2=fold
         alam=MAX(tmplam,0.1_dp*alam)
      END DO
   END SUBROUTINE lnsrchdp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE lnsrchsp(xold,fold,g,p,x,f,stpmax,check,func,errorstatus)
      IMPLICIT NONE
      INTEGER :: errorstatus
      REAL(sp), DIMENSION(:), INTENT(IN) :: xold,g
      REAL(sp), DIMENSION(:), INTENT(INOUT) :: p
      REAL(sp), INTENT(IN) :: fold,stpmax
      REAL(sp), DIMENSION(:), INTENT(OUT) :: x
      REAL(sp), INTENT(OUT) :: f
      LOGICAL(LGT), INTENT(OUT) :: check
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp) :: func
            REAL(sp), DIMENSION(:), INTENT(IN) :: x
         END FUNCTION func
      END INTERFACE
      REAL(sp), PARAMETER :: ALF=1.0e-4_sp,TOLX=epsilon(x)
      INTEGER(I4B) :: ndum
      REAL(sp) :: a,alam,alam2,alamin,b,disc,f2,fold2,pABS,rhs1,rhs2,slope,&
         tmplam
         
      ndum=assert_eq4(size(g),size(p),size(x),size(xold),'lnsrch')
      check=.false.
      pABS=vABSsp(p(:))
      IF (pABS > stpmax) p(:)=p(:)*stpmax/pABS
      slope=DOT_PRODUCT(g,p)
      alamin=TOLX/MAXVAL(ABS(p(:))/MAX(ABS(xold(:)),1.0_sp))
      alam=1.0
      DO
         x(:)=xold(:)+alam*p(:)
         f=func(x,errorstatus)
         IF (errorstatus/=0) RETURN
         IF (alam < alamin) THEN
            x(:)=xold(:)
            check=.true.
            RETURN
         ELSE IF (f <= fold+ALF*alam*slope) THEN
            RETURN
         ELSE
            IF (alam == 1.0) THEN
               tmplam=-slope/(2.0_sp*(f-fold-slope))
            ELSE
               rhs1=f-fold-alam*slope
               rhs2=f2-fold2-alam2*slope
               a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
               b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
                  (alam-alam2)
               IF (a == 0.0) THEN
                  tmplam=-slope/(2.0_sp*b)
               ELSE
                  disc=b*b-3.0_sp*a*slope
                  IF (disc < 0.0) THEN
                     WRITE(*,*) "Err>>lnsrch: roundoff problem"
                     STOP
                  END IF
                  tmplam=(-b+sqrt(disc))/(3.0_sp*a)
               END IF
               IF (tmplam > 0.5_sp*alam) tmplam=0.5_sp*alam
            END IF
         END IF
         alam2=alam
         f2=f
         fold2=fold
         alam=MAX(tmplam,0.1_sp*alam)
      END DO
   END SUBROUTINE lnsrchsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Lbfgsdp(p,ftol,gtol,xtol,iter,fret,func,dfunc,IsPrint1,ITMAX,errorstatus,memcorrections)
      !this is an interface to the F77 fortran code
      IMPLICIT NONE
      INTEGER, PARAMETER :: nmax=300000,mmax=50
      INTEGER :: ITMAX  !=100 !5000
      REAL(dp), DIMENSION(:), INTENT(INOUT) :: p
      REAL(dp), INTENT(IN) :: ftol,gtol,xtol
      REAL(dp), INTENT(OUT) :: fret
      INTEGER, INTENT(OUT) :: iter
      CHARACTER(len=60) :: task,csave
      INTEGER :: n,m,iprint,iwa(3*nmax),isave(44)
      REAL(dp), DIMENSION(SIZE(p)) :: g,pold
      REAL(dp) :: dsave(29),f,wa(2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax)
      REAL(dp) :: pgtol,factr,fold
      REAL(dp), DIMENSION(nmax) :: nbd,l,u
      REAL(dp) :: dpmax,dV,gg
      INTEGER :: i,errorstatus
      INTEGER, OPTIONAL :: memcorrections
      LOGICAL :: Converged,lsave(4)
      LOGICAL, OPTIONAL :: IsPrint1
      LOGICAL :: IsPrint
      
      EXTERNAL setulb
      
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), DIMENSION(:), INTENT(IN) :: x
            REAL(dp) :: func
         END FUNCTION func
         FUNCTION dfunc(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(dp), DIMENSION(:), INTENT(IN) :: x
            REAL(dp), DIMENSION(SIZE(x)) :: dfunc
         END FUNCTION dfunc
!          SUBROUTINE setulb(n,m,p,l,u,nbd,f,g,factr,pgtol,wa, &
!             iwa,task,iprint,csave,lsave,isave,dsave)
!             INTEGER, INTENT(IN) :: n,m,iprint
!             INTEGER, DIMENSION(:), INTENT(INOUT) :: iwa,isave
!             DOUBLE PRECISION, INTENT(IN) :: f,factr,pgtol
!             DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: p,wa,dsave
!             DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: l,u,g,nbd
!             CHARACTER(len=60), INTENT(INOUT) :: task,csave
!             LOGICAL, DIMENSION(:), INTENT(INOUT) :: lsave
!          END SUBROUTINE setulb
      END INTERFACE
      
      IF (PRESENT(IsPrint1)) THEN
         IsPrint=IsPrint1
      ELSE
         IsPrint=.FALSE.
      END IF

      iprint=-1 !freq of printing
      factr=1.e7_dp
      pgtol=1.e-8_dp
      
      n=SIZE(p) !dimension of the system
      IF (PRESENT(memcorrections)) THEN
         m=memcorrections
      ELSE
         m=5 !# limited memory corrections stored (earlier with m=25 AC got problems)
      END IF
      nbd=0
      l=1._dp
      u=1._dp
      
      task="START"
      iter=0
      f=func(p,errorstatus)
      IF (errorstatus/=0) RETURN
      IF (IsPrint) WRITE(UNIT=*,FMT='("...",I3," fp:",ES14.6)') iter,f

      pold=p
      fold=f
      CALL setulb(n,m,p,l,u,nbd,f,g,factr,pgtol,wa, &
         iwa,task,iprint,csave,lsave,isave,dsave)
      !IF (iter>2) Converged=(MAXVAL(ABS(p-pold))<1.e-10_dp)
      !IF (iter>2) Converged=Converged .OR. (ABS(f-fold)<1.e-10_dp)
      IF (task(1:2)=="FG") THEN
         f=func(p,errorstatus)
         IF (errorstatus/=0) RETURN
         fret=f
         g=dfunc(p,errorstatus)
         IF (errorstatus/=0) RETURN
         gg=DOT_PRODUCT(g,g)
         dpmax=MAXVAL(ABS(p-pold))
         dV=f-fold
      ELSEIF (task(1:5)=="NEW_X") THEN
         iter=iter+1
         IF (IsPrint .AND. MOD(iter,10)==1) WRITE(UNIT=6,FMT='("...",I3," fp: ",ES14.6, &
            " dpmax:",ES14.6, &
            " gg:",ES14.6," dV:",ES14.6)') iter,f,dpmax,gg,dV
         IF (dpmax<xtol .AND. gg<gtol .AND. ABS(dV)<ftol) RETURN
      END IF
         
      DO i=2,ITMAX
         pold=p
         fold=f
         CALL setulb(n,m,p,l,u,nbd,f,g,factr,pgtol,wa, &
            iwa,task,iprint,csave,lsave,isave,dsave)
         !IF (iter>2) Converged=(MAXVAL(ABS(p-pold))<1.e-10_dp)
         !IF (iter>2) Converged=Converged .OR. (ABS(f-fold)<1.e-10_dp)
         IF (task(1:2)=="FG") THEN
            f=func(p,errorstatus)
            IF (errorstatus/=0) RETURN
            fret=f
            g=dfunc(p,errorstatus)
            IF (errorstatus/=0) RETURN
            gg=DOT_PRODUCT(g,g)
            dpmax=MAXVAL(ABS(p-pold))
            dV=f-fold
         ELSEIF (task(1:5)=="NEW_X") THEN
            iter=iter+1
            !IF (IsPrint .AND. MOD(iter,50)==1) &
            !   WRITE(UNIT=6,FMT='("...",I3," fp: ",ES14.6, &
            !   " dpmax:",ES14.6, &
            !   " gg:",ES14.6," dV:",ES14.6)') iter,f,dpmax,gg,dV
            IF (dpmax<xtol .AND. gg<gtol .AND. ABS(dV)<ftol) THEN
               IF (IsPrint) WRITE(UNIT=6,FMT='("...",I3," fp: ",ES14.6, &
                 " dpmax:",ES14.6, &
                 " gg:",ES14.6," dV:",ES14.6)') iter,f,dpmax,gg,dV
               RETURN
            END IF
         END IF
      END DO
      IF (IsPrint) WRITE(UNIT=6,FMT='("...",I3," fp: ",ES14.6, &
         " dpmax:",ES14.6, &
         " gg:",ES14.6," dV:",ES14.6)') iter,f,dpmax,gg,dV
      IF (IsPrint) WRITE(6,*) "$Err>>L-BFGS: maximum iterations exceeded"
   END SUBROUTINE Lbfgsdp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Lbfgssp(p,ftol,gtol,xtol,iter,fret,func,dfunc,IsPrint1,ITMAX,errorstatus)
      !this is an interface to the F77 fortran code
      IMPLICIT NONE
      INTEGER, PARAMETER :: nmax=300000,mmax=50
      INTEGER :: ITMAX  !=100 !5000
      REAL(sp), DIMENSION(:), INTENT(INOUT) :: p
      REAL(sp), INTENT(IN) :: ftol,gtol,xtol
      REAL(sp), INTENT(OUT) :: fret
      INTEGER, INTENT(OUT) :: iter
      CHARACTER(len=60) :: task,csave
      INTEGER :: n,m,iprint,iwa(3*nmax),isave(44)
      REAL(sp), DIMENSION(SIZE(p)) :: g,pold
      REAL(sp) :: dsave(29),f,wa(2*mmax*nmax+4*nmax+12*mmax*mmax+12*mmax)
      REAL(sp) :: pgtol,factr,fold
      REAL(sp), DIMENSION(nmax) :: nbd,l,u
      REAL(sp) :: dpmax,dV,gg
      INTEGER :: i,errorstatus
      LOGICAL :: Converged,lsave(4)
      LOGICAL, OPTIONAL :: IsPrint1
      LOGICAL :: IsPrint
      
      INTERFACE
         FUNCTION func(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), DIMENSION(:), INTENT(IN) :: x
            REAL(sp) :: func
         END FUNCTION func
         FUNCTION dfunc(x,errorstatus)
            USE VARIABLE_TYPE
            IMPLICIT NONE
            INTEGER :: errorstatus
            REAL(sp), DIMENSION(:), INTENT(IN) :: x
            REAL(sp), DIMENSION(SIZE(x)) :: dfunc
         END FUNCTION dfunc
!          SUBROUTINE setulb(n,m,p,l,u,nbd,f,g,factr,pgtol,wa, &
!             iwa,task,iprint,csave,lsave,isave,dsave)
!             INTEGER, INTENT(IN) :: n,m,iprint
!             INTEGER, DIMENSION(:), INTENT(INOUT) :: iwa,isave
!             DOUBLE PRECISION, INTENT(IN) :: f,factr,pgtol
!             DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: p,wa,dsave
!             DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: l,u,g,nbd
!             CHARACTER(len=60), INTENT(INOUT) :: task,csave
!             LOGICAL, DIMENSION(:), INTENT(INOUT) :: lsave
!          END SUBROUTINE setulb
      END INTERFACE
      
      IF (PRESENT(IsPrint1)) THEN
         IsPrint=IsPrint1
      ELSE
         IsPrint=.FALSE.
      END IF

      iprint=-1 !freq of printing
      factr=1.e7_sp
      pgtol=1.e-8_sp
      
      n=SIZE(p) !dimension of the system
      m=25 !# limited memory corrections stored
      nbd=0
      l=1._sp
      u=1._sp
      
      task="START"
      iter=0
      f=func(p,errorstatus)
      IF (errorstatus/=0) RETURN
      IF (IsPrint) WRITE(UNIT=*,FMT='("...",I3," fp:",ES14.6)') iter,f

      pold=p
      fold=f
      CALL setulb(n,m,p,l,u,nbd,f,g,factr,pgtol,wa, &
         iwa,task,iprint,csave,lsave,isave,dsave)
      !IF (iter>2) Converged=(MAXVAL(ABS(p-pold))<1.e-10_sp)
      !IF (iter>2) Converged=Converged .OR. (ABS(f-fold)<1.e-10_sp)
      IF (task(1:2)=="FG") THEN
         f=func(p,errorstatus)
         IF (errorstatus/=0) RETURN
         fret=f
         g=dfunc(p,errorstatus)
         IF (errorstatus/=0) RETURN
         gg=DOT_PRODUCT(g,g)
         dpmax=MAXVAL(ABS(p-pold))
         dV=f-fold
      ELSEIF (task(1:5)=="NEW_X") THEN
         iter=iter+1
         IF (IsPrint .AND. MOD(iter,10)==1) WRITE(UNIT=6,FMT='("...",I3," fp: ",ES14.6, &
            " dpmax:",ES14.6, &
            " gg:",ES14.6," dV:",ES14.6)') iter,f,dpmax,gg,dV
         IF (dpmax<xtol .AND. gg<gtol .AND. ABS(dV)<ftol) RETURN
      END IF
         
      DO i=2,ITMAX
         pold=p
         fold=f
         CALL setulb(n,m,p,l,u,nbd,f,g,factr,pgtol,wa, &
            iwa,task,iprint,csave,lsave,isave,dsave)
         !IF (iter>2) Converged=(MAXVAL(ABS(p-pold))<1.e-10_sp)
         !IF (iter>2) Converged=Converged .OR. (ABS(f-fold)<1.e-10_sp)
         IF (task(1:2)=="FG") THEN
            f=func(p,errorstatus)
            IF (errorstatus/=0) RETURN
            fret=f
            g=dfunc(p,errorstatus)
            IF (errorstatus/=0) RETURN
            gg=DOT_PRODUCT(g,g)
            dpmax=MAXVAL(ABS(p-pold))
            dV=f-fold
         ELSEIF (task(1:5)=="NEW_X") THEN
            iter=iter+1
            !IF (IsPrint .AND. MOD(iter,50)==1) &
            !   WRITE(UNIT=6,FMT='("...",I3," fp: ",ES14.6, &
            !   " dpmax:",ES14.6, &
            !   " gg:",ES14.6," dV:",ES14.6)') iter,f,dpmax,gg,dV
            IF (dpmax<xtol .AND. gg<gtol .AND. ABS(dV)<ftol) THEN
               IF (IsPrint) WRITE(UNIT=6,FMT='("...",I3," fp: ",ES14.6, &
                 " dpmax:",ES14.6, &
                 " gg:",ES14.6," dV:",ES14.6)') iter,f,dpmax,gg,dV
               RETURN
            END IF
         END IF
      END DO
      IF (IsPrint) WRITE(UNIT=6,FMT='("...",I3," fp: ",ES14.6, &
         " dpmax:",ES14.6, &
         " gg:",ES14.6," dV:",ES14.6)') iter,f,dpmax,gg,dV
      WRITE(6,*) "$Err>>L-BFGS: maximum iterations exceeded"
   END SUBROUTINE Lbfgssp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE OptimPackage
