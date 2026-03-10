MODULE nrprocedures 
   USE VARIABLE_TYPE
   IMPLICIT NONE
   INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
   INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
   INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
   INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
   INTEGER(I4B), PARAMETER :: NPAR_POLY=8
   INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
   
   INTERFACE reallocate
      MODULE PROCEDURE reallocate_rv,reallocate_rm,&
         reallocate_iv,reallocate_im,reallocate_hv
   END INTERFACE
   INTERFACE arth
      MODULE PROCEDURE arth_r, arth_d, arth_i
   END INTERFACE
   INTERFACE gammln
      MODULE PROCEDURE gammln_s,gammln_v
   END INTERFACE gammln
   INTERFACE gammp
      MODULE PROCEDURE gammp_s,gammp_v
   END INTERFACE gammp
   INTERFACE gcf
      MODULE PROCEDURE gcf_s,gcf_v
   END INTERFACE gcf
   INTERFACE gser
      MODULE PROCEDURE gser_s,gser_v
   END INTERFACE gser
   INTERFACE array_copy
      MODULE PROCEDURE array_copy_r, array_copy_d, array_copy_i
   END INTERFACE
   INTERFACE swap
      MODULE PROCEDURE swap_i,swap_r,swap_d,swap_rv,swap_c, &
         swap_cv,swap_cm,swap_z,swap_zv,swap_zm !, &
         !masked_swap_rs,masked_swap_rv,masked_swap_rm
   END INTERFACE
   INTERFACE assert
      MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
   END INTERFACE
   INTERFACE assert_eq
      MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
   END INTERFACE
   INTERFACE erf
      MODULE PROCEDURE erf_s,erf_v
   END INTERFACE erf
   INTERFACE sort_pick
      MODULE PROCEDURE sort_pick_abhijit,sort_pick_abhijit_int
   END INTERFACE sort_pick
   INTERFACE sort_heap
      MODULE PROCEDURE sort_heap_abhijit,sort_heap_abhijit_int
   END INTERFACE sort_heap
   INTERFACE sort_quick
      MODULE PROCEDURE sort_quick_abhijit,sort_quick_abhijit_int
   END INTERFACE sort_quick
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION reallocate_rv(p,n)
      REAL(SP), DIMENSION(:), POINTER :: p, reallocate_rv
      INTEGER(I4B), INTENT(IN) :: n
      INTEGER(I4B) :: nold,ierr
      
      allocate(reallocate_rv(n),stat=ierr)
      IF (ierr /= 0) &
         WRITE(*,*) 'reallocate_rv: problem in attempt to allocate memory'
      IF (.not. associated(p)) RETURN
      nold=size(p)
      reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
      deallocate(p)
   END FUNCTION reallocate_rv
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION reallocate_iv(p,n)
      INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
      INTEGER(I4B), INTENT(IN) :: n
      INTEGER(I4B) :: nold,ierr
      
      allocate(reallocate_iv(n),stat=ierr)
      IF (ierr /= 0) &
         WRITE(*,*) 'reallocate_iv: problem in attempt to allocate memory'
      IF (.not. associated(p)) RETURN
      nold=size(p)
      reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
      deallocate(p)
   END FUNCTION reallocate_iv
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION reallocate_hv(p,n)
      CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
      INTEGER(I4B), INTENT(IN) :: n
      INTEGER(I4B) :: nold,ierr
      
      allocate(reallocate_hv(n),stat=ierr)
      IF (ierr /= 0) &
         WRITE(*,*) 'reallocate_hv: problem in attempt to allocate memory'
      IF (.not. associated(p)) RETURN
      nold=size(p)
      reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
      deallocate(p)
   END FUNCTION reallocate_hv
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION reallocate_rm(p,n,m)
      IMPLICIT NONE
      REAL(SP), DIMENSION(:,:), POINTER :: p, reallocate_rm
      INTEGER(I4B), INTENT(IN) :: n,m
      INTEGER(I4B) :: nold,mold,ierr
      
      allocate(reallocate_rm(n,m),stat=ierr)
      IF (ierr /= 0) &
         WRITE(*,*) 'reallocate_rm: problem in attempt to allocate memory'
      IF (.not. associated(p)) RETURN
      nold=size(p,1)
      mold=size(p,2)
      reallocate_rm(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
      deallocate(p)
   END FUNCTION reallocate_rm
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION reallocate_im(p,n,m)
      INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
      INTEGER(I4B), INTENT(IN) :: n,m
      INTEGER(I4B) :: nold,mold,ierr
      
      allocate(reallocate_im(n,m),stat=ierr)
      IF (ierr /= 0) &
         WRITE(*,*) 'reallocate_im: problem in attempt to allocate memory'
      IF (.not. associated(p)) RETURN
      nold=size(p,1)
      mold=size(p,2)
      reallocate_im(1:min(nold,n),1:min(mold,m))=&
         p(1:min(nold,n),1:min(mold,m))
      deallocate(p)
   END FUNCTION reallocate_im
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION arth_r(first,increment,n)
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: first,increment
      INTEGER(I4B), INTENT(IN) :: n
      REAL(SP), DIMENSION(n) :: arth_r
      INTEGER(I4B) :: k,k2
      REAL(SP) :: temp
      
      IF (n > 0) arth_r(1)=first
      IF (n <= NPAR_ARTH) THEN
         DO k=2,n
            arth_r(k)=arth_r(k-1)+increment
         END DO
      ELSE
         DO k=2,NPAR2_ARTH
            arth_r(k)=arth_r(k-1)+increment
         END DO
         temp=increment*NPAR2_ARTH
         k=NPAR2_ARTH
         DO
            IF (k >= n) EXIT
            k2=k+k
            arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
            temp=temp+temp
            k=k2
         END DO
      END IF
   END FUNCTION arth_r
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION arth_d(first,increment,n)
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: first,increment
      INTEGER(I4B), INTENT(IN) :: n
      REAL(DP), DIMENSION(n) :: arth_d
      INTEGER(I4B) :: k,k2
      REAL(DP) :: temp
      
      IF (n > 0) arth_d(1)=first
      IF (n <= NPAR_ARTH) THEN
         DO k=2,n
            arth_d(k)=arth_d(k-1)+increment
         END DO
      ELSE
         DO k=2,NPAR2_ARTH
            arth_d(k)=arth_d(k-1)+increment
         END DO
         temp=increment*NPAR2_ARTH
         k=NPAR2_ARTH
         DO
            IF (k >= n) EXIT
            k2=k+k
            arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
            temp=temp+temp
            k=k2
         END DO
      END IF
   END FUNCTION arth_d
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION arth_i(first,increment,n)
      IMPLICIT NONE
      INTEGER(I4B), INTENT(IN) :: first,increment,n
      INTEGER(I4B), DIMENSION(n) :: arth_i
      INTEGER(I4B) :: k,k2,temp
      
      IF (n > 0) arth_i(1)=first
      IF (n <= NPAR_ARTH) THEN
         DO k=2,n
            arth_i(k)=arth_i(k-1)+increment
         END DO
      ELSE
         DO k=2,NPAR2_ARTH
            arth_i(k)=arth_i(k-1)+increment
         END DO
         temp=increment*NPAR2_ARTH
         k=NPAR2_ARTH
         DO
            IF (k >= n) EXIT
            k2=k+k
            arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
            temp=temp+temp
            k=k2
         END DO
      END IF
   END FUNCTION arth_i
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION vABS(v)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), INTENT(IN) :: v
      REAL(dp) :: vABS
      
      vABS=sqrt(DOT_PRODUCT(v,v))
   END FUNCTION vABS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION vABSsp(v)
      IMPLICIT NONE
      REAL(sp), DIMENSION(:), INTENT(IN) :: v
      REAL(sp) :: vABSsp
      
      vABSsp=sqrt(DOT_PRODUCT(v,v))
   END FUNCTION vABSsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION outerprod(a,b)
      IMPLICIT NONE
      REAL(dp), DIMENSION(:), INTENT(IN) :: a,b
      REAL(dp), DIMENSION(size(a),size(b)) :: outerprod
      
      outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
   END FUNCTION outerprod
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION outerprodsp(a,b)
      IMPLICIT NONE
      REAL(sp), DIMENSION(:), INTENT(IN) :: a,b
      REAL(sp), DIMENSION(size(a),size(b)) :: outerprodsp
      
      outerprodsp = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
   END FUNCTION outerprodsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE unit_matrix(mat)
      REAL(dp), DIMENSION(:,:), INTENT(OUT) :: mat
      INTEGER(I4B) :: i,n
      
      n=min(size(mat,1),size(mat,2))
      mat(:,:)=0.0_dp
      DO i=1,n
         mat(i,i)=1.0_dp
      END DO
   END SUBROUTINE unit_matrix
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE unit_matrixsp(mat)
      REAL(sp), DIMENSION(:,:), INTENT(OUT) :: mat
      INTEGER(I4B) :: i,n
      
      n=min(size(mat,1),size(mat,2))
      mat(:,:)=0.0_dp
      DO i=1,n
         mat(i,i)=1.0_dp
      END DO
   END SUBROUTINE unit_matrixsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE mov3(a,b,c,d,e,f)
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: d,e,f
      REAL(dp), INTENT(OUT) :: a,b,c
      
      a=d
      b=e
      c=f
   END SUBROUTINE mov3
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE mov3sp(a,b,c,d,e,f)
      IMPLICIT NONE
      REAL(sp), INTENT(IN) :: d,e,f
      REAL(sp), INTENT(OUT) :: a,b,c
      
      a=d
      b=e
      c=f
   END SUBROUTINE mov3sp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE shft(a,b,c,d)
      IMPLICIT NONE
      REAL(dp), INTENT(OUT) :: a
      REAL(dp), INTENT(INOUT) :: b,c
      REAL(dp), INTENT(IN) :: d
      
      a=b
      b=c
      c=d
   END SUBROUTINE shft
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE shftsp(a,b,c,d)
      IMPLICIT NONE
      REAL(sp), INTENT(OUT) :: a
      REAL(sp), INTENT(INOUT) :: b,c
      REAL(sp), INTENT(IN) :: d
      
      a=b
      b=c
      c=d
   END SUBROUTINE shftsp
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE swap_i(a,b)
      INTEGER(I4B), INTENT(INOUT) :: a,b
      INTEGER(I4B) :: dum
      
      dum=a
      a=b
      b=dum
   END SUBROUTINE swap_i
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE swap_r(a,b)
      REAL(SP), INTENT(INOUT) :: a,b
      REAL(SP) :: dum
      
      dum=a
      a=b
      b=dum
   END SUBROUTINE swap_r
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE swap_rv(a,b)
      REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
      REAL(SP), DIMENSION(SIZE(a)) :: dum
      
      dum=a
      a=b
      b=dum
   END SUBROUTINE swap_rv
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE swap_d(a,b)
      REAL(DP), INTENT(INOUT) :: a,b
      REAL(DP) :: dum
      
      dum=a
      a=b
      b=dum
   END SUBROUTINE swap_d
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE swap_c(a,b)
      COMPLEX(SPC), INTENT(INOUT) :: a,b
      COMPLEX(SPC) :: dum
      
      dum=a
      a=b
      b=dum
   END SUBROUTINE swap_c
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE swap_cv(a,b)
      COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
      COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
      
      dum=a
      a=b
      b=dum
   END SUBROUTINE swap_cv
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE swap_cm(a,b)
      COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
      COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
      
      dum=a
      a=b
      b=dum
   END SUBROUTINE swap_cm
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE swap_z(a,b)
      COMPLEX(DPC), INTENT(INOUT) :: a,b
      COMPLEX(DPC) :: dum
      
      dum=a
      a=b
      b=dum
   END SUBROUTINE swap_z
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE swap_zv(a,b)
      COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
      COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
      
      dum=a
      a=b
      b=dum
   END SUBROUTINE swap_zv
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE swap_zm(a,b)
      COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
      COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
      
      dum=a
      a=b
      b=dum
   END SUBROUTINE swap_zm
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION assert_eq2(n1,n2,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: n1,n2
      INTEGER :: assert_eq2
      
      IF (n1 == n2) THEN
         assert_eq2=n1
      ELSE
         WRITE (*,*) '$Err: an assert_eq failed with this tag:', &
            string
         STOP
      END IF
   END FUNCTION assert_eq2
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION assert_eq3(n1,n2,n3,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: n1,n2,n3
      INTEGER :: assert_eq3

      if (n1 == n2 .and. n2 == n3) then
         assert_eq3=n1
      else
         write (*,*) 'nrerror: an assert_eq failed with this tag:', &
              string
         STOP 'program terminated by assert_eq3'
      end if
   END FUNCTION assert_eq3
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION assert_eq4(n1,n2,n3,n4,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, INTENT(IN) :: n1,n2,n3,n4
      INTEGER :: assert_eq4
      
      IF (n1 == n2 .and. n2 == n3 .and. n3 == n4) THEN
         assert_eq4=n1
      ELSE
         WRITE (*,*) '$Err: an assert_eq failed with this tag:', &
            string
         STOP
      END IF
   END FUNCTION assert_eq4
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION assert_eqn(nn,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      INTEGER, DIMENSION(:), INTENT(IN) :: nn
      INTEGER :: assert_eqn

      if (all(nn(2:) == nn(1))) then
         assert_eqn=nn(1)
      else
         write (*,*) 'nrerror: an assert_eq failed with this tag:', &
              string
         STOP 'program terminated by assert_eqn'
      end if
   END FUNCTION assert_eqn
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION gammln_s(xx)
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: xx
      REAL(SP) :: gammln_s
      REAL(DP) :: tmp,x
      REAL(DP) :: stp = 2.5066282746310005_dp
      REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
         -86.50532032941677_dp,24.01409824083091_dp,&
         -1.231739572450155_dp,0.1208650973866179e-2_dp,&
         -0.5395239384953e-5_dp/)
         
      CALL assert(xx > 0.0, 'gammln_s arg')
      x=xx
      tmp=x+5.5_dp
      tmp=(x+0.5_dp)*log(tmp)-tmp
      gammln_s=tmp+log(stp*(1.000000000190015_dp+&
         sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
   END FUNCTION gammln_s
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION gammln_v(xx)
      IMPLICIT NONE
      INTEGER(I4B) :: i
      REAL(SP), DIMENSION(:), INTENT(IN) :: xx
      REAL(SP), DIMENSION(size(xx)) :: gammln_v
      REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
      REAL(DP) :: stp = 2.5066282746310005_dp
      REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
         -86.50532032941677_dp,24.01409824083091_dp,&
         -1.231739572450155_dp,0.1208650973866179e-2_dp,&
         -0.5395239384953e-5_dp/)
         
      IF (size(xx) == 0) RETURN
      CALL assert(all(xx > 0.0), 'gammln_v arg')
      x=xx
      tmp=x+5.5_dp
      tmp=(x+0.5_dp)*log(tmp)-tmp
      ser=1.000000000190015_dp
      y=x
      DO i=1,size(coef)
         y=y+1.0_dp
         ser=ser+coef(i)/y
      END DO
      gammln_v=tmp+log(stp*ser/x)
   END FUNCTION gammln_v
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE array_copy_r(src,dest,n_copied,n_not_copied)
      REAL(SP), DIMENSION(:), INTENT(IN) :: src
      REAL(SP), DIMENSION(:), INTENT(OUT) :: dest
      INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
      
      n_copied=min(size(src),size(dest))
      n_not_copied=size(src)-n_copied
      dest(1:n_copied)=src(1:n_copied)
   END SUBROUTINE array_copy_r
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE array_copy_d(src,dest,n_copied,n_not_copied)
      REAL(DP), DIMENSION(:), INTENT(IN) :: src
      REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
      INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
      
      n_copied=min(size(src),size(dest))
      n_not_copied=size(src)-n_copied
      dest(1:n_copied)=src(1:n_copied)
   END SUBROUTINE array_copy_d
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
      INTEGER(I4B), DIMENSION(:), INTENT(IN) :: src
      INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: dest
      INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
      
      n_copied=min(size(src),size(dest))
      n_not_copied=size(src)-n_copied
      dest(1:n_copied)=src(1:n_copied)
   END SUBROUTINE array_copy_i
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE assert1(n1,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1
      IF (.not. n1) THEN
         write (*,*) '$Err: an assertion failed with this tag:', &
            string
         STOP
      END IF
   END SUBROUTINE assert1
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE assert2(n1,n2,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1,n2
      IF (.not. (n1 .and. n2)) THEN
         write (*,*) '$Err: an assertion failed with this tag:', &
            string
         STOP
      END IF
   END SUBROUTINE assert2
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE assert3(n1,n2,n3,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1,n2,n3
      IF (.not. (n1 .and. n2 .and. n3)) THEN
         write (*,*) '$Err: an assertion failed with this tag:', &
            string
         STOP
      END IF
   END SUBROUTINE assert3
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE assert4(n1,n2,n3,n4,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, INTENT(IN) :: n1,n2,n3,n4
      IF (.not. (n1 .and. n2 .and. n3 .and. n4)) THEN
         write (*,*) '$Err: an assertion failed with this tag:', &
            string
         STOP
      END IF
   END SUBROUTINE assert4
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE assert_v(n,string)
      CHARACTER(LEN=*), INTENT(IN) :: string
      LOGICAL, DIMENSION(:), INTENT(IN) :: n
      IF (.not. all(n)) THEN
         write (*,*) '$Err: an assertion failed with this tag:', &
            string
         STOP
         
      END IF
   END SUBROUTINE assert_v
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  
   SUBROUTINE sort_pick_abhijit(arr,brr)
      !Trivial sort - use when SIZE(arr)<20
      IMPLICIT NONE
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
      INTEGER, DIMENSION(:), INTENT(INOUT) :: brr
      INTEGER :: i,j,n,b
      REAL(DP) :: a
      
      n=SIZE(arr)
      
      DO j=2,n
         a=arr(j)
         b=brr(j)
         DO i=j-1,1,-1
            IF (arr(i) <= a) EXIT
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
         END DO
         arr(i+1)=a
         brr(i+1)=b
      END DO
   END SUBROUTINE sort_pick_abhijit
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
   SUBROUTINE sort_heap_abhijit(arr,brr)
      !Heap sort - use when SIZE(arr)<1000
      IMPLICIT NONE
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
      INTEGER, DIMENSION(:), INTENT(INOUT) :: brr
      INTEGER(I4B) :: i,n
      
      n=SIZE(arr)
      DO i=n/2,1,-1
         CALL sift_down(i,n)
      END DO
      DO i=n,2,-1
         CALL swap(arr(1),arr(i))
         CALL swap(brr(1),brr(i))
         CALL sift_down(1,i-1)
      END DO
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE sift_down(l,r)
         INTEGER(I4B), INTENT(IN) :: l,r
         INTEGER(I4B) :: j,jold
         REAL(DP) :: a,b
         
         a=arr(l)
         b=brr(l)
         jold=l
         j=l+l
         DO
            IF (j > r) EXIT
            IF (j < r) THEN
               IF (arr(j) < arr(j+1)) j=j+1
            END IF
            IF (a >= arr(j)) EXIT
            arr(jold)=arr(j)
            brr(jold)=brr(j)
            jold=j
            j=j+j
         END DO
         arr(jold)=a
         brr(jold)=b
      END SUBROUTINE sift_down
   END SUBROUTINE sort_heap_abhijit
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE sort_quick_abhijit(arr,brr)
      !Quick sort - use when SIZE(arr)>1000
      IMPLICIT NONE
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
      INTEGER, DIMENSION(:), INTENT(INOUT) :: brr
      INTEGER, PARAMETER :: NN=15, NSTACK=50
      REAL(DP) :: a
      INTEGER :: b,n,k,i,j,jstack,l,r
      INTEGER, DIMENSION(NSTACK) :: istack
      
      n=SIZE(arr)
      jstack=0
      l=1
      r=n
      DO
         IF (r-l < NN) THEN
            DO j=l+1,r
               a=arr(j)
               b=brr(j)
               DO i=j-1,l,-1
                  IF (arr(i) <= a) EXIT
                  arr(i+1)=arr(i)
                  brr(i+1)=brr(i)
               END DO
               arr(i+1)=a
               brr(i+1)=b
            END DO
            IF (jstack == 0) RETURN
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
         ELSE
            k=(l+r)/2
            CALL swap(arr(k),arr(l+1))
            CALL swap(brr(k),brr(l+1))
            IF (arr(l)>arr(r)) THEN
               CALL swap(arr(l),arr(r))
               CALL swap(brr(l),brr(r))
            END IF
            IF (arr(l+1)>arr(r)) THEN
               CALL swap(arr(l+1),arr(r))
               CALL swap(brr(l+1),brr(r))
            END IF
            IF (arr(l)>arr(l+1)) THEN
               CALL swap(arr(l),arr(l+1))
               CALL swap(brr(l),brr(l+1))
            END IF
            i=l+1
            j=r
            a=arr(l+1)
            b=brr(l+1)
            DO
               DO
                  i=i+1
                  IF (arr(i) >= a) EXIT
               END DO
               DO
                  j=j-1
                  IF (arr(j) <= a) EXIT
               END DO
               IF (j < i) EXIT
               CALL swap(arr(i),arr(j))
               CALL swap(brr(i),brr(j))
            END DO
            arr(l+1)=arr(j)
            brr(l+1)=brr(j)
            arr(j)=a
            brr(j)=b
            jstack=jstack+2
            IF (jstack > NSTACK) THEN
               WRITE(*,*) "$Err>> Stack size small in sort_quick_abhijit"
               STOP
            END IF
            IF (r-i+1 >= j-l) THEN
               istack(jstack)=r
               istack(jstack-1)=i
               r=j-1
            ELSE
               istack(jstack)=j-1
               istack(jstack-1)=l
               l=i
            END IF
         END IF
      END DO
   END SUBROUTINE sort_quick_abhijit
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   
   
   SUBROUTINE sort_pick_abhijit_int(arr,brr)
      !Trivial sort - use when SIZE(arr)<20
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(INOUT) :: arr
      INTEGER, DIMENSION(:), INTENT(INOUT) :: brr
      INTEGER :: i,j,n,b
      REAL(DP) :: a
      
      n=SIZE(arr)
      
      DO j=2,n
         a=arr(j)
         b=brr(j)
         DO i=j-1,1,-1
            IF (arr(i) <= a) EXIT
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
         END DO
         arr(i+1)=a
         brr(i+1)=b
      END DO
   END SUBROUTINE sort_pick_abhijit_int
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
   SUBROUTINE sort_heap_abhijit_int(arr,brr)
      !Heap sort - use when SIZE(arr)<1000
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(INOUT) :: arr
      INTEGER, DIMENSION(:), INTENT(INOUT) :: brr
      INTEGER(I4B) :: i,n
      
      n=SIZE(arr)
      DO i=n/2,1,-1
         CALL sift_down(i,n)
      END DO
      DO i=n,2,-1
         CALL swap(arr(1),arr(i))
         CALL swap(brr(1),brr(i))
         CALL sift_down(1,i-1)
      END DO
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE sift_down(l,r)
         INTEGER(I4B), INTENT(IN) :: l,r
         INTEGER(I4B) :: j,jold
         REAL(DP) :: a,b
         
         a=arr(l)
         b=brr(l)
         jold=l
         j=l+l
         DO
            IF (j > r) EXIT
            IF (j < r) THEN
               IF (arr(j) < arr(j+1)) j=j+1
            END IF
            IF (a >= arr(j)) EXIT
            arr(jold)=arr(j)
            brr(jold)=brr(j)
            jold=j
            j=j+j
         END DO
         arr(jold)=a
         brr(jold)=b
      END SUBROUTINE sift_down
   END SUBROUTINE sort_heap_abhijit_int
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE sort_quick_abhijit_int(arr,brr)
      !Quick sort - use when SIZE(arr)>1000
      IMPLICIT NONE
      INTEGER, DIMENSION(:), INTENT(INOUT) :: arr
      INTEGER, DIMENSION(:), INTENT(INOUT) :: brr
      INTEGER, PARAMETER :: NN=15, NSTACK=50
      REAL(DP) :: a
      INTEGER :: b,n,k,i,j,jstack,l,r
      INTEGER, DIMENSION(NSTACK) :: istack
      
      n=SIZE(arr)
      jstack=0
      l=1
      r=n
      DO
         IF (r-l < NN) THEN
            DO j=l+1,r
               a=arr(j)
               b=brr(j)
               DO i=j-1,l,-1
                  IF (arr(i) <= a) EXIT
                  arr(i+1)=arr(i)
                  brr(i+1)=brr(i)
               END DO
               arr(i+1)=a
               brr(i+1)=b
            END DO
            IF (jstack == 0) RETURN
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
         ELSE
            k=(l+r)/2
            CALL swap(arr(k),arr(l+1))
            CALL swap(brr(k),brr(l+1))
            IF (arr(l)>arr(r)) THEN
               CALL swap(arr(l),arr(r))
               CALL swap(brr(l),brr(r))
            END IF
            IF (arr(l+1)>arr(r)) THEN
               CALL swap(arr(l+1),arr(r))
               CALL swap(brr(l+1),brr(r))
            END IF
            IF (arr(l)>arr(l+1)) THEN
               CALL swap(arr(l),arr(l+1))
               CALL swap(brr(l),brr(l+1))
            END IF
            i=l+1
            j=r
            a=arr(l+1)
            b=brr(l+1)
            DO
               DO
                  i=i+1
                  IF (arr(i) >= a) EXIT
               END DO
               DO
                  j=j-1
                  IF (arr(j) <= a) EXIT
               END DO
               IF (j < i) EXIT
               CALL swap(arr(i),arr(j))
               CALL swap(brr(i),brr(j))
            END DO
            arr(l+1)=arr(j)
            brr(l+1)=brr(j)
            arr(j)=a
            brr(j)=b
            jstack=jstack+2
            IF (jstack > NSTACK) THEN
               WRITE(*,*) "$Err>> Stack size small in sort_quick_abhijit_int"
               STOP
            END IF
            IF (r-i+1 >= j-l) THEN
               istack(jstack)=r
               istack(jstack-1)=i
               r=j-1
            ELSE
               istack(jstack)=j-1
               istack(jstack-1)=l
               l=i
            END IF
         END IF
      END DO
   END SUBROUTINE sort_quick_abhijit_int
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   
   
   FUNCTION gser_s(a,x,gln)
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: a,x
      REAL(SP), OPTIONAL, INTENT(OUT) :: gln
      REAL(SP) :: gser_s
      INTEGER(I4B), PARAMETER :: ITMAX=100
      REAL(SP), PARAMETER :: EPS=epsilon(x)
      INTEGER(I4B) :: n
      REAL(SP) :: ap,del,summ
      
      if (x == 0.0) then
         gser_s=0.0
         RETURN
      end if
      ap=a
      summ=1.0_sp/a
      del=summ
      do n=1,ITMAX
         ap=ap+1.0_sp
         del=del*x/ap
         summ=summ+del
         if (abs(del) < abs(summ)*EPS) exit
      end do
      if (n > ITMAX) then
         write(*,*) "nrerror: a too large, ITMAX too small in gser_s"
         stop
      end if
      if (present(gln)) then
         gln=gammln(a)
         gser_s=summ*exp(-x+a*log(x)-gln)
      else
         gser_s=summ*exp(-x+a*log(x)-gammln(a))
      end if
   END FUNCTION gser_s
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION gser_v(a,x,gln)
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
      REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
      REAL(SP), DIMENSION(size(a)) :: gser_v
      INTEGER(I4B), PARAMETER :: ITMAX=100
      REAL(SP), PARAMETER :: EPS=epsilon(x)
      INTEGER(I4B) :: n
      REAL(SP), DIMENSION(size(a)) :: ap,del,summ
      LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
      
      n=assert_eq(size(a),size(x),'gser_v')
      zero=(x == 0.0)
      where (zero) gser_v=0.0
      ap=a
      summ=1.0_sp/a
      del=summ
      converged=zero
      do n=1,ITMAX
         where (.not. converged)
            ap=ap+1.0_sp
            del=del*x/ap
            summ=summ+del
            converged = (abs(del) < abs(summ)*EPS)
         end where
         if (all(converged)) exit
      end do
      if (n > ITMAX) then 
         write(*,*) "nrerror: a too large, ITMAX too small in gser_v"
         stop
      end if
      if (present(gln)) then
         if (size(gln) < size(a)) then
            write(*,*) "nrerror: gser: Not enough space for gln"
            stop
         end if
         gln=gammln(a)
         where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gln)
      else
         where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gammln(a))
      end if
   END FUNCTION gser_v
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION gcf_s(a,x,gln)
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: a,x
      REAL(SP), OPTIONAL, INTENT(OUT) :: gln
      REAL(SP) :: gcf_s
      INTEGER(I4B), PARAMETER :: ITMAX=100
      REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
      INTEGER(I4B) :: i
      REAL(SP) :: an,b,c,d,del,h
      
      if (x == 0.0) then
         gcf_s=1.0
         RETURN
      end if
      b=x+1.0_sp-a
      c=1.0_sp/FPMIN
      d=1.0_sp/b
      h=d
      do i=1,ITMAX
         an=-i*(i-a)
         b=b+2.0_sp
         d=an*d+b
         if (abs(d) < FPMIN) d=FPMIN
         c=b+an/c
         if (abs(c) < FPMIN) c=FPMIN
         d=1.0_sp/d
         del=d*c
         h=h*del
         if (abs(del-1.0_sp) <= EPS) exit
      end do
      if (i > ITMAX) then
         WRITE(*,*) 'a too large, ITMAX too small in gcf_s'
         STOP
      end if
      if (present(gln)) then
         gln=gammln(a)
         gcf_s=exp(-x+a*log(x)-gln)*h
      else
         gcf_s=exp(-x+a*log(x)-gammln(a))*h
      end if
   END FUNCTION gcf_s
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION gcf_v(a,x,gln)
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
      REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
      REAL(SP), DIMENSION(size(a)) :: gcf_v
      INTEGER(I4B), PARAMETER :: ITMAX=100
      REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
      INTEGER(I4B) :: i
      REAL(SP), DIMENSION(size(a)) :: an,b,c,d,del,h
      LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
      
      i=assert_eq(size(a),size(x),'gcf_v')
      zero=(x == 0.0)
      where (zero)
         gcf_v=1.0
      elsewhere
         b=x+1.0_sp-a
         c=1.0_sp/FPMIN
         d=1.0_sp/b
         h=d
      end where
      converged=zero
      do i=1,ITMAX
         where (.not. converged)
            an=-i*(i-a)
            b=b+2.0_sp
            d=an*d+b
            d=merge(FPMIN,d, abs(d)<FPMIN )
            c=b+an/c
            c=merge(FPMIN,c, abs(c)<FPMIN )
            d=1.0_sp/d
            del=d*c
            h=h*del
            converged = (abs(del-1.0_sp)<=EPS)
         end where
         if (all(converged)) exit
      end do
      if (i > ITMAX) then
         write(*,*) 'a too large, ITMAX too small in gcf_v'
         stop
      end if
      if (present(gln)) then
         if (size(gln) < size(a)) then
            write(*,*) 'gser: Not enough space for gln'
            stop
         end if
         gln=gammln(a)
         where (.not. zero) gcf_v=exp(-x+a*log(x)-gln)*h
      else
         where (.not. zero) gcf_v=exp(-x+a*log(x)-gammln(a))*h
      end if
   END FUNCTION gcf_v
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION gammp_s(a,x)
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: a,x
      REAL(SP) :: gammp_s
      
      call assert( x >= 0.0,  a > 0.0, 'gammp_s args')
      if (x<a+1.0_sp) then
         gammp_s=gser(a,x)
      else
         gammp_s=1.0_sp-gcf(a,x)
      end if
   END FUNCTION gammp_s
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION gammp_v(a,x)
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
      REAL(SP), DIMENSION(size(x)) :: gammp_v
      LOGICAL(LGT), DIMENSION(size(x)) :: mask
      INTEGER(I4B) :: ndum
      
      ndum=assert_eq(size(a),size(x),'gammp_v')
      call assert( all(x >= 0.0),  all(a > 0.0), 'gammp_v args')
      mask = (x<a+1.0_sp)
      gammp_v=merge(gser(a,merge(x,0.0_sp,mask)), &
         1.0_sp-gcf(a,merge(x,0.0_sp,.not. mask)),mask)
   END FUNCTION gammp_v
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION erf_s(x)
      IMPLICIT NONE
      REAL(sp), INTENT(IN) :: x
      REAL(sp) :: erf_s
      
      erf_s=gammp(0.5_sp,x**2)
      if (x < 0.0) erf_s=-erf_s
   END FUNCTION erf_s
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION erf_v(x)
      IMPLICIT NONE
      REAL(sp), DIMENSION(:), INTENT(IN) :: x
      REAL(sp), DIMENSION(size(x)) :: erf_v
      
      erf_v=gammp(spread(0.5_sp,1,size(x)),x**2)
      where (x < 0.0) erf_v=-erf_v
   END FUNCTION erf_v
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE nrprocedures
