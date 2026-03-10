MODULE ran_state
   USE VARIABLE_TYPE
   USE nrprocedures
   IMPLICIT NONE
   INTEGER, PARAMETER :: K4B=selected_int_kind(9)
   INTEGER(K4B), PARAMETER :: hg=huge(1_K4B), hgm=-hg, hgng=hgm-1
   INTEGER(K4B), SAVE :: lenran=0, seq=0
   INTEGER(K4B), SAVE :: iran0,jran0,kran0,nran0,mran0,rans
   INTEGER(K4B), DIMENSION(:,:), POINTER, SAVE :: ranseeds
   INTEGER(K4B), DIMENSION(:), POINTER, SAVE :: iran,jran,kran, &
      nran,mran,ranv
   REAL(SP), SAVE :: amm
   
   INTERFACE ran_hash
      MODULE PROCEDURE ran_hash_s, ran_hash_v
   END INTERFACE
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran_init(length)
      IMPLICIT NONE
      INTEGER(K4B), INTENT(IN) :: length
      INTEGER(K4B) :: new,j,hgt
      IF (length < lenran) RETURN
      hgt=hg
      IF (hg /= 2147483647) WRITE(*,*) 'Err>>ran_init: arith assump 1 fails'
      IF (hgng >= 0)        WRITE(*,*) 'Err>>ran_init: arith assump 2 fails'
      IF (hgt+1 /= hgng)    WRITE(*,*) 'Err>>ran_init: arith assump 3 fails'
      IF (NOT(hg) >= 0)     WRITE(*,*) 'Err>>ran_init: arith assump 4 fails'
      IF (NOT(hgng) < 0)    WRITE(*,*) 'Err>>ran_init: arith assump 5 fails'
      IF (hg+hgng >= 0)     WRITE(*,*) 'Err>>ran_init: arith assump 6 fails'
      IF (NOT(-1_k4b) < 0)  WRITE(*,*) 'Err>>ran_init: arith assump 7 fails'
      IF (NOT(0_k4b) >= 0)  WRITE(*,*) 'Err>>ran_init: arith assump 8 fails'
      IF (NOT(1_k4b) >= 0)  WRITE(*,*) 'Err>>ran_init: arith assump 9 fails'
      IF (lenran > 0) THEN
         ranseeds=>reallocate(ranseeds,length,5)
         ranv=>reallocate(ranv,length-1)
         new=lenran+1
      ELSE
         ALLOCATE(ranseeds(length,5))
         ALLOCATE(ranv(length-1))
         new=1
         amm=nearest(1.0_sp,-1.0_sp)/hgng
         IF (amm*hgng >= 1.0 .or. amm*hgng <= 0.0) &
            WRITE(*,*) 'Err>>ran_init: arth assump 10 fails'
      END IF
      ranseeds(new:,1)=seq
      ranseeds(new:,2:5)=spread(arth(new,1,SIZE(ranseeds(new:,1))),2,4)
      DO j=1,4
         CALL ran_hash(ranseeds(new:,j),ranseeds(new:,j+1))
      END DO
      where (ranseeds(new:,1:3) < 0) &
         ranseeds(new:,1:3)=NOT(ranseeds(new:,1:3))
      where (ranseeds(new:,4:5) == 0) ranseeds(new:,4:5)=1
      IF (new == 1) THEN
         iran0=ranseeds(1,1)
         jran0=ranseeds(1,2)
         kran0=ranseeds(1,3)
         mran0=ranseeds(1,4)
         nran0=ranseeds(1,5)
         rans=nran0
      END IF
      IF (length > 1) THEN
         iran => ranseeds(2:,1)
         jran => ranseeds(2:,2)
         kran => ranseeds(2:,3)
         mran => ranseeds(2:,4)
         nran => ranseeds(2:,5)
         ranv = nran
      END IF
      lenran=length
   END SUBROUTINE ran_init
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE rand_deallocate
      IF (lenran > 0) THEN
         DEALLOCATE(ranseeds,ranv)
         nullIFy(ranseeds,ranv,iran,jran,kran,mran,nran)
         lenran = 0
      END IF
   END SUBROUTINE rand_deallocate
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran_seed(sequence,SIZE,put,get)
      IMPLICIT NONE
      INTEGER, OPTIONAL, INTENT(IN) :: sequence
      INTEGER, OPTIONAL, INTENT(OUT) :: SIZE
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: put
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: get
      
      IF (present(SIZE)) THEN
         SIZE=5*lenran
      ELSE IF (present(put)) THEN
         IF (lenran == 0) RETURN
         ranseeds=reshape(put,shape(ranseeds))
         where (ranseeds(:,1:3) < 0) ranseeds(:,1:3)=NOT(ranseeds(:,1:3))
         where (ranseeds(:,4:5) == 0) ranseeds(:,4:5)=1
         iran0=ranseeds(1,1)
         jran0=ranseeds(1,2)
         kran0=ranseeds(1,3)
         mran0=ranseeds(1,4)
         nran0=ranseeds(1,5)
      ELSE IF (present(get)) THEN
         IF (lenran == 0) RETURN
         ranseeds(1,1:5)=(/ iran0,jran0,kran0,mran0,nran0 /)
         get=reshape(ranseeds,shape(get))
      ELSE IF (present(sequence)) THEN
         CALL rand_deallocate
         seq=sequence
      END IF
   END SUBROUTINE ran_seed
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran_hash_s(il,ir)
      IMPLICIT NONE
      INTEGER(K4B), INTENT(INOUT) :: il,ir
      INTEGER(K4B) :: is,j
      DO j=1,4
         is=ir
         ir=IEOR(ir,ISHFT(ir,5))+1422217823
         ir=IEOR(ir,ISHFT(ir,-16))+1842055030
         ir=IEOR(ir,ISHFT(ir,9))+80567781
         ir=IEOR(il,ir)
         il=is
      END DO
   END SUBROUTINE ran_hash_s
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran_hash_v(il,ir)
      IMPLICIT NONE
      INTEGER(K4B), DIMENSION(:), INTENT(INOUT) :: il,ir
      INTEGER(K4B), DIMENSION(SIZE(il)) :: is
      INTEGER(K4B) :: j
      DO j=1,4
         is=ir
         ir=IEOR(ir,ISHFT(ir,5))+1422217823
         ir=IEOR(ir,ISHFT(ir,-16))+1842055030
         ir=IEOR(ir,ISHFT(ir,9))+80567781
         ir=IEOR(il,ir)
         il=is
      END DO
   END SUBROUTINE ran_hash_v
END MODULE ran_state

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MODULE Ecuyer_random
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur
! L'Ecuyer's 1996 random number generator.
! Fortran version by Alan.Miller @ vic.cmis.csiro.au
! N.B. This version is compatible with Lahey's ELF90
! http://www.ozemail.com.au/~milleraj
! Latest revision - 30 March 1999
!L'Ecuyer's 1996 Tausworthe random number generator, and lfsr113.f90  L'Ecuyer's 1999 Tausworthe random 
!number generator. The first has a cycle of 2^88 while the second is a little slower but has a cycle of 
!2^113. Both are translations from C. N.B. These both assume that default integers are 32-bit. 
   USE VARIABLE_TYPE, ONLY : dp
   IMPLICIT NONE
   INTEGER, PARAMETER :: dprec = SELECTED_REAL_KIND(14, 60)

   ! These are unsigned integers in the C version
   INTEGER, SAVE :: s1 = 1234, s2 = -4567, s3 = 7890
   
   CONTAINS

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE InitializeTaus88(seed)
      IMPLICIT NONE
      INTEGER :: seed,seed1,i
      REAL(dp) :: random_numb
      
      seed1=ABS(seed)
      CALL init_seeds(seed1,seed1+440,seed1+123)
      DO i=1,100000
         random_numb=taus88()
      END DO
   END SUBROUTINE InitializeTaus88
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE init_seeds(i1, i2, i3)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i1, i2, i3

      s1 = i1
      s2 = i2
      s3 = i3
      IF (IAND(s1,-2) == 0) s1 = i1 - 1023
      IF (IAND(s2,-8) == 0) s2 = i2 - 1023
      IF (IAND(s3,-16) == 0) s3 = i3 - 1023

      RETURN
   END SUBROUTINE init_seeds
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION taus88() RESULT(random_numb)
   ! Generates a random number between 0 and 1.  Translated from C function in:
   ! Reference:
   ! L'Ecuyer, P. (1996) `Maximally equidistributed combined Tausworthe generators', Math. of Comput., 65, 203-213.

   ! The cycle length is claimed to be about 2^(88) or about 3E+26. Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).
      IMPLICIT NONE
      REAL (dprec) :: random_numb
      INTEGER   :: b
      
      ! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation; to the left if j > 0, otherwise to the right.
      b  = ISHFT( IEOR( ISHFT(s1,13), s1), -19)
      s1 = IEOR( ISHFT( IAND(s1,-2), 12), b)
      b  = ISHFT( IEOR( ISHFT(s2,2), s2), -25)
      s2 = IEOR( ISHFT( IAND(s2,-8), 4), b)
      b  = ISHFT( IEOR( ISHFT(s3,3), s3), -11)
      s3 = IEOR( ISHFT( IAND(s3,-16), 17), b)
      random_numb = IEOR( IEOR(s1,s2), s3) * 2.3283064365E-10_dprec + 0.5_dprec
      RETURN
   
   END FUNCTION taus88
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


END MODULE Ecuyer_random
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MODULE randomgen
   USE VARIABLE_TYPE
   USE nrprocedures
   IMPLICIT NONE
   
   INTERFACE ran1
      MODULE PROCEDURE ran1_s,ran1_v
   END INTERFACE ran1
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION ran(idum)
      IMPLICIT NONE
      INTEGER, PARAMETER :: K4B=selected_int_kind(9)
      INTEGER(K4B), INTENT(INOUT) :: idum
      REAL :: ran
      INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
      REAL, SAVE :: am
      INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
      
      IF (idum <= 0 .or. iy < 0) THEN
         am=nearest(1.0,-1.0)/IM
         iy=ior(IEOR(888889999,abs(idum)),1)
         ix=IEOR(777755555,abs(idum))
         idum=abs(idum)+1
      END IF
      ix=IEOR(ix,ISHFT(ix,13))
      ix=IEOR(ix,ISHFT(ix,-17))
      ix=IEOR(ix,ISHFT(ix,5))
      k=iy/IQ
      iy=IA*(iy-k*IQ)-IR*k
      IF (iy < 0) iy=iy+IM
      ran=am*ior(iAND(IM,IEOR(ix,iy)),1)
   END FUNCTION ran
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran0_s(harvest)
      USE ran_state, ONLY: K4B,amm,lenran,ran_init,iran0,jran0,kran0,nran0,rans
      IMPLICIT NONE
      REAL(SP), INTENT(OUT) :: harvest
      
      IF (lenran < 1) CALL ran_init(1)
      rans=iran0-kran0
      IF (rans < 0) rans=rans+2147483579_k4b
      iran0=jran0
      jran0=kran0
      kran0=rans
      nran0=IEOR(nran0,ISHFT(nran0,13))
      nran0=IEOR(nran0,ISHFT(nran0,-17))
      nran0=IEOR(nran0,ISHFT(nran0,5))
      rans=IEOR(nran0,rans)
      harvest=amm*MERGE(rans,NOT(rans), rans<0 )
   END SUBROUTINE ran0_s
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran0_v(harvest)
      USE ran_state, ONLY: K4B,amm,lenran,ran_init,iran,jran,kran,nran,ranv
      IMPLICIT NONE
      
      REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
      INTEGER(K4B) :: n
      n=SIZE(harvest)
      IF (lenran < n+1) CALL ran_init(n+1)
      ranv(1:n)=iran(1:n)-kran(1:n)
      where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
      iran(1:n)=jran(1:n)
      jran(1:n)=kran(1:n)
      kran(1:n)=ranv(1:n)
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),13))
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),-17))
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),5))
      ranv(1:n)=IEOR(nran(1:n),ranv(1:n))
      harvest=amm*MERGE(ranv(1:n),NOT(ranv(1:n)), ranv(1:n)<0 )
   END SUBROUTINE ran0_v
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran1_s(harvest)
      USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
         iran0,jran0,kran0,nran0,mran0,rans
      IMPLICIT NONE
      REAL(SP), INTENT(OUT) :: harvest
      IF (lenran < 1) CALL ran_init(1)
      rans=iran0-kran0
      IF (rans < 0) rans=rans+2147483579_k4b
      iran0=jran0
      jran0=kran0
      kran0=rans
      nran0=IEOR(nran0,ISHFT(nran0,13))
      nran0=IEOR(nran0,ISHFT(nran0,-17))
      nran0=IEOR(nran0,ISHFT(nran0,5))
      IF (nran0 == 1) nran0=270369_k4b
      mran0=IEOR(mran0,ISHFT(mran0,5))
      mran0=IEOR(mran0,ISHFT(mran0,-13))
      mran0=IEOR(mran0,ISHFT(mran0,6))
      rans=IEOR(nran0,rans)+mran0
      harvest=amm*MERGE(rans,NOT(rans), rans<0 )
   END SUBROUTINE ran1_s
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran1_v(harvest)
      USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
         iran,jran,kran,nran,mran,ranv
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
      INTEGER(K4B) :: n
      
      n=SIZE(harvest)
      IF (lenran < n+1) CALL ran_init(n+1)
      ranv(1:n)=iran(1:n)-kran(1:n)
      where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
      iran(1:n)=jran(1:n)
      jran(1:n)=kran(1:n)
      kran(1:n)=ranv(1:n)
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),13))
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),-17))
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),5))
      where (nran(1:n) == 1) nran(1:n)=270369_k4b
      mran(1:n)=IEOR(mran(1:n),ISHFT(mran(1:n),5))
      mran(1:n)=IEOR(mran(1:n),ISHFT(mran(1:n),-13))
      mran(1:n)=IEOR(mran(1:n),ISHFT(mran(1:n),6))
      ranv(1:n)=IEOR(nran(1:n),ranv(1:n))+mran(1:n)
      harvest=amm*MERGE(ranv(1:n),NOT(ranv(1:n)), ranv(1:n)<0 )
   END SUBROUTINE ran1_v
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran2_s(harvest)
      USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
         iran0,jran0,kran0,nran0,mran0,rans
      IMPLICIT NONE
      REAL(SP), INTENT(OUT) :: harvest
      
      IF (lenran < 1) CALL ran_init(1)
      rans=iran0-kran0
      IF (rans < 0) rans=rans+2147483579_k4b
      iran0=jran0
      jran0=kran0
      kran0=rans
      nran0=IEOR(nran0,ISHFT(nran0,13))
      nran0=IEOR(nran0,ISHFT(nran0,-17))
      nran0=IEOR(nran0,ISHFT(nran0,5))
      rans=iAND(mran0,65535)
      mran0=ISHFT(3533*ISHFT(mran0,-16)+rans,16)+ &
         3533*rans+820265819_k4b
      rans=IEOR(nran0,kran0)+mran0
      harvest=amm*MERGE(rans,NOT(rans), rans<0 )
   END SUBROUTINE ran2_s
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran2_v(harvest)
      USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
         iran,jran,kran,nran,mran,ranv
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
      INTEGER(K4B) :: n
      
      n=SIZE(harvest)
      IF (lenran < n+1) CALL ran_init(n+1)
      ranv(1:n)=iran(1:n)-kran(1:n)
      where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
      iran(1:n)=jran(1:n)
      jran(1:n)=kran(1:n)
      kran(1:n)=ranv(1:n)
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),13))
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),-17))
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),5))
      ranv(1:n)=iAND(mran(1:n),65535)
      mran(1:n)=ISHFT(3533*ISHFT(mran(1:n),-16)+ranv(1:n),16)+ &
         3533*ranv(1:n)+820265819_k4b
      ranv(1:n)=IEOR(nran(1:n),kran(1:n))+mran(1:n)
      harvest=amm*MERGE(ranv(1:n),NOT(ranv(1:n)), ranv(1:n)<0 )
   END SUBROUTINE ran2_v
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran3_s(harvest)
      USE ran_state, ONLY: K4B,amm,lenran,ran_init,ran_hash,mran0,nran0,rans
      IMPLICIT NONE
      REAL(SP), INTENT(OUT) :: harvest
      INTEGER(K4B) :: temp
      
      IF (lenran < 1) CALL ran_init(1)
      nran0=IEOR(nran0,ISHFT(nran0,13))
      nran0=IEOR(nran0,ISHFT(nran0,-17))
      nran0=IEOR(nran0,ISHFT(nran0,5))
      IF (nran0 == 1) nran0=270369_k4b
      rans=nran0
      mran0=IEOR(mran0,ISHFT(mran0,5))
      mran0=IEOR(mran0,ISHFT(mran0,-13))
      mran0=IEOR(mran0,ISHFT(mran0,6))
      temp=mran0
      CALL ran_hash(temp,rans)
      harvest=amm*MERGE(rans,NOT(rans), rans<0 )
   END SUBROUTINE ran3_s
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ran3_v(harvest)
      USE ran_state, ONLY: K4B,amm,lenran,ran_init,ran_hash,mran,nran,ranv
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
      INTEGER(K4B), DIMENSION(SIZE(harvest)) :: temp
      INTEGER(K4B) :: n
      
      n=SIZE(harvest)
      IF (lenran < n+1) CALL ran_init(n+1)
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),13))
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),-17))
      nran(1:n)=IEOR(nran(1:n),ISHFT(nran(1:n),5))
      where (nran(1:n) == 1) nran(1:n)=270369_k4b
      ranv(1:n)=nran(1:n)
      mran(1:n)=IEOR(mran(1:n),ISHFT(mran(1:n),5))
      mran(1:n)=IEOR(mran(1:n),ISHFT(mran(1:n),-13))
      mran(1:n)=IEOR(mran(1:n),ISHFT(mran(1:n),6))
      temp=mran(1:n)
      CALL ran_hash(temp,ranv(1:n))
      harvest=amm*MERGE(ranv(1:n),NOT(ranv(1:n)), ranv(1:n)<0 )
   END SUBROUTINE ran3_v
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION poidev(xm)
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: xm
      REAL(SP) :: poidev
      REAL(SP) :: em,harvest,t,y
      REAL(SP), SAVE :: alxm,g,oldm=-1.0_sp,sq
      
      IF (xm < 12.0) THEN
         IF (xm /= oldm) THEN
            oldm=xm
            g=exp(-xm)
         END IF
         em=-1
         t=1.0
         DO
            em=em+1.0_sp
            CALL ran1(harvest)
            t=t*harvest
            IF (t <= g) exit
         END DO
      ELSE
         IF (xm /= oldm) THEN
            oldm=xm
            sq=SQRT(2.0_sp*xm)
            alxm=LOG(xm)
            g=xm*alxm-gammln(xm+1.0_sp)
         END IF
         DO
            DO
               CALL ran1(harvest)
               y=tan(PI*harvest)
               em=sq*y+xm
               IF (em >= 0.0) exit
            END DO
            em=int(em)
            t=0.9_sp*(1.0_sp+y**2)*exp(em*alxm-gammln(em+1.0_sp)-g)
            CALL ran1(harvest)
            IF (harvest <= t) exit
         END DO
      END IF
      poidev=em
   END FUNCTION poidev
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE gasdev_s(harvest)
      IMPLICIT NONE
      REAL(SP), INTENT(OUT) :: harvest
      REAL(SP) :: rsq,v1,v2
      REAL(SP), SAVE :: g
      LOGICAL, SAVE :: gaus_stored=.FALSE.
      
      IF (gaus_stored) THEN
         harvest=g
         gaus_stored=.FALSE.
      ELSE
         DO
            CALL ran1(v1)
            CALL ran1(v2)
            v1=2.0_sp*v1-1.0_sp
            v2=2.0_sp*v2-1.0_sp
            rsq=v1**2+v2**2
            IF (rsq > 0.0 .AND. rsq < 1.0) exit
         END DO
         rsq=SQRT(-2.0_sp*LOG(rsq)/rsq)
         harvest=v1*rsq
         g=v2*rsq
         gaus_stored=.TRUE.
      END IF
   END SUBROUTINE gasdev_s
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE gasdev_v(harvest,harvest1)
      !harvest1 added to make it bivariate
      IMPLICIT NONE
      REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
      REAL(SP), DIMENSION(:), OPTIONAL :: harvest1
      REAL(SP), DIMENSION(SIZE(harvest)) :: rsq,v1,v2
      REAL(SP), ALLOCATABLE, DIMENSION(:), SAVE :: g
      INTEGER(I4B) :: n,ng,nn,m
      INTEGER(I4B), SAVE :: last_ALLOCATED=0
      LOGICAL, SAVE :: gaus_stored=.FALSE.
      LOGICAL, DIMENSION(SIZE(harvest)) :: mask
      
      n=SIZE(harvest)
      IF (PRESENT(harvest1)) THEN
         IF (SIZE(harvest)/=SIZE(harvest1)) THEN
            WRITE(*,*) "$Err>> bivariate outputs of unequal size"
            STOP
         END IF
      END IF
      
      IF (n /= last_ALLOCATEd) THEN
         IF (last_ALLOCATED /= 0) DEALLOCATE(g)
         ALLOCATE(g(n))
         last_ALLOCATED=n
         gaus_stored=.FALSE.
      END IF
      IF (gaus_stored) THEN
         harvest=g
         gaus_stored=.FALSE.
      ELSE
         ng=1
         DO
            IF (ng > n) exit
            CALL ran1(v1(ng:n))
            CALL ran1(v2(ng:n))
            v1(ng:n)=2.0_sp*v1(ng:n)-1.0_sp
            v2(ng:n)=2.0_sp*v2(ng:n)-1.0_sp
            rsq(ng:n)=v1(ng:n)**2+v2(ng:n)**2
            mask(ng:n)=(rsq(ng:n)>0.0 .AND. rsq(ng:n)<1.0)
            CALL array_copy(PACK(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
            v2(ng:ng+nn-1)=PACK(v2(ng:n),mask(ng:n))
            rsq(ng:ng+nn-1)=PACK(rsq(ng:n),mask(ng:n))
            ng=ng+nn
         END DO
         rsq=SQRT(-2.0_sp*LOG(rsq)/rsq)
         harvest=v1*rsq
         g=v2*rsq
         IF (PRESENT(harvest1)) THEN !bivariate
            harvest1=g
         END IF
         gaus_stored=.TRUE.
      END IF
   END SUBROUTINE gasdev_v
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION bnldev(pp,n)
      IMPLICIT NONE
      REAL(SP), INTENT(IN) :: pp
      INTEGER(I4B), INTENT(IN) :: n
      REAL(SP) :: bnldev
      INTEGER(I4B) :: j
      INTEGER(I4B), SAVE :: nold=-1
      REAL(SP) :: am,em,g,h,p,sq,t,y,arr(24)
      REAL(SP), SAVE :: pc,pLOG,pcLOG,en,oldg,pold=-1.0
      
      p=MERGE(pp,1.0_sp-pp, pp <= 0.5_sp )
      am=n*p
      IF (n < 25) THEN
         CALL ran1(arr(1:n))
         bnldev=count(arr(1:n)<p)
      ELSE IF (am < 1.0) THEN
         g=exp(-am)
         t=1.0
         DO j=0,n
            CALL ran1(h)
            t=t*h
            IF (t < g) exit
         END DO
         bnldev=MERGE(j,n, j <= n)
      ELSE
         IF (n /= nold) THEN
            en=n
            oldg=gammln(en+1.0_sp)
            nold=n
         END IF
         IF (p /= pold) THEN
            pc=1.0_sp-p
            pLOG=LOG(p)
            pcLOG=LOG(pc)
            pold=p
         END IF
         sq=SQRT(2.0_sp*am*pc)
         DO
            CALL ran1(h)
            y=tan(PI*h)
            em=sq*y+am
            IF (em < 0.0 .or. em >= en+1.0_sp) cycle
            em=int(em)
            t=1.2_sp*sq*(1.0_sp+y**2)*exp(oldg-gammln(em+1.0_sp)-&
               gammln(en-em+1.0_sp)+em*pLOG+(en-em)*pcLOG)
            CALL ran1(h)
            IF (h <= t) exit
         END DO
         bnldev=em
      END IF
      IF (p /= pp) bnldev=n-bnldev
   END FUNCTION bnldev
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE randomgen
