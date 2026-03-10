MODULE NMode
   USE VARIABLE_TYPE
   USE PotentialPackage
   USE NEB_package
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE HTST(AL1,AL2,ntransitionstates,TSEnergy,prefactor)
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL1,AL2,AL3
      INTEGER :: ntransitionstates,errorstatus
      REAL(dp) :: TSEnergy,springconstant
      REAL(dp), OPTIONAL :: prefactor
      
      !we assume that VL of AL1 and AL2 are associated
      
      springconstant=0.5_dp
      ntransitionstates=0
      TSEnergy=0._dp
      IF (ASSOCIATED(AL3)) NULLIFY(AL3)
      
      CALL NEB(AL1,AL2,nimages=15,imode=1101,omode=13,interpmode=1, &
         springconstant=springconstant,ITMAX=50,ntransitionstates=ntransitionstates, &
         TSEnergy=TSEnergy,TS=AL3,errorstatus=errorstatus)
      !prefactor=1.e13_dp
      IF (PRESENT(prefactor)) THEN
         prefactor=0._dp
         IF (errorstatus==0 .AND. ntransitionstates==1) THEN
            CALL PotentialInitialize(AL3)
            CALL WritePotential(AL3)
            CALL AddVerletList(AL3,NAtoms=AL3%NAtoms)
            CALL Vineyard(AL1,AL3,prefactor)
         ELSE
            TSEnergy=-9999._dp
         END IF
      END IF
      CALL Delete(AL3) !delete transition state
   END SUBROUTINE HTST
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Vineyard(AL1,AL2,prefactor)
   !calculates the pre-exponential factor given 
   !the min energy state AL1 and transition state AL2
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL1,AL2
      REAL(dp) :: prefactor,lnfreqprod1,lnfreqprod2
      INTEGER :: NAtoms,nnegative1,nnegative2,nmodes
      REAL(dp), DIMENSION(:), ALLOCATABLE :: eigen1,eigen2
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: eigenvector
      
      NAtoms=AL1%NAtoms
      
      !IF (NAtoms/=AL2%NAtoms) THEN
      !   WRITE(6,*) "$Err>> Number of atoms is not same in initial and transition state"
      !   STOP
      !END IF
      
      !initialize
      ALLOCATE(eigen1(3*NAtoms))
      ALLOCATE(eigen2(3*NAtoms))
      
      !find the normal modes for initial state (AL1)
      nmodes=0
      ALLOCATE(eigenvector(1,3*NAtoms))
      CALL NormalModeAnalysis(AL1,NAtoms,eigen1,nnegative1,lnfreqprod1,eigenvector,0)
      IF (nnegative1/=0) THEN
         WRITE(6,*) TRIM(TxtHeader)//"Err>> Number of negative eigenvalues for minimum energy state:",nnegative1
         !STOP
      END IF
      WRITE(6,*) TRIM(TxtHeader)//"Number of negative eigenvalues @ min:",nnegative1
      WRITE(6,*) TRIM(TxtHeader)//"Log product eigenvalues:",lnfreqprod1
      CALL FLUSH(6)
      

      !find the normal modes for transition state (AL2)
      CALL NormalModeAnalysis(AL2,NAtoms,eigen2,nnegative2,lnfreqprod2,eigenvector,0)
      IF (nnegative2/=1) THEN
         WRITE(6,*) TRIM(TxtHeader)//"Err>> Number of negative eigenvalues for transition state:",nnegative2
         !STOP
      END IF
      WRITE(6,*) TRIM(TxtHeader)//"Number of negative eigenvalues @ TS:",nnegative2
      WRITE(6,*) TRIM(TxtHeader)//"Log product eigenvalues:",lnfreqprod2
      
      !find the prefactor
      prefactor=EXP(lnfreqprod1-lnfreqprod2) !units of s^-1
      WRITE(UNIT=6,FMT='("   Prefactor:",ES15.5,"s-1")') prefactor
      IF (nnegative2-nnegative1/=1) THEN !something is odd
         WRITE(6,*) TRIM(TxtHeader)//" Err>> Number of transistion states odd",nnegative2,nnegative1
         WRITE(6,*) TRIM(TxtHeader)//"   Prefactor:",prefactor,"s-1"
         prefactor=0._dp
      END IF
      
      !finalize
      DEALLOCATE(eigen1,eigen2)
      DEALLOCATE(eigenvector)
   END SUBROUTINE Vineyard
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE NormalModeAnalysis(AL,NAtoms,eigenvalues,nnegative,lnfreqproduct,eigenvector,nmodes)
   !finds the normal modes and returns the lnfreqproduct
   !eigenvalues are always computed
   !nmodes is the number of largest modes that need to be returned
      IMPLICIT NONE
      INTEGER, PARAMETER :: maxnmove=1000 !size of curvature is 3maxnmove*(3maxnmove+1)
      INTEGER, PARAMETER :: maxndimn=(3*maxnmove*(3*maxnmove+1))/2
      INTEGER, DIMENSION(maxnmove) :: AtomMap
      REAL(dp), DIMENSION(maxndimn) :: curvature
      DOUBLE PRECISION, DIMENSION(9*maxnmove) :: workarray
      TYPE(SystemContainer), POINTER :: AL
      REAL(dp) :: lnfreqproduct,fconv,eig
      DOUBLE PRECISION, DIMENSION(:) :: eigenvalues
      INTEGER :: idimn,ndimn,info,i,ListType,NAtoms,nnegative,nmove
      INTEGER, OPTIONAL :: nmodes
      DOUBLE PRECISION, DIMENSION(:,:), OPTIONAL :: eigenvector !3*NAtoms*nmodes elements
      
      CALL GetHessianVL(AL,NAtoms,curvature,ndimn,AtomMap,nmove)
      !nmove is number of moving atoms, 3Nmove is order of Hessian, ndimn=(3Nmove*(3Nmove+1))/2
      
!SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
!*  DSPEV computes all the eigenvalues and, optionally, eigenvector of a real symmetric matrix A in packed storage.
!*  JOBZ    (input) CHARACTER*1 - 'N':  Compute eigenvalues only; 'V':  Compute eigenvalues and eigenvector.
!*  UPLO    (input) CHARACTER*1 - 'U':  Upper triangle of A is stored; 'L':  Lower triangle of A is stored.
!*  N       (input) INTEGER - The order of the matrix A.  N >= 0.
!*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!*  On entry, the upper or lower triangle of the symmetric matrix A, packed columnwise in a linear array.  The j-th column of A
!* is stored in the array AP as follows:
!*  if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!*  On exit, AP is overwritten by values generated during the reduction to tridiagonal form. 
!If UPLO = 'U', the diagonal and first superdiagonal of the tridiagonal matrix T overwrite
! the corresponding elements of A, and if UPLO = 'L', the diagonal and first subdiagonal of T overwrite the corresponding elements of A.
!*  W       (output) DOUBLE PRECISION array, dimension (N) If INFO = 0, the eigenvalues in ascending order.
!*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
! If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal eigenvector of the matrix A, with the i-th column of Z
! holding the eigenvector associated with W(i). If JOBZ = 'N', then Z is not referenced.
!*  LDZ     (input) INTEGER - The leading dimension of the array Z.  LDZ >= 1, and if JOBZ = 'V', LDZ >= max(1,N).
!*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
!*  INFO    (output) INTEGER = 0:  successful exit. < 0:  if INFO = -i, the i-th argument had an illegal value.
! > 0:  if INFO = i, the algorithm failed to converge; i off-diagonal elements of an intermediate tridiagonal form did not converge to zero.
      !IF (nmodes>0) THEN
      !   CALL DSPEV('V','U',3*nmove,curvature,eigenvalues(1:3*nmove), &
      !      eigenvector(nmodes,1:3*nmove),nmodes,workarray(1:9*nmove),info)
      !ELSEIF (nmodes==0) THEN
!         CALL DSPEV('N','U',3*nmove,curvature,eigenvalues(1:3*nmove), &
!            eigenvector(1,1:3*nmove),1,workarray(1:9*nmove),info)
      !END IF
      
      !convert diagonal elements to frequencies
      nnegative=0
      !nproduct=0
      fconv=SQRT(fev/famu)/fang !conversion factor to obtain s^-1 units
      lnfreqproduct=0._dp
      DO idimn=1,3*nmove
         eig=eigenvalues(idimn)
         IF (eig<0.) nnegative=nnegative+1
         eig=SQRT(ABS(eig))*SIGN(1._dp,eig)/TWOPI
         eig=eig*fconv
         eigenvalues(idimn)=eig
         IF (eig>0._dp) lnfreqproduct=lnfreqproduct+LOG(eig)
      END DO
      
   END SUBROUTINE NormalModeAnalysis
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE NMode
