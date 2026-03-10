MODULE EAMGenerateModule
   
!Art Voter's makeampot.f


! makapot - Self-contained example of how to create the 3 interpolation arrays
!           needed for an EAM potential,
!           and write it out in the clsman trpman form,   afv 3/5/09

!           The functions phi, rho, and f are functions that the user must
!           write to define the shapes of phi, rho, and the embedding function f.

!           Also be sure to set the desired values for xstart and xstop for each fn.

!           See the portapot package to understand the file formats better.

!           NOTE: each of the functions must be declared external in the main,
!           so that its name can be passed into maktrp and called from there.
!           If you make different function names, be sure to declare it external, 
!           and put your new name into the calling list for maktrp as well.
!           The function pairsub must have just one parameter in the call:  r.
!
!           The length unit is angstroms, and the energy returned from phi and f
!           has to be in hartress. (you can convert from eV using 27.21 eV/hartree, 
!           which will match the internal conversion in clsman)

!           Note that you can use "mkstrp" instead of "maktrp" to have cutoff smoothing
!           done automatically.  You can find mkstrp, and other routines you might
!           find useful, in the trpsubs.f package, which has been included in its
!           totality here.

!           This run will make both binary and text versions of the files. 
!           They will be named "x.phi", x.rho, x.f ...
!           The text versions will be named "x.phit", x.rhot, ...

!           To convert the text files to binary for a particular machine,
!           which is necessary to use them as input for clsman or tad, 
!           use the "potcong" program, which is in the portapot package.

!  based on auercolessi.f, written in 2006

   INTEGER, PARAMETER, PRIVATE :: dp = KIND(1.0d0)
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE EAMGenerate()
      IMPLICIT NONE
      INTEGER :: npt,iflo,ifhi,iunit,itxt,j
      REAL(dp) :: xstart,xstop,dx,E1,E2,E3,E4,r1,r2,r3,r4,delr
      DOUBLE PRECISION, DIMENSION(100000) :: work

      npt=2000 !W - Dudarev
      npt=3000 !Fe- Mendelev

! phi array

      !Tungsten - Dudarev
      xstart=0.1d0   ! not sure what to use here.
      xstop=4.2689d0   ! from Dudarev et al paper, Tungsten phi cutoff
      !Iron - Mendelev
      xstart=0.1
      xstop=5.307d0
      iflo=1
      ifhi=0

      iunit=200
      !open(unit=iunit,file='x.phit',form='unformatted')
      open(unit=iunit,file="x.phit")
      !call maktrp(work,npt,xstart,xstop, phi ,iflo,ifhi,iunit)
      r1=xstart
      dx=(xstop-xstart)/REAL(npt,dp)
      delr=1./dx
      r2=xstart+dx
      r3=xstart+2.*dx
      r4=xstart+3.*dx
      WRITE(iunit,*) npt,xstart,xstop,delr 
      
      DO j=1,npt,4
         E1=phi(r1) !eV
         E2=phi(r2) !ev 
         E3=phi(r3) !ev 
         E4=phi(r4) !ev 
         WRITE(UNIT=iunit,FMT='(4E20.12)') E1,E2,E3,E4
         !WRITE(6,*) r1,E1
         !WRITE(6,*) r2,E2
         !WRITE(6,*) r3,E3
         !WRITE(6,*) r4,E4
         r1=r1+4.*dx
         r2=r2+4.*dx
         r3=r3+4.*dx
         r4=r4+4.*dx
      END DO
      close(unit=iunit)

! rho array
      !Tungsten - Dudarev
      xstart=0.1d0   ! not sure what to use here.
      xstop=4.2689d0    ! Here we are using the longer of the two (they have 3.3216)
      !Iron - Mendelev
      xstart=0.1
      xstop=5.307d0
      iflo=1
      ifhi=0

      iunit=71
      !open(unit=iunit,file='x.rhot',form='unformatted')
      open(unit=iunit,file='x.rhot')
      r1=xstart
      dx=(xstop-xstart)/REAL(npt,dp)
      delr=1./dx
      r2=xstart+dx
      r3=xstart+2.*dx
      r4=xstart+3.*dx
      WRITE(iunit,*) npt, xstart,xstop,delr 
      DO j=1,npt,4
         E1=rho(r1) !eV
         E2=rho(r2) !ev 
         E3=rho(r3) !ev 
         E4=rho(r4) !ev 
         WRITE(UNIT=iunit,FMT='(4E20.12)') E1,E2,E3,E4
         !WRITE(6,*) r1,E1
         !WRITE(6,*) r2,E2
         !WRITE(6,*) r3,E3
         !WRITE(6,*) r4,E4
         r1=r1+4.*dx
         r2=r2+4.*dx
         r3=r3+4.*dx
         r4=r4+4.*dx
      END DO
      close(unit=iunit)

! F array

      !xstart=0.0d0
      !xstop=4000.0d0   ! not much idea yet how big to make this
      xstart=0.d0
      xstop=1000.d0
      iflo=0
      ifhi=1

      iunit=91
      !open(unit=iunit,file='x.ft',form='unformatted')
      open(unit=iunit,file='x.ft')
      !call maktrp(work,npt,xstart,xstop, f ,iflo,ifhi,iunit)
      r1=xstart
      dx=(xstop-xstart)/REAL(npt,dp)
      delr=1./dx
      r2=xstart+dx
      r3=xstart+2.*dx
      r4=xstart+3.*dx
      WRITE(iunit,*) npt, xstart,xstop,delr 
      DO j=1,npt,4
         E1=f(r1) !eV
         E2=f(r2) !ev 
         E3=f(r3) !ev 
         E4=f(r4) !ev 
         WRITE(UNIT=iunit,FMT='(4E20.12)') E1,E2,E3,E4
         WRITE(6,*) r1,E1
         WRITE(6,*) r2,E2
         WRITE(6,*) r3,E3
         WRITE(6,*) r4,E4
         r1=r1+4.*dx
         r2=r2+4.*dx
         r3=r3+4.*dx
         r4=r4+4.*dx
      END DO
      close(unit=iunit)
      STOP


! Now make text (instead of binary) versions of these three files:

      !itxt=19    ! unit for this new file
      !open(unit=itxt,file='x.phit',status='unknown',form='formatted')
      !iunit=81
      !call trptxo(iunit,itxt)


      !itxt=20    ! unit for this new file
      !open(unit=itxt,file='x.rhot',status='unknown',form='formatted')
      !iunit=71
      !call trptxo(iunit,itxt)


      !itxt=21    ! unit for this new file
      !open(unit=itxt,file='x.ft',status='unknown',form='formatted')
      !iunit=91
      !call trptxo(iunit,itxt)



      !write(6,*) 'EAM potential files ready, both binary and text'
      !write(6,*) 'You must separately create a .doc file'
      WRITE(6,*) "EAM tables have been generated"

   END SUBROUTINE EAMGenerate



!====================================================
!HERE are the functions the user must re-write:
      !INCLUDE "PotentialEAMTungsten.f90" !Tungsten potential - radiation damage
      INCLUDE "PotentialEAMFe.f90" !Fe potential - radiation damage
      
      
      FUNCTION phi(r)
      IMPLICIT NONE
      DOUBLE PRECISION :: phi,phireturned,r

      call phit(r,phireturned)

      phi=phireturned

      RETURN
      END FUNCTION phi


      FUNCTION rho(r)
      IMPLICIT NONE
      DOUBLE PRECISION :: rho,rhoreturned,r

      call rhot(r,rhoreturned)
      rho=rhoreturned

      return
      END FUNCTION rho


      FUNCTION f(rhobar)
      IMPLICIT NONE
      DOUBLE PRECISION :: f,freturned,rhobar

      call ft(rhobar,freturned)
      f=freturned  ! already in hartree

      return
      END FUNCTION f

!c====================================================

  
END MODULE EAMGenerateModule
