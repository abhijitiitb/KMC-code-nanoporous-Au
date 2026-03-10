! Contains subroutine files needed by AFV's makeampot.f to generate CLSMAN potential files

!PROGRAM MAIN


!   IMPLICIT NONE
!   INTEGER, PARAMETER :: dp = KIND(1.0d0)
!   REAL(dp) :: E
  
!    CALL ft(0.246_dp,E); WRITE(*,*) E
!    CALL rhot(1.35_dp,E); WRITE(*,*) E
!    CALL phit(2.536_dp,E); WRITE(*,*) E
! END PROGRAM MAIN

!Potential of 
![1] Mendelev et al, Philosophical Maganize, 83, 35, 3977-3994 - potential 2
![2] Mendelev et al, Philosophical Maganize, 83, 35, 3977-3994 - potential 4
![3] Dudarev et al, J. Phys.: Condens. Matter 17 (2005) 7097–7118 - case I
![4] Dudarev et al, J. Phys.: Condens. Matter 17 (2005) 7097–7118 - case II

   !####################################
   SUBROUTINE ft(r,E)
      IMPLICIT NONE
      INTEGER, PARAMETER :: dp = KIND(1.0d0)
      REAL(dp), PARAMETER :: A=-0.00035387096579929_dp ![1]
      !REAL(dp), PARAMETER :: A=-0.00034906178363530_dp ![2]
      !REAL(dp), PARAMETER :: A=3.527586256672234_dp, B=1.642855167616477_dp, &
      !   rc, nu= ![3]
      !REAL(dp), PARAMETER :: A=4.100199340884814_dp, B=1.565647547483517_dp, &
      !   rc, nu= ![3]
      REAL(dp) :: rho
      REAL(dp), INTENT(IN) :: r
      REAL(dp), INTENT(OUT) :: E

      rho=r
      E=A*rho**2-SQRT(rho) ![1-2]
      !E=-A*SQRT(r)-B*(SQRT(rc)-SQRT(r))**2*Heaviside(rc-r)/(nu+SQRT(rc)-SQRT(r)) [3-4]
   END SUBROUTINE ft
   !####################################
   SUBROUTINE rhot(r,E)
      IMPLICIT NONE
      INTEGER, PARAMETER :: dp = KIND(1.0d0)
      REAL(dp), PARAMETER :: f(3)=(/11.686859407970_dp,-0.014710740098830_dp,0.47193527075943_dp/) ![1-2]
      REAL(dp), PARAMETER :: rf(3)=(/2.4_dp,3.2_dp,4.2_dp/)
      REAL(dp), DIMENSION(3) :: drf,drf3
      REAL(dp), INTENT(IN) :: r
      REAL(dp), INTENT(OUT) :: E

      drf=rf-r
      drf3=drf**3*Heaviside(drf)
      E=DOT_PRODUCT(f,drf3)
      CONTAINS
      !####################################
      FUNCTION Heaviside(r)
         IMPLICIT NONE
         INTEGER, PARAMETER :: dp = KIND(1.0d0)
         REAL(dp), DIMENSION(:), INTENT(IN) :: r
         REAL(dp), DIMENSION(SIZE(r)) :: Heaviside

         Heaviside=0._dp
         WHERE (r>0._dp) Heaviside=1._dp
      END FUNCTION Heaviside
      !####################################
   END SUBROUTINE rhot
   !####################################
   SUBROUTINE phit(r,E)
   !unit is electron volt
      IMPLICIT NONE
      INTEGER, PARAMETER :: dp = KIND(1.0d0)
      REAL(dp) :: r1, r2, rs, rb, x, phix
      REAL(dp) :: B0, B1, B2, B3
      REAL(dp), PARAMETER :: v(15)=(/0.000000_dp,-24.028204854115_dp,11.300691696477_dp, &
         5.3144495820462_dp,-4.6659532856049_dp, 5.9637758529194_dp, &
         1.7710262006061_dp,0.85913830768731_dp,-2.1845362968261_dp, &
         2.6424377007466_dp,-1.0358345370208_dp,0.33548264951582_dp, &
         -0.046448582149334_dp,-0.0070294963048689_dp,0.000000000_dp/) ![1]
      REAL(dp), PARAMETER :: rv(15)=(/2.1_dp,2.2_dp,2.3_dp, &
         2.4_dp,2.5_dp,2.6_dp,2.7_dp,2.8_dp,3.0_dp, &
         3.3_dp,3.7_dp,4.2_dp,4.7_dp,5.3_dp,6.0_dp/) ![1-2]
      !REAL(dp), PARAMETER :: v(15)=(/195.92322853994_dp,17.516698453315_dp, &
      !   1.4926525164290_dp,6.4129476125197_dp,-6.8157461860553_dp, &
      !   9.6582581963600_dp,-5.3412764419_dp,1.7996558048346_dp, &
      !   -1.4788966636288_dp,1.8530435283665_dp,-0.64164344859316_dp, &
      !   0.24463630025168_dp,-0.057721650527383_dp,0.023358616514826_dp, &
      !   -0.0097064921265079_dp/) ![2] - 15th value is doubtful -- it appears there is a typo in the main paper
      REAL(dp), DIMENSION(15) :: drv,drv3
      REAL(dp), INTENT(IN) :: r
      REAL(dp), INTENT(OUT) :: E
      
      r1=1._dp
      r2=2._dp
      rb=0.529_dp
      rs=0.88534_dp*rb/(SQRT(2._dp)*26._dp**(1./3.))
      

      IF (r<=r1) THEN
         x=r/rs
         phix=0.1818_dp*exp(-3.2_dp*x)+0.5099_dp*exp(-0.9423_dp*x)+ &
            0.2802_dp*exp(-0.4029_dp*x)+0.02817_dp*exp(-0.2016_dp*x) ![1-4]
         E=26._dp**2*14.3992_dp*phix/r ![1-4]
      ELSEIF (r>r1 .AND. r<r2) THEN
         B0=6.4265260576348_dp ![1]
         B1=1.7900488524286_dp ![1]
         B2=-4.5108316729807_dp ![1]
         B3=1.0866199373306_dp ![1]
         !B0=14.996917289290_dp ![2]
         !B1=-20.533174190155_dp ![2]
         !B2=14.002591780752_dp ![2]
         !B3=-3.6473736591143_dp ![2]
         E=exp(B0+B1*r+B2*r*r+B3*r*r*r) ![1-2]
      ELSE ! (r>r2)
         drv=rv-r
         drv3=drv**3*Heaviside(drv)
         E=DOT_PRODUCT(v,drv3)
      END IF
      
      CONTAINS
      !####################################
      FUNCTION Heaviside(r)
         IMPLICIT NONE
         INTEGER, PARAMETER :: dp = KIND(1.0d0)
         REAL(dp), DIMENSION(:), INTENT(IN) :: r
         REAL(dp), DIMENSION(SIZE(r)) :: Heaviside

         Heaviside=0._dp
         WHERE (r>0._dp) Heaviside=1._dp
      END FUNCTION Heaviside
      !####################################
   END SUBROUTINE phit
