cMODULE MEAMInitialize
c   USE VARIABLE_TYPE
c   IMPLICIT NONE
c   REAL(dp), PRIVATE :: Eu,Ee,scr1,scr2
c   CHARACTER(len=10), PRIVATE :: reference
c   TYPE(SystemContainer), POINTER :: AL1
c   CONTAINS
c   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE Define_EAM_Potential(Ncomp,LineTable,Ierr)

c!  Define Pairwise Potential and Make Tables

      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      common /lj_1/ PhiTab(350002,5,5),DPhiTab(350002,5,5)
      common /lj_2/ Rmin,Rsqmin,Rcutoff,OverDeltaRsq
      common /eam1/ Ec(10,10),Re(10),B(10,10),A(10),beta(10,4),t(10,3)
      common /eam2/ RhoZ(10),Rho0_Bar(10),Z(10),Z2(10),A2(10)
      common /eam3/ Cmin(9,9,9),Cmax(9,9,9),Dcutoff,OverDeltaR,Eu_d(10)
      common /eam4/ Re_Alloy(10,10),alpha(10,10),omega(10,10)
      common /ener/ ene_pot(200000),ene_kin(200000),Amass(5)
      common /scrn/ scr_2nn
      common /comp/ Icomp(5)
      dimension data_ele(14)
      character Icomp*2,element*2,card*80,card2*80
      character lattice*6,lattice_element(10)*6
      character alloy*5,alloy_1_name*5,alloy_2_name*5
      character alloy3*8,ternary_name*8

c!  ===========================================================================
c!     Open EAM_DB, search for the elements and Read MEAM parameters
c!  ===========================================================================

      open(unit=2,file='eam_db.tdb',status='old')
      do n = 1, Ncomp
         read(2,'(a)') card
         read(2,'(a)') card2
    1    continue
         read(2,'(a)',end=4444) card
         read(2,'(a)') card2
         read(card,'(a2)') element
         if(element.ne.Icomp(n)) goto 1
         if(element.eq.'XX') goto 4444
         read(card,2000) lattice,amass(n),(data_ele(k), k=1,6)
 2000    format(3x,a6,d12.3,d9.2,d10.3,d11.4,3d9.2)
         read(card2,2001) Eu_d(n),(data_ele(k), k=7,14)
 2001    format(2x,d10.3,2d9.2,3d10.2,d9.2,x,2f4.2)
c
         Ec(n,n) = Data_Ele(1)
         Re(n)   = Data_Ele(2)
         B(n,n)  = Data_Ele(3)
         A(n)    = Data_Ele(4)
         beta(n,4) = Data_Ele(5)
         beta(n,1) = Data_Ele(6)
         beta(n,2) = Data_Ele(7)
         beta(n,3) = Data_Ele(8)
         t(n,1) = Data_Ele(9)
         t(n,2) = Data_Ele(10)
         t(n,3) = Data_Ele(11)
         RhoZ(n)   = Data_Ele(12)
         Cmin(n,n,n) = Data_Ele(13)
         Cmax(n,n,n) = Data_Ele(14)
         lattice_element(n) = lattice
c
         rewind(2)
c
c
c ---------------This Part should be further extended -------------------------
c
         if(lattice.eq.'FCC_A1') then
c
            a2(n)   = dsqrt(2.d0)
            Z(n)    = 12.0d0
            scr_2nn = screen_on_2nn(lattice,Cmax(n,n,n),Cmin(n,n,n))
            Z2(n) = 6.0d0 * scr_2nn
            ZR2 = Z2(n) / Z(n)
            alat = dsqrt(2.d0) * Re(n)
            omega(n,n) = alat * alat * alat / 4.d0
c
         elseif(lattice.eq.'BCC_A2') then
c
            a2(n)   = dsqrt(4.d0/3.d0)
            Z(n)    = 8.d0
            scr_2nn = screen_on_2nn(lattice,Cmax(n,n,n),Cmin(n,n,n))
            Z2(n) = 6.0d0 * scr_2nn
            ZR2 = Z2(n) / Z(n)
            alat = 2.d0 * Re(n) / dsqrt(3.d0)
            omega(n,n) = alat * alat * alat / 2.d0
c
         elseif(lattice.eq.'HCP_A3') then
c
            a2(n)   = dsqrt(2.d0)
            Z(n)    = 12.0d0
            scr_2nn = screen_on_2nn(lattice,Cmax(n,n,n),Cmin(n,n,n))
            Z2(n) = 6.0d0 * scr_2nn
            ZR2 = Z2(n) / Z(n)
            alat = Re(n)
            omega(n,n) = alat * alat * alat / a2(n)
c
         elseif(lattice.eq.'DIA_A4') then
c
            a2(n)   = 4.d0 / dsqrt(6.d0)
            Z(n)    = 4.0d0
            scr_2nn = screen_on_2nn(lattice,Cmax(n,n,n),Cmin(n,n,n))
            Z2(n) = 12.0d0 * scr_2nn
            ZR2 = Z2(n) / Z(n)
            alat = 4.d0 * Re(n) / dsqrt(3.d0)
            omega(n,n) = alat * alat * alat / 8.d0
c
         elseif(lattice.eq.'DIMER') then
c  This is an arbitrary treatment, especially those for lat and omega(n,n)
            a2(n)   = dsqrt(2.d0)
            Z(n)    = 1.0d0
            scr_2nn = 0.d0
            Z2(n) = scr_2nn
            ZR2 = Z2(n) / Z(n)
            alat = Re(n)
            omega(n,n) = alat * alat * alat
c
         else
            goto 4445
         endif
c
         alpha(n,n) = dsqrt(9.d0*B(n,n)*omega(n,n)/Ec(n,n))
         if(lattice.eq.'DIMER') alpha(n,n) = B(n,n)
         alphaN = alpha(n,n) 
         write(*,'(1x,a2,a9,f7.4)') Icomp(n),': alpha =',alphaN
         a2N = a2(n)
         R2factor = a2(n) - 1.d0
         Rho0_Bar(n) = RhoZ(n) * Z(n)
     &       + Z2(n) * RhoZ(n) * dexp ( - beta(n,4) * R2factor )
c
c  add t(3) term in the case of HCP_A3 or DIA_A4 lattice
c
         if(lattice.eq.'HCP_A3' .or. lattice.eq.'DIA_A4') then
            if(lattice.eq.'HCP_A3') then
               Rho3 = 0.5d0 * RhoZ(n) 
     &    - scr_2nn * dsqrt(2.d0) * RhoZ(n) * dexp(-beta(n,3)*R2factor)
               Rho3P = 4.d0 * Rho3 * Rho3 / 3.d0
            elseif(lattice.eq.'DIA_A4') then
               Rho3 = RhoZ(n)
               Rho3P = 32.d0 * Rho3 * Rho3 / 9.d0
            endif
c
            Rho0 = Rho0_Bar(n)
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.d0) RhoFactor = 1.d-20
            GAMMA = t(n,3) * Rho3P / RhoFactor
            Rho0_Bar(n) = Rho0 * 2.d0 / ( 1.d0 + dexp(-GAMMA) )
c           Rho0_Bar(n) = Rho0 * dsqrt( 1.d0 + GAMMA )
         endif
c
c  add t(1), t(2) and t(3) term in the case of DIMER elements (Gases)
c
         if(lattice.eq.'DIMER') then
            Rho1 = RhoZ(n)
            Rho1P = Rho1 * Rho1
            Rho2 = RhoZ(n)
            Rho2P = 2.d0 / 3.d0 * Rho2 * Rho2
            Rho3 = RhoZ(n)
            Rho3P = 2.d0 / 5.d0 * Rho3 * Rho3
c
            Rho0 = Rho0_Bar(n)
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.d0) RhoFactor = 1.d-20
            GAMMA = (t(n,1)*Rho1P+t(n,2)*Rho2P+t(n,3)*Rho3P) / RhoFactor
            Rho0_Bar(n) = Rho0 * 2.d0 / ( 1.d0 + dexp(-GAMMA) )
c           Rho0_Bar(n) = Rho0 * dsqrt( 1.d0 + GAMMA )
         endif
c
c -----------------------------------------------------------------------------
c
c     Make a table for pair potential and its derivative between each pair
c
         DeltaR = ( Rcutoff - RMin ) / dble( LineTable - 1 )
         OverDeltaR = 1.d0 / DeltaR
c
         do k=1,LineTable
            Dij = RMin + (k-1) * DeltaR
            call Define_Pair_Potential
     &          (N,lattice,alphaN,ZR2,a2N,Dij,Phi,DPhi)
            PhiTab(k,n,n)  = Phi
            DPhiTab(k,n,n) = DPhi
         enddo
      enddo
c
c  ===========================================================================
c     Search for the alloy parameters
c  ===========================================================================
c
      if(Ncomp.gt.1) then
      do i = 1, Ncomp-1
      do j = i+1, Ncomp
         alloy_1_name(1:2) = Icomp(i)         
         alloy_1_name(3:3) = '-'
         alloy_1_name(4:5) = Icomp(j)         
         alloy_2_name(1:2) = Icomp(j)         
         alloy_2_name(3:3) = '-'
         alloy_2_name(4:5) = Icomp(i)         
c
         read(2,'(a)') card
         read(2,'(a)') card2
    2    continue
         read(2,'(a)',end=4446) card
         read(card,'(a5)') alloy
         if(alloy.ne.alloy_1_name .and. alloy.ne.alloy_2_name) goto 2
         read(card,2010) lattice,delta_Ec,Re_input,B_ij,d_input
     &                 , C_iji, C_jij, C_iij, C_ijj
         read(2,2012) Cx_iji, Cx_jij, Cx_iij, Cx_ijj
c2010    format(7x,a6,5d12.4)
 2010    format(7x,a6,3d12.4,1X,5f5.2)
 2012    format(55x,4f5.2)
c
         rewind(2)
c
         if(C_iji.eq.0.0) then
            Cmin(i,j,i) = Cmin(i,i,i)
         else
            Cmin(i,j,i) = C_iji
         endif
         if(C_jij.eq.0.0) then
            Cmin(j,i,j) = Cmin(j,j,j)
         else
            Cmin(j,i,j) = C_jij
         endif
         if(Cx_iji.ne.0.0) Cmax(i,j,i) = Cx_iji
         if(Cx_jij.ne.0.0) Cmax(j,i,j) = Cx_jij
         if(Cx_iij.ne.0.0) then
            Cmax(i,i,j) = Cx_iij
            Cmax(j,i,i) = Cx_iij
         endif
         if(Cx_ijj.ne.0.0) then
            Cmax(i,j,j) = Cx_ijj
            Cmax(j,j,i) = Cx_ijj
         endif
c
c ---------------This Part should be further extended -------------------------
c
         if(lattice.eq.'FCC_B1' .or. lattice.eq.'NaCl') then
c
c          for fcc_B1/NaCl
c
            a2N      = dsqrt(2.d0)
            Ec(i,j)    = .5d0 * (Ec(i,i) + Ec(j,j)) - delta_Ec
            if(Re_input.eq.0.0d0) then
               omega(i,j) = .5d0 * (omega(i,i) + omega(j,j))
                   tmp_ij = 8.d0 * omega(i,j)
               Re_ij      = tmp_ij**(1.d0/3.d0) / 2.d0
            else
               Re_ij      =  Re_input
                   tmp_ij =  2.d0 * Re_ij
               omega(i,j) = .125d0 * tmp_ij * tmp_ij * tmp_ij 
            endif
            if(B_ij.eq.0.0d0) then
               B(i,j) = .5d0 * B(i,i) + .5d0 * B(j,j)
            else
               B(i,j) =  B_ij
            endif
            alpha_ij  = dsqrt(9.d0*B(i,j)*omega(i,j)/Ec(i,j))

            if(d_input.eq.0.0d0) then
               d_ij   = .5d0 * Eu_d(i) + .5d0 * Eu_d(j)
            elseif(d_input.lt.0.0d0) then
               d_ij   =  0.0d0  
            else
               d_ij   =  d_input
            endif

            Z_ij      = 6.0d0
              scr_2nn = screen_on_2nn(lattice,Cmax(i,j,i),Cmin(i,j,i))
            Z2_i      = 12.0d0 * scr_2nn
              scr_2nn = screen_on_2nn(lattice,Cmax(j,i,j),Cmin(j,i,j))
            Z2_j      = 12.0d0 * scr_2nn

            if(C_iij.eq.0.0) then
               tmp_ij = 0.5d0*(dsqrt(Cmin(i,i,i))+dsqrt(Cmin(j,j,j)))
               Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               Cmin(i,i,j) = C_iij
            endif
            Cmin(j,i,i) = Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij = 0.5d0*(dsqrt(Cmin(i,i,i))+dsqrt(Cmin(j,j,j)))
               Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               Cmin(i,j,j) = C_ijj
            endif
            Cmin(j,j,i) = Cmin(i,j,j)
c
         elseif(lattice.eq.'MC_BCC') then
c
c          for BCC M-C Bianry 
c
            a2N     = dsqrt(2.d0)
            Ec(i,j) = .25d0 * Ec(i,i) + .75d0 * Ec(j,j) - delta_Ec
            if(Re_input.eq.0.0d0) then
               omega(i,j) = .25d0 * omega(i,i) + .75d0 * omega(j,j)
                   tmp_ij = omega(i,j)
               Re_ij      = tmp_ij**(1.d0/3.d0)
            else
               Re_ij      =  Re_input
               omega(i,j) =  Re_ij * Re_ij * Re_ij 
            endif
            if(B_ij.eq.0.0d0) then
               B(i,j) = .25d0 * B(i,i) + .75d0 * B(j,j)
            else
               B(i,j) =  B_ij
            endif
            alpha_ij = dsqrt(9.d0*B(i,j)*omega(i,j)/Ec(i,j))

            if(d_input.eq.0.0d0) then
               d_ij   = .25d0 * Eu_d(i) + .75d0 * Eu_d(j)
            elseif(d_input.lt.0.0d0) then
               d_ij   =  0.0d0  
            else
               d_ij   =  d_input
            endif

            Z_ij      = 6.0d0
              scr_2nn = screen_on_2nn(lattice,Cmax(i,j,i),Cmin(i,j,i))
            Z2_i      = 12.0d0 * scr_2nn
c  Here, Cmin(j,i,j) should be equal to Cmin(j,j,j) .. 2002. 6. 14, B.-J. Lee
              scr_2nn = screen_on_2nn(lattice,Cmax(j,j,j),Cmin(j,j,j))
            Z2_j      = 12.0d0 * scr_2nn

            if(C_iij.eq.0.0) then
               tmp_ij =.5d0*dsqrt(Cmin(i,i,i))+.5d0*dsqrt(Cmin(j,j,j))
               Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               Cmin(i,i,j) = C_iij
            endif
            Cmin(j,i,i) = Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij =.5d0*dsqrt(Cmin(i,i,i))+.5d0*dsqrt(Cmin(j,j,j))
               Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               Cmin(i,j,j) = C_ijj
            endif
            Cmin(j,j,i) = Cmin(i,j,j)
c
         elseif(lattice.eq.'BCC_B2') then
c
c        for bcc_B2
c
            a2N     = dsqrt(4.d0/3.d0)
            Ec(i,j) = .5d0 * (Ec(i,i) + Ec(j,j)) - delta_Ec
            if(Re_input.eq.0.0d0) then
               omega(i,j) = .5d0 * (omega(i,i) + omega(j,j))
                   tmp_ij = 2.d0 * omega(i,j)
               Re_ij      = tmp_ij**(1.d0/3.d0) / a2N
            else
               Re_ij      =  Re_input
                   tmp_ij =  a2N * Re_ij
               omega(i,j) = .5d0 * tmp_ij * tmp_ij * tmp_ij 
            endif
            if(B_ij.eq.0.0d0) then
               B(i,j) = .5d0 * B(i,i) + .5d0 * B(j,j)
            else
               B(i,j) =  B_ij
            endif
            alpha_ij = dsqrt(9.d0*B(i,j)*omega(i,j)/Ec(i,j))

            if(d_input.eq.0.0d0) then
               d_ij   = .5d0 * Eu_d(i) + .5d0 * Eu_d(j)
            elseif(d_input.lt.0.0d0) then
               d_ij   =  0.0d0  
            else
               d_ij   =  d_input
            endif

            Z_ij      = 8.0d0
              scr_2nn = screen_on_2nn(lattice,Cmax(i,j,i),Cmin(i,j,i))
            Z2_i      = 6.0d0 * scr_2nn
              scr_2nn = screen_on_2nn(lattice,Cmax(j,i,j),Cmin(j,i,j))
            Z2_j      = 6.0d0 * scr_2nn

            if(C_iij.eq.0.0) then
               tmp_ij = 0.5d0*(dsqrt(Cmin(i,i,i))+dsqrt(Cmin(j,j,j)))
               Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               Cmin(i,i,j) = C_iij
            endif
            Cmin(j,i,i) = Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij = 0.5d0*(dsqrt(Cmin(i,i,i))+dsqrt(Cmin(j,j,j)))
               Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               Cmin(i,j,j) = C_ijj
            endif
            Cmin(j,j,i) = Cmin(i,j,j)
c
         elseif(lattice.eq.'L12AB3') then
c
c        for fcc_L12 (AB3 type)
c
            a2N     = dsqrt(2.d0)
            Ec(i,j) = .25d0 * Ec(i,i) + .75d0 * Ec(j,j) - delta_Ec
            if(Re_input.eq.0.0d0) then
               omega(i,j) = .25d0 * omega(i,i) + .75d0 * omega(j,j)
                   tmp_ij = 4.d0 * omega(i,j)
               Re_ij      = tmp_ij**(1.d0/3.d0) / a2N
            else
               Re_ij      =  Re_input
                   tmp_ij =  a2N * Re_ij
               omega(i,j) = .25d0 * tmp_ij * tmp_ij * tmp_ij 
            endif
            if(B_ij.eq.0.0d0) then
               B(i,j) = .25d0 * B(i,i) + .75d0 * B(j,j)
            else
               B(i,j) =  B_ij
            endif
            alpha_ij = dsqrt(9.d0*B(i,j)*omega(i,j)/Ec(i,j))

            if(d_input.eq.0.0d0) then
               d_ij   = .25d0 * Eu_d(i) + .75d0 * Eu_d(j)
            elseif(d_input.lt.0.0d0) then
               d_ij   =  0.0d0  
            else
               d_ij   =  d_input
            endif

            Z_ij      = 12.0d0
              scr_2nn = screen_on_2nn(lattice,Cmax(i,j,i),Cmin(i,j,i))
            Z2_i      = 6.0d0 * scr_2nn
c  Here, Cmin(j,i,j) should be equal to Cmin(j,j,j) .. 2002. 6. 14, B.-J. Lee
              scr_2nn = screen_on_2nn(lattice,Cmax(j,j,j),Cmin(j,j,j))
            Z2_j      = 6.0d0 * scr_2nn

            if(C_iij.eq.0.0) then
               tmp_ij =.5d0*dsqrt(Cmin(i,i,i))+.5d0*dsqrt(Cmin(j,j,j))
               Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               Cmin(i,i,j) = C_iij
            endif
            Cmin(j,i,i) = Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij =.5d0*dsqrt(Cmin(i,i,i))+.5d0*dsqrt(Cmin(j,j,j))
               Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               Cmin(i,j,j) = C_ijj
            endif
            Cmin(j,j,i) = Cmin(i,j,j)
c
         elseif(lattice.eq.'L12A3B') then
c
c        for fcc_L12 (A3B type)
c
            a2N     = dsqrt(2.d0)
            Ec(i,j) = .75d0 * Ec(i,i) + .25d0 * Ec(j,j) - delta_Ec
            if(Re_input.eq.0.0d0) then
               omega(i,j) = .75d0 * omega(i,i) + .25d0 * omega(j,j)
                   tmp_ij = 4.d0 * omega(i,j)
               Re_ij      = tmp_ij**(1.d0/3.d0) / a2N
            else
               Re_ij      =  Re_input
                   tmp_ij =  a2N * Re_ij
               omega(i,j) = .25d0 * tmp_ij * tmp_ij * tmp_ij 
            endif
            if(B_ij.eq.0.0d0) then
               B(i,j) = .75d0 * B(i,i) + .25d0 * B(j,j)
            else
               B(i,j) =  B_ij
            endif
            alpha_ij = dsqrt(9.d0*B(i,j)*omega(i,j)/Ec(i,j))

            if(d_input.eq.0.0d0) then
               d_ij   = .75d0 * Eu_d(i) + .25d0 * Eu_d(j)
            elseif(d_input.lt.0.0d0) then
               d_ij   =  0.0d0  
            else
               d_ij   =  d_input
            endif

            Z_ij      = 12.0d0
c  Here, Cmin(i,j,i) should be equal to Cmin(i,i,i) .. 2002. 6. 14, B.-J. Lee
              scr_2nn = screen_on_2nn(lattice,Cmax(i,i,i),Cmin(i,i,i))
            Z2_i      = 6.0d0 * scr_2nn
              scr_2nn = screen_on_2nn(lattice,Cmax(j,i,j),Cmin(j,i,j))
            Z2_j      = 6.0d0 * scr_2nn

            if(C_iij.eq.0.0) then
               tmp_ij =.5d0*dsqrt(Cmin(i,i,i))+.5d0*dsqrt(Cmin(j,j,j))
               Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               Cmin(i,i,j) = C_iij
            endif
            Cmin(j,i,i) = Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij =.5d0*dsqrt(Cmin(i,i,i))+.5d0*dsqrt(Cmin(j,j,j))
               Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               Cmin(i,j,j) = C_ijj
            endif
            Cmin(j,j,i) = Cmin(i,j,j)
c
         elseif(lattice.eq.'ZnS_B3') then
c
c          for ZnS_B3
c
            a2N     = 4.d0 / dsqrt(6.d0)
            Ec(i,j)    = .5d0 * (Ec(i,i) + Ec(j,j)) - delta_Ec
            if(Re_input.eq.0.0d0) then
               omega(i,j) = .5d0 * (omega(i,i) + omega(j,j))
                   tmp_ij = 8.d0 * omega(i,j)
               Re_ij      = tmp_ij**(1.d0/3.d0) * dsqrt(3.d0) / 4.d0
            else
               Re_ij      =  Re_input
                   tmp_ij =  4.d0 * Re_ij / dsqrt(3.d0)
               omega(i,j) = .125d0 * tmp_ij * tmp_ij * tmp_ij 
            endif
            if(B_ij.eq.0.0d0) then
               B(i,j) = .5d0 * B(i,i) + .5d0 * B(j,j)
            else
               B(i,j) =  B_ij
            endif
            alpha_ij  = dsqrt(9.d0*B(i,j)*omega(i,j)/Ec(i,j))

            if(d_input.eq.0.0d0) then
               d_ij   = .5d0 * Eu_d(i) + .5d0 * Eu_d(j)
            elseif(d_input.lt.0.0d0) then
               d_ij   =  0.0d0  
            else
               d_ij   =  d_input
            endif

            Z_ij      = 4.0d0
              scr_2nn = screen_on_2nn(lattice,Cmax(i,j,i),Cmin(i,j,i))
            Z2_i      = 12.0d0 * scr_2nn
              scr_2nn = screen_on_2nn(lattice,Cmax(j,i,j),Cmin(j,i,j))
            Z2_j      = 12.0d0 * scr_2nn

            if(C_iij.eq.0.0) then
               tmp_ij = 0.5d0*(dsqrt(Cmin(i,i,i))+dsqrt(Cmin(j,j,j)))
               Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               Cmin(i,i,j) = C_iij
            endif
            Cmin(j,i,i) = Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij = 0.5d0*(dsqrt(Cmin(i,i,i))+dsqrt(Cmin(j,j,j)))
               Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               Cmin(i,j,j) = C_ijj
            endif
            Cmin(j,j,i) = Cmin(i,j,j)
c
         elseif(lattice.eq.'DIMER') then
c
c          for Diatomic Gases
c
            a2N     = dsqrt(2.d0)
            Ec(i,j) = .5d0 * (Ec(i,i) + Ec(j,j)) - delta_Ec
            if(Re_input.eq.0.0d0) then
               omega(i,j) = .5d0 * (omega(i,i) + omega(j,j))
               Re_ij      = omega(i,j)**(1.d0/3.d0)
            else
               Re_ij      =  Re_input
               omega(i,j) = Re_ij * Re_ij * Re_ij
            endif
            if(B_ij.eq.0.0d0) then
               B(i,j) = .5d0 * B(i,i) + .5d0 * B(j,j)
            else
               B(i,j) =  B_ij
            endif
            alpha_ij  = dsqrt(9.d0*B(i,j)*omega(i,j)/Ec(i,j))

            if(d_input.eq.0.0d0) then
               d_ij   = .5d0 * Eu_d(i) + .5d0 * Eu_d(j)
            elseif(d_input.lt.0.0d0) then
               d_ij   =  0.0d0  
            else
               d_ij   =  d_input
            endif

            Z_ij      = 1.0d0
              scr_2nn = 0.d0
            Z2_i      = scr_2nn
            Z2_j      = scr_2nn

            if(C_iij.eq.0.0) then
               tmp_ij = 0.5d0*(dsqrt(Cmin(i,i,i))+dsqrt(Cmin(j,j,j)))
               Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               Cmin(i,i,j) = C_iij
            endif
            Cmin(j,i,i) = Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij = 0.5d0*(dsqrt(Cmin(i,i,i))+dsqrt(Cmin(j,j,j)))
               Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               Cmin(i,j,j) = C_ijj
            endif
            Cmin(j,j,i) = Cmin(i,j,j)
c
         else
            goto 4445
         endif
c
         Re_Alloy(i,j) = Re_ij
         alpha(i,j) = alpha_ij
         write(*,'(1x,a5,a7,f7.4,a10,f7.4)') 
     &   alloy_1_name,':  Re =',Re_ij,',  alpha =',alpha_ij
c
c -----------------------------------------------------------------------------
c
c  Make a table for pair potential and its derivative between different element
c
         do k=1,LineTable+1
            Dij = RMin + (k-1) * DeltaR
c
c        Calculation of Eu(R) and dEdR for reference structure at r = Dij
c
            Rfactor = Dij / Re_ij - 1.d0
            astar   = alpha_ij * Rfactor
            astar3  = astar * astar * astar
            Eu = - Ec(i,j) * (1.d0 + astar + d_ij*astar3) * dexp(-astar)
            dEdR = Ec(i,j) * astar * alpha_ij / Re_ij * dexp(-astar)
     &           * (1.d0 + d_ij * astar * astar - 3.d0 * d_ij * astar)
c
c        Calculation of Phi and DPhi at r = a2N * Dij for each element
c
            if(lattice(1:4).ne.'MC_B' .and. lattice(1:3).ne.'DIM') then
            rr = a2N * Dij
            ZR2 = Z2(i) / Z(i)
            alphaN = alpha(i,i) 
            scr_2nn = 
     &        screen_on_2nn(lattice_element(i),Cmax(i,i,i),Cmin(i,i,i))
            call Define_Pair_Potential
     &          (i,lattice_element(i),alphaN,ZR2,a2(i),rr,Phi_i,DPhi_i)
            ZR2 = Z2(j) / Z(j)
            alphaN = alpha(j,j) 
            scr_2nn = 
     &        screen_on_2nn(lattice_element(j),Cmax(j,j,j),Cmin(j,j,j))
            call Define_Pair_Potential
     &          (j,lattice_element(j),alphaN,ZR2,a2(j),rr,Phi_j,DPhi_j)
            endif
c
            R12factor = dexp( - beta(j,4) * (Dij / Re(j) - 1.d0) )
            R21factor = dexp( - beta(i,4) * (a2N*Dij/Re(i) - 1.d0) )
            R11factor = dexp( - beta(i,4) * (Dij / Re(i) - 1.d0) )
            R22factor = dexp( - beta(j,4) * (a2N*Dij/Re(j) - 1.d0) )
c
         if(lattice.eq.'FCC_B1' .or. lattice.eq.'NaCl' .or.
     &      lattice.eq.'BCC_B2') then
c
c           calculation of Rho_Bar, F, dFdR for reference structure at r = Dij
c
            Rho_Bar_i = Z_ij * RhoZ(j) * R12factor
     &                + Z2_i * RhoZ(i) * R21factor
            Rho_Bar_j = Z_ij * RhoZ(i) * R11factor
     &                + Z2_j * RhoZ(j) * R22factor
c
            if(rho_bar_i.eq.0.d0) rho_bar_i = 1.0d-20
            if(rho_bar_j.eq.0.d0) rho_bar_j = 1.0d-20
            rr0_i = Rho_Bar_i / Rho0_Bar(i)
            rr0_j = Rho_Bar_j / Rho0_Bar(j)
            F_i = A(i) * Ec(i,i) * rr0_i * dlog(rr0_i)
            F_j = A(j) * Ec(j,j) * rr0_j * dlog(rr0_j)
c
c           dFdR = - beta(4)/Re * A * Ec * rr0 * ( dlog(rr0) + 1.d0 )
c
            dRho_i = beta(j,4) / Re(j) * Z_ij * RhoZ(j) * R12factor
     &       + a2N * beta(i,4) / Re(i) * Z2_i * RhoZ(i) * R21factor
            dRho_j = beta(i,4) / Re(i) * Z_ij * RhoZ(i) * R11factor
     &       + a2N * beta(j,4) / Re(j) * Z2_j * RhoZ(j) * R22factor
            dFdR_i = - A(i) * Ec(i,i) / Rho0_Bar(i) 
     &             * ( dlog(rr0_i) + 1.d0 ) * dRho_i
            dFdR_j = - A(j) * Ec(j,j) / Rho0_Bar(j)
     &             * ( dlog(rr0_j) + 1.d0 ) * dRho_j
c
c           Pair potential and its derivative between i and j atoms at r = Dij
c
            PhiTab(k,i,j) = 1.d0 / Z_ij *
     &        (2.d0*Eu - F_i - F_j - .5d0*(Z2_i*Phi_i + Z2_j*Phi_j))
            PhiTab(k,j,i) = PhiTab(k,i,j)
            DPhiTab(k,i,j) = (2.d0*dEdR -dFdR_i -dFdR_j) / Z_ij / Dij
     &          - .5d0 * a2N * a2N * (Z2_i*DPhi_i+Z2_j*DPhi_j) / Z_ij
            DPhiTab(k,j,i) = DPhiTab(k,i,j)
c
         elseif(lattice.eq.'L12AB3') then
c
c           calculation of Rho_Bar and F for reference structure at r = Dij
c
            Rho_Bar_i = Z_ij * RhoZ(j) * R12factor
     &                + Z2_i * RhoZ(i) * R21factor

            Rho0 = 4.d0 * RhoZ(i) * R11factor
     &           + 8.d0 * RhoZ(j) * R12factor
     &           + Z2_j * RhoZ(j) * R22factor
            Rho2 = RhoZ(j) * dexp( - beta(j,2) * (Dij/Re(j) - 1.d0) )
     &           - RhoZ(i) * dexp( - beta(i,2) * (Dij/Re(i) - 1.d0) )
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.d0) RhoFactor = 1.d-20
            Rho2P = 8.d0 * Rho2 * Rho2 / 3.d0
            GAMMA = t(j,2) * Rho2P / RhoFactor
            Rho_Bar_j = Rho0 * 2.d0 / ( 1.d0 + dexp(-GAMMA) )
c
            if(rho_bar_i.eq.0.d0) rho_bar_i = 1.0d-20
            if(rho_bar_j.eq.0.d0) rho_bar_j = 1.0d-20
            rr0_i = Rho_Bar_i / Rho0_Bar(i)
            rr0_j = Rho_Bar_j / Rho0_Bar(j)
            F_i = A(i) * Ec(i,i) * rr0_i * dlog(rr0_i)
            F_j = A(j) * Ec(j,j) * rr0_j * dlog(rr0_j)
c
c           Pair potential and its derivative between i and j atoms at r = Dij
c
            PhiTab(k,i,j) = 1.d0 / Z_ij *
     &        ( 4.d0 * Eu - F_i - 3.d0*F_j - Z_ij * PhiTab(k,j,j) 
     &          - .5d0 * (Z2_i * Phi_i + 3.d0 * Z2_j * Phi_j) )
            PhiTab(k,j,i) = PhiTab(k,i,j)
            if(k.gt.2) then
               DPhiTab(k-1,i,j) = .5d0*(PhiTab(k,i,j) - PhiTab(k-2,i,j))
     &                          / DeltaR / (Dij-DeltaR)
               DPhiTab(k-1,j,i) = DPhiTab(k-1,i,j)
            endif
c
         elseif(lattice.eq.'L12A3B') then
c
c           calculation of Rho_Bar and F for reference structure at r = Dij
c
            Rho_Bar_j = Z_ij * RhoZ(i) * R11factor
     &                + Z2_j * RhoZ(j) * R22factor

            Rho0 = 4.d0 * RhoZ(j) * R12factor
     &           + 8.d0 * RhoZ(i) * R11factor
     &           + Z2_i * RhoZ(i) * R21factor
            Rho2 = RhoZ(i) * dexp( - beta(i,2) * (Dij/Re(i) - 1.d0) )
     &           - RhoZ(j) * dexp( - beta(j,2) * (Dij/Re(j) - 1.d0) )
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.d0) RhoFactor = 1.d-20
            Rho2P = 8.d0 * Rho2 * Rho2 / 3.d0
            GAMMA = t(i,2) * Rho2P / RhoFactor
            Rho_Bar_i = Rho0 * 2.d0 / ( 1.d0 + dexp(-GAMMA) )
c
            if(rho_bar_i.eq.0.d0) rho_bar_i = 1.0d-20
            if(rho_bar_j.eq.0.d0) rho_bar_j = 1.0d-20
            rr0_i = Rho_Bar_i / Rho0_Bar(i)
            rr0_j = Rho_Bar_j / Rho0_Bar(j)
            F_i = A(i) * Ec(i,i) * rr0_i * dlog(rr0_i)
            F_j = A(j) * Ec(j,j) * rr0_j * dlog(rr0_j)
c
c           Pair potential and its derivative between i and j atoms at r = Dij
c
            PhiTab(k,i,j) = 1.d0 / Z_ij *
     &        ( 4.d0 * Eu - 3.d0*F_i - F_j - Z_ij * PhiTab(k,i,i) 
     &          - .5d0 * (Z2_j * Phi_j + 3.d0 * Z2_i * Phi_i) )
            PhiTab(k,j,i) = PhiTab(k,i,j)
            if(k.gt.2) then
               DPhiTab(k-1,i,j) = .5d0*(PhiTab(k,i,j) - PhiTab(k-2,i,j))
     &                          / DeltaR / (Dij-DeltaR)
               DPhiTab(k-1,j,i) = DPhiTab(k-1,i,j)
            endif
c
         elseif(lattice.eq.'MC_BCC') then
c
c           calculation of Rho_Bar and F for reference structure at r = Dij
c
            Rho_Bar_i = Z_ij * RhoZ(j) * R12factor

            Rho0 = 2.d0 * RhoZ(i) * R11factor
     &           + 4.d0 * RhoZ(j) * R12factor
            Rho2 = RhoZ(j) * dexp( - beta(j,2) * (Dij/Re(j) - 1.d0) )
     &           - RhoZ(i) * dexp( - beta(i,2) * (Dij/Re(i) - 1.d0) )
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.d0) RhoFactor = 1.d-20
            Rho2P = 8.d0 * Rho2 * Rho2 / 3.d0
            GAMMA = t(j,2) * Rho2P / RhoFactor
            Rho_Bar_j = Rho0 * 2.d0 / ( 1.d0 + dexp(-GAMMA) )
c
            if(rho_bar_i.eq.0.d0) rho_bar_i = 1.0d-20
            if(rho_bar_j.eq.0.d0) rho_bar_j = 1.0d-20
            rr0_i = Rho_Bar_i / Rho0_Bar(i)
            rr0_j = Rho_Bar_j / Rho0_Bar(j)
            F_i = A(i) * Ec(i,i) * rr0_i * dlog(rr0_i)
            F_j = A(j) * Ec(j,j) * rr0_j * dlog(rr0_j)
c
c           Pair potential and its derivative between i and j atoms at r = Dij
c
            PhiTab(k,i,j) = 1.d0 / Z_ij *
     &        ( 4.d0 * Eu - F_i - 3.d0 * F_j - Z_ij * PhiTab(k,j,j) )
            PhiTab(k,j,i) = PhiTab(k,i,j)
            if(k.gt.2) then
               DPhiTab(k-1,i,j) = .5d0*(PhiTab(k,i,j) - PhiTab(k-2,i,j))
     &                          / DeltaR / (Dij-DeltaR)
               DPhiTab(k-1,j,i) = DPhiTab(k-1,i,j)
            endif
c
         elseif(lattice.eq.'ZnS_B3') then
c
c           calculation of Rho_Bar and F for reference structure at r = Dij
c
            Rho0 = Z_ij * RhoZ(j) * R12factor
     &           + Z2_i * RhoZ(i) * R21factor
            Rho3 = RhoZ(j) * dexp( - beta(j,3) * (Dij/Re(j) - 1.d0) )
            Rho3P = 32.d0 * Rho3 * Rho3 / 9.d0
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.d0) RhoFactor = 1.d-20
            GAMMA = t(i,3) * Rho3P / RhoFactor
            Rho_Bar_i = Rho0 * 2.d0 / ( 1.d0 + dexp(-GAMMA) )

            Rho0 = Z_ij * RhoZ(i) * R11factor
     &           + Z2_j * RhoZ(j) * R22factor
            Rho3 = RhoZ(i) * dexp( - beta(i,3) * (Dij/Re(i) - 1.d0) )
            Rho3P = 32.d0 * Rho3 * Rho3 / 9.d0
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.d0) RhoFactor = 1.d-20
            GAMMA = t(j,3) * Rho3P / RhoFactor
            Rho_Bar_j = Rho0 * 2.d0 / ( 1.d0 + dexp(-GAMMA) )
c
            if(rho_bar_i.eq.0.d0) rho_bar_i = 1.0d-20
            if(rho_bar_j.eq.0.d0) rho_bar_j = 1.0d-20
            rr0_i = Rho_Bar_i / Rho0_Bar(i)
            rr0_j = Rho_Bar_j / Rho0_Bar(j)
            F_i = A(i) * Ec(i,i) * rr0_i * dlog(rr0_i)
            F_j = A(j) * Ec(j,j) * rr0_j * dlog(rr0_j)
c
c           Pair potential and its derivative between i and j atoms at r = Dij
c
            PhiTab(k,i,j) = 1.d0 / Z_ij *
     &        (2.d0*Eu - F_i - F_j - .5d0*(Z2_i*Phi_i + Z2_j*Phi_j))
            PhiTab(k,j,i) = PhiTab(k,i,j)
            if(k.gt.2) then
               DPhiTab(k-1,i,j) = .5d0*(PhiTab(k,i,j) - PhiTab(k-2,i,j))
     &                          / DeltaR / (Dij-DeltaR)
               DPhiTab(k-1,j,i) = DPhiTab(k-1,i,j)
            endif
c
         elseif(lattice.eq.'DIMER') then
c
c           calculation of Rho_Bar and F for reference structure at r = Dij
c
            Rho0 = Z_ij * RhoZ(j) * R12factor
            Rho1 = RhoZ(j) * dexp( - beta(j,1) * (Dij/Re(j) - 1.d0) )
            Rho2 = RhoZ(j) * dexp( - beta(j,2) * (Dij/Re(j) - 1.d0) )
            Rho3 = RhoZ(j) * dexp( - beta(j,3) * (Dij/Re(j) - 1.d0) )
            Rho1P = Rho1 * Rho1
            Rho2P = 2.d0 * Rho2 * Rho2 / 3.d0
            Rho3P = 2.d0 * Rho3 * Rho3 / 5.d0
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.d0) RhoFactor = 1.d-20
            GAMMA = (t(i,1)*Rho1P+t(i,2)*Rho2P+t(i,3)*Rho3P) / RhoFactor
            Rho_Bar_i = Rho0 * 2.d0 / ( 1.d0 + dexp(-GAMMA) )
c           Rho_Bar_i = Rho0 * dsqrt ( 1.d0 + GAMMA )

            Rho0 = Z_ij * RhoZ(i) * R11factor
            Rho1 = RhoZ(i) * dexp( - beta(i,1) * (Dij/Re(i) - 1.d0) )
            Rho2 = RhoZ(i) * dexp( - beta(i,2) * (Dij/Re(i) - 1.d0) )
            Rho3 = RhoZ(i) * dexp( - beta(i,3) * (Dij/Re(i) - 1.d0) )
            Rho1P = Rho1 * Rho1
            Rho2P = 2.d0 * Rho2 * Rho2 / 3.d0
            Rho3P = 2.d0 * Rho3 * Rho3 / 5.d0
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.d0) RhoFactor = 1.d-20
            GAMMA = (t(j,1)*Rho1P+t(j,2)*Rho2P+t(j,3)*Rho3P) / RhoFactor
            Rho_Bar_j = Rho0 * 2.d0 / ( 1.d0 + dexp(-GAMMA) )
c           Rho_Bar_j = Rho0 * dsqrt ( 1.d0 + GAMMA )
c
            if(rho_bar_i.eq.0.d0) rho_bar_i = 1.0d-20
            if(rho_bar_j.eq.0.d0) rho_bar_j = 1.0d-20
            rr0_i = Rho_Bar_i / Rho0_Bar(i)
            rr0_j = Rho_Bar_j / Rho0_Bar(j)
            F_i = A(i) * Ec(i,i) * rr0_i * dlog(rr0_i)
            F_j = A(j) * Ec(j,j) * rr0_j * dlog(rr0_j)
c
c           Pair potential and its derivative between i and j atoms at r = Dij
c
            PhiTab(k,i,j) = (2.d0*Eu - F_i - F_j) / Z_ij
            PhiTab(k,j,i) = PhiTab(k,i,j)
            if(k.gt.2) then
               DPhiTab(k-1,i,j) = .5d0*(PhiTab(k,i,j) - PhiTab(k-2,i,j))
     &                          / DeltaR / (Dij-DeltaR)
               DPhiTab(k-1,j,i) = DPhiTab(k-1,i,j)
            endif
c
c        elseif(lattice.eq.'HCP_A3' .or. lattice.eq.'PtCo') then
c     Description for PtCo type ordered HCP structure should be given
c     2000. 5. 4.  B.-J. Lee
c
         endif
         enddo
      enddo
      enddo
      endif
c     
      if(Ncomp.gt.2) then
         ternary_name(1:2) = Icomp(1)         
         ternary_name(3:3) = '-'
         ternary_name(4:5) = Icomp(2)         
         ternary_name(6:6) = '-'
         ternary_name(7:8) = Icomp(3)         
         
         do i = 1, Ncomp-2
         do j = i+1, Ncomp-1
         do k = j+1, Ncomp
c
            read(2,'(a)') card
            read(2,'(a)') card2
    3       continue
            read(2,'(a)',end=4447) card
            read(card,'(a8)') alloy3
            if(alloy3.ne.ternary_name) goto 3

            read(card,2011) Cx_ikj, Cx_ijk, Cx_jik, C_ikj, C_ijk, C_jik
 2011       format(40x,6f5.2)
 
            if(C_ikj.eq.0.0) then
               tmp_ij = 0.5d0*(dsqrt(Cmin(i,i,i))+dsqrt(Cmin(j,j,j)))
               Cmin(i,k,j) = tmp_ij * tmp_ij
            else
               Cmin(i,k,j) = C_ikj
            endif
            Cmin(j,k,i) = Cmin(i,k,j)

            if(C_ijk.eq.0.0) then
               tmp_ij = 0.5d0*(dsqrt(Cmin(i,i,i))+dsqrt(Cmin(k,k,k)))
               Cmin(i,j,k) = tmp_ij * tmp_ij
            else
               Cmin(i,j,k) = C_ijk
            endif
            Cmin(k,j,i) = Cmin(i,j,k)

            if(C_jik.eq.0.0) then
               tmp_ij = 0.5d0*(dsqrt(Cmin(j,j,j))+dsqrt(Cmin(k,k,k)))
               Cmin(j,i,k) = tmp_ij * tmp_ij
            else
               Cmin(j,i,k) = C_jik
            endif
            Cmin(k,i,j) = Cmin(j,i,k)

            if(Cx_ikj.ne.0.0) then
               Cmax(i,k,j) = Cx_ikj
               Cmax(j,k,i) = Cx_ikj
            endif
            if(Cx_ijk.ne.0.0) then
               Cmax(i,j,k) = Cx_ijk
               Cmax(k,j,i) = Cx_ijk
            endif
            if(Cx_jik.ne.0.0) then
               Cmax(j,i,k) = Cx_jik
               Cmax(k,i,j) = Cx_jik
            endif
c
         enddo
         enddo
         enddo
      endif
c
      Ierr = 0
      close(unit=2)
      return
c
 4444 continue
      Ierr = 1
      write(*,'(/,a,a,/)') ' **** ERROR : Non-Registered Element ',
     &      Icomp(n)
      return
 4445 continue
      Ierr = 1
      write(*,'(/,a,a,/)') ' **** ERROR : Non-Registered Structure ',
     &      lattice
      return
 4446 continue
      Ierr = 1
      write(*,'(/,a,a,/)') ' **** ERROR : Non-Registered Alloy ',
     &      alloy_1_name
      return
 4447 continue
      Ierr = 1
      write(*,'(/,a,a,/)') ' **** ERROR : Non-Registered Alloy ',
     &      ternary_name
      return
      end


      function screen_on_2nn(lattice,Cmax,Cmin)
c
c  Gives the value of screening function on 2nd NN atoms as a fuction of Cmin 
c  for the given structure
c
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      character lattice*6
c
      if(lattice(1:6).eq.'FCC_A1' .or. lattice(1:3).eq.'L12'
     &  .or. lattice(1:3).eq.'HCP') then
         if(Cmax.eq.2.80d0) then
            if(Cmin.eq.0.04d0) then
               screen_on_2nn =  2.02613642d-01
            elseif(Cmin.eq.0.09d0) then
               screen_on_2nn =  1.76993632d-01
            elseif(Cmin.eq.0.16d0) then
               screen_on_2nn =  1.42574615d-01
            elseif(Cmin.eq.0.25d0) then
               screen_on_2nn =  1.01972180161d-01
            elseif(Cmin.eq.0.36d0) then
               screen_on_2nn =  6.0225031562d-02
            elseif(Cmin.eq.0.49d0) then
               screen_on_2nn =  2.5236627552d-02
            elseif(Cmin.eq.0.52d0) then
               screen_on_2nn =  1.956036545d-02
            elseif(Cmin.eq.0.64d0) then
               screen_on_2nn =  5.163452059d-03
            elseif(Cmin.eq.0.74d0) then
               screen_on_2nn =  9.15416424d-04
            elseif(Cmin.eq.0.78d0) then
               screen_on_2nn =  3.47474570d-04
            elseif(Cmin.eq.0.80d0) then
               screen_on_2nn =  1.95639888d-04
            elseif(Cmin.eq.0.81d0) then
               screen_on_2nn =  1.42742161d-04
            else
               screen_on_2nn =  0.d0
            endif
         elseif(Cmax.eq.1.44d0) then
            if(Cmin.eq.0.04d0) then
               screen_on_2nn =  9.24561664010d-01
            elseif(Cmin.eq.0.09d0) then
               screen_on_2nn =  9.13211457257d-01
            elseif(Cmin.eq.0.16d0) then
               screen_on_2nn =  8.93607056678d-01
            elseif(Cmin.eq.0.25d0) then
               screen_on_2nn =  8.59899566947d-01
            elseif(Cmin.eq.0.36d0) then
               screen_on_2nn =  7.99723000590d-01
            elseif(Cmin.eq.0.49d0) then
               screen_on_2nn =  6.86002760884d-01
            elseif(Cmin.eq.0.52d0) then
               screen_on_2nn =  6.50575052285d-01
            elseif(Cmin.eq.0.64d0) then
               screen_on_2nn =  4.64061512672d-01
            elseif(Cmin.eq.0.74d0) then
               screen_on_2nn =  2.57220007402d-01
            elseif(Cmin.eq.0.78d0) then
               screen_on_2nn =  1.71959709240d-01
            elseif(Cmin.eq.0.80d0) then
               screen_on_2nn =  1.32300505821d-01
            elseif(Cmin.eq.0.81d0) then
               screen_on_2nn =  1.13752424573d-01
            elseif(Cmin.eq.0.90d0) then
               screen_on_2nn =  9.56257003220d-03
            elseif(Cmin.eq.0.95d0) then
               screen_on_2nn =  2.24318162818d-04
            else
               screen_on_2nn =  0.d0
            endif
         endif
      elseif(lattice(1:3).eq.'BCC') then
         if(Cmax.eq.2.80d0) then
         if(Cmin.eq.0.04d0) then
            screen_on_2nn =  0.9449061210945d0
         elseif(Cmin.eq.0.09d0) then
            screen_on_2nn =  0.940836817464d0
         elseif(Cmin.eq.0.16d0) then
            screen_on_2nn =  0.934499543653d0
         elseif(Cmin.eq.0.25d0) then
            screen_on_2nn =  0.9250795156595d0
         elseif(Cmin.eq.0.36d0) then
            screen_on_2nn =  0.91120732155d0
         elseif(Cmin.eq.0.49d0) then
            screen_on_2nn =  0.890549742627d0
         elseif(Cmin.eq.0.64d0) then
            screen_on_2nn =  0.85901541d0
         elseif(Cmin.eq.0.74d0) then
            screen_on_2nn =  0.831883d0
         elseif(Cmin.eq.0.78d0) then
            screen_on_2nn =  0.819328578d0
         elseif(Cmin.eq.0.8d0) then
            screen_on_2nn =  0.81264d0
         elseif(Cmin.eq.0.81d0) then
            screen_on_2nn =  0.809187293d0
         elseif(Cmin.eq.1.21d0) then
            screen_on_2nn =  0.588685d0
         elseif(Cmin.eq.1.44d0) then
            screen_on_2nn =  0.360517d0
         elseif(Cmin.eq.1.46d0) then
            screen_on_2nn =  0.33725192772d0
         elseif(Cmin.eq.1.48d0) then
            screen_on_2nn =  0.31366443673d0
         elseif(Cmin.eq.1.50d0) then
            screen_on_2nn =  0.289850789d0
         elseif(Cmin.eq.1.52d0) then
            screen_on_2nn =  0.2659238719d0
         elseif(Cmin.eq.1.69d0) then
            screen_on_2nn =  0.080808587d0
         else
            screen_on_2nn =  0.d0
         endif
         elseif(Cmax.eq.1.44d0) then
            screen_on_2nn =  1.d0
         endif
      elseif(lattice(1:6).eq.'FCC_B1' .or. lattice(1:4).eq.'NaCl'
     &  .or. lattice(1:2).eq.'SC' .or. lattice(1:6).eq.'MC_BCC') then
         if(Cmax.eq.2.80d0) then
            if(Cmin.eq.0.04d0) then
               screen_on_2nn =  0.450126251261d0
            elseif(Cmin.eq.0.09d0) then
               screen_on_2nn =  0.420706111242d0
            elseif(Cmin.eq.0.16d0) then
               screen_on_2nn =  0.377590538937d0
            elseif(Cmin.eq.0.25d0) then
               screen_on_2nn =  0.319330831837d0
            elseif(Cmin.eq.0.36d0) then
               screen_on_2nn =  0.245407888141d0
            elseif(Cmin.eq.0.49d0) then
               screen_on_2nn =  0.158860402718d0
            elseif(Cmin.eq.0.64d0) then
               screen_on_2nn =  7.1857164283d-02
            elseif(Cmin.eq.0.74d0) then
               screen_on_2nn =  3.025584941d-02
            elseif(Cmin.eq.0.78d0) then
               screen_on_2nn =  1.864066979d-02
            elseif(Cmin.eq.0.8d0) then
               screen_on_2nn =  1.3987132961d-02
            elseif(Cmin.eq.0.81d0) then
               screen_on_2nn =  1.1947475102d-02
            else
               screen_on_2nn =  0.d0
            endif
         elseif(Cmax.eq.2.00d0) then
            if(Cmin.eq.0.04d0) then
               screen_on_2nn =  0.755283957216d0
            elseif(Cmin.eq.0.09d0) then
               screen_on_2nn =  0.731653578617d0
            elseif(Cmin.eq.0.16d0) then
               screen_on_2nn =  0.694098984589d0
            elseif(Cmin.eq.0.25d0) then
               screen_on_2nn =  0.637001633253d0
            elseif(Cmin.eq.0.36d0) then
               screen_on_2nn =  0.551507187967d0
            elseif(Cmin.eq.0.49d0) then
               screen_on_2nn =  0.425493343617d0
            elseif(Cmin.eq.0.64d0) then
               screen_on_2nn =  0.250825155120d0
            elseif(Cmin.eq.0.74d0) then
               screen_on_2nn =  0.132430188269d0
            elseif(Cmin.eq.0.78d0) then
               screen_on_2nn =  9.057870368d-02
            elseif(Cmin.eq.0.8d0) then
               screen_on_2nn =  7.1857164283d-02
            elseif(Cmin.eq.0.81d0) then
               screen_on_2nn =  6.3168287352d-02
            else
               screen_on_2nn =  0.d0
            endif
         elseif(Cmax.eq.1.44d0) then
            if(Cmin.eq.0.04d0) then
               screen_on_2nn =  0.961541036051d0
            elseif(Cmin.eq.0.09d0) then
               screen_on_2nn =  0.955620979917d0
            elseif(Cmin.eq.0.16d0) then
               screen_on_2nn =  0.945307916331d0
            elseif(Cmin.eq.0.25d0) then
               screen_on_2nn =  0.927307698095d0
            elseif(Cmin.eq.0.36d0) then
               screen_on_2nn =  0.894272330216d0
            elseif(Cmin.eq.0.49d0) then
               screen_on_2nn =  0.828252836327d0
            elseif(Cmin.eq.0.64d0) then
               screen_on_2nn =  0.681220604996d0
            elseif(Cmin.eq.0.74d0) then
               screen_on_2nn =  0.507168618313d0
            elseif(Cmin.eq.0.78d0) then
               screen_on_2nn =  0.414680249397d0
            elseif(Cmin.eq.0.8d0) then
               screen_on_2nn =  0.363731364912d0
            elseif(Cmin.eq.0.81d0) then
               screen_on_2nn =  0.337272033488d0
            elseif(Cmin.eq.0.85d0) then
               screen_on_2nn =  0.227571090763d0
            elseif(Cmin.eq.0.90d0) then
               screen_on_2nn =  9.77883941590d-02
            elseif(Cmin.eq.0.95d0) then
               screen_on_2nn =  1.49772548492d-02
            else
               screen_on_2nn =  0.d0
            endif
         elseif(Cmax.eq.1.00d0) then
            if(Cmin.le.1.00d0) then
               screen_on_2nn =  1.d0
            else
               screen_on_2nn =  0.d0
            endif
         endif
      elseif(lattice(1:3).eq.'DIA' .or. lattice(1:3).eq.'ZnS') then
         if(Cmin.eq.0.04d0) then
            screen_on_2nn =  0.268061866522d0
         elseif(Cmin.eq.0.09d0) then
            screen_on_2nn =  0.231514106451d0
         elseif(Cmin.eq.0.16d0) then
            screen_on_2nn =  0.179693727623d0
         elseif(Cmin.eq.0.25d0) then
            screen_on_2nn =  0.114354647869d0
         elseif(Cmin.eq.0.36d0) then
            screen_on_2nn =  4.4310323188d-02
         elseif(Cmin.eq.0.49d0) then
            screen_on_2nn =  2.95974168d-04
         else
            screen_on_2nn =  0.d0
         endif
      else
         screen_on_2nn =  0.d0
      endif
c
      return
      end


      subroutine Define_Pair_Potential
     &          (N,lattice,alphaN,ZR2,a2N,Dij,Phi,DPhi)
c
c  Compute phi(R) and dphi(R) for element according to 2nd NN approximation
c
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      character lattice*6
c
      call Compute_Phi(N,lattice,alphaN,1.d0,Dij,P_Phi,P_DPhi)
           Phi  = P_Phi
           DPhi = P_DPhi
c
      if(ZR2.eq.0.d0) return
c
      call Compute_Phi(N,lattice,alphaN,a2N,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2 * P_Phi
           DPhi = DPhi - ZR2 * P_DPhi
      call Compute_Phi(N,lattice,alphaN,a2N*a2N,Dij,P_Phi,P_DPhi)
           Phi  = Phi  + ZR2*ZR2 * P_Phi
           DPhi = DPhi + ZR2*ZR2 * P_DPhi
      call Compute_Phi(N,lattice,alphaN,a2N**3,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2*ZR2*ZR2 * P_Phi
           DPhi = DPhi - ZR2*ZR2*ZR2 * P_DPhi
      call Compute_Phi(N,lattice,alphaN,a2N**4,Dij,P_Phi,P_DPhi)
           Phi  = Phi  + ZR2*ZR2*ZR2*ZR2 * P_Phi
           DPhi = DPhi + ZR2*ZR2*ZR2*ZR2 * P_DPhi
      call Compute_Phi(N,lattice,alphaN,a2N**5,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2*ZR2*ZR2*ZR2*ZR2 * P_Phi
           DPhi = DPhi - ZR2*ZR2*ZR2*ZR2*ZR2 * P_DPhi
      call Compute_Phi(N,lattice,alphaN,a2N**6,Dij,P_Phi,P_DPhi)
           Phi  = Phi  + ZR2*ZR2*ZR2*ZR2*ZR2*ZR2 * P_Phi
           DPhi = DPhi + ZR2*ZR2*ZR2*ZR2*ZR2*ZR2 * P_DPhi
      call Compute_Phi(N,lattice,alphaN,a2N**7,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2**7 * P_Phi
           DPhi = DPhi - ZR2**7 * P_DPhi
      call Compute_Phi(N,lattice,alphaN,a2N**8,Dij,P_Phi,P_DPhi)
           Phi  = Phi  + ZR2**8 * P_Phi
           DPhi = DPhi + ZR2**8 * P_DPhi
      call Compute_Phi(N,lattice,alphaN,a2N**9,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2**9 * P_Phi
           DPhi = DPhi - ZR2**9 * P_DPhi
      call Compute_Phi(N,lattice,alphaN,a2N**10,Dij,P_Phi,P_DPhi)
           Phi  = Phi  + ZR2**10 * P_Phi
           DPhi = DPhi + ZR2**10 * P_DPhi
      call Compute_Phi(N,lattice,alphaN,a2N**11,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2**11 * P_Phi
           DPhi = DPhi - ZR2**11 * P_DPhi
c
      return
      end


      subroutine Compute_Phi(i,lattice,alphaN,aa,Dij,Phi,DPhi)
c
c  Compute Eu(r) and F(rho(r)) for calculation of phi-potential
c
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      common /eam1/ Ec(10,10),Re(10),B(10,10),A(10),beta(10,4),t(10,3)
      common /eam2/ RhoZ(10),Rho0_Bar(10),Z(10),Z2(10),A2(10)
      common /eam3/ Cmin(9,9,9),Cmax(9,9,9),Dcutoff,OverDeltaR,Eu_d(10)
      common /lj_2/ Rmin,Rsqmin,Rcutoff,OverDeltaRsq
      common /scrn/ scr_2nn
      character lattice*6
c
c  calculation of Eu(R) and dEdR for reference structure at r = aa*Dij
c
      Rfactor = aa * Dij / Re(i) - 1.d0
      astar  = alphaN * Rfactor
      astar3 = astar * astar * astar
      Eu = - Ec(i,i) * (1.d0 + astar + Eu_d(i) * astar3) * dexp(-astar)
      dEdR = Ec(i,i) * astar * alphaN / Re(i) * dexp(-astar) * aa
     &     * (1.d0 + Eu_d(i) * astar * astar - 3.d0 * Eu_d(i) * astar)
c
c  for calculation of Rho_Bar, F and dFdR for reference structure at r = aa*Dij
c
      R2factor = aa * a2(i) * Dij / Re(i) - 1.d0
      Rho_Bar =  Z(i) * RhoZ(i) * dexp(-beta(i,4)*Rfactor)
     &        + Z2(i) * RhoZ(i) * dexp(-beta(i,4)*R2factor)
      dRho0 = - beta(i,4) / Re(i) * aa 
     &      * ( Z(i) * RhoZ(i) * dexp(-beta(i,4)*Rfactor)
     &      +  a2(i) * Z2(i) * RhoZ(i) * dexp(-beta(i,4)*R2factor) )
      if(Rho_Bar.lt.1.0d-20) then
         Phi   = 0.d0 
         DPhi  = 0.d0 
         return
      endif
c
      G_GAMMA  = 1.0d0
      dRho_App = 0.0d0
c
c  add t(3) term in the case of HCP_A3 or DIA_A4 lattice
c
      if(lattice.eq.'HCP_A3' .or. lattice.eq.'DIA_A4') then
         if(lattice.eq.'HCP_A3') then
            Rho3  = 0.5d0 * RhoZ(i) * dexp(-beta(i,3)*Rfactor)
     &  - scr_2nn * dsqrt(2.d0) * RhoZ(i) * dexp(-beta(i,3)*R2factor)
            Rho3P = 4.d0 * Rho3 * Rho3 / 3.d0
            dRho3 = - beta(i,3) / Re(i) * aa * RhoZ(i)
     &            * ( 0.5d0 * dexp(-beta(i,3)*Rfactor)
     &  - scr_2nn * a2(i) * dsqrt(2.d0) * dexp(-beta(i,3)*R2factor) )
         elseif(lattice.eq.'DIA_A4') then
            Rho3  = RhoZ(i) * dexp(-beta(i,3)*Rfactor)
            Rho3P = 32.d0 * Rho3 * Rho3 / 9.d0
            dRho3 = - beta(i,3) / Re(i) * aa * Rho3
         endif
c
         Rho0 = Rho_Bar
         RhoFactor = Rho0 * Rho0
         if(RhoFactor.eq.0.d0) RhoFactor = 1.d-20
         GAMMA = t(i,3) * Rho3P / RhoFactor
            if(GAMMA.lt.-500.d0) then
               Phi   = 0.d0 
               DPhi  = 0.d0 
               return
            endif
         EXP_GAMMA = dexp(-GAMMA) 
         G_GAMMA = 2.d0 / ( 1.d0 + EXP_GAMMA )
         Rho_Bar = Rho0 * G_GAMMA
c        Rho_Bar = Rho0 * dsqrt ( 1.d0 + GAMMA )
c
         G_prime = 2.d0*EXP_GAMMA / (1.d0+EXP_GAMMA) / (1.d0+EXP_GAMMA)
         dRho_App = Rho0 * G_prime * 2.d0 * GAMMA * Rho0 / Rho3
     &            * ( dRho3 / Rho0 - Rho3 / Rho0 / Rho0 * dRho0 )
c        dRho_App = Rho0 / dsqrt( 1.d0 + GAMMA ) * GAMMA
c    &            * ( dRho3 / Rho3 - dRho0 / Rho0 )
      endif
c
c  add t(1), t(2) and t(3) terms in the case of DIMER species
c
      if(lattice.eq.'DIMER') then
         Rho1  = RhoZ(i) * dexp(-beta(i,1)*Rfactor)
         Rho2  = RhoZ(i) * dexp(-beta(i,2)*Rfactor)
         Rho3  = RhoZ(i) * dexp(-beta(i,3)*Rfactor)
         Rho1P = Rho1 * Rho1
         Rho2P = 2.d0 * Rho2 * Rho2 / 3.d0
         Rho3P = 2.d0 * Rho3 * Rho3 / 5.d0
         dRho1 = - beta(i,1) / Re(i) * aa * Rho1
         dRho2 = - beta(i,2) / Re(i) * aa * Rho2
         dRho3 = - beta(i,3) / Re(i) * aa * Rho3
c
         Rho0 = Rho_Bar
         RhoFactor = Rho0 * Rho0
         if(RhoFactor.eq.0.d0) RhoFactor = 1.d-20
         GAMMA = (t(i,1)*Rho1P + t(i,2)*Rho2P + t(i,3)*Rho3P) /RhoFactor
            if(GAMMA.lt.-500.d0) then
               Phi   = 0.d0 
               DPhi  = 0.d0 
               return
            endif
         EXP_GAMMA = dexp(-GAMMA) 
         G_GAMMA = 2.d0 / ( 1.d0 + EXP_GAMMA )
         Rho_Bar = Rho0 * G_GAMMA
c        Rho_Bar = Rho0 * dsqrt ( 1.d0 + GAMMA )
c
         G_prime = 2.d0*EXP_GAMMA / (1.d0+EXP_GAMMA) / (1.d0+EXP_GAMMA)
         dRho_App = Rho0 * G_prime * 2.d0 *
     1   ( t(i,1)*Rho1/Rho0 * (dRho1*Rho0-Rho1*dRho0) / RhoFactor 
     2   + t(i,2)*2.d0/3.d0*Rho2/Rho0*(dRho2*Rho0-Rho2*dRho0) /RhoFactor
     3   + t(i,3)*.4d0*Rho3/Rho0 * (dRho3*Rho0-Rho3*dRho0) / RhoFactor )
c        dRho_App = Rho0 / dsqrt( 1.d0 + GAMMA ) * 
c    1   ( t(i,1)*Rho1/Rho0 * (dRho1*Rho0-Rho1*dRho0) / RhoFactor 
c    2   + t(i,2)*2.d0/3.d0*Rho2/Rho0*(dRho2*Rho0-Rho2*dRho0) /RhoFactor
c    3   + t(i,3)*.4d0*Rho3/Rho0 * (dRho3*Rho0-Rho3*dRho0) / RhoFactor )
      endif
c
      if(rho_bar.eq.0.d0) rho_bar = 1.0d-20
      rr0 = Rho_Bar / Rho0_Bar(i)
      F = A(i) * Ec(i,i) * rr0 * dlog(rr0)

c     dFdR = A * Ec * ( 1 + dlog(rr0) ) / Rho0_Bar * dRho_dR
      dRho_dR = dRho0 * G_GAMMA + dRho_App
c     dRho_dR = dRho0 * dsqrt(1.d0 + GAMMA) + dRho_App
      dFdR = A(i) * Ec(i,i) * (1.d0+dlog(rr0)) / Rho0_Bar(i) * dRho_dR
c
c  Calculation of Phi and DPhi at r = aa * Dij
c
      Phi  = 2.d0 * ( Eu - F ) / Z(i)
      DPhi = 2.d0 * ( dEdR - dFdR ) / Z(i) / Dij
c
      return
      end

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c!END MODULE MEAMInitialize

