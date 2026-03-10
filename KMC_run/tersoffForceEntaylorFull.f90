!Computes the Tersoff energy and forces based on Tersoff, PRB, 37, 12, 6991 (1998)
!Computes the forces using interpolation with full list
      
      ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(iatom) !range for ith atom
      r1=ra+1
      r2=ra+r0

      s1=AtomSpecies(iatom)
      
      DO j=r1,r2 !this will take care of the pair and many-body terms in one go
      
         jatom=VLList(j)
         s2=AtomSpecies(jatom)
         
         smin=MIN(s1,s2)
         smax=MAX(s1,s2)
         drij=VLdr(3*j-2:3*j)  !vector between atoms i and j
         drmagij=VLdrmag(j)   ! spacing between atoms i and j
         
         !xxxxxxxxxxxxxxxxxxxx  PAIR POTENTIAL TERM  xxxxxxxxxxxxxxxxxxx
         sectionij=((2*NSpecies-smin)*(smin-1))/2+smax
         Rangefcij=Rangefc(sectionij)
         IF (drmagij<Rangefcij .OR. iatom<jatom) THEN !below cut-off or valid pair
         
            Aij=A(sectionij)
            Bij=B(sectionij)
            lam1ij=lam1(sectionij)
            lam2ij=lam2(sectionij)
            lam3ij=lam3(sectionij)
            betaij=beta(sectionij)
            Rij=R(sectionij)
            RDij=RD(sectionij)
            cij=c(sectionij)
            dij=d(sectionij)
            hij=h(sectionij)
            nij=n(sectionij)
            
            !energy
            frepl=Aij*EXP(-lam1ij*drmagij)
            fattr=-Bij*EXP(-lam2ij*drmagij)
            
            drt=PIO2*(MAX(drmagij,Rij-RDij)-Rij)/RDij
            fcut=0.5_dp-0.5_dp*SIN(drt)
            
            !force
            dfrepl=-lam1ij*frepl !dertivative dfrepl/drij
            dfattr=-lam2ij*fattr !dfr/drij
            dfcut=-0.5_dp*PIO2*COS(drt)/RDij !dfc/drij
         
            b_ij=0._dp
            b_ji=0._dp
            
            INCLUDE "tersoff-bondorder.f90"
            
            EnergyTersoff=EnergyTersoff+fcut*(frepl+0.5_dp*(b_ij+b_ji)*fattr)
            fterm = (dfcut*(frepl+0.5_dp*(b_ij+b_ji)*fattr) + fcut*(dfrepl + 0.5_dp*(b_ij+b_ji)*dfattr))*drij_rixi
            ForceTersoff_i = ForceTersoff_i - fterm - fcut*(0.5_dp*(db_ij_drixi+db_ji_drixi)*fattr)
            ForceTersoff_j = ForceTersoff_j + fterm - fcut*(0.5_dp*(db_ij_drjxi+db_ji_drjxi)*fattr)
            
            write (*,*) frepl,fattr,fcut,b_ij,b_ji,drmagij
         END IF
      END DO
write (*,*) EnergyTersoff
