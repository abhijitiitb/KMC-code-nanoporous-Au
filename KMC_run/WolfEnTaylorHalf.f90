      !Computes the forces using interpolation
      
      !Initialize Verlet list variables for atom iatom
      ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(iatom) !range for ith atom
      r1=ra+1
      r2=ra+r0

      z1=AtomCharge(iatom)

      DO j=r1,r2
         jatom=VLList(j)
         z2=AtomCharge(jatom)
         drmag=VLdrmag(j)

         !xxxxxxxxxxxxxxxxxxxx  PAIR POTENTIAL TERM  xxxxxxxxxxxxxxxxxxx
         EnergyPP=EnergyPP+ z1*z2*(ERFC(alfa*drmag)/drmag-earc)

      END DO
      
      
