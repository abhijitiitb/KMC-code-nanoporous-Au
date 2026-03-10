      
      !Initialize Verlet list variables for atom iatom
      ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(iatom) !range for ith atom
      r1=ra+1
      r2=ra+r0

      z1=AtomCharge(iatom)
      
      DO j=r1,r2
         
         jatom=VLList(j)
         drmag=VLdrmag(j)
         z2=AtomCharge(jatom)
            
         !Pair potential term
         ear=ERFC(alfa*drmag)/drmag
         EnergyPP=EnergyPP+ z1*z2*(ear-earc)
         
         !slope1=z1*z2*(ear/drmag/drmag+2._dp*alfa/SQRT(PI)*EXP(-alfa*alfa*drmag*drmag)/drmag)
         tgrad=slope
         !slope2=-z1*z2*(earc/Rc/Rc+2._dp*alfa/SQRT(PI)*EXP(-alfa*alfa*Rc*Rc)/Rc)
         WRITE(6,*) "Error here"
         STOP
         tgrad=tgrad+slope
         tgrad=tgrad/drmag
         
         grad=tgrad*VLdr(3*j-2:3*j)
         AtomForce(3*iatom-2:3*iatom)=AtomForce(3*iatom-2:3*iatom)+grad
         AtomForce(3*jatom-2:3*jatom)=AtomForce(3*jatom-2:3*jatom)-grad
         
      END DO
