      IF (MDSiml%LangevinCoeff==0._dp) THEN
         WRITE(*,*) "$Err>> Langevin coefficient has not been setup"
         STOP
      END IF
      
      IF (MDSiml%dt>1.e-12_dp .OR. MDSiml%dt<1.e-18_dp) THEN
         WRITE(*,*) "$Err>> MD time step has been improperly setup"
         STOP
      END IF
      
      NULLIFY(sim)
      NULLIFY(vf)
      NULLIFY(vc)
      NULLIFY(ap)
      NULLIFY(ac)
      NULLIFY(sigr)
      NULLIFY(sigv)
      NULLIFY(g1)
      NULLIFY(g2)
      NULLIFY(crv)
      NULLIFY(c)
      
      CALL MakeSize(sim,3*NAtoms)
      CALL MakeSize(vf,3*NAtoms)
      CALL MakeSize(vc,3*NAtoms)
      CALL MakeSize(ap,3*NAtoms)
      CALL MakeSize(ac,3*NAtoms)
      CALL MakeSize(sigr,3*NAtoms)
      CALL MakeSize(sigv,3*NAtoms)
      CALL MakeSize(g1,3*NAtoms)
      CALL MakeSize(g2,3*NAtoms)
      CALL MakeSize(crv,3*NAtoms)
      CALL MakeSize(c,3*NAtoms)
      
      ldt=MDSiml%LangevinCoeff*dt !dimensionless langevin coeff (s * 1/s)
      ldtinv=1._dp/ldt
      dt2=dt*dt
      m=AL%AtomMass
      
      TemperatureRamp=MDSiml%TemperatureRamp
      
      IF (MDSiml%LangevinCoeff>0._dp) THEN
         c0=exp(-ldt)
         c1=(1.0_dp-c0)/ldt
         c2=(1.0_dp-c1)/ldt
      ELSE
         c0=1._dp
         c1=1._dp
         c2=0.5_dp
      END IF
      
      sim=SQRT(kboltzmann*MDSiml%Temperature/AL%AtomMass) !in velocity units sqrt(eV/amu)
      
      !variances & correlations from A&T p.262
      a1=2._dp-ldtinv*(3._dp-4._dp*c0+c0*c0)
      sigr=fv*dt*sim*SQRT(a1/ldt) !std dev for position in Ang.
      sigv=sim*SQRT(1._dp-c0*c0) !Eq. 9.23b in SQRT(eV/amu)
      !crv=a2/SQRT(ldt)/SQRT(a1*(1._dp-c0*c0)) !Eq. 9.23c 
      
      IF (a1<0._dp) THEN
         WRITE(6,*) "$Err>> sigr in Langevin dynamics is a complex number"
         STOP
      END IF
      
      a2=(1._dp-c0)*(1._dp-c0)
      crv=fv*dt*sim*sim*a2/ldt/sigr/sigv
         !crv ends up being a scalar quantity
         
      IF (MAXVAL(crv)>1.0_dp .OR. MINVAL(crv)<0.0_dp) THEN
         WRITE(*,*) "$Err>> Correlation coefficient beyond range [0,1]"
         STOP
      END IF
      c=SQRT(1.0_dp-crv*crv)
