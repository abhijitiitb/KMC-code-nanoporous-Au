      
      !Initialize Verlet list variables for atom i
      ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(iatom) !range for ith atom
      r1=ra+1
      r2=ra+r0

      !s1=AL%AtomSpecies(i)
      s1=AtomSpecies(iatom)
      
      DO j=r1,r2
         
         drmag=VLdrmag(j)
         
         jatom=VLList(j)
         s2=AtomSpecies(jatom)
         smin=MIN(s1,s2) !this is automatically s1 because of the DO loop
         smax=MAX(s1,s2)
         section=((2*NSpecies-smin)*(smin-1))/2+smax
            
         !Pair potential term
         dxinv=IncrPPinv(section)
         mincutoff=RangePP(2*section-1)
      
         tmp=(drmag-mincutoff)*dxinv
         
         k0=FLOOR(tmp)
         
         tgphi=0._dp
         
         IF (k0<2) THEN
         
            IF (CatchOutOfBound) THEN
               WRITE(6,*) TRIM(TxtHeader)//"$Err>> EAM index below PP array lower bound"
               WRITE(UNIT=*,FMT='("...too low drmag:",ES14.5)') drmag
               WRITE(6,*) "Index in interp table=", k0
               WRITE(UNIT=6,FMT='("Main atom coord:",3f15.3)') AtomCoord(3*iatom-2:3*iatom)
               WRITE(6,*) "DRMAG="
               DO pos=r1,r2
                  jatom=VLList(pos)
                  WRITE(UNIT=6,FMT='("Rad:",f15.3," Interac atom coord:",3f15.3)') VLdrmag(pos), &
                     AtomCoord(3*jatom-2:3*jatom)
               END DO
               errorstatus=1
               STOP
               !RETURN
            ELSE
               k=(section-1)*MaxEAMTableSize+1
               !linear extrapolation
               EnergyPP=EnergyPP+ PP(k) + tmp*(PP(k+1)-PP(k))
               tgphi=(PP(k+1)-PP(k))*dxinv/drmag !dphi/drij/drmag
               !WRITE(6,*) TRIM(TxtHeader)//"S3-1:",tgphi,section,drmag
            END IF
            
         ELSEIF (drmag<RangePP(2*section)) THEN
         
            p=MIN(tmp-REAL(k0,dp),1._dp)
            k=MIN(k0,SizePP(section)-1)+(section-1)*MaxEAMTableSize
            EnergyPPTerm= PP(k)+0.5_dp*p*(PP(k+1)-PP(k-1)  +p*(PP(k+1)-PP(k)-PP(k)+PP(k-1)))
            IF (.NOT. VLListDomainAtom(j)) EnergyPPTerm=EnergyPPTerm*0.5_dp !only half contribution added
            EnergyPP=EnergyPP+EnergyPPTerm
            
            k=MIN(k0,SizePP(section)-2)+(section-1)*MaxEAMTableSize
            gi=PP(k+1)-PP(k-1)
            gi1=PP(k+2)-PP(k)
            slope=0.5_dp*(gi+p*(gi1-gi))*dxinv
            tgphi=slope/drmag
         
         END IF
         
         !Density term -- role of jatom on iatom
         section=s2
         dxinv=IncrDTinv(section)
         mincutoff=RangeDT(2*section-1)
         tmp=(drmag-mincutoff)*dxinv
         k0=FLOOR(tmp)
         
         tgrho=0._dp
         
         IF (k0<2) THEN
            !IF (CatchOutOfBound) THEN !this must have already been checked in Step1.f90
            
            k=(section-1)*MaxEAMTableSize+1
            !linear extrapolation
            tgrho=(DT(k+1)-DT(k))*dxinv/drmag !drho/drij/drmag
            !WRITE(6,*) TRIM(TxtHeader)//"S3-2:",tgrho,section,drmag
         ELSEIF (drmag<RangeDT(2*section)) THEN
            p=MIN(tmp-REAL(k0,dp),1._dp)
         
            k=MIN(k0,SizeDT(section)-2)+(section-1)*MaxEAMTableSize
            gi=DT(k+1)-DT(k-1)
            gi1=DT(k+2)-DT(k)
            slope=0.5_dp*(gi+p*(gi1-gi))*dxinv
            tgrho=slope/drmag
         
         !use of jatom here requires some discussion. if jatom lies inside the domain then the forces are anyways correct
         !however if jatom is outside the domain, then we have not computed its density. because of the periodic boundary conditions
         !we can take density at jatom atom as the same as the jatom in the domain
         !grad=tgrad*VLdr(3*j-2:3*j)
         END IF
         
         tgrad=tgphi+DensityArray1(jatom)*tgrho
         
         !Density term -- role of iatom on jatom
         section=s1
         dxinv=IncrDTinv(section)
         mincutoff=RangeDT(2*section-1)
         tmp=(drmag-mincutoff)*dxinv
         k0=FLOOR(tmp)
         tgrho=0._dp
         IF (k0<2) THEN
            !IF (CatchOutOfBound) THEN !this must have already been checked in Step1.f90
            
            k=(section-1)*MaxEAMTableSize+1
            !linear extrapolation
            tgrho=(DT(k+1)-DT(k))*dxinv/drmag !drho/drij/drmag
            !WRITE(6,*) TRIM(TxtHeader)//"S3-3:",tgrho,section,drmag
         ELSEIF (drmag<RangeDT(2*section)) THEN
            p=MIN(tmp-REAL(k0,dp),1._dp)
         
            k=MIN(k0,SizeDT(section)-2)+(section-1)*MaxEAMTableSize
            gi=DT(k+1)-DT(k-1)
            gi1=DT(k+2)-DT(k)
            slope=0.5_dp*(gi+p*(gi1-gi))*dxinv
            tgrho=slope/drmag
         END IF
         
         tgrad=tgrad+DensityArray1(iatom)*tgrho
         
         grad=tgrad*VLdr(3*j-2:3*j)
         
         AtomForce(3*iatom-2:3*iatom)=AtomForce(3*iatom-2:3*iatom)+grad
         IF (VLListDomainAtom(j)) AtomForce(3*jatom-2:3*jatom)=AtomForce(3*jatom-2:3*jatom)-grad !then atom j lies in the domain
         
      END DO
