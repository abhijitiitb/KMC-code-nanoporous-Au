      
      !Initialize Verlet list variables for atom iatom
      ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(iatom) !range for ith atom
      r1=ra+1
      r2=ra+r0

      s1=AtomSpecies(iatom)
      
      DO j=r1,r2
         
         drmag=VLdrmag(j)
         jatom=VLList(j)
         s2=AtomSpecies(jatom)
         smin=MIN(s1,s2)
         smax=MAX(s1,s2)
            
         !Pair potential term
         section=((2*NSpecies-smin)*(smin-1))/2+smax
         dxinv=IncrPPinv(section)
         mincutoff=RangePP(2*section-1)
         tmp=(drmag-mincutoff)*dxinv
         k0=FLOOR(tmp)
         
         IF (k0<2) THEN
            WRITE(6,*) "$Err>> EAM index below PP array lower bound"
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
            RETURN
         END IF
         
         p=MIN(tmp-REAL(k0,dp),1._dp)
         k=MIN(k0,SizePP(section)-1)+(section-1)*MaxPairPotentialTableSize
         
         EnergyPPTerm= PP(k)+0.5_dp*p*(PP(k+1)-PP(k-1)  +p*(PP(k+1)-PP(k)-PP(k)+PP(k-1)))
         IF (.NOT. VLListDomainAtom(j)) EnergyPPTerm=EnergyPPTerm*0.5_dp !only half contribution added
         EnergyPP=EnergyPP+EnergyPPTerm
            
         k=MIN(k0,SizePP(section)-2)+(section-1)*MaxPairPotentialTableSize
         gi=PP(k+1)-PP(k-1)
         gi1=PP(k+2)-PP(k)
         slope=0.5_dp*(gi+p*(gi1-gi))*dxinv
         tgphi=slope/drmag
         
         tgrad=tgphi
         grad=tgrad*VLdr(3*j-2:3*j)
         AtomForce(3*iatom-2:3*iatom)=AtomForce(3*iatom-2:3*iatom)+grad
         IF (VLListDomainAtom(j)) AtomForce(3*jatom-2:3*jatom)=AtomForce(3*jatom-2:3*jatom)-grad !atom lies in domain
         
      END DO
