      !Computes the forces using interpolation
      
      !Initialize Verlet list variables for atom i
      ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(iatom) !range for ith atom
      r1=ra+1
      r2=ra+r0

      s1=AtomSpecies(iatom)

      !First loop to calculate the density and the pair potential term
      Density=0._dp
      
      DO j=r1,r2
      
         jatom=VLList(j)
         drmag=VLdrmag(j)
         
         s2=AtomSpecies(jatom)
         
         !xxxxxxxxxxxxxxxxxxxx  PAIR POTENTIAL TERM  xxxxxxxxxxxxxxxxxxx
         smin=MIN(s1,s2)
         smax=MAX(s1,s2)
         section=((2*NSpecies-smin)*(smin-1))/2+smax
         dxinv=IncrPPinv(section)
         mincutoff=RangePP(2*section-1)
         
         tmp=(drmag-mincutoff)*dxinv
         k0=FLOOR(tmp) !table index for reading density and pair potential, use this instead of abi
         IF (k0<2) THEN
            WRITE(6,*) TRIM(TxtHeader)//"$Err>> EAM index below PP array lower bound"
            WRITE(UNIT=6,FMT='("...too low drmag:",ES14.5)') drmag
            WRITE(6,*) "Index in interp table=", k
            WRITE(UNIT=6,FMT='("Main atom index:",I4,"           Main atom coord:",3f10.3)') &
               iatom,AtomCoord(3*iatom-2:3*iatom)
            WRITE(6,*) "DRMAG="
            DO pos=r1,r2
               jatom=VLList(pos)
               WRITE(UNIT=6,FMT='("Atom index:",I4,"  Rad:",f10.3,"    Interac atom coord:",3f10.3)') &
                 jatom,VLdrmag(pos),AtomCoord(3*jatom-2:3*jatom)
            END DO
            errorstatus=1  !pair potential is below range
            RETURN
         END IF
         
         p=MIN(tmp-REAL(k0,dp),1._dp)
         k=MIN(k0,SizePP(section)-1)+(section-1)*MaxEAMTableSize
         EnergyPP=EnergyPP+PP(k)+ 0.5_dp*p*(PP(k+1)-PP(k-1)  +p*(PP(k+1)-PP(k)-PP(k)+PP(k-1)))
         
         k=MIN(k0,SizePP(section)-2)+(section-1)*MaxEAMTableSize
         gi=PP(k+1)-PP(k-1)
         gi1=PP(k+2)-PP(k)
         tgrad=0.5_dp*(gi+p*(gi1-gi))*dxinv !gradient of the pair potential between iatom and jatom
         AtomForce(3*iatom-2:3*iatom)=AtomForce(3*iatom-2:3*iatom)+(tgrad/drmag)*VLdr(3*j-2:3*j)
         
         !xxxxxxxxxxxxxxxxxxxx  DENSITY TERM  xxxxxxxxxxxxxxxxxxx
         section=s2
         dxinv=IncrDTinv(section)
         mincutoff=RangeDT(2*section-1)
         tmp=(drmag-mincutoff)*dxinv
         k0=FLOOR(tmp) !table index for reading density and pair potential, use this instead of abi
         IF (k0<2) THEN
            WRITE(*,*) TRIM(TxtHeader)//"$Err>> EAM index below DT array lower bound"
            WRITE(UNIT=*,FMT='("...too low drmag:",ES14.5)') drmag
            WRITE(*,*) "Index in interp table=", k
            WRITE(UNIT=6,FMT='("Main atom coord:",3f15.3)') AtomCoord(3*iatom-2:3*iatom)
            WRITE(6,*) "DRMAG="
            DO pos=r1,r2
               jatom=VLList(pos)
               WRITE(UNIT=6,FMT='("Rad:",f15.3," Interac atom coord:",3f15.3)') VLdrmag(pos), &
                 AtomCoord(3*jatom-2:3*jatom)
            END DO
            errorstatus=2
            RETURN
         END IF
         
         p=MIN(tmp-REAL(k0,dp),1._dp)
         k=MIN(k0,SizeDT(section)-1)+(section-1)*MaxEAMTableSize
         Density=Density+DT(k)+0.5_dp*p*(DT(k+1)-DT(k-1)  +p*(DT(k+1)-DT(k)-DT(k)+DT(k-1)))
         
      END DO
      
      !Compute the embedding energy
      section=s1
      dxinv=IncrEEinv(section)
      mincutoff=RangeEE(2*section-1)
      tmp=(Density-mincutoff)*dxinv
      k0=FLOOR(tmp)
      
      IF (k0<=1) CYCLE  !free atom  -- the steps that follow are not needed for a free atom
      
      p=MIN(tmp-REAL(k0,dp),1._dp)
      k=MIN(k0,SizeEE(section)-1)+(section-1)*MaxEAMTableSize
      EnergyEmbedding=EnergyEmbedding+EE(k)+0.5_dp*p*(EE(k+1)-EE(k-1)  +p*(EE(k+1)-EE(k)-EE(k)+EE(k-1)))
      
      k=MIN(k0,SizeEE(section)-2)+(section-1)*MaxEAMTableSize
      gi=EE(k+1)-EE(k-1)
      gi1=EE(k+2)-EE(k)
      ro1=0.5_dp*(gi+p*(gi1-gi))*dxinv  !dF/drho for atom i
   
      !Compute the forces
      DO j=r1,r2
         jatom=VLList(j)
         s2=AtomSpecies(jatom)
         drmag=VLdrmag(j)
         
         !xxxxxxxxxxxxxxxxxxxx  DENSITY TERM  xxxxxxxxxxxxxxxxxxx
         section=s2
         dxinv=IncrDTinv(section)
         mincutoff=RangeDT(2*section-1)
         tmp=(drmag-mincutoff)*dxinv
         k0=FLOOR(tmp)
         
         p=MIN(tmp-REAL(k0,dp),1._dp)
         k=MIN(k0,SizeDT(section)-2)+(section-1)*MaxEAMTableSize
         gi=DT(k+1)-DT(k-1)
         gi1=DT(k+2)-DT(k)
         tgrad=0.5_dp*(gi+p*(gi1-gi))*dxinv !gradient of the pair potential between i and jatom
         grad=(ro1*tgrad/drmag)*VLdr(3*j-2:3*j)
         AtomForce(3*iatom-2:3*iatom)=AtomForce(3*iatom-2:3*iatom)+grad
         !IF (VLListDomainAtom(j)) -- in the case of full list if an atom i can see an image of another atom j then 
         ! atom j can see the image of atom i - so we directly include it
         AtomForce(3*jatom-2:3*jatom)=AtomForce(3*jatom-2:3*jatom)-grad !atom lies in domain
         !WRITE(*,*) "force is incorrect"
         !WRITE(*,*) "see example 8c when 1x1x1 unit cell is used"
         !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      END DO
