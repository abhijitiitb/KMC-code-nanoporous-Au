      !Computes the forces using interpolation
      
      !Initialize Verlet list variables for atom i
      ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(iatom) !range for ith atom
      r1=ra+1
      r2=ra+r0

      s1=AtomSpecies(iatom)

      !First loop to calculate the density and the pair potential term
      Density=0._dp
      
      DO j=r1,r2 !over all neighbors
         jatom=VLList(j)
         s2=AtomSpecies(jatom)
         drmag=VLdrmag(j)
            
         !xxxxxxxxxxxxxxxxxxxx  PAIR POTENTIAL TERM  xxxxxxxxxxxxxxxxxxx
         smin=MIN(s1,s2) !this is automatically s1 because of the DO loop
         smax=MAX(s1,s2)
         section=((2*NSpecies-smin)*(smin-1))/2+smax
         dxinv=IncrPPinv(section)
         mincutoff=RangePP(2*section-1)
         tmp=(drmag-mincutoff)*dxinv
         k=FLOOR(tmp) !table index for reading density and pair potential, use this instead of abi
         IF (k<2) THEN
            WRITE(6,*) TRIM(TxtHeader)//"$Err>> EAM index below PP array lower bound"
            WRITE(UNIT=6,FMT='("...too low drmag:",ES14.5)') drmag
            WRITE(6,*) TRIM(TxtHeader)//"Index in interp table=", k
            WRITE(UNIT=6,FMT='("Main atom index:",I4,"           Main atom coord:",3f10.3)') &
               iatom,AtomCoord(3*iatom-2:3*iatom)
            WRITE(6,*) TRIM(TxtHeader)//"DRMAG="
            DO pos=r1,r2
               jatom=VLList(pos)
               WRITE(UNIT=6,FMT='("Atom index:",I4,"  Rad:",f10.3,"    Interac atom coord:",3f10.3)') &
                 jatom,VLdrmag(pos),AtomCoord(3*jatom-2:3*jatom)
            END DO
            STOP
         END IF
         
         p=MIN(tmp-REAL(k,dp),1._dp)
         k=MIN(k,SizePP(section)-1)+(section-1)*MaxEAMTableSize
         EnergyPP=EnergyPP+PP(k)+ 0.5_dp*p*(PP(k+1)-PP(k-1)  +p*(PP(k+1)-PP(k)-PP(k)+PP(k-1)))
         
         !xxxxxxxxxxxxxxxxxxxx  DENSITY TERM  xxxxxxxxxxxxxxxxxxx
         section=s2
         dxinv=IncrDTinv(section)
         mincutoff=RangeDT(2*section-1)
         tmp=(drmag-mincutoff)*dxinv
         k=FLOOR(tmp)
         IF (k<2) THEN
            WRITE(6,*) TRIM(TxtHeader)//"$Err>> EAM index below DT array lower bound"
            WRITE(UNIT=6,FMT='("...too low drmag:",ES14.5)') drmag
            WRITE(6,*) TRIM(TxtHeader)//"Index in interp table=", k
            WRITE(UNIT=6,FMT='("Main atom coord:",3f15.3)') AtomCoord(3*iatom-2:3*iatom)
            WRITE(6,*) TRIM(TxtHeader)//"DRMAG="
            DO pos=r1,r2
               jatom=VLList(pos)
               WRITE(UNIT=6,FMT='("Rad:",f15.3," Interac atom coord:",3f15.3)') VLdrmag(pos), &
                 AtomCoord(3*jatom-2:3*jatom)
            END DO
            STOP
         END IF
            
         p=MIN(tmp-REAL(k,dp),1._dp)
         k=MIN(k,SizeDT(section)-1)+(section-1)*MaxEAMTableSize
         Density=Density+DT(k)+0.5_dp*p*(DT(k+1)-DT(k-1)  +p*(DT(k+1)-DT(k)-DT(k)+DT(k-1)))
      END DO
      
      !Compute the embedding energy
      section=s1
      dxinv=IncrEEinv(section)
      mincutoff=RangeEE(2*section-1)
      tmp=(Density-mincutoff)*dxinv
      k=FLOOR(tmp)
      IF (k>SizeEE(section)-2) THEN
         WRITE(6,*) "$Err[0]>> Increase range of density array"
         WRITE(6,*) "...embedding energy term"
         WRITE(UNIT=6,FMT='("... density encountered:",ES14.5)') Density
         STOP
      END IF
         
      IF (k<=1) CYCLE  !free atom
      
      p=MIN(tmp-REAL(k,dp),1._dp)
      EnergyEmbedding=EnergyEmbedding+EE(k)+0.5_dp*p*(EE(k+1)-EE(k-1)  +p*(EE(k+1)-EE(k)-EE(k)+EE(k-1)))
