      !Computes the forces using interpolation
      
      !Initialize Verlet list variables for atom iatom
      ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(iatom) !range for ith atom
      r1=ra+1
      r2=ra+r0

      s1=AtomSpecies(iatom)
      
      DO j=r1,r2
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
         EnergyPPTerm= PP(k)+0.5_dp*p*(PP(k+1)-PP(k-1)  +p*(PP(k+1)-PP(k)-PP(k)+PP(k-1)))
         IF (.NOT. VLListDomainAtom(j)) EnergyPPTerm=EnergyPPTerm*0.5_dp !only half contribution added
         EnergyPP=EnergyPP+EnergyPPTerm
                  
         !xxxxxxxxxxxxxxxxxxxx  DENSITY TERM FOR ATOM i xxxxxxxxxxxxxxxxxxx
         section=s2
         dxinv=IncrDTinv(section)
         mincutoff=RangeDT(2*section-1)
         tmp=(drmag-mincutoff)*dxinv
         k=FLOOR(tmp) !table index for reading density and pair potential, use this instead of abi
         IF (k<2) THEN
            WRITE(6,*) "$Err>> EAM index below DT array lower bound"
            WRITE(UNIT=6,FMT='("...too low drmag:",ES14.5)') drmag
            WRITE(6,*) "Index in interp table=", k
            WRITE(UNIT=6,FMT='("Main atom coord:",3f15.3)') AtomCoord(3*iatom-2:3*iatom)
            WRITE(6,*) "DRMAG="
            DO pos=r1,r2
               jatom=VLList(pos)
               WRITE(UNIT=6,FMT='("Rad:",f15.3," Interac atom coord:",3f15.3)') VLdrmag(pos), &
                 AtomCoord(3*jatom-2:3*jatom)
            END DO
            STOP
         END IF
             
         p=MIN(tmp-REAL(k,dp),1._dp)
         k=MIN(k,SizeDT(section)-1)+(section-1)*MaxEAMTableSize
         tmp=DT(k)+ 0.5_dp*p*(DT(k+1)-DT(k-1)   +p*(DT(k+1)-DT(k)-DT(k)+DT(k-1)))
         DensityArray(iatom)=DensityArray(iatom)+ tmp
         
         !xxxxxxxxxxxxxxxxxxxx  DENSITY TERM FOR ATOM i1 xxxxxxxxxxxxxxxxxxx
         IF (VLListDomainAtom(j)) THEN !atom i1 is a domain atom
            section=s1
            dxinv=IncrDTinv(section)
            mincutoff=RangeDT(2*section-1)
            tmp=(drmag-mincutoff)*dxinv
            k=FLOOR(tmp) !table index for reading density and pair potential, use this instead of abi
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
            tmp=DT(k)+ 0.5_dp*p*(DT(k+1)-DT(k-1)   +p*(DT(k+1)-DT(k)-DT(k)+DT(k-1)))
            DensityArray(jatom)=DensityArray(jatom)+ tmp
         END IF
      END DO
      
