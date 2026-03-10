      ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(iatom) !range for ith atom
      r1=ra+1
      r2=ra+r0

      s1=AtomSpecies(iatom)

      DO j=r1,r2
         
         drmag=VLdrmag(j)
         jatom=VLList(j)
         
         s2=AtomSpecies(jatom)
                  
         !xxxxxxxxxxxxxxxxxxxx  DENSITY TERM  xxxxxxxxxxxxxxxxxxx
         !contribution of jatom to the density at iatom
         section=s2
         dxinv=IncrDTinv(section)
         mincutoff=RangeDT(2*section-1)
         
         tmp=(drmag-mincutoff)*dxinv
         k=FLOOR(tmp)
!COMMENT 1: AC--  k=0 corresponds to mincutoff - so we should add 1 to k at this point 
!to bring it to the correct table (array)  position. Now I have added 1
!The force term requires k=0 to be k=0 and not k=1 so I have removed the +1 added - 13/07/2012 
         IF (k<2) THEN ! (k<1) THEN !see Comment 1-- k<2 was required when 1 not added
            IF (CatchOutOfBound) THEN
               WRITE(6,*) TRIM(TxtHeader)//"$Err>> EAM index below DT array lower bound"
               WRITE(UNIT=*,FMT='("...too low drmag:",ES14.5)') drmag
               WRITE(6,*) "Index in interp table=", k
               WRITE(UNIT=6,FMT='("Main atom coord:",3f15.3)') AtomCoord(3*iatom-2:3*iatom)
               WRITE(6,*) "DRMAG="
               DO pos=r1,r2
                  jatom=VLList(pos)
                  WRITE(UNIT=6,FMT='("Rad:",f15.3," Interac atom coord:",3f15.3)') VLdrmag(pos), &
                    AtomCoord(3*jatom-2:3*jatom)
               END DO
               errorstatus=2
               STOP
               !RETURN
            ELSE  !we shall (linear) extrapolate -- note that k<2
               !use tmp to find the function value - 11/07/12
               !note we don't need to increment tmp by 1 - 11/07/12
               !let t denote the (real/fractional) array position corresponding to tmp then --
               ! (f(t)-f(1))/(t-1)=f(2)-f(1) - where f(1) is fn value at mincutoff
               !or f(t)=f(1)+tmp*(f(2)-f(1)) - 11/07/12
               !more generally f(t)=f(1)+(tmp+1+W-1-W)*(f(2)-f(1)) -- where W=(section-1)*MaxEAMTableSize
               k=(section-1)*MaxEAMTableSize+1
               DensityArray(iatom)=DensityArray(iatom)+ DT(k) + tmp*(DT(k+1)-DT(k))
               !WRITE(6,*) TRIM(TxtHeader)//"S1-1:",tmp,DT(k) + tmp*(DT(k+1)-DT(k)),DT(k),DT(k+1),section,drmag
            END IF
         ELSE !use the tables

            p=MIN(tmp-REAL(k,dp),1._dp)
            k=MIN(k,SizeDT(section)-1)+(section-1)*MaxEAMTableSize !+1 !1 added -- see Comment 1
            tmp=DT(k)+ 0.5_dp*p*(DT(k+1)-DT(k-1)   +p*(DT(k+1)-DT(k)-DT(k)+DT(k-1)))
            DensityArray(iatom)=DensityArray(iatom)+ tmp
            
         END IF
         !xxxxxxxxxxxxxxxxxx
         !contribution of iatom to the density at jatom
         IF (VLListDomainAtom(j)) THEN !then atom jatom lies in the domain
            section=s1
            dxinv=IncrDTinv(section)
            mincutoff=RangeDT(2*section-1)
         
            tmp=(drmag-mincutoff)*dxinv
            k=FLOOR(tmp)
            IF (k<2) THEN !(k<1) THEN !AC-- see Comment 1 - k<2 required when 1 was not added - 11/07/12
               IF (CatchOutOfBound) THEN
                  WRITE(6,*) TRIM(TxtHeader)//"$Err>> EAM index below DT array lower bound"
                  WRITE(UNIT=*,FMT='("...too low drmag:",ES14.5)') drmag
                  WRITE(6,*) "Index in interp table=", k
                  WRITE(UNIT=6,FMT='("Main atom coord:",3f15.3)') AtomCoord(3*iatom-2:3*iatom)
                  WRITE(6,*) "DRMAG="
                  DO pos=r1,r2
                     jatom=VLList(pos)
                     WRITE(UNIT=6,FMT='("Rad:",f15.3," Interac atom coord:",3f15.3)') VLdrmag(pos), &
                        AtomCoord(3*jatom-2:3*jatom)
                  END DO
                  errorstatus=2
                  STOP
                  !RETURN
               ELSE  !we shall (linear) extrapolate -- note that k<2
                  !use tmp to find the function value - 11/07/12
                  !note we don't need to increment tmp by 1 - 11/07/12
                  !let t denote the (real/fractional) array position corresponding to tmp then --
                  ! (f(t)-f(1))/(t-1)=f(2)-f(1) - where f(1) is fn value at mincutoff
                  !or f(t)=f(1)+tmp*(f(2)-f(1)) - 11/07/12
                  k=(section-1)*MaxEAMTableSize+1
                  DensityArray(jatom)=DensityArray(jatom)+ DT(k) + tmp*(DT(k+1)-DT(k))
               !WRITE(6,*) TRIM(TxtHeader)//"S1-2:",tmp,DT(k) + tmp*(DT(k+1)-DT(k)),DT(k),DT(k+1),section,drmag
                  
               ENDIF   
            ELSEIF (drmag<RangeDT(2*section)) THEN !use the tables

               p=MIN(tmp-REAL(k,dp),1._dp)
               k=MIN(k,SizeDT(section)-1)+(section-1)*MaxEAMTableSize !+1 !AC-- 1 added -- see comment 1-- 11/07/12
               tmp=DT(k)+ 0.5_dp*p*(DT(k+1)-DT(k-1)   +p*(DT(k+1)-DT(k)-DT(k)+DT(k-1)))
               DensityArray(jatom)=DensityArray(jatom)+ tmp
            END IF
         END IF
      END DO
