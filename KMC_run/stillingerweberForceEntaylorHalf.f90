!Computes the forces using interpolation with half list

!Issue that needs to be resolved: the triplets might not be complete
!iatom's neighborlist contains jatom and katom
!but jatom's neighborlist does not contains iatom -- so a triplet such as iatom-jatom-latom cannot be formed
! where latom is present in the neighborlist of jatom
      STOP
      
      ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(iatom) !range for ith atom
      r1=ra+1
      r2=ra+r0

      s1=AtomSpecies(iatom)
      
      !first loop to calculate force for pair porential term
      !WRITE(6,*) "iatom:",iatom,s1
!WRITE(6,*) "Atom position i:",AL%AtomCoord(3*iatom-2:3*iatom)
      DO j=r1,r2 !this will take care of the pair and many-body terms in one go
      
         jatom=VLList(j)
         s2=AtomSpecies(jatom)
!WRITE(6,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
!WRITE(6,*) "jatom:",jatom,s2
         
         smin=MIN(s1,s2)
         smax=MAX(s1,s2)
         drij=VLdr(3*j-2:3*j)  !vector between atoms i and j
         drmagij=VLdrmag(j)   ! spacing between atoms i and j
         
!WRITE(6,*) "drij:",drij,drmagij
!WRITE(6,*) "Atom position j:",AL%AtomCoord(3*jatom-2:3*jatom)
         
         !xxxxxxxxxxxxxxxxxxxx  PAIR POTENTIAL TERM  xxxxxxxxxxxxxxxxxxx
         sectionij=((2*NSpecies-smin)*(smin-1))/2+smax
         
         !IF (drmagij>RangePP(2*sectionij)) CYCLE !beyond cut-off
         
         dxinv=IncrPPinv(sectionij)
         mincutoff=RangePP(2*sectionij-1)
         tmp=(drmagij-mincutoff)*dxinv
         k0=FLOOR(tmp) !table index for reading density and pair potential, use this instead of abi
         
         IF (k0<2) THEN
            WRITE(6,*) "$Err>> SW index below PP array lower bound"
            WRITE(UNIT=6,FMT='("...too low drmagij:",ES14.5)') drmagij
            WRITE(6,*) "Index in interp table=", k0
            WRITE(UNIT=6,FMT='("Main atom index:",I4,"           Main atom coord:",3f10.3)') &
               iatom,AtomCoord(3*iatom-2:3*iatom)
            WRITE(6,*) "DRMAG="
            DO pos=r1,r2
               jatom=VLList(pos)
               WRITE(UNIT=6,FMT='("Atom index:",I4,"  Rad:",f10.3," Interac atom coord:",3f10.3)') &
                 jatom,VLdrmag(pos),AtomCoord(3*jatom-2:3*jatom)
            END DO
            errorstatus=1
            RETURN
         END IF
         
         p=MIN(tmp-REAL(k0,dp),1.0_dp) !remaining part

         k=MIN(k0,SizePP(sectionij)-1)+(sectionij-1)*MaxSWTableSize
         EnergySWTerm=PP(k)+ 0.5_dp*p*(PP(k+1)-PP(k-1)  +p*(PP(k+1)-PP(k)-PP(k)+PP(k-1)))
         !IF (.NOT. VLListDomainAtom(j)) EnergySWTerm=EnergySWTerm*0.5_dp !only half contribution
         EnergySW=EnergySW+EnergySWTerm
         !write(*,*) "pair:",EnergySWTerm,drmagij

         k=MIN(k0,SizePP(sectionij)-2)+(sectionij-1)*MaxSWTableSize
         gi=PP(k+1)-PP(k-1)
         gi1=PP(k+2)-PP(k)
         tgrad=0.5_dp*(gi+p*(gi1-gi))*dxinv !gradient of the pair potential between i and jatom

         AtomForce(3*iatom-2:3*iatom)=AtomForce(3*iatom-2:3*iatom)+(tgrad/drmagij)*drij
         !IF (VLListDomainAtom(j)) 
         AtomForce(3*jatom-2:3*jatom)=AtomForce(3*jatom-2:3*jatom)-(tgrad/drmagij)*drij

         
         !force calculation for triplet term iatom-jatom-katom
         epsij=epssqrt(sectionij)
         Lfnij=lam(sectionij) !this is (lam1*lam2)**0.25 in Eq
!WRITE(6,*) "epsij,lam:",epsij,Lfnij

         !get F3ij
         dxinv=IncrF3inv(sectionij)
         mincutoff=RangeF3(2*sectionij-1)
         tmp=(drmagij-mincutoff)*dxinv
         k0=FLOOR(tmp)
         
         IF (k0<2) THEN
            WRITE(6,*) "$Err>> SW index below F3 array lower bound"
            WRITE(UNIT=6,FMT='("...too low drmagij:",ES14.5)') drmagij
            WRITE(6,*) "Index in interp table=", k0
            WRITE(UNIT=6,FMT='("Main atom index:",I4,"           Main atom coord:",3f10.3)') &
               iatom,AtomCoord(3*iatom-2:3*iatom)
            WRITE(6,*) "DRMAG="
            DO pos=r1,r2
               jatom=VLList(pos)
               WRITE(UNIT=6,FMT='("Atom index:",I4,"  Rad:",f10.3," Interac atom coord:",3f10.3)') &
                  jatom,VLdrmag(pos),AtomCoord(3*jatom-2:3*jatom)
            END DO
            errorstatus=1
            RETURN
         END IF
         
         p=MIN(tmp-REAL(k0,dp),1.0_dp) !remaining part

         k=MIN(k0,SizeF3(sectionij)-2)+(sectionij-1)*MaxSWTableSize
            
         F3ij=F3(k)+ 0.5_dp*p*(F3(k+1)-F3(k-1)  +p*(F3(k+1)-F3(k)-F3(k)+F3(k-1)))
         
         k=MIN(k0,SizeF3(sectionij)-2)+(sectionij-1)*MaxSWTableSize
         gi=F3(k+1)-F3(k-1)
         gi1=F3(k+2)-F3(k)
         F3ijgrad=0.5_dp*(gi+p*(gi1-gi))*dxinv

            
         DO j1=j+1,r2 !triplet terms
         
            katom=VLList(j1)
            s3=AtomSpecies(katom)
            
            drik=VLdr(3*j1-2:3*j1)            ! gives vector between atoms i and k
            drmagik=VLdrmag(j1)      ! gives spacing between the same
            
            !xxxxxxxxxxxxxxxxxxxxxxxxxxSETTING UP F3 ENERGY AND FORCE TERMSxxxxxxxxxxxxxxxx
            smin=MIN(s1,s3)
            smax=MAX(s1,s3)
            sectionik=((2*NSpecies-smin)*(smin-1))/2+smax
            
            !IF (drmagik>RangePP(2*sectionik)) CYCLE !beyond cut-off
         
            !WRITE(6,*) "iatom,jatom,katom:",iatom,jatom,katom
            
            drjk=drik-drij
            drmagjk=SQRT(drjk(1)*drjk(1)+drjk(2)*drjk(2)+drjk(3)*drjk(3))
            
            !IF (drmagjk>RangePP(2*sectionjk)) CYCLE !beyond cut-off
            
            !WRITE(6,*) "drmags:",drmagij,drmagjk,drmagik
            
            costheta_ijk=-(drij(1)*drjk(1)+drij(2)*drjk(2)+drij(3)*drjk(3))/drmagij/drmagjk
            costheta_jik=(drij(1)*drik(1)+drij(2)*drik(2)+drij(3)*drik(3))/drmagij/drmagik
            costheta_ikj=(drik(1)*drjk(1)+drik(2)*drjk(2)+drik(3)*drjk(3))/drmagik/drmagjk
            phi_ijk=costheta_ijk+1._dp/3._dp
            phi_jik=costheta_jik+1._dp/3._dp
            phi_ikj=costheta_ikj+1._dp/3._dp
            
            epsik=epssqrt(sectionik)
            Lfnik=lam(sectionik)
            !get F3ik
            dxinv=IncrF3inv(sectionik)
            mincutoff=RangeF3(2*sectionik-1)
            tmp=(drmagik-mincutoff)*dxinv
            k0=FLOOR(tmp)
            IF (k0<2) THEN
               WRITE(6,*) "$Err>> SW index below F3 array lower bound"
               WRITE(UNIT=6,FMT='("...too low drmagij:",ES14.5)') drmagik
               WRITE(6,*) "Index in interp table=", k0
               errorstatus=1
               RETURN
            END IF
         
            p=MIN(tmp-REAL(k0,dp),1.0_dp) !remaining part

            k=MIN(k0,SizeF3(sectionik)-1)+(sectionik-1)*MaxSWTableSize
            
            F3ik=F3(k)+ 0.5_dp*p*(F3(k+1)-F3(k-1)  +p*(F3(k+1)-F3(k)-F3(k)+F3(k-1)))
            
            k=MIN(k0,SizeF3(sectionik)-2)+(sectionik-1)*MaxSWTableSize
            gi=F3(k+1)-F3(k-1)
            gi1=F3(k+2)-F3(k)
            F3ikgrad=0.5_dp*(gi+p*(gi1-gi))*dxinv
            
            !j-k terms
            smin=MIN(s2,s3)
            smax=MAX(s2,s3)
            sectionjk=((2*NSpecies-smin)*(smin-1))/2+smax
            epsjk=epssqrt(sectionjk)
            Lfnjk=lam(sectionjk)
            !get F3jk
            dxinv=IncrF3inv(sectionjk)
            mincutoff=RangeF3(2*sectionjk-1)
            tmp=(drmagjk-mincutoff)*dxinv
            k0=FLOOR(tmp)
            IF (k0<2) THEN
               WRITE(6,*) "$Err>> SW index below F3 array lower bound"
               WRITE(UNIT=6,FMT='("...too low drmagij:",ES14.5)') drmagjk
               WRITE(6,*) "Index in interp table=", k0
               errorstatus=1
               RETURN
            END IF
         
            p=MIN(tmp-REAL(k0,dp),1.0_dp) !remaining part

            k=MIN(k0,SizeF3(sectionjk)-1)+(sectionjk-1)*MaxSWTableSize
            
            F3jk=F3(k)+ 0.5_dp*p*(F3(k+1)-F3(k-1)  +p*(F3(k+1)-F3(k)-F3(k)+F3(k-1)))
            
            k=MIN(k0,SizeF3(sectionjk)-2)+(sectionjk-1)*MaxSWTableSize
            gi=F3(k+1)-F3(k-1)
            gi1=F3(k+2)-F3(k)
            F3jkgrad=0.5_dp*(gi+p*(gi1-gi))*dxinv
            
            !xxxxxxxxxxxxxxxxCOMPLETED SETTING UP F3 ENERGY AND FORCE TERMSxxxxxxxxxxxxxxxx
            
            !Energy and force calculation for triplet - jik
            EnergySW=EnergySW+epsij*epsik*Lfnij*Lfnik*F3ij*F3ik*phi_jik*phi_jik
            EnergySW=EnergySW+epsij*epsjk*Lfnij*Lfnjk*F3ij*F3jk*phi_ijk*phi_ijk
            EnergySW=EnergySW+epsik*epsjk*Lfnik*Lfnjk*F3ik*F3jk*phi_ikj*phi_ikj
         !write(*,*) "triplet:",epsij*epsik*Lfnij*Lfnik*F3ij*F3ik*phi_jik*phi_jik
         !write(*,*) "triplet:",epsij*epsjk*Lfnij*Lfnjk*F3ij*F3jk*phi_ijk*phi_ijk
         !write(*,*) "triplet:",epsik*epsjk*Lfnik*Lfnjk*F3ik*F3jk*phi_ikj*phi_ikj
            
            !Calculation of forces for many-body term
            !force acting on atom i (top experssion to be used for j and k as well)
            drij1=drij/drmagij
            drik1=drik/drmagik
            drjk1=drjk/drmagjk
            
            !force acting on atom i 
            dhdr= &
               !first term involving h(j,i,k) -- with atom i at center
               epsij*epsik*Lfnij*Lfnik* &
                 ( phi_jik*phi_jik*( -F3ijgrad*drij1*F3ik -F3ikgrad*drik1*F3ij )+ &
                   2._dp*F3ij*F3ik*phi_jik* &
                   (-(drij1/drmagik+drik1/drmagij) + costheta_jik*(drij1/drmagij+drik1/drmagik)) )+ &
               !second term involving h(i,j,k) -- with atom j at center
               epsij*epsjk*Lfnij*Lfnjk* &
                 ( -F3ijgrad*drij1*F3jk*phi_ijk*phi_ijk+ &
                   2._dp*F3ij*F3jk*phi_ijk*(drjk1/drmagij + costheta_ijk*drij1/drmagij) )+ &
               !third term involving h(i,k,j) -- with atom k at center
               epsik*epsjk*Lfnik*Lfnjk* &
                 ( -F3ikgrad*drik1*F3jk*phi_ikj*phi_ikj+ &
                   2._dp*F3ik*F3jk*phi_ikj*(-drjk1/drmagik + costheta_ikj*drik1/drmagik) )
            AtomForce(3*iatom-2:3*iatom)=AtomForce(3*iatom-2:3*iatom)-dhdr
            
            !force acting on atom j
            !IF (VLListDomainAtom(j)) 
            dhdr= &
               !first term involving h(i,j,k) -- with atom j at center
               epsij*epsjk*Lfnij*Lfnjk* &
                 ( (F3ijgrad*drij1*F3jk -F3ij*F3jkgrad*drjk1)*phi_ijk*phi_ijk+ &
                   2._dp*F3ij*F3jk*phi_ijk* &
                  ((drij1/drmagjk-drjk1/drmagij) + costheta_ijk*(-drij1/drmagij+drjk1/drmagjk)) ) + &
               !second term involving h(j,i,k) -- with atom i at center
               epsij*epsik*Lfnij*Lfnik* &
                 ( F3ijgrad*drij1*F3ik*phi_jik*phi_jik+ &
                   2._dp*F3ij*F3ik*phi_jik*(drik1/drmagij - costheta_jik*drij1/drmagij) )+ &
               !third term involving h(i,k,j) -- with atom k at center
               epsik*epsjk*Lfnik*Lfnjk* &
                 ( -F3jkgrad*F3ik*drjk1*phi_ikj*phi_ikj+ &
                   2._dp*F3ik*F3jk*phi_ikj*(-drik1/drmagjk + costheta_ikj*drjk1/drmagjk) )
            AtomForce(3*jatom-2:3*jatom)=AtomForce(3*jatom-2:3*jatom)-dhdr
            
            !force acting on atom k
            !IF (VLListDomainAtom(j1)) 
            dhdr= &
               !first term involving h(i,k,j) -- with atom k at center
               epsjk*epsik*Lfnjk*Lfnik* &
               ( (F3jkgrad*drjk1+F3ik*F3ikgrad*drik1*F3jk)*phi_ikj*phi_ikj+ &
                 2._dp*F3jk*F3ik*phi_ikj* &
               (drjk1/drmagik+drik1/drmagjk - costheta_ikj*(drjk1/drmagjk+drik1/drmagik)) )+ &
               !second term involving h(i,j,k) -- with atom j at center
               epsij*epsjk*Lfnij*Lfnjk* &
               ( F3ij*F3jkgrad*drjk1*phi_ijk*phi_ijk+ &
                 2._dp*F3ij*F3jk*phi_ijk*(-drij1/drmagjk - costheta_ijk*drjk1/drmagjk) )+ &
               !third term involving h(j,i,k) -- with atom i at center
               epsij*epsik*Lfnij*Lfnik* &
               ( F3ij*F3ikgrad*drik1*phi_jik*phi_jik+ &
               2._dp*F3ij*F3ik*phi_jik*(drij1/drmagik - costheta_jik*drik1/drmagik) )
            AtomForce(3*katom-2:3*katom)=AtomForce(3*katom-2:3*katom)-dhdr

            
         END DO
      END DO
      !WRITE(6,*) EnergySW
      !STOP
