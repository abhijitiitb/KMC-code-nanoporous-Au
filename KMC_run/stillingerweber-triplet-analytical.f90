!finds the triplet contribution from iatom-jatom-katom
            drik=VLdr(3*j1-2:3*j1)            ! gives vector between atoms i and k
            drmagik=VLdrmag(j1)      ! gives spacing between the same
            
            !xxxxxxxxxxxxxxxxxxxxxxxxxxSETTING UP F3 ENERGY AND FORCE TERMSxxxxxxxxxxxxxxxx
            smin=MIN(s1,s3)
            smax=MAX(s1,s3)
            sectionik=((2*NSpecies-smin)*(smin-1))/2+smax
            
            !IF (drmagik>RangePP(2*sectionik)) CYCLE !beyond cut-off
         
            !WRITE(6,*) "iatom,jatom,katom:",iatom,jatom,katom
            
            costheta_jik=(drij(1)*drik(1)+drij(2)*drik(2)+drij(3)*drik(3))/drmagij/drmagik
            phi_jik=costheta_jik+1._dp/3._dp
            
!            epsik=epssqrt(sectionik)
!            Lfnik=lam(sectionik)
            !get F3ik
!            dxinv=IncrF3inv(sectionik)
!            mincutoff=RangeF3(2*sectionik-1)
!            tmp=(drmagik-mincutoff)*dxinv
!            k0=FLOOR(tmp)
!            p=MIN(tmp-REAL(k0,dp),1.0_dp) !remaining part
            
!            IF (k0<2) THEN
!               WRITE(6,*) "$Err>> SW index below F3 array lower bound"
!               WRITE(UNIT=6,FMT='("...too low drmagij:",ES14.5)') drmagik
!               WRITE(6,*) "Index in interp table=", k0
!               errorstatus=1
!               RETURN
               
!               k=(sectionik-1)*MaxSWTableSize+1
!               F3ik=F3(k) + tmp*(F3(k+1)-F3(k))
!               F3ikgrad=(F3(k+1)-F3(k))*dxinv
!            END IF
         

!               k=MIN(k0,SizeF3(sectionik)-1)+(sectionik-1)*MaxSWTableSize
            
!               F3ik=F3(k)+ 0.5_dp*p*(F3(k+1)-F3(k-1)  +p*(F3(k+1)-F3(k)-F3(k)+F3(k-1)))
            
!               k=MIN(k0,SizeF3(sectionik)-2)+(sectionik-1)*MaxSWTableSize
!               gi=F3(k+1)-F3(k-1)
!               gi1=F3(k+2)-F3(k)
!               F3ikgrad=0.5_dp*(gi+p*(gi1-gi))*dxinv
!            END IF
            
            !xxxxxxxxxxxxxxxxCOMPLETED SETTING UP F3 ENERGY AND FORCE TERMSxxxxxxxxxxxxxxxx
            !drjk(1)=drik(1)-drij(1)
            !drjk(2)=drik(2)-drij(2)
            !drjk(3)=drik(3)-drij(3)
            !drmagjk=drjk(1)*drjk(1)+drjk(2)*drjk(2)+drjk(3)*drjk(3)
            
            !costheta_ijk=-(drij(1)*drjk(1)+drij(2)*drjk(2)+drij(3)*drjk(3))/drmagij/drmagjk
            !phi_ijk=costheta_ijk+1._dp/3._dp
            
            !costheta_ikj=(drik(1)*drjk(1)+drik(2)*drjk(2)+drik(3)*drjk(3))/drmagij/drmagik
            !phi_ikj=costheta_ikj+1._dp/3._dp
            
            !rjk=drmagjk/sigma
            rik=drmagik/sigma(sectionik)
            expogmrij=exp(gamma(sectionij)/(rij-ac(sectionij)))
            expogmrik=exp(gamma(sectionik)/(rik-ac(sectionik)))
            
            epslnsqrt=sqrt(epsln(sectionij)*epsln(sectionik))
            !expogmrjk=exp(gamma/(rjk-ac))
            
            
            
!            IF (rik<ac .and. rij<ac .and. rjk<ac) THEN
            !Energy and force calculation for triplet - jik
            h_jik=0._dp
            !h_ijk=0._dp
            !h_ikj=0._dp
            
            IF (rij<ac(sectionij) .and. rik<ac(sectionik)) THEN
            
            h_jik=epslnsqrt*lmd(sectionij)*lmd(sectionik)*expogmrij*expogmrik*phi_jik*phi_jik
            
            !IF (rjk<ac .and. rij<ac) h_ijk=epsln*lmd*expogmrij*expogmrjk*phi_ijk*phi_ijk
            
            !IF (rik<ac .and. rjk<ac) h_ikj=epsln*lmd*expogmrik*expogmrjk*phi_ikj*phi_ikj
            
            EnergySWv3=h_jik !+h_ijk+h_ikj
            
            
            
!            write(*,*) phi_jik
!            stop
!            ELSE
!            EnergySWv3=0.0_dp
!            END IF
            
            EnergySW=EnergySW+EnergySWv3       !epsij*epsik*Lfnij*Lfnik*F3ij*F3ik*phi_jik*phi_jik
         !   IF (iatom==4) THEN
         !      write(*,*) rij
         !      write(*,*) drij
         !      write(*,*) rik
         !      write(*,*) drik
         !      write(*,*) costheta_jik
         !      write(*,*) EnergySW
            END IF                                                     
            
            
            !Calculation of forces for many-body term
            !force acting on atom i (top experssion to be used for j and k as well)
            IF (rik<ac(sectionik) .and. rij<ac(sectionij)) THEN
               drij1=drij/drmagij
               drik1=drik/drmagik
            
            !force acting on atom i 
!            dhdr= &
!               !term involving h(j,i,k) -- with atom i at center
!               epsij*epsik*Lfnij*Lfnik* &
!                 ( phi_jik*phi_jik*( -F3ijgrad*drij1*F3ik -F3ikgrad*drik1*F3ij )+ &
!                   2._dp*F3ij*F3ik*phi_jik* &
!                   (-(drij1/drmagik+drik1/drmagij) + costheta_jik*(drij1/drmagij+drik1/drmagik)) )
               F3ijgrad=(-gamma(sectionij)*expogmrij)/sigma(sectionij)/(rij-ac(sectionij))**2
               F3ikgrad=(-gamma(sectionik)*expogmrik)/sigma(sectionik)/(rik-ac(sectionik))**2 
                  
               dhdr=epslnsqrt*lmd(sectionij)*lmd(sectionik)* &
                 ( phi_jik*phi_jik*( -F3ijgrad*drij1*expogmrik -F3ikgrad*drik1*expogmrij )+ &
                   2._dp*expogmrij*expogmrik*phi_jik* &
                   (-(drij1/drmagik+drik1/drmagij) + costheta_jik*(drij1/drmagij+drik1/drmagik)) )     
               AtomForce(3*iatom-2:3*iatom)=AtomForce(3*iatom-2:3*iatom)-dhdr
            
            !force acting on atom j
            !IF (VLListDomainAtom(j)) 
!            dhdr= &
!               !term involving h(j,i,k) -- with atom i at center
!               epsij*epsik*Lfnij*Lfnik* &
!                 ( F3ijgrad*drij1*F3ik*phi_jik*phi_jik+ &
!                   2._dp*F3ij*F3ik*phi_jik*(drik1/drmagij - costheta_jik*drij1/drmagij) )
                   
               dhdr=epslnsqrt*lmd(sectionij)*lmd(sectionik)* &
                 ( F3ijgrad*drij1*expogmrik*phi_jik*phi_jik+ &
                   2._dp*expogmrij*expogmrik*phi_jik*(drik1/drmagij - costheta_jik*drij1/drmagij) )
               AtomForce(3*jatom-2:3*jatom)=AtomForce(3*jatom-2:3*jatom)-dhdr
            
            !force acting on atom k
            !IF (VLListDomainAtom(j1)) 
!            dhdr= &
!               !term involving h(j,i,k) -- with atom i at center
!               epsij*epsik*Lfnij*Lfnik* &
!               ( F3ij*F3ikgrad*drik1*phi_jik*phi_jik+ &
!              2._dp*F3ij*F3ik*phi_jik*(drij1/drmagik - costheta_jik*drik1/drmagik) )
               
               dhdr=epslnsqrt*lmd(sectionij)*lmd(sectionik)* &
               ( expogmrij*F3ikgrad*drik1*phi_jik*phi_jik+ &
               2._dp*expogmrij*expogmrik*phi_jik*(drij1/drmagik - costheta_jik*drik1/drmagik) )
               AtomForce(3*katom-2:3*katom)=AtomForce(3*katom-2:3*katom)-dhdr
            END IF



