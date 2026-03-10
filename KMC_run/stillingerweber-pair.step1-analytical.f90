!finds the v2-term
!         dxinv=IncrPPinv(sectionij)
!         mincutoff=RangePP(2*sectionij-1)
!         tmp=(drmagij-mincutoff)*dxinv
!         k0=FLOOR(tmp) !table index for reading density and pair potential, use this instead of abi
!         p=MIN(tmp-REAL(k0,dp),1.0_dp) !remaining part
         
!         IF (k0<2) THEN
!            WRITE(6,*) "$Err>> SW index below PP array lower bound"
!            WRITE(UNIT=6,FMT='("...too low drmagij:",ES14.5)') drmagij
!            WRITE(6,*) "Index in interp table=", k0
!            WRITE(UNIT=6,FMT='("Main atom index:",I4,"           Main atom coord:",3f10.3)') &
!              iatom,AtomCoord(3*iatom-2:3*iatom)
!           WRITE(6,*) "DRMAG="
!            DO pos=r1,r2
!               jatom=VLList(pos)
!               WRITE(UNIT=6,FMT='("Atom index:",I4,"  Rad:",f10.3," Interac atom coord:",3f10.3)') &
!                 jatom,VLdrmag(pos),AtomCoord(3*jatom-2:3*jatom)
!            END DO
!            errorstatus=1
!            RETURN
!            k=(sectionij-1)*MaxSWTableSize+1
!            EnergySW=EnergySW+ PP(k) + tmp*(PP(k+1)-PP(k))
!            tgrad=(PP(k+1)-PP(k))*dxinv
!            AtomForce(3*iatom-2:3*iatom)=AtomForce(3*iatom-2:3*iatom)+(tgrad/drmagij)*drij
!            AtomForce(3*jatom-2:3*jatom)=AtomForce(3*jatom-2:3*jatom)-(tgrad/drmagij)*drij
         
!         END IF
         
!            k=MIN(k0,SizePP(sectionij)-1)+(sectionij-1)*MaxSWTableSize
!            EnergySWTerm=PP(k)+ 0.5_dp*p*(PP(k+1)-PP(k-1)  +p*(PP(k+1)-PP(k)-PP(k)+PP(k-1)))
            !IF (.NOT. VLListDomainAtom(j)) EnergySWTerm=EnergySWTerm*0.5_dp !only half contribution
            
            
!            write(*,*) rij
            IF (rij<ac(sectionij)) THEN
               rij2=rij*rij
               rij4inv=1._dp/rij2/rij2
               exporij=exp(1._dp/(rij-ac(sectionij)))
         
               EnergySWTerm=epsln(sectionij)*A(sectionij)*(B(sectionij)*rij4inv-1._dp)*exporij
               EnergySW=EnergySW+EnergySWTerm
               
               
               !write(*,*) "pair:",EnergySWTerm,drmagij
            
            

!            k=MIN(k0,SizePP(sectionij)-2)+(sectionij-1)*MaxSWTableSize
!            gi=PP(k+1)-PP(k-1)
!            gi1=PP(k+2)-PP(k)
!            tgrad=0.5_dp*(gi+p*(gi1-gi))*dxinv !gradient of the pair potential between i and jatom

!            AtomForce(3*iatom-2:3*iatom)=AtomForce(3*iatom-2:3*iatom)+(tgrad/drmagij)*drij
!            !IF (VLListDomainAtom(j)) 
!            AtomForce(3*jatom-2:3*jatom)=AtomForce(3*jatom-2:3*jatom)-(tgrad/drmagij)*drij
!         END IF
            
               forceoniatom=(epsln(sectionij)*A(sectionij)*exporij)*(-4._dp*B(sectionij)*rij4inv/rij) &
                                 - (EnergySWTerm)/(rij-ac(sectionij))/(rij-ac(sectionij))
               forceoniatom=forceoniatom/sigma(sectionij)/drmagij
               AtomForce(3*iatom-2:3*iatom)= AtomForce(3*iatom-2:3*iatom)+forceoniatom*drij
            
               AtomForce(3*jatom-2:3*jatom)= AtomForce(3*jatom-2:3*jatom)-forceoniatom*drij
            END IF























