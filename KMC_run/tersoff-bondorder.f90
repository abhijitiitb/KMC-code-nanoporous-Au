!find atoms katom that can form a triplet
			   zetaij=0._dp
            zetaji=0._dp
            dzetaij_i=0._dp
            dzetaij_j=0._dp
            dzetaji_i=0._dp
            dzetaji_j=0._dp
            
            DO j1=r1,r2 !neighbors of i for calculation of bij and bji
         
               katom=VLList(j1)
               s3=AtomSpecies(katom)
            
               !bji terms with theta_jik obtained
               drik=VLdr(j1)
               drmagik=VLdrmag(j1)
               drjk=drik-drij
               drmagjk=SQRT(drjk(1)*drjk(1)+drjk(2)*drjk(2)+drjk(3)*drjk(3))
               
               write(*,*) "reached1"
               write(*,*) drmagik,drmagjk
               !!derivative of r's wrt r_ixi annd r_jxi
               drij_rjxi = drij/drmagij
               drij_rixi = -1*drij_rjxi
               drik_rixi = -1*drik/drmagik
               drjk_rjxi = -1*drjk/drmagjk
               
               !bij terms with theta_ijk obtained
               write(*,*) "Rangefcij"
               write(*,*) Rangefcij
               IF (drmagik<Rangefcij) THEN !add contribution
                  write(*,*) "reached2"
                  drtik=PIO2*(MAX(drmagik,Rij-RDij)-Rij)/RDij
                  fcutik=0.5_dp-0.5_dp*SIN(drtik)
                  costheta=(drik(1)*drij(1)+drik(2)*drij(2)+drik(3)*drij(3))/drmagik/drmagij
                  gtheta=1+cij*cij/dij/dij-cij*cij/(dij*dij+(hij-costheta)*(hij-costheta))
                  lterm=lam3ij*(drmagij-drmagik)
                  expr=EXP(lterm*lterm*lterm)
                  zetaji=zetaji+fcutik*gtheta*expr
                  write(*,*) "b_ij"
                  write(*,*) fcutik,costheta,gtheta,lterm,zetaij
                  
                  !derivative of zetaij wrt rixi and rjxi
                  dfcutik=-0.5_dp*PIO2*COS(drtik)/RDij !dfc_rik/drik
                  dcostheta_i=-costheta*(drij_rixi/drmagij + drik_rixi/drmagik) + drij_rixi/drmagik + drik_rixi/drmagij!array
                  dcostheta_j=-costheta*(drij_rjxi/drmagij) - drik_rixi/drmagij!array
                  gterm=2*cij*cij*(hij-costheta)/((dij*dij + ((hij-costheta)**2))**2) !real
                  dg_i=gterm*dcostheta_i!array
                  dg_j=gterm*dcostheta_j!array
                  zterm=expr*fcutik*gtheta*3*lam3ij*lterm*lterm!real
                  dzetaij_i=dzetaij_i + expr*(dfcutik*drik_rixi*gtheta + fcutik*dg_i) + zterm*(drij_rixi-drik_rixi)!array
                  dzetaij_j=dzetaij_j + expr*fcutik*dg_j + zterm*drij_rjxi!array
               END IF
               
               !bji terms with theta_jik obtained
               IF (drmagjk<Rangefcij) THEN !add contribution
                  write(*,*) "reached3"
				      drtjk=PIO2*(MAX(drmagjk,Rij-RDij)-Rij)/RDij
				      fcutjk=0.5_dp-0.5_dp*SIN(drtjk)
                  costheta=-(drjk(1)*drij(1)+drjk(2)*drij(2)+drjk(3)*drij(3))/drmagjk/drmagij
                  gtheta=1+cij*cij/dij/dij-cij*cij/(dij*dij+(hij-costheta)*(hij-costheta))
                  lterm=lam3ij*(drmagij-drmagjk)
                  expr=EXP(lterm*lterm*lterm)
                  zetaij=zetaij+fcutjk*gtheta*expr
                  write(*,*) "b_ji"
                  write(*,*) fcutjk,costheta,gtheta,lterm,zetaij
                  
                  !derivative of zetaji wrt rixi and rjxi
                  dfcutjk=-0.5_dp*PIO2*COS(drtjk)/RDij !dfc_rjk/drjk
                  dcostheta_i=-costheta*(drij_rixi/drmagij) - drjk_rjxi/drmagij
                  dcostheta_j=-costheta*(drij_rjxi/drmagij + drjk_rjxi/drmagjk) - drij_rixi/drmagjk + drjk_rjxi/drmagij
                  gterm=2*cij*cij*(hij-costheta)/((dij*dij + ((hij-costheta)**2))**2) 
                  dg_i=gterm*dcostheta_i
                  dg_j=gterm*dcostheta_j
                  zterm=expr*fcutjk*gtheta*3*lam3ij*lterm*lterm
                  dzetaji_i=dzetaji_i+expr*fcutjk*dg_i + zterm*drij_rixi
                  dzetaji_j=dzetaji_j + expr*(dfcutjk*drjk_rjxi*gtheta + fcutjk*dg_j) + zterm*(drij_rjxi-drjk_rjxi)
		  
               END IF
            END DO
         
            raj=(iatom-1)*MaxAtomPerAtom !start index minus 1
            r0j=VLListRange(jatom) !range for ith atom
            r1j=ra+1
            r2j=ra+r0
            DO j1=r1j,r2j !neighbors of j for calculation of bij and bji
         
               katom=VLList(j1)
               s3=AtomSpecies(katom)
            
               drjk=VLdr(j1)
               drmagjk=VLdrmag(j1)
               drik=drjk+drij
               drmagik=SQRT(drik(1)*drik(1)+drik(2)*drik(2)+drik(3)*drik(3))
               
               !IF (drmagik> .AND. drmagjk<Rangefcij) THEN !if drmagik< then katom would already have been a neighbor of iatom and considered earlier
                  
               !END IF
            END DO
            
            bzn_ij = (betaij*zetaij)**nij
            bzn_ji = (betaij*zetaji)**nij
            b_ij=(1._dp+bzn_ij)**(-0.5_dp/nij)!
            b_ji=(1._dp+bzn_ji)**(-0.5_dp/nij)!
            
            bterm_ij = -0.5_dp*((1._dp + bzn_ij)**(-1._dp-0.5_dp*nij))*bzn_ij/zetaij
            bterm_ji = -0.5_dp*((1._dp + bzn_ji)**(-1._dp-0.5_dp*nij))*bzn_ji/zetaji
            
            db_ij_drixi = bterm_ij*dzetaij_i
            db_ij_drjxi = bterm_ij*dzetaij_j
            db_ji_drixi = bterm_ji*dzetaji_i
            db_ji_drjxi = bterm_ji*dzetaji_j
            
