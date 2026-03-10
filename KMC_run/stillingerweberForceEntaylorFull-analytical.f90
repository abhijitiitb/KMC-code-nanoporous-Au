!Computes the forces using interpolation with full list
      
      ra=(iatom-1)*MaxAtomPerAtom !start index minus 1
      r0=VLListRange(iatom) !range for ith atom
      
      r1=ra+1
      r2=ra+r0

      s1=AtomSpecies(iatom)
      
         !IF (iatom==700) write(6,*)"Si ",AL%AtomCoord(3*iatom-2:3*iatom)
      !first loop to calculate force for pair porential term
      !WRITE(6,*) "iatom:",iatom,s1
!WRITE(6,*) "Atom position i:",AL%AtomCoord(3*iatom-2:3*iatom)
      
      DO j=r1,r2 !this will take care of the pair and many-body terms in one go
         jatom=VLList(j)
         !IF (iatom==700) write(6,*)"Si ",AL%AtomCoord(3*jatom-2:3*jatom)
         s2=AtomSpecies(jatom)
!WRITE(6,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
!WRITE(6,*) "jatom:",jatom,s2
         
         smin=MIN(s1,s2)
         smax=MAX(s1,s2)
         sectionij=((2*NSpecies-smin)*(smin-1))/2+smax
         drij=VLdr(3*j-2:3*j)  !vector between atoms i and j
         drmagij=VLdrmag(j)   ! spacing between atoms i and j
         rij=drmagij/sigma(sectionij)
!         write(*,*) iatom  
!         write(*,*) jatom 
!         write(*,*) drij 
!         write(*,*) drmagij
!         stop
      
!WRITE(6,*) "drij:",drij,drmagij
!WRITE(6,*) "Atom position j:",AL%AtomCoord(3*jatom-2:3*jatom)
         
         !xxxxxxxxxxxxxxxxxxxx  PAIR POTENTIAL TERM  xxxxxxxxxxxxxxxxxxx
         
         !IF (drmagij>RangePP(2*sectionij)) CYCLE !beyond cut-off
         
         IF (iatom<jatom) THEN !we require iatom<jatom

            INCLUDE "stillingerweber-pair.step1-analytical.f90"
           ! INCLUDE "stillingerweber-pair.step1array.f90"
         
         END IF
         
         !force calculation for triplet term iatom-jatom-katom
        ! INCLUDE "stillingerweber-pair.step2.f90"
        ! INCLUDE "stillingerweber-pair.step2array.f90"
         
         DO j1=r1,r2 !triplet term involving hjik
         
            katom=VLList(j1)
            s3=AtomSpecies(katom)
            
            IF (jatom<katom) THEN
          !  IF (jatom < katom) THEN
            
               INCLUDE "stillingerweber-triplet-analytical.f90"          
              ! INCLUDE "stillingerweber-tripletarray.f90"
            
            END IF
            
         END DO
         
      END DO
      !WRITE(6,*) EnergySW
      !STOP
