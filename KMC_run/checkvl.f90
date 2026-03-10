
         ri0=VL%ListRange(i)
         ri1=VL%MaxAtomPerAtom*(i-1)+1
         ri2=ri1-1+ri0
         
         !without using Verlet list
         DO j=i+1,NAtoms !j is greater than i
         
            !IF (i==j) CYCLE
            !dr=PBCdistance(AL,j,i)
            dr=AL%AtomCoord(3*j-2:3*j)-AL%AtomCoord(3*i-2:3*i)
            dr=dr-BoxSize*NINT(dr/BoxSize)
            drmag=SQRT(DOT_PRODUCT(dr,dr))
            
            IF (ANY(VL%List(ri1:ri2)==j) .AND. drmag>CutOff) THEN
               WRITE(6,*) "Err>> Atom should not be added to VL, still it is in the list"
               STOP
            END IF
            
            IF (drmag<=CutOff) THEN
         
               !first check the list of atom i
               Found=ANY(VL%List(ri1:ri2)==j)
               IF (.NOT. Found) THEN
                  WRITE(UNIT=6,FMT='("$ChkErr>> a) Unable to find atom index",I4)') j
                  WRITE(UNIT=6,FMT='(" ... in VL of atom index",I4)') i
                  WRITE(UNIT=6,FMT='(" >> distance:",ES15.8)') drmag
                  WRITE(UNIT=6,FMT='(" >> cutoff distance:",ES15.8)') CutOff
                  WRITE(UNIT=6,FMT='(" >> Coord  atom:",3ES15.8)') &
                     AL%AtomCoord(3*j-2:3*j)
                  WRITE(UNIT=6,FMT='(" >> Cell :",i5)') LL%AtomCellIndx(j)
                  r1=ri1
                  r2=ri2
                  CALL PrintList(i)
                  IF (drmag<CutOff-VL%Buffer) STOP
                  CheckVL=.FALSE.
                  STOP
               END IF
               counter(i)=counter(i)+1
               
               !now lets check if the distance is correct
               Found=.FALSE.
               DO k=ri1,ri2
                  IF (VL%List(k)==j) THEN
                     Found=.TRUE.
                     EXIT
                  END IF
               END DO
                      
               IF (ABS(drmag-VL%drmag(k))>1.e-8_dp &
                  .OR. ANY(ABS(dr-VL%dr(3*k-2:3*k))>1.e-8_dp)) THEN
                  WRITE(6,*) "$ChkErr>> a) Incorrect position or"// &
                     " distance stored in Verlet list"
                  WRITE(UNIT=6,FMT='("Problem with atom index",I4)') j
                  WRITE(UNIT=6,FMT='(" ... in VL of atom index",I4)') i
                  WRITE(UNIT=6,FMT='(" ... found drmag:",ES10.3)') drmag
                  WRITE(UNIT=6,FMT='(" ... stored drmag:",ES10.3)') VL%drmag(k)
                  WRITE(UNIT=6,FMT='(" ... coord1",3ES10.3)') dr
                  WRITE(UNIT=6,FMT='(" ... coord2",3ES10.3)') VL%dr(3*k-2:3*k)
                  !r1=ri1
                  !r2=ri2
                  CALL PrintList(i)
                  CheckVL=.FALSE.
                  !IF (drmag<CutOff-VL%Buffer) STOP
                  IF (VL%drmag(k)>0._dp .AND. ABS(drmag-VL%drmag(k))>0.1_dp) STOP
               END IF
               
               !Now repeat this check for the list of atom j and find i
               IF (.NOT. halflist) THEN
                  rj0=VL%ListRange(j)
                  rj1=VL%MaxAtomPerAtom*(j-1)+1
                  rj2=rj1-1+rj0
                  !this atom should be in the list
                  Found=ANY(VL%List(rj1:rj2)==i)
                  IF (.NOT. Found) THEN
                     WRITE(UNIT=6,FMT='("$ChkErr>> b) Unable to find atom index",I4)') i
                     WRITE(UNIT=6,FMT='(" ... in VL of atom index",I4)') j
                     WRITE(UNIT=6,FMT='(" >> distance:",ES10.3)') drmag
                     WRITE(UNIT=6,FMT='(" >> cutoff distance:",ES15.8)') CutOff
                     WRITE(UNIT=6,FMT='(" >> Coord  atom:",3ES10.3)') &
                       AL%AtomCoord(3*i-2:3*i)
                     WRITE(UNIT=6,FMT='(" >> Cell :",i5)') LL%AtomCellIndx(i)
                     !r1=rj1
                     !r2=rj2
                     CALL PrintList(j)
                     IF (drmag<CutOff-VL%Buffer) STOP
                     CheckVL=.FALSE.
                  END IF
                  counter(j)=counter(j)+1
               
                  !now lets check if the distance is correct
                  Found=.FALSE.
                  DO k=rj1,rj2
                     IF (VL%List(k)==i) THEN
                        Found=.TRUE.
                        EXIT
                     END IF
                  END DO
                  
                  IF (ABS(drmag-VL%drmag(k))>1.e-8_dp &
                   .OR. ANY(ABS(-dr-VL%dr(3*k-2:3*k))>1.e-8_dp)) THEN
                     WRITE(*,*) "$ChkErr>> b) Incorrect position or"// &
                       " distance stored in Verlet list"
                     WRITE(UNIT=6,FMT='("Problem with atom index",I4)') i
                     WRITE(UNIT=6,FMT='(" ... in VL of atom index",I4)') j
                     WRITE(UNIT=6,FMT='(" ... found drmag:",ES10.3)') drmag
                     WRITE(UNIT=6,FMT='(" ... stored drmag:",ES10.3)') VL%drmag(k)
                     WRITE(UNIT=6,FMT='(" ... coord1",3ES10.3)') dr
                     WRITE(UNIT=6,FMT='(" ... coord2",3ES10.3)') VL%dr(3*k-2:3*k)
                     !r1=rj1
                     !r2=rj2
                     CALL PrintList(j)
                     CheckVL=.FALSE.
                     !IF (drmag<CutOff-VL%Buffer) STOP
                     IF (VL%drmag(k)>0._dp .AND. ABS(drmag-VL%drmag(k))>0.1_dp) STOP
                  END IF
                  
               END IF
               
            END IF
            
         END DO
         
         
         IF (HalfList .AND. Counter(i)/=VL%ListRange(i)) THEN
            WRITE(6,*) "Verlet list has incorrect number of atoms"
            WRITE(6,*) "Atoms:",i
            WRITE(6,*) "Found # atoms:",Counter(i)
            WRITE(6,*) "Supposed # atoms:",VL%ListRange(i)
            CALL PrintList(i)
            stop
         END IF
