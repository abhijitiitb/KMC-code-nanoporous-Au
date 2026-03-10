   
   SUBROUTINE VLAddFullAAShortCutOff()
   !full list, same cells A and A (maincellindx==neighcellindx)
   !CellOutsideBox is FALSE(by defn since cutoff is short), PotentialCutOff is short compared to BoxSize
   
      IMPLICIT NONE
      INTEGER :: i,j,istart,iend
      INTEGER :: position,atomindx1,atomindx2
      INTEGER :: atomnumberneigh
      REAL(dp) :: dr(3),dr2mag,drmag
      
      istart=(maincellindx-1)*MaxAtomPerCell+1
      iend=istart+LLListRange(maincellindx)-1
      
      
      DO i=istart,iend-1
         
         atomindx1=LLList(i) !atom index from first cell
         
         DO j=i+1,iend

            atomindx2=LLList(j) !atom index from other cell

            dr=AtomCoord(3*atomindx2-2:3*atomindx2)-AtomCoord(3*atomindx1-2:3*atomindx1)
            dr=dr-BoxSize*NINT(dr/BoxSize)
            dr2mag=DOT_PRODUCT(dr,dr)
            
            IF (dr2mag<=CutOff2) THEN

!write(*,*) "Full list, short interaction AA:",atomindx1,atomindx2
               drmag=SQRT(dr2mag)
               
               atomnumberneigh=VLListRange(atomindx1)+1

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindx1-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindx2
               VLdr(3*position-2:3*position)=dr
               VLdrmag(position)=drmag
               VLListRange(atomindx1)=atomnumberneigh 
               
               
               atomnumberneigh=VLListRange(atomindx2)+1 !VL%ListRange()=0 to begin with

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindx2-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindx1
               VLdr(3*position-2:3*position)=-dr
               VLdrmag(position)=drmag
               VLListRange(atomindx2)=atomnumberneigh 
               
            END IF

         END DO

      END DO

   END SUBROUTINE VLAddFullAAShortCutOff
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE VLAddFullABShortCutOff()
   !full list, two different cells A and B (maincellindx/=neighcellindx)
   !CellOutsideBox is FALSE (by defn since CutOff is short), PotentialCutOff is short compared to box size

      IMPLICIT NONE
      INTEGER :: i,j,istart,iend,jstart,jend
      INTEGER :: position,atomindx1,atomindx2
      INTEGER :: atomnumberneigh
      REAL(dp) :: dr(3),dr2mag,drmag
      
      istart=(maincellindx-1)*MaxAtomPerCell+1
      iend=istart+LLListRange(maincellindx)-1

      jstart=(neighcellindx-1)*MaxAtomPerCell+1
      jend=jstart+LLListRange(neighcellindx)-1
      
      
      DO i=istart,iend
         
         atomindx1=LLList(i)
         
         DO j=jstart,jend

            atomindx2=LLList(j)

            dr=AtomCoord(3*atomindx2-2:3*atomindx2)-AtomCoord(3*atomindx1-2:3*atomindx1)
            dr=dr-BoxSize*NINT(dr/BoxSize)
            dr2mag=DOT_PRODUCT(dr,dr)
            
            IF (dr2mag<=CutOff2) THEN

!write(*,*) "Full list, short interaction AB:",atomindx1,atomindx2
               drmag=SQRT(dr2mag)
               
               atomnumberneigh=VLListRange(atomindx1)+1 !VL%ListRange()=0 to begin with

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindx1-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindx2
               VLdr(3*position-2:3*position)=dr
               VLdrmag(position)=drmag
               VLListRange(atomindx1)=atomnumberneigh
               

               atomnumberneigh=VLListRange(atomindx2)+1

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(*,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindx2-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindx1
               VLdr(3*position-2:3*position)=-dr
               VLdrmag(position)=drmag
               VLListRange(atomindx2)=atomnumberneigh

            END IF

         END DO

      END DO

   END SUBROUTINE VLAddFullABShortCutOff
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE VLAddHalfAAShortCutOff()
   !half list, same cells A and A (maincellindx==neighcellindx)
   !CellOutsideBox is FALSE (by defn since cut off is short), Potential cutoff is short compared to box size
      IMPLICIT NONE
      INTEGER :: i,j,istart,iend
      INTEGER :: position,atomindx1,atomindx2,atomindxmin,atomindxmax
      INTEGER :: atomnumberneigh
      REAL(dp) :: dr(3),dr2mag,drmag
      
      istart=(maincellindx-1)*MaxAtomPerCell+1
      iend=istart+LLListRange(maincellindx)-1
      
      DO i=istart,iend-1
         
         atomindx1=LLList(i) !atom index from first cell
         
         DO j=i+1,iend

            atomindx2=LLList(j) !atom index from other cell
            
            atomindxmin=MIN(atomindx1,atomindx2)
            atomindxmax=MAX(atomindx1,atomindx2)

            dr=AtomCoord(3*atomindxmax-2:3*atomindxmax)-AtomCoord(3*atomindxmin-2:3*atomindxmin)
            dr=dr-BoxSize*NINT(dr/BoxSize)
            dr2mag=DOT_PRODUCT(dr,dr)
            
            IF (dr2mag<=CutOff2) THEN
!write(*,*) "half-short-AA",atomindx1,atomindx2
               drmag=SQRT(dr2mag)

               atomnumberneigh=VLListRange(atomindxmin)+1

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindxmin-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindxmax
               VLdr(3*position-2:3*position)=dr
               VLdrmag(position)=drmag
               VLListRange(atomindxmin)=atomnumberneigh 
               
            END IF

         END DO

      END DO
   
   END SUBROUTINE VLAddHalfAAShortCutOff
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE VLAddHalfABShortCutOff()
   !half list, two different cells A and B (maincellindx/=neighcellindx)
   !CellOutsideBox is FALSE (by defn since cutoff is short), Potential cutoff is short compared to box size

      IMPLICIT NONE
      INTEGER :: i,j,istart,iend,jstart,jend
      INTEGER :: position,atomindx1,atomindx2,atomindxmin,atomindxmax
      INTEGER :: atomnumberneigh
      REAL(dp) :: dr(3),dr2mag,drmag
      
      istart=(maincellindx-1)*MaxAtomPerCell+1
      iend=istart+LLListRange(maincellindx)-1

      jstart=(neighcellindx-1)*MaxAtomPerCell+1
      jend=jstart+LLListRange(neighcellindx)-1
      
      
      DO i=istart,iend
         
         atomindx1=LLList(i)
         
         DO j=jstart,jend

            atomindx2=LLList(j)
            
            atomindxmin=MIN(atomindx1,atomindx2)
            atomindxmax=MAX(atomindx1,atomindx2)

            dr=AtomCoord(3*atomindxmax-2:3*atomindxmax)-AtomCoord(3*atomindxmin-2:3*atomindxmin)
            dr=dr-BoxSize*NINT(dr/BoxSize)
            dr2mag=DOT_PRODUCT(dr,dr)
            
            IF (dr2mag<=CutOff2) THEN
!write(*,*) "half-short-AB",atomindx1,atomindx2
               drmag=SQRT(dr2mag)

               atomnumberneigh=VLListRange(atomindxmin)+1

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindxmin-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindxmax
               VLdr(3*position-2:3*position)=dr
               VLdrmag(position)=drmag
               VLListRange(atomindxmin)=atomnumberneigh 

            END IF

         END DO

      END DO

   END SUBROUTINE VLAddHalfABShortCutOff
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      
   SUBROUTINE VLAddFullAACellInBoxLongCutOff()
   !full list, same cells A and A (maincellindx==neighcellindx)
   !CellOutsideBox is FALSE, potential cutoff is long in comparison to box size
   
      IMPLICIT NONE
      INTEGER :: i,j,istart,iend
      INTEGER :: position,atomindx1,atomindx2
      INTEGER :: atomnumberneigh
      REAL(dp) :: dr(3),dr2mag,drmag
      
      istart=(maincellindx-1)*MaxAtomPerCell+1
      iend=istart+LLListRange(maincellindx)-1
      
      
      DO i=istart,iend-1
         
         atomindx1=LLList(i)
         
         DO j=i+1,iend

            atomindx2=LLList(j) !atom index from other cell

            dr=AtomCoord(3*atomindx2-2:3*atomindx2)-AtomCoord(3*atomindx1-2:3*atomindx1)
            dr=dr-BoxSize*(NINT(dr/BoxSize)*AllowPBC)
            !dr=dr+REAL(BoxPosition)*BoxSize  !since BoxPosition=0
            dr2mag=DOT_PRODUCT(dr,dr)
            
            IF (dr2mag<=CutOff2) THEN

               drmag=SQRT(dr2mag)
               
               atomnumberneigh=VLListRange(atomindx1)+1

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindx1-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindx2
               VLdr(3*position-2:3*position)=dr
               VLdrmag(position)=drmag
               VLListRange(atomindx1)=atomnumberneigh
               !IF (CellOutsideBox) VLListDomainAtom(position)=.FALSE. !atomindx2 does not belong to the domain
               
               
               atomnumberneigh=VLListRange(atomindx2)+1

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindx2-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindx1
               VLdr(3*position-2:3*position)=-dr
               VLdrmag(position)=drmag
               VLListRange(atomindx2)=atomnumberneigh
               !IF (CellOutsideBox) VLListDomainAtom(position)=.FALSE. !atomindx2 does not belong to the domain
            END IF

         END DO

      END DO

   END SUBROUTINE VLAddFullAACellInBoxLongCutOff
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE VLAddFullABCellInBoxLongCutOff()
   !full list, two different cells A and B (maincellindx/=neighcellindx)
   !CellOutsideBox is FALSE, potential cut off is long
      IMPLICIT NONE
      INTEGER :: i,j,istart,iend,jstart,jend
      INTEGER :: position,atomindx1,atomindx2
      INTEGER :: atomnumberneigh
      REAL(dp) :: dr(3),dr2mag,drmag
      
      istart=(maincellindx-1)*MaxAtomPerCell+1
      iend=istart+LLListRange(maincellindx)-1

      jstart=(neighcellindx-1)*MaxAtomPerCell+1
      jend=jstart+LLListRange(neighcellindx)-1
      
      
      DO i=istart,iend
         
         atomindx1=LLList(i)
         
         DO j=jstart,jend

            atomindx2=LLList(j)

            dr=AtomCoord(3*atomindx2-2:3*atomindx2)-AtomCoord(3*atomindx1-2:3*atomindx1)
            dr=dr-BoxSize*(NINT(dr/BoxSize)*AllowPBC)
            !dr=dr+REAL(BoxPosition)*BoxSize !since BoxSize==0
            dr2mag=DOT_PRODUCT(dr,dr)
            
            !dr2mag can be 0 when maincellindx==neighcellindx and atom is overlapping
            IF (dr2mag<=CutOff2) THEN

               drmag=SQRT(dr2mag)
               
               atomnumberneigh=VLListRange(atomindx1)+1
               
               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindx1-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindx2
               VLdr(3*position-2:3*position)=dr
               VLdrmag(position)=drmag
               VLListRange(atomindx1)=atomnumberneigh 
               !IF (CellOutsideBox) VLListDomainAtom(position)=.FALSE. !atomindx2 does not belong to the domain
               
               
               atomnumberneigh=VLListRange(atomindx2)+1
               
               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindx2-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindx1
               VLdr(3*position-2:3*position)=-dr
               VLdrmag(position)=drmag
               VLListRange(atomindx2)=atomnumberneigh 
               
            END IF
            
         END DO
         
      END DO
      
   END SUBROUTINE VLAddFullABCellInBoxLongCutOff
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE VLAddFullCellOutBoxLongCutOff()
   !full list, same cells A and A or two different cells A and B (maincellindx/=neighcellindx)
   !CellOutsideBox is TRUE, potential cut off is long
      IMPLICIT NONE
      INTEGER :: i,j,istart,iend,jstart,jend
      INTEGER :: position,atomindx1,atomindx2
      INTEGER :: atomnumberneigh
      REAL(dp) :: dr(3),dr2mag,drmag
      
      istart=(maincellindx-1)*MaxAtomPerCell+1
      iend=istart+LLListRange(maincellindx)-1

      jstart=(neighcellindx-1)*MaxAtomPerCell+1
      jend=jstart+LLListRange(neighcellindx)-1
      
      
      DO i=istart,iend
         
         atomindx1=LLList(i)
         
         DO j=jstart,jend

            atomindx2=LLList(j) !atom index outside box

            dr=AtomCoord(3*atomindx2-2:3*atomindx2)-AtomCoord(3*atomindx1-2:3*atomindx1)
            dr=dr-BoxSize*(NINT(dr/BoxSize)*AllowPBC)
            dr=dr+REAL(BoxPosition)*BoxSize
            dr2mag=DOT_PRODUCT(dr,dr)
            
            IF (dr2mag<=CutOff2) THEN

!write(*,*) "FullCellOutLongCutOff:",atomindx1,atomindx2
               drmag=SQRT(dr2mag)
               
               atomnumberneigh=VLListRange(atomindx1)+1

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindx1-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindx2
               VLdr(3*position-2:3*position)=dr
               VLdrmag(position)=drmag
               VLListRange(atomindx1)=atomnumberneigh
               VLListDomainAtom(position)=.FALSE.
!IF (atomindx1==1) write(*,*) "ooo",atomindx2
               
            END IF

         END DO

      END DO
     
   END SUBROUTINE VLAddFullCellOutBoxLongCutOff
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE VLAddHalfAACellInBoxLongCutOff()
   !half list, same cells A and A (maincellindx==neighcellindx)
   !CellOutsideBox is FALSE, potential cut off is long
      IMPLICIT NONE
      INTEGER :: i,j,istart,iend
      INTEGER :: position,atomindx1,atomindx2,atomindxmin,atomindxmax
      INTEGER :: atomnumberneigh
      REAL(dp) :: dr(3),dr2mag,drmag
      
      istart=(maincellindx-1)*MaxAtomPerCell+1
      iend=istart+LLListRange(maincellindx)-1
      
      DO i=istart,iend-1
         
         atomindx1=LLList(i) !atom index from first cell
         
         DO j=i+1,iend

            atomindx2=LLList(j) !atom index from other cell
            
            atomindxmin=MIN(atomindx1,atomindx2)
            atomindxmax=MAX(atomindx1,atomindx2)
            
            dr=AtomCoord(3*atomindxmax-2:3*atomindxmax)-AtomCoord(3*atomindxmin-2:3*atomindxmin)
            dr=dr-BoxSize*(NINT(dr/BoxSize)*AllowPBC)
            !dr=dr+REAL(BoxPosition)*BoxSize !BoxPosition==0
            dr2mag=DOT_PRODUCT(dr,dr)
            
            IF (dr2mag<=CutOff2) THEN
!write(*,*) "Half-AA-CellInBox-LOng",atomindx1,atomindx2
               drmag=SQRT(dr2mag)
                  
               atomnumberneigh=VLListRange(atomindxmin)+1

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindxmin-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindxmax
               VLdr(3*position-2:3*position)=dr
               VLdrmag(position)=drmag
               VLListRange(atomindxmin)=atomnumberneigh

            END IF

         END DO

      END DO
   
   END SUBROUTINE VLAddHalfAACellInBoxLongCutOff
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE VLAddHalfABCellInBoxLongCutOff()
   !half list, two different cells A and B (maincellindx/=neighcellindx)
   !CellOutsideBox is FALSE, potential cut off is long
      IMPLICIT NONE
      INTEGER :: i,j,istart,iend,jstart,jend
      INTEGER :: position,atomindx1,atomindx2,atomindxmin,atomindxmax
      INTEGER :: atomnumberneigh
      REAL(dp) :: dr(3),dr2mag,drmag
      
      istart=(maincellindx-1)*MaxAtomPerCell+1
      iend=istart+LLListRange(maincellindx)-1

      jstart=(neighcellindx-1)*MaxAtomPerCell+1
      jend=jstart+LLListRange(neighcellindx)-1
      
      
      DO i=istart,iend
         
         atomindx1=LLList(i) !atom index from first cell
         
         DO j=jstart,jend

            atomindx2=LLList(j) !atom index from other cell
            
            atomindxmin=MIN(atomindx1,atomindx2)
            atomindxmax=MAX(atomindx1,atomindx2)

            dr=AtomCoord(3*atomindxmax-2:3*atomindxmax)-AtomCoord(3*atomindxmin-2:3*atomindxmin)
            dr=dr-BoxSize*(NINT(dr/BoxSize)*AllowPBC)
            !dr=dr+REAL(BoxPosition)*BoxSize !since BoxPosition==0
            dr2mag=DOT_PRODUCT(dr,dr)
            
            IF (dr2mag<=CutOff2) THEN

!write(*,*) "Half-AB-CellInBox-LOng",atomindx1,atomindx2
               drmag=SQRT(dr2mag)
               
               atomnumberneigh=VLListRange(atomindxmin)+1

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindxmin-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindxmax
               VLdr(3*position-2:3*position)=dr
               VLdrmag(position)=drmag
               VLListRange(atomindxmin)=atomnumberneigh

            END IF

         END DO

      END DO
   
   END SUBROUTINE VLAddHalfABCellInBoxLongCutOff
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   SUBROUTINE VLAddHalfCellOutBoxLongCutOff()
   !half list, same cells A and A (maincellindx==neighcellindx) or different cells A and B
   !CellOutsideBox is TRUE, potential cut off is long
      IMPLICIT NONE
      INTEGER :: i,j,istart,iend,jstart,jend
      INTEGER :: position,atomindx1,atomindx2
      INTEGER :: atomnumberneigh
      REAL(dp) :: dr(3),dr2mag,drmag
      
      istart=(maincellindx-1)*MaxAtomPerCell+1
      iend=istart+LLListRange(maincellindx)-1

      jstart=(neighcellindx-1)*MaxAtomPerCell+1
      jend=jstart+LLListRange(neighcellindx)-1
      
      
      DO i=istart,iend
         
         atomindx1=LLList(i)
         
         DO j=jstart,jend !IMPORTANT: Make i interact with all atoms

            atomindx2=LLList(j) !atom index from cell outside the box

            dr=AtomCoord(3*atomindx2-2:3*atomindx2)-AtomCoord(3*atomindx1-2:3*atomindx1)
            dr=dr-BoxSize*(NINT(dr/BoxSize)*AllowPBC)
!IF ((maincellindx==1 .OR. neighcellindx==141) .AND. (maincellindx==141 .OR. neighcellindx==1)) THEN
! write(*,*) "Atom1,2:",atomindx1,atomindx2
! write(*,*) "drstg1:",dr
! END IF
            dr=dr+REAL(BoxPosition,dp)*BoxSize
!IF ((maincellindx==1 .OR. neighcellindx==141) .AND. (maincellindx==141 .OR. neighcellindx==1)) THEN
!write(*,*) "drstg2:",dr
! END IF
            dr2mag=DOT_PRODUCT(dr,dr)
!IF ((maincellindx==1 .OR. neighcellindx==141) .AND. (maincellindx==141 .OR. neighcellindx==1)) THEN
!write(*,*) "drmag:",sqrt(dr2mag)
!END IF
            
            IF (dr2mag<=CutOff2) THEN

!write(*,*) "Half-CellOutBox-LOng",atomindx1,atomindx2,maincellindx,neighcellindx
               drmag=SQRT(dr2mag)
               
               !cell is outside the box so atom has to be added no matter what
               atomnumberneigh=VLListRange(atomindx1)+1

               IF (atomnumberneigh>MaxAtomPerAtom) THEN
                  WRITE(6,*) "$Err>> Max atoms per atom in VL is small"
                  WRITE(UNIT=6,FMT='("... current # neighbors:",I3)') atomnumberneigh
                  WRITE(UNIT=6,FMT='("... maximum # neighbors:",I3)') MaxAtomPerAtom
                  STOP
               END IF

               position=(atomindx1-1)*MaxAtomPerAtom+atomnumberneigh
               VLList(position)=atomindx2 !add atom2
               VLdr(3*position-2:3*position)=dr
               VLdrmag(position)=drmag
               VLListRange(atomindx1)=atomnumberneigh
               VLListDomainAtom(position)=.FALSE.
               
            END IF
            
         END DO

      END DO
   
   END SUBROUTINE VLAddHalfCellOutBoxLongCutOff
   
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
