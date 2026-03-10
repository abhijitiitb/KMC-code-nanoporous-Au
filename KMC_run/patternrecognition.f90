MODULE PatternRecognition

   USE VARIABLE_TYPE
   USE db_manipulate
   USE NeighborList
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE CompareALToALNearestNeighbors(AL1,AL2,cutoffoptimized,cutoffunoptimized,IsMatch,errorstatus)
   !compares whether the nearest neighbors are same
   !REQUIRES VerletList
   !AL1 is the optimized reference state, while AL2 is (unoptimized)
   !this is not a fool-proof subroutine, 
   ! -- chances are that it will think a transition has occurred, but actually no transition occurred
   ! -- is it possible there are chances that it will miss transitions (not clear to me so far ..)
   !however, this is lower cost compared to the CompareALToAL with Optimization
   !AL1 is the reference state and may comprise of more than one image
   !AL2 is the state which has be compared to the reference state and may comprise of more than one image
   
   !subroutine is to be used in two steps:
   !first, when calling this for the first time setup cutoff, e.g., 
   !  CALL CompareALToALNearestNeighbors(cutoffoptimized=3.2,cutoffunoptimized=4.05)
   !later, you dont need to provide cutoff
   
   !this would work well in case of crystalline systems. in the case of amorphous materials defining
   !cutoffoptimized,cutoffunoptimized might be tricky
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER, OPTIONAL :: AL1,AL2
      TYPE(SystemContainer), POINTER :: ALCurr1,ALCurr2
      TYPE(VerletListContainer), POINTER :: VL1,VL2
      REAL(dp), OPTIONAL :: cutoffoptimized,cutoffunoptimized
      REAL(dp), SAVE :: cutoffoptimizedsaved=0._dp,cutoffunoptimizedsaved=0._dp !initialized as zero
      REAL(dp), DIMENSION(:), POINTER :: VLdrmag1,VLdrmag2
      INTEGER, DIMENSION(:), POINTER :: VLList1,VLList2,VLListRange1,VLListRange2
      REAL(dp), DIMENSION(:), POINTER :: AtomCoordOrig1,AtomCoordOrig2
      LOGICAL, DIMENSION(:), POINTER :: AtomIsMoving1
      LOGICAL, OPTIONAL :: IsMatch
      INTEGER :: j,j1,iatom,jatom,r1,weightedcount1,weightedcount2,errorstatus,MaxAtomPerAtom
      INTEGER :: nimages1,nimages2,image1,image2
      LOGICAL :: Found
      
      errorstatus=0 !perfect run
      
      IF (PRESENT(cutoffoptimized)) cutoffoptimizedsaved=cutoffoptimized
      IF (PRESENT(cutoffunoptimized)) cutoffunoptimizedsaved=cutoffunoptimized
      IF (cutoffoptimizedsaved==0._dp .OR. cutoffunoptimizedsaved==0._dp) THEN
         WRITE(6,*) "$Err>> Cutoff is zero in CompareALToALNearestNeighbors"
         STOP
      END IF
      
      IF (PRESENT(IsMatch)) THEN !do matching
         !check if AL present
         IF (.NOT. PRESENT(AL1) .OR. .NOT. PRESENT(AL2)) THEN
            WRITE(6,*) "$Err>> One or more atomlist missing in CompareALToALNearestNeighbors"
            errorstatus=1
            RETURN
         END IF
         !check if Verlet list is present
         IF (.NOT. ASSOCIATED(AL1%VL)) THEN
            WRITE(6,*) "$Err>> VL of one or more atomlist missing in CompareALToALNearestNeighbors"
            errorstatus=2
            RETURN
         END IF
         IF (.NOT. ASSOCIATED(AL2%VL)) THEN
            WRITE(6,*) "$Err>> VL of one or more atomlist missing in CompareALToALNearestNeighbors"
            errorstatus=3
            RETURN
         END IF
         !now we can compare
         AtomIsMoving1=>AL1%AtomIsMoving
         
         !nimages1=CountImages(AL1) !reference state
         !nimages2=CountImages(AL2) !state that needs to be compared to reference
         
         nimages1=1 !presently we shall force comparison between first images
         nimages2=1
         
         IF (nimages1==1 .AND. nimages2==1) THEN
            CALL CompareALToALNearestNeighborsDriver()
         ELSE !now we need to do some jugglery involving VL based comparison
            !goal is to modify the atom coordinates in the first image of the chain of states
            !then update the VL and check whether there is match -- if there is a match we are good to go
            AtomCoordOrig1=>AL1%AtomCoord
            AtomCoordOrig2=>AL2%AtomCoord
            ALCurr1=>AL1
            ALCurr2=>AL2
            DO image1=1,nimages1
               AL1%AtomCoord=>ALCurr1%AtomCoord
               CALL AddVerletList(AL1,NAtoms=AL1%NAtoms)
               DO image2=1,nimages2
                  AL2%AtomCoord=>ALCurr2%AtomCoord
                  CALL AddVerletList(AL2,NAtoms=AL2%NAtoms)
                  CALL CompareALToALNearestNeighborsDriver()
                  IF (IsMatch) EXIT
                  ALCurr2=>ALCurr2%NextNeigh
               END DO
               IF (IsMatch) EXIT
               ALCurr1=>ALCurr1%NextNeigh
            END DO
            AL1%AtomCoord=>AtomCoordOrig1
            CALL AddVerletList(AL1,NAtoms=AL1%NAtoms)
            AL2%AtomCoord=>AtomCoordOrig2
            CALL AddVerletList(AL2,NAtoms=AL2%NAtoms)
         END IF
      END IF
      
      CONTAINS
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      SUBROUTINE CompareALToALNearestNeighborsDriver()
         VL1=>AL1%VL
         VL2=>AL2%VL
         VLdrmag1=>VL1%drmag
         VLdrmag2=>VL2%drmag
         VLListRange1=>VL1%ListRange
         VLListRange2=>VL2%ListRange
         VLList1=>VL1%List
         VLList2=>VL2%List
         MaxAtomPerAtom=VL1%MaxAtomPerAtom
         IsMatch=.TRUE.
         DO iatom=1,AL1%NAtoms
            IF (.NOT. ALL(AtomIsMoving1(3*iatom-2:3*iatom))) CYCLE !we only consider moving atoms
            r1=(iatom-1)*MaxAtomPerAtom
            !weightedcount1=0
            DO j=r1+1,r1+VLListRange1(iatom)
               IF (VLdrmag1(j)<cutoffoptimizedsaved) THEN
                  Found=.FALSE.
                  jatom=VLList1(j)
                  DO j1=r1+1,r1+VLListRange2(iatom)
                     IF (jatom==VLList2(j1)) THEN
                        Found=.TRUE.
                        IF (VLdrmag2(j1)>cutoffunoptimizedsaved) THEN !it is outside cutoff
                           WRITE(6,*) "outside cutoff",VLdrmag2(j1),jatom
                           IsMatch=.FALSE.
                           RETURN
                        END IF
                        EXIT
                     END IF
                  END DO
                  IF (.NOT. Found) THEN
                     WRITE(*,*) "could not find"
                     IsMatch=.FALSE.
                     RETURN
                  END IF
               END IF
            END DO
         END DO
      END SUBROUTINE CompareALToALNearestNeighborsDriver
      !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   END SUBROUTINE CompareALToALNearestNeighbors
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

END MODULE PatternRecognition
