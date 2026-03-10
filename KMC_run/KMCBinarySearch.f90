MODULE KMC_BinarySearch
!performs binary searches with the KMC algorithm
!requires a fixed array of rates which will be updated externally
!the module maintains a database to update the sum of rates 
!so that a process can be efficiently selected with binary search

!algorithm used - create a hierarchy of arrays to update sum of rates
   USE KMC_VARIABLE_TYPE
   USE Ecuyer_random
   IMPLICIT NONE
   
   TYPE BinaryTreeLayer
      REAL(dp), DIMENSION(:), POINTER :: RateArray
      INTEGER :: n !size of the rate array
      INTEGER :: LayerNo=0 !position of the layer - top layer is layer 1
      TYPE(BinaryTreeLayer), POINTER :: NextNeigh=>NULL(),PrevNeigh=>NULL()
   END TYPE BinaryTreeLayer
   TYPE(BinaryTreeLayer), POINTER :: BinaryTree=>NULL()
   INTEGER :: NLayers=0 !number of layers in the binary tree
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE BinaryTreeInitialize(rate,n)
   !rate: the rate array to be used for KMC search
   !n: size of the rate array
      IMPLICIT NONE
      INTEGER :: n,nrates_in_new_layer,ilayer
      REAL(dp), DIMENSION(n), TARGET :: rate
      TYPE(BinaryTreeLayer), POINTER :: BT=>NULL()
      LOGICAL :: NotSatisfied
      
      !create the binary tree
      nrates_in_new_layer=n
      ALLOCATE(BinaryTree)
      BinaryTree%RateArray=>rate
      BinaryTree%n=n
      NLayers=1
      
      NotSatisfied=nrates_in_new_layer>1
      DO WHILE (NotSatisfied) !create a new layer while not satisfied
         nrates_in_new_layer=CEILING(REAL(nrates_in_new_layer)/2.)
         ALLOCATE(BinaryTree%PrevNeigh)
         BinaryTree%PrevNeigh%NextNeigh=>BinaryTree
         BinaryTree=>BinaryTree%PrevNeigh
         ALLOCATE(BinaryTree%RateArray(nrates_in_new_layer))
         BinaryTree%n=nrates_in_new_layer
         NotSatisfied=nrates_in_new_layer>1
         NLayers=NLayers+1
      END DO
  !    WRITE(6,*) NLayers
      
      !number the layers
      BT=>BinaryTree
      DO ilayer=1,NLayers
         BT%LayerNo=ilayer
         IF (ilayer<NLayers) BT=>BT%NextNeigh
      END DO
      
      !set up the rates (note BT is pointing to the last layer)
      CALL BinaryTreeRefresh()
   END SUBROUTINE BinaryTreeInitialize
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE BinaryTreeRefresh()
   !refresh sum of rates for entire BTree to clear away truncation error accumulated over time
      IMPLICIT NONE
      TYPE(BinaryTreeLayer), POINTER :: BT
      REAL(dp), DIMENSION(:), POINTER :: RateArray,RateArrayParent
      INTEGER :: idata,ilayer
      
      !Layer 1 is always the latest one so it does not need to be updated
      BT=>BinaryTree
      DO WHILE (ASSOCIATED(BT%NextNeigh))
         BT=>BT%NextNeigh
      END DO
      !now that we reached last layer, update previous layers
      DO ilayer=NLayers-1,1,-1
         RateArray=>BT%RateArray
         RateArrayParent=>BT%PrevNeigh%RateArray
         DO idata=1,BT%n/2
            RateArrayParent(idata)=RateArray(2*idata-1)+RateArray(2*idata)
    !      WRITE(*,*) RateArray(2*idata-1),RateArray(2*idata),">",RateArrayParent(idata)
         END DO
         IF (MOD(BT%n,2)==1) RateArrayParent(BT%n/2+1)=RateArray(BT%n) !account odd number of rates in RateArray
         BT=>BT%PrevNeigh
         !STOP
      END DO
   END SUBROUTINE BinaryTreeRefresh
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintBinaryTree(BT,iunit)
      IMPLICIT NONE
      TYPE(BinaryTreeLayer), POINTER :: BT,BT1
      REAL(dp), DIMENSION(:), POINTER :: RateArray
      INTEGER :: ilayer,iunit,idata
      
      BT1=>BT
      DO ilayer=1,NLayers
         WRITE(UNIT=iunit,FMT='("Layer # ",i3)') BT1%LayerNo
         WRITE(iunit,*) "  rate information ..."
         RateArray=>BT1%RateArray
         DO idata=1,BT1%n
            WRITE(UNIT=iunit,FMT='(F10.3," ")',ADVANCE="NO") RateArray(idata)
         END DO
         BT1=>BT1%NextNeigh
         WRITE(iunit,*) ""
      END DO
   END SUBROUTINE PrintBinaryTree
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE BinaryTreeSearch(BT,position)
   !perform binary search and provide the position in the rate array that was selected
      IMPLICIT NONE
      INTEGER :: position
      TYPE(BinaryTreeLayer), POINTER :: BT,BT1
      REAL(dp), DIMENSION(:), POINTER :: RateArrayChild
      REAL(dp) :: r
      INTEGER :: ilayer
      
      position=1
      r=taus88()*BT%RateArray(1) !random number
      BT1=>BT
!      WRITE(6,*) "Random rate:",r,BT1%RateArray(1)

      DO ilayer=2,NLayers
         RateArrayChild=>BT1%NextNeigh%RateArray
   !      WRITE(6,*) "Current layer : ",ilayer,"left rate:",RateArrayChild(2*position-1)
         IF (RateArrayChild(2*position-1)>=r) THEN !select left child
            position=2*position-1
         ELSE !select right child provided it exist
            r=r-RateArrayChild(2*position-1)
            position=2*position
         END IF
   !      WRITE(6,*) "selected position in next layer:",position," search rate ",r
         IF (position>BT1%NextNeigh%n) THEN
   !         WRITE(6,*) "Err>> Rate not found in array ... stopping"
            STOP
         END IF
         BT1=>BT1%NextNeigh
      END DO
   END SUBROUTINE BinaryTreeSearch
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE BinaryTreeLocalUpdate(position,m)
   !the subroutine updates the sum of rates 
   !given the positions in the rate array where the rates are affected
   !subroutine updates the entire binary tree
      IMPLICIT NONE
      INTEGER :: m,Parentm,Childm !size of the position array
      TYPE(BinaryTreeLayer), POINTER :: BT
      INTEGER, DIMENSION(m) :: position
      REAL(dp), DIMENSION(:), POINTER :: RateArrayParent,RateArrayChild
      INTEGER, DIMENSION(1000) :: Parentposition=0,Childposition=0 !assuming there will be no more than 1000 processes to be updated
      INTEGER :: n,idata,pos,Parentpos,ilayer
      
      !Layer 1 is always the latest one so it does not need to be updated
      BT=>BinaryTree
      DO WHILE (ASSOCIATED(BT%NextNeigh))
         BT=>BT%NextNeigh
      END DO
      
      Childposition(1:m)=position(1:m)
      Childm=m
      !now that we reached last layer, update previous layers
      DO ilayer=NLayers-1,1,-1
         RateArrayChild=>BT%RateArray
         RateArrayParent=>BT%PrevNeigh%RateArray
         Parentm=0 !position in parent
         n=BT%n
         DO idata=1,Childm
            pos=Childposition(idata)
            Parentpos=(pos+1)/2 !this is how the data is filled
            Parentm=Parentm+1
            Parentposition(Parentm)=Parentpos
            IF (2*Parentpos<=n) THEN
               RateArrayParent(Parentpos)=RateArrayChild(2*Parentpos-1)+RateArrayChild(2*Parentpos)
            ELSE
               RateArrayParent(Parentpos)=RateArrayChild(2*Parentpos-1)
            END IF
         END DO
         Childposition(1:Parentm)=Parentposition(1:Parentm)
         Childm=Parentm
         BT=>BT%PrevNeigh
      END DO
   END SUBROUTINE BinaryTreeLocalUpdate
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE KMC_BinarySearch
