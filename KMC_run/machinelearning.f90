MODULE machinelearning
!contains tools for performing supervised and unsupervised machine learning
   USE VARIABLE_TYPE
   IMPLICIT NONE
   
   CONTAINS
   
   !==========================================================================================
   !Classification algorithms (supervised algorithms predicting categorical labels)
   !Maximum entropy classifier (aka logistic regression, multinomial logistic regression): Note that logistic regression is an algorithm for classification, despite its name. (The name comes from the fact that logistic regression uses an extension of a linear regression model to model the probability of an input being in a particular class.)
   !Naive Bayes classifier
   !Decision trees, decision lists
   !Support vector machines
   !Kernel estimation and K-nearest-neighbor algorithms
   !Perceptrons
   !Neural networks (multi-level perceptrons)
   !==========================================================================================
   !Clustering algorithms (unsupervised algorithms predicting categorical labels)
   !Categorical mixture models
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE kmeansClustering
   !K-means clustering
   END SUBROUTINE kmeansClustering
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE GaussianMixtureModel
   !Guassian mixture model using expectation-maximization algorithm
      IMPLICIT NONE
   END SUBROUTINE GaussianMixtureModel
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE kNearestNeighborClustering
   END SUBROUTINE kNearestNeighborClustering
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION HammingDistance(v1,v2,n,w)
   !find the HammingDistance between two arrays
      IMPLICIT NONE
      REAL(dp) :: HammingDistance
      INTEGER :: n,i
      REAL(dp), DIMENSION(n) :: v1,v2
      REAL(dp), DIMENSION(n) :: w !weight associated with each element of the array
      
      STOP
      HammingDistance=0._dp
      DO i=1,n
         IF (v1(i)/=v2(i)) HammingDistance=HammingDistance+ABS(w(i))
      END DO
   END FUNCTION HammingDistance
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   !Hierarchical clustering (agglomerative or divisive)
   !Kernel principal component analysis (Kernel PCA)
   !==========================================================================================
   ![edit]Regression algorithms (predicting real-valued labels)
   !Supervised:
   !Linear regression and extensions
   !Neural networks
   !Gaussian process regression (kriging)
   !Unsupervised:
   !Principal components analysis (PCA)
   !Independent component analysis (ICA)
   !==========================================================================================
   ![edit]Categorical sequence labeling algorithms (predicting sequences of categorical labels)
   !Supervised:
   !Hidden Markov models (HMMs)
   !Maximum entropy Markov models (MEMMs)
   !Conditional random fields (CRFs)
   !Unsupervised:
   !Hidden Markov models (HMMs)
   !==========================================================================================
   ![edit]Real-valued sequence labeling algorithms (predicting sequences of real-valued labels)
   !Kalman filters
   !Particle filters
   !==========================================================================================
   ![edit]Parsing algorithms (predicting tree structured labels)
   !Supervised and unsupervised:
   !Probabilistic context free grammars (PCFGs)
   !==========================================================================================
   ![edit]General algorithms for predicting arbitrarily-structured labels
   !Bayesian networks
   !Markov random fields
   !==========================================================================================
   ![edit]Ensemble learning algorithms (supervised meta-algorithms for combining multiple learning algorithms together)
   !Bootstrap aggregating ("bagging")
   !Boosting
   !Ensemble averaging
   !Mixture of experts, hierarchical mixture of experts
   !==========================================================================================
   SUBROUTINE RecognizeClusterAtoms(AtomIndex,AtomCoordList,natoms,rcutoff,IsClusterAtom,nClusterAtoms)
   !AtomCoord is the coordinates of the atom which is a part of the cluster
   !AtomCoordList is the coordinates of all natoms atoms from which some of the atoms are part of
   ! the cluster
   !natoms is the number of atoms in AtomCoordList
   !rcutoff is the cutoff distance for identifying neighbors for an atom
   !ClusterAtoms will contain the atom indices (from AtomCoordList) which belong to the cluster
   !nClusterAtoms is the number of atoms in the cluster
      IMPLICIT NONE
      INTEGER :: natoms,nClusterAtoms,AtomIndex
      REAL(dp), DIMENSION(3) :: AtomCoord
      REAL(dp), DIMENSION(3*natoms) :: AtomCoordList
      REAL(dp) :: rcutoff
      LOGICAL, DIMENSION(natoms) :: IsClusterAtom
      LOGICAL :: NotSatisfied
      
      !AtomCoord=AtomCoordList(3*AtomIndex-2:3*AtomIndex)
      !IsClusterAtom=.FALSE.
      !IsClusterAtom(AtomIndex)=.TRUE.
      !NotSatisfied=.TRUE.
      !NewClusterAtoms(1)=AtomIndex
      !nNewClusterAtoms=1
      !DO WHILE (NotSatisfied)
      !   nNewClusterAtoms1=nNewClusterAtoms
      !   DO i=1,natoms
      !      CALL 
      !   END DO
      !END DO
   END SUBROUTINE RecognizeClusterAtoms
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE machinelearning
