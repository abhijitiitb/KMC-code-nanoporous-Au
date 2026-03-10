MODULE StatisticsPackage
   USE VARIABLE_TYPE
   IMPLICIT NONE
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION RSquared(y,f,n)
   !y is target, f is model output, n is number of data points
      IMPLICIT NONE
      INTEGER :: n
      REAL(dp), DIMENSION(n) :: y,f
      REAL(dp) :: RSquared,meany
      
      RSquared=0._dp
      meany=SUM(y(1:n))/REAL(n,dp)
      RSquared=1._dp-SUM((y(1:n)-f(1:n))*(y(1:n)-f(1:n)))/ &
         SUM((y(1:n)-meany)*(y(1:n)-meany))
   END FUNCTION RSquared
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE StatisticsPackage
