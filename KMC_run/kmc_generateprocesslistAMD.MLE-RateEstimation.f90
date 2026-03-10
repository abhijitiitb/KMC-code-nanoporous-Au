   FUNCTION TimeNPoisson(ntransitions,tstop)
   !finds the total time involving ntransitions when each successful transition took less than tstop time
   !time is given in dimensionless units, i.e., k*actual_time
      IMPLICIT NONE
      REAL(dp) :: TimeNPoisson,tstop,dt,rand1
      INTEGER :: ntransitions,itransition
      
      TimeNPoisson=0._dp
      DO itransition=1,ntransitions
         dt=tstop+0.1
         DO WHILE (dt>tstop)
            rand1=taus88()
            dt=-log(rand1)
         END DO
         TimeNPoisson=TimeNPoisson+dt
      END DO
   END FUNCTION TimeNPoisson
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION NTransitionTail(ntransitions,tstop,MaxTotalTime,nattempts)
   !finds how many attempts resulted in total time after NTransitions to fall
   !beyond time MaxTotalTime
      IMPLICIT NONE
      INTEGER :: iattempt,nattempts,ntransitions,NTransitionTail
      REAL(dp) :: ndt,tstop,MaxTotalTime
      
      NTransitionTail=0
      DO iattempt=1,nattempts
         ndt=TimeNPoisson(ntransitions,tstop)
         IF (ndt>MaxTotalTime) THEN
            NTransitionTail=NTransitionTail+1
         END IF
      END DO
   END FUNCTION NTransitionTail
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE FindTStop(ntransitions,MaxTotalTime,tstop,tstop0,nattempts_per_tstop,probability)
   !finds the tstop value that will result in the probability of total time greater than
   !MaxTotalTime specified by the probability that is provided
      IMPLICIT NONE
      INTEGER, PARAMETER :: ndata=10
      INTEGER :: ntransitions,nattempts_per_tstop,n0,n1,ntail,nA(ndata),nB(ndata),nC(ndata),iter
      REAL(dp) :: tstop0,MaxTotalTime,probability,tstop,tstopA,tstopB,tstopC
      REAL(dp) :: meanA,stdevA,stderrA,meanB,stdevB,stderrB,targetnumber
      REAL(dp) :: meanC,stdevC,stderrC
      LOGICAL :: NotSatisfied
      
      !Target number of events is 10
      targetnumber=10._dp
      nattempts_per_tstop=CEILING(targetnumber/probability) !this is how many we need -- probability is the probability of having time > MaxTotalTime
      WRITE(356,*) nattempts_per_tstop; CALL FLUSH(356)
      
      WRITE(356,*) "Number of transitions:",ntransitions
      WRITE(356,*) "Tail starts from:",MaxTotalTime
      WRITE(356,*) "Initial tstop guess:",tstop0
      
      !find lower bracket
      NotSatisfied=.TRUE.
      tstopA=tstop0
      WRITE(356,*) "Finding lower tstop bracket ..."
      DO WHILE (NotSatisfied)
         WRITE(356,*) "Current tstopA:",tstopA
         DO iter=1,ndata
            nA(iter)=NTransitionTail(ntransitions,tstopA,MaxTotalTime,nattempts_per_tstop)
            WRITE(356,*) "nA(",iter,")",nA(iter); CALL FLUSH(356)
         END DO
         meanA=sum(nA)/REAL(ndata,dp)
         stdevA=sum((nA-meanA)*(nA-meanA))/REAL(ndata-1,dp)
         stdevA=SQRT(stdevA)
         stderrA=stdevA/SQRT(REAL(ndata,dp))
         WRITE(356,*) "A statistics:",meanA,"+/-",stderrA; CALL FLUSH(356)
         IF (meanA+stderrA<targetnumber) THEN
            NotSatisfied=.FALSE.
         ELSE
            tstopA=tstopA/1.05_dp
         END IF
      END DO
      
      !find revised lower and upper bracket
      NotSatisfied=.TRUE.
      tstopB=tstopA*1.02 !this should give larger value
      WRITE(356,*) "Finding upper tstop bracket ..."
      DO WHILE (NotSatisfied)
         WRITE(356,*) "Current tstopB:",tstopB
         DO iter=1,10
            nB(iter)=NTransitionTail(ntransitions,tstopB,MaxTotalTime,nattempts_per_tstop)
            WRITE(356,*) "nB(",iter,")",nB(iter); CALL FLUSH(356)
         END DO
         meanB=sum(nB)/REAL(ndata,dp)
         stdevB=sum((nB-meanB)*(nB-meanB))/REAL(ndata-1,dp)
         stdevB=SQRT(stdevB)
         stderrB=stdevB/SQRT(REAL(ndata,dp))
         WRITE(356,*) "B statistics:",meanB,"+/-",stderrB; CALL FLUSH(356)
         IF (meanB-2._dp*stderrB>targetnumber) THEN
            NotSatisfied=.FALSE.
         ELSE
            IF (meanB+2._dp*stderrB<targetnumber) THEN
               tstopA=tstopB
               meanA=meanB
            END IF
            tstopB=tstopB*1.02_dp
         END IF
      END DO
      WRITE(356,*) "Found lower tstop bracket:",tstopA," with mean:",meanA
      WRITE(356,*) "Found upper tstop bracket:",tstopB," with mean:",meanB
      
      NotSatisfied=.TRUE.
      DO WHILE (NotSatisfied) !do binary search
         tstopC=tstopA+(tstopB-tstopA)*(targetnumber-REAL(meanA,dp))/(REAL(meanB,dp)-REAL(meanA,dp))
         DO iter=1,ndata
            nC(iter)=NTransitionTail(ntransitions,tstopC,MaxTotalTime,nattempts_per_tstop)
            WRITE(356,*) "nC(",iter,")",nC(iter); CALL FLUSH(356)
         END DO
         meanC=sum(nC)/REAL(ndata,dp)
         stdevC=sum((nC-meanC)*(nC-meanC))/REAL(ndata-1,dp)
         stdevC=SQRT(stdevC)
         stderrC=stdevC/SQRT(REAL(ndata,dp))
         WRITE(356,*) "C statistics:",meanC,"+/-",stderrC; CALL FLUSH(356)
         WRITE(356,*) "Time brackets:",tstopA,tstopB,tstopC
         WRITE(356,*) "Count brackets:",meanA,meanB,meanC
         
         NotSatisfied=ABS(meanC-targetnumber)>1.5_dp*stderrC
         IF (meanC>targetnumber) THEN
            WRITE(356,*) "Replace B with C"
            meanB=meanC
            stderrB=stderrC
            tstopB=tstopC
         ELSE
            WRITE(356,*) "Replace A with C"
            meanA=meanC
            stderrA=stderrC
            tstopA=tstopC
         END IF
      END DO
      WRITE(356,*) "Found the optimum stop time:",tstopC
      tstop=tstopC !output
   END SUBROUTINE FindTStop
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE FindTStopSet()
      IMPLICIT NONE
      INTEGER :: ntransitions,nattempts_per_tstop,ntransitions_set(31),i,n
      REAL(dp) :: tstopguess,tstop0,tstop
      
      tstopguess=1.2_dp
      !n=31
      !OPEN(UNIT=355,FILE="TstopResults")
      !OPEN(UNIT=356,FILE="TstopDetailedResults")
      !ntransitions_set=(/10,15,20,25,30,40,50,60,70,80,100,120,150,180,210,240, &
      !  270,300,350,400,450,500,600,700,800,900,1000,1200,1400,1600,1800/) !31 different transitions
      
      !Set 1
      n=3
      tstopguess=1.2_dp
      OPEN(UNIT=355,FILE="TstopResults.1")
      OPEN(UNIT=356,FILE="TstopDetailedResults.1")
      ntransitions_set(1:n)=(/2500,4000,8000/)
      
      !Set 2
      !n=3
      !OPEN(UNIT=355,FILE="TstopResults.2")
      !OPEN(UNIT=356,FILE="TstopDetailedResults.2")
      !ntransitions_set(1:n)=(/12000,18000,25000/)
      
      !Set 3
      !n=1
      !OPEN(UNIT=355,FILE="TstopResults.3")
      !OPEN(UNIT=356,FILE="TstopDetailedResults.3")
      !ntransitions_set(1:n)=(/32000/)
      
      !Set 4
      !n=1
      !OPEN(UNIT=355,FILE="TstopResults.4")
      !OPEN(UNIT=356,FILE="TstopDetailedResults.4")
      !ntransitions_set(1:n)=(/60000/)
      
      !Set 5
      !n=1
      !OPEN(UNIT=355,FILE="TstopResults.5")
      !OPEN(UNIT=356,FILE="TstopDetailedResults.5")
      !ntransitions_set(1:n)=(/100000/)
      WRITE(355,*) "Ntransitions   tstopfound"
      DO i=1,n
         ntransitions=ntransitions_set(i)
         tstopguess=-0.0071_dp+1.1635_dp*log10(REAL(ntransitions))
         CALL FindTStop(ntransitions=ntransitions,MaxTotalTime=REAL(ntransitions,dp),tstop=tstop, &
            tstop0=tstopguess, &
            nattempts_per_tstop=nattempts_per_tstop,probability=1.e-5_dp)
         WRITE(355,*) ntransitions,tstop
         CALL FLUSH(355)
         tstopguess=tstop !guess for next ntransitions
      END DO
      CLOSE(355)
      CLOSE(356)
   END SUBROUTINE FindTStopSet
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE InitializeTStop
   !sets up arrays for calculating the total rate constant estimates
      IMPLICIT NONE
      INTEGER :: idata,m
      
      OPEN(UNIT=135,FILE="BEStstop.input")
      READ(135,*) NTarget
      WRITE(6,*) "Initializing the target stop times"
      WRITE(6,*) "----------------------------"
      WRITE(6,*) "       m        Stop time"
      WRITE(6,*) "----------------------------"
      DO idata=1,NTarget
         READ(135,*) m
         TargetNBasinEscapes(idata)=m
         TargetStopTime(idata)=InterpTStop(m)
         WRITE(6,*) m,TargetStopTime(idata)
      END DO
      CLOSE(135)
      WRITE(6,*) "----------------------------"
   END SUBROUTINE InitializeTStop
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION InterpTStop(m)
   !obtains a tstop value from the array TargetStopTime
      IMPLICIT NONE
      INTEGER, DIMENSION(37) :: marr=(/10,15,20,25,30,40,50,60,70,80,100,120,150,180,210,240,270,300,350,400, &
         450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,3000,4000,6000,15000,70000/)
      REAL(dp), DIMENSION(37) :: tsarr=(/1.2800,1.4300,1.5400,1.6300, &
         1.7100,1.8400,1.9450,2.0310,2.1080,2.1680,2.2820,2.3700, &
         2.4880,2.5742,2.6600,2.7200,2.8010,2.8490,2.9248,2.9980, &
         3.0690,3.1236,3.2290,3.3100,3.3840,3.4510,3.5110,3.6090,3.6970,3.7745, &
         3.8375,3.8929,4.1340,4.3069,4.5359,5.0902,5.9518/)
      REAL(dp) :: InterpTStop
      INTEGER :: m,i
      
      IF (m>marr(37)) THEN
         WRITE(6,*) "Err>> The value of m in TargetNBasinEscapes exceeds 70000"
         STOP
      END IF
      IF (m<marr(1)) THEN
         WRITE(6,*) "Err>> The value of m in TargetNBasinEscapes is less than 10"
         STOP
      END IF
      
      DO i=2,37
         IF (m<marr(i)) EXIT
      END DO
      InterpTStop=tsarr(i)-(tsarr(i)-tsarr(i-1))*REAL(marr(i)-m)/REAL(marr(i)-marr(i-1))
   END FUNCTION InterpTStop
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
