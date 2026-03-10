MODULE KMCProcessSearchParallel
!Principal Author: Abhijit Chatterjee, Department of Chemical Engineering, IIT Kanpur 
!Provided by user find the processes possible in the current system - performs parallel calculations
   USE VARIABLE_TYPE
   USE KMCNEBProcessSearch, ONLY : AddKMCProcessInitialState,AddKMCProcessFinalState,ProcessNEB, &
      PrcInitialCoord,PrcSaddleCoord,PrcFinalCoord,ProcessAtom
   USE mpimanager
   USE MDPackage
   USE mdmpi_domaindecomposition
   USE mpi
   IMPLICIT NONE
   !INCLUDE "mpif.h"
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PerformNEBParallel(nfiles)
   !executed by the master processor
      IMPLICIT NONE
      INTEGER :: ifile1
      INTEGER :: master,istatus,jobid
      INTEGER :: ifile,nfiles,nremain,proc,i
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status=0
      LOGICAL :: NotSatisfied

      !CALL MPI_INIT(ierr)
      !CALL MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)
      !CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
      
      master=0
      jobid=1 !NEB type calculations is jobid 1
      
      !CALL SetNumProcs(num_procs1-1)
      !CALL ResetProcStatus()
      CALL ReserveProc(num_procs,jobid=1,istatus=istatus)
      IF (istatus==1) THEN
         WRITE(6,*) "Unable to obtain any processors for mpi_neb"
         STOP
      END IF
         
      ifile=0
      NotSatisfied= ifile<nfiles
      DO WHILE (NotSatisfied)
            
         proc=GetFreeProc(jobid=jobid)
         
         IF (proc==0) THEN !we have no free processors
            CALL MPI_RECV(proc,1,MPI_INTEGER,MPI_ANY_SOURCE,TAG,MPI_COMM_WORLD,status,ierr)
            CALL SetProcStatus(proc=proc,procstatus=0,jobid=jobid)
         END IF
         
         IF (proc>0) THEN !submit job
            ifile=ifile+1
            NotSatisfied= ifile<nfiles !this will be the last job when  ifile==nfiles
            CALL SetProcStatus(proc=proc,procstatus=1,jobid=jobid)
            TAG=0
            CALL MPI_SEND(jobid,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(ifile,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,ierr)
         END IF
      END DO
      
      !get the remaining jobs back
      nremain=GetNumProcs(jobid=jobid,procstatus=1)
      DO i=1,nremain
         CALL MPI_RECV(proc,1,MPI_INTEGER,MPI_ANY_SOURCE,TAG,MPI_COMM_WORLD,status,ierr)
         CALL SetProcStatus(proc=proc,procstatus=0,jobid=jobid)
      END DO
      !DO proc=1,num_procs !send negative number as a signal to stop
      !   CALL MPI_SEND(-1,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,ierr)
      !END DO
      
      !CALL MPI_FINALIZE(ierr)
      !IF (ierr/=0) THEN
      !   WRITE(6,*) "$Err>> MPI_MD_DOMAINDECOMPOSITION failed to finalize MPI"
      !END IF
   END SUBROUTINE PerformNEBParallel
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PerformMDProcessSearchParallel(stateindex,nsearch,Temperature,DetectTransitionFrequency, &
      MDTimePerBCMD,TotalMDTime,nsuccess,StoppingOption,ObtainRate)
   !uses MD to search for processes using parallel computation
   !called by MASTER processor
   !stateindex is initial state
   !nsearch is number of search attempts
   !irank is the rank of the processor
   !nsuccess contains the number of successful searches
   !StoppingOption= 0 (perform nsearches with stop time given by MDTimePerBCMD)
   !StoppingOption= 1 (perform nsearches with stop time given by MDTimePerBCMD(if it is present), but later modify 
   !   MDTimePerBCMD when a better value of MDTimePerBCMD is known. this is useful when the rates 
   !   from a state are completely unknown)
      IMPLICIT NONE
      INTEGER :: isearch1
      REAL(dp) :: Temperature,ForcedBCMDTime
      INTEGER :: master,istatus,jobid,isearchcompleted
      INTEGER :: isearch,nsearch,nremain,proc,i,iseed_amd
      INTEGER :: DetectTransitionFrequency,stateindex
      REAL(sp) :: ran1scl
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status=0
      LOGICAL :: NotSatisfied,success
      REAL(dp) :: MDTimePerBCMD,TotalMDTime,dTime !,TotalMDTimeWithStopOption1
      INTEGER :: nsuccess,StoppingOption,StoppingOption1,msgtype,ObtainRate
      CHARACTER(len=100) :: filename
      CHARACTER(len=NINT2CHAR) :: charstateindx,charisearch
      REAL(dp) :: Rate

      !CALL MPI_INIT(ierr)
      !CALL MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)
      !CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)
      
      master=0
      jobid=2 !MD process search calculations is jobid 2
      StoppingOption1=StoppingOption
      
      !CALL SetNumProcs(num_procs1-1)
      !CALL ResetProcStatus()
      nsuccess=0
      TotalMDTime=0._dp
      ForcedBCMDTime=MDTimePerBCMD
      !TotalMDTimeWithStopOption1=0._dp
      !OPEN(UNIT=314,FILE="MDTransitionSearchStatus")
      !WRITE(314,*) MDTimePerBCMD !all processors have to use this time
      !CLOSE(314)
      
      CALL ReserveProc(nsearch,jobid=jobid,istatus=istatus)
      IF (istatus==1) THEN
         WRITE(6,*) "Unable to obtain any processors for mpi_md"
         STOP
      END IF
      isearch=0
      NotSatisfied= isearch<nsearch
      

!OPEN(UNIT=377,FILE="par.out")
      DO WHILE (NotSatisfied)
            
         proc=GetFreeProc(jobid=jobid)
         
         IF (proc==0) THEN !we have no free processors
            CALL MPI_RECV(proc,1,MPI_INTEGER,MPI_ANY_SOURCE,TAG,MPI_COMM_WORLD,status,ierr) !wait for someone to reply back
            CALL MPI_RECV(msgtype,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,status,ierr) !understand what type of msg is being sent
            IF (msgtype==1) THEN !msgtype==1 implies the job is done and MD time information will be supplied
               CALL MPI_RECV(success,1,MPI_LOGICAL,proc,TAG,MPI_COMM_WORLD,status,ierr)
               IF (success) nsuccess=nsuccess+1 !MD successfully found a transition
               CALL MPI_RECV(dTime,1,MPI_DOUBLE_PRECISION,proc,TAG,MPI_COMM_WORLD,status,ierr)
            write(378,*) "processor",proc,dTime; CALL FLUSH(378)
               TotalMDTime=TotalMDTime+dTime
               CALL MPI_RECV(isearchcompleted,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,status,ierr)
!WRITE(377,*) "MD completed - proc, isearch:",proc,isearchcompleted
               CALL SetProcStatus(proc=proc,procstatus=0,jobid=jobid)
               ForcedBCMDTime=0._dp

            ELSEIF (msgtype==2) THEN !processor is asking for stopping time information
               IF (StoppingOption1==1) THEN !MD finds transition only before ForcedBCMD
                  CALL MPI_SEND(ForcedBCMDTime,1,MPI_DOUBLE_PRECISION,proc,TAG,MPI_COMM_WORLD,ierr)
               ELSEIF (StoppingOption1==0) THEN !let MD find a transition
                  CALL MPI_SEND(MDTimePerBCMD,1,MPI_DOUBLE_PRECISION,proc,TAG,MPI_COMM_WORLD,ierr)
               ELSE
                  WRITE(6,*) "Err>> StoppingOption1 is having wierd value"
                  STOP
               END IF
               proc=0 !no processor free yet
            ELSE
               WRITE(6,*) "Err>> Message from compute processor was not understood"
               STOP
            END IF

            !IF (StoppingOption1==1) THEN
            !   charstateindx=INT2CHAR(stateindex)
            !   charisearch=INT2CHAR(isearchcompleted)
            !   filename="./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".param"
            !   OPEN(UNIT=314,FILE=filename)
            !   READ(314,*) ForcedBCMDTime !time elapsed
            !   CLOSE(314)
            !   !IF (StoppingOption1==1) TotalMDTimeWithStopOption1=ForcedBCMDTime
               
            !   OPEN(314,FILE="MDTransitionSearchStatus1") !presently the MD does not read the file regularly. we need to implement this
            !   WRITE(314,*) 0._dp !ForcedBCMDTime*2._dp
            !   CLOSE(314)
            !   CALL SYSTEM("mv MDTransitionSearchStatus1 MDTransitionSearchStatus")
            !   StoppingOption1=0 !I am removing this since this serves no purpose
            !END IF
         END IF
           
         IF (proc>0) THEN !submit job
            isearch=isearch+1
            NotSatisfied= isearch<nsearch !this will be the last job when  isearch==nsearchs
            CALL SetProcStatus(proc=proc,procstatus=1,jobid=jobid)
!WRITE(377,*) "Allocating process and isearch:",proc,isearch; CALL FLUSH(377)
            
            !MPI is working oddly in some situations. By having a send and recieve the error
            !was observed to be removed
            TAG=0
            CALL MPI_SEND(jobid,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,ierr)
            IF (ierr/=0) THEN
               WRITE(6,*) "Err>> Error in sending message"
               STOP
            END IF
            
            CALL MPI_SEND(isearch,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(stateindex,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(Temperature,1,MPI_DOUBLE_PRECISION,proc,TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(DetectTransitionFrequency,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(ObtainRate,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,ierr)
            CALL ran1(ran1scl)
            iseed_amd=INT(ran1scl*10000000)
            CALL MPI_SEND(iseed_amd,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,ierr)
         END IF
      
      END DO

      !get the remaining jobs back
      nremain=GetNumProcs(jobid=jobid,procstatus=1)
      DO WHILE (nremain>0)
         CALL MPI_RECV(proc,1,MPI_INTEGER,MPI_ANY_SOURCE,TAG,MPI_COMM_WORLD,status,ierr)
         CALL MPI_RECV(msgtype,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,status,ierr) !understand what type of msg is being sent
         IF (msgtype==1) THEN !msgtype==1 implies the job is done and MD time information will be supplied
            CALL MPI_RECV(success,1,MPI_LOGICAL,proc,TAG,MPI_COMM_WORLD,status,ierr)
            IF (success) nsuccess=nsuccess+1 !MD successfully found a transition
            CALL MPI_RECV(dTime,1,MPI_DOUBLE_PRECISION,proc,TAG,MPI_COMM_WORLD,status,ierr)
            TotalMDTime=TotalMDTime+dTime
            CALL MPI_RECV(isearchcompleted,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,status,ierr)
!WRITE(377,*) "MD completed - proc, isearch:",proc,isearchcompleted
            CALL SetProcStatus(proc=proc,procstatus=0,jobid=jobid)
            nremain=nremain-1
            ForcedBCMDTime=0._dp !stop all MD right now
!WRITE(378,*) "Received answer from ",proc
         ELSEIF (msgtype==2) THEN !processor is asking for stopping time information
            IF (StoppingOption1==1) THEN
               CALL MPI_SEND(ForcedBCMDTime,1,MPI_DOUBLE_PRECISION,proc,TAG,MPI_COMM_WORLD,ierr)
            ELSEIF (StoppingOption1==0) THEN
               CALL MPI_SEND(MDTimePerBCMD,1,MPI_DOUBLE_PRECISION,proc,TAG,MPI_COMM_WORLD,ierr)
            ELSE
               WRITE(6,*) "Err>> StoppingOption1 is having wierd value"
               STOP
            END IF
         ELSE
            WRITE(6,*) "Err>> Message from compute processor was not understood"
            STOP
         END IF
         
         !IF (StoppingOption1==1) THEN !the number of processors is less than number 
         !   !of searches performed and stop time has to be found
         !   charstateindx=INT2CHAR(stateindex)
         !   charisearch=INT2CHAR(isearchcompleted)
         !   filename="./GlobalAMD/MD."//TRIM(charstateindx)//"."//TRIM(charisearch)//".param"
         !   OPEN(UNIT=314,FILE=filename)
         !   READ(314,*) ForcedBCMDTime !time elapsed
         !   CLOSE(314)
         !   !IF (StoppingOption1==1) TotalMDTimeWithStopOption1=ForcedBCMDTime
            
         !   OPEN(314,FILE="MDTransitionSearchStatus1") !presently the MD does not read the file regularly. 
         !      !we need to implement this
         !   WRITE(314,*) 0._dp !ForcedBCMDTime
         !   CLOSE(314)
         !   CALL SYSTEM("mv MDTransitionSearchStatus1 MDTransitionSearchStatus")
         !   StoppingOption1=0  !I am removing this since this serves no purpose
         !END IF
      END DO
!CLOSE(377)
      
      !IF (TotalMDTimeWithStopOption1>0._dp) TotalMDTime=TotalMDTimeWithStopOption1 !this is analogous to parrep
      CALL ReleaseProc(jobid)
      !DO proc=1,num_procs !send negative number as a signal to stop
      !   CALL MPI_SEND(-1,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,ierr)
      !END DO
      
      !CALL MPI_FINALIZE(ierr)
      !IF (ierr/=0) THEN
      !   WRITE(6,*) "$Err>> MPI_MD_DOMAINDECOMPOSITION failed to finalize MPI"
      !END IF
      
   END SUBROUTINE PerformMDProcessSearchParallel
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PerformLocalMDProcessSearchParallel(nfiles,nsearch)!state,NKMCMove,Time,nsearch,Temperature, &
      !DetectTransitionFrequency,WriteTrajectoryFile)
   !Used for searching processes while using LEKMC-MD with domains
   
      IMPLICIT NONE
      INTEGER :: nfiles
      INTEGER, DIMENSION(nfiles) :: nsearch
      
   END SUBROUTINE PerformLocalMDProcessSearchParallel
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SpatiallyParallelMD(MDSiml,nmd,Domains,WriteXYZTrajectoryFile,WriteTrajectoryFrequency)
      IMPLICIT NONE
      TYPE(MDContainer), POINTER :: MDSiml
      INTEGER, DIMENSION(3) :: Domains
      INTEGER :: num_procs1,irank,WriteTrajectoryFrequency,jobid,jobid1,nproc_available
      INTEGER :: nmd,proc,iproc,istatus
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status=0
      CHARACTER(len=100) :: WriteXYZTrajectoryFile
      REAL(dp) :: dt
      
      jobid=4 !Spatially parallel MD is jobid 4
      num_procs1=PRODUCT(Domains)
      nproc_available=GetNumProcs(jobid=0,procstatus=0)
      IF (nproc_available<num_procs1) THEN
         WRITE(6,*) "MPI err>> Number of processors should be same or more than number ..."
         WRITE(6,*) "... of domains while performing spatially parallel MD"
         WRITE(UNIT=6,FMT='("Number of processors available for MD:",I5)') nproc_available
         WRITE(UNIT=6,FMT='("Number of domains present:",I5)') num_procs1
         STOP
      END IF
      CALL ReserveProc(num_procs1,jobid=jobid,istatus=istatus)
      IF (istatus==1) THEN
         WRITE(6,*) "Unable to obtain any processors for mpi_md"
         STOP
      END IF
      
      !Secure the processors
      WRITE(6,*) "Reserving processors for spatially parallel MD ..."
      DO iproc=1,num_procs1
         proc=GetFreeProc(jobid=jobid)
         WRITE(6,*) "Reserved processor: ",proc
         CALL FLUSH(6)
         IF (proc>0) THEN
            TAG=0
            CALL MPI_SEND(jobid,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,ierr)
            IF (ierr/=0) THEN
               WRITE(6,*) "Err>> Error in sending message"
               STOP
            END IF
            CALL MPI_RECV(jobid1,1,MPI_INTEGER,proc,TAG,MPI_COMM_WORLD,status,ierr)
            IF (ierr/=0) THEN
               WRITE(6,*) "Err>> Error in receiving message"
               STOP
            END IF
            IF (jobid1/=jobid) THEN
               WRITE(6,*) "MPI is not working"
               WRITE(6,*) "Found jobid as",jobid
               STOP
            END IF
         ELSE
            WRITE(6,*) "Err>> Unable to secure processor for spatially parallized MD"
            STOP
         END IF
         CALL SetProcStatus(proc=proc,procstatus=1,jobid=jobid)
      END DO
      
      CALL PrintProcStatus()
      
      !Perform MD
      dt=MDSiml%dt
      CALL md_mpi_master(MDSiml,nmd,Domains,num_procs1,0,dt,WriteXYZTrajectoryFile,WriteTrajectoryFrequency)
      
      CALL ReleaseProc(jobid)
      
   END SUBROUTINE SpatiallyParallelMD
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SlaveTaskManager(irank)
      IMPLICIT NONE
      INTEGER :: irank,jobid,status=0,i
      INTEGER :: ifile1  !for NEB
      INTEGER :: isearch1,stateindex,DetectTransitionFrequency,iseed_amd !for MD
      INTEGER :: msgtype,ObtainRate
      REAL(dp) :: Temperature,MDTime !for MD
      REAL(sp) :: ran1scl
      REAL :: ActivationEnergy,MaxDisplacement
      LOGICAL :: NotSatisfied,success
      
      NotSatisfied=.TRUE.
      
      DO WHILE (NotSatisfied)
         
         TAG=0
         CALL MPI_RECV(jobid,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,status,ierr)
         IF (ierr/=0) THEN
            WRITE(6,*) "Err>> Error in recieving message"
            STOP
         END IF
         !CALL MPI_SEND(jobid,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,ierr)
         IF (jobid==1) THEN  !NEB based search for KMC processes
      
            CALL MPI_RECV(ifile1,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,status,ierr)
            CALL ProcessNEB(MinDisplacement=0.,ActivationEnergy=ActivationEnergy, &
               MaxDisplacement=MaxDisplacement,ifile=ifile1,ReuseVL=.FALSE.)
            WRITE(6,*) "I have completed my job",irank," with ifile ",ifile1
            CALL MPI_SEND(irank,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,ierr) !report job is done
      
         ELSEIF (jobid==2) THEN  !MD based search for KMC processes

            CALL MPI_RECV(isearch1,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(stateindex,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(Temperature,1,MPI_DOUBLE_PRECISION,0,TAG,MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(DetectTransitionFrequency,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,status,ierr)
            CALL MPI_RECV(ObtainRate,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,status,ierr)
            !ObtainRate=1
            CALL MPI_RECV(iseed_amd,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,status,ierr)
            CALL ran_seed(iseed_amd)
            DO i=1,10000
               CALL ran1(ran1scl)
            END DO
            
            CALL SearchProcessesFromState(stateindex,Temperature,DetectTransitionFrequency, &
               isearch1,success,MDTime,ObtainRate)
            
            WRITE(6,*) "I have completed my job",irank," with isearch ",isearch1
            TAG=0
            CALL MPI_SEND(irank,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,ierr) !report job is done
            msgtype=1 !tell job is done
            CALL MPI_SEND(msgtype,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(success,1,MPI_LOGICAL,0,TAG,MPI_COMM_WORLD,ierr) !report job is done
            CALL MPI_SEND(MDTime,1,MPI_DOUBLE_PRECISION,0,TAG,MPI_COMM_WORLD,ierr) !report job is done
            CALL MPI_SEND(isearch1,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,ierr) !report job is done
      
         ELSEIF (jobid==3) THEN  !Dimer based search for KMC processes
         
         ELSEIF (jobid==4) THEN  !spatially parallel MD
         
            !CALL md_mpi()
            
         ELSE !signal to stop
      
            WRITE(6,*) "Shutting down processor:",irank
            NotSatisfied=.FALSE.
      
         END IF
!        CALL FLUSH(6)
      END DO
      
   END SUBROUTINE SlaveTaskManager
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE KMCProcessSearchParallel 
