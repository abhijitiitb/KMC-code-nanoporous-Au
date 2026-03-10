MODULE mpimanager
!manages the nodes for job allocation
!NProcs is the number of processors
!ProcStatus is 0 or 1 depending on whether it is free or not
!ProcJobType can 0 implying it has not been assigned a job type, if a non-zero value is present
  !it implies that the job 
   IMPLICIT NONE
   INTEGER, PARAMETER :: MaxProcs=1000
   INTEGER, DIMENSION(MaxProcs) :: ProcessorStatus=0,ProcessorJobType=0
   INTEGER :: num_procs=0
   INTEGER :: ierr,irank,TAG=4
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetNumProcs(nprocs)
   !sets the total number of processors that are available
      IMPLICIT NONE
      INTEGER :: nprocs
      
      IF (nprocs>MaxProcs) THEN
         WRITE(6,*) "Err>> Number of processors exceeds MaxProcs, increase MaxProcs"
         STOP
      END IF
      
      num_procs=nprocs !these are the number of processors available
   END SUBROUTINE SetNumProcs
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetNumProcs(jobid,procstatus)
   !number of processor with jobid=jobid and status=procstatus
      IMPLICIT NONE
      INTEGER :: jobid,procstatus,i,GetNumProcs
      
      GetNumProcs=0
      
      IF (num_procs==0) THEN
         WRITE(6,*) "$Err>> In mpimanager number processors is set to zero"
         STOP
      END IF
      
      DO i=1,num_procs
         IF (ProcessorJobType(i)==jobid .AND. ProcessorStatus(i)==procstatus) GetNumProcs=GetNumProcs+1
      END DO
   END FUNCTION GetNumProcs
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ResetProcStatus()
      IMPLICIT NONE
      
      ProcessorJobType(1:num_procs)=0
      ProcessorStatus(1:num_procs)=0
      WRITE(6,*) "Number of slave processors:",num_procs
      
   END SUBROUTINE ResetProcStatus
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetMPINumberBusyProcs(jobid)
      IMPLICIT NONE
      INTEGER :: GetMPINumberBusyProcs,jobid,i
      
      GetMPINumberBusyProcs=0
      DO i=1,num_procs
         IF (ProcessorJobType(i)==jobid .AND. ProcessorStatus(i)/=0) THEN
            GetMPINumberBusyProcs=GetMPINumberBusyProcs+1
         END IF
      END DO
   END FUNCTION GetMPINumberBusyProcs
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   FUNCTION GetFreeProc(jobid)
   !obtains a free processor of JobType jobid
      IMPLICIT NONE
      INTEGER :: GetFreeProc,jobid,i
      
      GetFreeProc=0
      
      IF (num_procs==0) THEN
         WRITE(6,*) "$Err>> In mpimanager number processors is set to zero"
         STOP
      END IF
      
      DO i=1,num_procs
         IF (ProcessorJobType(i)==jobid .AND. ProcessorStatus(i)==0) THEN
            GetFreeProc=i
            RETURN
         END IF
      END DO
   END FUNCTION GetFreeProc
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE SetProcStatus(proc,procstatus,jobid)
   !sets a processor to processor status and assigns the jobid
      IMPLICIT NONE
      INTEGER :: proc,procstatus,jobid
      
      ProcessorStatus(proc)=procstatus
      ProcessorJobType(proc)=jobid
   END SUBROUTINE SetProcStatus
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReserveProc(nresv,jobid,istatus)
   !reserves nresv processor of job type jobid
      IMPLICIT NONE
      INTEGER :: nresv,jobid,istatus,i,proc,nresv1
      
      istatus=0
      nresv1=0
      DO i=1,nresv
         proc=GetFreeProc(jobid=0) !proc with no allocation
         IF (proc==0) THEN !failed to obtain a free processor
            IF (nresv1==0) THEN
               istatus=1 !failed completely
            ELSEIF (nresv1>0) THEN
               istatus=2
            END IF
            WRITE(UNIT=6,FMT='(">>Failed to reserve required number of processors ",I5," for jobid ",I3)') & 
               nresv,jobid
            WRITE(UNIT=6,FMT='(">>Instead reserving only ",I5," processors for jobid ",I3)') nresv1,jobid
            RETURN
         ELSE
            CALL SetProcStatus(proc,procstatus=0,jobid=jobid)
            nresv1=nresv1+1
         END IF
      END DO
   END SUBROUTINE ReserveProc
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE ReleaseProc(jobid)
   !finds all processors with jobid given by jobid and sets their jobid to 0
      IMPLICIT NONE
      INTEGER :: jobid,istatus,i
      
      istatus=1
      DO i=1,num_procs
         IF (ProcessorJobType(i)==jobid) THEN
            ProcessorStatus(i)=0
            ProcessorJobType(i)=0
         END IF
      END DO
      istatus=0
      
   END SUBROUTINE ReleaseProc
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE PrintProcStatus(proc)
   !Prints the status of the processor
      IMPLICIT NONE
      INTEGER, OPTIONAL :: proc
      INTEGER :: proc1
      
      WRITE(6,*) "Printing status of the processors ..."
      WRITE(6,*) "Processor      Status    Jobid"
      WRITE(6,*) "---------      ------    -----"
      WRITE(6,*) "   0              1      Master"
      IF (PRESENT(proc)) THEN
         WRITE(UNIT=6,FMT='(I4,I5)',ADVANCE="NO") proc,ProcessorStatus(proc)
         IF (ProcessorStatus(proc)>0) THEN
            WRITE(UNIT=6,FMT='(I6)') ProcessorJobType(proc)
         ELSE
            WRITE(6,*) ""
         END IF
      ELSE
         DO proc1=1,num_procs
            WRITE(UNIT=6,FMT='(I5,"          ",I5)',ADVANCE="NO") proc1,ProcessorStatus(proc1)
            IF (ProcessorStatus(proc1)>0 .OR. ProcessorJobType(proc1)/=0) THEN
               WRITE(UNIT=6,FMT='("   ",I6)') ProcessorJobType(proc1)
            ELSE
               WRITE(UNIT=6,FMT='("        -")')
            END IF
         END DO
      END IF
   END SUBROUTINE PrintProcStatus
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE mpimanager
