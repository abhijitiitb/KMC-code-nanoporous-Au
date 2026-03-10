
      IF (.NOT. ASSOCIATED(AL)) THEN
         WRITE(6,*) "$Err>> AL is not associated in WriteXYZ"
         STOP
      END IF
      
      IF (PRESENT(fileopen)) THEN
         fileopen1=fileopen
      ELSE
         fileopen1=.TRUE. !file needs to be created and opened if TRUE
      END IF
      
      IF (fileopen1) THEN
         iunit1=301
         IF (.NOT. (PRESENT(filename))) THEN
            WRITE(6,*) "$Err>> Filename should have been provided"
            STOP
         END IF
      ELSE
         IF (PRESENT(iunit)) THEN
            iunit1=iunit
         ELSE
            WRITE(6,*) "$Err>> A unit for the file should be specified"
            STOP
         END IF
      END IF
      
      IF (PRESENT(fileclose)) THEN
         fileclose1=fileclose
      ELSE
         fileclose1=.TRUE.
      END IF
      
      IF (PRESENT(iprint)) THEN
         iprint1=iprint
      ELSE
         iprint1=1
      END IF
      
      IF (fileopen1) OPEN(UNIT=iunit1,FILE=TRIM(filename))
