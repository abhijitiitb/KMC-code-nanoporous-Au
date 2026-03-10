
      ReadClsman=(errorstatus==0)
      IF (ReadClsmanPrevImage) THEN !at least some images were read so far
         IF (.NOT. ReadClsman) THEN
            ReadClsman=.TRUE. !it was successful since we were able to read images
            WRITE(UNIT=6,FMT='(" ...number of images read ",I4)') nimages
            CLOSE(UnitTmp)
            RETURN
         END IF
      ELSE
         IF (.NOT. ReadClsman) THEN
            CLOSE(UnitTmp)
            RETURN
         END IF
      END IF
