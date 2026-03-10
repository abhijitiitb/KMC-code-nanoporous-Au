
      ReadXYZ=(errorstatus==0)
      IF (ReadXYZPrevImage) THEN
         IF (.NOT. ReadXYZ) THEN
            ReadXYZ=.TRUE. !we were able to read some images
            WRITE(UNIT=6,FMT='(" ...number of images read ",I4)') nimages
            CLOSE(UnitTmp)
            RETURN
         END IF
      ELSE
         IF (.NOT. ReadXYZ) THEN
            CLOSE(UnitTmp)
            RETURN
         END IF
      END IF
