
      s1=AtomSpecies(iatom)
      
      !Compute the embedding energy
      section=s1
      dxinv=IncrEEinv(section)
      mincutoff=RangeEE(2*section-1)
      tmp=(DensityArray(iatom)-mincutoff)*dxinv
      k=FLOOR(tmp)
      IF (k>SizeEE(section)-2) THEN
         WRITE(6,*) TRIM(TxtHeader)//"$Err[0]>> Increase range of density array"
         WRITE(6,*) TRIM(TxtHeader)//"...embedding energy term"
         WRITE(UNIT=6,FMT='("... density encountered:",ES14.5)') Density
         STOP
      END IF
         
      IF (k<=1) CYCLE   !free atom
      
      p=MIN(tmp-REAL(k,dp),1._dp)
      k=MIN(k,SizeEE(section)-1)+(section-1)*MaxEAMTableSize
      EnergyEmbedding=EnergyEmbedding+EE(k)+0.5_dp*p*(EE(k+1)-EE(k-1)  +p*(EE(k+1)-EE(k)-EE(k)+EE(k-1)))
        !embedding energy completed
