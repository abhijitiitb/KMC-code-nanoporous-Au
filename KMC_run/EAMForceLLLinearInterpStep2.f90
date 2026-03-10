 
!Compute the embedding energy
section=s1
dxinv=IncrEEinv(section)
mincutoff=RangeEE(2*section-1)

tmp=(DensityArray(iatom)-mincutoff)*dxinv
k0=FLOOR(tmp)
      
IF (k0>1) THEN  !k0<=1 implies free atom  -- the steps that follow are not needed for a free atom

IF (k0>SizeEE(section)-2) THEN
   IF (CatchOutOfBound) THEN
      WRITE(6,*) "$Err[0]>> Increase range of density array"
      WRITE(UNIT=6,FMT='("... density encountered:",ES14.5)') DensityArray(iatom)
      errorstatus=2
      k0=2
      STOP
   ELSE
      !perform linear extrapolation -- added 11/07/12 by AC
      !use tmp to find the function value - 11/07/12
      !let t=tmp+1 denote the (real/fractional) array position corresponding to tmp then --
      ! (f(t)-f(e))/(t-e)=f(e)-f(e-1) - where f(e) is fn value at cutoff
      !or f(t)=f(e)+(tmp+1)*(f(e)-f(e-1)) - 11/07/12
      k=(section-1)*MaxEAMTableSize+SizeEE(section) !correct position for section value
      !move tmp to the correct location in the array given by tmp+1.+(section-1)*MaxEAMTableSize
      EnergyEmbedding=EnergyEmbedding+ EE(k)+(tmp+1._dp+REAL((section-1)*MaxEAMTableSize,dp)-k)* &
         (EE(k)-EE(k-1))
      DensityArray1(iatom)=(EE(k)-EE(k-1))*dxinv !dF/drho
      !WRITE(6,*) TRIM(TxtHeader)//"S2-1:",(EE(k)-EE(k-1))*dxinv,EE(k),EE(k-1),k,section,DensityArray(iatom)
   END IF
ELSE
   
   p=MIN(tmp-REAL(k0,dp),1._dp)
   k=MIN(k0,SizeEE(section)-1)+(section-1)*MaxEAMTableSize
   EnergyEmbedding=EnergyEmbedding+EE(k)+0.5_dp*p*(EE(k+1)-EE(k-1)  +p*(EE(k+1)-EE(k)-EE(k)+EE(k-1)))

   !compute the derivative
   k=MIN(k0,SizeEE(section)-2)+(section-1)*MaxEAMTableSize
   gi=EE(k+1)-EE(k-1)
   gi1=EE(k+2)-EE(k)
   DensityArray1(iatom)=0.5_dp*(gi+p*(gi1-gi))*dxinv  !dF/drho for atom i

END IF

END IF
