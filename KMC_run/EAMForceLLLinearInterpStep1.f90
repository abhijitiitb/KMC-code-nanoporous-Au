!for the pair iatom-jatom find the density terms

dr=AtomCoord(3*jatom-2:3*jatom)-AtomCoord(3*iatom-2:3*iatom)
dr=dr-BoxSize*NINT(dr/BoxSize) !PBC
drmag=SQRT( dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) )

IF (drmag<=cutoff) THEN

s2=AtomSpecies(jatom)

!xxxxxxxxxxxxxxxxxxxx  DENSITY TERM  xxxxxxxxxxxxxxxxxxx
!contribution of jatom to the density at iatom
section=s2 !density of species 2 at location of iatom
dxinv=IncrDTinv(section)
mincutoff=RangeDT(2*section-1)

tmp=(drmag-mincutoff)*dxinv
k0=FLOOR(tmp) !table index for reading density and pair potential, use this instead of abi
IF (k0<2) THEN
   IF (CatchOutOfBound) THEN
      WRITE(6,*) TRIM(TxtHeader)//"$Err>> EAM index below DT array lower bound"
      WRITE(UNIT=6,FMT='("...too low drmag:",2ES14.5)') drmag,mincutoff
      WRITE(6,*) "Index in interp table=", k0
      WRITE(UNIT=6,FMT='("Main atom coord:",3f15.3)') AtomCoord(3*iatom-2:3*iatom)
      errorstatus=2
      k0=2
      CALL FLUSH(6)
      STOP
   ELSE
      !we shall (linear) extrapolate -- note that k<2
      !use tmp to find the function value - 11/07/12
      !note we don't need to increment tmp by 1 - 11/07/12
      !let t denote the (real/fractional) array position corresponding to tmp then --
      ! (f(t)-f(1))/(t-1)=f(2)-f(1) - where f(1) is fn value at mincutoff
      !or f(t)=f(1)+tmp*(f(2)-f(1)) - 11/07/12
      !more generally f(t)=f(1)+(tmp+1+W-1-W)*(f(2)-f(1)) -- where W=(section-1)*MaxEAMTableSize
      k=(section-1)*MaxEAMTableSize+1
      DensityArray(iatom)=DensityArray(iatom)+ DT(k) + tmp*(DT(k+1)-DT(k))
      !WRITE(6,*) TRIM(TxtHeader)//"S1-1:",tmp,DT(k) + tmp*(DT(k+1)-DT(k)),DT(k),DT(k+1),section,drmag
   END IF
ELSE

   p=MIN(tmp-REAL(k0,dp),1._dp)
   k=MIN(k0,SizeDT(section)-1)+(section-1)*MaxEAMTableSize
   DensityArray(iatom)=DensityArray(iatom)+DT(k)+0.5_dp*p*(DT(k+1)-DT(k-1)  +p*(DT(k+1)-DT(k)-DT(k)+DT(k-1)))
   
END IF

END IF

