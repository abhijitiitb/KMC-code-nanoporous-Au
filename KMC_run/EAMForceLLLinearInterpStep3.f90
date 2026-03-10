!for the pair iatom-jatom find the pair and density terms

dr=AtomCoord(3*jatom-2:3*jatom)-AtomCoord(3*iatom-2:3*iatom)
dr=dr-BoxSize*NINT(dr/BoxSize) !PBC
drmag=SQRT( dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3) )

IF (drmag<=cutoff) THEN

s2=AtomSpecies(jatom)

!xxxxxxxxxxxxxxxxxxxx  PAIR POTENTIAL TERM  xxxxxxxxxxxxxxxxxxx
smin=MIN(s1,s2)
smax=MAX(s1,s2)

section=((2*NSpecies-smin)*(smin-1))/2+smax
dxinv=IncrPPinv(section)
mincutoff=RangePP(2*section-1)

tmp=(drmag-mincutoff)*dxinv
k0=FLOOR(tmp) !table index for reading density and pair potential, use this instead of abi
IF (k0<2) THEN
   IF (CatchOutOfBound) THEN
      WRITE(6,*) TRIM(TxtHeader)//"$Err>> EAM index below PP array lower bound"
      WRITE(UNIT=6,FMT='("...too low drmag:",ES14.5)') drmag
      WRITE(6,*) "Index in interp table=", k
      WRITE(UNIT=6,FMT='("Main atom index:",I4,"           Main atom coord:",3f10.3)') &
         iatom,AtomCoord(3*iatom-2:3*iatom)
      errorstatus=1  !pair potential is below range
      k0=2
      STOP
      !RETURN
   ELSE
      k=(section-1)*MaxEAMTableSize+1
      !linear extrapolation
      EnergyPP=EnergyPP+ PP(k) + tmp*(PP(k+1)-PP(k))
      tgphi=(PP(k+1)-PP(k))*dxinv/drmag !dphi/drij/drmag
      !WRITE(6,*) TRIM(TxtHeader)//"S3-1:",tgphi,section,drmag
   END IF
ELSE
   p=MIN(tmp-REAL(k0,dp),1._dp)
   k=MIN(k0,SizePP(section)-1)+(section-1)*MaxEAMTableSize
   EnergyPP=EnergyPP+PP(k)+ 0.5_dp*p*(PP(k+1)-PP(k-1)  +p*(PP(k+1)-PP(k)-PP(k)+PP(k-1)))
   !note it is assumed that jatom is part of the Domain i.e.,VLListDomainAtom(j)=.TRUE.

   !gradient from the phi term
   k=MIN(k0,SizePP(section)-2)+(section-1)*MaxEAMTableSize
   gi=PP(k+1)-PP(k-1)
   gi1=PP(k+2)-PP(k)
   slope=0.5_dp*(gi+p*(gi1-gi))*dxinv !gradient of the pair potential between iatom and jatom
   tgphi=slope/drmag
END IF

!Density term -- role of jatom on iatom
section=s2
dxinv=IncrDTinv(section)
mincutoff=RangeDT(2*section-1)
tmp=(drmag-mincutoff)*dxinv
k0=FLOOR(tmp)

tgrho=0._dp
IF (k0<2) THEN
   !IF (CatchOutOfBound) THEN !this must have already been checked in Step1.f90
   
   k=(section-1)*MaxEAMTableSize+1
   !linear extrapolation
   tgrho=(DT(k+1)-DT(k))*dxinv/drmag !drho/drij/drmag
   !WRITE(6,*) TRIM(TxtHeader)//"S3-2:",tgrho,section,drmag
ELSEIF (drmag<RangeDT(2*section)) THEN
   p=MIN(tmp-REAL(k0,dp),1._dp)
   
   k=MIN(k0,SizeDT(section)-2)+(section-1)*MaxEAMTableSize
   gi=DT(k+1)-DT(k-1)
   gi1=DT(k+2)-DT(k)
   slope=0.5_dp*(gi+p*(gi1-gi))*dxinv
   tgrho=slope/drmag
END IF

tgrad=tgphi+DensityArray1(jatom)*tgrho
   !use of jatom here requires some discussion. if jatom lies inside the domain then the forces are anyways correct
   !however if jatom is outside the domain, then we have not computed its density. because of the periodic boundary conditions
   !we can take density at jatom atom as the same as the jatom in the domain
   !grad=tgrad*VLdr(3*j-2:3*j)
   
!Density term -- role of iatom on jatom
section=s1
dxinv=IncrDTinv(section)
mincutoff=RangeDT(2*section-1)
tmp=(drmag-mincutoff)*dxinv
k0=FLOOR(tmp)

tgrho=0._dp
IF (k0<2) THEN
   !IF (CatchOutOfBound) THEN !this must have already been checked in Step1.f90
   
   k=(section-1)*MaxEAMTableSize+1
   !linear extrapolation
   tgrho=(DT(k+1)-DT(k))*dxinv/drmag !drho/drij/drmag
   !WRITE(6,*) TRIM(TxtHeader)//"S3-3:",tgrho,section,drmag
ELSEIF (drmag<RangeDT(2*section)) THEN
   p=MIN(tmp-REAL(k0,dp),1._dp)

   k=MIN(k0,SizeDT(section)-2)+(section-1)*MaxEAMTableSize
   gi=DT(k+1)-DT(k-1)
   gi1=DT(k+2)-DT(k)
   slope=0.5_dp*(gi+p*(gi1-gi))*dxinv
   tgrho=slope/drmag
END IF

tgrad=tgrad+DensityArray1(iatom)*tgrho

grad=tgrad*dr

AtomForce(3*iatom-2:3*iatom)=AtomForce(3*iatom-2:3*iatom)+grad

END IF
