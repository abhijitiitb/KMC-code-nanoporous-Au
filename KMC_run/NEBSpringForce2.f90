!NEB improved tangent implementation
IF (CI_NEB1 .AND. TSImage==image) THEN
   !do nothing -- since the TS cannot see the springs
ELSE
   t(1:3*NAtoms)=AtomCoordNext(1:3*NAtoms)-AtomCoord(1:3*NAtoms)
   normt=SpringConst*SQRT(DOT_PRODUCT(t(1:3*NAtoms),t(1:3*NAtoms)))
   t(1:3*NAtoms)=AtomCoord(1:3*NAtoms)-AtomCoordPrev(1:3*NAtoms)
   normt=normt-SpringConst*SQRT(DOT_PRODUCT(t(1:3*NAtoms),t(1:3*NAtoms)))
   
   ChOSForce(1:3*NAtoms)=ChOSForce(1:3*NAtoms)+normt*ChOSTangent(1:3*NAtoms)
END IF
