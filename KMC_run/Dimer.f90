MODULE DimerMethod
   !works only with the conjugate-gradient minimizer
   !the minimizer is added here itself since the dimer uses a different type of CG
   !t:tangent along the dimer
   !n:unit vector along the perpendicular force
   !n0: unit vector from the cg minimizer
   !old: appended to indicate value from previous iteration
   
   USE OptimizationDimer
   
   IMPLICIT NONE
   
   CONTAINS
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   SUBROUTINE Dimer(AL,AL1,rmode,tmode,iprint) !returns a saddle AL for the state AL1
      !imode: 1=Modified Newton method
      !imode: 2=CG method
      !imode: 3=LBFGS method
      
      IMPLICIT NONE
      TYPE(SystemContainer), POINTER :: AL,AL1
      TYPE(ChOSContainer), POINTER :: chos
      INTEGER :: NAtoms
      INTEGER, OPTIONAL :: rmode,tmode,iprint
      INTEGER :: iprint1
      
      IF (PRESENT(iprint)) THEN
         iprint1=iprint
      ELSE
         iprint1=0
      END IF
      
      NULLIFY(chos) !ensure ChOS is empty
      NULLIFY(AL1)
      CALL AddChOSDimer(s1=AL,chos=chos) !builds a random chain of 2 states
      CALL OptimizeDimer(chos,rmode,tmode) !align dimer
      
      CALL Delete(chos)
   END SUBROUTINE Dimer
   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END MODULE DimerMethod
