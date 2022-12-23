MODULE radin_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: radin
  PUBLIC :: radlg

CONTAINS

  ! ==================================================================
  SUBROUTINE radin(mesh,c,func,asum)
    ! ==--------------------------------------------------------------==
    ! == SIMPSONS RULE INTEGRATION FOR HERMAN SKILLMAN MESH           ==
    ! == MESH - # OF MESH POINTS                                      ==
    ! == C    - 0.8853418/Z**(1/3.)                                   ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: mesh
    REAL(real_8)                             :: c, func(mesh), asum

    INTEGER                                  :: i, i1, j, k, nblock
    REAL(real_8)                             :: a1, a2e, a2es, a2o, h

    a1=0.0_real_8
    a2e=0.0_real_8
    asum=0.0_real_8
    h=0.0025_real_8*c
    nblock=mesh/40
    i=1
    func(1)=0.0_real_8
    DO j=1,nblock
       DO k=1,20
          i=i+2
          i1=i-1
          a2es=a2e
          a2o=func(i1)/12._real_8
          a2e=func(i)/12._real_8
          a1=a1+5._real_8*a2es+8._real_8*a2o-a2e
          func(i1)=asum+a1*h
          a1=a1-a2es+8._real_8*a2o+5._real_8*a2e
          func(i)=asum+a1*h
       ENDDO
       asum=func(i)
       a1=0._real_8
       h=h+h
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE radin
  ! ==================================================================
  SUBROUTINE radlg(mesh,func,rab,asum)
    ! ==--------------------------------------------------------------==
    ! == SIMPSONS RULE INTEGRATOR FOR FUNCTION STORED ON THE          ==
    ! == RADIAL LOGARITHMIC MESH                                      ==
    ! == RAB(MESH) : LOGARITHMIC RADIAL MESH INFORMATION              ==
    ! == FUNC(MESH): FUNCTION TO BE INTEGRATED                        ==
    ! ==--------------------------------------------------------------==
    ! == ROUTINE ASSUMES THAT MESH IS AN ODD NUMBER SO RUN CHECK      ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: mesh
    REAL(real_8)                             :: func(mesh), rab(mesh), asum

    INTEGER                                  :: i
    REAL(real_8)                             :: f1, f2, f3, r12

    IF ( mesh - ( mesh / 2 ) * 2 .NE. 1 ) THEN
       IF (paral%io_parent)&
            WRITE(6,*) '***ERROR IN SUBROUTINE RADLG'
       IF (paral%io_parent)&
            WRITE(6,*) 'ROUTINE ASSUMES MESH IS ODD BUT MESH =',mesh
       CALL stopgm('RADLG',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    asum = 0.0_real_8
    r12 = 1.0_real_8 / 12.0_real_8
    f3  = func(1) * rab(1) * r12
    func(1) = 0.0_real_8
    ! ==--------------------------------------------------------------==
    DO i = 2,mesh-1,2
       f1 = f3
       f2 = func(i) * rab(i) * r12
       f3 = func(i+1) * rab(i+1) * r12
       asum = asum + 5.0_real_8*f1 + 8.0_real_8*f2 - 1.0_real_8*f3
       func(i) = asum
       asum = asum - 1.0_real_8*f1 + 8.0_real_8*f2 + 5.0_real_8*f3
       func(i+1) = asum
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE radlg
  ! ==================================================================

END MODULE radin_utils
