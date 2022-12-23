MODULE zdiis_utils
  USE kinds,                           ONLY: real_8
  USE linres,                          ONLY: lr01,&
                                             lr02
  USE mp_interface,                    ONLY: mp_sum
  USE odiis_utils,                     ONLY: solve,&
                                             updis
  USE parac,                           ONLY: parai
  USE system,                          ONLY: maxdis,&
                                             ncpw
  USE td_dav_utils,                    ONLY: dotps,&
                                             setvpp
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: zdiis
  !public :: applpre

CONTAINS

  ! ==================================================================
  SUBROUTINE zdiis(c0,c2,vpp,eigv,nstate,pme,gde,&
       svar2,reinit,orbital)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vpp(:), eigv(*)
    INTEGER                                  :: nstate
    COMPLEX(real_8) :: c2(ncpw%ngw,nstate), c0(ncpw%ngw,nstate), &
      pme(ncpw%ngw*nstate,*), gde(ncpw%ngw*nstate,*)
    REAL(real_8)                             :: svar2
    LOGICAL                                  :: reinit
    CHARACTER(len=*)                         :: orbital

    INTEGER                                  :: i, isub, j, nsize
    INTEGER, SAVE                            :: istate = 0, ndiis, nowv
    REAL(real_8)                             :: bc_array(maxdis+1,maxdis+1), &
                                                vc(maxdis+1), vppmin
    REAL(real_8), SAVE                       :: diism(maxdis,maxdis)

! ==--------------------------------------------------------------==

    CALL tiset('     ZDIIS',isub)
    IF (reinit.OR.istate.EQ.0) THEN
       ! Reinitialization of cntl%diis procedure
       istate=0
       reinit=.FALSE.
       ndiis=0
       ! INIT VPP
       CALL setvpp(vpp,ncpw%ngw,lr02%lr_hthrs,vppmin)
    ENDIF
    IF (istate.EQ.0) CALL updis(ndiis,nowv,nsize,lr01%mldiis,0)
    istate=2
    ! Perform an electronic cntl%diis step
    CALL updis(ndiis,nowv,nsize,lr01%mldiis,1)
    ! Update cntl%diis buffers
    CALL dcopy(2*ncpw%ngw*nstate,c0,1,pme(1,nowv),1)
    ! Scale Gradient with Preconditioner
    CALL dscal(2*ncpw%ngw*nstate,-1.0_real_8,c2,1)
    CALL applpre(c2,vpp,eigv,nstate,orbital)
    CALL dcopy(2*ncpw%ngw*nstate,c2,1,gde(1,nowv),1)
    ! Update cntl%diis matrix
    DO i=1,nsize-1
       diism(i,nowv)=dotps(ncpw%ngw,nstate,gde(:,i:i+nstate-1),gde(:,nowv:nowv+nstate-1))
    ENDDO
    CALL mp_sum(diism(:,nowv),nsize-1,parai%allgrp)
    DO i=1,nsize-1
       diism(nowv,i)=diism(i,nowv)
    ENDDO
    ! Set up cntl%diis Matrix
    CALL zeroing(bc_array)!,(maxdis+1)*(maxdis+1))
    DO i=1,nsize-1
       DO j=1,nsize-1
          bc_array(i,j)=diism(i,j)
       ENDDO
    ENDDO
    DO i=1,nsize-1
       vc(i)=0._real_8
       bc_array(i,nsize)=-1._real_8
       bc_array(nsize,i)=-1._real_8
    ENDDO
    vc(nsize)=-1._real_8
    bc_array(nsize,nsize)=0._real_8
    ! Solve System of Linear Equations
    CALL solve(bc_array,maxdis+1,nsize,vc)
    ! Compute Interpolated Coefficient Vectors
    CALL zeroing(c0)!,SIZE(c0))
    DO i=1,nsize-1
       CALL daxpy(2*ncpw%ngw*nstate,vc(i),pme(1,i),1,c0,1)
    ENDDO
    ! Estimate New Parameter Vectors 
    IF (nsize-1.EQ.1) THEN
       CALL daxpy(2*ncpw%ngw*nstate,-svar2*vc(1),gde,1,c0,1)
    ELSE
       DO i=1,nsize-1
          CALL daxpy(2*ncpw%ngw*nstate,-vc(i),gde(1,i),1,c0,1)
       ENDDO
    ENDIF
    CALL tihalt('     ZDIIS',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zdiis
  ! ==================================================================
  SUBROUTINE applpre(c2,vpp,eigv,nstate,orbital)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: vpp(ncpw%ngw)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: eigv(nstate,nstate)
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate)
    CHARACTER(len=*)                         :: orbital

    INTEGER                                  :: ig, is
    REAL(real_8)                             :: deij

    DO is=1,nstate
       IF (INDEX(orbital,"CANON").NE.0) THEN
          deij=-eigv(is,1)
       ELSE
          deij=eigv(is,is)
       ENDIF
       DO ig=1,ncpw%ngw
          c2(ig,is)=c2(ig,is)/(vpp(ig)+deij)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE applpre
  ! ==================================================================

END MODULE zdiis_utils
